#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/vector_map.hpp"
#include "common/utils/template_utils.hpp"
#include "common/utils/string_utils.hpp"

namespace mtg {
namespace graph {
namespace align {


class ILabeledDBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::QueryGenerator QueryGenerator;

    // undefined target column
    static constexpr uint64_t kNTarget = std::numeric_limits<uint64_t>::max() - 1;

    ILabeledDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config);

    virtual ~ILabeledDBGAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    typedef std::pair<std::vector<node_index> /* forward */,
                      std::vector<node_index> /* reverse complement */ > Mapping;
    typedef std::pair<sdsl::bit_vector /* forward */,
                      sdsl::bit_vector /* reverse complement */ > Signature;
    typedef std::vector<Mapping> BatchMapping;
    typedef std::vector<VectorMap<uint64_t, Signature>> BatchLabels;

    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    std::pair<BatchMapping, BatchLabels>
    map_and_label_query_batch(const QueryGenerator &generate_query) const;
};

template <class BaseSeeder>
class LabeledSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::node_index node_index;
    typedef typename BaseSeeder::Seed Seed;

    template <typename... Args>
    LabeledSeeder(uint64_t target_column, Args&&... args)
          : BaseSeeder(std::forward<Args>(args)...), target_column_(target_column) {}

    virtual ~LabeledSeeder() {}

    virtual std::vector<Seed> get_seeds() const override {
        std::vector<Seed> seeds = BaseSeeder::get_seeds();
        for (Seed &seed : seeds) {
            assert(seed.get_offset() || target_column_ != ILabeledDBGAligner::kNTarget);
            seed.target_column = !seed.get_offset()
                ? target_column_
                : ILabeledDBGAligner::kNTarget;
        }

        return seeds;
    }

  private:
    uint64_t target_column_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class BaseSeeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledDBGAligner : public ILabeledDBGAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    template <typename... Args>
    LabeledDBGAligner(Args&&... args) : ILabeledDBGAligner(std::forward<Args>(args)...) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

  protected:
    typedef LabeledSeeder<BaseSeeder> Seeder;

    Seeder build_seeder(uint64_t target_column,
                        std::string_view query,
                        bool is_reverse_complement,
                        std::vector<node_index>&& base_nodes,
                        const sdsl::bit_vector &signature) const;
};

template <typename NodeType>
class LabeledColumnExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::score_t score_t;
    typedef typename IExtender<NodeType>::node_index node_index;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~LabeledColumnExtender() {}

    virtual void initialize(const DBGAlignment &path) override;

  protected:
    typedef std::vector<std::pair<NodeType, char>> Edges;
    typedef typename DefaultColumnExtender<NodeType>::Column Column;
    typedef typename DefaultColumnExtender<NodeType>::Scores Scores;
    typedef typename DefaultColumnExtender<NodeType>::AlignNode AlignNode;
    typedef typename DefaultColumnExtender<NodeType>::AlignNodeHash AlignNodeHash;

    virtual Edges get_outgoing(const AlignNode &node) const override;

    virtual void add_scores_to_column(Column &column,
                                      Scores&& scores,
                                      const AlignNode &node) override {
        const AlignNode &prev = std::get<6>(scores);
        assert(align_node_to_target_.count(prev));

        if (!align_node_to_target_.count(node)) {
            align_node_to_target_[node] = align_node_to_target_[prev];
#ifndef NDEBUG
        } else {
            assert(align_node_to_target_[prev] == align_node_to_target_[node]
                || target_columns_.at(align_node_to_target_[prev])
                    == ILabeledDBGAligner::kNTarget);
#endif
        }

        DefaultColumnExtender<NodeType>::add_scores_to_column(
            column, std::move(scores), node
        );
    }

    virtual void backtrack(score_t min_path_score,
                           AlignNode best_node,
                           tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
                           std::vector<DBGAlignment> &extensions) const override {
        size_t old_size = extensions.size();
        DefaultColumnExtender<NodeType>::backtrack(min_path_score, best_node,
                                                   prev_starts, extensions);
        assert(extensions.size() - old_size <= 1);
        if (extensions.size() > old_size) {
            extensions.back().target_column = std::min(
                extensions.back().target_column,
                ILabeledDBGAligner::kNTarget
            );
            assert(extensions.back().get_offset()
                || align_node_to_target_.count(best_node));
            assert(extensions.back().target_column == ILabeledDBGAligner::kNTarget
                || extensions.back().target_column
                    == target_columns_[align_node_to_target_[best_node]]);
            if (extensions.back().target_column == ILabeledDBGAligner::kNTarget) {
                extensions.back().target_column
                    = target_columns_[align_node_to_target_[best_node]];
            }
            assert(extensions.back().get_offset()
                || extensions.back().target_column != ILabeledDBGAligner::kNTarget);
        }
    }

  private:
    const AnnotatedDBG &anno_graph_;
    mutable std::vector<uint64_t> target_columns_;
    mutable tsl::hopscotch_map<NodeType, tsl::hopscotch_map<size_t, Edges>> cached_edge_sets_;
    mutable tsl::hopscotch_map<AlignNode, size_t, AlignNodeHash> align_node_to_target_;
};


template <class BaseSeeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    typedef SeedAndExtendAlignerCore<AlignmentCompare> AlignerCore;

    size_t i = 0;
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        assert(config_.num_alternative_paths);
        AlignerCore main_aligner_core(graph_, config_, query, is_reverse_complement);
        tsl::hopscotch_map<uint64_t, AlignerCore> labeled_aligner_cores;

        DBGQueryAlignment &paths = main_aligner_core.get_paths();
        std::string_view this_query = paths.get_query(is_reverse_complement);
        std::string_view reverse = paths.get_query(true);
        assert(this_query == query);

        Extender extender(anno_graph_, config_, this_query);
        Extender extender_rc(anno_graph_, config_, reverse);

        // const auto &[query_nodes_pair, target_columns] = mapped_batch;
        // assert(target_columns[i].size());

        // const auto &[nodes, nodes_rc] = query_nodes_pair[i];

        // for (const auto &[target_column, signature_pair] : target_columns[i]) {
        {
            uint64_t target_column = ILabeledDBGAligner::kNTarget;
            std::string target_label;
            if (config_.label.size()) {
                const auto &label_encoder = anno_graph_.get_annotation().get_label_encoder();
                if (config_.label == "HEADER") {
                    std::string test(header);
                    test = utils::split_string(test, ".")[0];
                    if (label_encoder.label_exists(test)) {
                        common::logger->trace("Found label {}", test);
                        target_column = label_encoder.encode(test);
                        target_label = test;
                    } else {
                        common::logger->trace("Label {} not found", test);
                    }
                } else {
                    target_column = label_encoder.encode(config_.label);
                    target_label = config_.label;
                }
            }
            common::logger->trace("Target column {}", target_column);
            assert(target_column <= ILabeledDBGAligner::kNTarget);
            AlignerCore &aligner_core = labeled_aligner_cores.emplace(
                target_column,
                AlignerCore { graph_, config_, paths.get_query_ptr(false),
                              paths.get_query_ptr(true)
                }
            ).first.value();

            std::vector<node_index> nodes = map_sequence_to_nodes(graph_, query);
            std::vector<node_index> nodes_rc;

            // const auto &[signature, signature_rc] = signature_pair;
            sdsl::bit_vector signature;
            sdsl::bit_vector signature_rc;

            if (target_label.size()) {
                for (auto &[label, t_signature] : anno_graph_.get_top_label_signatures(query, anno_graph_.get_annotation().num_labels())) {
                    if (label == target_label) {
                        signature = t_signature;
                        break;
                    }
                }
            }

            if (signature.empty()) {
                callback(header, std::move(paths));
                ++i;
                return;
            }

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                    || config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);
                std::string dummy(query);
                nodes_rc = nodes;
                reverse_complement_seq_path(graph_, dummy, nodes_rc);
                assert(dummy == aligner_core.get_paths().get_query(true));
                assert(nodes_rc.size() == nodes.size());

                if (target_label.size()) {
                    for (auto &[label, t_signature] : anno_graph_.get_top_label_signatures(dummy, anno_graph_.get_annotation().num_labels())) {
                        if (config_.label.size() && label == target_label) {
                            signature_rc = t_signature;
                            break;
                        }
                    }
                }

                if (signature_rc.empty()) {
                    callback(header, std::move(paths));
                    ++i;
                    return;
                }
            }


            Seeder seeder = build_seeder(
                target_column,
                this_query, // use this_query since paths stores a copy
                is_reverse_complement,
                std::vector<node_index>(nodes),
                signature
            );

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
                assert(!is_reverse_complement);

                auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                         const auto &callback) {
#ifndef NDEBUG
                    for (const auto &a : rev_comp_seeds) {
                        assert(a.get_offset()
                            || a.target_column < ILabeledDBGAligner::kNTarget);
                    }
#endif
                    callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
                };

                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                Seeder seeder_rc = build_seeder(target_column,
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_both_directions(seeder, seeder_rc,
                                                   extender, extender_rc,
                                                   build_rev_comp_alignment_core);

            } else if (config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);

                Seeder seeder_rc = build_seeder(target_column,
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_best_direction(seeder, seeder_rc, extender, extender_rc);

            } else {
                aligner_core.align_one_direction(is_reverse_complement, seeder, extender);
            }
        }

        for (auto it = labeled_aligner_cores.begin(); it != labeled_aligner_cores.end(); ++it) {
            AlignerCore &aligner_core = it.value();
            aligner_core.flush();
            main_aligner_core.align_aggregate([&](const auto &callback, const auto &) {
                for (DBGAlignment &path : aligner_core.get_paths()) {
                    assert(path.target_column <= ILabeledDBGAligner::kNTarget);
                    if (path.target_column < ILabeledDBGAligner::kNTarget)
                        callback(std::move(path));
                }
            });
        }

        auto &aggregator = main_aligner_core.get_aggregator().data();
        if (aggregator.size() > config_.num_top_labels) {
            std::vector<std::pair<uint64_t, score_t>> scored_labels;
            scored_labels.reserve(aggregator.size());
            for (const auto &[target, path_queue] : aggregator) {
                scored_labels.emplace_back(target, path_queue.maximum().get_score());
            }

            std::sort(scored_labels.begin(), scored_labels.end(), utils::GreaterSecond());
            auto start = scored_labels.begin() + config_.num_top_labels - 1;
            auto it = std::find_if(start + 1, scored_labels.end(), [&](const auto &a) {
                return a.second < start->second;
            });
            for ( ; it != scored_labels.end(); ++it) {
                aggregator.erase(it->first);
            }
        }

        main_aligner_core.flush();
        callback(header, std::move(paths));
        ++i;
    });
}

template <class BaseSeeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::build_seeder(uint64_t target_column,
               std::string_view query,
               bool is_reverse_complement,
               std::vector<node_index>&& base_nodes,
               const sdsl::bit_vector &signature) const -> Seeder {
    assert(base_nodes.size() == signature.size());
    for (size_t i = 0; i < base_nodes.size(); ++i) {
        if (!signature[i])
            base_nodes[i] = DeBruijnGraph::npos;
    }

    return Seeder(target_column, graph_, query, is_reverse_complement,
                  std::move(base_nodes), config_);
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
