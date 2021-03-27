#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"

namespace mtg {
namespace graph {

class AnnotatedDBG;

namespace align {


class ILabeledDBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef Vector<uint64_t> Targets;

    ILabeledDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config);

    virtual ~ILabeledDBGAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    typedef std::pair<std::vector<node_index> /* forward */,
                      std::vector<node_index> /* reverse complement */ > Mapping;
    typedef std::pair<sdsl::bit_vector /* forward */,
                      sdsl::bit_vector /* reverse complement */ > Signature;
    typedef std::vector<Mapping> BatchMapping;
    typedef std::vector<std::pair<Targets, Signature>> QueryLabels;
    typedef std::vector<QueryLabels> BatchLabels;

    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    std::pair<BatchMapping, BatchLabels>
    map_and_label_query_batch(const QueryGenerator &generate_query) const;
};

bool check_targets(const AnnotatedDBG &anno_graph,
                   const Alignment<DeBruijnGraph::node_index> &path);

class LabeledSeedFilter : public ISeedFilter {
  public:
    LabeledSeedFilter(size_t k) : k_(k) {}
    Vector<uint64_t> labels_to_keep(const DBGAlignment &seed);
    void update_seed_filter(const LabeledNodeRangeGenerator &generator);

  private:
    size_t k_;
    tsl::hopscotch_map<node_index, tsl::hopscotch_map<uint64_t, std::pair<size_t, size_t>>> visited_nodes_;
};

template <class BaseSeeder>
class LabeledSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::node_index node_index;
    typedef typename BaseSeeder::Seed Seed;
    typedef ILabeledDBGAligner::Targets Targets;

    template <typename... Args>
    LabeledSeeder(Targets&& target_columns, Args&&... args)
          : BaseSeeder(std::forward<Args>(args)...),
            target_columns_(std::move(target_columns)) {}

    virtual ~LabeledSeeder() {}

    virtual std::vector<Seed> get_seeds() const override {
        std::vector<Seed> seeds = BaseSeeder::get_seeds();
        if (target_columns_.size()) {
            for (Seed &seed : seeds) {
                if (!seed.get_offset())
                    seed.target_columns = target_columns_;
            }
        }

        return seeds;
    }

  private:
    Targets target_columns_;
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

    Seeder build_seeder(Vector<uint64_t>&& target_columns,
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
    typedef ILabeledDBGAligner::Targets Targets;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~LabeledColumnExtender() {}

    virtual void initialize(const DBGAlignment &path) override;

  protected:
    typedef std::vector<std::pair<NodeType, char>> Edges;
    typedef std::vector<std::tuple<NodeType, char, size_t>> LabeledEdges;
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

        if (!align_node_to_target_.count(node))
            align_node_to_target_[node] = align_node_to_target_[prev];

        DefaultColumnExtender<NodeType>::add_scores_to_column(
            column, std::move(scores), node
        );
    }

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &,
                                      const AlignNode &node) const override {
        assert(align_node_to_target_.count(node));
        size_t target_columns_idx = align_node_to_target_[node];

        auto find = backtrack_start_counter_.find(target_columns_idx);
        assert(find == backtrack_start_counter_.end()
            || find->second <= this->config_.num_alternative_paths);

        return find != backtrack_start_counter_.end()
            && find->second == this->config_.num_alternative_paths;
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
            if (extensions.back().get_offset()) {
                extensions.pop_back();
                return;
            }

            size_t target_columns_idx = align_node_to_target_.find(best_node)->second;
            assert(target_columns_idx < target_columns_.size());

            ++backtrack_start_counter_[target_columns_idx];
            extensions.back().target_columns
                = *(target_columns_.begin() + target_columns_idx);

            assert(check_targets(anno_graph_, extensions.back()));
        }
    }

  private:
    const AnnotatedDBG &anno_graph_;
    mutable VectorSet<Targets, utils::VectorHash> target_columns_;
    mutable tsl::hopscotch_map<NodeType, tsl::hopscotch_map<size_t, LabeledEdges>> cached_edge_sets_;
    mutable tsl::hopscotch_map<AlignNode, size_t, AlignNodeHash> align_node_to_target_;
    mutable tsl::hopscotch_map<size_t, size_t> backtrack_start_counter_;

    AlignNode get_next_align_node(NodeType node, char c, size_t dist_from_origin) const;
};


template <class BaseSeeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    typedef SeedAndExtendAlignerCore<AlignmentCompare> AlignerCore;
    auto mapped_batch = map_and_label_query_batch(generate_query);

    size_t i = 0;
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        const auto &[query_nodes_pair, target_columns] = mapped_batch;
        assert(config_.num_alternative_paths);

        LabeledSeedFilter seed_filter(this->graph_.get_k());
        AlignerCore aligner_core(graph_, config_, seed_filter, query, is_reverse_complement);
        DBGQueryAlignment &paths = aligner_core.get_paths();

        std::string_view this_query = paths.get_query(is_reverse_complement);
        std::string_view reverse = paths.get_query(true);
        assert(this_query == query);

        Extender extender(anno_graph_, config_, this_query);
        Extender extender_rc(anno_graph_, config_, reverse);

        const auto &[nodes, nodes_rc] = query_nodes_pair[i];

        for (const auto &[target_columns, signature_pair] : target_columns[i]) {
            const auto &[signature, signature_rc] = signature_pair;

            Seeder seeder = build_seeder(
                Targets(target_columns),
                this_query, // use this_query since paths stores a copy
                is_reverse_complement,
                std::vector<node_index>(nodes),
                signature
            );

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
                assert(!is_reverse_complement);

                auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                         const auto &callback) {
                    callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
                };

                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                Seeder seeder_rc = build_seeder(Targets(target_columns),
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_both_directions(seeder, seeder_rc,
                                                   extender, extender_rc,
                                                   build_rev_comp_alignment_core);

            } else if (config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);

                Seeder seeder_rc = build_seeder(Targets(target_columns),
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_best_direction(seeder, seeder_rc, extender, extender_rc);

            } else {
                aligner_core.align_one_direction(is_reverse_complement, seeder, extender);
            }
        }

        auto &aggregator = aligner_core.get_aggregator().data();
        aggregator.erase(std::numeric_limits<uint64_t>::max());
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

        aligner_core.flush();
        callback(header, std::move(paths));
        ++i;
    });
}

template <class BaseSeeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::build_seeder(Vector<uint64_t>&& target_columns,
               std::string_view query,
               bool is_reverse_complement,
               std::vector<node_index>&& base_nodes,
               const sdsl::bit_vector &signature) const -> Seeder {
    assert(base_nodes.size() == signature.size());
    for (size_t i = 0; i < base_nodes.size(); ++i) {
        if (!signature[i])
            base_nodes[i] = DeBruijnGraph::npos;
    }

    return Seeder(std::move(target_columns), graph_, query, is_reverse_complement,
                  std::move(base_nodes), config_);
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
