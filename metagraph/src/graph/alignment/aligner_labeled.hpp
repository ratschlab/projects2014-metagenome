#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include "dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/vector_map.hpp"
#include "common/hashers/hash.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class Seeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class LabeledDBGAligner : public SeedAndExtendAligner<Seeder, Extender> {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    LabeledDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : anno_graph_(anno_graph), config_(config) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual ~LabeledDBGAligner() {}

    const DeBruijnGraph& get_graph() const override { return anno_graph_.get_graph(); }
    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    Seeder build_seeder() const override { return Seeder(get_graph(), config_); }
    Extender build_extender() const override { return Extender(anno_graph_, config_); }

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths. Each path has the property that there exists
    // at least one label which is shared by all nodes
    void align_aggregate(DBGQueryAlignment &paths,
                         const AlignmentGenerator &alignment_generator) const override;

  private:
    const AnnotatedDBG& anno_graph_;
    const DBGAlignerConfig config_;
};

template <typename NodeType>
class LabeledColumnExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef typename Extender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename Extender<NodeType>::score_t score_t;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config),
            anno_graph_(anno_graph) {
        assert(this->get_config().check_config_scores());
    }

    virtual ~LabeledColumnExtender() {}

    virtual void initialize(const DBGAlignment &path) override;

  protected:
    virtual std::deque<std::pair<NodeType, char>>
    fork_extension(NodeType node,
                   std::function<void(DBGAlignment&&, NodeType)> callback,
                   score_t min_path_score) override;

  private:
    const AnnotatedDBG &anno_graph_;
    std::vector<uint64_t> target_columns_;

    LabeledColumnExtender fork_extender(std::vector<uint64_t>&& new_target_labels) const;
};


template <class Seeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<Seeder, Extender, AlignmentCompare>
::align_aggregate(DBGQueryAlignment &paths,
                  const AlignmentGenerator &alignment_generator) const {
    AlignmentAggregator<node_index, AlignmentCompare> path_queue(
        paths.get_query(), paths.get_query_reverse_complement(), config_
    );

    alignment_generator(
        [&](DBGAlignment&& alignment) { path_queue.add_alignment(std::move(alignment)); },
        [&](const DBGAlignment &seed) { return path_queue.get_min_path_score(seed); }
    );

    path_queue.call_alignments([&](auto&& alignment) {
        assert(alignment.is_valid(get_graph(), &config_));
        paths.emplace_back(std::move(alignment));
    });
}

template <typename NodeType>
inline void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    DefaultColumnExtender<NodeType>::initialize(path);

    VectorMap<AnnotatedDBG::row_index, size_t> row_map;
    row_map.reserve(path.size());
    for (NodeType i : path) {
        row_map[anno_graph_.graph_to_anno_index(i)]++;
    }

    auto rows = const_cast<std::vector<std::pair<AnnotatedDBG::row_index, size_t>>&&>(
        row_map.values_container()
    );

    // pick all labels which this path could have
    auto label_counts = anno_graph_.get_top_labels(
        rows, anno_graph_.get_annotation().num_labels()
    );

    // encode labels and store their codes
    target_columns_.reserve(label_counts.size());
    target_columns_.clear();

    size_t max_count = 0;
    const auto &label_encoder = anno_graph_.get_annotation().get_label_encoder();
    for (const auto &[label, count] : label_counts) {
        max_count = std::max(max_count, count);
        target_columns_.push_back(label_encoder.encode(label));
    }

    assert(std::is_sorted(target_columns_.begin(), target_columns_.end()));

    mtg::common::logger->trace("Seed has {} labels", target_columns_.size());

    if (max_count != path.size()) {
        mtg::common::logger->warn(
            "Seed has no common labels. Maximum label multiplicity is {} / {}",
            max_count, path.size()
        );
    }
}

template <typename NodeType>
inline LabeledColumnExtender<NodeType> LabeledColumnExtender<NodeType>
::fork_extender(std::vector<uint64_t>&& new_target_labels) const {
    auto fork = *this;
    fork.target_columns_ = std::move(new_target_labels);
    return fork;
}

template <typename NodeType>
inline std::deque<std::pair<NodeType, char>> LabeledColumnExtender<NodeType>
::fork_extension(NodeType node,
                 std::function<void(DBGAlignment&&, NodeType)> callback,
                 score_t min_path_score) {
    if (target_columns_.empty())
        return {};

    const auto &mat = anno_graph_.get_annotation().get_matrix();
    auto base_edges = DefaultColumnExtender<NodeType>::fork_extension(
        node, callback, min_path_score
    );

    std::deque<std::pair<NodeType, char>> edges;
    for (const auto &edge : base_edges) {
        const auto &[out_node, c] = edge;
        AnnotatedDBG::row_index out_row = anno_graph_.graph_to_anno_index(out_node);

        if (target_columns_.size() == 1) {
            if (mat.get(out_row, target_columns_[0]))
                edges.push_back(edge);

        } else {
            auto row = mat.get_row(out_row);
            assert(std::is_sorted(row.begin(), row.end()));

            std::vector<AnnotatedDBG::row_index> intersection;
            intersection.reserve(std::min(row.size(), target_columns_.size()));
            std::set_intersection(target_columns_.begin(), target_columns_.end(),
                                  row.begin(), row.end(), std::back_inserter(intersection));

            if (intersection.size() == target_columns_.size()) {
                edges.push_back(edge);

            } else if (!intersection.empty()) {
                // discard the labels which are in the intersection
                std::vector<uint64_t> diff;
                diff.reserve(target_columns_.size() - intersection.size());
                std::set_difference(target_columns_.begin(), target_columns_.end(),
                                    intersection.begin(), intersection.end(),
                                    std::back_inserter(diff));
                std::swap(diff, target_columns_);

                // assign intersection labels to the fork
                auto fork = fork_extender(std::move(intersection));
                fork.update_columns(node, { edge }, min_path_score);
                fork.extend_main([&](DBGAlignment&& extension, NodeType start_node) {
                    if (start_node)
                        callback(std::move(extension), start_node);
                }, min_path_score);
            }
        }
    }

    return edges;
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
