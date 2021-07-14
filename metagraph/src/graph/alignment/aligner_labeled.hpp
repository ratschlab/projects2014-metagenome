#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_set.hpp"
#include "common/hashers/hash.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

namespace mtg {
namespace graph {
namespace align {


class DynamicLabeledGraph {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;
    typedef VectorSet<Vector<Column>, utils::VectorHash> Storage;
    typedef Storage::const_iterator const_iterator;

    DynamicLabeledGraph(const AnnotatedDBG &anno_graph) : anno_graph_(anno_graph) {
        targets_set_.emplace(); // insert empty vector
    }

    const DeBruijnGraph& get_graph() const { return anno_graph_.get_graph(); }
    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }

    void flush();
    void add_node(node_index node);
    void add_path(const std::vector<node_index> &path, std::string sequence);

    // get all sequences coordinates (k-mers) generating the node
    // coordinates from different labels are shifted by offsets
    std::vector<size_t> get_coords(node_index node) const;

    const_iterator begin() const { return targets_set_.begin(); }
    const_iterator end() const { return targets_set_.end(); }
    const_iterator find(node_index node) const {
        auto it = targets_.find(node);
        if (it == targets_.end() || it->second == std::numeric_limits<size_t>::max()) {
            return targets_set_.end();
        } else {
            return targets_set_.begin() + it->second;
        }
    }

  private:
    const AnnotatedDBG &anno_graph_;

    Storage targets_set_;
    // map nodes to indexes in targets_set_
    tsl::hopscotch_map<node_index, size_t> targets_;

    std::vector<Row> added_rows_;
    std::vector<node_index> added_nodes_;
};


template <class AlignmentCompare = LocalAlignmentLess>
class ILabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    ILabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(anno_graph.get_graph(), config),
            labeled_graph_(anno_graph) {}

    virtual ~ILabeledAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const IDBGAligner::AlignmentCallback &callback) const override {
        ISeedAndExtendAligner<AlignmentCompare>::align_batch(seq_batch,
            [&](std::string_view header, IDBGAligner::DBGQueryAlignment&& alignments) {
                auto it = std::remove_if(
                    alignments.begin(), alignments.end(),
                    [](const auto &a) { return a.target_columns.empty(); }
                );
                alignments.erase(it, alignments.end());

                callback(header, std::move(alignments));
            }
        );
    }

  protected:
    mutable DynamicLabeledGraph labeled_graph_;
};


template <typename NodeType = DeBruijnGraph::node_index>
class LabeledBacktrackingExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef DefaultColumnExtender<DeBruijnGraph::node_index> BaseExtender;
    typedef typename BaseExtender::score_t score_t;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;
    typedef typename BaseExtender::ScoreVec ScoreVec;
    typedef AlignmentAggregator<node_index, LocalAlignmentLess> Aggregator;

    LabeledBacktrackingExtender(DynamicLabeledGraph &anno_graph,
                                const DBGAlignerConfig &config,
                                const Aggregator &aggregator,
                                std::string_view query)
          : BaseExtender(anno_graph.get_graph(), config, query),
            labeled_graph_(anno_graph),
            aggregator_(aggregator),
            extensions_(anno_graph.get_graph(), aggregator_.get_query(false),
                        aggregator_.get_query(true), this->config_) {}

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual std::vector<DBGAlignment> extend(score_t min_path_score) override;

    virtual void init_backtrack() override {
        labeled_graph_.flush();
        diff_target_sets_.clear();
    }

    virtual bool terminate_backtrack_start(const std::vector<DBGAlignment> &) const override { return false; }

    virtual bool terminate_backtrack() const override { return target_intersection_.empty(); }

    virtual bool skip_backtrack_start(size_t i) override;

    virtual bool update_seed_filter(node_index node,
                                    size_t query_start,
                                    const score_t *s_begin,
                                    const score_t *s_end) override;

    virtual void call_outgoing(NodeType node,
                               size_t max_prefetch_distance,
                               const std::function<void(NodeType, char)> &callback) override;

    virtual void call_alignments(score_t cur_cell_score,
                                 score_t end_score,
                                 score_t min_path_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> &trace,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 const std::function<void(DBGAlignment&&)> &callback) override;

  private:
    DynamicLabeledGraph &labeled_graph_;
    const Aggregator &aggregator_;
    Aggregator extensions_;
    Vector<uint64_t> target_intersection_;
    size_t last_path_size_;
    tsl::hopscotch_map<size_t, Vector<uint64_t>> diff_target_sets_;
};


template <class Extender = LabeledBacktrackingExtender<>,
          class Seeder = UniMEMSeeder<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ILabeledAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    LabeledAligner(Args&&... args)
          : ILabeledAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  protected:
    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query,
                   const typename Extender::Aggregator &aggregator) const override {
        return std::make_shared<Extender>(this->labeled_graph_, this->get_config(), aggregator, query);
    }

    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<DeBruijnGraph::node_index> &nodes) const override {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
