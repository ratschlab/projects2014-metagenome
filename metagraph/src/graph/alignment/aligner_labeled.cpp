#include "aligner_labeled.hpp"

#include <tsl/hopscotch_set.h>

#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef DeBruijnGraph::node_index node_index;
using MIM = annot::matrix::MultiIntMatrix;


bool check_targets(const DeBruijnGraph &graph,
                   const AnnotatedDBG &anno_graph,
                   const Alignment<node_index> &path) {
    std::string query = path.get_sequence();
    if (dynamic_cast<const RCDBG*>(&graph))
        ::reverse_complement(query.begin(), query.end());

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<annot::binmat::BinaryMatrix::Column> ref_targets;
    for (const std::string &label : anno_graph.get_labels(query, 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (annot::binmat::BinaryMatrix::Column target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}

void DynamicLabeledGraph::flush() {
    auto it = added_nodes_.begin();
#ifndef NDEBUG
    auto kt = added_rows_.begin();
#endif
    for (const auto &labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows_)) {
        assert(it != added_nodes_.end());
        auto jt = targets_set_.emplace(labels).first;
        assert(labels == *jt);
        assert(targets_[*it].first == *kt);
        targets_[*it].second = jt - targets_set_.begin();
        ++it;
#ifndef NDEBUG
        ++jt;
#endif
    }
    assert(it == added_nodes_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

std::vector<size_t> DynamicLabeledGraph::get_coords(node_index node) const {
    std::vector<size_t> coordinates;

    if (const auto *multi_int
            = dynamic_cast<const MIM *>(&anno_graph_.get_annotation().get_matrix())) {
        Row row = AnnotatedDBG::graph_to_anno_index(anno_graph_.get_graph().get_base_node(node));
        for (const auto &[j, tuple] : multi_int->get_row_tuples(row)) {
            for (uint64_t coord : tuple) {
                // TODO: make sure the offsets are correct (query max_int in multi_int)
                // TODO: if this takes up a significant amount of time, preallocate
                //       the entire vector beforehand
                coordinates.push_back(j * 1e15 + coord);
            }
        }
        assert(std::is_sorted(coordinates.begin(), coordinates.end()));
    }

    return coordinates;
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::call_outgoing(NodeType node,
                size_t max_prefetch_distance,
                const std::function<void(NodeType, char /* last char */)> &callback) {
    auto cached_labels = labeled_graph_[node];
    if (this->config_.label_every_n && cached_labels) {
        max_prefetch_distance = std::min(max_prefetch_distance, this->config_.label_every_n);
        std::vector<NodeType> nodes { node };
        std::string seq(this->graph_->get_k(), '#');
        for (size_t i = 0; i < max_prefetch_distance; ++i) {
            size_t outdegree = 0;
            this->graph_->call_outgoing_kmers(nodes.back(), [&](NodeType next, char c) {
                if (c != boss::BOSS::kSentinel) {
                    if (!outdegree) {
                        nodes.push_back(next);
                        seq += c;
                    }

                    ++outdegree;
                }
            });

            if (outdegree != 1)
                break;
        }

        labeled_graph_.add_path(nodes, seq);
        labeled_graph_.flush();
        cached_labels = labeled_graph_[node];
    }

    std::function<void(NodeType, char)> call = callback;

    auto coords = labeled_graph_.get_coords(node);
    assert(std::is_sorted(coords.begin(), coords.end()));

    if (coords.size()) {
        // coordinate consistency
        for (size_t &j : coords) {
            ++j;
        }

        call = [&,bc{std::move(coords)}](NodeType next, char c) {
            auto next_coords = labeled_graph_.get_coords(next);
            assert(std::is_sorted(next_coords.begin(), next_coords.end()));
            if (utils::count_intersection(bc.begin(), bc.end(),
                                          next_coords.begin(), next_coords.end())) {
                callback(next, c);
            }
        };
    } else if (cached_labels) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        call = [&](NodeType next, char c) {
            auto next_labels = labeled_graph_[next];
            // If labels at the next node are not cached, always take the edge.
            // In this case, the label consistency will be checked later.
            // If they are cached, the existence of at least one common label is checked.
            if (!next_labels || utils::count_intersection(cached_labels->get().begin(),
                                                          cached_labels->get().end(),
                                                          next_labels->get().begin(),
                                                          next_labels->get().end())) {
                callback(next, c);
            }
        };
    }

    DefaultColumnExtender<NodeType>::call_outgoing(node, max_prefetch_distance, call);
}

void DynamicLabeledGraph::add_path(const std::vector<node_index> &path,
                                   std::string query) {
    assert(anno_graph_.get_graph().get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    if (path.empty())
        return;

    const DeBruijnGraph &graph = anno_graph_.get_graph();
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph.get_base_graph());
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    if (!canonical && graph.get_mode() == DeBruijnGraph::CANONICAL && query.front() == '#')
        query = graph.get_node_sequence(path[0]) + query.substr(graph.get_k());

    auto [base_path, reversed] = graph.get_base_path(path, query);
    for (size_t i = 0; i < base_path.size(); ++i) {
        if (base_path[i] != DeBruijnGraph::npos) {
            if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i])))
                continue; // skip dummy nodes

            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            if (targets_.emplace(path[i], std::make_pair(row, nannot)).second) {
                added_rows_.push_back(row);
                added_nodes_.push_back(path[i]);
            }
        }
    }
}

void DynamicLabeledGraph::add_node(node_index node) {
    add_path({ node }, std::string(anno_graph_.get_graph().get_k(), '#'));
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::update_seed_filter(node_index node,
                                                               size_t query_start,
                                                               const score_t *s_begin,
                                                               const score_t *s_end) {
    if (SeedFilteringExtender<NodeType>::update_seed_filter(node, query_start, s_begin, s_end)) {
        labeled_graph_.add_node(node);
        return true;
    } else {
        return false;
    }
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::skip_backtrack_start(size_t i) {
    target_intersection_.clear();

    if (this->prev_starts.emplace(i).second) {
        // if backtracking hasn't been started from here yet, get its labels

        auto target_find = diff_target_sets_.find(i);
        if (target_find != diff_target_sets_.end()) {
            // extract a subset of the labels if this node was previously traversed
            target_intersection_ = target_find->second;

        } else if (auto labels = labeled_graph_[std::get<3>(this->table[i])]) {
            if (this->seed_->target_columns.size()) {
                // if the seed had labels, intersect with those
                std::set_intersection(labels->get().begin(),
                                      labels->get().end(),
                                      this->seed_->target_columns.begin(),
                                      this->seed_->target_columns.end(),
                                      std::back_inserter(target_intersection_));
            } else {
                // otherwise take the full label set
                target_intersection_ = *labels;
            }
        }

        // we already have the labels for the first node in the path
        last_path_size_ = 1;
    }

    // skip backtracking from this node if no labels could be determined for it
    return target_intersection_.empty();
}

template <typename NodeType>
void set_target_coordinates(const DynamicLabeledGraph &labeled_graph,
                            const MIM &multi_int,
                            Alignment<NodeType> &alignment) {
    using Column = DynamicLabeledGraph::Column;
    const std::vector<NodeType> &path = alignment.get_nodes();
    const Vector<uint64_t> &target_columns = alignment.target_columns;
    auto &target_coordinates = alignment.target_coordinates;
    target_coordinates.resize(target_columns.size());

    typedef tsl::hopscotch_map<int64_t, std::vector<size_t>> RelativeCoordsMap;
    tsl::hopscotch_map<Column, RelativeCoordsMap> row_coordinates;
    for (Column target : target_columns) {
        row_coordinates.emplace(target, RelativeCoordsMap{});
    }

    auto tuples = multi_int.get_row_tuples(labeled_graph.get_anno_rows(path));
    for (size_t i = 0; i < tuples.size(); ++i) {
        for (const auto &[j, coords] : tuples[i]) {
            auto find = row_coordinates.find(j);
            if (find != row_coordinates.end()) {
                for (int64_t coord : coords) {
                    find.value()[coord - i].emplace_back(i);
                }
            }
        }
    }

    for (auto it = row_coordinates.begin(); it != row_coordinates.end(); ++it) {
        std::vector<std::pair<size_t, std::pair<uint64_t, uint64_t>>> &cur_target_coords
            = target_coordinates[it->first];
        for (auto jt = it.value().begin(); jt != it.value().end(); ++jt) {
            int64_t start = jt->first;
            auto &relative_coords = jt.value();
            if (relative_coords.size()) {
                std::sort(relative_coords.begin(), relative_coords.end());
                relative_coords.erase(std::unique(relative_coords.begin(),
                                                  relative_coords.end()),
                                      relative_coords.end());
                std::pair<size_t, size_t> cur_range(relative_coords[0], relative_coords[0]);
                for (size_t i = 1; i < relative_coords.size(); ++i) {
                    if (relative_coords[i] == cur_range.second + 1) {
                        ++cur_range.second;
                    } else {
                        cur_target_coords.emplace_back(
                            cur_range.first,
                            std::make_pair(cur_range.first + start,
                                           cur_range.second + start)
                        );
                        cur_range.first = relative_coords[i];
                        cur_range.second = relative_coords[i];
                    }
                }
                cur_target_coords.emplace_back(
                    cur_range.first,
                    std::make_pair(cur_range.first + start, cur_range.second + start)
                );
            }
        }

        std::sort(cur_target_coords.begin(), cur_target_coords.end(),
                  [](const auto &a, const auto &b) {
            return std::make_pair(b.second.second - b.second.first, a.first)
                < std::make_pair(a.second.second - a.second.first, b.first);
        });

        auto first_sub_alignment = std::find_if(cur_target_coords.begin(),
                                                cur_target_coords.end(),
                                                [&path](const auto &a) {
            return a.second.second - a.second.first + 1 < path.size();
        });
        cur_target_coords.erase(first_sub_alignment, cur_target_coords.end());
    }
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::call_alignments(score_t cur_cell_score,
                  score_t end_score,
                  score_t min_path_score,
                  const std::vector<node_index> &path,
                  const std::vector<size_t> &trace,
                  const Cigar &ops,
                  size_t clipping,
                  size_t offset,
                  std::string_view window,
                  const std::string &match,
                  const std::function<void(DBGAlignment&&)> & /* callback */) {
    assert(path.size());
    assert(ops.size());
    assert(target_intersection_.size());
    assert(trace.size() >= this->graph_->get_k());

    size_t label_path_end = trace.size() - this->graph_->get_k() + 1;
    if (label_path_end > last_path_size_) {
        for (size_t i = last_path_size_; i < label_path_end; ++i) {
            assert(static_cast<size_t>(i) < path.size());
            if (auto labels = labeled_graph_[path[i]]) {
                Vector<Column> inter;
                std::set_intersection(target_intersection_.begin(),
                                      target_intersection_.end(),
                                      labels->get().begin(), labels->get().end(),
                                      std::back_inserter(inter));

                if (labels->get().size() > inter.size() && this->prev_starts.count(trace[i])) {
                    Vector<Column> diff;
                    auto prev_labels = diff_target_sets_.find(trace[i]);
                    if (prev_labels == diff_target_sets_.end()) {
                        std::set_difference(labels->get().begin(), labels->get().end(),
                                            target_intersection_.begin(),
                                            target_intersection_.end(),
                                            std::back_inserter(diff));
                        if (diff.size()) {
                            diff_target_sets_.emplace(trace[i], std::move(diff));
                            this->prev_starts.erase(trace[i]);
                        }
                    } else {
                        std::set_difference(prev_labels->second.begin(),
                                            prev_labels->second.end(),
                                            target_intersection_.begin(),
                                            target_intersection_.end(),
                                            std::back_inserter(diff));
                        std::swap(prev_labels.value(), diff);
                        if (diff.size())
                            this->prev_starts.erase(trace[i]);
                    }
                }

                std::swap(target_intersection_, inter);

                if (target_intersection_.empty())
                    return;
            }
        }

        last_path_size_ = label_path_end;
    }

    if (ops.back().first == Cigar::MATCH
            && window.size() >= this->config_.min_seed_length
            && end_score - cur_cell_score >= min_path_score) {
        score_t target_score = std::max(
            aggregator_.get_min_path_score(target_intersection_),
            extensions_.get_min_path_score(target_intersection_)
        );

        if (end_score - cur_cell_score >= target_score) {
            DBGAlignment alignment = this->construct_alignment(
                ops, clipping, window, path, match, end_score - cur_cell_score, offset
            );

            assert(!alignment.get_offset());
            alignment.target_columns = target_intersection_;
            assert(check_targets(*this->graph_, labeled_graph_.get_anno_graph(), alignment));

            extensions_.add_alignment(std::move(alignment));
        }
    }
}

template <typename NodeType>
auto LabeledBacktrackingExtender<NodeType>
::extend(score_t min_path_score, bool fixed_seed) -> std::vector<DBGAlignment> {
    extensions_.clear();
    DefaultColumnExtender<NodeType>::extend(min_path_score, fixed_seed);

    const auto *multi_int = dynamic_cast<const MIM *>(
        &labeled_graph_.get_anno_graph().get_annotation().get_matrix()
    );

    std::vector<DBGAlignment> extensions;
    extensions_.call_alignments(
        [&](DBGAlignment&& alignment) {
            if (multi_int)
                set_target_coordinates(labeled_graph_, *multi_int, alignment);

            extensions.emplace_back(std::move(alignment));
        },
        [&]() { return extensions.size() && extensions.back().get_score() < min_path_score; }
    );

    return extensions;
}

template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
