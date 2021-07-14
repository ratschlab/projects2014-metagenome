#include "aligner_labeled.hpp"

#include <tsl/hopscotch_set.h>

#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef DeBruijnGraph::node_index node_index;


bool check_targets(const DeBruijnGraph &graph,
                   const AnnotatedDBG &anno_graph,
                   const Alignment<node_index> &path) {
    std::string query = path.get_sequence();
    if (dynamic_cast<const RCDBG*>(&graph))
        ::reverse_complement(query.begin(), query.end());

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<uint64_t> ref_targets;
    for (const std::string &label : anno_graph.get_labels(query, 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (uint64_t target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}


template <class Callback>
void process_seq_path(const DeBruijnGraph &graph,
                      std::string_view query,
                      const std::vector<node_index> &query_nodes,
                      const Callback &callback) {
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph.get_base_graph());
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    auto run_callback = [&](node_index node, size_t i) {
        if (!boss || boss->get_W(dbg_succ->kmer_to_boss_index(node)))
            callback(AnnotatedDBG::graph_to_anno_index(node), i);
    };

    if (canonical) {
        if (query_nodes.size()) {
            auto first = std::find_if(query_nodes.begin(), query_nodes.end(),
                                      [](auto i) -> bool { return i; });
            if (first == query_nodes.end())
                return;

            size_t start = first - query_nodes.begin();

            if (canonical->get_base_node(*first) == *first) {
                for (size_t i = start; i < query_nodes.size(); ++i) {
                    if (query_nodes[i] != DeBruijnGraph::npos)
                        run_callback(canonical->get_base_node(query_nodes[i]), i);
                }
            } else {
                for (size_t i = query_nodes.size(); i > start; --i) {
                    if (query_nodes[i - 1] != DeBruijnGraph::npos)
                        run_callback(canonical->get_base_node(query_nodes[i - 1]), i - 1);
                }
            }
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos)
                run_callback(query_nodes[i], i);
        }
    } else {
        size_t i = 0;
        if (query.front() == '#') {
            std::string map_query
                = graph.get_node_sequence(query_nodes[0]).substr(0, graph.get_k());
            map_query += query.substr(graph.get_k());
            graph.map_to_nodes(map_query, [&](node_index node) {
                if (node != DeBruijnGraph::npos)
                    run_callback(node, i);

                ++i;
            });
        } else {
            graph.map_to_nodes(query, [&](node_index node) {
                if (node != DeBruijnGraph::npos)
                    run_callback(node, i);

                ++i;
            });
        }
        assert(i == query_nodes.size());
    }
}

void DynamicLabeledGraph::flush() {
    auto it = added_nodes_.begin();
    for (const auto &labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows_)) {
        assert(it != added_nodes_.end());
        auto jt = targets_set_.emplace(labels).first;
        assert(labels == *jt);
        targets_[*it] = jt - targets_set_.begin();
        ++it;
    }
    assert(it == added_nodes_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

std::vector<size_t> DynamicLabeledGraph::get_coords(node_index /* node */) const {
    return {};
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::call_outgoing(NodeType node,
                size_t max_prefetch_distance,
                const std::function<void(NodeType, char)> &callback) {
    auto it = anno_graph_.find(node);
    if (this->config_.label_every_n && it == anno_graph_.end()) {
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

        anno_graph_.add_path(nodes, seq);
        anno_graph_.flush();
        it = anno_graph_.find(node);
    }

    std::function<void(NodeType, char)> call = callback;

    auto coords = anno_graph_.get_coords(node);
    if (coords.size()) {
        for (size_t &j : coords) {
            ++j;
        }

        call = [this,oc=call,bc=coords](NodeType next, char c) {
            auto next_coords = anno_graph_.get_coords(next);
            if (utils::count_intersection(bc.begin(), bc.end(),
                                          next_coords.begin(), next_coords.end())) {
                oc(next, c);
            }
        };
    }

    if (it != anno_graph_.end()) {
        call = [this,it,oc=call](NodeType next, char c) {
            auto next_it = anno_graph_.find(next);
            if (next_it == anno_graph_.end()
                    || utils::count_intersection(it->begin(), it->end(),
                                                 next_it->begin(), next_it->end())) {
                oc(next, c);
            }
        };
    }

    DefaultColumnExtender<NodeType>::call_outgoing(node, max_prefetch_distance, call);
}

void DynamicLabeledGraph::add_path(const std::vector<node_index> &path,
                                   std::string_view spelling) {
    process_seq_path(get_graph(), spelling, path, [&](auto row, size_t i) {
        if (targets_.emplace(path[i], std::numeric_limits<size_t>::max()).second) {
            added_rows_.push_back(row);
            added_nodes_.push_back(path[i]);
        }
    });
}

void DynamicLabeledGraph::add_node(node_index node) {
    std::string dummy(get_graph().get_k(), '#');
    add_path(std::vector<node_index>{ node }, dummy);
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::update_seed_filter(node_index node,
                                                               size_t query_start,
                                                               const score_t *s_begin,
                                                               const score_t *s_end) {
    if (SeedFilteringExtender<NodeType>::update_seed_filter(node, query_start, s_begin, s_end)) {
        anno_graph_.add_node(node);
        return true;
    } else {
        return false;
    }
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::skip_backtrack_start(size_t i) {
    target_intersection_.clear();
    if (this->prev_starts.emplace(i).second) {
        auto target_find = diff_target_sets_.find(i);
        if (target_find != diff_target_sets_.end()) {
            target_intersection_ = target_find->second;

        } else {
            NodeType node = std::get<3>(this->table[i]);
            auto find = anno_graph_.find(node);
            if (find == anno_graph_.end() || find->empty())
                return true;

            target_intersection_ = *find;
            if (this->seed_->target_columns.size()) {
                Vector<uint64_t> inter;
                std::set_intersection(target_intersection_.begin(),
                                      target_intersection_.end(),
                                      this->seed_->target_columns.begin(),
                                      this->seed_->target_columns.end(),
                                      std::back_inserter(inter));
                if (inter.empty())
                    return true;

                std::swap(inter, target_intersection_);
            }
        }

        last_path_size_ = 1;

        return false;
    } else {
        return true;
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
            auto find = anno_graph_.find(path[i]);
            if (find != anno_graph_.end()) {
                Vector<uint64_t> inter;
                std::set_intersection(target_intersection_.begin(),
                                      target_intersection_.end(),
                                      find->begin(), find->end(),
                                      std::back_inserter(inter));

                if (find->size() > inter.size() && this->prev_starts.count(trace[i])) {
                    Vector<uint64_t> diff;
                    auto prev_find = diff_target_sets_.find(trace[i]);
                    if (prev_find == diff_target_sets_.end()) {
                        std::set_difference(find->begin(), find->end(),
                                            target_intersection_.begin(),
                                            target_intersection_.end(),
                                            std::back_inserter(diff));
                        if (diff.size()) {
                            diff_target_sets_.emplace(trace[i], std::move(diff));
                            this->prev_starts.erase(trace[i]);
                        }
                    } else {
                        std::set_difference(prev_find->second.begin(),
                                            prev_find->second.end(),
                                            target_intersection_.begin(),
                                            target_intersection_.end(),
                                            std::back_inserter(diff));
                        std::swap(prev_find.value(), diff);
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
            assert(check_targets(*this->graph_, anno_graph_.get_anno_graph(), alignment));

            extensions_.add_alignment(std::move(alignment));
        }
    }
}

template <typename NodeType>
auto LabeledBacktrackingExtender<NodeType>
::extend(score_t min_path_score) -> std::vector<DBGAlignment> {
    extensions_.clear();
    std::vector<DBGAlignment> extensions;

    DefaultColumnExtender<NodeType>::extend(min_path_score);
    extensions_.call_alignments(
        [&](DBGAlignment&& alignment) { extensions.emplace_back(std::move(alignment)); },
        [&]() { return extensions.size() && extensions.back().get_score() < min_path_score; }
    );

    return extensions;
}

template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
