#include "aligner_labeled.hpp"

#include <tsl/hopscotch_set.h>

#include "aligner_prefix_suffix.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/utils/template_utils.hpp"

namespace mtg {
namespace graph {
namespace align {


bool check_targets(const AnnotatedDBG &anno_graph,
                   const Alignment<DeBruijnGraph::node_index> &path) {
    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<uint64_t> ref_targets;
    for (const std::string &label : anno_graph.get_labels(path.get_sequence(), 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (uint64_t target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}


void
process_seq_path(const DeBruijnGraph &graph,
                 std::string_view query,
                 const std::vector<DeBruijnGraph::node_index> &query_nodes,
                 const std::function<void(AnnotatedDBG::row_index, size_t)> &callback) {
    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&graph);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!dbg_succ && canonical)
        dbg_succ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());

    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    auto run_callback = [&](DeBruijnGraph::node_index node, size_t i) {
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
            graph.map_to_nodes(map_query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    run_callback(node, i);

                ++i;
            });
        } else {
            graph.map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    run_callback(node, i);

                ++i;
            });
        }
        assert(i == query_nodes.size());
    }
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>::initialize(const DBGAlignment &seed) {
    DefaultColumnExtender<NodeType>::initialize(seed);
    min_scores_.clear();
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::populate_min_scores(const Vector<uint64_t> &targets) {
    for (uint64_t target : targets) {
        score_t score = aggregator_.get_min_path_score(target);
        if (!min_scores_.count(target)) {
            min_scores_.emplace(target, score);
        } else {
            min_scores_[target] = std::max(min_scores_[target], score);
        }
    }
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::update_seed_filter(size_t j) {
    bool result = DefaultColumnExtender<NodeType>::update_seed_filter(j);

    const auto &[S, E, F, OS, OE, OF, next, i_prev, c, offset, max_pos, begin] = this->table[j];
    auto [it, inserted] = targets_.emplace(next, 0);
    if (inserted) {
        process_seq_path(this->graph_, std::string(this->graph_.get_k(), '#'),
                         { next }, [&](auto row, size_t) {
            added_rows.push_back(row);
            added_nodes.push_back(next);
        });
    } else {
        populate_min_scores(*(targets_set_.begin() + it->second));
    }

    return result;
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>::init_backtrack() {
    auto it = added_nodes.begin();
    for (const auto &labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows)) {
        assert(it != added_nodes.end());
        populate_min_scores(labels);
        auto jt = targets_set_.emplace(labels).first;
        targets_[*it] = jt - targets_set_.begin();
        ++it;
    }
    assert(it == added_nodes.end());

    added_rows.clear();
    added_nodes.clear();
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::process_extension(DBGAlignment&& alignment,
                    const std::vector<size_t> &trace,
                    tsl::hopscotch_set<size_t> &prev_starts,
                    const std::function<void(DBGAlignment&&)> &callback) {
    assert(alignment.is_valid(this->graph_, &this->config_));
    assert(trace.size() == alignment.size());

    if (alignment.empty() || alignment.get_offset()) {
        prev_starts.insert(trace.begin(), trace.end());
        return;
    }

    size_t target_id = targets_[alignment.back()];
    if (!target_id)
        return;

    Vector<uint64_t> target_intersection = *(targets_set_.begin() + target_id);

    AlignmentSuffix<node_index> suffix(alignment, this->config_, this->graph_);
    while (!suffix.get_added_offset())
        ++suffix;

    auto suffix_shift = [&suffix]() {
        assert(!suffix.reof());
        --suffix;
        while (!suffix.reof() && suffix.get_front_op() == Cigar::INSERTION) {
            --suffix;
        }
    };

    suffix_shift();

    auto it = trace.begin();

    bool skipped = false;
    prev_starts.emplace(*it);
    ++it;

    auto is_valid = [&]() -> bool {
        return suffix.get_front_op() == Cigar::MATCH
            && suffix.get_score() >= this->config_.min_cell_score
            && target_intersection.size();
    };

    AlignmentSuffix<node_index> last_valid(suffix);
    Vector<uint64_t> last_valid_target_intersection;

    auto aln_from_suffix = [&]() -> bool {
        DBGAlignment aln_suffix(last_valid);
        aln_suffix.target_columns = last_valid_target_intersection;
        assert(check_targets(anno_graph_, aln_suffix));
        auto target_it = std::remove_if(
            aln_suffix.target_columns.begin(),
            aln_suffix.target_columns.end(),
            [&](uint64_t target) {
                if (aln_suffix.get_score() > min_scores_[target]) {
                    min_scores_[target] = aln_suffix.get_score();
                    return false;
                }

                return true;
            }
        );
        aln_suffix.target_columns.erase(target_it, aln_suffix.target_columns.end());

        if (aln_suffix.target_columns.size()) {
            callback(std::move(aln_suffix));
            return true;
        }

        return false;
    };

    for (size_t i = alignment.size() - 1; i > 0; --i) {
        const auto &cur_targets = *(targets_set_.begin() + targets_[alignment[i - 1]]);
        Vector<uint64_t> inter;
        std::set_intersection(target_intersection.begin(), target_intersection.end(),
                              cur_targets.begin(), cur_targets.end(),
                              std::back_inserter(inter));

        if (inter.empty())
            break;

        assert(it != trace.end());
        if ((!skipped && inter.size() == target_intersection.size())
                || inter.size() == cur_targets.size()) {
            prev_starts.emplace(*it);
        } else if (!skipped) {
            Vector<uint64_t> diff;
            std::set_difference(cur_targets.begin(), cur_targets.end(),
                                target_intersection.begin(), target_intersection.end(),
                                std::back_inserter(diff));
            if (diff.empty()) {
                prev_starts.emplace(*it);
            } else {
                skipped = true;
            }
        }
        ++it;

        if (is_valid()) {
            last_valid = suffix;
            last_valid_target_intersection = target_intersection;
        }

        if (inter.size() < target_intersection.size())
            aln_from_suffix();

        suffix_shift();
        std::swap(target_intersection, inter);
    }

    if (is_valid()) {
        std::swap(last_valid, suffix);
        std::swap(last_valid_target_intersection, target_intersection);
    }

    aln_from_suffix();
}

template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
