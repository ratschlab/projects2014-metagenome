#include "aligner_extender_methods.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef DBGAlignerConfig::score_t score_t;
constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;


template <typename NodeType>
DefaultColumnExtender<NodeType>::DefaultColumnExtender(const DeBruijnGraph &graph,
                                                       const DBGAlignerConfig &config,
                                                       std::string_view query)
      : graph_(graph), config_(config), query_(query) {
    assert(config_.check_config_scores());
    partial_sums_.reserve(query_.size() + 1);
    partial_sums_.resize(query_.size(), 0);
    std::transform(query_.begin(), query_.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query_) == partial_sums_.front());
    assert(config_.get_row(query_.back())[query_.back()] == partial_sums_.back());
    partial_sums_.push_back(0);

    for (char c : graph_.alphabet()) {
        auto &p_score_row = profile_score_.emplace(c, query_.size() + 9).first.value();
        auto &p_op_row = profile_op_.emplace(c, query_.size() + 9).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = kCharToOp[c];

        // the first cell in a DP table row is one position before the last matched
        // character, so we need to shift the indices of profile_score_ and profile_op_
        std::transform(query_.begin(), query_.end(), p_score_row.begin() + 1,
                       [&row](char q) { return row[q]; });

        std::transform(query_.begin(), query_.end(), p_op_row.begin() + 1,
                       [&op_row](char q) { return op_row[q]; });
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &seed) {
    assert(seed.size());
    assert(seed.get_cigar().size());
    assert(seed.get_cigar().back().first == Cigar::MATCH
        || seed.get_cigar().back().first == Cigar::MISMATCH);
    assert(seed.is_valid(graph_, &config_));

    auto it = conv_checker.find(seed.back());

    if (it != conv_checker.end()
            && it->second.at(seed.get_query().size() + seed.get_clipping()) < seed.get_score())
        it = conv_checker.end();

    if (it == conv_checker.end()) {
        seed_ = &seed;
        table.clear();
    } else {
        DEBUG_LOG("Skipping seed: {}", seed);
        seed_ = nullptr;
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::process_extension(DBGAlignment&& extension,
                    const std::vector<size_t> &trace,
                    tsl::hopscotch_set<size_t> &prev_starts,
                    const std::function<void(DBGAlignment&&)> &callback) {
    prev_starts.insert(trace.begin(), trace.end());
    callback(std::move(extension));
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>::update_seed_filter(size_t i) {
    const auto &[S, E, F, OS, OE, OF, next, i_prev, c, offset, max_pos, begin] = table[i];
    auto it = conv_checker.find(next);
    if (it == conv_checker.end()) {
        ScoreVec &s_merged = conv_checker.emplace(next, ScoreVec(query_.size() + 1, ninf)).first.value();
        std::copy(S.begin(), S.end(), s_merged.begin() + start + begin);
        return true;
    } else {
        ScoreVec &s_merged = it.value();
        bool converged = true;
        for (size_t j = begin; j < S.size() + begin; ++j) {
            if (S[j - begin] > s_merged[start + j]) {
                converged = false;
                s_merged[j + start] = S[j - begin];
            }
        }

        return !converged;
    }
}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>::get_extensions(score_t min_path_score)
        -> std::vector<DBGAlignment> {
    std::vector<DBGAlignment> extensions;

    if (!seed_)
        return extensions;


    std::string_view window(seed_->get_query().data(),
                            query_.data() + query_.size() - seed_->get_query().data());

    start = seed_->get_clipping();

    table.emplace_back(ScoreVec(window.size() + 1, 0),
                       ScoreVec(window.size() + 1, ninf),
                       ScoreVec(window.size() + 1, ninf),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       seed_->front(), -1, '\0', seed_->get_offset() - 1, 0, 0);

    score_t xdrop_cutoff = std::max(-config_.xdrop, ninf + 1);
    assert(config_.xdrop > 0);
    assert(xdrop_cutoff < 0);

    using Ref = std::tuple<score_t,
                           size_t /* table idx */,
                           ssize_t /* max pos */,
                           ssize_t /* offset */>;
    Ref best_score { 0, 0, 0, seed_->get_offset() - 1 };

    struct RefCmp {
        bool operator()(const Ref &a, const Ref &b) const {
            const auto &[a_score, a_id, a_max_pos, a_offset] = a;
            const auto &[b_score, b_id, b_max_pos, b_offset] = b;
            return std::make_pair(a_score, std::abs(b_max_pos - b_offset))
                < std::make_pair(b_score, std::abs(a_max_pos - a_offset));
        }
    };

    std::priority_queue<Ref, std::vector<Ref>, RefCmp> queue;
    queue.emplace(best_score);

    assert(partial_sums_.at(start) == config_.match_score(window));

    while (queue.size()) {
        size_t i = std::get<1>(queue.top());
        queue.pop();

        std::vector<std::pair<NodeType, char>> outgoing;
        size_t next_offset = -1;

        size_t begin = 0;
        size_t end = window.size() + 1;

        {
            const auto &[S, E, F, OS, OE, OF, node, i_prev, c, offset, max_pos, trim] = table[i];
            next_offset = offset + 1;

            if (static_cast<double>(table.size()) / window.size() >= config_.max_nodes_per_seq_char)
                continue;

            auto in_range = [xdrop_cutoff](score_t s) { return s >= xdrop_cutoff; };

            begin = std::find_if(S.begin(), S.end(), in_range) - S.begin() + trim;
            end = std::find_if(S.rbegin(), S.rend(), in_range).base() - S.begin() + trim;

            if (end != window.size() + 1)
                ++end;

            if (end <= begin)
                continue;

            bool has_extension = false;
            size_t loop_end = std::min(end, S.size() + trim);
            for (size_t j = begin; j < loop_end; ++j) {
                assert(partial_sums_.at(start + j) == config_.match_score(window.substr(j)));
                score_t ext_score = S[j - trim] + partial_sums_.at(start + j);
                if ((config_.num_alternative_paths == 1 && ext_score > std::get<0>(best_score))
                        || ext_score >= min_path_score) {
                    has_extension = true;
                    break;
                }
            }

            if (!has_extension)
                continue;

            if (next_offset - seed_->get_offset() < seed_->get_sequence().size()) {
                if (next_offset < graph_.get_k()) {
                    outgoing.emplace_back(seed_->front(), seed_->get_sequence()[next_offset - seed_->get_offset()]);
                } else {
                    outgoing.emplace_back((*seed_)[next_offset - graph_.get_k() + 1],
                                          seed_->get_sequence()[next_offset - seed_->get_offset()]);
                    assert(graph_.traverse(node, outgoing.back().second) == outgoing.back().first);
                }
            } else {
                graph_.call_outgoing_kmers(node, [&](NodeType next, char c) {
                    if (c != boss::BOSS::kSentinel)
                        outgoing.emplace_back(next, c);
                });
            }
        }

        for (const auto &[next, c] : outgoing) {
            assert(std::get<0>(best_score) > xdrop_cutoff);

            table.emplace_back(ScoreVec(end - begin, ninf),
                               ScoreVec(end - begin, ninf),
                               ScoreVec(end - begin, ninf),
                               OpVec(end - begin, Cigar::CLIPPED),
                               OpVec(end - begin, Cigar::CLIPPED),
                               OpVec(end - begin, Cigar::CLIPPED),
                               next, i, c, next_offset, begin, begin);

            const auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev,
                         node_prev, i_prev, c_prev, offset_prev, max_pos_prev, trim_prev] = table[i];

            auto &[S, E, F, OS, OE, OF, node_cur, i_cur, c_stored, offset, max_pos, trim] = table.back();

            if (next != node_prev && !trim && !trim_prev)
                S[0] = S_prev[0];

            assert(i_cur == i);
            assert(node_cur == next);
            assert(c_stored == c);
            assert(offset == offset_prev + 1);
            assert(c == graph_.get_node_sequence(node_cur)[std::min(graph_.get_k() - 1, offset)]);

            const int8_t *profile_scores = profile_score_[c].data() + start;
            const Cigar::Operator *profile_ops = profile_op_[c].data() + start;

            // match scores
            size_t loop_begin = std::max({ begin, trim_prev });
            size_t loop_end = std::min({ S_prev.size() + trim_prev, end - 1 });

            auto update_match = [S_prev_v=S_prev.data() - trim_prev,
                                 S_v=S.data() - begin + 1,
                                 OS_v=OS.data() - begin + 1,
                                 profile_scores,profile_ops,xdrop_cutoff](size_t j) {
                score_t match = S_prev_v[j] + profile_scores[j + 1];
                if (match > std::max(S_v[j], xdrop_cutoff - 1)) {
                    S_v[j] = match;
                    OS_v[j] = profile_ops[j + 1];
                }
            };

            auto update_del = [S_prev_v=S_prev.data() - trim_prev,
                               F_prev_v=F_prev.data() - trim_prev,
                               F_v=F.data() - begin,
                               S_v=S.data() - begin,
                               OS_v=OS.data() - begin,
                               OF_v=OF.data() - begin,
                               this,xdrop_cutoff](size_t j) {
                score_t del_open = S_prev_v[j] + config_.gap_opening_penalty;
                score_t del_extend = F_prev_v[j] + config_.gap_extension_penalty;
                score_t del_score = std::max(del_open, del_extend);

                if (del_score >= xdrop_cutoff) {
                    F_v[j] = del_score;
                    OF_v[j] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
                    if (del_score > S_v[j]) {
                        S_v[j] = del_score;
                        OS_v[j] = Cigar::DELETION;
                    }
                }
            };

            if (loop_begin > begin)
                update_match(begin);

            #pragma omp simd
            for (size_t j = loop_begin; j < loop_end; ++j) {
                update_match(j);
                update_del(j);
            }

            if (end - 1 > loop_end)
                update_del(end - 1);

            // update insertion scores
            // since each element is dependent on the previous one, this can't
            // be vectorized easily
            #pragma omp simd
            for (size_t j = begin + 1; j < end; ++j) {
                score_t ins_open = S[j - 1 - begin] + config_.gap_opening_penalty;
                score_t ins_extend = E[j - 1 - begin] + config_.gap_extension_penalty;
                score_t ins_score = std::max(ins_open, ins_extend);

                if (ins_score >= xdrop_cutoff) {
                    E[j - begin] = ins_score;
                    OE[j - begin] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

                    if (ins_score > S[j - begin]) {
                        S[j - begin] = ins_score;
                        OS[j - begin] = Cigar::INSERTION;
                    }
                }
            }

            while (S.back() >= xdrop_cutoff && S.size() + trim < window.size() + 1) {
                score_t ins_open = S.back() + config_.gap_opening_penalty;
                score_t ins_extend = E.back() + config_.gap_extension_penalty;
                score_t ins_score = std::max(ins_open, ins_extend);

                if (ins_score >= xdrop_cutoff) {
                    E.push_back(ins_score);
                    S.push_back(ins_score);
                    F.push_back(ninf);
                    OE.push_back(ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH);
                    OS.push_back(Cigar::INSERTION);
                    OF.push_back(Cigar::CLIPPED);
                } else {
                    break;
                }
            }

            if (offset + 1 >= begin && offset + 1 < end) {
                #pragma omp simd
                for (size_t j = begin; j < end; ++j) {
                    if (std::make_pair(S[j - trim], std::abs(static_cast<ssize_t>(max_pos) - static_cast<ssize_t>(offset + 1)))
                            > std::make_pair(S[max_pos - trim], std::abs(static_cast<ssize_t>(j) - static_cast<ssize_t>(offset + 1)))) {
                        max_pos = j;
                    }
                }
            } else {
                max_pos = std::max_element(S.begin(), S.end()) - S.begin() + trim;
            }

            if (S[max_pos - trim] >= xdrop_cutoff) {
                Ref next_score { S[max_pos - trim], table.size() - 1, max_pos, offset };

                if (S[max_pos - trim] - xdrop_cutoff > config_.xdrop)
                    xdrop_cutoff = S[max_pos - trim] - config_.xdrop;

                if (S[max_pos - trim] > std::get<0>(best_score))
                    best_score = next_score;

                if (update_seed_filter(table.size() - 1))
                    queue.emplace(next_score);

            } else {
                table.pop_back();
            }
        }
    }

    if (table.size()) {
        init_backtrack();
        tsl::hopscotch_set<size_t> prev_starts;

        std::vector<size_t> indices;
        indices.reserve(table.size());
        for (size_t i = 0; i < table.size(); ++i) {
            const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[i];
            if (OS[max_pos - trim] == Cigar::MATCH)
                indices.emplace_back(i);
        }

        std::sort(indices.begin(), indices.end(), [this](size_t i, size_t j) {
            const auto &a = table[i];
            const auto &b = table[j];
            return std::get<0>(a)[std::get<10>(a) - std::get<11>(a)]
                > std::get<0>(b)[std::get<10>(b) - std::get<11>(b)];
        });

        for (size_t i = 0; i < indices.size(); ++i) {
            if (skip_backtrack_start(extensions))
                break;

            size_t j = indices[i];
            if (!prev_starts.emplace(j).second)
                continue;

            std::vector<NodeType> path;
            std::vector<size_t> trace;
            Cigar ops;
            std::string seq;

            size_t pos = 0;
            size_t end_pos = 0;
            score_t score = 0;
            size_t align_offset = seed_->get_offset();

            {
                const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];
                if (S[max_pos - trim] < min_path_score)
                    break;

                if (offset < graph_.get_k() - 1)
                    continue;

                if (OS[max_pos - trim] != Cigar::MATCH)
                    continue;

                pos = max_pos;
                end_pos = max_pos;
                score = S[pos - trim];
            }

            auto get_extension = [&](Cigar cigar,
                                     std::vector<size_t> trace,
                                     std::vector<NodeType> final_path,
                                     std::string match) {
                assert(trace.size() == final_path.size());
                assert(final_path.size());
                cigar.append(Cigar::CLIPPED, pos);

                std::reverse(cigar.begin(), cigar.end());
                std::reverse(final_path.begin(), final_path.end());
                std::reverse(match.begin(), match.end());

                Alignment<NodeType> extension(
                    window.substr(pos, end_pos - pos),
                    std::move(final_path), std::move(match), score, std::move(cigar),
                    0, seed_->get_orientation(), align_offset
                );
                assert(extension.is_valid(graph_, &config_));

                trace.erase(trace.end() - extension.trim_offset(), trace.end());
                extension.extend_query_begin(query_.data());
                extension.extend_query_end(query_.data() + query_.size());
                assert(extension.is_valid(graph_, &config_));

                return std::make_pair(extension, trace);
            };

            while (j) {
                assert(j != static_cast<size_t>(-1));
                const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];
                assert(c == graph_.get_node_sequence(node)[std::min(graph_.get_k() - 1, offset)]);

                Cigar::Operator last_op = OS[pos - trim];
                align_offset = std::min(offset, graph_.get_k() - 1);

                if (S[pos - trim] == 0) {
                    auto [extension, trimmed_trace] = get_extension(ops, trace, path, seq);
                    process_extension(std::move(extension), std::move(trimmed_trace),
                                      prev_starts,
                                      [&](DBGAlignment&& extension) {
                        extensions.emplace_back(std::move(extension));
                    });
                }

                if (pos == max_pos)
                    prev_starts.emplace(j);

                switch (last_op) {
                    case Cigar::MATCH:
                    case Cigar::MISMATCH: {
                        if (offset >= graph_.get_k() - 1) {
                            path.emplace_back(node);
                            trace.emplace_back(j);
                        }

                        ops.append(last_op);
                        seq += c;
                        --pos;
                        assert(j_prev != static_cast<size_t>(-1));
                        j = j_prev;
                    } break;
                    case Cigar::INSERTION: {
                        while (last_op == Cigar::INSERTION) {
                            ops.append(last_op);
                            last_op = OE[pos - trim];
                            --pos;
                        }
                        assert(last_op == Cigar::MATCH);
                    } break;
                    case Cigar::DELETION: {
                        while (last_op == Cigar::DELETION && j) {
                            const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];

                            last_op = OF[pos - trim];
                            if (offset >= graph_.get_k() - 1) {
                                path.emplace_back(node);
                                trace.emplace_back(j);
                            }

                            seq += c;
                            ops.append(Cigar::DELETION);
                            assert(j_prev != static_cast<size_t>(-1));
                            j = j_prev;
                        }
                        assert(last_op == Cigar::MATCH || !j);
                    } break;
                    case Cigar::CLIPPED: { j = 0; } break;
                }
            }

            auto [extension, trimmed_trace]
                = get_extension(std::move(ops), std::move(trace),
                                std::move(path), std::move(seq));
            process_extension(std::move(extension), std::move(trimmed_trace),
                              prev_starts, [&](DBGAlignment&& extension) {
                extensions.emplace_back(std::move(extension));
            });
        }
    }

    return extensions;
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::skip_backtrack_start(const std::vector<DBGAlignment> &extensions) const {
    return extensions.size() >= config_.num_alternative_paths;
}

template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
