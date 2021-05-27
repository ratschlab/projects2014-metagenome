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
      : graph_(graph), config_(config), query_(query),
        start_node_(graph_.max_index() + 1, '$', 0, 0) {
    // std::cerr << "query\n";
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
        reset();
        xdrop_cutoff_ = ninf;
        table.clear();
    } else {
        DEBUG_LOG("Skipping seed: {}", seed);
        seed_ = nullptr;
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

    // std::cerr << "start\t" << *seed_ << "\t" << window << "\n";
    // std::cerr << "start\t" << *seed_ << "\n";

    size_t start = seed_->get_clipping();
    const score_t *partial = partial_sums_.data() + start;
    assert(*partial == config_.match_score(window));

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

    while (queue.size()) {
        size_t i = std::get<1>(queue.top());
        queue.pop();

        std::vector<std::pair<NodeType, char>> outgoing;
        size_t next_offset = -1;

        size_t begin = 0;
        size_t loop_begin = 0;
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

            if (end == begin)
                continue;

            loop_begin = std::max(begin, (size_t)1);

            bool has_extension = false;
            for (size_t j = begin; j < end; ++j) {
                assert(partial[j] == config_.match_score(window.substr(j)));
                score_t ext_score = S[j - trim] + partial[j];
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

            // std::cerr << "\t" << c;
            // for (const auto &[next, cc] : outgoing) {
            //     std::cerr << "," << cc;
            // }
            // std::cerr << "\t" << S[max_pos] << "," << max_pos << "," << next_offset;
            // // std::cerr << c << "," << outgoing.size() << "," << begin << "," << end;
            // for (score_t s : S) {
            //     if (s == ninf) {
            //         std::cerr << "\te";
            //     } else {
            //         std::cerr << "\t" << s;
            //     }
            // }
            // std::cerr << "\n";
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

            // std::cerr << "\tnext\t" << node_cur << " " << table.size() << " " << std::get<0>(best_score) << " " << offset << " " << S.size() << "\n";

            if (next != node_prev && !trim && !trim_prev)
                S[0] = S_prev[0];

            assert(i_cur == i);
            assert(node_cur == next);
            assert(c_stored == c);
            assert(offset == offset_prev + 1);
            assert(c == graph_.get_node_sequence(node_cur)[std::min(graph_.get_k() - 1, offset)]);

            auto update_del = [&](size_t j, score_t &s_score, Cigar::Operator &s_op) {
                if (j < trim_prev || j - trim_prev >= S_prev.size())
                    return;

                score_t del_open = S_prev[j - trim_prev] + config_.gap_opening_penalty;
                score_t del_extend = F_prev[j - trim_prev] + config_.gap_extension_penalty;
                score_t del_score = std::max(del_open, del_extend);

                if (del_score >= xdrop_cutoff) {
                    F[j - trim] = del_score;
                    OF[j - trim] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
                    if (del_score > s_score) {
                        s_score = del_score;
                        s_op = Cigar::DELETION;
                    }
                }
            };

            if (!trim && begin == 0 && next != node_prev)
                update_del(0, S[0], OS[0]);

            const int8_t *profile_scores = profile_score_[c].data() + start;
            const Cigar::Operator *profile_ops = profile_op_[c].data() + start;

            // update match and delete scores
            #pragma omp simd
            for (size_t j = loop_begin; j < end; ++j) {
                score_t s_score = ninf;
                Cigar::Operator s_op = Cigar::CLIPPED;

                // deletion
                update_del(j, s_score, s_op);

                if (j - 1 >= trim_prev && j - 1 - trim_prev < S_prev.size()) {
                    score_t match = S_prev[j - 1 - trim_prev] + profile_scores[j];
                    if (match > std::max(s_score, xdrop_cutoff - 1)) {
                        s_score = match;
                        s_op = profile_ops[j];
                    }
                }

                S[j - trim] = s_score;
                OS[j - trim] = s_op;
            }

            // update insertion scores
            // since each element is dependent on the previous one, this can't
            // be vectorized easily
            #pragma omp simd
            for (size_t j = loop_begin; j < end; ++j) {
                if (j - 1 >= trim && j - 1 - trim < S.size()) {
                    score_t ins_open = S[j - 1 - trim] + config_.gap_opening_penalty;
                    score_t ins_extend = E[j - 1 - trim] + config_.gap_extension_penalty;
                    score_t ins_score = std::max(ins_open, ins_extend);

                    if (ins_score >= xdrop_cutoff) {
                        E[j - trim] = ins_score;
                        OE[j - trim] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

                        if (ins_score > S[j - trim]) {
                            S[j - trim] = ins_score;
                            OS[j - trim] = Cigar::INSERTION;
                        }
                    }
                }
            }

            #pragma omp simd
            for (size_t j = begin; j < end; ++j) {
                if (std::make_pair(S[j - trim], std::abs(static_cast<ssize_t>(j) - static_cast<ssize_t>(offset + 1)))
                        > std::make_pair(S[max_pos - trim], std::abs(static_cast<ssize_t>(max_pos) - static_cast<ssize_t>(offset + 1)))) {
                    max_pos = j;
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

            // std::cerr << c << "," << outgoing.size() << "," << begin << "," << end;
            // for (score_t s : S) {
            //     if (s == ninf) {
            //         std::cerr << "\te";
            //     } else {
            //         std::cerr << "\t" << s;
            //     }
            // }
            // std::cerr << "\n";

            // TODO: find closest to diagonal
            max_pos = std::max_element(S.begin(), S.end()) - S.begin() + trim;
            // std::cerr << "\t" << max_pos << "\t" << S[max_pos] << "\n";

            if (S[max_pos - trim] >= xdrop_cutoff) {
                Ref next_score { S[max_pos - trim], table.size() - 1, max_pos, offset };

                if (S[max_pos - trim] - xdrop_cutoff > config_.xdrop)
                    xdrop_cutoff = S[max_pos - trim] - config_.xdrop;

                if (S[max_pos - trim] > std::get<0>(best_score)) {
                    best_score = next_score;
                    queue.emplace(next_score);
                } else {
                    auto it = conv_checker.find(next);
                    if (it == conv_checker.end()) {
                        ScoreVec &s_merged = conv_checker.emplace(next, ScoreVec(query_.size() + 1, ninf)).first.value();
                        std::copy(S.begin(), S.end(), s_merged.begin() + start + begin);
                        queue.emplace(next_score);
                    } else {
                        ScoreVec &s_merged = it.value();
                        bool converged = true;
                        for (size_t j = begin; j < end; ++j) {
                            if (S[j - trim] > s_merged[start + j]) {
                                converged = false;
                                s_merged[j + start] = S[j - trim];
                            }
                        }

                        if (!converged)
                            queue.emplace(next_score);
                    }
                }

            } else {
                table.pop_back();
            }
        }
    }

    if (table.size()) {
        std::vector<size_t> indices;
        tsl::hopscotch_set<size_t> prev_starts;
        if (config_.num_alternative_paths == 1) {
            indices.push_back(std::get<1>(best_score));
        } else {
            indices.resize(table.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [this](size_t i, size_t j) {
                const auto &a = table[i];
                const auto &b = table[j];
                return std::get<0>(a)[std::get<10>(a) - std::get<11>(a)]
                    > std::get<0>(b)[std::get<10>(b) - std::get<11>(b)];
            });
        }

        for (size_t i = 0; i < indices.size(); ++i) {
            if (extensions.size() >= config_.num_alternative_paths)
                break;

            size_t j = indices[i];
            if (!prev_starts.emplace(j).second)
                continue;

            std::vector<NodeType> path;
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

            auto add_extension = [&](Cigar cigar,
                                     std::vector<NodeType> final_path,
                                     std::string match) {
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

                extension.trim_offset();
                extension.extend_query_begin(query_.data());
                extension.extend_query_end(query_.data() + query_.size());
                assert(extension.is_valid(graph_, &config_));

                // std::cerr << "\text\t" << extension << "\n";
                // std::cerr << "\text\t" << extension.get_cigar().to_string() << "\n";

                extensions.emplace_back(std::move(extension));
            };

            while (j) {
                assert(j != static_cast<size_t>(-1));
                const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];
                assert(c == graph_.get_node_sequence(node)[std::min(graph_.get_k() - 1, offset)]);

                Cigar::Operator last_op = OS[pos - trim];
                align_offset = std::min(offset, graph_.get_k() - 1);

                if (S[pos - trim] == 0)
                    add_extension(ops, path, seq);

                if (pos == max_pos)
                    prev_starts.emplace(j);

                switch (last_op) {
                    case Cigar::MATCH:
                    case Cigar::MISMATCH: {
                        // std::cerr << "\t" << j << "," << pos << "," << node << "," << Cigar::opt_to_char(last_op);
                        if (offset >= graph_.get_k() - 1)
                            path.emplace_back(node);

                        ops.append(last_op);
                        seq += c;
                        --pos;
                        assert(j_prev != static_cast<size_t>(-1));
                        j = j_prev;
                    } break;
                    case Cigar::INSERTION: {
                        // std::cerr << "\t" << j << "," << pos << "," << node << "," << Cigar::opt_to_char(last_op);
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
                            // std::cerr << "\t" << j << "," << pos << "," << node << "," << Cigar::opt_to_char(last_op);

                            last_op = OF[pos - trim];
                            if (offset >= graph_.get_k() - 1)
                                path.emplace_back(node);

                            seq += c;
                            ops.append(Cigar::DELETION);
                            assert(j_prev != static_cast<size_t>(-1));
                            j = j_prev;
                        }
                        assert(last_op == Cigar::MATCH || !j);
                    } break;
                    case Cigar::CLIPPED: {
                        // std::cerr << "\t" << j << "," << pos << "," << node << "," << Cigar::opt_to_char(last_op);
                        j = 0;
                    } break;
                }
            }

            // std::cerr << "\n";

            add_extension(std::move(ops), std::move(path), std::move(seq));
        }
    }

    return extensions;
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::skip_backtrack_start(const std::vector<DBGAlignment> &, const AlignNode &) const {
    return true;
}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>
::backtrack(score_t, AlignNode,
            tsl::hopscotch_set<AlignNode, AlignNodeHash> &,
            std::vector<DBGAlignment> &) const -> std::vector<AlignNode> {
    return {};
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &) const {
    return;
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>::has_converged(const Column &, const Scores &) {
    return true;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::sanitize(Scores &, score_t) {}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>::get_outgoing(const AlignNode &) const
        -> std::vector<AlignNode> { return {}; }

template <typename NodeType>
auto DefaultColumnExtender<NodeType>
::get_band(const AlignNode &, const Column &, score_t) -> std::pair<size_t, size_t> { return {}; }

template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
