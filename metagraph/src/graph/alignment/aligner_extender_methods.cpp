#include "aligner_extender_methods.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"

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
        start_node_(graph_.max_index() + 1, '$', 0, 0),
        seed_filter_(graph_.get_k()) {
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

    if (seed_filter_.labels_to_keep(seed).size()) {
        seed_ = &seed;
        reset();
        xdrop_cutoff_ = ninf;
        table.clear();
    } else {
        // std::cerr << "skipping\t" << seed << "\n";
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

    // std::cerr << "starting\t" << *seed_ << "\t" << window.size() << "\n";

    size_t start = seed_->get_query().data() - query_.data();
    const score_t *partial = partial_sums_.data() + start;
    assert(*partial == config_.match_score(window));

    table.emplace_back(ScoreVec(window.size() + 1, ninf),
                       ScoreVec(window.size() + 1, ninf),
                       ScoreVec(window.size() + 1, ninf),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       OpVec(window.size() + 1, Cigar::CLIPPED),
                       graph_.max_index() + 1, -1, '\0', seed_->get_offset() - 1, 0);

    score_t xdrop_cutoff = -config_.xdrop;
    assert(config_.xdrop > 0);
    assert(xdrop_cutoff < 0);

    {
        auto &[S, E, F, OS, OE, OF, node, i_prev, c, offset, max_pos] = table[0];
        S[0] = 0;

        if (S.size() > 1) {
            S[1] = config_.gap_opening_penalty;
            E[1] = config_.gap_opening_penalty;
            OE[1] = Cigar::MATCH;
        }

        for (size_t j = 2; j < S.size(); ++j) {
            S[j] = config_.gap_opening_penalty + (j - 1) * config_.gap_extension_penalty;
            E[j] = S[j];
            OS[j] = Cigar::INSERTION;
            OE[j] = Cigar::INSERTION;

            if (S[j] < xdrop_cutoff)
                break;
        }
    }

    using Ref = std::pair<score_t, size_t>;
    Ref best_score { 0, 0 };

    std::priority_queue<Ref> queue;
    queue.emplace(best_score);

    size_t num_termina = 0;
    constexpr size_t max_num_termina = 10;



    while (queue.size() && num_termina < max_num_termina) {
        size_t i = queue.top().second;
        queue.pop();

        std::vector<std::pair<NodeType, char>> outgoing;
        size_t next_offset = -1;

        size_t begin = 0;
        size_t end = window.size() + 1;

        {
            const auto &[S, E, F, OS, OE, OF, node, i_prev, c, offset, max_pos] = table[i];
            next_offset = offset + 1;

            if (static_cast<double>(next_offset) / window.size() >= 2) {
                ++num_termina;
                continue;
            }

            begin = std::find_if(S.begin(), S.end(),
                                 [&](score_t s) { return s >= xdrop_cutoff; }) - S.begin();
            end = std::find_if(S.begin() + begin, S.end(),
                               [&](score_t s) { return s < xdrop_cutoff; }) - S.begin();

            if (end != S.size())
                ++end;

            if (end == begin) {
                ++num_termina;
                continue;
            }

            seed_filter_.update_seed_filter([&](const auto &callback) {
                callback(node, 0, begin, end);
            });

            if (config_.num_alternative_paths == 1) {
                bool has_extension = false;
                for (size_t j = begin; j < end; ++j) {
                    assert(partial[j] == config_.match_score(window.substr(j)));
                    if (S[j] + partial[j] >= best_score.first) {
                        has_extension = true;
                        break;
                    }
                }

                if (!has_extension) {
                    ++num_termina;
                    continue;
                }
            }

            if (next_offset - seed_->get_offset() < seed_->get_sequence().size()) {
                if (next_offset < graph_.get_k()) {
                    outgoing.emplace_back(seed_->front(), seed_->get_sequence()[next_offset - seed_->get_offset()]);
                } else {
                    outgoing.emplace_back((*seed_)[next_offset - graph_.get_k() + 1],
                                          seed_->get_sequence()[next_offset - seed_->get_offset()]);
                    assert(graph_.traverse(node, outgoing.back().second) == outgoing.back().first);
                }
            } else {
                assert(node != graph_.max_index() + 1);
                graph_.call_outgoing_kmers(node, [&](NodeType next, char c) {
                    outgoing.emplace_back(next, c);
                });
            }
        }

        for (const auto &[next, c] : outgoing) {
            assert(best_score.first > xdrop_cutoff);

            table.emplace_back(ScoreVec(window.size() + 1, ninf),
                               ScoreVec(window.size() + 1, ninf),
                               ScoreVec(window.size() + 1, ninf),
                               OpVec(window.size() + 1, Cigar::CLIPPED),
                               OpVec(window.size() + 1, Cigar::CLIPPED),
                               OpVec(window.size() + 1, Cigar::CLIPPED),
                               next, i, c, next_offset, 0);

            const auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev,
                         node_prev, i_prev, c_prev, offset_prev, max_pos_prev] = table[i];

            auto &[S, E, F, OS, OE, OF, node_cur, i_cur, c_stored, offset, max_pos] = table.back();

            // std::cerr << "\tcheck\t" << best_score.first << "\t" << table.size() << "," << offset << " "
            //           << max_pos_prev << "/" << window.size() << "\n";

            assert(i_cur == i);
            assert(node_cur == next);
            assert(c_stored == c);
            assert(offset == offset_prev + 1);
            assert(c == graph_.get_node_sequence(node_cur)[std::min(graph_.get_k() - 1, offset)]);

            auto update_del = [&](size_t j, score_t &s_score, Cigar::Operator &s_op) {
                score_t del_open = S_prev[j] + config_.gap_opening_penalty;
                score_t del_extend = F_prev[j] + config_.gap_extension_penalty;
                score_t del_score = std::max(del_open, del_extend);

                if (del_score >= xdrop_cutoff) {
                    F[j] = del_score;
                    OF[j] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
                    if (del_score > s_score) {
                        s_score = del_score;
                        s_op = Cigar::DELETION;
                    }
                }
            };

            if (begin == 0)
                update_del(0, S[0], OS[0]);

            const int8_t *profile_scores = profile_score_[c].data() + start;
            const Cigar::Operator *profile_ops = profile_op_[c].data() + start;
            size_t loop_begin = std::max(begin, (size_t)1);

            // update match and delete scores
            #pragma omp simd
            for (size_t j = loop_begin; j < end; ++j) {
                score_t s_score = ninf;
                Cigar::Operator s_op = Cigar::CLIPPED;

                // deletion
                update_del(j, s_score, s_op);

                score_t match = S_prev[j - 1] + profile_scores[j];
                if (match > std::max(s_score, xdrop_cutoff - 1)) {
                    s_score = match;
                    s_op = profile_ops[j];
                }

                S[j] = s_score;
                OS[j] = s_op;
            }

            // update insertion scores
            // since each element is dependent on the previous one, this can't
            // be vectorized easily
            #pragma omp simd
            for (size_t j = loop_begin; j < end; ++j) {
                score_t ins_open = S[j - 1] + config_.gap_opening_penalty;
                score_t ins_extend = E[j - 1] + config_.gap_extension_penalty;
                score_t ins_score = std::max(ins_open, ins_extend);

                if (ins_score >= xdrop_cutoff) {
                    E[j] = ins_score;
                    OE[j] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

                    if (ins_score > S[j]) {
                        S[j] = ins_score;
                        OS[j] = Cigar::INSERTION;
                    }
                }
            }

            max_pos = std::max_element(S.begin() + begin, S.begin() + end) - S.begin();

            if (S[max_pos] >= xdrop_cutoff) {
                Ref next_score { S[max_pos], table.size() - 1 };
                queue.emplace(next_score);

                if (S[max_pos] - xdrop_cutoff > config_.xdrop)
                    xdrop_cutoff = S[max_pos] - config_.xdrop;

                if (S[max_pos] > best_score.first) {
                    num_termina = 0;
                    best_score = next_score;
                }

            } else {
                table.pop_back();
                ++num_termina;
            }
        }
    }

    if (table.size()) {
        std::vector<size_t> indices;
        if (config_.num_alternative_paths == 1) {
            indices.push_back(best_score.second);
        } else {
            indices.resize(table.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [this](size_t i, size_t j) {
                const auto &a = table[i];
                const auto &b = table[j];
                return std::get<0>(a)[std::get<10>(a)] > std::get<0>(b)[std::get<10>(b)];
            });
        }

        for (size_t i = 0; i < indices.size(); ++i) {
            if (extensions.size() >= config_.num_alternative_paths)
                break;

            std::vector<NodeType> path;
            Cigar ops;
            std::string seq;

            size_t pos = 0;
            size_t end_pos = 0;
            score_t score = 0;

            size_t j = indices[i];
            {
                const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos] = table[j];
                if (S[max_pos] < min_path_score)
                    break;

                if (offset < graph_.get_k() - 1)
                    continue;

                if (OS[max_pos] != Cigar::MATCH && OS[max_pos] != Cigar::MISMATCH)
                    continue;

                pos = max_pos;
                end_pos = max_pos;
                score = S[pos];
            }

            while (pos || j) {
                assert(j != static_cast<size_t>(-1));
                const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos] = table[j];
                assert(c == graph_.get_node_sequence(node)[std::min(graph_.get_k() - 1, offset)]);

                Cigar::Operator last_op = OS[pos];
                switch (last_op) {
                    case Cigar::MATCH:
                    case Cigar::MISMATCH: {
                        if (offset >= graph_.get_k() - 1)
                            path.emplace_back(node);

                        ops.append(last_op);
                        seq += c;
                        --pos;
                        assert(j_prev != static_cast<size_t>(-1));
                        j = j_prev;
                    } break;
                    case Cigar::INSERTION: {
                        while (last_op == Cigar::INSERTION) {
                            ops.append(last_op);
                            last_op = OE[pos--];
                        }
                        assert(last_op == Cigar::MATCH);
                    } break;
                    case Cigar::DELETION: {
                        while (last_op == Cigar::DELETION && (pos || j)) {
                            const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos] = table[j];

                            last_op = OF[pos];
                            if (offset >= graph_.get_k() - 1)
                                path.emplace_back(node);

                            seq += c;
                            ops.append(Cigar::DELETION);
                            assert(j_prev != static_cast<size_t>(-1));
                            j = j_prev;
                        }
                        assert(last_op == Cigar::MATCH || (!pos && !j));
                    } break;
                    case Cigar::CLIPPED: { assert(false); } break;
                }
            }

            assert(path.size());

            std::reverse(ops.begin(), ops.end());
            std::reverse(path.begin(), path.end());
            std::reverse(seq.begin(), seq.end());

            Alignment<NodeType> extension(
                window.substr(0, end_pos),
                std::move(path), std::move(seq), score, std::move(ops),
                0, seed_->get_orientation(), seed_->get_offset()
            );
            assert(extension.is_valid(graph_, &config_));

            std::ignore = config_;
            extension.trim_offset();
            extension.extend_query_begin(query_.data());
            extension.extend_query_end(query_.data() + query_.size());
            assert(extension.is_valid(graph_, &config_));

            extensions.emplace_back(std::move(extension));
            // std::cerr << "ext\t" << extensions.back() << "\n";
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
