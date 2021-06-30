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
      : SeedFilteringExtender<NodeType>(query),
        graph_(&graph), config_(config), query_(query) {
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

    for (char c : graph_->alphabet()) {
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
bool SeedFilteringExtender<NodeType>::set_seed(const DBGAlignment &seed) {
    assert(seed.size());
    assert(seed.get_cigar().size());
    assert(seed.get_cigar().back().first == Cigar::MATCH
        || seed.get_cigar().back().first == Cigar::MISMATCH);

    seed_ = nullptr;

    auto it = conv_checker_.find(seed.back());

    if (it != conv_checker_.end()
            && it->second.at(seed.get_query().size() + seed.get_clipping() - 1)
                < seed.get_score()) {
        it = conv_checker_.end();
    }

    if (it == conv_checker_.end()) {
        seed_ = &seed;
    } else {
        DEBUG_LOG("Skipping seed: {}", seed);
    }

    return seed_;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::process_extension(DBGAlignment&& extension,
                    const std::vector<size_t> &trace,
                    tsl::hopscotch_set<size_t> &prev_starts,
                    score_t min_path_score,
                    const std::function<void(DBGAlignment&&)> &callback) {
    prev_starts.insert(trace.begin(), trace.end());
    if (extension.get_score() >= min_path_score)
        callback(std::move(extension));
}

template <typename NodeType>
bool SeedFilteringExtender<NodeType>::update_seed_filter(node_index node,
                                                         size_t query_start,
                                                         const score_t *s_begin,
                                                         const score_t *s_end) {
    assert(s_end >= s_begin);
    assert(query_start + (s_end - s_begin) <= query_size_);
    auto it = conv_checker_.find(node);
    if (it == conv_checker_.end()) {
        ScoreVec &s_merged = conv_checker_.emplace(
            node, ScoreVec(query_size_, ninf)
        ).first.value();
        std::copy(s_begin, s_end, s_merged.begin() + query_start);
        return true;
    } else {
        score_t *s_merged = it.value().data() + query_start;
        bool converged = true;
        size_t size = s_end - s_begin;
        for (size_t j = 0; j < size; ++j) {
            if (s_begin[j] > s_merged[j]) {
                converged = false;
                s_merged[j] = s_begin[j];
            }
        }

        return !converged;
    }
}

template <typename NodeType>
bool SeedFilteringExtender<NodeType>::filter_nodes(node_index node,
                                                   size_t query_start,
                                                   size_t query_end) {
    assert(query_end >= query_start);
    assert(query_end <= query_size_);
    constexpr score_t mscore = -ninf;

    auto it = conv_checker_.find(node);
    if (it == conv_checker_.end()) {
        ScoreVec &s_merged = conv_checker_.emplace(
            node, ScoreVec(query_size_, ninf)
        ).first.value();
        std::fill(s_merged.begin() + query_start, s_merged.begin() + query_end, mscore);
        return true;
    } else {
        score_t *s_merged = it.value().data();
        bool converged = true;
        for (size_t j = query_start; j < query_end; ++j) {
            if (mscore > s_merged[j]) {
                converged = false;
                s_merged[j] = mscore;
            }
        }

        return !converged;
    }
}

template <class ScoreVec, class OpVec>
void update_column(size_t prev_end,
                   const score_t *S_prev_v,
                   const score_t *F_prev_v,
                   ScoreVec &S_v,
                   ScoreVec &F_v,
                   ScoreVec &E_v,
                   OpVec &OS_v,
                   OpVec &OF_v,
                   OpVec &OE_v,
                   const score_t *profile_scores,
                   const Cigar::Operator *profile_ops,
                   score_t xdrop_cutoff,
                   const DBGAlignerConfig &config_) {
    {
        size_t j = 0;
        score_t del_open = S_prev_v[j] + config_.gap_opening_penalty;
        score_t del_extend = F_prev_v[j] + config_.gap_extension_penalty;
        score_t del_score = std::max(del_open, del_extend);

        if (del_score >= xdrop_cutoff) {
            F_v[j] = del_score;
            S_v[j] = del_score;
            OF_v[j] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
            OS_v[j] = Cigar::DELETION;
        }
    }

#ifndef __AVX2__
    for (size_t j = 1; j < prev_end; ++j) {
        score_t match = S_prev_v[j - 1] + profile_scores[j];
        Cigar::Operator op = profile_ops[j];

        score_t del_open = S_prev_v[j] + config_.gap_opening_penalty;
        score_t del_extend = F_prev_v[j] + config_.gap_extension_penalty;
        score_t del_score = std::max(del_open, del_extend);

        if (del_score >= xdrop_cutoff) {
            F_v[j] = del_score;
            OF_v[j] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
            if (del_score > match) {
                match = del_score;
                op = Cigar::DELETION;
            }
        }

        if (match >= xdrop_cutoff) {
            S_v[j] = match;
            OS_v[j] = op;
        }
    }
#else
    for (size_t j = 1; j < prev_end; j += 8) {
        // match = S_prev_v[j - 1] + profile_scores[j];
        __m256i match = _mm256_loadu_si256((__m256i*)&S_prev_v[j - 1]);
        __m256i profile_v = _mm256_loadu_si256((__m256i*)&profile_scores[j]);
        __m128i ops = mm_loadu_si64(&profile_ops[j]);

        match = _mm256_add_epi32(match, profile_v);


        // del_open = S_prev_v[j] + config_.gap_opening_penalty;
        // del_extend = F_prev_v[j] + config_.gap_extension_penalty;
        // del_score = std::max(del_open, del_extend);
        __m256i del_open = _mm256_add_epi32(
            _mm256_loadu_si256((__m256i*)&S_prev_v[j]),
            _mm256_set1_epi32(config_.gap_opening_penalty)
        );
        __m256i del_extend = _mm256_add_epi32(
            _mm256_loadu_si256((__m256i*)&F_prev_v[j]),
            _mm256_set1_epi32(config_.gap_extension_penalty)
        );
        __m256i del_score = _mm256_max_epi32(del_open, del_extend);

        __m128i del_op = _mm_blendv_epi8(
            _mm_set1_epi8(Cigar::MATCH),
            _mm_set1_epi8(Cigar::DELETION),
            mm256_cvtepi32_epi8(_mm256_cmpgt_epi32(del_extend, del_open))
        );

        __m256i xdrop_v = _mm256_set1_epi32(xdrop_cutoff - 1);

        // del_score >= xdrop_cutoff
        __m256i del_mask = _mm256_cmpgt_epi32(del_score, xdrop_v);
        __m128i del_mask_small = mm256_cvtepi32_epi8(del_mask);

        del_score = _mm256_blendv_epi8(_mm256_set1_epi32(ninf), del_score, del_mask);
        del_op = _mm_blendv_epi8(_mm_set1_epi8(Cigar::CLIPPED), del_op, del_mask_small);

        _mm256_storeu_si256((__m256i*)&F_v[j], del_score);
        mm_storeu_si64((int8_t*)&OF_v[j], del_op);

        // match = max(match, del_score)
        __m256i match_mask = _mm256_cmpgt_epi32(del_score, match);
        __m128i match_mask_small = mm256_cvtepi32_epi8(match_mask);
        match = _mm256_max_epi32(match, del_score);
        ops = _mm_blendv_epi8(ops, _mm_set1_epi8(Cigar::DELETION), match_mask_small);

        // match >= xdrop_cutoff
        __m256i mask = _mm256_cmpgt_epi32(match, xdrop_v);
        __m128i mask_small = mm256_cvtepi32_epi8(mask);

        match = _mm256_blendv_epi8(_mm256_set1_epi32(ninf), match, mask);
        ops = _mm_blendv_epi8(_mm_set1_epi8(Cigar::CLIPPED), ops, mask_small);

        _mm256_storeu_si256((__m256i*)&S_v[j], match);
        mm_storeu_si64((int8_t*)&OS_v[j], ops);
    }
    _mm256_zeroupper();
#endif

    if (S_v.size() > prev_end) {
        size_t j = S_v.size() - 1;
        score_t match = S_prev_v[j - 1] + profile_scores[j];
        if (match >= xdrop_cutoff) {
            S_v[j] = match;
            OS_v[j] = profile_ops[j];
        }
    }

    // update insertion scores
    // since each element is dependent on the previous one, this can't
    // be vectorized easily
    for (size_t j = 1; j < S_v.size(); ++j) {
        score_t ins_open = S_v[j - 1] + config_.gap_opening_penalty;
        score_t ins_extend = E_v[j - 1] + config_.gap_extension_penalty;
        score_t ins_score = std::max(ins_open, ins_extend);

        if (ins_score >= xdrop_cutoff) {
            E_v[j] = ins_score;
            OE_v[j] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

            if (ins_score > S_v[j]) {
                S_v[j] = ins_score;
                OS_v[j] = Cigar::INSERTION;
            }
        }
    }
}

template <class ScoreVec, class OpVec>
void extend_ins(ScoreVec &S,
                ScoreVec &E,
                ScoreVec &F,
                OpVec &OS,
                OpVec &OE,
                OpVec &OF,
                size_t max_size,
                score_t xdrop_cutoff,
                const DBGAlignerConfig &config_) {
    while (S.back() >= xdrop_cutoff && S.size() < max_size) {
        score_t ins_open = S.back() + config_.gap_opening_penalty;
        score_t ins_extend = E.back() + config_.gap_extension_penalty;
        score_t ins_score = std::max(ins_open, ins_extend);

        if (ins_score >= xdrop_cutoff) {
            S.push_back(ins_score);
            E.push_back(ins_score);
            F.push_back(ninf);
            OS.push_back(Cigar::INSERTION);
            OE.push_back(ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH);
            OF.push_back(Cigar::CLIPPED);
        } else {
            break;
        }
    }
}

template <typename NodeType, class Table, class Skipper, class Processor>
std::vector<Alignment<NodeType>> backtrack(const Table &table,
                                           const Skipper &skip_backtrack_start,
                                           const Processor &process_extension,
                                           const DeBruijnGraph *graph_,
                                           score_t min_path_score,
                                           size_t seed_offset,
                                           bool orientation,
                                           std::string_view query_,
                                           std::string_view window) {
    typedef Alignment<NodeType> DBGAlignment;
    std::vector<DBGAlignment> extensions;
    tsl::hopscotch_set<size_t> prev_starts;

    std::vector<size_t> indices;
    indices.reserve(table.size());
    for (size_t i = 0; i < table.size(); ++i) {
        const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[i];
        if (OS[max_pos - trim] == Cigar::MATCH)
            indices.emplace_back(i);
    }

    std::sort(indices.begin(), indices.end(), [&table](size_t i, size_t j) {
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

        std::vector<DeBruijnGraph::node_index> path;
        std::vector<size_t> trace;
        Cigar ops;
        std::string seq;

        size_t pos = 0;
        size_t end_pos = 0;
        score_t score = 0;
        size_t align_offset = seed_offset;

        {
            const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);

            if (S[max_pos - trim] < min_path_score)
                break;

            if (offset < graph_->get_k() - 1)
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
                0, orientation, align_offset
            );

            trace.erase(trace.end() - extension.trim_offset(), trace.end());
            extension.extend_query_begin(query_.data());
            extension.extend_query_end(query_.data() + query_.size());

            return std::make_pair(extension, trace);
        };

        while (j) {
            assert(j != static_cast<size_t>(-1));
            const auto &[S, E, F, OS, OE, OF, node, j_prev, c, offset, max_pos, trim] = table[j];
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);
            assert(c == graph_->get_node_sequence(node)[std::min(graph_->get_k() - 1, offset)]);

            Cigar::Operator last_op = OS[pos - trim];

            if (S[pos - trim] == 0 && ops.size() && (ops.back().first == Cigar::MATCH || ops.back().first == Cigar::MISMATCH)) {
                auto [extension, trimmed_trace] = get_extension(ops, trace, path, seq);
                process_extension(std::move(extension), std::move(trimmed_trace),
                                  prev_starts, min_path_score,
                                  [&](DBGAlignment&& extension) {
                    extensions.emplace_back(std::move(extension));
                });
            }

            align_offset = std::min(offset, graph_->get_k() - 1);

            if (pos == max_pos)
                prev_starts.emplace(j);

            switch (last_op) {
                case Cigar::MATCH:
                case Cigar::MISMATCH: {
                    if (offset >= graph_->get_k() - 1) {
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
                        assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);

                        last_op = OF[pos - trim];
                        if (offset >= graph_->get_k() - 1) {
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

        if (ops.size() && (ops.back().first == Cigar::MATCH || ops.back().first == Cigar::MISMATCH)) {
            auto [extension, trimmed_trace]
                = get_extension(std::move(ops), std::move(trace),
                                std::move(path), std::move(seq));
            process_extension(std::move(extension), std::move(trimmed_trace),
                              prev_starts, min_path_score, [&](DBGAlignment&& extension) {
                extensions.emplace_back(std::move(extension));
            });
        }
    }

    return extensions;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::call_outgoing(NodeType node,
                size_t /* max_prefetch_distance */,
                const std::function<void(NodeType, char)> &callback) {
    graph_->call_outgoing_kmers(node, [&](NodeType next, char c) {
        if (c != boss::BOSS::kSentinel)
            callback(next, c);
    });
}

template <class Column, typename NodeType>
Column alloc_column(size_t size,
                    NodeType node,
                    size_t i_prev,
                    char c,
                    size_t offset,
                    size_t max_pos,
                    size_t trim) {
    Column column;
    auto &[S, E, F, OS, OE, OF, node_, i_prev_, c_, offset_, max_pos_, trim_] = column;

    S.reserve(size + 9);
    E.reserve(size + 9);
    F.reserve(size + 9);
    OS.reserve(size + 9);
    OE.reserve(size + 9);
    OF.reserve(size + 9);

    S.resize(size);
    E.resize(size);
    F.resize(size);
    OS.resize(size);
    OE.resize(size);
    OF.resize(size);

    std::fill(S.data(), S.data() + S.capacity(), ninf);
    std::fill(E.data(), E.data() + E.capacity(), ninf);
    std::fill(F.data(), F.data() + F.capacity(), ninf);
    std::fill(OS.data(), OS.data() + OS.capacity(), Cigar::CLIPPED);
    std::fill(OE.data(), OE.data() + OE.capacity(), Cigar::CLIPPED);
    std::fill(OF.data(), OF.data() + OF.capacity(), Cigar::CLIPPED);

    node_ = node;
    i_prev_ = i_prev;
    c_ = c;
    offset_ = offset;
    max_pos_ = max_pos;
    trim_ = trim;

    return column;
}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>
::extend(score_t min_path_score) -> std::vector<DBGAlignment> {
    assert(this->seed_);

    typedef typename SeedFilteringExtender<NodeType>::ScoreVec ScoreVec;
    typedef AlignedVector<Cigar::Operator> OpVec;

    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec, OpVec, OpVec, OpVec,
                       NodeType, size_t, char, size_t, size_t, size_t> Column;
    typedef std::vector<Column> Table;
    Table table;

    size_t start = this->seed_->get_clipping();

    std::string_view window(this->seed_->get_query().data(),
                            query_.data() + query_.size() - this->seed_->get_query().data());

    table.emplace_back(alloc_column<Column, NodeType>(
        window.size() + 1, this->seed_->front(), -1, '\0', this->seed_->get_offset() - 1, 0, 0
    ));
    {
        auto &S = std::get<0>(table.back());
        std::fill(S.begin(), S.end(), 0);
    }

    score_t xdrop_cutoff = std::max(-config_.xdrop, ninf + 1);
    assert(config_.xdrop > 0);
    assert(xdrop_cutoff < 0);

    using Ref = std::tuple<score_t,
                           size_t /* table idx */,
                           ssize_t /* max pos */,
                           ssize_t /* offset */>;
    Ref best_score { 0, 0, 0, this->seed_->get_offset() - 1 };

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

        bool nofork = true;
        while (nofork) {
            nofork = false;

            std::vector<std::pair<NodeType, char>> outgoing;
            size_t next_offset = -1;

            size_t prev_begin = 0;
            size_t prev_end = window.size() + 1;

            {
                const auto &[S, E, F, OS, OE, OF, node, i_prev, c, offset, max_pos, trim] = table[i];
                next_offset = offset + 1;

                if (static_cast<double>(table.size()) / window.size() >= config_.max_nodes_per_seq_char)
                    continue;

                auto in_range = [xdrop_cutoff](score_t s) { return s >= xdrop_cutoff; };

                prev_begin = std::find_if(S.begin(), S.end(), in_range) - S.begin() + trim;
                prev_end = std::find_if(S.rbegin(), S.rend(), in_range).base() - S.begin() + trim;

                if (prev_end <= prev_begin)
                    continue;

                bool has_extension = false;
                for (size_t j = prev_begin; j < prev_end; ++j) {
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

                if (next_offset - this->seed_->get_offset() < this->seed_->get_sequence().size()) {
                    if (next_offset < graph_->get_k()) {
                        outgoing.emplace_back(this->seed_->front(), this->seed_->get_sequence()[next_offset - this->seed_->get_offset()]);
                    } else {
                        outgoing.emplace_back((*this->seed_)[next_offset - graph_->get_k() + 1],
                                              this->seed_->get_sequence()[next_offset - this->seed_->get_offset()]);
                        assert(graph_->traverse(node, outgoing.back().second) == outgoing.back().first);
                    }
                } else {
                    call_outgoing(node, window.size() + 1 - offset - S.size(),
                                  [&](NodeType next, char c) { outgoing.emplace_back(next, c); });
                }
            }

            size_t begin = prev_begin;
            size_t end = std::min(prev_end, window.size()) + 1;

            for (const auto &[next, c] : outgoing) {
                assert(std::get<0>(best_score) > xdrop_cutoff);

                table.emplace_back(alloc_column<Column, NodeType>(
                    end - begin, next, i, c, next_offset, begin, begin
                ));

                const auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev,
                             node_prev, i_prev, c_prev, offset_prev, max_pos_prev, trim_prev] = table[i];

                auto &[S, E, F, OS, OE, OF, node_cur, i_cur, c_stored, offset, max_pos, trim] = table.back();

                if (next != node_prev && !trim && !trim_prev)
                    S[0] = S_prev[0];

                assert(i_cur == i);
                assert(node_cur == next);
                assert(c_stored == c);
                assert(offset == offset_prev + 1);
                assert(c == graph_->get_node_sequence(node_cur)[std::min(graph_->get_k() - 1, offset)]);

                update_column(prev_end - trim,
                              S_prev.data() + trim - trim_prev,
                              F_prev.data() + trim - trim_prev,
                              S, F, E, OS, OF, OE,
                              profile_score_[c].data() + start + trim,
                              profile_op_[c].data() + start + trim,
                              xdrop_cutoff, config_);

                extend_ins(S, E, F, OS, OE, OF, window.size() + 1 - trim, xdrop_cutoff, config_);

                assert(max_pos >= trim);
                assert(max_pos - trim < S.size());

                ssize_t s_offset = static_cast<ssize_t>(offset + 1);
                for (size_t j = 0; j < S.size(); ++j) {
                    if (std::make_pair(S[j], std::abs(static_cast<ssize_t>(max_pos) - s_offset))
                            > std::make_pair(S[max_pos - begin], std::abs(static_cast<ssize_t>(j + begin) - s_offset))) {
                        max_pos = j + begin;
                    }
                }
                assert(max_pos >= trim);
                assert(max_pos - trim < S.size());

                score_t max_val = S[max_pos - trim];
                if (max_val >= xdrop_cutoff) {
                    Ref next_score { max_val, table.size() - 1, max_pos, offset };

                    if (max_val - xdrop_cutoff > config_.xdrop)
                        xdrop_cutoff = max_val - config_.xdrop;

                    if (max_val > std::get<0>(best_score))
                        best_score = next_score;

                    size_t vec_offset = start + begin;
                    score_t *s_begin = S.data();
                    score_t *s_end = S.data() + S.size();

                    // skip the first index since it corresponds to the position
                    // before the query start
                    if (!begin) {
                        ++s_begin;
                    } else {
                        --vec_offset;
                    }

                    assert(s_begin <= s_end);
                    assert(vec_offset + (s_end - s_begin) <= query_.size());

                    if (this->update_seed_filter(next, vec_offset, s_begin, s_end)) {
                        if (outgoing.size() == 1 && max_val >= std::get<0>(queue.top())) {
                            nofork = true;
                            i = table.size() - 1;
                        } else {
                            queue.emplace(next_score);
                        }
                    }

                } else {
                    table.pop_back();
                }
            }
        }
    }

    if (table.size()) {
        init_backtrack();
        return backtrack<NodeType>(table,
            [this](const std::vector<DBGAlignment> &extensions) {
                return skip_backtrack_start(extensions);
            },
            [this](DBGAlignment&& extension,
                   const std::vector<size_t> &trace,
                   tsl::hopscotch_set<size_t> &prev_starts,
                   score_t min_path_score,
                   const std::function<void(DBGAlignment&&)> &callback) {
                assert(extension.is_valid(*graph_, &config_));
                process_extension(std::move(extension), trace, prev_starts,
                                  min_path_score, callback);
            },
            graph_, min_path_score,
            this->seed_->get_offset(), this->seed_->get_orientation(),
            query_, window
        );
    } else {
        return {};
    }
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
