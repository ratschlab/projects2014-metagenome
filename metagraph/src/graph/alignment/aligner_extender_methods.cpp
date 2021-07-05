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

// to ensure that SIMD operations on arrays don't read out of bounds
constexpr size_t kPadding = 5;

#define WRAP_MEMBER(x) [this](auto&&... args) { return x(std::forward<decltype(args)>(args)...); }


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
        auto &p_score_row = profile_score_.emplace(c, query_.size() + kPadding).first.value();
        auto &p_op_row = profile_op_.emplace(c, query_.size() + kPadding).first.value();

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

    if (it != conv_checker_.end()) {
        size_t pos = seed.get_query().size() + seed.get_clipping() - 1;
        const auto &[start, vec] = it->second;
        if (pos < start || pos - start >= vec.size() || vec[pos - start] < seed.get_score())
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

    size_t size = s_end - s_begin;

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end()) {
        conv_checker_.emplace(node, ScoreVec(query_start, { s_begin, s_end }));
        return true;

    } else {
        auto &[start, vec] = it.value();
        if (query_start + size <= start) {
            vec.insert(vec.begin(), start - query_start, ninf);
            std::copy(s_begin, s_end, vec.begin());
            start = query_start;
            return true;

        } else if (query_start >= start + vec.size()) {
            vec.reserve(query_start + size - start);
            vec.insert(vec.end(), query_start - start - vec.size(), ninf);
            vec.insert(vec.end(), s_begin, s_end);
            return true;

        } else {
            // overlap
            if (query_start < start) {
                vec.insert(vec.begin(), start - query_start, ninf);
                start = query_start;
            }

            if (query_start + size > start + vec.size())
                vec.resize(query_start + size - start, ninf);

            bool converged = true;
            score_t *v = vec.data() + query_start - start;
            for (size_t j = 0; j < size; ++j) {
                if (s_begin[j] > v[j]) {
                    converged = false;
                    v[j] = s_begin[j];
                }
            }

            return !converged;
        }
    }
}

template <typename NodeType>
bool SeedFilteringExtender<NodeType>::filter_nodes(node_index node,
                                                   size_t query_start,
                                                   size_t query_end) {
    assert(query_end >= query_start);
    assert(query_end <= query_size_);
    constexpr score_t mscore = -ninf;
    size_t size = query_end - query_start;

    auto it = conv_checker_.find(node);
    if (it == conv_checker_.end()) {
        conv_checker_.emplace(
            node, ScoreVec(query_start, AlignedVector<score_t>(size, mscore))
        );
        return true;

    } else {
        auto &[start, vec] = it.value();
        if (query_start + size <= start) {
            vec.insert(vec.begin(), start - query_start, ninf);
            std::fill(vec.begin(), vec.begin() + size, mscore);
            start = query_start;
            return true;

        } else if (query_start >= start + vec.size()) {
            vec.reserve(query_start + size - start);
            vec.insert(vec.end(), query_start - start - vec.size(), ninf);
            vec.insert(vec.end(), size, mscore);
            return true;

        } else {
            // overlap
            if (query_start < start) {
                vec.insert(vec.begin(), start - query_start, ninf);
                start = query_start;
            }

            if (query_start + size > start + vec.size())
                vec.resize(query_start + size - start, ninf);

            bool converged = true;
            score_t *v = vec.data() + query_start - start;
            for (size_t j = 0; j < size; ++j) {
                if (mscore > v[j]) {
                    converged = false;
                    v[j] = mscore;
                }
            }

            return !converged;
        }
    }
}

template <class ScoreVec>
void update_column(size_t prev_end,
                   const score_t *S_prev_v,
                   const score_t *F_prev_v,
                   ScoreVec &S_v,
                   ScoreVec &F_v,
                   const score_t *profile_scores,
                   score_t xdrop_cutoff,
                   const DBGAlignerConfig &config_) {
#ifndef __SSE4_1__
    for (size_t j = 0; j < prev_end; ++j) {
        score_t match = j ? (S_prev_v[j - 1] + profile_scores[j]) : ninf;
        score_t del_score = std::max(
            S_prev_v[j] + config_.gap_opening_penalty,
            F_prev_v[j] + config_.gap_extension_penalty
        );

        if (del_score >= xdrop_cutoff)
            F_v[j] = del_score;

        match = std::max(del_score, match);

        if (match >= xdrop_cutoff)
            S_v[j] = match;
    }
#else
    const __m128i gap_open = _mm_set1_epi32(config_.gap_opening_penalty);
    const __m128i gap_extend = _mm_set1_epi32(config_.gap_extension_penalty);
    const __m128i xdrop_v = _mm_set1_epi32(xdrop_cutoff - 1);
    const __m128i ninf_v = _mm_set1_epi32(ninf);
    const __m128i prev_end_v = _mm_set1_epi32(prev_end);
    __m128i j_v = _mm_set_epi32(3, 2, 1, 0);
    for (size_t j = 0; j < prev_end; j += 4) {
        // match = j ? S_prev_v[j - 1] + profile_scores[j] : ninf;
        __m128i match;
        if (j) {
            match = _mm_add_epi32(
                _mm_loadu_si128((__m128i*)&S_prev_v[j - 1]),
                _mm_loadu_si128((__m128i*)&profile_scores[j])
            );
        } else {
            // rotate elements to the right, then insert ninf in first cell
            match = _mm_shuffle_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), 0b10010000);
            match = _mm_add_epi32(match, _mm_loadu_si128((__m128i*)&profile_scores[j]));
            match = _mm_insert_epi32(match, ninf, 0);
        }

        // del_score = std::max(del_open, del_extend);
        __m128i del_score = _mm_max_epi32(
            _mm_add_epi32(_mm_loadu_si128((__m128i*)&S_prev_v[j]), gap_open),
            _mm_add_epi32(_mm_loadu_si128((__m128i*)&F_prev_v[j]), gap_extend)
        );

        // j < prev_end
        __m128i bound = _mm_cmpgt_epi32(prev_end_v, j_v);
        j_v = _mm_add_epi32(j_v, _mm_set1_epi32(4));

        // if (del_score >= xdrop_cutoff) F_v[j] = del_score
        __m128i del_mask = _mm_cmpgt_epi32(del_score, xdrop_v);
        del_mask = _mm_and_si128(del_mask, bound);
        del_score = _mm_blendv_epi8(ninf_v, del_score, del_mask);
        _mm_store_si128((__m128i*)&F_v[j], del_score);

        // match = max(match, del_score)
        match = _mm_max_epi32(match, del_score);

        // if (match >= xdrop_cutoff) S_v[j] = match
        __m128i mask = _mm_cmpgt_epi32(match, xdrop_v);
        mask = _mm_and_si128(mask, bound);
        match = _mm_blendv_epi8(ninf_v, match, mask);
        _mm_store_si128((__m128i*)&S_v[j], match);
    }
#endif

    if (S_v.size() > prev_end) {
        size_t j = S_v.size() - 1;
        score_t match = S_prev_v[j - 1] + profile_scores[j];
        if (match >= xdrop_cutoff)
            S_v[j] = match;
    }
}

template <class ScoreVec>
void update_ins(ScoreVec &S,
                ScoreVec &E,
                score_t xdrop_cutoff,
                const DBGAlignerConfig &config_) {
    // update insertion scores
    // elements are dependent on the previous one, so this can't be vectorized easily
    for (size_t j = 1; j < S.size(); ++j) {
        score_t ins_score = std::max(
            S[j - 1] + config_.gap_opening_penalty,
            E[j - 1] + config_.gap_extension_penalty
        );

        if (ins_score >= xdrop_cutoff) {
            E[j] = ins_score;
            S[j] = std::max(S[j], ins_score);
        }
    }
}

template <class ScoreVec>
void extend_ins(ScoreVec &S,
                ScoreVec &E,
                ScoreVec &F,
                size_t max_size,
                score_t xdrop_cutoff,
                const DBGAlignerConfig &config_) {
    while (S.back() >= xdrop_cutoff && S.size() < max_size) {
        score_t ins_score = std::max(
            S.back() + config_.gap_opening_penalty,
            E.back() + config_.gap_extension_penalty
        );

        if (ins_score >= xdrop_cutoff) {
            S.push_back(ins_score);
            E.push_back(ins_score);
            F.push_back(ninf);
        } else {
            break;
        }
    }
}

template <typename NodeType, class Table, class Skipper, class Processor, class ProfileScores, class ProfileOps>
std::vector<Alignment<NodeType>> backtrack(const Table &table,
                                           const Skipper &skip_backtrack_start,
                                           const Processor &process_extension,
                                           const DeBruijnGraph *graph_,
                                           const DBGAlignerConfig &config_,
                                           const ProfileScores &profile_scores,
                                           const ProfileOps &profile_ops,
                                           score_t min_path_score,
                                           const Alignment<NodeType> &seed,
                                           std::string_view query_,
                                           std::string_view window) {
    typedef Alignment<NodeType> DBGAlignment;
    std::vector<DBGAlignment> extensions;
    tsl::hopscotch_set<size_t> prev_starts;

    size_t seed_clipping = seed.get_clipping();
    size_t seed_offset = seed.get_offset();
    bool orientation = seed.get_orientation();

    std::vector<size_t> indices;
    indices.reserve(table.size());
    for (size_t i = 1; i < table.size(); ++i) {
        const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[i];
        const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

        if (max_pos < trim_p + 1)
            continue;

        size_t pos = max_pos - trim;
        size_t pos_p = max_pos - trim_p - 1;
        if (S[pos] >= min_path_score
                && offset >= graph_->get_k() - 1
                && S[pos] == S_p[pos_p] + profile_scores.find(c)->second[seed_clipping + max_pos]
                && profile_ops.find(c)->second[seed_clipping + max_pos] == Cigar::MATCH) {
            indices.emplace_back(i);
        }
    }

    std::sort(indices.begin(), indices.end(), [&table](size_t i, size_t j) {
        const auto &[S_a, E_a, F_a, node_a, j_prev_a, c_a, offset_a, max_pos_a, trim_a] = table[i];
        const auto &[S_b, E_b, F_b, node_b, j_prev_b, c_b, offset_b, max_pos_b, trim_b] = table[j];
        ssize_t off_diag_a = std::abs(static_cast<ssize_t>(max_pos_a) - static_cast<ssize_t>(offset_a));
        ssize_t off_diag_b = std::abs(static_cast<ssize_t>(max_pos_b) - static_cast<ssize_t>(offset_b));
        return std::make_tuple(S_a[max_pos_a - trim_a], off_diag_b, j)
            > std::make_tuple(S_b[max_pos_b - trim_b], off_diag_a, i);
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
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);

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

            return std::make_pair(std::move(extension), std::move(trace));
        };

        while (j) {
            assert(j != static_cast<size_t>(-1));
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

            assert(pos >= trim);
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);
            assert(c == graph_->get_node_sequence(node)[std::min(graph_->get_k() - 1, offset)]);

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

            if (S[pos - trim] == ninf) {
                j = 0;
            } else if (pos && pos >= trim_p + 1
                    && S[pos - trim] == S_p[pos - trim_p - 1]
                        + profile_scores.find(c)->second[seed_clipping + pos]) {
                // match/mismatch
                if (offset >= graph_->get_k() - 1) {
                    path.emplace_back(node);
                    trace.emplace_back(j);
                }

                seq += c;
                ops.append(profile_ops.find(c)->second[seed_clipping + pos]);
                --pos;
                assert(j_prev != static_cast<size_t>(-1));
                j = j_prev;

            } else if (S[pos - trim] == F[pos - trim] && ops.size() && ops.back().first != Cigar::INSERTION) {
                // deletion
                Cigar::Operator last_op = Cigar::DELETION;
                while (last_op == Cigar::DELETION && j) {
                    const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim] = table[j];
                    const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p] = table[j_prev];

                    assert(pos >= trim_p);

                    assert(F[pos - trim] == F_p[pos - trim_p] + config_.gap_extension_penalty
                        || F[pos - trim] == S_p[pos - trim_p] + config_.gap_opening_penalty);

                    last_op = F[pos - trim] == F_p[pos - trim_p] + config_.gap_extension_penalty
                        ? Cigar::DELETION
                        : Cigar::MATCH;

                    if (offset >= graph_->get_k() - 1) {
                        path.emplace_back(node);
                        trace.emplace_back(j);
                    }

                    seq += c;
                    ops.append(Cigar::DELETION);
                    assert(j_prev != static_cast<size_t>(-1));
                    j = j_prev;
                }
            } else if (pos && S[pos - trim] == E[pos - trim] && ops.size() && ops.back().first != Cigar::DELETION) {
                // insertion
                Cigar::Operator last_op = Cigar::INSERTION;
                while (last_op == Cigar::INSERTION) {
                    ops.append(last_op);

                    assert(E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        || E[pos - trim] == S[pos - trim - 1] + config_.gap_opening_penalty);

                    last_op = E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        ? Cigar::INSERTION
                        : Cigar::MATCH;

                    --pos;
                }
            } else {
                assert(false && "One of the above should apply");
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
    auto &[S, E, F, node_, i_prev_, c_, offset_, max_pos_, trim_] = column;

    S.reserve(size + kPadding);
    E.reserve(size + kPadding);
    F.reserve(size + kPadding);

    S.resize(size, ninf);
    E.resize(size, ninf);
    F.resize(size, ninf);

    std::fill(S.data() + S.size(), S.data() + S.capacity(), ninf);
    std::fill(E.data() + E.size(), E.data() + E.capacity(), ninf);
    std::fill(F.data() + F.size(), F.data() + F.capacity(), ninf);

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

    typedef AlignedVector<score_t> ScoreVec;
    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec,
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
            return std::make_tuple(a_score, std::abs(b_max_pos - b_offset), b_id)
                < std::make_tuple(b_score, std::abs(a_max_pos - a_offset), a_id);
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
                const auto &[S, E, F, node, i_prev, c, offset, max_pos, trim] = table[i];
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

                const auto &[S_prev, E_prev, F_prev,
                             node_prev, i_prev, c_prev, offset_prev, max_pos_prev, trim_prev] = table[i];

                auto &[S, E, F, node_cur, i_cur, c_stored, offset, max_pos, trim] = table.back();

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
                              S, F,
                              profile_score_[c].data() + start + trim,
                              xdrop_cutoff, config_);

                update_ins(S, E, xdrop_cutoff, config_);

                extend_ins(S, E, F, window.size() + 1 - trim, xdrop_cutoff, config_);

                assert(max_pos >= trim);
                assert(max_pos - trim < S.size());

                ssize_t s_offset = static_cast<ssize_t>(offset + 1);
                for (size_t j = 0; j < S.size(); ++j) {
                    if (std::make_pair(S[j], std::abs(static_cast<ssize_t>(max_pos) - static_cast<ssize_t>(s_offset)))
                            > std::make_pair(S[max_pos - begin], std::abs(static_cast<ssize_t>(j + begin) - static_cast<ssize_t>(s_offset)))) {
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
        return backtrack<NodeType>(table, WRAP_MEMBER(skip_backtrack_start),
                                   WRAP_MEMBER(process_extension), graph_, config_,
                                   profile_score_, profile_op_, min_path_score,
                                   *this->seed_, query_, window);
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
