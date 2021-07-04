#ifndef __ALIGNER_AGGREGATOR_HPP__
#define __ALIGNER_AGGREGATOR_HPP__

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/vector_map.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType, class AlignmentCompare>
class AlignmentAggregator {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef std::vector<Alignment<NodeType>> Chain;
    typedef typename DBGAlignment::score_t score_t;
    typedef boost::container::priority_deque<DBGAlignment,
                                             std::vector<DBGAlignment>,
                                             AlignmentCompare> PathQueue;

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(DBGAlignment&& alignment);

    score_t get_min_path_score() const;
    score_t get_max_path_score() const;

    score_t get_min_path_score(const DBGAlignment &) const {
        return get_min_path_score();
    }

    score_t get_max_path_score(const DBGAlignment &) const {
        return get_max_path_score();
    }


    const DBGAlignment& maximum() const { return path_queue_.maximum(); }
    void pop_maximum() { path_queue_.pop_maximum(); }

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                         const std::function<bool()> &terminate = []() { return false; });

    void call_alignment_chains(const std::function<void(Chain&&, score_t)> &callback,
                               const std::function<bool()> &terminate = []() { return false; });

    void
    construct_alignment_chain(std::string_view query,
                              std::vector<DBGAlignment>&& chain,
                              score_t score,
                              typename std::vector<DBGAlignment>::iterator begin,
                              typename std::vector<DBGAlignment>::iterator end,
                              std::vector<score_t> &best_score,
                              const std::function<void(Chain&&, score_t)> &callback);

    size_t size() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    void clear() { path_queue_.clear(); }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    const DeBruijnGraph &graph_;
    PathQueue path_queue_;
    AlignmentCompare cmp_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    if (std::find(path_queue_.begin(), path_queue_.end(), alignment) != path_queue_.end())
        return;

    if (config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths) {
        path_queue_.emplace(std::move(alignment));
    } else if (!cmp_(alignment, path_queue_.minimum())) {
        path_queue_.update(path_queue_.begin(), std::move(alignment));
    }
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score() const -> score_t {
    return config_.chain_alignments || path_queue_.size() < config_.num_alternative_paths
        ? config_.min_path_score
        : std::max(static_cast<score_t>(path_queue_.maximum().get_score() * config_.fraction_of_top),
                   path_queue_.minimum().get_score());
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_max_path_score() const -> score_t {
    return path_queue_.size() ? path_queue_.maximum().get_score() : config_.min_path_score;
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                  const std::function<bool()> &terminate) {
    while (!terminate() && path_queue_.size()) {
        callback(DBGAlignment(path_queue_.maximum()));
        path_queue_.pop_maximum();
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignment_chains(const std::function<void(Chain&&, score_t)> &callback,
                        const std::function<bool()> &terminate) {
    if (path_queue_.empty() || terminate())
        return;

    std::vector<DBGAlignment> alignments[2];
    for (const auto &alignment : path_queue_) {
        alignments[alignment.get_orientation()].push_back(alignment);
    }

    if (alignments[0].empty() && alignments[1].empty())
        return;

    path_queue_.clear();

    for (auto &aln : alignments) {
        std::sort(aln.begin(), aln.end(), [](const auto &a, const auto &b) {
            return a.get_clipping() + a.get_query().size()
                < b.get_clipping() + b.get_query().size();
        });
    }

    std::vector<score_t> best_score(query_.size() + 1, 0);
    Chain best_chain;
    score_t best_chain_score = 0;
    auto update_chain = [&](Chain&& chain, score_t score) {
        if (chain.empty())
            return;

        if (score > best_chain_score) {
            best_chain = std::move(chain);
            best_chain_score = score;
        }
    };

    for (auto it = alignments[0].begin(); it != alignments[0].end() && !terminate(); ++it) {
        const char *chain_end = it->get_query().data() + it->get_query().size();

        assert(!it->get_orientation());
        assert(query_.data() + it->get_clipping() == it->get_query().data());
        assert(chain_end + it->get_end_clipping() == query_.data() + query_.size());

        if (it->get_score() > best_score[chain_end - query_.data()]) {
            best_score[chain_end - query_.data()] = it->get_score();
            construct_alignment_chain(query_, { *it }, it->get_score(), it + 1,
                                      alignments[0].end(), best_score, update_chain);
        }
    }

    std::vector<score_t> best_score_rc(rc_query_.size() + 1, 0);
    for (auto it = alignments[1].begin(); it != alignments[1].end() && !terminate(); ++it) {
        const char *chain_end = it->get_query().data() + it->get_query().size();

        assert(it->get_orientation());
        assert(rc_query_.data() + it->get_clipping() == it->get_query().data());
        assert(chain_end + it->get_end_clipping() == rc_query_.data() + rc_query_.size());

        if (it->get_score() > best_score_rc[chain_end - rc_query_.data()]) {
            best_score_rc[chain_end - rc_query_.data()] = it->get_score();
            construct_alignment_chain(rc_query_, { *it }, it->get_score(), it + 1,
                                      alignments[1].end(), best_score_rc, update_chain);
        }
    }

    callback(std::move(best_chain), best_chain_score);
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::construct_alignment_chain(std::string_view query,
                            std::vector<DBGAlignment>&& chain,
                            score_t score,
                            typename std::vector<DBGAlignment>::iterator begin,
                            typename std::vector<DBGAlignment>::iterator end,
                            std::vector<score_t> &best_score,
                            const std::function<void(Chain&&, score_t)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.back().get_query().data();
    const char *chain_end = chain_begin + chain.back().get_query().size();
    if (begin == end || chain_end == query_.data() + query_.size()) {
        callback(std::move(chain), score);
        return;
    }

    size_t k = graph_.get_k();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        const char *next_begin = it->get_query().data();

        assert(chain_begin - chain.back().get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.back().get_orientation());

        const char *next_end = next_begin + it->get_query().size();

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        std::vector<DBGAlignment> next_chain;
        score_t next_score = 0;

        if (next_begin >= chain_end) {
            // Score assumes that k-1 dummy sink k-mers are added
            // e.g.,
            // k = 4
            // gap = 2
            // ATGCTATGCA
            //             ACGTACGACT
            assert(k >= 2);
            next_score = score + it->get_score()
                + config_.gap_opening_penalty
                + (k - 2) * config_.gap_extension_penalty;

            if (next_score > best_score[next_end - query.data()]) {
                best_score[next_end - query.data()] = next_score;
                next_chain = chain;
                next_chain.push_back(*it);
            }
        } else {
            // alignments overlap
            // num_added = k - 1 - matching_overlap
            // e.g.,
            // k = 4
            // overlap = 4
            // matching overlap = 2
            // ATGCTATGCA
            //       ACCAACGACT
            // first trim front of the incoming alignment
            // ATGCTATGCA
            //           ACGACT
            //       TGCA
            //        GCAA - added
            //         CAAC
            //          AACG
            //           ACGA
            assert(chain.back().get_end_clipping());
            size_t chain_last_match = (chain.back().get_cigar().end() - 2)->second;

            DBGAlignment aln(*it);
            size_t overlap = std::min(
                chain_last_match,
                aln.trim_query_prefix(chain_end - it->get_query().data(), graph_, config_)
            );

            if (aln.empty())
                continue;

            assert(aln.get_query().data()
                == chain.back().get_query().data() + chain.back().get_query().size());

            if ((aln.get_cigar().begin() + static_cast<bool>(aln.get_clipping()))->first != Cigar::MATCH)
                continue;

            next_score = score + aln.get_score();

            if (k - 1 <= overlap && aln.get_offset() == k - 1) {
                // they overlap on a node, so they can be joined
                assert(graph_.traverse(chain.back().back(), aln.get_sequence()[0])
                    == aln.front());

                if (next_score > best_score[next_end - query.data()]) {
                    best_score[next_end - query.data()] = next_score;

                    next_chain = chain;
                    next_chain.back().trim_end_clipping();
                    next_chain.back().append(std::move(aln));
                    assert(next_chain.back().is_valid(graph_, &config_));
                }
            } else {
                // they can't be joined since the overlap is too small

                next_score += config_.gap_opening_penalty
                    + (k - overlap - 2) * config_.gap_extension_penalty;

                if (next_score > best_score[next_end - query.data()]) {
                    best_score[next_end - query.data()] = next_score;
                    next_chain = chain;
                    aln.trim_offset();
                    aln.extend_query_begin(query.data());
                    next_chain.emplace_back(std::move(aln));
                }
            }

        }

        if (next_chain.size()) {
            called = true;
            construct_alignment_chain(query, std::move(next_chain), next_score,
                                      it + 1, end, best_score, callback);
        }
    }

    if (!called)
        callback(std::move(chain), score);
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
