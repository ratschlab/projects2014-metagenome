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
    struct SharedPtrCmp {
        bool operator()(const std::shared_ptr<Alignment<NodeType>> &a,
                        const std::shared_ptr<Alignment<NodeType>> &b) const {
            return base_cmp_(*a, *b);
        }

        AlignmentCompare base_cmp_;
    };

  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef std::vector<Alignment<NodeType>> Chain;
    typedef typename DBGAlignment::score_t score_t;
    typedef boost::container::priority_deque<std::shared_ptr<DBGAlignment>,
                                             std::vector<std::shared_ptr<DBGAlignment>>,
                                             SharedPtrCmp> PathQueue;

    AlignmentAggregator(const DeBruijnGraph &graph,
                        std::string_view query,
                        std::string_view rc_query,
                        const DBGAlignerConfig &config)
          : query_(query), rc_query_(rc_query), config_(config), graph_(graph) {
        assert(config_.num_alternative_paths);
    }

    void add_alignment(DBGAlignment&& alignment);

    score_t get_min_path_score(const Vector<uint64_t> &targets) const;
    score_t get_max_path_score(const Vector<uint64_t> &targets) const;

    score_t get_min_path_score(uint64_t target = std::numeric_limits<uint64_t>::max()) const;
    score_t get_max_path_score(uint64_t target = std::numeric_limits<uint64_t>::max()) const;

    score_t get_min_path_score(const DBGAlignment &seed) const {
        return get_min_path_score(seed.target_columns);
    }

    score_t get_max_path_score(const DBGAlignment &seed) const {
        return get_max_path_score(seed.target_columns);
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

    size_t size() const {
        size_t size = 0;
        for (const auto &[target, queue] : path_queue_)
            size += queue.size();

        return size;
    }

    size_t num_targets() const { return path_queue_.size(); }

    bool empty() const { return path_queue_.empty(); }

    VectorMap<uint64_t, PathQueue>& data() { return path_queue_; }

    void clear() { path_queue_.clear(); }

    std::string_view get_query(bool is_reverse_complement) const {
        return is_reverse_complement ? rc_query_ : query_;
    }

  private:
    std::string_view query_;
    std::string_view rc_query_;
    const DBGAlignerConfig &config_;
    const DeBruijnGraph &graph_;
    VectorMap<uint64_t, PathQueue> path_queue_;
    SharedPtrCmp cmp_;
};


template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::add_alignment(DBGAlignment&& alignment) {
    auto packaged_alignment = std::make_shared<DBGAlignment>(std::move(alignment));

    auto add_to_target = [&](uint64_t target) {
        auto &cur_queue = path_queue_[target];

        for (const auto &aln : cur_queue) {
            if (*packaged_alignment == *aln)
                return;
        }

        if (config_.chain_alignments || cur_queue.size() < config_.num_alternative_paths) {
            cur_queue.emplace(packaged_alignment);
        } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
            cur_queue.update(cur_queue.begin(), packaged_alignment);
        }
    };

    add_to_target(std::numeric_limits<uint64_t>::max());
    std::for_each(packaged_alignment->target_columns.begin(),
                  packaged_alignment->target_columns.end(),
                  add_to_target);
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(const Vector<uint64_t> &targets) const -> score_t {
    score_t global_min = !config_.chain_alignments
        ? get_max_path_score() * config_.fraction_of_top
        : std::numeric_limits<score_t>::min();

    if (targets.empty())
        return std::max(global_min, get_min_path_score());

    score_t min_score = std::numeric_limits<score_t>::max();
    for (uint64_t target : targets) {
        if (min_score < global_min)
            break;

        min_score = std::min(min_score, get_min_path_score(target));
    }

    return std::max(global_min, min_score);
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_max_path_score(const Vector<uint64_t> &targets) const -> score_t {
    if (targets.empty())
        return get_max_path_score();

    score_t max_score = std::numeric_limits<score_t>::min();
    for (uint64_t target : targets) {
        max_score = std::max(max_score, get_max_path_score(target));
    }

    return max_score;
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_min_path_score(uint64_t target) const -> score_t {
    auto find = path_queue_.find(target);
    return config_.chain_alignments || find == path_queue_.end() || find->second.size() < config_.num_alternative_paths
        ? config_.min_path_score
        : find->second.minimum()->get_score();
}

template <typename NodeType, class AlignmentCompare>
inline auto AlignmentAggregator<NodeType, AlignmentCompare>
::get_max_path_score(uint64_t target) const -> score_t {
    auto find = path_queue_.find(target);
    return find == path_queue_.end() ? config_.min_path_score
                                     : find->second.maximum()->get_score();
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignments(const std::function<void(DBGAlignment&&)> &callback,
                  const std::function<bool()> &terminate) {
    typedef std::pair<uint64_t, PathQueue> queue_value;
    auto queues = const_cast<std::vector<queue_value>&&>(path_queue_.values_container());

    if (queues.empty())
        return;

    auto cmp = [this](const queue_value &a, const queue_value &b) {
        if (b.second.empty())
            return false;

        if (a.second.empty() || cmp_(a.second.maximum(), b.second.maximum()))
            return true;

        if (cmp_(b.second.maximum(), a.second.maximum()))
            return false;

        return a.second.maximum().get() < b.second.maximum().get();
    };

    std::make_heap(queues.begin(), queues.end(), cmp);

    auto begin = queues.begin();
    auto end = queues.end();
    std::shared_ptr<DBGAlignment> last_alignment;
    while (!terminate() && queues.size() && queues[0].second.size()) {
        if (queues[0].second.maximum().get() != last_alignment.get()) {
            if (last_alignment) {
                assert(last_alignment->size());
                callback(std::move(*last_alignment));
                *last_alignment = DBGAlignment();
            }

            last_alignment = queues[0].second.maximum();
        }

        queues[0].second.pop_maximum();
        std::pop_heap(queues.begin(), queues.end(), cmp);
        if (queues.back().second.empty()) {
            queues.pop_back();
            if (queues.empty())
                break;

            begin = queues.begin();
            end = queues.end();
        }

        if (--end == begin && begin->second.size()) {
            end = queues.end();
            std::make_heap(begin, end, cmp);
        }
    }

    if (last_alignment) {
        assert(last_alignment->size());
        callback(std::move(*last_alignment));
    }

    path_queue_.clear();
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::call_alignment_chains(const std::function<void(Chain&&, score_t)> &callback,
                        const std::function<bool()> &terminate) {
    if (path_queue_.empty() || terminate())
        return;

    const auto &all_queue_ = path_queue_.find(std::numeric_limits<uint64_t>::max())->second;
    if (all_queue_.empty())
        return;

    std::vector<DBGAlignment> alignments[2];
    for (const auto &alignment : all_queue_) {
        alignments[alignment->get_orientation()].push_back(*alignment);
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

        Vector<uint64_t> target_columns;
        if (it->target_columns.size()) {
            std::set_intersection(it->target_columns.begin(), it->target_columns.end(),
                                  chain.back().target_columns.begin(),
                                  chain.back().target_columns.end(),
                                  std::back_inserter(target_columns));
            if (target_columns.empty())
                continue;
        }

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
                next_chain.back().target_columns = std::move(target_columns);
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
                    next_chain.back().target_columns = std::move(target_columns);
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
                    next_chain.back().target_columns = std::move(target_columns);
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
