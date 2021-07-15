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

    template <typename Type, typename Sequence, typename Compare>
    class PriorityDeque : public boost::container::priority_deque<Type, Sequence, Compare> {
      public:
        void sort_interval_heap() {
            boost::heap::sort_interval_heap(this->sequence().begin(), this->sequence().end(),
                                            this->compare());
        }

        Sequence& data() { return this->sequence(); }
    };


  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;
    typedef PriorityDeque<std::shared_ptr<DBGAlignment>,
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
                         const std::function<bool()> &terminate);

    void call_alignments(const std::function<void(DBGAlignment&&)> &callback);

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

    void construct_alignment_chains();
    void construct_alignment_chain(std::string_view query,
                                   DBGAlignment&& chain,
                                   typename std::vector<DBGAlignment>::iterator begin,
                                   typename std::vector<DBGAlignment>::iterator end,
                                   std::vector<score_t> &best_score,
                                   const std::function<void(DBGAlignment&&)> &callback);
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
    if (config_.chain_alignments)
        construct_alignment_chains();

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
::call_alignments(const std::function<void(DBGAlignment&&)> &callback) {
    call_alignments(callback, []() { return false; });
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::construct_alignment_chains() {
    if (path_queue_.empty())
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

    auto push_to_queue = [&](DBGAlignment&& chain) {
        auto packaged_alignment = std::make_shared<DBGAlignment>(std::move(chain));

        auto add_to_target = [&](uint64_t target) {
            auto &cur_queue = path_queue_[target];

            for (const auto &aln : cur_queue) {
                if (*packaged_alignment == *aln)
                    return;
            }

            if (cur_queue.size() < config_.num_alternative_paths) {
                cur_queue.emplace(packaged_alignment);
            } else if (!cmp_(packaged_alignment, cur_queue.minimum())) {
                cur_queue.update(cur_queue.begin(), packaged_alignment);
            }
        };

        add_to_target(std::numeric_limits<uint64_t>::max());
        std::for_each(packaged_alignment->target_columns.begin(),
                      packaged_alignment->target_columns.end(),
                      add_to_target);
    };

    for (bool orientation : { false, true }) {
        auto &aln = alignments[orientation];

        // sort by endpoint (using beginning point and scores as tie-breakers)
        std::sort(aln.begin(), aln.end(), [](const auto &a, const auto &b) {
            return std::make_tuple(a.get_clipping() + a.get_query().size(),
                                   a.get_clipping(),
                                   b.get_score(),
                                   a.get_sequence().size())
                < std::make_tuple(b.get_clipping() + b.get_query().size(),
                                  b.get_clipping(),
                                  a.get_score(),
                                  b.get_sequence().size());
        });

        // recursively construct chains
        std::string_view this_query = get_query(orientation);
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = aln.begin(); it != aln.end(); ++it) {
            size_t end_pos = it->get_query().data() + it->get_query().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain(this_query, DBGAlignment(*it),
                                          it + 1, aln.end(),
                                          best_score, push_to_queue);
            }
        }
    }
}

template <typename NodeType, class AlignmentCompare>
inline void AlignmentAggregator<NodeType, AlignmentCompare>
::construct_alignment_chain(std::string_view query,
                            DBGAlignment&& chain,
                            typename std::vector<DBGAlignment>::iterator begin,
                            typename std::vector<DBGAlignment>::iterator end,
                            std::vector<score_t> &best_score,
                            const std::function<void(DBGAlignment&&)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query().data();
    const char *chain_end = chain_begin + chain.get_query().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    size_t k = graph_.get_k();
    score_t score = chain.get_score();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        if (it->get_offset())
            continue;

        const char *next_begin = it->get_query().data();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        const char *next_end = next_begin + it->get_query().size();

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        Vector<uint64_t> target_columns;
        if (chain.target_columns.size()) {
            std::set_intersection(it->target_columns.begin(), it->target_columns.end(),
                                  chain.target_columns.begin(),
                                  chain.target_columns.end(),
                                  std::back_inserter(target_columns));

            if (target_columns.empty())
                continue;
        }

        DBGAlignment aln(*it);

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, graph_, config_);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query().data(), graph_, config_)
            );

            if (aln.empty() || aln.get_sequence().size() < graph_.get_k()
                    || aln.get_cigar().front().first != Cigar::MATCH)
                continue;

            assert(aln.get_query().data() == chain.get_query().data() + chain.get_query().size());

            if (overlap < k - 1)
                aln.insert_gap_prefix(-overlap, graph_, config_);
        }

        score_t next_score = score + aln.get_score();
        if (next_score > best_score[next_end - query.data()]) {
            best_score[next_end - query.data()] = next_score;

            DBGAlignment next_chain(chain);
            next_chain.trim_end_clipping();
            next_chain.append(std::move(aln));
            assert(next_chain.get_score() == next_score);
            assert(next_chain.is_valid(graph_, &config_));
            if (next_chain.size()) {
                if (target_columns.size() == chain.target_columns.size())
                    called = true;

                next_chain.target_columns = std::move(target_columns);
                construct_alignment_chain(query, std::move(next_chain), it + 1,
                                          end, best_score, callback);
            }
        }
    }

    if (!called)
        callback(std::move(chain));
}


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_AGGREGATOR_HPP__
