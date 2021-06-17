#ifndef __DBG_EXTENDER_METHODS_HPP__
#define __DBG_EXTENDER_METHODS_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "aligner_alignment.hpp"
#include "common/aligned_vector.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename NodeType = uint64_t>
class IExtender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~IExtender() {}

    std::vector<DBGAlignment>
    get_extensions(const DBGAlignment &seed,
                   score_t min_path_score = std::numeric_limits<score_t>::min()) {
        return set_seed(seed) ? extend(min_path_score) : std::vector<DBGAlignment>{};
    }

    virtual void set_graph(const DeBruijnGraph &graph) = 0;

  protected:
    virtual const DBGAlignment& get_seed() const = 0;
    virtual bool set_seed(const DBGAlignment &seed) = 0;

    virtual std::vector<DBGAlignment> extend(score_t min_path_score) = 0;
};

template <typename NodeType = uint64_t>
class SeedFilteringExtender : public IExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;

    SeedFilteringExtender(std::string_view query) : query_size_(query.size()) {}

    virtual ~SeedFilteringExtender() {}

    virtual void set_graph(const DeBruijnGraph &) override { conv_checker_.clear(); }

  protected:
    const DBGAlignment *seed_ = nullptr;
    size_t query_size_;

    typedef AlignedVector<score_t> ScoreVec;
    tsl::hopscotch_map<NodeType, ScoreVec> conv_checker_;

    virtual const DBGAlignment& get_seed() const override final { return *seed_; }
    virtual bool set_seed(const DBGAlignment &seed) override;

    virtual bool update_seed_filter(node_index node,
                                    size_t query_start,
                                    const score_t *s_begin,
                                    const score_t *s_end);

    virtual bool filter_nodes(node_index node, size_t query_start, size_t query_end);
};


template <typename NodeType = uint64_t>
class DefaultColumnExtender : public SeedFilteringExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual void set_graph(const DeBruijnGraph &graph) override {
        SeedFilteringExtender<NodeType>::set_graph(graph);
        graph_ = &graph;
    }

  protected:
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;

    virtual std::vector<DBGAlignment> extend(score_t min_path_score) override;

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &extensions) const;

    virtual void process_extension(DBGAlignment&& extension,
                                   const std::vector<size_t> &trace,
                                   tsl::hopscotch_set<size_t> &prev_starts,
                                   score_t min_path_score,
                                   const std::function<void(DBGAlignment&&)> &callback);

    virtual void init_backtrack() {}

    virtual void call_outgoing(NodeType node,
                               size_t max_prefetch_distance,
                               const std::function<void(NodeType, char)> &callback);

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<score_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
