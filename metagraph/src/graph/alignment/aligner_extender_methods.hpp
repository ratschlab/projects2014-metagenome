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

    virtual std::vector<DBGAlignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &seed) = 0;

    virtual void set_graph(const DeBruijnGraph &graph) = 0;

  protected:
    virtual const DBGAlignment& get_seed() const = 0;
};


template <typename NodeType = uint64_t>
class DefaultColumnExtender : public IExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual std::vector<DBGAlignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &seed) override;

    virtual void set_graph(const DeBruijnGraph &graph) override {
        graph_ = &graph;
        table.clear();
        conv_checker.clear();
    }

  protected:
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;

    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<Cigar::Operator> OpVec;

    typedef std::vector<std::tuple<ScoreVec, ScoreVec, ScoreVec, OpVec, OpVec, OpVec,
                                   NodeType, size_t, char, size_t, size_t, size_t>> Table;

    Table table;
    tsl::hopscotch_map<NodeType, ScoreVec> conv_checker;

    // the initial seed
    const DBGAlignment *seed_ = nullptr;

    size_t start;

    virtual const DBGAlignment& get_seed() const override { return *seed_; }

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &extensions) const;

    virtual void process_extension(DBGAlignment&& extension,
                                   const std::vector<size_t> &trace,
                                   tsl::hopscotch_set<size_t> &prev_starts,
                                   const std::function<void(DBGAlignment&&)> &callback);

    virtual bool update_seed_filter(size_t j);

    virtual void init_backtrack() {}

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
