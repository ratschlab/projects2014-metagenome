#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__

#include "aligner_helper.hpp"
#include "common/vectors/bitmap.hpp"


template <typename NodeType = typename DeBruijnGraph::node_index>
class Seeder {
  public:
    typedef Alignment<NodeType> Seed;

    virtual ~Seeder() {}

    virtual void initialize(std::string_view /* query string */,
                            bool /* orientation */) {}

    virtual void call_seeds(std::function<void(Seed&&)> callback) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
    virtual const std::string_view get_query() const = 0;
};


template <typename NodeType>
class ExactMapSeeder;


template <typename NodeType = typename DeBruijnGraph::node_index>
class SuffixSeeder : public Seeder<NodeType> {
  public:
    typedef typename Seeder<NodeType>::Seed Seed;

    SuffixSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    void initialize(std::string_view query, bool orientation);

    void call_seeds(std::function<void(Seed&&)> callback) const;

    const DeBruijnGraph& get_graph() const { return base_seeder_->get_graph(); }
    virtual const std::string_view get_query() const { return base_seeder_->get_query(); }
    const DBGAlignerConfig& get_config() const { return base_seeder_->get_config(); }

  private:
    std::unique_ptr<ExactMapSeeder<NodeType>> base_seeder_;
    std::vector<std::vector<NodeType>> alt_query_nodes;

    std::vector<NodeType>& get_query_nodes() { return base_seeder_->get_query_nodes(); }
    std::vector<uint8_t>& get_offsets() { return base_seeder_->get_offsets(); }
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactMapSeeder : public Seeder<NodeType> {
  friend SuffixSeeder<NodeType>;

  public:
    typedef typename Seeder<NodeType>::Seed Seed;

    ExactMapSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) { assert(config_.check_config_scores()); }

    virtual ~ExactMapSeeder() {}

    virtual const DeBruijnGraph& get_graph() const override final { return graph_; }
    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

    const std::string_view get_query() const { return query_; }
    const std::vector<NodeType>& get_query_nodes() const { return query_nodes_; }
    const std::vector<uint8_t>& get_offsets() const { return offsets_; }
    const std::vector<score_t>& get_partial_sums() const { return partial_sum_; }
    bool get_orientation() const { return orientation_; }

    virtual void initialize(std::string_view query, bool orientation);

  protected:
    std::vector<NodeType>& get_query_nodes() { return query_nodes_; }
    std::vector<uint8_t>& get_offsets() { return offsets_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;
    std::vector<NodeType> query_nodes_;
    std::vector<uint8_t> offsets_;
    bool orientation_;
    std::vector<score_t> partial_sum_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactSeeder : public ExactMapSeeder<NodeType> {
  public:
    typedef typename Seeder<NodeType>::Seed Seed;

    ExactSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : ExactMapSeeder<NodeType>(graph, config) {}

    void call_seeds(std::function<void(Seed&&)> callback) const;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class MEMSeeder : public ExactMapSeeder<NodeType> {
  friend SuffixSeeder<NodeType>;

  public:
    typedef typename Seeder<NodeType>::Seed Seed;

    MEMSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : ExactMapSeeder<NodeType>(graph, config) {}

    virtual ~MEMSeeder() {}

    void initialize(std::string_view query, bool orientation);

    void call_seeds(std::function<void(Seed&&)> callback) const;

    const bitmap& get_mem_terminator() const {
        assert(is_mem_terminus_.get());
        return *is_mem_terminus_;
    }

  protected:
    MEMSeeder(const DeBruijnGraph &graph,
              const DBGAlignerConfig &config,
              std::unique_ptr<bitmap>&& is_mem_terminus)
          : ExactMapSeeder<NodeType>(graph, config),
            is_mem_terminus_(std::move(is_mem_terminus)) {
        assert(is_mem_terminus_->size() == graph.max_index() + 1);
    }

    std::vector<uint8_t>& get_query_node_flags() { return query_node_flags_; }

  private:
    const std::unique_ptr<bitmap> is_mem_terminus_;
    std::vector<uint8_t> query_node_flags_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class UniMEMSeeder : public MEMSeeder<NodeType> {
  public:
    typedef typename Seeder<NodeType>::Seed Seed;

    UniMEMSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
        : MEMSeeder<NodeType>(graph, config,
                              std::make_unique<bitmap_lazy>(
                                      [&](auto i) {
                                          return graph.has_multiple_outgoing(i)
                                                  || graph.indegree(i) > 1;
                                      },
                                      graph.max_index() + 1)) {}
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class Extender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~Extender() {}

    virtual void
    operator()(const DBGAlignment &path,
               std::string_view query,
               std::function<void(DBGAlignment&&, NodeType)> callback,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &path) = 0;
    virtual void initialize_query(const std::string_view query) = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;

  protected:
    virtual void set_graph(const DeBruijnGraph &graph) = 0;
    virtual void reset() = 0;
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class DefaultColumnExtender : public Extender<NodeType> {
  public:
    typedef typename Extender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename Extender<NodeType>::node_index node_index;
    typedef typename Extender<NodeType>::score_t score_t;

    DefaultColumnExtender(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(&graph), config_(config) { assert(config_.check_config_scores()); }

    virtual void
    operator()(const DBGAlignment &path,
               std::string_view query,
               std::function<void(DBGAlignment&&, NodeType)> callback,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &) override {}

    virtual void initialize_query(const std::string_view query) override;

    virtual const DeBruijnGraph& get_graph() const override {
        assert(graph_);
        return *graph_;
    }

    virtual const DBGAlignerConfig& get_config() const override { return config_; }

  protected:
    virtual void set_graph(const DeBruijnGraph &graph) override { graph_ = &graph; }
    virtual void reset() override { dp_table.clear(); }

  private:
    typedef std::pair<NodeType, score_t> ColumnRef;
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;

    DPTable<NodeType> dp_table;

    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;
};


#endif // __DBG_ALIGNER_METHODS_HPP__
