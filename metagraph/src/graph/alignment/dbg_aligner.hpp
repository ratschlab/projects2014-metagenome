#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "aligner_aggregator.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {

class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;
    typedef std::function<void(std::string_view /* header */,
                               std::string_view /* seq */,
                               bool /* orientation */)> QueryCallback;
    typedef std::function<void(const QueryCallback&)> QueryGenerator;
    typedef std::function<void(std::string_view, DBGQueryAlignment&&)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const = 0;
    virtual DBGQueryAlignment align(const std::string_view query,
                                    bool orientation = false) const;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};

template <class ThisSeeder = ExactSeeder<>,
          class ThisExtender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class SeedAndExtendAligner : public IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~SeedAndExtendAligner() {}

    virtual void
    align_batch(const QueryGenerator &generate_query,
                const AlignmentCallback &callback) const override {
        generate_query([&](std::string_view header, std::string_view query, bool orientation) {
            assert(get_config().num_alternative_paths);
            if (get_graph().is_canonical_mode()) {
                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                callback(header, align_both_directions(query));
            } else {
                callback(header, align_one_direction(query, orientation));
            }
        });
    }

    virtual const DeBruijnGraph& get_graph() const override = 0;
    virtual const DBGAlignerConfig& get_config() const override = 0;

  protected:
    typedef const std::function<void(const std::function<void(Seeder<node_index>&&, Extender<node_index>&&)>&)> AlignmentCoreGenerator;
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    // Generate seeds, then extend them
    void align_core(const std::string_view query,
                    const AlignmentCoreGenerator &seed_generator,
                    const std::function<void(DBGAlignment&&)> &callback,
                    const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const;

    virtual std::shared_ptr<Seeder<node_index>> build_seeder() const = 0;
    virtual std::shared_ptr<Extender<node_index>> build_extender() const = 0;

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    DBGQueryAlignment align_one_direction(const std::string_view query,
                                          bool orientation) const;

    // Align both forwards and backwards from a given seed. Procedure
    // 1. Given each seed, extend forward to produce an alignment A
    // 2. Reverse complement the alignment to get A', treated like a new seed
    // 3. Extend A' forwards
    // 4. Reverse complement A' to get the final alignment A''
    DBGQueryAlignment align_both_directions(const std::string_view query) const;

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    virtual void align_aggregate(DBGQueryAlignment &paths,
                                 const AlignmentGenerator &alignment_generator) const;

    virtual AlignmentCoreGenerator build_alignment_core_generator(const std::string_view query,
                                                                  bool orientation) const;

    virtual AlignmentCoreGenerator
    build_alignment_core_generator_from_seeds(const std::string_view query,
                                              bool orientation,
                                              std::vector<DBGAlignment>&& seeds) const;
};


template <class ThisSeeder = ExactSeeder<>,
          class ThisExtender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class DBGAligner : public SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare> {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  protected:
    typedef typename SeedAndExtendAligner<ThisSeeder, ThisExtender>::AlignmentGenerator AlignmentGenerator;

  private:
    std::shared_ptr<Seeder<node_index>> build_seeder() const override {
        return std::make_shared<ThisSeeder>(graph_, config_);
    }

    std::shared_ptr<Extender<node_index>> build_extender() const override {
        return std::make_shared<ThisExtender>(graph_, config_);
    }

    const DeBruijnGraph& graph_;
    const DBGAlignerConfig config_;
};


template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline void SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::align_core(const std::string_view query,
             const AlignmentCoreGenerator &alignment_core_generator,
             const std::function<void(DBGAlignment&&)> &callback,
             const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const {
    alignment_core_generator([&](Seeder<node_index>&& seeder, Extender<node_index>&& extend) {
        std::vector<DBGAlignment> seeds;
        seeder.call_seeds([&](DBGAlignment&& seed) {
            assert(seed.is_valid(get_graph(), &get_config()));
            seeds.emplace_back(std::move(seed));
        });

        for (auto &seed : seeds) {
#ifndef NDEBUG
            mtg::common::logger->trace("Seed: {}", seed);
#endif
            score_t min_path_score = get_min_path_score(seed);

            if (seed.get_query_end() == query.data() + query.size()) {
                if (seed.get_score() >= min_path_score) {
                    seed.trim_offset();
                    assert(seed.is_valid(get_graph(), &get_config()));
#ifndef NDEBUG
                    mtg::common::logger->trace("Alignment: {}", seed);
#endif
                    callback(std::move(seed));
                }

                continue;
            }

            bool extended = false;
            extend.initialize(seed);
            extend([&](DBGAlignment&& extension, auto start_node) {
                if (!start_node && !extended) {
                    // no good extension found
                    if (seed.get_score() >= min_path_score) {
                        seed.extend_query_end(query.data() + query.size());
                        seed.trim_offset();
                        assert(seed.is_valid(get_graph(), &get_config()));
#ifndef NDEBUG
                        mtg::common::logger->trace("Alignment: {}", seed);
#endif
                        callback(std::move(seed));
                    }
                    extended = true;
                    return;
                }

                assert(extension.is_valid(get_graph(), &get_config()));
                extension.extend_query_end(query.data() + query.size());

                if (extension.get_clipping() || start_node != seed.back()) {
                    // if the extension starts at a different position
                    // from the seed end, then it's a new alignment
                    extension.extend_query_begin(query.data());
                    extension.trim_offset();
                    assert(extension.is_valid(get_graph(), &get_config()));
#ifndef NDEBUG
                    mtg::common::logger->trace("Alignment: {}", extension);
#endif
                    callback(std::move(extension));
                    return;
                }

                assert(extension.get_offset() == get_graph().get_k() - 1);
                auto next_path = seed;
                next_path.append(std::move(extension));
                next_path.trim_offset();
                assert(next_path.is_valid(get_graph(), &get_config()));

#ifndef NDEBUG
                mtg::common::logger->trace("Alignment: {}", next_path);
#endif
                callback(std::move(next_path));
                extended = true;
            }, min_path_score);

            // if !extended, then the seed was not extended because of early cutoff
        }
    });
}

template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline auto SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::align_one_direction(const std::string_view query,
                      bool orientation) const -> DBGQueryAlignment {
    DBGQueryAlignment paths(query);

    if (orientation) {
        std::swap(const_cast<std::string&>(paths.get_query()),
                  const_cast<std::string&>(paths.get_query_reverse_complement()));
    }

    const auto &query_alignment = orientation ? paths.get_query_reverse_complement()
                                              : paths.get_query();
    assert(query_alignment == query);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(query_alignment,
                   build_alignment_core_generator(query_alignment, orientation),
                   alignment_callback,
                   get_min_path_score);
    });

    return paths;
}

template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline auto SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::align_both_directions(const std::string_view query) const -> DBGQueryAlignment {
    DBGQueryAlignment paths(query);

    std::vector<DBGAlignment> reverse_seeds;

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
#ifndef NDEBUG
        mtg::common::logger->trace("Aligning forwards");
#endif

        // First get forward alignments
        align_core(
            paths.get_query(),
            build_alignment_core_generator(paths.get_query(), false),
            [&](DBGAlignment&& path) {
                score_t min_path_score = get_min_path_score(path);

                // If the alignment starts from the beginning of the query,
                // there's no sequence left for aligning backwards.
                if (!path.get_clipping()) {
                    if (path.get_score() >= min_path_score)
                        alignment_callback(std::move(path));

                    return;
                }

                auto rev = path;
                rev.reverse_complement(get_graph(), paths.get_query_reverse_complement());
                if (rev.empty()) {
#ifndef NDEBUG
                    mtg::common::logger->trace("Alignment cannot be reversed, returning");
#endif
                    if (path.get_score() >= min_path_score)
                        alignment_callback(std::move(path));

                    return;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(rev.get_end_clipping());
                rev.trim_end_clipping();
                assert(rev.is_valid(get_graph(), &get_config()));

                // Pass the reverse complement of the forward alignment
                // as a seed for extension
                reverse_seeds.emplace_back(std::move(rev));
            },
            [&](const auto&) {
                // ignore the min path score for the forward alignment,
                // since it may have a score that is too low before it is
                // extended backwards
                return get_config().min_cell_score;
            }
        );

#ifndef NDEBUG
        mtg::common::logger->trace("Aligning backwards");
#endif

        // Then use the reverse complements of the forward alignments as seeds
        align_core(
            paths.get_query_reverse_complement(),
            build_alignment_core_generator_from_seeds(paths.get_query_reverse_complement(),
                                                      true, std::move(reverse_seeds)),
            [&](DBGAlignment&& path) {
                // If the path originated from a backwards alignment (forward alignment
                // of a reverse complement) and did not skip the first characters
                // (so it is unable to be reversed), change back to the forward orientation
                if (path.get_orientation()) {
                    auto forward_path = path;
                    forward_path.reverse_complement(get_graph(), paths.get_query());
                    if (!forward_path.empty()) {
                        path = std::move(forward_path);
#ifndef NDEBUG
                    } else {
                        mtg::common::logger->trace("Backwards alignment cannot be reversed, returning");
#endif
                    }
                }

                assert(path.is_valid(get_graph(), &get_config()));
                alignment_callback(std::move(path));
            },
            get_min_path_score
        );
    });

    return paths;
}

template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline auto SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::build_alignment_core_generator(const std::string_view query,
                                 bool orientation) const -> AlignmentCoreGenerator {
    return [this,query,orientation](const auto &callback) {
        auto seeder = build_seeder();
        seeder->initialize(query, orientation);
        auto extender = build_extender();
        extender->initialize_query(query);

        callback(std::move(*seeder), std::move(*extender));
    };
}

template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline auto SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::build_alignment_core_generator_from_seeds(const std::string_view query,
                                            bool orientation,
                                            std::vector<DBGAlignment>&& seeds) const
        -> AlignmentCoreGenerator {
    return [this,query,orientation,s=std::move(seeds)](const auto &callback) mutable {
        ManualSeeder<node_index> seeder(get_graph(), get_config(), std::move(s));
        seeder.initialize(query, orientation);
        auto extender = build_extender();
        extender->initialize_query(query);

        callback(std::move(seeder), std::move(*extender));
    };
}

template <class ThisSeeder, class ThisExtender, class AlignmentCompare>
inline void SeedAndExtendAligner<ThisSeeder, ThisExtender, AlignmentCompare>
::align_aggregate(DBGQueryAlignment &paths,
                  const AlignmentGenerator &alignment_generator) const {
    AlignmentAggregator<node_index, AlignmentCompare> path_queue(
        paths.get_query(), paths.get_query_reverse_complement(), get_config()
    );

    alignment_generator(
        [&](DBGAlignment&& alignment) { path_queue.add_alignment(std::move(alignment)); },
        [&](const DBGAlignment &seed) { return path_queue.get_min_path_score(seed); }
    );

    path_queue.call_alignments([&](auto&& alignment) {
        assert(alignment.is_valid(get_graph(), &get_config()));
        paths.emplace_back(std::move(alignment));
    });
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
