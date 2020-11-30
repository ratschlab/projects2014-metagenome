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

    virtual ~IDBGAligner() {}

    virtual DBGQueryAlignment align(const std::string_view query) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};

template <class Seeder = ExactMapSeeder<>,
          class Extender = DefaultColumnExtender<>>
class SeedAndExtendAligner : public IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~SeedAndExtendAligner() {}

    virtual DBGQueryAlignment align(const std::string_view query) const override final {
        assert(get_config().num_alternative_paths);
        if (get_graph().is_canonical_mode()) {
            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            return align_both_directions(query);
        } else {
            return get_config().forward_and_reverse_complement
                ? align_forward_and_reverse_complement(query)
                : align_one_direction(query, false);
        }
    }

    virtual const DeBruijnGraph& get_graph() const override = 0;
    virtual const DBGAlignerConfig& get_config() const override = 0;

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&)> SeedGenerator;
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    // Generate seeds, then extend them
    void align(const std::string_view query,
               const SeedGenerator &seed_generator,
               const std::function<void(DBGAlignment&&)> &callback,
               const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const;

    virtual Seeder build_seeder() const = 0;
    virtual Extender build_extender() const = 0;

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

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    DBGQueryAlignment align_forward_and_reverse_complement(const std::string_view query) const;

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    virtual void align_aggregate(DBGQueryAlignment &paths,
                                 const AlignmentGenerator &alignment_generator) const = 0;
};


template <class Seeder = ExactMapSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class DBGAligner : public SeedAndExtendAligner<Seeder, Extender> {
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

    virtual ~DBGAligner() {}

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&)> SeedGenerator;
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

  private:
    Seeder build_seeder() const override { return Seeder(graph_, config_); }
    Extender build_extender() const override { return Extender(graph_, config_); }

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    void align_aggregate(DBGQueryAlignment &paths,
                         const AlignmentGenerator &alignment_generator) const override;

    const DeBruijnGraph& graph_;
    const DBGAlignerConfig config_;
};


template <class Seeder, class Extender>
inline void SeedAndExtendAligner<Seeder, Extender>
::align(const std::string_view query,
        const SeedGenerator &seed_generator,
        const std::function<void(DBGAlignment&&)> &callback,
        const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const {
    auto extend = build_extender();
    extend.initialize_query(query);

    std::vector<DBGAlignment> seeds;
    seed_generator([&](DBGAlignment&& seed) {
        assert(seed.is_valid(get_graph(), &get_config()));
        seeds.emplace_back(std::move(seed));
    });

    for (auto &seed : seeds) {
        mtg::common::logger->trace("Seed: {}", seed);
        score_t min_path_score = get_min_path_score(seed);

        if (seed.get_query_end() == query.data() + query.size()) {
            if (seed.get_score() >= min_path_score) {
                seed.trim_offset();
                assert(seed.is_valid(get_graph(), &get_config()));
                mtg::common::logger->trace("Alignment: {}", seed);
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
                    mtg::common::logger->trace("Alignment: {}", seed);
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
                mtg::common::logger->trace("Alignment: {}", extension);
                callback(std::move(extension));
                return;
            }

            assert(extension.get_offset() == get_graph().get_k() - 1);
            auto next_path = seed;
            next_path.append(std::move(extension));
            next_path.trim_offset();
            assert(next_path.is_valid(get_graph(), &get_config()));

            mtg::common::logger->trace("Alignment: {}", next_path);
            callback(std::move(next_path));
            extended = true;
        }, min_path_score);

        // if !extended, then the seed was not extended because of early cutoff
    }
}

template <class Seeder, class Extender>
inline auto SeedAndExtendAligner<Seeder, Extender>
::align_one_direction(const std::string_view query,
                      bool orientation) const -> DBGQueryAlignment {
    auto seeder = build_seeder();
    DBGQueryAlignment paths(query);

    if (orientation) {
        std::swap(const_cast<std::string&>(paths.get_query()),
                  const_cast<std::string&>(paths.get_query_reverse_complement()));
    }

    const auto &query_alignment = orientation ? paths.get_query_reverse_complement()
                                              : paths.get_query();
    assert(query_alignment == query);

    seeder.initialize(query_alignment, orientation);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align(query_alignment, [&](const auto &callback) {
            seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
        }, alignment_callback, get_min_path_score);

    });

    return paths;
}

template <class Seeder, class Extender>
inline auto SeedAndExtendAligner<Seeder, Extender>
::align_both_directions(const std::string_view query) const -> DBGQueryAlignment {
    auto seeder = build_seeder();
    DBGQueryAlignment paths(query);

    seeder.initialize(paths.get_query(), false);
    std::vector<DBGAlignment> reverse_seeds;

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        mtg::common::logger->trace("Aligning forwards");

        // First get forward alignments
        align(paths.get_query(),
            [&](const auto &forward_seed_callback) {
                seeder.call_seeds([&](DBGAlignment&& seed) {
                    forward_seed_callback(std::move(seed));
                });
            },
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
                    mtg::common::logger->trace("Alignment cannot be reversed, returning");
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

        mtg::common::logger->trace("Aligning backwards");

        // Then use the reverse complements of the forward alignments as seeds
        align(paths.get_query_reverse_complement(),
            [&](const auto &reverse_seed_callback) {
                for (auto&& path : reverse_seeds) {
                    reverse_seed_callback(std::move(path));
                }
            },
            [&](DBGAlignment&& path) {
                // If the path originated from a backwards alignment (forward alignment
                // of a reverse complement) and did not skip the first characters
                // (so it is unable to be reversed), change back to the forward orientation
                if (path.get_orientation()) {
                    auto forward_path = path;
                    forward_path.reverse_complement(seeder.get_graph(), paths.get_query());
                    if (!forward_path.empty()) {
                        path = std::move(forward_path);
                    } else {
                        mtg::common::logger->trace("Backwards alignment cannot be reversed, returning");
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

template <class Seeder, class Extender>
inline auto SeedAndExtendAligner<Seeder, Extender>
::align_forward_and_reverse_complement(const std::string_view query) const -> DBGQueryAlignment {
    auto seeder = build_seeder();
    DBGQueryAlignment paths(query);
    const auto &forward = paths.get_query();
    const auto &reverse = paths.get_query_reverse_complement();

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        seeder.initialize(paths.get_query(), false);
        align(forward, [&](const auto &callback) {
            seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
        }, alignment_callback, get_min_path_score);

        seeder.initialize(paths.get_query_reverse_complement(), true);
        align(reverse, [&](const auto &callback) {
            seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
        }, alignment_callback, get_min_path_score);

    });

    return paths;
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void DBGAligner<Seeder, Extender, AlignmentCompare>
::align_aggregate(DBGQueryAlignment &paths,
                  const AlignmentGenerator &alignment_generator) const {
    AlignmentAggregator<node_index, AlignmentCompare> path_queue(
        paths.get_query(), paths.get_query_reverse_complement(), graph_, config_
    );

    alignment_generator(
        [&](DBGAlignment&& alignment) {
            alignment.trim_offset(&graph_);
            path_queue.add_alignment(std::move(alignment));
        },
        [&](const DBGAlignment &seed) { return path_queue.get_min_path_score(seed); }
    );

    path_queue.call_alignments([&](auto&& alignment) {
        assert(alignment.is_valid(graph_, &config_));
        paths.emplace_back(std::move(alignment));
    });
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
