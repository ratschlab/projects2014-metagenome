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
                               bool /* orientation of seq */)> QueryCallback;
    typedef std::function<void(const QueryCallback&)> QueryGenerator;
    typedef std::function<void(std::string_view, DBGQueryAlignment&&)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const = 0;
    virtual DBGQueryAlignment align(const std::string_view query,
                                    bool query_orientation = false) const;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};

template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class SeedAndExtendAligner : public IDBGAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;

    virtual ~SeedAndExtendAligner() {}

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override;

    virtual const DeBruijnGraph& get_graph() const override = 0;
    virtual const DBGAlignerConfig& get_config() const override = 0;

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    // Generate seeds, then extend them
    void align_core(const std::string_view query,
                    ISeeder<node_index>&& seeder,
                    IExtender<node_index>&& extender,
                    const std::function<void(DBGAlignment&&)> &callback,
                    const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const;

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    void align_one_direction(DBGQueryAlignment &paths,
                             bool orientation_to_align,
                             ISeeder<node_index>&& seeder) const;

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    void align_best_direction(DBGQueryAlignment &paths,
                              ISeeder<node_index>&& seeder,
                              ISeeder<node_index>&& seeder_rc) const;

    // Align both forwards and backwards from a given seed. Procedure
    // 1. Given each seed, extend forward to produce an alignment A
    // 2. Reverse complement the alignment to get A', treated like a new seed
    // 3. Extend A' forwards
    // 4. Reverse complement A' to get the final alignment A''
    void align_both_directions(DBGQueryAlignment &paths,
                               ISeeder<node_index>&& seeder) const;

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    virtual void align_aggregate(DBGQueryAlignment &paths,
                                 const AlignmentGenerator &alignment_generator) const;
};


template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class DBGAligner : public SeedAndExtendAligner<Seeder, Extender, AlignmentCompare> {
  public:
    typedef SeedAndExtendAligner<Seeder, Extender, AlignmentCompare> BaseAligner;
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;

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
    typedef typename BaseAligner::AlignmentGenerator AlignmentGenerator;

  private:
    const DeBruijnGraph& graph_;
    const DBGAlignerConfig config_;
};


template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    generate_query([&](std::string_view header, std::string_view query, bool query_orientation) {
        DBGQueryAlignment paths(query, query_orientation);
        std::string_view forward = paths.get_query();
        std::string_view reverse = paths.get_query_reverse_complement();

        std::string_view this_query = query_orientation ? reverse : forward;
        assert(this_query == query);

        assert(get_config().num_alternative_paths);
        Seeder seeder(get_graph(),
                      this_query, // use this_query since paths stores a copy
                      query_orientation,
                      map_sequence_to_nodes(get_graph(), query),
                      get_config());

        if (get_graph().is_canonical_mode()) {
            assert(!query_orientation);
            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            align_both_directions(paths, std::move(seeder));
        } else if (get_config().forward_and_reverse_complement) {
            assert(!query_orientation);

            align_best_direction(
                paths,
                std::move(seeder),
                Seeder(get_graph(),
                       reverse,
                       !query_orientation,
                       map_sequence_to_nodes(get_graph(), reverse),
                       get_config())
            );
        } else {
            align_one_direction(paths, query_orientation, std::move(seeder));
        }

        callback(header, std::move(paths));
    });
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
::align_core(const std::string_view query,
             ISeeder<node_index>&& seeder,
             IExtender<node_index>&& extender,
             const std::function<void(DBGAlignment&&)> &callback,
             const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const {
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
        extender.initialize(seed);
        extender([&](DBGAlignment&& extension, auto start_node) {
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
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
::align_one_direction(DBGQueryAlignment &paths,
                      bool orientation_to_align,
                      ISeeder<node_index>&& seeder) const {
    std::string_view query = orientation_to_align ? paths.get_query_reverse_complement()
                                                  : paths.get_query();

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(query,
                   std::move(seeder),
                   Extender(get_graph(), get_config(), query),
                   alignment_callback,
                   get_min_path_score);
    });
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
::align_best_direction(DBGQueryAlignment &paths,
                       ISeeder<node_index>&& seeder,
                       ISeeder<node_index>&& seeder_rc) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query_reverse_complement();

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(forward,
                   std::move(seeder),
                   Extender(get_graph(), get_config(), forward),
                   alignment_callback,
                   get_min_path_score);

        align_core(reverse,
                   std::move(seeder_rc),
                   Extender(get_graph(), get_config(), reverse),
                   alignment_callback,
                   get_min_path_score);

    });
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
::align_both_directions(DBGQueryAlignment &paths, ISeeder<node_index>&& seeder) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query_reverse_complement();

    std::vector<DBGAlignment> reverse_seeds;

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
#ifndef NDEBUG
        mtg::common::logger->trace("Aligning forwards");
#endif

        // First get forward alignments
        align_core(
            forward,
            std::move(seeder),
            Extender(get_graph(), get_config(), forward),
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
                rev.reverse_complement(get_graph(), reverse);
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
            reverse,
            ManualSeeder<node_index>(std::move(reverse_seeds)),
            Extender(get_graph(), get_config(), reverse),
            [&](DBGAlignment&& path) {
                // If the path originated from a backwards alignment (forward alignment
                // of a reverse complement) and did not skip the first characters
                // (so it is unable to be reversed), change back to the forward orientation
                if (path.get_orientation()) {
                    auto forward_path = path;
                    forward_path.reverse_complement(get_graph(), forward);
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
}

template <class Seeder, class Extender, class AlignmentCompare>
inline void SeedAndExtendAligner<Seeder, Extender, AlignmentCompare>
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
