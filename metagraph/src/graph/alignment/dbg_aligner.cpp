#include "dbg_aligner.hpp"

#include <tsl/hopscotch_map.h>

#include "common/algorithms.hpp"

namespace mtg {
namespace graph {
namespace align {


class SeedFilter {
  public:
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::node_index node_index;
    typedef std::function<void(node_index, uint64_t, size_t, size_t)> LabeledNodeRangeCallback;
    typedef std::function<void(const LabeledNodeRangeCallback&)> LabeledNodeRangeGenerator;

    SeedFilter(size_t k) : k_(k) {}
    Vector<uint64_t> labels_to_keep(const DBGAlignment &seed);
    void update_seed_filter(const LabeledNodeRangeGenerator &generator);

  private:
    size_t k_;
    tsl::hopscotch_map<node_index, std::pair<size_t, size_t>> visited_nodes_;
};


template <class AlignmentCompare = LocalAlignmentLess>
class SeedAndExtendAlignerCore {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;

    typedef std::function<void(DBGAlignment&&)> LocalAlignmentCallback;
    typedef std::function<score_t(const DBGAlignment&)> MinScoreComputer;
    typedef std::function<void(const LocalAlignmentCallback&,
                               const MinScoreComputer&)> AlignmentGenerator;

    typedef std::function<void(const ISeeder<node_index>&)> SeederCallback;
    typedef std::function<void(std::vector<DBGAlignment>&& /* rev_comp_seeds */,
                               const SeederCallback&)> SeederGenerator;

    template <typename... Args>
    SeedAndExtendAlignerCore(const DeBruijnGraph &graph,
                             const DBGAlignerConfig &config,
                             std::shared_ptr<SeedFilter> seed_filter,
                             Args&&... args)
          : graph_(graph), config_(config),
            seed_filter_(seed_filter),
            paths_(std::forward<Args>(args)...),
            aggregator_(paths_.get_query(false), paths_.get_query(true), config_) {}

    void flush(const std::function<bool(const DBGAlignment&)> &skip
                   = [](const auto &) { return false; },
               const std::function<bool()> &terminate = []() { return false; }) {
        aggregator_.call_alignments([&](auto&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            if (!skip(alignment))
                paths_.emplace_back(std::move(alignment));
        }, terminate);
    }

    DBGQueryAlignment& get_paths() { return paths_; }

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    void align_one_direction(bool orientation_to_align,
                             const ISeeder<node_index> &seeder,
                             IExtender<node_index> &extender);

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    void align_best_direction(const ISeeder<node_index> &seeder,
                              const ISeeder<node_index> &seeder_rc,
                              IExtender<node_index> &extender,
                              IExtender<node_index> &extender_rc);

    // Align the forward and reverse complement of the query sequence in both
    // directions and return the overall best alignment. e.g., for the forward query
    // 1. Find all seeds of its reverse complement
    // 2. Given a seed, extend forwards to get alignment A
    // 3. Reverse complement the alignment to get A', treat it like a new seed
    // 4. Extend A' forwards to get the final alignment A''
    void align_both_directions(const ISeeder<node_index> &forward_seeder,
                               const ISeeder<node_index> &reverse_seeder,
                               IExtender<node_index> &forward_extender,
                               IExtender<node_index> &reverse_extender,
                               const SeederGenerator &rev_comp_core_generator);

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    void align_aggregate(const AlignmentGenerator &alignment_generator);

    AlignmentAggregator<node_index, AlignmentCompare>& get_aggregator() {
        return aggregator_;
    }

  protected:
    // Generate seeds, then extend them
    void align_core(std::string_view query,
                    const ISeeder<node_index> &seeder,
                    IExtender<node_index> &extender,
                    const LocalAlignmentCallback &callback,
                    const MinScoreComputer &get_min_path_score);

    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;

    std::shared_ptr<SeedFilter> seed_filter_;
    DBGQueryAlignment paths_;
    AlignmentAggregator<node_index, AlignmentCompare> aggregator_;
};

IDBGAligner::DBGQueryAlignment IDBGAligner::align(std::string_view query,
                                                  bool is_reverse_complement) const {
    DBGQueryAlignment result(query);
    std::string empty_header;
    align_batch(
        [&](const QueryCallback &callback) {
            callback(empty_header, query, is_reverse_complement);
        },
        [&](std::string_view, DBGQueryAlignment&& alignment) {
            result = std::move(alignment);
        }
    );

    return result;
}

void IDBGAligner
::align_batch(const std::vector<std::pair<std::string, std::string>> &seq_batch,
              const AlignmentCallback &callback) const {
    align_batch([&](const QueryCallback &query_callback) {
        for (const auto &[header, seq] : seq_batch) {
            query_callback(header, seq, false /* orientation of seq */);
        }
    }, callback);
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        SeedAndExtendAlignerCore<AlignmentCompare> aligner_core(
            get_graph(), get_config(),
            std::make_shared<SeedFilter>(get_graph().get_k()),
            query, is_reverse_complement
        );
        auto &paths = aligner_core.get_paths();
        std::string_view this_query = paths.get_query(is_reverse_complement);
        std::string_view reverse = paths.get_query(!is_reverse_complement);
        assert(this_query == query);

        std::vector<node_index> nodes = map_sequence_to_nodes(get_graph(), query);
        std::vector<node_index> nodes_rc;
        if (get_graph().get_mode() == DeBruijnGraph::CANONICAL
                || get_config().forward_and_reverse_complement) {
            assert(!is_reverse_complement);
            std::string dummy(query);
            nodes_rc = nodes;
            reverse_complement_seq_path(get_graph(), dummy, nodes_rc);
            assert(dummy == paths.get_query(true));
            assert(nodes_rc.size() == nodes.size());
        }

        std::shared_ptr<ISeeder<node_index>> seeder = build_seeder(
            this_query, is_reverse_complement, std::move(nodes)
        );

        std::shared_ptr<ISeeder<node_index>> seeder_rc;

        if (get_config().forward_and_reverse_complement
                || get_graph().get_mode() == DeBruijnGraph::CANONICAL)
            seeder_rc = build_seeder(reverse, !is_reverse_complement, std::move(nodes_rc));

        auto extender = build_extender(this_query, aligner_core.get_aggregator());

        if (get_graph().get_mode() == DeBruijnGraph::CANONICAL) {
            auto extender_rc = build_extender(reverse, aligner_core.get_aggregator());

            auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                     const auto &callback) {
                callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
            };

            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            aligner_core.align_both_directions(*seeder, *seeder_rc,
                                               *extender, *extender_rc,
                                               build_rev_comp_alignment_core);

        } else if (get_config().forward_and_reverse_complement) {
            auto extender_rc = build_extender(reverse, aligner_core.get_aggregator());
            aligner_core.align_best_direction(*seeder, *seeder_rc, *extender, *extender_rc);

        } else {
            aligner_core.align_one_direction(is_reverse_complement, *seeder, *extender);
        }

        aligner_core.flush();

        callback(header, std::move(paths));
    });
}

inline Vector<uint64_t> SeedFilter::labels_to_keep(const DBGAlignment &seed) {
    size_t found_count = 0;
    std::pair<size_t, size_t> idx_range {
        seed.get_clipping(), seed.get_clipping() + k_ - seed.get_offset()
    };
    for (node_index node : seed) {
        auto emplace = visited_nodes_.emplace(node, idx_range);
        auto &range = emplace.first.value();
        if (emplace.second) {
        } else if (range.first > idx_range.first || range.second < idx_range.second) {
            DEBUG_LOG("Node: {}; Prev_range: [{},{})", node, range.first, range.second);
            range.first = std::min(range.first, idx_range.first);
            range.second = std::max(range.second, idx_range.second);
            DEBUG_LOG("Node: {}; cur_range: [{},{})", node, range.first, range.second);
        } else {
            ++found_count;
        }

        if (idx_range.second - idx_range.first == k_)
            ++idx_range.first;

        ++idx_range.second;
    }

    if (found_count == seed.size())
        return {};

    return { std::numeric_limits<uint64_t>::max() };
}

inline void SeedFilter::update_seed_filter(const LabeledNodeRangeGenerator &generator) {
    generator([&](node_index node, uint64_t, size_t begin, size_t end) {
        auto emplace = visited_nodes_.emplace(node, std::make_pair(begin, end));
        auto &range = emplace.first.value();
        if (!emplace.second) {
            range.first = std::min(range.first, begin);
            range.second = std::max(range.second, end);
        }
    });
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_core(std::string_view query,
             const ISeeder<node_index> &seeder,
             IExtender<node_index> &extender,
             const LocalAlignmentCallback &callback,
             const MinScoreComputer &get_min_path_score) {
    bool filter_seeds = dynamic_cast<const ExactSeeder<node_index>*>(&seeder);
    constexpr uint64_t nlabel = std::numeric_limits<uint64_t>::max();

    std::vector<DBGAlignment> seeds = seeder.get_seeds();

    if (filter_seeds
            && std::any_of(seeds.begin(), seeds.end(),
                           [&](const auto &a) { return a.target_columns.size(); })) {
        constexpr auto get_end_ptr = [](const DBGAlignment &a) {
            return a.get_query().data() + a.get_query().size();
        };

        std::sort(seeds.begin(), seeds.end(), [&](const auto &a, const auto &b) {
            return get_end_ptr(a) < get_end_ptr(b);
        });

        for (size_t i = 0; i + 1 < seeds.size(); ++i) {
            DBGAlignment &seed = seeds[i];
            if (get_end_ptr(seeds[i + 1]) - get_end_ptr(seed) == 1
                    && seed.target_columns.size()
                    && seeds[i + 1].target_columns.size()) {
                uint64_t cnt = utils::count_intersection(
                    seed.target_columns.begin(), seed.target_columns.end(),
                    seeds[i + 1].target_columns.begin(), seeds[i + 1].target_columns.end()
                );

                if (!cnt || (cnt == seed.target_columns.size()
                                && graph_.get_mode() == DeBruijnGraph::CANONICAL)) {
                    if (seed.get_score() >= get_min_path_score(seed)) {
                        seed.extend_query_end(query.data() + query.size());
                        seed.trim_offset();
                        assert(seed.is_valid(graph_, &config_));
                        DEBUG_LOG("Alignment (seed): {}", seed);
                        callback(std::move(seed));
                    }

                    DEBUG_LOG("Skipping seed: {}", seed);
                    seed = DBGAlignment();
                }
            }
        }
    }

    std::sort(seeds.begin(), seeds.end(), LocalAlignmentGreater());

    for (size_t i = 0; i < seeds.size(); ++i) {
        DBGAlignment &seed = seeds[i];
        if (seed.empty())
            continue;

        score_t min_path_score = get_min_path_score(seed);

        // check if this seed has been explored before in an alignment and discard
        // it if so
        if (filter_seeds) {
            seed.target_columns = seed_filter_->labels_to_keep(seed);
            if (seed.target_columns.empty()) {
                DEBUG_LOG("Skipping seed: {}", seed);
                continue;
            }

            if (seed.target_columns.size() == 1 && seed.target_columns[0] == nlabel)
                seed.target_columns.clear();
        }

        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

        extender.initialize(seed);

        auto extensions = extender.get_extensions(min_path_score);

        // if the ManualSeeder is not used, then add nodes to the visited_nodes_
        // table to allow for seed filtration
        if (filter_seeds) {
            tsl::hopscotch_set<uint64_t> targets(seed.target_columns.begin(),
                                                 seed.target_columns.end());
            targets.insert(nlabel);
            if (seed.target_columns.empty()) {
                for (const DBGAlignment &extension : extensions) {
                    targets.insert(extension.target_columns.begin(),
                                   extension.target_columns.end());
                }
            }

            seed_filter_->update_seed_filter([&](const auto &callback) {
                extender.call_visited_nodes([&](node_index node, size_t begin, size_t end) {
                    for (uint64_t target : targets) {
                        callback(node, target, begin, end);
                    }
                });
            });
        }

        if (extensions.empty() && seed.get_score() >= min_path_score) {
            seed.extend_query_end(query.data() + query.size());
            seed.trim_offset();
            assert(seed.is_valid(graph_, &config_));
            DEBUG_LOG("Alignment (seed): {}", seed);
            callback(std::move(seed));
        }

        for (auto&& extension : extensions) {
            assert(extension.is_valid(graph_, &config_));
            callback(std::move(extension));
        }
    }
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_one_direction(bool orientation_to_align,
                      const ISeeder<node_index> &seeder,
                      IExtender<node_index> &extender) {
    std::string_view query = paths_.get_query(orientation_to_align);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        align_core(query, seeder, extender, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_best_direction(const ISeeder<node_index> &seeder,
                       const ISeeder<node_index> &seeder_rc,
                       IExtender<node_index> &extender,
                       IExtender<node_index> &extender_rc) {
    std::string_view forward = paths_.get_query();
    std::string_view reverse = paths_.get_query(true);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        align_core(forward, seeder, extender, alignment_callback, get_min_path_score);
        align_core(reverse, seeder_rc, extender_rc, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_both_directions(const ISeeder<node_index> &forward_seeder,
                        const ISeeder<node_index> &reverse_seeder,
                        IExtender<node_index> &forward_extender,
                        IExtender<node_index> &reverse_extender,
                        const SeederGenerator &rev_comp_core_generator) {
    std::string_view forward = paths_.get_query();
    std::string_view reverse = paths_.get_query(true);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        auto get_forward_alignments = [&](std::string_view query,
                                          std::string_view query_rc,
                                          const ISeeder<node_index> &seeder,
                                          IExtender<node_index> &extender) {
            std::vector<DBGAlignment> rc_of_alignments;

            DEBUG_LOG("Extending in forwards direction");
            align_core(query, seeder, extender, [&](DBGAlignment&& path) {
                score_t min_path_score = get_min_path_score(path);

                if (path.get_score() >= min_path_score)
                    alignment_callback(DBGAlignment(path));

                if (!path.get_clipping() || path.get_offset())
                    return;

                auto rev = path;
                rev.reverse_complement(graph_, query_rc);
                if (rev.empty()) {
                    DEBUG_LOG("This local alignment cannot be reversed, skipping");
                    return;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(rev.get_end_clipping());
                rev.trim_end_clipping();
                assert(rev.is_valid(graph_, &config_));

                // Pass the reverse complement of the forward alignment
                // as a seed for extension
                rc_of_alignments.emplace_back(std::move(rev));
            }, [&](const auto&) { return config_.min_cell_score; });

            return rc_of_alignments;
        };

        std::vector<DBGAlignment> rc_of_reverse = get_forward_alignments(
            reverse, forward, reverse_seeder, reverse_extender
        );
        std::vector<DBGAlignment> rc_of_forward = get_forward_alignments(
            forward, reverse, forward_seeder, forward_extender
        );

        DEBUG_LOG("Extending in reverse direction");
        rev_comp_core_generator(std::move(rc_of_reverse), [&](const auto &seeder) {
            align_core(forward, seeder, forward_extender,
                       alignment_callback, get_min_path_score);
        });

        rev_comp_core_generator(std::move(rc_of_forward), [&](const auto &seeder) {
            align_core(reverse, seeder, reverse_extender,
                       alignment_callback, get_min_path_score);
        });
    });
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_aggregate(const AlignmentGenerator &alignment_generator) {
    alignment_generator(
        [&](DBGAlignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            aggregator_.add_alignment(std::move(alignment));
        },
        [&](const DBGAlignment &seed) { return aggregator_.get_min_path_score(seed); }
    );
}

template class ISeedAndExtendAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
