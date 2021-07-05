#include "dbg_aligner.hpp"

#include "common/algorithms.hpp"
#include "graph/representation/rc_dbg.hpp"

namespace mtg {
namespace graph {
namespace align {


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
                             Args&&... args)
          : graph_(graph), config_(config), rc_dbg_(graph_),
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

    // Align the forward and reverse complement of the query sequence in both
    // directions and return the overall best alignment. e.g., for the forward query
    // 1. Find all seeds of its reverse complement
    // 2. Given a seed, extend forwards to get alignment A
    // 3. Reverse complement the alignment to get A', treat it like a new seed
    // 4. Extend A' forwards to get the final alignment A''
    void align_both_directions(const ISeeder<node_index> &forward_seeder,
                               const ISeeder<node_index> &reverse_seeder,
                               IExtender<node_index> &forward_extender,
                               IExtender<node_index> &reverse_extender);

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
    RCDBG rc_dbg_;

    DBGQueryAlignment paths_;
    AlignmentAggregator<node_index, AlignmentCompare> aggregator_;
};

IDBGAligner::DBGQueryAlignment IDBGAligner::align(std::string_view query,
                                                  bool is_reverse_complement) const {
    DBGQueryAlignment result(query);
    align_batch({ Query{ std::string{}, query, is_reverse_complement} },
        [&](std::string_view, DBGQueryAlignment&& alignment) {
            result = std::move(alignment);
        }
    );

    return result;
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
              const AlignmentCallback &callback) const {
    for (const auto &[header, query, is_reverse_complement] : seq_batch) {
        SeedAndExtendAlignerCore<AlignmentCompare> core(
            graph_, config_, query, is_reverse_complement
        );
        auto &paths = core.get_paths();
        std::string_view this_query = paths.get_query(is_reverse_complement);
        assert(this_query == query);

        std::vector<node_index> nodes = map_sequence_to_nodes(graph_, query);

        auto seeder = build_seeder(this_query, is_reverse_complement, nodes);
        auto extender = build_extender(this_query, core.get_aggregator());

#if ! _PROTEIN_GRAPH
        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            std::vector<node_index> nodes_rc(nodes);
            std::string dummy(query);
            reverse_complement_seq_path(graph_, dummy, nodes_rc);
            assert(dummy == paths.get_query(!is_reverse_complement));
            assert(nodes_rc.size() == nodes.size());

            std::string_view reverse = paths.get_query(!is_reverse_complement);

            auto seeder_rc = build_seeder(reverse, !is_reverse_complement, nodes_rc);
            auto extender_rc = build_extender(reverse, core.get_aggregator());

            core.align_both_directions(*seeder, *seeder_rc, *extender, *extender_rc);

        } else {
            core.align_one_direction(is_reverse_complement, *seeder, *extender);
        }
#else
        core.align_one_direction(is_reverse_complement, *seeder, *extender);
#endif

        core.flush();

        callback(header, std::move(paths));
    };
}

template <class AlignmentCompare>
inline void SeedAndExtendAlignerCore<AlignmentCompare>
::align_core(std::string_view query,
             const ISeeder<node_index> &seeder,
             IExtender<node_index> &extender,
             const LocalAlignmentCallback &callback,
             const MinScoreComputer &get_min_path_score) {
    for (DBGAlignment &seed : seeder.get_seeds()) {
        if (seed.empty())
            continue;

        score_t min_path_score = get_min_path_score(seed);

        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

        auto extensions = extender.get_extensions(seed, min_path_score);

        if (extensions.empty() && seed.get_score() >= min_path_score) {
            seed.extend_query_end(query.data() + query.size());
            seed.trim_offset();
            DEBUG_LOG("Alignment (seed): {}", seed);
            callback(std::move(seed));
        }

        for (auto&& extension : extensions) {
            DEBUG_LOG("Alignment (extension): {}", extension);
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
::align_both_directions(const ISeeder<node_index> &forward_seeder,
                        const ISeeder<node_index> &reverse_seeder,
                        IExtender<node_index> &forward_extender,
                        IExtender<node_index> &reverse_extender) {
#if _PROTEIN_GRAPH
    throw std::runtime_error("Only alignment in one direction supported for Protein graphs");
#endif

    std::string_view forward = paths_.get_query();
    std::string_view reverse = paths_.get_query(true);

    bool use_rcdbg = graph_.get_mode() != DeBruijnGraph::CANONICAL && config_.forward_and_reverse_complement;
    const DeBruijnGraph &rc_graph = use_rcdbg ? rc_dbg_ : graph_;

    auto is_reversible = [this](const DBGAlignment &alignment) {
        return graph_.get_mode() == DeBruijnGraph::CANONICAL
            && alignment.get_orientation()
            && !alignment.get_offset();
    };

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        auto get_forward_alignments = [&](std::string_view query,
                                          std::string_view query_rc,
                                          const ISeeder<node_index> &seeder,
                                          IExtender<node_index> &extender) {
            size_t farthest_reach = 0;
            score_t max_score = config_.min_cell_score;
            std::vector<DBGAlignment> rc_of_alignments;

            DEBUG_LOG("Extending in forwards direction");
            align_core(query, seeder, extender,
                [&](DBGAlignment&& path) {
                    score_t min_path_score = get_min_path_score(path);

                    farthest_reach = std::max(farthest_reach, path.get_query().size() + path.get_clipping());
                    max_score = std::max(max_score, path.get_score());

                    if (path.get_score() >= min_path_score) {
                        if (is_reversible(path)) {
                            DBGAlignment out_path = path;
                            out_path.reverse_complement(graph_, query_rc);
                            assert(out_path.size());
                            alignment_callback(std::move(out_path));
                        } else {
                            alignment_callback(DBGAlignment(path));
                        }
                    }

                    if (!path.get_clipping() || path.get_offset())
                        return;

                    DBGAlignment rev = path;
                    rev.reverse_complement(rc_graph, query_rc);
                    assert(rev.is_valid(rc_graph, &config_));

                    if (rev.empty()) {
                        DEBUG_LOG("This local alignment cannot be reversed, skipping");
                        return;
                    }

                    // Remove any character skipping from the end so that the
                    // alignment can proceed
                    assert(rev.get_end_clipping());
                    rev.trim_end_clipping();

                    // Pass the reverse complement of the forward alignment
                    // as a seed for extension
                    rc_of_alignments.emplace_back(std::move(rev));
                },
                [&](const DBGAlignment &seed) {
                    return seed.get_clipping() <= farthest_reach
                        && config_.fraction_of_top > 0
                            ? max_score * config_.fraction_of_top
                            : config_.min_cell_score;
                }
            );

            std::sort(rc_of_alignments.begin(), rc_of_alignments.end(),
                      LocalAlignmentGreater());

            return rc_of_alignments;
        };

        ManualSeeder<node_index> rc_of_reverse(get_forward_alignments(
            reverse, forward, reverse_seeder, reverse_extender
        ));

        ManualSeeder<node_index> rc_of_forward(get_forward_alignments(
            forward, reverse, forward_seeder, forward_extender
        ));

        auto finish_alignment = [&](std::string_view query,
                                    std::string_view query_rc,
                                    const ManualSeeder<node_index> &seeder,
                                    IExtender<node_index> &extender) {
            if (use_rcdbg)
                extender.set_graph(rc_dbg_);

            align_core(query_rc, seeder, extender,
                [&](DBGAlignment&& path) {
                    if (use_rcdbg || is_reversible(path)) {
                        path.reverse_complement(rc_graph, query);
                        if (path.empty())
                            return;
                    }

                    assert(path.is_valid(graph_, &config_));
                    alignment_callback(std::move(path));
                },
                get_min_path_score
            );
        };

        if (rc_of_forward.data().size() && (rc_of_reverse.data().empty()
                || rc_of_forward.data()[0].get_score() >= rc_of_reverse.data()[0].get_score())) {
            finish_alignment(forward, reverse, rc_of_forward, reverse_extender);
            finish_alignment(reverse, forward, rc_of_reverse, forward_extender);
        } else {
            finish_alignment(reverse, forward, rc_of_reverse, forward_extender);
            finish_alignment(forward, reverse, rc_of_forward, reverse_extender);
        }
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
        [&](const DBGAlignment &seed) {
            return std::max(
                aggregator_.get_min_path_score(seed),
                static_cast<score_t>(aggregator_.get_max_path_score() * config_.fraction_of_top)
            );
        }
    );
}

template class ISeedAndExtendAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
