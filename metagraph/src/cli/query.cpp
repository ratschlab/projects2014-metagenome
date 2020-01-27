#include "query.hpp"

#include <ips4o.hpp>
#include <fmt/format.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/int_vector_algorithm.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "align.hpp"

static const size_t kRowBatchSize = 100'000;
const bool kPrefilterWithBloom = true;

using mg::common::logger;


void execute_query(const std::string &seq_name,
                   const std::string &sequence,
                   bool count_labels,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream) {
    std::string output;
    output.reserve(1'000);

    if (count_labels) {
        auto top_labels = anno_graph.get_top_labels(sequence,
                                                    num_top_labels,
                                                    discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        output += seq_name;

        for (const auto &[label, count] : top_labels) {
            output += "\t<";
            output += label;
            output += ">:";
            output += fmt::format_int(count).c_str();
        }

        output += '\n';

    } else {
        auto labels_discovered = anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return;

        output += seq_name;
        output += '\t';
        output += utils::join_strings(labels_discovered, anno_labels_delimiter);
        output += '\n';
    }

    output_stream << output;
}

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads) {
    const auto *full_dbg = &anno_graph.get_graph();
    if (!full_dbg)
        throw std::runtime_error("Error: batch queries are supported only for de Bruijn graphs");

    const auto &full_annotation = anno_graph.get_annotation();

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(full_dbg);
    if (kPrefilterWithBloom && dbg_succ) {
        if (dbg_succ->get_bloom_filter())
            logger->trace("[Query graph construction] Started indexing k-mers pre-filtered with Bloom filter");

        call_sequences([&graph,&dbg_succ](const std::string &sequence) {
            graph->add_sequence(sequence, get_missing_kmer_skipper(
                dbg_succ->get_bloom_filter(),
                sequence
            ));
        });
    } else {
        call_sequences([&graph](const std::string &sequence) {
            graph->add_sequence(sequence);
        });
    }

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    graph->call_sequences(
        [&](const std::string &contig, const auto &path) { contigs.emplace_back(contig, path); },
        full_dbg->is_canonical_mode()
    );

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    if (full_dbg->is_canonical_mode()) {
        // construct graph storing all distinct k-mers in query
        graph = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), true);

        for (const auto &pair : contigs) {
            graph->add_sequence(pair.first);
        }

        logger->trace("[Query graph construction] k-mers reindexed in canonical mode in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    // map contigs onto the full graph
    auto index_in_full_graph
        = std::make_shared<std::vector<uint64_t>>(graph->max_index() + 1, 0);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
    for (size_t i = 0; i < contigs.size(); ++i) {

        auto contig = std::move(contigs[i].first);
        auto path = std::move(contigs[i].second);

        if (graph->is_canonical_mode()) {
            size_t j = 0;
            graph->map_to_nodes(contig, [&](auto node) { path[j++] = node; });
            assert(j == path.size());
        }

        size_t j = 0;

        full_dbg->map_to_nodes(contig,
            [&](auto node_in_full) { (*index_in_full_graph)[path[j++]] = node_in_full; }
        );

        assert(j == path.size());
    }

    logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
    timer.reset();

    contigs.clear();

    assert(!(*index_in_full_graph)[0]);

    if (discovery_fraction > 0) {
        sdsl::bit_vector mask(graph->max_index() + 1, false);

        call_sequences([&](const std::string &sequence) {
            if (sequence.length() < graph->get_k())
                return;

            const size_t num_kmers = sequence.length() - graph->get_k() + 1;
            const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
            const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
            size_t num_kmers_discovered = 0;
            size_t num_kmers_missing = 0;

            std::vector<DeBruijnGraph::node_index> nodes;
            nodes.reserve(num_kmers);

            graph->map_to_nodes(sequence,
                [&](auto node) {
                    if ((*index_in_full_graph)[node]) {
                        num_kmers_discovered++;
                        nodes.push_back(node);
                    } else {
                        num_kmers_missing++;
                    }
                },
                [&]() { return num_kmers_missing > max_kmers_missing
                                || num_kmers_discovered >= min_kmers_discovered; }
            );

            if (num_kmers_missing <= max_kmers_missing) {
                for (auto node : nodes) { mask[node] = true; }
            }
        });

        // correcting the mask
        call_zeros(mask, [&](auto i) { (*index_in_full_graph)[i] = 0; });

        logger->trace("[Query graph construction] Reduced k-mer dictionary in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    assert(index_in_full_graph.get());

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_query;
    from_full_to_query.reserve(index_in_full_graph->size());

    for (uint64_t node = 0; node < index_in_full_graph->size(); ++node) {
        if ((*index_in_full_graph)[node]) {
            from_full_to_query.emplace_back(
                AnnotatedDBG::graph_to_anno_index((*index_in_full_graph)[node]),
                AnnotatedDBG::graph_to_anno_index(node)
            );
        }
    }

    ips4o::parallel::sort(from_full_to_query.begin(), from_full_to_query.end(),
                          utils::LessFirst(), num_threads);

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = std::make_unique<annotate::RowCompressed<>>(
        graph->max_index(),
        full_annotation.get_label_encoder().get_labels(),
        [&](annotate::RowCompressed<>::CallRow call_row) {

            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (uint64_t batch_begin = 0;
                                batch_begin < from_full_to_query.size();
                                                batch_begin += kRowBatchSize) {

                const uint64_t batch_end
                    = std::min(batch_begin + kRowBatchSize,
                               static_cast<uint64_t>(from_full_to_query.size()));

                std::vector<uint64_t> row_indexes;
                row_indexes.reserve(batch_end - batch_begin);

                for (uint64_t i = batch_begin; i < batch_end; ++i) {
                    assert(from_full_to_query[i].first < full_annotation.num_objects());

                    row_indexes.push_back(from_full_to_query[i].first);
                }

                auto rows = full_annotation.get_label_codes(row_indexes);

                assert(rows.size() == batch_end - batch_begin);

                for (uint64_t i = batch_begin; i < batch_end; ++i) {
                    call_row(from_full_to_query[i].second,
                             std::move(rows[i - batch_begin]));
                }
            }
        }
    );

    logger->trace("[Query graph construction] Query annotation constructed in {} sec",
                  timer.elapsed());
    timer.reset();

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(graph,
        [=](auto i) -> bool { return (*index_in_full_graph)[i]; }
    );

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(masked_graph, std::move(annotation));
}


int query_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() == 1);

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1);

    Timer timer;

    std::unique_ptr<IDBGAligner> aligner;
    if (config->align_sequences && !config->fast)
        aligner.reset(build_aligner(*graph, *config).release());

    // iterate over input files
    for (const auto &file : files) {
        logger->trace("Parsing sequences from file '{}'", file);

        Timer curr_timer;

        size_t seq_count = 0;

        const auto *graph_to_query = anno_graph.get();

        // Graph constructed from a batch of queried sequences
        // Used only in fast mode
        std::unique_ptr<AnnotatedDBG> query_graph;
        std::vector<std::pair<std::string, std::string>> named_queries;

        if (aligner) {
            read_fasta_file_critical(
                file,
                [&](kseq_t *seq) {
                    std::string query(seq->seq.s);
                    auto alignments = aligner->align(query);
                    named_queries.emplace_back(
                        alignments.size() ? std::move(alignments.front().get_sequence())
                                          : std::move(query),
                        seq->name.s
                    );
                }
            );
        } else {
            read_fasta_file_critical(
                file,
                [&](kseq_t *seq) { named_queries.emplace_back(seq->seq.s, seq->name.s); },
                config->forward_and_reverse
            );
        }

        if (config->fast) {
            query_graph = construct_query_graph(*anno_graph,
                [&](auto call_sequence) {
                    for (const auto &[query, name] : named_queries) {
                        call_sequence(query);
                    }
                },
                config->count_labels ? 0 : config->discovery_fraction,
                get_num_threads()
            );

            graph_to_query = query_graph.get();

            logger->trace("Query graph constructed for '{}' in {} sec",
                          file, curr_timer.elapsed());
        }

        for (const auto &[query, name] : named_queries) {
            thread_pool.enqueue(execute_query,
                fmt::format_int(seq_count++).str() + "\t" + name,
                std::string(query),
                config->count_labels,
                config->print_signature,
                config->suppress_unlabeled,
                config->num_top_labels,
                config->discovery_fraction,
                config->anno_labels_delimiter,
                std::ref(*graph_to_query),
                std::ref(std::cout)
            );
        }

        // wait while all threads finish processing the current file
        thread_pool.join();

        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}
