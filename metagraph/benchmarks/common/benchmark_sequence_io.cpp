#include <cmath>
#include <random>

#include <benchmark/benchmark.h>

#include "common/threads/threading.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "../benchmark_graph_helpers.hpp"

const std::string file_prefix = "/tmp/bm_mg_outfile.fasta.gz";


template <int num_seqs>
static void BM_WriteRandomSequences(benchmark::State& state) {
    const std::string alphabet = "ATGCN";
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist4(0, 4);
    std::uniform_int_distribution<std::mt19937::result_type> dist1000(10, 1000);

    std::vector<std::string> sequences;
    sequences.reserve(num_seqs);
    for (size_t i = 0; i < num_seqs; ++i) {
        sequences.emplace_back(dist1000(rng), 'A');
        std::for_each(sequences.back().begin(),
                      sequences.back().end(),
                      [&](char &c) { c = alphabet[dist4(rng)]; });
    }

    set_num_threads(state.range(0));
    for (auto _ : state) {
        FastaWriter writer(file_prefix, "", false, state.range(0) - 1);
        for (const auto &sequence : sequences) {
            writer.write(sequence);
        }
    }
    set_num_threads(1);
}

BENCHMARK_TEMPLATE(BM_WriteRandomSequences, 10000)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);

BENCHMARK_TEMPLATE(BM_WriteRandomSequences, 100000)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);


static void BM_WriteContigs(benchmark::State& state) {
    std::string graph_file = "../tests/data/transcripts_1000.fa";
    auto graph = mg::bm::build_graph(graph_file);
    std::vector<std::string> contigs;
    graph->call_sequences([&](const auto &contig, auto&&) { contigs.push_back(contig); });

    set_num_threads(state.range(0));
    for (auto _ : state) {
        FastaWriter writer(file_prefix, "", false, state.range(0) - 1);
        for (const auto &contig : contigs) {
            writer.write(contig);
        }
    }
    set_num_threads(1);
}

BENCHMARK(BM_WriteContigs)
    ->Unit(benchmark::kMillisecond)
    ->DenseRange(1, 2, 1);
