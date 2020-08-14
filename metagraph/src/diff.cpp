#include <cmath>
#include <random>

#include <gflags/gflags.h>
#include <progress_bar.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/select_support_mcl.hpp>

#include "annotation//binary_matrix/row_diff/row_diff.hpp"
#include "annotation//binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/annotation_converters.hpp"
#include "annotation//binary_matrix/row_flat/flat_matrix.hpp"
#include "cli/load/load_annotation.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

DEFINE_string(action, "store", "Action to perform [load, store]");
DEFINE_string(format, "diff", "Format to compress in: diff or zip");
DEFINE_string(o, "/tmp/graph_750_." + FLAGS_format, "Output file");
DEFINE_string(graph, "/tmp/diffs/graph_subset_750_primary.dbg", "Input graph file");
DEFINE_string(annotation,
              "/tmp/graph_subset_750_rowflat.flat.annodbg",
              "Input annotation file");
DEFINE_int32(threads, 4, "Numnber of threads to use");
DEFINE_string(kmer, "ACTA", "Kmer to read the annotation for");

// Note: if testing with many chunks use 'ulimit -n <max_files>' to increase the maxium
// number of files the system allows you to open. On mac the default is only 256!

using namespace mtg;

std::shared_ptr<graph::DBGSuccinct> load_graph(const std::string &name) {
    std::shared_ptr<graph::DBGSuccinct> result = std::make_shared<graph::DBGSuccinct>(20);
    if (!result->load(name)) {
        std::cerr << "Cannot load graph " << name << std::endl;
        std::exit(1);
    }
    std::cout << "Loaded graph with k-mer length " << result->get_k() << ", type "
              << result->get_state() << " number of nodes: " << result->num_nodes()
              << std::endl;
    return result;
}

std::unique_ptr<annot::MultiLabelEncoded<std::string>>
load_annotation(const std::string &anno_file) {
    std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation
            = cli::initialize_annotation(anno_file);

    if (!annotation->load(anno_file)) {
        std::cerr << "Can't load annotation from " << anno_file << std::endl;
        exit(1);
    }
    std::cout << "Loaded annotation " << anno_file << std::endl;
    return annotation;
}

void store_diff(const graph::DBGSuccinct &graph,
                annot::RowCompressed<std::string> &&annotation) __attribute__((optnone)) {
    //_Pragma("clang optimize off")
    uint64_t nnodes = graph.num_nodes();
    std::vector<std::vector<uint64_t>> tdiffs(nnodes + 1);
    sdsl::bit_vector terminal(nnodes + 1);

    std::unique_ptr<annot::RowDiffAnnotator> annotator
            = annot::convert_to_row_diff(graph, std::move(annotation), 4);

    annotator->serialize(FLAGS_o);

    //_Pragma("clang optimize on")
}

std::string to_string(const Vector<uint64_t> &v) {
    std::string result;
    for (const auto el : v) {
        result += (std::to_string(el) + " ");
    }
    return result;
}

void store_zip(const graph::AnnotatedDBG &annotated_dbg) {
    std::ofstream out(FLAGS_o);
    const graph::DeBruijnGraph &graph = annotated_dbg.get_graph();
    ProgressBar progress_bar(graph.num_nodes(), "Building Diff-anno");
    graph.call_sequences([&](const std::string &, const std::vector<uint64_t> &path) {
        std::vector<uint64_t> anno_ids;

        for (const uint64_t node_id : path) {
            anno_ids.push_back(graph::AnnotatedSequenceGraph::graph_to_anno_index(node_id));
        }
        std::vector<Vector<uint64_t>> rows
                = annotated_dbg.get_annotation().get_matrix().get_rows(anno_ids);

        for (uint32_t i = 0; i < rows.size(); ++i) {
            sdsl::bit_vector row(annotated_dbg.get_annotation().get_matrix().num_columns());
            for (const auto v : rows[i]) {
                row[v] = 1;
            }
            row.serialize(out);
        }

        progress_bar += path.size();
    });

    out.close();
}

int main(int argc, char* argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::trace);

    std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation
            = load_annotation(FLAGS_annotation);
    std::shared_ptr<graph::DBGSuccinct> graph = load_graph(FLAGS_graph);

    if (FLAGS_action == "store") {
        if (FLAGS_format == "diff") {
            store_diff(*graph, *annotation);
        } else if (FLAGS_format == "zip") {
            graph::AnnotatedDBG annotated_dbg(graph, std::move(annotation), false);
            store_zip(annotated_dbg);
        } else {
            std::cerr << "Invalid -format, only 'zip' and 'diff' supported.";
            std::exit(1);
        }
    } else if (FLAGS_action == "load") {
        annot::binmat::RowDiff anno;
        anno.load(FLAGS_annotation);
        uint64_t node_id = graph->kmer_to_node(FLAGS_kmer);
        if (node_id == 0) {
            std::cout << "Kmer " << FLAGS_kmer << " doesn't exist.\n";
            std::exit(1);
        }
        std::cout << "Annotation for k-mer " << FLAGS_kmer << " is "
                  << to_string(anno.get_row(node_id));
    }
}
