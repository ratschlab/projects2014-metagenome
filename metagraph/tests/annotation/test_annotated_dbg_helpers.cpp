#include "test_annotated_dbg_helpers.hpp"

#include "../graph/all/test_dbg_helpers.hpp"

#include "annotation/annotation_converters.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph;
using namespace mtg::annot;

template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels,
                                               DeBruijnGraph::Mode mode) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences, mode);
    auto canonical = std::dynamic_pointer_cast<CanonicalDBG>(graph);

    uint64_t max_index = canonical ? canonical->get_graph().max_index()
                                   : graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<ColumnCompressed<>>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    if (!std::is_same<Annotation, ColumnCompressed<>>::value)
        anno_graph = std::make_unique<AnnotatedDBG>(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                convert<Annotation>(
                    std::move(dynamic_cast<ColumnCompressed<>&>(*anno_graph->annotator_)
                ))
            )
        );

    return anno_graph;
}

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode);


} // namespace test
} // namespace mtg
