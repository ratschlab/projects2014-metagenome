#ifndef __TEST_ANNOTATED_DBG_HELPERS_HPP__
#define __TEST_ANNOTATED_DBG_HELPERS_HPP__

#include "gtest/gtest.h"

#include <memory>
#include <vector>
#include <string>

#define protected public
#define private public
#include "annotated_dbg.hpp"
#include "annotate_column_compressed.hpp"
#include "static_annotators_def.hpp"

#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "masked_graph.hpp"

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels);

MaskedDeBruijnGraph build_masked_graph(const AnnotatedDBG &anno_graph,
                                       const std::vector<std::string> &ingroup,
                                       const std::vector<std::string> &outgroup,
                                       double outlabel_mixture = 0.0,
                                       double lazy_evaluation_density_cutoff = 0.05);

MaskedDeBruijnGraph build_masked_graph_lazy(const AnnotatedDBG &anno_graph,
                                            const std::vector<std::string> &ingroup,
                                            const std::vector<std::string> &outgroup,
                                            double outlabel_mixture = 0.0,
                                            double lazy_evaluation_density_cutoff = 0.05);

typedef ::testing::Types<std::pair<DBGBitmap, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashString, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annotate::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annotate::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashString, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annotate::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annotate::RowFlatAnnotator>
                        > GraphAnnotationPairTypes;

typedef ::testing::Types<std::pair<DBGBitmap, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, annotate::ColumnCompressed<>>,
                         std::pair<DBGBitmap, annotate::RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, annotate::RowFlatAnnotator>
                        > GraphNoNAnnotationPairTypes;

typedef ::testing::Types<std::pair<DBGHashString, annotate::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annotate::ColumnCompressed<>>,
                         std::pair<DBGHashString, annotate::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annotate::RowFlatAnnotator>
                        > GraphWithNAnnotationPairTypes;

#endif // __TEST_ANNOTATED_DBG_HELPERS_HPP__
