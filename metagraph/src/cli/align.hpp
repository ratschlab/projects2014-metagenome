#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

class Config;
class DeBruijnGraph;
class AnnotatedDBG;
class IDBGAligner;

std::unique_ptr<IDBGAligner>
build_aligner(const DeBruijnGraph &graph, const Config &config);

std::unique_ptr<IDBGAligner>
build_masked_aligner(const AnnotatedDBG &anno_graph, const Config &config);

int align_to_graph(Config *config);

#endif // __ALIGN_GRAPH_HPP__
