#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <sdsl/enc_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/logger.hpp"
#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {
namespace binmat {

const size_t RD_PATH_RESERVE_SIZE = 2;

/**
 * Sparsified representation of the underlying #BinaryMatrix that stores diffs between
 * successive nodes, rather than the full annotation.
 * The successor of a node (that is the node to diff against) is determined by a path in
 * an external graph structure, a #graph::DBGSuccinct, which has the property that rows
 * on a path are likely identical or very similar.
 *
 * RowDiff sparsification can be applied to any BinaryMatrix instance.
 * The row-diff binary matrix is defined by three data structures:
 *   1. #diffs_ the underlying sparsified (diffed) #BinaryMatrix
 *   2. #terminal_ rows marked as terminal are stored in full
 *   3. #graph_ the graph that was used to determine adjacent rows for sparsification
 * Retrieving data from RowDiff requires the associated #graph_. In order to get the
 * annotation for  i-th row, we start traversing the node corresponding to i in #graph_
 * and accumulate the values in #diffs until we hit a terminal node, which is stored in
 * full.
 */
//NOTE: Clang aggressively abuses the clause in the C++ standard (14.7.1/11) that allows
// virtual methods in template classes to not be instantiated if unused and mistakenly
// does not instantiate the virtual methods in this class, so I had to move definitions
// to the header (gcc works fine)
template <class BaseMatrix>
class RowDiff : public BinaryMatrix {
  public:
    typedef bit_vector_small anchor_bv_type;

    RowDiff() {}

    RowDiff(const graph::DBGSuccinct *graph,
               BaseMatrix &&diffs,
               const std::string &anchors_filename,
               bool load_anchors = true)
        : graph_(graph), diffs_(std::move(diffs)), anchors_filename_(anchors_filename) {
        if (load_anchors)
            load_terminal(anchors_filename_, &terminal_);
    }

    uint64_t num_columns() const override { return diffs_.num_columns(); }

    const std::string& anchors_filename() const { return anchors_filename_; }

    const graph::DBGSuccinct* graph() const { return graph_; }

    /**
     * Returns the number of set bits in the matrix.
     */
    uint64_t num_relations() const override { return diffs_.num_relations(); }

    uint64_t num_rows() const override { return diffs_.num_rows(); }
    void set_graph(const graph::DBGSuccinct *graph) { graph_ = graph; }

    bool get(Row row, Column column) const override;

    /**
     * Returns the given column.
     */
    std::vector<Row> get_column(Column column) const override;

    SetBitPositions get_row(Row row) const override;

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &row_ids) const override;

    bool load(std::istream &f) override;
    void serialize(std::ostream &f) const override ;

    void serialize(const std::string &filename) const;
    bool load(const std::string &filename);

    const anchor_bv_type& terminal() const { return terminal_; }

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    Vector<uint64_t> get_diff(uint64_t node_id) const { return diffs_.get_row(node_id); }

  private:
    static void load_terminal(const std::string &filename, anchor_bv_type *terminal);

    static void merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2);

    const graph::DBGSuccinct *graph_ = nullptr;

    BaseMatrix diffs_;
    anchor_bv_type terminal_;

    std::string anchors_filename_;
};

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

template <class BaseMatrix>
BinaryMatrix::SetBitPositions RowDiff<BaseMatrix>::get_row(Row row) const {
    Vector<uint64_t> result = get_diff(row);

    uint64_t boss_edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(row));
    const graph::boss::BOSS &boss = graph_->get_boss();

    while (!terminal_[row]) {
        graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
        assert(boss_edge > 1 && w != 0);

        // fwd always selects the last outgoing edge for a given node
        boss_edge = boss.fwd(boss_edge, w % boss.alph_size);
        row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                graph_->boss_to_kmer_index(boss_edge));
        merge(&result, get_diff(row));
    };
    return result;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
RowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    rd_ids.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());

    const graph::boss::BOSS &boss = graph_->get_boss();

    for (size_t i = 0; i < row_ids.size(); ++i) {
        Row row = row_ids[i];

        graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(
                graph::AnnotatedSequenceGraph::anno_to_graph_index(row));

        while (true) {
            row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                    graph_->boss_to_kmer_index(boss_edge));

            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            rd_paths_trunc[i].push_back(it.value());
            // If the node had been reached before, interrupt the diff path here.
            if (!is_new)
                break;

            rd_ids.push_back(row);

            if (terminal()[row])
                break;

            graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
            assert(boss_edge > 1 && w != 0);
            // fwd always selects the last outgoing edge for a given node
            boss_edge = boss.fwd(boss_edge, w % boss.alph_size);
        }
    }

    node_to_rd = VectorMap<Row, size_t>();

    std::vector<SetBitPositions> rd_rows = diffs_.get_rows(rd_ids);

    rd_ids = std::vector<Row>();

    // reconstruct annotation rows from row-diff
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        SetBitPositions &result = rows[i];
        // propagate back to reconstruct full annotations for predecessors
        for (auto it = rd_paths_trunc[i].rbegin(); it != rd_paths_trunc[i].rend(); ++it) {
            merge(&result, rd_rows[*it]);
            rd_rows[*it] = result;
        }
    }

    return rows;
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(std::istream &f) {
    uint64_t len;
    f.read(reinterpret_cast<char *>(&len), sizeof(uint64_t));
    anchors_filename_ = std::string(len, '\0');
    f.read(anchors_filename_.data(), len);
    common::logger->trace("Loading terminal nodes from {}", anchors_filename_);
    std::ifstream fterm(anchors_filename_, ios::binary);
    terminal_.load(fterm);
    fterm.close();

    return diffs_.load(f);
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(std::ostream &f) const {
    uint64_t len = anchors_filename_.size();
    f.write(reinterpret_cast<char *>(&len), sizeof(uint64_t));
    f.write(anchors_filename_.c_str(), len);
    diffs_.serialize(f);
};

template <class BaseMatrix>
void RowDiff<BaseMatrix>::load_terminal(const std::string &filename,
                                        anchor_bv_type *terminal) {
    std::ifstream f(filename, ios::binary);
    if (!f.good()) {
        common::logger->error("Could not open anchor file {}", filename);
        std::exit(1);
    }
    terminal->load(f);
    f.close();
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::merge(Vector<uint64_t> *result, const Vector<uint64_t> &diff2) {
    assert(std::is_sorted(result->begin(), result->end())
                   && std::is_sorted(diff2.begin(), diff2.end()));
    if (diff2.empty()) {
        return;
    }
    Vector<uint64_t> diff1;
    std::swap(*result, diff1);
    result->reserve(std::max(diff1.size(), diff2.size()));
    uint64_t idx1 = 0;
    uint64_t idx2 = 0;
    while (idx1 < diff1.size() && idx2 < diff2.size()) {
        if (diff1[idx1] == diff2[idx2]) {
            idx1++;
            idx2++;
        } else if (diff1[idx1] < diff2[idx2]) {
            result->push_back(diff1[idx1++]);
        } else {
            result->push_back(diff2[idx2++]);
        }
    }
    while (idx1 < diff1.size()) {
        result->push_back(diff1[idx1++]);
    }
    while (idx2 < diff2.size()) {
        result->push_back(diff2[idx2++]);
    }
}

} // namespace binmat
} // namespace annot
} // namespace mtg
