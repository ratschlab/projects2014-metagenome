#include "row_diff.hpp"

#include <fstream>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(const std::string &filename) const {
    std::ofstream f(filename, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(const std::string &filename) {
    std::ifstream f(filename, ios::binary);
    bool result = load(f);
    f.close();
    return result;
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::load_anchor(const std::string& filename) {
    if (!std::filesystem::exists(filename)) {
        common::logger->error("Can't read anchor file: {}", filename);
        std::exit(1);
    }
    std::ifstream f(filename, ios::binary);
    if (!f.good()) {
        common::logger->error("Could not open anchor file {}", filename);
        std::exit(1);
    }
    anchor_.load(f);
    f.close();
}

template
class RowDiff<ColumnMajor>;
template
class RowDiff<BRWT>;

} // namespace binmat
} // namespace annot
} // namespace mtg
