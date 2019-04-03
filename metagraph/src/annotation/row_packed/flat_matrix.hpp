#ifndef __FLAT_MATRIX_HPP__
#define __FLAT_MATRIX_HPP__

#include <vector>

#include "binary_matrix.hpp"
#include "bit_vector.hpp"


// TODO: fix this and make consistent with other annotators
// (Empty annotatino can have multiple objects)
template <typename BitVector = bit_vector_sd>
class RowConcatenated : public BinaryMatrix {
  public:
    using RowCallback = std::function<void(const std::vector<Column> &)>;
    using ValueCallback = std::function<void(Row, Column)>;

    RowConcatenated() {}
    RowConcatenated(const std::function<void(RowCallback)> &call_rows,
                    uint64_t num_columns,
                    uint64_t num_rows,
                    uint64_t num_set_bits);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return compressed_rows_->size() / num_columns_; }

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    uint64_t select1(uint64_t i) {
        return compressed_rows_->select1(i);
    }
    void dump() {
        for(uint64_t i = 1; i <= compressed_rows_->num_set_bits(); i++) {
            uint64_t idx = compressed_rows_->select1(i);
            std::cout << idx << std::endl;
            assert(get(idx/num_columns(), idx%num_columns()));
        }
    }

    // number of ones in the matrix
    uint64_t num_relations() const { return compressed_rows_->num_set_bits(); }

    const BitVector& data() const { return *compressed_rows_; }

  private:
    std::unique_ptr<BitVector> compressed_rows_;
    uint64_t num_columns_ = 0;
};

#endif // __FLAT_MATRIX_HPP__
