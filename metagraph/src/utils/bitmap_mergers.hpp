#ifndef __BIT_VECTOR_UTILS_HPP__
#define __BIT_VECTOR_UTILS_HPP__

#include <vector>
#include <string>
#include <memory>
#include <functional>

#include "bit_vector.hpp"
#include "serialization.hpp"


namespace utils {

// Read indices of set bits from a vector of VectorStreams
class RowsFromColumnsTransformer {
  public:
    explicit RowsFromColumnsTransformer(const std::vector<const bit_vector*> &columns);
    explicit RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector>> &columns);
    explicit RowsFromColumnsTransformer(const std::vector<std::shared_ptr<bit_vector>> &columns);

    template <typename BitVector>
    explicit RowsFromColumnsTransformer(const std::vector<BitVector> &columns) {
        std::vector<const bit_vector*> cols;
        for (const auto &col : columns) {
            cols.push_back(&col);
        }
        initialize(cols);
    }

    RowsFromColumnsTransformer(const bit_vector &columns_concatenated,
                               uint64_t column_size);

    RowsFromColumnsTransformer(RowsFromColumnsTransformer&&) = default;
    RowsFromColumnsTransformer& operator=(RowsFromColumnsTransformer&&) = default;

    using ValueCallback = std::function<void(uint64_t /*row*/,
                                             uint64_t /*column*/)>;
    void call_next(ValueCallback callback);

    uint64_t rows() const { return num_rows_; }
    uint64_t columns() const { return streams_.size(); }

    // get the number of set bits in the vectors left
    uint64_t values_left() const { return num_set_bits_left_; }

    template <typename Vector = std::vector<uint64_t>>
    void call_rows(const std::function<void(const Vector &)> &callback);

  private:
    void initialize(const std::vector<const bit_vector*> &columns);
    void init_heap();

    std::vector<std::unique_ptr<VectorInStream>> streams_;
    uint64_t num_set_bits_left_ = 0;
    uint64_t num_rows_;

    // store pair of kmer index and label index
    struct kmer_label_pair {
        uint64_t row_id;
        uint64_t col_id;

        bool operator<(const kmer_label_pair &other) const {
            return row_id > other.row_id || (row_id == other.row_id
                                                && col_id > other.col_id);
        }
    };

    std::priority_queue<kmer_label_pair,
                        std::vector<kmer_label_pair>> index_heap_;
};

class RowsFromColumnsIterator {
  public:
    RowsFromColumnsIterator(std::unique_ptr<RowsFromColumnsTransformer> transformer) {
        transformer_ = std::move(transformer);
    }

    std::tuple<uint64_t, uint64_t> next_set_bit();
    uint64_t values_left() { return transformer_->values_left(); }

    template <typename Vector = std::vector<uint64_t>>
    Vector next_row();

  private:
    std::unique_ptr<RowsFromColumnsTransformer> transformer_;

    uint64_t i_ = 0;
    uint64_t row_ = 0;
    uint64_t column_;
};

template <class BitVectorType = bit_vector_stat>
std::vector<std::unique_ptr<bit_vector>>
transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix);

} // namespace utils

#endif // __BIT_VECTOR_UTILS_HPP__
