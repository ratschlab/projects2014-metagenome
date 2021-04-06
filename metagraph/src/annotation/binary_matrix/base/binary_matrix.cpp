#include "binary_matrix.hpp"

#include "common/serialization.hpp"


namespace mtg {
namespace annot {
namespace binmat {

std::vector<sdsl::bit_vector>
BinaryMatrix::has_column(const std::vector<Row> &row_ids, const SetBitPositions &columns) const {
    std::vector<sdsl::bit_vector> row_masks;
    row_masks.reserve(columns.size());

    for (Column column : columns) {
        sdsl::bit_vector &row_mask = row_masks.emplace_back(row_ids.size(), false);
        for (size_t i = 0; i < row_ids.size(); ++i) {
            if (get(row_ids[i], column))
                row_mask[i] = true;
        }
    }

    return row_masks;
}

std::vector<BinaryMatrix::SetBitPositions>
BinaryMatrix::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows[i] = get_row(row_ids[i]);
    }

    return rows;
}

std::vector<BinaryMatrix::Column>
BinaryMatrix::slice_rows(const std::vector<Row> &row_ids) const {
    std::vector<BinaryMatrix::Column> slice;

    for (uint64_t i : row_ids) {
        for (uint64_t j : get_row(i)) {
            slice.push_back(j);
        }
        slice.push_back(std::numeric_limits<Column>::max());
    }

    return slice;
}

std::vector<std::pair<BinaryMatrix::Column, size_t /* count */>>
BinaryMatrix::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
                       size_t min_count,
                       size_t count_cap) const {
    assert(count_cap >= min_count);

    if (!count_cap)
        return {};

    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<size_t> col_counts(num_columns(), 0);
    size_t max_matched = 0;
    size_t total_checked = 0;

    for (auto [i, count] : index_counts) {
        if (max_matched + (total_sum_count - total_checked) < min_count)
            break;

        for (size_t j : get_row(i)) {
            assert(j < col_counts.size());

            col_counts[j] += count;
            max_matched = std::max(max_matched, col_counts[j]);
        }

        total_checked += count;
    }

    if (max_matched < min_count)
        return {};

    std::vector<std::pair<uint64_t, size_t>> result;
    result.reserve(col_counts.size());

    for (size_t j = 0; j < num_columns(); ++j) {
        if (col_counts[j] >= min_count) {
            result.emplace_back(j, std::min(col_counts[j], count_cap));
        }
    }

    return result;
}


template <typename RowType>
StreamRows<RowType>::StreamRows(const std::string &filename, size_t offset) {
    std::ifstream instream(filename, std::ios::binary);

    if (!instream.good() || !instream.seekg(offset).good())
        throw std::ifstream::failure("Cannot read rows from file " + filename);

    (void)load_number(instream);
    (void)load_number(instream);

    inbuf_ = sdsl::int_vector_buffer<>(filename,
                                       std::ios::in | std::ios::binary,
                                       1024 * 1024,
                                       0,
                                       false,
                                       instream.tellg());
}

template <typename RowType>
RowType* StreamRows<RowType>::next_row() {
    row_.clear();

    while (i_ < inbuf_.size()) {
        auto value = inbuf_[i_++];
        if (value - 1 > std::numeric_limits<typename RowType::value_type>::max())
            throw std::ifstream::failure("Integer overflow: trying to read too"
                                         " large column index: " + std::to_string(value - 1));
        if (value) {
            row_.push_back(value - 1);
        } else {
            return &row_;
        }
    }
    return nullptr;
}

template class StreamRows<BinaryMatrix::SetBitPositions>;


void append_row_major(const std::string &filename,
                      const std::function<void(BinaryMatrix::RowCallback)> &call_rows,
                      uint64_t num_cols) {
    std::ofstream outstream(filename, std::ios::binary | std::ios::app);

    uint64_t num_rows = 0;

    // write dummy num_rows value to fill in later
    const uint64_t header_offs = outstream.tellp();
    serialize_number(outstream, 0);
    serialize_number(outstream, num_cols);
    const uint64_t iv_offs = outstream.tellp();
    outstream.close();

    {
        auto outbuf = sdsl::int_vector_buffer<>(filename,
                                                std::ios::out | std::ios::binary,
                                                1024 * 1024,
                                                sdsl::bits::hi(num_cols) + 1,
                                                false,
                                                iv_offs);

        call_rows([&](const auto &row) {
            for (auto val : row) {
                outbuf.push_back(val + 1);
            }
            outbuf.push_back(0);
            num_rows++;
        });

        outbuf.close();
    }

    outstream.open(filename, std::ios::in | std::ios::out | std::ios::binary);
    outstream.seekp(header_offs);
    serialize_number(outstream, num_rows);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
