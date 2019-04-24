#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <cassert>
#include <vector>
#include <functional>
#include <filesystem>
#include <progress_bar.hpp>

#include "annotate.hpp"
#include "static_annotators_def.hpp"


namespace annotate {

typedef LabelEncoder<std::string> LEncoder;


template <typename Label>
class RowCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(RowCompressed<Label>&& annotation);


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(ColumnCompressed<Label>&& annotation);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_simple_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_threads = 1);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_greedy_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t num_threads = 1);

template <class StaticAnnotation>
void relax_BRWT(StaticAnnotation *annotation,
                size_t relax_max_arity,
                size_t num_threads = 1);

template <class StaticAnnotation, typename Label = std::string>
std::unique_ptr<StaticAnnotation>
convert(const std::string &filename) {
    using MatrixType = typename StaticAnnotation::binary_matrix_type;

    uint64_t num_rows;
    uint64_t num_relations;
    RowCompressed<Label>::stream_counts(filename, &num_rows, &num_relations);

    auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);

    constexpr size_t num_passes = std::is_same<MatrixType, Rainbowfish>::value ? 2u : 1u;
    ProgressBar progress_bar(num_rows*num_passes, "Processing rows");

    auto callback = [&](auto callback) {
        auto annotator = std::make_unique<typename RowCompressed<Label>::StreamRows>(filename);
        for (uint64_t r = 0; r < num_rows; ++r) {
            auto row = annotator->next_row();
            std::sort(row->begin(), row->end());
            callback(*row);
            ++progress_bar;
        }
    };

    std::unique_ptr<MatrixType> matrix;
    if constexpr (std::is_same<MatrixType, RowConcatenated<> >::value) {
        matrix = std::make_unique<MatrixType>(callback, label_encoder->size(), num_rows, num_relations);
    } else if constexpr (std::is_same<MatrixType, Rainbowfish>::value) {
        matrix = std::make_unique<MatrixType>(callback, label_encoder->size());
    } else if constexpr (std::is_same<MatrixType, BinRelWT>::value) {
        matrix = std::make_unique<MatrixType>(callback, num_relations, label_encoder->size());
    } else if constexpr (std::is_same<MatrixType, BinRelWT_sdsl>::value) {
        matrix = std::make_unique<MatrixType>(callback, num_relations, label_encoder->size());
    } else {
        static_assert(utils::dependent_false<StaticAnnotation>::value);
    }

    return std::make_unique<StaticAnnotation>(std::move(matrix), *label_encoder);
}

void merge_rows(const std::vector<const LEncoder*> &label_encoders,
                std::function<const std::vector<uint64_t>(uint64_t)> get_next_row,
                uint64_t num_rows,
                std::function<void(const LEncoder& /* merged encoder */,
                                   std::function<void(BinaryMatrix::RowCallback)> /* call row */)> callback);

template <class ToAnnotation, typename Label>
uint64_t merge(const std::vector<const MultiLabelEncoded<uint64_t, Label>*> &annotators,
               const std::vector<std::string> &filenames,
               const std::string &outfile) {
    assert((annotators.size() || filenames.size()) && "nothing to merge");

    uint64_t num_rows;
    uint64_t num_relations;

    if (annotators.size()) {
        num_rows = annotators.at(0)->num_objects();
    } else {
        RowCompressed<Label>::stream_counts(filenames.at(0), &num_rows, &num_relations);
    }
    assert(num_rows);

    std::vector<const LEncoder*> label_encoders;
    std::vector<std::unique_ptr<IterateRows>> annotator_row_iterators;
    for (auto annotator : annotators) {
        label_encoders.push_back(&annotator->label_encoder_);
        annotator_row_iterators.push_back(annotator->iterator());
    }

    std::vector<std::unique_ptr<const LEncoder> > loaded_label_encoders;
    std::vector<std::unique_ptr<class RowCompressed<Label>::StreamRows> > streams;
    for (auto filename : filenames) {
        if (utils::ends_with(filename, RowCompressed<Label>::kExtension)) {

            auto label_encoder = RowCompressed<Label>::load_label_encoder(filename);
            label_encoders.push_back(label_encoder.get());
            loaded_label_encoders.push_back(std::move(label_encoder));

            auto annotator = std::make_unique<class RowCompressed<Label>::StreamRows>(filename);
            streams.push_back(std::move(annotator));

        } else {
            throw std::runtime_error("streaming only supported for rowcompressed annotator");
        }
    }

    merge_rows(
        label_encoders,
        [&](uint64_t annotator_idx) -> const std::vector<uint64_t> {
            if (annotator_idx < annotators.size()) {
                return annotator_row_iterators.at(annotator_idx)->next_row();
            } else {
                return *streams[annotator_idx-annotators.size()]->next_row();
            }
        },
        num_rows,
        [&](const LEncoder &merged_label_enc,
                std::function<void(BinaryMatrix::RowCallback)> callback) {
            return RowCompressed<Label>::write_rows(outfile,
                                                    merged_label_enc,
                                                    callback);
        }
    );

    if constexpr (!std::is_same<RowCompressed<Label>, ToAnnotation>::value) {
        auto out_annotator = convert<ToAnnotation, Label>(outfile);
        out_annotator->serialize(outfile);
    }

    return num_rows;
}

} // namespace annotate

#endif // __ANNOTATION_CONVERTERS_HPP__
