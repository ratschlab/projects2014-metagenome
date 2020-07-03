#ifndef __ANNOTATION_MATRIX_HPP__
#define __ANNOTATION_MATRIX_HPP__

#include <memory>
#include <vector>

#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace anno {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelEncoded<Label> {
  public:
    typedef BinaryMatrixType binary_matrix_type;
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;

    explicit StaticBinRelAnnotator() : matrix_(new BinaryMatrixType()) {}

    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder);

    bool has_label(Index i, const Label &label) const override;
    bool has_labels(Index i, const VLabels &labels) const override;

    void serialize(const std::string &filename) const override;
    bool merge_load(const std::vector<std::string> &filenames) override;
    // Dump columns to separate files in human-readable format
    bool dump_columns(const std::string &prefix, uint64_t num_threads = 1) const;

    uint64_t num_objects() const override;
    uint64_t num_relations() const override;

    void set(Index, const VLabels &) override { except_dyn(); }
    void add_labels(const std::vector<Index> &, const VLabels &) override { except_dyn(); }
    void insert_rows(const std::vector<Index> &) override { except_dyn(); }

    const BinaryMatrixType& get_matrix() const override { return *matrix_; };

    std::string file_extension() const override;

    static const std::string kExtension;

  private:
    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    using MultiLabelEncoded<Label>::label_encoder_;
};

} // namespace anno
} // namespace mtg

#endif // __ANNOTATION_MATRIX_HPP__
