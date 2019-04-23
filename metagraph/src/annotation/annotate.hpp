#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <tuple>


namespace annotate {

// class GenomeAnnotation {
//   public:
//     typedef uint64_t Index;
//     enum Label {
//         OUTGOING_EDGE_LABELS = 0,
//         SEQUENCE_NAME,
//     };

//     virtual const Annotation& get(Index i, Label j) const = 0;
//     virtual void set(Index i, Label j, const Annotation &label) = 0;
//     virtual void merge(const GenomeAnnotation &annotation) = 0;

//     virtual bool load(const std::string &filename) = 0;
//     virtual void serialize(const std::string &filename) const = 0;
// };


template <typename Index, typename LabelType>
class AnnotationCategory {
  public:
    virtual ~AnnotationCategory() {}

    virtual LabelType get(Index i) const = 0;
    virtual void add(Index i, const LabelType &label) = 0;
    virtual void set(Index i, const LabelType &label) = 0;

    virtual void serialize(const std::string &filename) const = 0;
    virtual bool load(const std::string &filename) { return merge_load({ filename }); }
    virtual bool merge_load(const std::vector<std::string> &filenames) = 0;
};

// Graph annotation
// An annotated graph is a graph with labeled edges.
// A labeled edge is an edge carrying certain labels.
template <typename IndexType, typename LabelType>
class MultiLabelAnnotation
      : public AnnotationCategory<IndexType, std::vector<LabelType>> {

  public:
    typedef IndexType Index;
    typedef LabelType Label;
    typedef std::vector<Label> VLabels;

    virtual ~MultiLabelAnnotation() {}

    /***************** Inherited member functions ****************/

    virtual VLabels get(Index i) const override final {
        return get_labels(i);
    }
    virtual void set(Index i, const VLabels &labels) override final {
        set_labels(i, labels);
    }
    virtual void add(Index i, const VLabels &labels) override final {
        add_labels(i, labels);
    }

    /******************* General functionality *******************/

    virtual void set_labels(Index i, const VLabels &labels) = 0;
    virtual VLabels get_labels(Index i) const = 0;

    virtual void add_label(Index i, const Label &label) = 0;
    virtual void add_labels(Index i, const VLabels &labels) = 0;
    virtual void add_labels(const std::vector<Index> &indices,
                            const VLabels &labels) = 0;

    virtual bool has_label(Index i, const Label &label) const = 0;
    virtual bool has_labels(Index i, const VLabels &labels) const = 0;

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const std::unordered_map<Label, Label> &dict) = 0;

    /*********************** Special queries **********************/

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    virtual VLabels get_labels(const std::vector<Index> &indices,
                               double presence_ratio) const = 0;

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    // Skip labels with count frequency smaller than |min_label_frequency|.
    virtual std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1),
                   double min_label_frequency = 0.0) const = 0;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const = 0;
    virtual size_t num_labels() const = 0;
    virtual uint64_t num_relations() const = 0;
};


// A dictionary to encode annotation labels
template <typename Label = std::string>
class LabelEncoder {
  public:
    /**
     * If the label passed does not exist, insert
     * that label and return its code.
     */
    size_t insert_and_encode(const Label &label);

    /**
     * Return the code of the label passed.
     * Throw exception if it does not exist.
     */
    size_t encode(const Label &label) const;

    /**
     * Throws an exception if a bad code is passed.
     */
    const Label& decode(size_t code) const { return decode_label_.at(code); }

    size_t size() const { return decode_label_.size(); }

    bool load(std::istream &instream);
    void serialize(std::ostream &outstream) const;

    void merge(const std::vector<const LabelEncoder<Label>*> &label_encoders);

    void clear() { encode_label_.clear(); decode_label_.clear(); }

  private:
    std::unordered_map<Label, uint64_t> encode_label_;
    std::vector<Label> decode_label_;
};

template <typename IndexType, typename LabelType>
class IterateRows;

template <typename IndexType, typename LabelType>
class IterateRowsByIndex;


template <typename IndexType, typename LabelType>
class MultiLabelEncoded
      : public MultiLabelAnnotation<IndexType, LabelType> {
    template <class A, typename L>
    friend uint64_t merge(const std::vector<const MultiLabelEncoded<uint64_t, L>*>&, const std::vector<std::string>&, const std::string&);

    template <typename I, typename L>
    friend class IterateRowsByIndex;
  public:
    using Index = typename MultiLabelAnnotation<IndexType, LabelType>::Index;
    using Label = typename MultiLabelAnnotation<IndexType, LabelType>::Label;
    using VLabels = typename MultiLabelAnnotation<IndexType, LabelType>::VLabels;

    virtual ~MultiLabelEncoded() {}

    /******************* General functionality *******************/

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const std::unordered_map<Label, Label> &dict) override;

    /*********************** Special queries **********************/

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    // Skip labels with count frequency smaller than |min_label_frequency|.
    virtual std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1),
                   double min_label_frequency = 0.0) const override final;

    virtual std::unique_ptr<IterateRows<IndexType, LabelType> > iterator() const { 
        return std::move(std::make_unique<IterateRowsByIndex<IndexType, LabelType> >(*this));
    };

  protected:
    // TODO: add |min_label_frequency| parameter: return only frequent labels
    virtual std::vector<uint64_t>
    count_labels(const std::vector<Index> &indices) const = 0;

    LabelEncoder<Label> label_encoder_;

    virtual std::vector<uint64_t> get_label_indexes(Index i) const {
        VLabels labels = this->get_labels(i);
        std::vector<uint64_t> indexes;
        for (const auto &label : labels) {
            indexes.push_back(label_encoder_.encode(label));
        }
        return indexes;
    }
};

template <typename IndexType, typename LabelType>
class IterateRows {
  public:
    virtual std::vector<uint64_t> next_row() = 0;
};

template <typename IndexType, typename LabelType>
class IterateRowsByIndex : public IterateRows<IndexType, LabelType> {
  public:
    IterateRowsByIndex(const MultiLabelEncoded<IndexType, LabelType>& annotator) : annotator_(annotator) {};
    virtual std::vector<uint64_t> next_row() { return annotator_.get_label_indexes(i_++); };
  private:
    typename MultiLabelEncoded<IndexType, LabelType>::Index i_ = 0;
    const MultiLabelEncoded<IndexType, LabelType> &annotator_;
};

template <class SetBitsIterator, typename IndexType, typename LabelType>
class IterateRowsBySetBits : public IterateRows<IndexType, LabelType> {
  public:
    IterateRowsBySetBits(std::unique_ptr<SetBitsIterator> set_bits_iterator) :
        set_bits_iterator_(std::move(set_bits_iterator)) {};
    
    virtual std::vector<uint64_t> next_row() {
        std::vector<uint64_t> indices;

        if (i_ > 0 && (row_ == i_)) {
            indices.push_back(column_);
        }

        if (!set_bits_iterator_->values_left() || row_ > i_) {
            i_++;
            return indices;
        }

        while (true) {
            if (!set_bits_iterator_->values_left())
                break;
            std::tie(row_, column_) = set_bits_iterator_->next_set_bit();
            if (row_ != i_)
                break;
            indices.push_back(column_);
        }
        i_++;
        return indices;
    }
  protected:
    std::unique_ptr<SetBitsIterator> set_bits_iterator_;

    typename MultiLabelEncoded<IndexType, LabelType>::Index i_ = 0;
    typename MultiLabelEncoded<IndexType, LabelType>::Index row_ = 0;
    uint64_t column_;
};

} // namespace annotate

#endif // __ANNOTATE_HPP__
