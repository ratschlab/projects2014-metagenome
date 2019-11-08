#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>
#include <vector>
#include <deque>
#include <functional>
#include <stdexcept>
#include <queue>
#include <utility>
#include <set>
#include <bitset>

#if _USE_FOLLY
#include <folly/FBVector.h>
#include <folly/small_vector.h>
    template <typename... Args>
    using Vector = folly::fbvector<Args...>;

    template <typename T, size_t NumReserved = 2, typename SizeType = uint32_t>
    using SmallVector = folly::small_vector<T, NumReserved, SizeType>;
#else
    template <typename... Args>
    using Vector = std::vector<Args...>;

    template <typename T, size_t NumReserved = 2, typename SizeType = uint32_t>
    using SmallVector = std::vector<T>;
#endif

// Branch prediction helper macros
#ifndef LIKELY
#define LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)
#endif

#include "serialization.hpp"


struct VectorHash {
    template <class Vector>
    std::size_t operator()(const Vector &vector) const {
        uint64_t hash = 0;
        for (uint64_t value : vector) {
            hash ^= value + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return static_cast<std::size_t>(hash);
    }
};


namespace utils {

    template <typename T>
    struct is_pair : std::false_type {};

    template <typename T1, typename T2>
    struct is_pair<std::pair<T1,T2>> : std::true_type {};

    static_assert(is_pair<std::pair<int,size_t>>::value);
    static_assert(is_pair<std::pair<uint64_t,size_t>>::value);

    static_assert(!is_pair<std::tuple<int,size_t>>::value);
    static_assert(!is_pair<int>::value);

    template <typename T, typename... Us>
    inline const T& get_first(const std::tuple<T, Us...> &tuple) { return std::get<0>(tuple); }

    template <typename T, typename... Us>
    inline T& get_first(std::tuple<T, Us...> &tuple) { return std::get<0>(tuple); }

    template <typename T, typename U>
    inline const T& get_first(const std::pair<T, U> &pair) { return pair.first; }

    template <typename T, typename U>
    inline T& get_first(std::pair<T, U> &pair) { return pair.first; }

    template <typename T>
    inline T& get_first(T &value) { return value; }

    template <typename T>
    struct LessFirst {
        bool operator()(const T &p1, const T &p2) const {
            return get_first(p1) < get_first(p2);
        }
    };

    template <typename T>
    struct EqualFirst {
        bool operator()(const T &p1, const T &p2) const {
            return get_first(p1) == get_first(p2);
        }
    };

    bool get_verbose();
    void set_verbose(bool verbose);

    /**
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    // For each occurrence of |label| in |array|, mark segment
    // of length |segment_length| starting at that position.
    // Example: [X]***[X]****[X]***
    //    --->  [111]0[111]00[111]0
    template <class Vector>
    inline std::vector<bool> drag_and_mark_segments(const Vector &array,
                                                    typename Vector::value_type label,
                                                    size_t segment_length) {
        std::vector<bool> mask(array.size(), false);
        size_t last_occurrence
            = std::find(array.data(), array.data() + array.size(), label)
                - array.data();

        for (size_t i = last_occurrence; i < array.size(); ++i) {
            if (array[i] == label)
                last_occurrence = i;

            if (i - last_occurrence < segment_length)
                mask[i] = true;
        }

        return mask;
    }

    inline uint32_t code_length(uint64_t x) { return sdsl::bits::hi(x) + 1; }

    inline uint64_t max_uint(uint8_t width) { return ~uint64_t(0) >> (64 - width); }

    template <class AIt, class BIt>
    uint64_t count_intersection(AIt first_begin, AIt first_end,
                                BIt second_begin, BIt second_end) {
        assert(std::is_sorted(first_begin, first_end));
        assert(std::is_sorted(second_begin, second_end));
        assert(std::set<typename AIt::value_type>(first_begin, first_end).size()
                    == static_cast<uint64_t>(std::distance(first_begin, first_end)));
        assert(std::set<typename BIt::value_type>(second_begin, second_end).size()
                    == static_cast<uint64_t>(std::distance(second_begin, second_end)));

        uint64_t count = 0;

        while (first_begin != first_end && second_begin != second_end) {
            first_begin = std::lower_bound(first_begin, first_end, *second_begin);

            if (first_begin == first_end)
                break;

            second_begin = std::lower_bound(second_begin, second_end, *first_begin);

            if (second_begin == second_end)
                break;

            if (*first_begin == *second_begin) {
                ++count;
                ++first_begin;
                ++second_begin;
            }
        }

        return count;
    }

    template <typename T>
    struct Hash {
        size_t operator()(const T &x) const {
            return hasher(reinterpret_cast<const std::bitset<sizeof(T) * 8>&>(x));
        }

        std::hash<std::bitset<sizeof(T) * 8>> hasher;
    };

    // Bitmap |new_indexes| marks positions of inserted values in the final vector
    template <class Vector>
    void insert(Vector *vector,
                const bitmap &new_indexes,
                const typename Vector::value_type &value) {
        assert(vector);
        assert(new_indexes.size() == vector->size() + new_indexes.num_set_bits());

        if (!new_indexes.num_set_bits())
            return;

        vector->resize(vector->size() + new_indexes.num_set_bits());

        // need to move everything to the right as we call ones from left to right
        std::move_backward(vector->begin(),
                           vector->end() - new_indexes.num_set_bits(),
                           vector->end());

        size_t i = new_indexes.num_set_bits();
        size_t curpos = 0;

        // call ones from left to right
        new_indexes.call_ones([&](auto new_index) {
            assert(new_index >= curpos);

            // move block
            while (curpos < new_index) {
                (*vector)[curpos++] = std::move((*vector)[i++]);
            }
            // insert new value
            (*vector)[curpos++] = value;
        });
    }

    // new_indexes - positions of inserted values in the final vector
    template <class Vector>
    void insert(Vector *vector,
                const std::vector<uint64_t> &new_indexes,
                const typename Vector::value_type &value) {
        assert(vector);
        assert(std::is_sorted(new_indexes.begin(), new_indexes.end()));
        assert(!new_indexes.size()
                    || new_indexes.back() < vector->size() + new_indexes.size());

        if (!new_indexes.size())
            return;

        vector->resize(vector->size() + new_indexes.size());

        uint64_t i = vector->size() - 1;
        uint64_t shift = new_indexes.size();

        for (auto it = new_indexes.rbegin(); it != new_indexes.rend(); ++it) {
            while (i > *it) {
                assert(i >= shift && "Invalid indexes for insertion");
                (*vector)[i] = std::move((*vector)[i - shift]);
                i--;
            }
            // insert new value
            shift--;
            (*vector)[i--] = value;
        }
    }

    template <class Array, class Mask>
    void erase(Array *vector, const Mask &erase_mask) {
        assert(vector);
        assert(vector->size() == erase_mask.size());

        size_t j = 0;
        for (size_t i = 0; i < erase_mask.size(); ++i) {
            if (!erase_mask[i])
                (*vector)[j++] = (*vector)[i];
        }
        vector->resize(j);
    }

    // Read indices of set bits from a vector of VectorStreams
    class RowsFromColumnsTransformer {
      public:
        // Files store serialized vectors as plain indexes of the set bits
        RowsFromColumnsTransformer(uint64_t num_rows,
                                   const std::vector<std::string> &files);

        template <typename BitVectorPtr>
        RowsFromColumnsTransformer(const std::vector<BitVectorPtr> &columns);

        RowsFromColumnsTransformer(const bit_vector_small &columns_concatenated,
                                   uint64_t column_size);

        using ValueCallback = std::function<void(uint64_t /*row*/,
                                                 uint64_t /*column*/)>;
        void call_next(ValueCallback callback);

        uint64_t rows() const { return num_rows_; }
        uint64_t columns() const { return streams_.size(); }

        // get the number of set bits in the vectors left
        uint64_t values_left() const { return num_set_bits_left_; }

      private:
        void init_heap();

        std::vector<std::unique_ptr<VectorStream>> streams_;
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
        RowsFromColumnsIterator(std::unique_ptr<utils::RowsFromColumnsTransformer> transformer) {
            transformer_ = std::move(transformer);
        }

        std::tuple<uint64_t, uint64_t> next_set_bit();
        uint64_t values_left() { return transformer_->values_left(); };

        std::vector<uint64_t> next_row();

      private:
        std::unique_ptr<utils::RowsFromColumnsTransformer> transformer_;

        uint64_t i_ = 0;
        uint64_t row_ = 0;
        uint64_t column_;
    };

    void call_rows(const std::function<void(const std::vector<uint64_t>&)> &callback,
                   RowsFromColumnsTransformer&& transformer);

    template <class BitVectorType = bit_vector_stat>
    std::vector<std::unique_ptr<bit_vector>>
    transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix);


    template <typename T>
    std::vector<T> arange(T first, size_t size) {
        std::vector<T> result(size);
        std::iota(result.begin(), result.end(), first);
        return result;
    }


    // indexes are distinct and sorted
    sdsl::bit_vector subvector(const bit_vector &col,
                               const std::vector<uint64_t> &indexes);

    std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                         uint64_t sample_size,
                                         std::mt19937 &gen);

    template<class T> struct dependent_false : std::false_type {};


    template <typename T>
    T get_quantile(const std::vector<std::pair<T, uint64_t>> &count_hist, double q) {
        assert(q >= 0.0);
        assert(q <= 1.0);

        assert(std::is_sorted(count_hist.begin(), count_hist.end(),
                              [](const auto &first, const auto &second) {
                                  return first.first < second.first;
                              }));

        uint64_t sum = 0;
        for (const auto &pair : count_hist) {
            sum += pair.second;
        }

        const double threshold = q * sum;
        uint64_t partial_sum = 0;

        for (const auto &pair : count_hist) {
            partial_sum += pair.second;
            if (partial_sum >= threshold)
                return pair.first;
        }

        assert(false);
        return count_hist.back().first;
    }

    template <typename T>
    class DequeStorage {
      public:
        typedef T value_type;
        using iterator = typename std::deque<T>::iterator;
        using const_iterator = typename std::deque<T>::const_iterator;

        template <class... Args>
        DequeStorage(Args&&... args) : deque_(std::forward<Args>(args)...),
                                       size_(deque_.size()) {}
        ~DequeStorage() {}

        DequeStorage& operator=(const DequeStorage &other) = default;
        DequeStorage& operator=(DequeStorage&& other) = default;

        inline void push_back(T&& value) {
            if (size_ == deque_.size())
                try_reserve(size_ * growth_factor, size_ + 1);

            deque_[size_++] = std::move(value);
        }

        inline void push_back(const T &value) {
            if (size_ == deque_.size())
                try_reserve(size_ * growth_factor, size_ + 1);

            deque_[size_++] = value;
        }

        template <class... Args>
        inline void emplace_back(Args&&... args) {
            if (size_ == deque_.size())
                try_reserve(size_ * growth_factor, size_ + 1);

            deque_[size_++] = T(std::forward<Args>(args)...);
        }

        inline iterator erase(iterator first, iterator last) {
            return __erase(first, last);
        }

        inline const_iterator erase(const_iterator first, const_iterator last) {
            return __erase(first, last);
        }

        inline void clear() {
            size_ = 0;
            deque_.clear();
        }

        inline void resize(size_t size) {
            reserve(size);
            size_ = size;
        }

        inline void reserve(size_t size) {
            if (size > deque_.size())
                deque_.resize(size);
        }

        inline void try_reserve(size_t size, size_t min_size = 0) {
            size = std::max(size, min_size);

            while (size > min_size) {
                try {
                    reserve(size);
                    return;
                } catch (const std::bad_alloc &exception) {
                    size = min_size + (size - min_size) * 2 / 3;
                }
            }
            reserve(min_size);
        }

        inline size_t size() const { return size_; }

        inline size_t capacity() const { return deque_.size(); }

        inline T& operator[](size_t pos) { return deque_[pos]; }
        inline const T& operator[](size_t pos) const { return deque_[pos]; }

        inline T& at(size_t pos) {
            if (pos > size_)
                throw std::out_of_range("Out of range error");
            return deque_[pos];
        }

        inline const T& at(size_t pos) const {
            if (pos > size_)
                throw std::out_of_range("Out of range error");
            return deque_[pos];
        }

        inline T& front() { return deque_.front(); }
        inline const T& front() const { return deque_.front(); }

        inline T& back() { return deque_[size_ - 1]; }
        inline const T& back() const { return deque_[size_ - 1]; }

        inline void shrink_to_fit() {
            deque_.resize(size_);
            deque_.shrink_to_fit();
        }

        inline auto begin() { return deque_.begin(); }
        inline auto end() { return deque_.begin() + size_; }

        inline auto begin() const { return deque_.begin(); }
        inline auto end() const { return deque_.begin() + size_; }

      private:
        template <typename Iterator>
        inline Iterator __erase(Iterator first, Iterator last) {
            assert(first <= last);

            if (last == end()) {
                size_ -= (last - first);
                return end();
            }

            size_t old_size = deque_.size();
            auto it = deque_.erase(first, last);
            size_ -= old_size - deque_.size();
            return it;
        }

        static constexpr double growth_factor = 3. / 2;
        std::deque<T> deque_;
        size_t size_;
    };

    class NoCleanup {
      public:
        template <class Array>
        static void cleanup(Array*) {}
    };

    // removes redundant dummy BOSS k-mers from a sorted list
    class DummyKmersCleaner {
      public:
        template <class Array>
        static void cleanup(Array *kmers) {
            using KMER = std::remove_reference_t<decltype(get_first(kmers->at(0)))>;

            assert(std::is_sorted(kmers->begin(), kmers->end(),
                                  utils::LessFirst<typename Array::value_type>()));
            assert(std::unique(kmers->begin(), kmers->end(),
                               utils::EqualFirst<typename Array::value_type>()) == kmers->end());

            if (kmers->size() < 2)
                return;

            // last k-mer is never redundant. Start with the next one.
            uint64_t last = kmers->size() - 1;

            typename KMER::CharType edge_label, node_last_char;

            std::vector<uint64_t> last_kmer(1llu << KMER::kBitsPerChar, kmers->size());

            last_kmer[get_first(kmers->at(last))[0]] = last;

            for (int64_t i = last - 1; i >= 0; --i) {
                const KMER &kmer = get_first(kmers->at(i));
                node_last_char = kmer[1];
                edge_label = kmer[0];

                // assert((edge_label || node_last_char)
                //             && "dummy k-mer cannot be both source and sink dummy");

                if (!edge_label) {
                    // sink dummy k-mer

                    // skip if redundant
                    if (node_last_char && KMER::compare_suffix(kmer, get_first(kmers->at(last)), 0))
                        continue;

                } else if (!node_last_char) {
                    // source dummy k-mer

                    // skip if redundant
                    if (last_kmer[edge_label] < kmers->size()
                            && KMER::compare_suffix(kmer, get_first(kmers->at(last_kmer[edge_label])), 1))
                        continue;
                }

                // the k-mer is either not dummy, or not redundant -> keep the k-mer
                kmers->at(--last) = kmers->at(i);
                last_kmer[edge_label] = last;
            }

            kmers->erase(kmers->begin(), kmers->begin() + last);
        }
    };

    template <typename, template <typename, typename...> class>
    struct is_instance : public std::false_type {};

    template <typename... Ts, template <typename, typename...> class U>
    struct is_instance<U<Ts...>, U> : public std::true_type {};

} // namespace utils

template <typename T>
std::set<T> convert_to_set(const std::vector<T> &vector);

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector);

#endif // __UTILS_HPP__
