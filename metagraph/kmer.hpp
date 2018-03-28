#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include <sdsl/uint256_t.hpp>

using sdsl::uint256_t;
typedef uint64_t KMerCharType;

#ifdef _PROTEIN_GRAPH
const int kBitsPerChar = 5;
#else
const int kBitsPerChar = 3;
#endif


class KMer {
    friend std::ostream& operator<<(std::ostream &os, const KMer &kmer);

  public:
    KMer() {}
    template <typename T>
    KMer(const T &arr, size_t k)
       : seq_(pack_kmer(arr, k)) {}

    template <class Map, class String>
    KMer(const String &seq, Map &&to_alphabet);

    KMer(KMer &&other) : seq_(other.seq_) {}
    KMer(const KMer &other) : seq_(other.seq_) {}
    explicit KMer(sdsl::uint256_t &&seq) : seq_(seq) {}
    explicit KMer(const sdsl::uint256_t &seq) : seq_(seq) {}

    KMer& operator=(KMer &&other) { seq_ = other.seq_; return *this; }
    KMer& operator=(const KMer &other) { seq_ = other.seq_; return *this; }

    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    KMerCharType operator[](size_t i) const;

    template <size_t bits_per_digit>
    uint64_t get_digit(size_t i) const;

    static bool compare_suffix(const KMer &k1,
                               const KMer &k2, size_t minus = 0);

    std::string to_string(const std::string &alphabet) const;

    template<typename T>
    static sdsl::uint256_t pack_kmer(const T &arr, size_t k);

    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     */
    static void update_kmer(size_t k,
                            KMerCharType edge_label,
                            KMerCharType last,
                            sdsl::uint256_t *kmer);
  private:
    sdsl::uint256_t seq_; // kmer sequence
};

template <typename T>
sdsl::uint256_t KMer::pack_kmer(const T &arr, size_t k) {
    sdsl::uint256_t result(0);
    assert(k * kBitsPerChar < 256 && k >= 2
            && "String must be between lengths 2 and 256 / kBitsPerChar");

    for (int i = k - 2; i >= 0; --i) {
        assert(static_cast<uint64_t>(arr[i] + 1) < (1llu << kBitsPerChar)
                 && "Alphabet size too big for the given number of bits");

        result = result << kBitsPerChar;
        result += arr[i] + 1;
    }
    result = result << kBitsPerChar;
    result += arr[k - 1] + 1;
    return result;
}

template <class Map, class String>
KMer::KMer(const String &seq, Map &&to_alphabet) {
    std::vector<uint8_t> arr(seq.size());
    std::transform(seq.begin(), seq.end(), arr.begin(), to_alphabet);
    *this = KMer(arr.data(), arr.size());
}

template <size_t digit_size>
uint64_t KMer::get_digit(size_t i) const {
    static_assert(digit_size <= 64, "too big digit");
    return static_cast<uint64_t>(seq_ >> (digit_size * i))
             % (1llu << digit_size);
}

#endif // __KMER_HPP__
