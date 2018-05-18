#include "kmer.hpp"

#include <sdsl/uint256_t.hpp>

const KMerBaseType kFirstCharMask = (1 << kBitsPerChar) - 1;


bool KMer::compare_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
             == k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}

KMerCharType KMer::operator[](size_t i) const {
    assert(get_digit<kBitsPerChar>(i) > 0);
    return get_digit<kBitsPerChar>(i) - 1;
}

std::string KMer::to_string(const std::string &alphabet) const {
    std::string seq;
    seq.reserve(sizeof(KMerBaseType) * 8 / kBitsPerChar + 1);

    KMerCharType cur;
    for (size_t i = 0; (cur = get_digit<kBitsPerChar>(i)); ++i) {
        seq.push_back(alphabet.at(cur - 1));
    }
    return seq;
}

std::ostream& operator<<(std::ostream &os, const KMer &kmer) {
    return os << sdsl::uint256_t(kmer.seq_);
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
void KMer::update_kmer(size_t k,
                       KMerCharType edge_label,
                       KMerCharType last,
                       KMerBaseType *kmer) {
    *kmer = *kmer >> kBitsPerChar;
    *kmer += KMerBaseType(last + 1) << static_cast<int>(kBitsPerChar * k);
    *kmer |= kFirstCharMask;
    *kmer -= kFirstCharMask;
    *kmer += edge_label + 1;
}
