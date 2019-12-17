#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#define private public
#define protected public

#include "kmer_boss.hpp"
#include "kmer_extractor.hpp"


template <typename T>
T encode_c(char c, const T *char_map) {
    assert(static_cast<size_t>(c) < 128);
    return char_map[static_cast<size_t>(c)];
}

template <typename T>
char decode_c(T a, const std::string &alphabet) {
    assert(a < alphabet.size());
    return alphabet[a];
}

template <typename T>
std::vector<T> encode_c(const std::string &sequence, const T *char_map) {
    std::vector<T> encoded;
    std::transform(sequence.begin(), sequence.end(),
                   std::back_inserter(encoded),
                   [&](char c) { return encode_c(c, char_map); });
    assert(encoded.size() == sequence.size());

    return encoded;
}

template <typename T>
std::string decode_c(const std::vector<T> &encoded, const std::string &alphabet) {
    std::string decoded;
    std::transform(encoded.begin(), encoded.end(),
                   std::back_inserter(decoded),
                   [&](T a) { return decode_c(a, alphabet); });
    assert(decoded.size() == encoded.size());

    return decoded;
}

typedef uint8_t TAlphabet;
template std::vector<TAlphabet> encode_c(const std::string &sequence, const TAlphabet *char_map);


 // Nucleotide
std::vector<TAlphabet> encode_nucleotide(const std::string &sequence) {
    const TAlphabet kCharToNucleotide[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    return encode_c(sequence, kCharToNucleotide);
}

template <typename G, int L>
std::string decode_nucleotide(const KMerBOSS<G, L> &kmer, size_t k) {
    return kmer.to_string(k, "ACGTN");
}


template <class KMER>
class KmerBOSS : public ::testing::Test { };
typedef ::testing::Types<uint64_t,
                         sdsl::uint128_t,
                         sdsl::uint256_t> IntTypes;
TYPED_TEST_CASE(KmerBOSS, IntTypes);

TYPED_TEST(KmerBOSS, nucleotide_alphabet_pack) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    const auto encoded = encode_nucleotide(sequence);

    for (uint64_t k = 2; k < sizeof(TypeParam) * 8 / 3; ++k) {
        if (k > encoded.size())
            break;

        ASSERT_LE(k, encoded.size());
        std::vector<KMerBOSS<TypeParam, 3>> kmers;
        KMerBOSS<TypeParam, 3> kmer_packed(encoded.data(), k);
        for (uint64_t i = 0; i + k <= encoded.size(); ++i) {
            kmers.emplace_back(std::vector<TAlphabet>(encoded.begin() + i,
                                                      encoded.begin() + i + k));
            auto decoded = decode_nucleotide<TypeParam, 3>(kmers.back(), k);
            EXPECT_EQ(sequence.substr(i, k), decoded) << k << " " << i;

            KMerBOSS<TypeParam, 3> kmer_alt(encoded.data() + i, k);
            EXPECT_EQ(kmers.back(), kmer_alt) << k << " " << i;

            EXPECT_EQ(kmers.back(), kmer_packed) << k << " " << i;

            if (i + k < encoded.size())
                kmer_packed.to_next(k, encoded[i + k], encoded[i + k - 1]);
        }
    }
}

//typedef uint64_t KMerBaseType;
const size_t kBitsPerChar = KmerExtractorBOSS::bits_per_char;

template <typename TypeParam>
using KMER = KMerBOSS<TypeParam, kBitsPerChar>;

#define kSizeOfKmer ( sizeof(TypeParam) )
//typedef KMerBOSS<KMerBaseType, kBitsPerChar> KMER;
//const size_t kSizeOfKmer = sizeof(KMerBaseType);

template <typename KMER>
std::string kmer_codec(const std::string &test_kmer) {
    std::vector<uint64_t> kmer(test_kmer.size());
    std::transform(test_kmer.begin(), test_kmer.end(), kmer.begin(),
        [](char c) {
            return c == KmerExtractorBOSS::alphabet[0]
                        ? 0
                        : KmerExtractorBOSS::encode(c);
        }
    );
    return KMER(kmer).to_string(test_kmer.length(), KmerExtractorBOSS::alphabet);
}

template <typename TypeParam>
void test_kmer_codec(const std::string &test_kmer,
                     const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_codec<KMER<TypeParam>>(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_codec<KMER<TypeParam>>(test_kmer));
}

TYPED_TEST(KmerBOSS, Invertible) {
    test_kmer_codec<TypeParam>("ATGG", "ATGG");
}

TYPED_TEST(KmerBOSS, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    //test bit shifting
    KMER<TypeParam> kmer_builtup(KmerExtractorBOSS::encode(
        std::string(long_seq.rbegin() + 1,
                    long_seq.rbegin() + 3)
    ));
    for (int i = long_seq.length() - 4; i >= 0; --i) {
        kmer_builtup.seq_ <<= static_cast<uint64_t>(kBitsPerChar);
        kmer_builtup.seq_ |= KmerExtractorBOSS::encode(long_seq[i]);
    }
    kmer_builtup.seq_ <<= static_cast<uint64_t>(kBitsPerChar);
    kmer_builtup.seq_ |= KmerExtractorBOSS::encode(long_seq[long_seq.length() - 1]);
    std::string dec = kmer_builtup.to_string(long_seq.length(), KmerExtractorBOSS::alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_codec<TypeParam>(long_seq, long_seq);
}

TYPED_TEST(KmerBOSS, UpdateKmer) {
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode("ATGC")),
        KMER<TypeParam>(KmerExtractorBOSS::encode("TGCT"))
    };
    KMER<TypeParam> updated = kmer[0];
    updated.to_next(4, KmerExtractorBOSS::encode('T'), KmerExtractorBOSS::encode('C'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractorBOSS::encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TYPED_TEST(KmerBOSS, NextPrevKmer) {
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode("ATGC")),
        KMER<TypeParam>(KmerExtractorBOSS::encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractorBOSS::encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, KmerExtractorBOSS::encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TYPED_TEST(KmerBOSS, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode(long_seq)),
        KMER<TypeParam>(KmerExtractorBOSS::encode(long_seq_alt))
    };
    kmer[0].to_next(long_seq.length(), KmerExtractorBOSS::encode('T'),
                                       KmerExtractorBOSS::encode(long_seq.back()));
    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), KmerExtractorBOSS::alphabet),
              kmer[0].to_string(long_seq.length(), KmerExtractorBOSS::alphabet));
}

TYPED_TEST(KmerBOSS, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = KmerExtractorBOSS::encode(long_seq0);
    KMER<TypeParam> kmer0(seq0.begin(), seq0.size());
    kmer0.to_next(long_seq0.length(), KmerExtractorBOSS::encode(long_seq1.back()),
                                      KmerExtractorBOSS::encode(long_seq0.back()));
    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               KmerExtractorBOSS::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(KmerExtractorBOSS::encode(long_seq1.back()));
    KMER<TypeParam> kmer1(seq0.begin() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               KmerExtractorBOSS::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TYPED_TEST(KmerBOSS, InvertibleEndDol) {
    test_kmer_codec<TypeParam>("ATG$", "ATG$");
}

TYPED_TEST(KmerBOSS, InvertibleStartDol) {
    test_kmer_codec<TypeParam>("$ATGG", "$ATGG");
}

TYPED_TEST(KmerBOSS, InvertibleBothDol) {
    test_kmer_codec<TypeParam>("$ATG$", "$ATG$");
}

TYPED_TEST(KmerBOSS, InvalidChars) {
#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_codec<TypeParam>("ATGH", "ATGN");
    test_kmer_codec<TypeParam>("ATGЯ", "ATGNN");
#elif _PROTEIN_GRAPH
    test_kmer_codec<TypeParam>("ATGH", "ATGH");
    test_kmer_codec<TypeParam>("ATGЯ", "ATGXX");
#elif _DNA_GRAPH
    ASSERT_DEATH(test_kmer_codec<TypeParam>("ATGH", "ATGN"), "");
    ASSERT_DEATH(test_kmer_codec<TypeParam>("ATGЯ", "ATGNN"), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

template <typename TypeParam>
void test_kmer_less(const std::string &k1,
                    const std::string &k2, bool truth) {
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode(k1)),
        KMER<TypeParam>(KmerExtractorBOSS::encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TYPED_TEST(KmerBOSS, LessEdge) {
    test_kmer_less<TypeParam>("ATGC", "ATGG", true);
}

TYPED_TEST(KmerBOSS, Less) {
    test_kmer_less<TypeParam>("ACTG", "GCTG", true);
}

TYPED_TEST(KmerBOSS, LessLong) {
    test_kmer_less<TypeParam>(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "C",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "T",
        true
    );

    test_kmer_less<TypeParam>(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "CA",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "TA",
        true
    );
}

template <typename TypeParam>
void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode(k1)),
        KMER<TypeParam>(KmerExtractorBOSS::encode(k2))
    };
    ASSERT_EQ(truth, KMER<TypeParam>::compare_suffix(kmer[0], kmer[1], 1));
}

TYPED_TEST(KmerBOSS, CompareSuffixTrue) {
    test_kmer_suffix<TypeParam>("ACTG", "GCTG", true);
}

TYPED_TEST(KmerBOSS, CompareSuffixFalse) {
    test_kmer_suffix<TypeParam>("ATTG", "ACTG", false);
}

TYPED_TEST(KmerBOSS, CompareSuffixTrueLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[0] = 'T';
    KMER<TypeParam> kmer[2] = {
        KMER<TypeParam>(KmerExtractorBOSS::encode(long_seq)),
        KMER<TypeParam>(KmerExtractorBOSS::encode(long_seq_alt))
    };
    ASSERT_TRUE(KMER<TypeParam>::compare_suffix(kmer[0], kmer[1], 1));

    //shift, then compare
    long_seq_alt[kSizeOfKmer * 8 / kBitsPerChar - 2] = 'T';

    kmer[0].seq_
        = kmer[0].seq_ >> static_cast<int>((kSizeOfKmer * 8 / kBitsPerChar - 2)
                                                * kBitsPerChar);

    kmer[1] = KMER<TypeParam>(KmerExtractorBOSS::encode(
        long_seq_alt.substr(kSizeOfKmer * 8 / kBitsPerChar - 2)
    ));

    ASSERT_TRUE(KMER<TypeParam>::compare_suffix(kmer[0], kmer[1], 1));
}

TYPED_TEST(KmerBOSS, CompareSuffixFalseLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[1] = 'T';

    test_kmer_suffix<TypeParam>(long_seq, long_seq_alt, false);
}

TYPED_TEST(KmerBOSS, SizeOfClass) {
    EXPECT_EQ(kSizeOfKmer, sizeof(KMER<TypeParam>));
}


TEST(KmerBOSS, TestPrint64) {
    size_t size = sizeof(uint64_t) * 8 / kBitsPerChar;
    KMER<uint64_t> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("0000000000000000000000000000000000000000000000001249249249249249", out);
#endif
}

TEST(KmerBOSS, TestPrint128) {
    size_t size = sizeof(sdsl::uint128_t) * 8 / kBitsPerChar;
    KMER<sdsl::uint128_t> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("0000000000000000000000000000000009249249249249249249249249249249", out);
#endif
}

TEST(KmerBOSS, TestPrint256) {
    size_t size = sizeof(sdsl::uint256_t) * 8 / kBitsPerChar;
    KMER<sdsl::uint256_t> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("1249249249249249249249249249249249249249249249249249249249249249", out);
#endif
}
