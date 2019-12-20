#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include "common/hash/rolling_hasher.hpp"
#include "common/hash/bloom_filter.hpp"


/**
 * BloomFilter for inserting and checking k-mers extracted from sequences.
 * When in canonical mode, all k-mers are converted to their canonical forms
 * before inserting and checking.
 */
template <class KmerHasher = RollingHash<>>
class KmerBloomFilter {
  public:
    typedef KmerHasher KmerHasherType;

    /**
     * Construction a KmerBloomFilter given k and the desired canonical mode.
     * See BloomFilter constructors for possible combinations of Bloom filter
     * parameters.
     */
    template <typename... Args>
    KmerBloomFilter(size_t k, bool canonical_mode, Args&&... args)
          : filter_(std::forward<Args>(args)...),
            canonical_mode_(canonical_mode),
            k_(k),
            hasher_(k_) {}

    /**
     * Insert the k-mers of a sequence into the Bloom filter.
     */
    void add_sequence(const char *begin, const char *end);

    /**
     * Insert the k-mers of a sequence into the Bloom filter.
     */
    void add_sequence(const std::string &sequence) {
        add_sequence(sequence.data(), sequence.data() + sequence.size());
    }

    /**
     * Insert the k-mers of the generated sequences into the Bloom filter.
     */
    typedef std::function<void(const std::string&)> CallString;
    void add_sequences(const std::function<void(const CallString&)> &generate_sequences);

    /**
     * Check the k-mers of a sequence for presence/absence in the Bloom filter
     */
    sdsl::bit_vector check_kmer_presence(const char *begin, const char *end) const;

    /**
     * Check the k-mers of a sequence for presence/absence in the Bloom filter
     */
    sdsl::bit_vector check_kmer_presence(const std::string &sequence) const {
        return check_kmer_presence(sequence.data(), sequence.data() + sequence.size());
    }

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }
    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); }

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    const KmerHasher& get_hasher() const { return hasher_; }

    const BloomFilter& get_filter() const { return filter_; }

  private:
    BloomFilter filter_;
    bool canonical_mode_;
    size_t k_;
    const KmerHasher hasher_;
};


/**
 * Construct a callback which can be called `end`-`begin`+1 times, where
 * the i-th call returns `true` if the i-th k-mer in sequence is invalid
 * and `false` otherwise.
 */
std::function<bool()> get_missing_kmer_skipper(const KmerBloomFilter<> *bloom_filter,
                                               const char *begin,
                                               const char *end);

#endif // __KMER_BLOOM_FILTER_HPP__
