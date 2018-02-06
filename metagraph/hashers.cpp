#include "hashers.hpp"
#include <cyclichash.h>

namespace annotate {

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                                      const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "ORing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] | b[i];
    }
    return merged;
}

std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                       const std::vector<uint64_t> &b) {

    assert(a.size() == b.size() && "ANDing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] & b[i];
    }
    return merged;
}

uint64_t popcount(const std::vector<uint64_t> &a) {
    uint64_t popcount = 0;
    for (auto value : a) {
        popcount += __builtin_popcountl(value);
    }
    return popcount;
}

bool equal(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "Checking different sizes");
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

bool test_bit(const std::vector<uint64_t> &a, const size_t col) {
    return (a[col >> 6] | (1llu << (col % 64))) == a[col >> 6];
}

void set_bit(std::vector<uint64_t> &a, const size_t col) {
    a[col >> 6] |= 1llu << (col % 64);
}

void print(const std::vector<uint64_t> &a) {
    for (auto it = a.begin(); it != a.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

//HASH ITERATORS

HashIterator::HashIterator(const std::string &sequence, const size_t num_hash, const size_t k)
  : seq_begin(sequence.c_str()),
    seq_cur(seq_begin),
    seq_end(sequence.c_str() + sequence.length()),
    hashes_(num_hash),
    k_(k) {
    assert(sequence.length() >= k_);
}

HashIterator::HashIterator(const size_t num_hash, const size_t k)
  : seq_begin(NULL),
    seq_cur(NULL),
    seq_end(NULL),
    hashes_(num_hash),
    k_(k) { }

MultiHash HashIterator::get_hash() {
    assert(hashes_.size());
    return MultiHash(operator*(), hashes_.size());
}

std::vector<MultiHash> HashIterator::generate_hashes() {
    std::vector<MultiHash> hashes;
    for (; *this != end(); operator++()) {
        hashes.emplace_back(operator*(), hashes_.size());
    }
    return hashes;
}

//MurmurHash
void MurmurHashIterator::compute_hashes() {
    for (size_t i = 0; i < hashes_.size(); ++i) {
        //use index as seed
        Murmur3Hasher(seq_cur, k_, i, &hashes_[i]);
    }
}

MurmurHashIterator& MurmurHashIterator::operator++() {
    if (seq_cur + k_ <= seq_end)
        compute_hashes();
    seq_cur++;
    return *this;
}

MurmurHashIterator::MurmurHashIterator(const std::string &kmer, const size_t num_hash)
  : HashIterator(num_hash, kmer.length()),
    cache_(kmer.begin(), kmer.end()) {
    for (size_t i = 0; i < hashes_.size(); ++i) {
        Murmur3Hasher(kmer.c_str(), kmer.length(), i, &hashes_[i]);
    }
    assert(kmer.length() == k_ || *this != end());
}

MurmurHashIterator& MurmurHashIterator::update(const char next) {
    cache_.push_back(next);
    cache_.pop_front();
    std::string kmer(cache_.begin(), cache_.end());
    for (size_t i = 0; i < hashes_.size(); ++i) {
        Murmur3Hasher(kmer.c_str(), kmer.length(), i, &hashes_[i]);
    }
    return *this;
}

MurmurHashIterator& MurmurHashIterator::reverse_update(const char prev) {
    cache_.push_front(prev);
    cache_.pop_back();
    std::string kmer(cache_.begin(), cache_.end());
    for (size_t i = 0; i < hashes_.size(); ++i) {
        Murmur3Hasher(kmer.c_str(), kmer.length(), i, &hashes_[i]);
    }
    return *this;
}

//CyclicHashIterator
void CyclicHashIterator::init(const std::string &sequence) {
    chashers_.reserve(hashes_.size());
    for (uint32_t j = 0; j < hashes_.size(); ++j) {
        chashers_.push_back(new CyclicHash<uint64_t>(k_, j, j + 1, 64lu));
        for (size_t i = 0; i < k_; ++i) {
            reinterpret_cast<CyclicHash<uint64_t>*>(chashers_.back())->eat(sequence[i]);
        }
        hashes_[j] = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_.back())->hashvalue;
    }
}

void CyclicHashIterator::compute_hashes() {
    assert(seq_cur + k_ <= seq_end);
    assert(seq_cur > seq_begin);
    for (size_t i = 0; i < hashes_.size(); ++i) {
        hashes_[i] = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->hashvalue;
    }
}

CyclicHashIterator& CyclicHashIterator::operator++() {
    if (seq_cur + k_ <= seq_end) {
        for (size_t i = 0; i < hashes_.size(); ++i) {
            reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->update(*(seq_cur - 1), *(seq_cur - 1 + k_));
        }
        compute_hashes();
    }
    seq_cur++;
    return *this;
}

CyclicHashIterator::CyclicHashIterator(const std::string &kmer, const size_t num_hash)
  : HashIterator(num_hash, kmer.length()),
    cache_(kmer.begin(), kmer.end()) {
    init(kmer);
    assert(kmer.length() == k_ || *this != end());
}

CyclicHashIterator::CyclicHashIterator(const std::string &sequence, const size_t num_hash, const size_t k)
  : HashIterator(sequence, num_hash, k) {
    init(sequence);
    seq_cur++;
    assert(sequence.length() == k_ || *this != end());
}

CyclicHashIterator::~CyclicHashIterator() {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        delete reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
    }
}

CyclicHashIterator& CyclicHashIterator::update(const char next) {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->update(cache_.front(), next);
        hashes_[i] = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->hashvalue;
    }
    cache_.pop_front();
    cache_.push_back(next);
    return *this;
}

CyclicHashIterator& CyclicHashIterator::reverse_update(const char prev) {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->reverse_update(prev, cache_.back());
        hashes_[i] = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->hashvalue;
    }
    cache_.pop_back();
    cache_.push_front(prev);
    return *this;
}

} // namespace annotate
