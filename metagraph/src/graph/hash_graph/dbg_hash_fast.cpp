#include "dbg_hash_fast.hpp"

#include <cassert>
#include <limits>

#include <tsl/ordered_set.h>
#include <libmaus2/util/NumberSerialisation.hpp>

#include "common/seq_tools/reverse_complement.hpp"
#include "serialization.hpp"
#include "bit_vector.hpp"
#include "utils/algorithms.hpp"
#include "utils/hash_utils.hpp"
#include "utils/string_utils.hpp"


template <typename KMER = KmerExtractor2Bit::Kmer64>
class DBGHashFastImpl : public DBGHashFast::DBGHashFastInterface {
    using Kmer = KMER;
    using KmerWord = typename KMER::WordType;
    using TAlphabet = KmerExtractor2Bit::TAlphabet;

#if _PROTEIN_GRAPH
    using Bits = uint32_t;
#elif _DNA_CASE_SENSITIVE_GRAPH
    using Bits = uint16_t;
#elif _DNA5_GRAPH
    using Bits = uint16_t;
#elif _DNA_GRAPH
    using Bits = unsigned char;
#else
    static_assert(false, "invalid or undefined alphabet");
#endif

    using KmerIndex = tsl::ordered_set<KmerWord,
                                     utils::Hash<KmerWord>,
                                     std::equal_to<KmerWord>,
                                     std::allocator<KmerWord>,
                                     std::vector<KmerWord, std::allocator<KmerWord>>,
                                     std::uint64_t>;

    using KmerIterator = typename KmerIndex::iterator;
    using KmerConstIterator = typename KmerIndex::const_iterator;

  public:
    explicit DBGHashFastImpl(size_t k,
                              bool canonical_mode,
                              bool packed_serialization,
                              size_t reserve = 0);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate) const;

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const;

    void call_nodes(const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &stop_early) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index) const;
    bool has_multiple_outgoing(node_index) const;

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index) const;

    node_index kmer_to_node(const std::string &kmer) const;

    std::string get_node_sequence(node_index node) const;

    size_t get_k() const { return k_; }
    bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const {
        uint64_t nnodes = 0;
        call_nodes([&](auto) { nnodes++; }, [](){ return false; });
        return nnodes;
    };
    uint64_t max_index() const { return kmers_.size() << kBitsPerChar; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const;

    const std::string& alphabet() const { return seq_encoder_.alphabet; }

    bool in_graph(node_index node) const;

  private:
    Vector<std::pair<Kmer, bool>> sequence_to_kmers(const std::string &sequence,
                                                    bool canonical = false) const {
        return seq_encoder_.sequence_to_kmers<Kmer>(sequence, k_, canonical);
    }

    KmerIterator find_kmer(const KmerWord &kmer_prefix) {
        KmerIterator find = kmers_.find(kmer_prefix);
        assert(find == kmers_.end()
                   || !in_graph(get_node_index(find))
                   || KMER(kmer_prefix) == get_kmer(get_node_index(find), find));
        return find;
    }

    KmerConstIterator find_kmer(const KmerWord &kmer_prefix) const {
        KmerConstIterator find = kmers_.find(kmer_prefix);
        assert(find == kmers_.end()
                   || !in_graph(get_node_index(find))
                   || KMER(kmer_prefix) == get_kmer(get_node_index(find), find));
        return find;
    }

    bool kmers_overlap(const Kmer &out_kmer, const Kmer &in_kmer) const {
        KMER overlap = out_kmer;
        overlap.to_next(k_, in_kmer[k_ - 1]);
        return overlap == in_kmer;
    }

    node_index get_node_index(const KmerConstIterator &iter) const;
    node_index get_node_index(const Kmer &kmer) const;

    size_t get_vec_index(const KmerConstIterator &iter) const {
        return iter - kmers_.begin();
    }

    KmerIterator get_iter(node_index node);
    KmerConstIterator get_const_iter(node_index node) const;

    const Kmer get_kmer(node_index node) const {
        return get_kmer(node, get_const_iter(node));
    }
    const Kmer get_kmer(node_index node,
                        KmerConstIterator kmer_iter) const;

    void print_internal_representation() const {
        for (auto it = kmers_.begin(); it != kmers_.end(); ++it) {
            const auto index = get_vec_index(it);
            std::cout << (uint64_t)(it - kmers_.begin()) << "," << get_node_index(it) << ":";
            std::cout << KMER(it.key()).to_string(k_ - 1, seq_encoder_.alphabet) << " ";
            std::cout << std::bitset<std::numeric_limits<Bits>::digits>(bits_[index]) << " " << std::endl;
        }
    }

    size_t k_;
    bool canonical_mode_;

    KmerIndex kmers_;
    std::vector<Bits> bits_;
    KmerExtractor2Bit seq_encoder_;

    bool packed_serialization_;

    const KmerWord kIgnoreLastCharMask;

    // bits layout:
    //     <--- leading 0's --->
    //     < 1-bit flag for "may contain source kmer" >
    //     < kBitsPerNullableChar-bits for "next_char + 1" >
    //     < alphabet_size length bit field for each possible k'th char in kmer>
    // e.g. 0 100 1001 for two kmers with the same (k-1) prefix ending in 00 and 11, which
    // are known not to be source kmers, and where the (k-1) prefix hash entry for the 11
    // outgoing edge is known to be in the previous slot in KmerIndex (next_char = 100 - 1 == 11)
    const Bits kLastCharNBits = seq_encoder_.alphabet.size();
    const Bits kLastCharMask = (Bits(1) << kLastCharNBits) - 1;
    const Bits kIncomingEdgesMask = kLastCharMask << seq_encoder_.alphabet.size();
    const size_t kBitsPerNullableChar = utils::code_length(seq_encoder_.alphabet.size());
    const Bits kNextCharMask = (Bits(1) << kBitsPerNullableChar) - 1;
    const Bits kMayContainSourceKmer = Bits(1) << (kLastCharNBits + kBitsPerNullableChar);

    static constexpr auto kExtension = DBGHashFast::kExtension;
    static constexpr auto kBitsPerChar = KMER::kBitsPerChar;
};

template <typename KMER>
DBGHashFastImpl<KMER>::DBGHashFastImpl(size_t k,
                                         bool canonical_mode,
                                         bool packed_serialization,
                                         size_t reserve)
      : k_(k),
        canonical_mode_(canonical_mode),
        packed_serialization_(packed_serialization),
        kIgnoreLastCharMask(~(KmerWord((1llu << kBitsPerChar) - 1) << int(kBitsPerChar * (k - 1)))) {

    kmers_.reserve(reserve);
    bits_.reserve(reserve);
};


template <typename KMER>
void DBGHashFastImpl<KMER>::add_sequence(const std::string &sequence,
                                         bit_vector_dyn *nodes_inserted) {
    assert(!nodes_inserted || nodes_inserted->size() == max_index() + 1);

    const auto &kmers = sequence_to_kmers(sequence);

    bool may_contain_source_kmer = true;
    bool next_iter = false;
    for (auto kmer_it = kmers.rbegin(); kmer_it != kmers.rend(); ++kmer_it) {
        const auto &kmer = kmer_it->first;

        may_contain_source_kmer = kmer_it + 1 == kmers.rend() || !kmers_overlap((kmer_it + 1)->first, kmer);

        Bits val = Bits(1) << kmer[k_ - 1];

        const auto &&key = KmerWord(kmer.data()) & kIgnoreLastCharMask;

        auto &&[iter, inserted] = kmers_.insert(key);
        if (inserted) {
            val |= next_iter ? (kmer[k_ - 1] + 1) << kLastCharNBits : 0;
            val |= may_contain_source_kmer ? kMayContainSourceKmer  : 0;
            bits_.push_back(val);
        }
        assert(iter != kmers_.end());
        assert(iter == find_kmer(key));

        if (!inserted) {
            const auto index = get_vec_index(iter);
            val |= may_contain_source_kmer ? kMayContainSourceKmer : 0;
            bits_[index] |= val;
        }

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);

        next_iter = inserted && !may_contain_source_kmer;
    }

    if (!canonical_mode_)
        return;

    may_contain_source_kmer = true;
    next_iter = false;

    auto rev_comp = sequence;
    reverse_complement(rev_comp.begin(), rev_comp.end());

    auto rev_kmers = sequence_to_kmers(rev_comp);

    for (auto kmer_it = rev_kmers.rbegin(); kmer_it != rev_kmers.rend(); ++kmer_it) {
        const auto &kmer = kmer_it->first;

        may_contain_source_kmer = kmer_it + 1 == kmers.rend() || !kmers_overlap((kmer_it + 1)->first, kmer);

        Bits val = Bits(1) << kmer[k_ - 1];

        const auto &&key = KmerWord(kmer.data()) & kIgnoreLastCharMask;

        auto &&[iter, inserted] = kmers_.insert(key);
        if (inserted) {
            val |= next_iter ? (kmer[k_ - 1] + 1) << kLastCharNBits : 0;
            val |= may_contain_source_kmer ? kMayContainSourceKmer : 0;
            bits_.push_back(val);
        }
        assert(iter != kmers_.end());
        assert(iter == find_kmer(key));

        if (!inserted) {
            const auto index = get_vec_index(iter);
            val |= may_contain_source_kmer ? kMayContainSourceKmer : 0;
            bits_[index] |= val;
        }

        if (iter != kmers_.end() && nodes_inserted)
            nodes_inserted->insert_bit(kmers_.size() - 1, true);

        next_iter = inserted && !may_contain_source_kmer;
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence.
// In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
template <typename KMER>
void DBGHashFastImpl<KMER>::map_to_nodes_sequentially(
                                std::string::const_iterator begin,
                                std::string::const_iterator end,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers({ begin, end })) {
        if (terminate())
            return;

        assert((find_kmer(KmerWord(kmer.data()) & kIgnoreLastCharMask) == kmers_.end())
               || get_node_index(kmer) == npos
               || kmer == get_kmer(get_node_index(kmer)));

        callback(is_valid ? get_node_index(kmer) : npos);
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
template <typename KMER>
void DBGHashFastImpl<KMER>::map_to_nodes(const std::string &sequence,
                                         const std::function<void(node_index)> &callback,
                                         const std::function<bool()> &terminate) const {
    for (const auto &[kmer, is_valid] : sequence_to_kmers(sequence, canonical_mode_)) {
        if (terminate())
            return;

        assert(!get_node_index(kmer) || kmer == get_kmer(get_node_index(kmer)));

        callback(is_valid ? get_node_index(kmer) : npos);
    }
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_outgoing_kmers(node_index node,
                                                 const OutgoingEdgeCallback &callback) const {
    assert(in_graph(node));

    KMER kmer = get_kmer(node);

    node_index next_kmer_base_index;
    KmerConstIterator next_kmer_prefix_it;

    const Bits val = bits_[get_vec_index(get_const_iter(node))];
    unsigned char next_char = ((val >> kLastCharNBits) & kNextCharMask) - 1;
    if (kmer[k_ - 1] == next_char) {

        next_kmer_base_index = ((((node - 1) >> kBitsPerChar) - 1) << kBitsPerChar) + 1;
        next_kmer_prefix_it = get_const_iter(next_kmer_base_index);

        assert(next_kmer_prefix_it == find_kmer(KmerWord(kmer.data()) >> kBitsPerChar));
    } else {
        next_kmer_prefix_it = find_kmer(KmerWord(kmer.data()) >> kBitsPerChar);

        if (next_kmer_prefix_it == kmers_.end())
            return;

        next_kmer_base_index = get_node_index(next_kmer_prefix_it);
    }

    const Bits next_val = bits_[get_vec_index(next_kmer_prefix_it)];

    for (size_t i = 0; i < kLastCharNBits; ++i) {
        if ((next_val >> i) & 1) {
            auto next_kmer_index = next_kmer_base_index + i;

            assert(in_graph(next_kmer_index));

            callback(next_kmer_index, seq_encoder_.decode(i));
        }
    }
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_incoming_kmers(node_index node,
                                                 const IncomingEdgeCallback &callback) const {
    assert(in_graph(node));


    const auto &kmer = get_kmer(node);

    for (char c : seq_encoder_.alphabet) {
        auto prev_kmer = kmer;
        prev_kmer.to_prev(k_, seq_encoder_.encode(c));

        auto prev = get_node_index(prev_kmer);
        if (prev != npos)
            callback(prev, c);
    }
    return;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::traverse(node_index node, char next_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_next(k_, seq_encoder_.encode(next_char));
    return get_node_index(kmer);
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    auto kmer = get_kmer(node);
    kmer.to_prev(k_, seq_encoder_.encode(prev_char));
    return get_node_index(kmer);
}

template <typename KMER>
void
DBGHashFastImpl<KMER>
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
}

template <typename KMER>
void
DBGHashFastImpl<KMER>
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

template <typename KMER>
size_t DBGHashFastImpl<KMER>::outdegree(node_index node) const {
    assert(in_graph(node));

    KMER kmer = get_kmer(node);

    kmer.to_next(k_, 0);

    const auto next_kmer_prefix_it
        = find_kmer(KmerWord(kmer.data()) & kIgnoreLastCharMask);

    if (next_kmer_prefix_it == kmers_.end())
        return 0;

    const auto val = bits_[get_vec_index(next_kmer_prefix_it)] & kLastCharMask;

#if _PROTEIN_GRAPH
    return __builtin_popcountl(val);
#elif _DNA_CASE_SENSITIVE_GRAPH
    return __builtin_popcount(val);
#elif _DNA5_GRAPH
    return __builtin_popcount(val);
#elif _DNA_GRAPH
    return __builtin_popcount(val);
#else
    static_assert(false, "invalid or undefined alphabet");
#endif
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::has_single_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) == 1;
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::has_multiple_outgoing(node_index node) const {
    assert(in_graph(node));

    return outdegree(node) > 1;
}

template <typename KMER>
size_t DBGHashFastImpl<KMER>::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;
    call_incoming_kmers(node, [&](auto, auto) { indegree++; });
    return indegree;
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::has_no_incoming(node_index node) const {
    assert(in_graph(node));


    if (!(bits_[get_const_iter(node) - kmers_.begin()] & kMayContainSourceKmer))
        return false;

    return indegree(node) == 0;
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::has_single_incoming(node_index node) const {
    assert(in_graph(node));

    return indegree(node) == 1;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::kmer_to_node(const std::string &kmer) const {
    assert(kmer.length() == k_);

    return get_node_index(seq_encoder_.encode(kmer));
}

template <typename KMER>
std::string DBGHashFastImpl<KMER>::get_node_sequence(node_index node) const {
    assert(in_graph(node));

    return seq_encoder_.kmer_to_sequence(get_kmer(node), k_);
}

class Serializer {
  public:
    explicit Serializer(std::ostream &os) : os_(os) {}

    template <class T>
    void operator()(const T &value) {
        os_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

  private:
    std::ostream &os_;
};

class Deserializer {
  public:
    explicit Deserializer(std::istream &is) : is_(is) {}

    template <class T>
    T operator()() {
        T value;
        is_.read(reinterpret_cast<char *>(&value), sizeof(T));
        return value;
    }

  private:
    std::istream &is_;
};

template <typename KMER>
void DBGHashFastImpl<KMER>::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    out.exceptions(out.badbit | out.failbit);

    serialize_number(out, k_);

    Serializer serializer(out);

    if (packed_serialization_) {
        serialize_number(out, kmers_.size());
        std::for_each(kmers_.begin(), kmers_.end(), serializer);
        serialize_number(out, bits_.size());
        std::for_each(bits_.begin(), bits_.end(), serializer);
    } else {
        serialize_number(out, std::numeric_limits<uint64_t>::max());
        kmers_.serialize(serializer);
        serialize_number(out, bits_.size());
        std::for_each(bits_.begin(), bits_.end(), serializer);
    }

    serialize_number(out, canonical_mode_);
}

template <typename KMER>
void DBGHashFastImpl<KMER>::serialize(const std::string &filename) const {
    std::ofstream out(utils::remove_suffix(filename, kExtension) + kExtension,
                      std::ios::binary);
    serialize(out);
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::load(std::istream &in) {
    if (!in.good())
        return false;

    in.exceptions(in.badbit | in.failbit | in.eofbit);

    kmers_.clear();
    bits_.clear();

    try {
        k_ = load_number(in);

        Deserializer deserializer(in);

        uint64_t tag = load_number(in);

        if (tag < std::numeric_limits<uint64_t>::max()) {
            packed_serialization_ = true;

            const auto size = tag;
            kmers_.reserve(size + 1);
            for (uint64_t i = 0; i < size; ++i) {
                kmers_.insert(deserializer.operator()<KmerWord>());
            }

            uint64_t size2 = load_number(in);
            bits_.reserve(size2 + 1);
            for (uint64_t i = 0; i < size2; ++i) {
                bits_.push_back(deserializer.operator()<Bits>());
            }

        } else {
            packed_serialization_ = false;
            kmers_ = KmerIndex::deserialize(deserializer, true);

            uint64_t size2 = load_number(in);
            bits_.reserve(size2 + 1);
            for (uint64_t i = 0; i < size2; ++i) {
                bits_.push_back(deserializer.operator()<Bits>());
            }
        }

        canonical_mode_ = load_number(in);

        return in.good();

    } catch (...) {
        return false;
    }
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (!dynamic_cast<const DBGHashFastImpl*>(&other))
        throw std::runtime_error("Not implemented");

    const auto &other_hash = *dynamic_cast<const DBGHashFastImpl*>(&other);

    if (this == &other_hash)
        return true;

    assert(k_ == other_hash.k_);
    assert(canonical_mode_ == other_hash.canonical_mode_);
    assert(kmers_.size() == other_hash.kmers_.size());

    return kmers_ == other_hash.kmers_;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::get_node_index(const KmerConstIterator &kmer_iter) const {
    return ((kmer_iter - kmers_.begin()) << kBitsPerChar) + 1;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::node_index
DBGHashFastImpl<KMER>::get_node_index(const Kmer &kmer) const {
    const auto find = find_kmer(KmerWord(kmer.data()) & kIgnoreLastCharMask);
    if (find == kmers_.end())
        return npos;

    if (!(Bits(1) & (bits_[get_vec_index(find)] >> kmer[k_ - 1])))
        return npos;

    node_index node = get_node_index(find) + kmer[k_ - 1];
    assert(get_kmer(node, find) == kmer);
    return node;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::KmerIterator
DBGHashFastImpl<KMER>::get_iter(node_index node) {
    KmerIterator iter = (kmers_.begin() + ((node - 1) >> kBitsPerChar));
    return iter;
}

template <typename KMER>
typename DBGHashFastImpl<KMER>::KmerConstIterator
DBGHashFastImpl<KMER>::get_const_iter(node_index node) const {
    KmerConstIterator iter = (kmers_.begin() + ((node - 1) >> kBitsPerChar));
    return iter;
}

template <typename KMER>
const KMER DBGHashFastImpl<KMER>
::get_kmer(node_index node,
           KmerConstIterator kmer_iter) const {
    assert(in_graph(node));

    KMER prefix = KMER(*kmer_iter);

    KmerWord c = (node - 1) & KmerWord((1llu << kBitsPerChar) - 1);

    return KMER(KmerWord(prefix.data()) | (c << static_cast<int>(kBitsPerChar * (k_ - 1))));
}

template <typename KMER>
bool DBGHashFastImpl<KMER>::in_graph(node_index node) const {
    assert(node > 0 && node <= max_index());

    auto it = get_const_iter(node);
    KmerWord c = (node - 1) & KmerWord((1llu << kBitsPerChar) - 1);
    return Bits(1) & (bits_[get_vec_index(it)] >> c);
}

template <typename KMER>
void DBGHashFastImpl<KMER>::call_nodes(const std::function<void(node_index)> &callback,
                               const std::function<bool()> &stop_early) const {
    for (auto iter = kmers_.begin(); iter != kmers_.end() && !stop_early(); ++iter) {
        const auto val = bits_[get_vec_index(iter)];
        for (size_t i = 0; i < kLastCharNBits; ++i) {
            if (val & (Bits(1) << i)) {
                callback(get_node_index(iter) + i);
            }
        }
    }
}

std::unique_ptr<DBGHashFast::DBGHashFastInterface>
DBGHashFast::initialize_graph(size_t k,
                                 bool canonical_mode,
                                 bool packed_serialization) {
    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer64>>(
            k, canonical_mode, packed_serialization
        );
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer128>>(
            k, canonical_mode, packed_serialization
        );
    } else {
        return std::make_unique<DBGHashFastImpl<KmerExtractor2Bit::Kmer256>>(
            k, canonical_mode, packed_serialization
        );
    }
}

DBGHashFast::DBGHashFast(size_t k,
                               bool canonical_mode,
                               bool packed_serialization) {
    hash_dbg_ = initialize_graph(k, canonical_mode, packed_serialization);
}

bool DBGHashFast::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        auto pos = in.tellg();
        auto k = load_number(in);
        in.seekg(pos, in.beg);

        // the actual value of |canonical| will be set in load
        hash_dbg_ = initialize_graph(k, false, false);
        return hash_dbg_->load(in) && in.good();
    } catch (...) {
        return false;
    }
}

bool DBGHashFast::load(const std::string &filename) {
    std::ifstream in(utils::remove_suffix(filename, kExtension) + kExtension,
                     std::ios::binary);
    return load(in);
}
