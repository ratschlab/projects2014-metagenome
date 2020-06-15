#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>
#include "merge.hpp"
#include "common/elias_fano_file_merger.hpp"
#include "common/logger.hpp"
#include "common/sorted_multiset.hpp"
#include "common/sorted_multiset_disk.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_set_disk.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "kmer/kmer_collector.hpp"
#include "kmer/kmer_to_int_converter.hpp"
#include "boss_chunk.hpp"


namespace mtg {
namespace succinct {

using mtg::common::logger;
using mtg::common::ChunkedWaitQueue;
using mtg::kmer::get_int_t;
using mtg::kmer::get_kmer_t;
using mtg::kmer::KmerCollector;
using mtg::kmer::KmerExtractorBOSS;
using utils::get_first;
using utils::get_first_type_t;
using TAlphabet = KmerExtractorBOSS::TAlphabet;

constexpr size_t ENCODER_BUFFER_SIZE = 100'000;


template <class Container, typename KMER>
inline void push_back(Container &kmers, const KMER &kmer) {
    if constexpr(utils::is_pair_v<typename Container::value_type>) {
        kmers.emplace_back(kmer, 0);
    } else {
        kmers.push_back(kmer);
    }
}

/**
 * Generates non-redundant dummy sink kmers (a1a2...ak->$) for the given #kmers_p.
 *
 * The method traverses #kmers_p and generates the corresponding dummy sink k-mer for each
 * kmer in #kmers_p. The dummy sink k-mer is added if it's not redundant, i.e. if there is
 * no kmer in #kmers_p with the same prefix. Since #kmers_p are ordered, the generated
 * dummy sink k-mers for a given edge label will also be ordered. The method keeps an
 * iterator for each edge label that points to the first k-mer with prefix greater or
 * equal to the current dummy k-mer. The current dummy k-mer is marked as redundant if its
 * corresponding iterator points to a kmer with identical prefix.
 * @tparam T the type of kmers being processed, typically either a KMer64/128/256 or an
 * std::pair<KMer64/128/256, int8/16/32> if counting kmers.
 * @param k node length (the actual k-mer has length k+1)
 * @param kmers_p list of sorted, non-redundant kmers of length k+1
 */
template <typename T>
void add_dummy_sink_kmers(size_t k, Vector<T> *kmers_p) {
    using KMER = get_first_type_t<T>;
    using KMER_INT = typename KMER::WordType;

    const size_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    Vector<T> &kmers = *kmers_p;

    // points to the current k-mer with the given first character
    std::vector<size_t> max_it(alphabet_size);
    std::vector<size_t> it(alphabet_size);
    for (TAlphabet c = 1; c < alphabet_size; ++c) {
        std::vector<KMER_INT> zeros(k + 1, 0);
        zeros[k - 1] = c;
        it[c] = std::lower_bound(kmers.data(), kmers.data() + kmers.size(),
                                 KMER(zeros, k + 1), // the $$...i->$ k-mer
                                 [](const T &a, const KMER &b) -> bool {
                                     return get_first(a) < b;
                                 }) - kmers.data();
        max_it[c - 1] = it[c];
    }
    max_it[alphabet_size - 1] = kmers.size();

    std::vector<KMER> last_dummy(alphabet_size, KMER(0));
    size_t size = kmers.size();
    for (size_t i = 1; i < size; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // none of the original k-mers is a dummy k-mer
        assert(kmer[1] != 0 && kmer[0] != 0);

        KMER dummy_sink = kmer;
        dummy_sink.to_next(k + 1, BOSS::kSentinelCode);

        TAlphabet last_char = kmer[0];

        while (it[last_char] < max_it[last_char]
                && KMER::less(get_first(kmers[it[last_char]]), dummy_sink)) {
            it[last_char]++;
        }
        if (last_dummy[last_char] != dummy_sink
            && (it[last_char] == max_it[last_char]
                || !KMER::compare_suffix(get_first(kmers[it[last_char]]), dummy_sink))) {
            push_back(kmers, dummy_sink);
            last_dummy[last_char] = dummy_sink;
        }
    }
}

/**
 * Adds non-redundant dummy source nodes with prefix length 1 for the given kmers.
 *
 * For each kmer in #kmers_p, the method creates the corresponding dummy sources with
 * sentinels of length 1 (aka dummy-1 sources). The method then checks if the dummy-1
 * source node is redundant, i.e. if there is another k-mer that is identical
 * except for the first character. For example, $ACGT is redundant with TACGT. To do this
 * efficiently, the method uses an iterator called dummy_it that points to the first kmer
 * with the suffix equal to or larger than the current kmer. Since #kmers_p are ordered,
 * the generated dummy-1 kmers will also be ordered as long as the last character
 * in the original k-mers stays the same. This means that dummy_it only needs to move
 * forward and must be reset only when the last character in the original kmer changes.
 *
 * For example, the k-mer ACG->T, generates the dummy-kmer $AC->G, while the k-mer AGG->T
 * generates the dummy k-mer $AG->T. Because ACG->T precedes AGG->T, so will their
 * generated dummy-1 kmers. To check if $AC->G is redundant, dummy_it is advanced until we
 * find or pass k-mers with the suffix AC. Then we check all kmers with the AC suffix, and
 * if we find one with the same edge label, such as TAC->G, then $AC->G is redundant and
 * will be skipped.
 */
template <typename T>
void add_dummy_source_kmers(size_t k, Vector<T> *kmers_p, size_t end) {
    using KMER = get_first_type_t<T>;
    Vector<T> &kmers = *kmers_p;

    // points to the first k-mer that may be redundant with the current dummy source k-mer
    size_t real_it = 0;

    for (size_t i = 1; i < end; ++i) { // starting at 1 to skip the $$...$$ k-mer
        const KMER &kmer = get_first(kmers[i]);
        // none of the original k-mers is a dummy k-mer
        assert(kmer[1] != 0 && kmer[0] != 0);

        if (kmer[k] != get_first(kmers[i - 1])[k]) {
            // the last (most significant) character changed -> reset the iterator
            real_it = 0;
        }

        if (KMER::compare_suffix(kmer, get_first(kmers[i - 1])))
            continue; // i would generate the same dummy-1 k-mer as i-1, skip

        KMER prev_kmer = kmer;
        prev_kmer.to_prev(k + 1, BOSS::kSentinelCode);

        while (real_it < end && KMER::less(get_first(kmers[real_it]), prev_kmer, 1)) {
            real_it++;
        }
        bool is_redundant = false;
        while (real_it < end
                && KMER::compare_suffix(get_first(kmers[real_it]), prev_kmer, 1)) {
            if (get_first(kmers[real_it])[0] == prev_kmer[0]) {  // edge labels match
                is_redundant = true;
                break;
            }
            real_it++;
        }
        if (!is_redundant) {
            push_back(kmers, prev_kmer);
        }
    }
}

template <typename T>
inline T
reverse_complement(size_t k, const T &v, const std::vector<TAlphabet> &complement_code) {
    using KMER = utils::get_first_type_t<T>;
    using INT = typename KMER::WordType;
    INT kmer = utils::get_first(v).data();
    constexpr uint64_t mask = KMER::kFirstCharMask;
    auto cc = [&complement_code](uint64_t v) { return complement_code[v]; };
    INT last_two_chars = cc(kmer & mask);
    kmer >>= KMER::kBitsPerChar;
    last_two_chars = (last_two_chars << KMER::kBitsPerChar) | cc(kmer & mask);
    kmer >>= KMER::kBitsPerChar;
    INT result = 0;
    for (uint32_t i = 2; i < k; ++i) {
        TAlphabet next_char = kmer & mask;
        assert(next_char >= 0 && next_char < complement_code.size());
        result = (result << KMER::kBitsPerChar) | cc(next_char);
        kmer >>= KMER::kBitsPerChar;
    }
    result = (result << 2 * KMER::kBitsPerChar) | last_two_chars;
    if constexpr (utils::is_pair_v<T>) {
        return T(KMER(result), v.second);
    } else {
        return KMER(result);
    }
}

template <typename T>
void add_reverse_complements(size_t k, size_t num_threads, Vector<T> *kmers) {
    logger->trace("Adding reverse complements...");
    // reserve space for the reverse complements and the ~10% dummy source/sink k-mers
    size_t size = kmers->size();
    if (size < 2) {
        return;
    }
    kmers->reserve(2 * size);

    T *kmer = kmers->data() + 1; // skip $$...$
    const T *end = kmers->data() + size;
    for (; kmer != end; ++kmer) {
        const T &rc = reverse_complement(k + 1, *kmer, KmerExtractorBOSS::kComplementCode);
        if (get_first(rc) != get_first(*kmer)) {
            kmers->push_back(std::move(rc));
        } else {
            if constexpr (utils::is_pair_v<T>) {
                kmer->second *= 2;
            }
        }
    }
    ips4o::parallel::sort(kmers->begin() + size, kmers->end(), utils::LessFirst(), num_threads);
    size_t block_size = static_cast<size_t>(std::sqrt(size));
    common::merge(kmers->data(), 0, size, kmers->size(), block_size);
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
/**
 * Adds dummy nodes for the given kmers.
 * The method first adds dummy sink kmers, then dummy sources with sentinels of length 1
 * (aka dummy-1 sources). The method will then gradually add dummy sources with sentinels
 * of length 2, 3, ... up to k-1.
 *
 * @param k the node length in the BOOS graph (so k-mer length is k+1)
 * @param num_threads number of threads available for sorting
 * @param both_strands for each k-mer, also add its reverse complement
 * @tparam T the type of kmers being processed, typically either a KMer64/128/256 or an
 * std::pair<KMer64/128/256, int8/16/32> if counting kmers.
 */
template <typename T>
void recover_dummy_nodes(size_t k, bool both_strands, size_t num_threads, Vector<T> *kmers_p) {
    using KMER = get_first_type_t<T>;
    Vector<T> &kmers = *kmers_p;

    if (both_strands) {
        add_reverse_complements(k, num_threads, &kmers);
    }
    size_t original_end = kmers.size();
    logger->trace("Total number of real k-mers: {}", original_end);
    add_dummy_sink_kmers(k, &kmers);

    size_t dummy_source_begin = kmers.size();
    add_dummy_source_kmers(k, &kmers, original_end);

    ips4o::parallel::sort(kmers.begin() + dummy_source_begin, kmers.end(),
                          utils::LessFirst(), num_threads);

    logger->trace("Number of dummy k-mers with dummy prefix of length 1: {}",
                  kmers.size() - dummy_source_begin);

    for (size_t c = 2; c < k + 1; ++c) {
        size_t dummy_end = kmers.size();

        for (size_t i = dummy_source_begin; i < dummy_end; ++i) {
            KMER kmer = get_first(kmers[i]);
            if (KMER::compare_suffix(kmer, get_first(kmers[i - 1])))
                continue; // i would generate the same dummy source k-mer as i-1, skip

            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            push_back(kmers, kmer);
        }
        dummy_source_begin = dummy_end;
        ips4o::parallel::sort(kmers.begin() + dummy_source_begin, kmers.end(),
                              utils::LessFirst(), num_threads);

        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}", c,
                      kmers.size() - dummy_source_begin);
    }

    ips4o::parallel::sort(kmers.begin(), kmers.end(),
                          utils::LessFirst(), num_threads);
}

template <typename T>
using Encoder = common::EliasFanoEncoderBuffered<T>;
template <typename T>
using Decoder = common::EliasFanoDecoder<T>;

/**
 * Splits #kmers by W (kmer[0]) and F (kmer[k]) into |ALPHABET\{$}|^2 chunks.
 */
template <typename T>
std::vector<std::string>
split(size_t k, const std::filesystem::path &dir, const ChunkedWaitQueue<T> &kmers) {
    using T_INT = get_int_t<T>;

    const uint8_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    size_t chunk_count = std::pow(alphabet_size - 1, 2);

    logger->trace("Splitting k-mers into {} chunks...", chunk_count);

    std::vector<Encoder<T_INT>> sinks;
    std::vector<std::string> names(chunk_count);
    for (size_t i = 0; i < names.size(); ++i) {
        names[i] = dir/("original_split_by_F_W_" + std::to_string(i));
        sinks.emplace_back(names[i], ENCODER_BUFFER_SIZE);
    }

    size_t num_parent_kmers = 0;
    for (auto &it = kmers.begin(); it != kmers.end(); ++it) {
        const T &kmer = *it;
        TAlphabet F = get_first(kmer)[k];
        TAlphabet W = get_first(kmer)[0];

        if (W == 0) {
            // erase the $$$...$ k-mer that was added at the beginning
            // it will be added as a source dummy k-mer later
            assert(F == 0);
            continue;
        }
        assert(F != 0);

        // subtract 1 -- the sentinel with code 0
        size_t idx = (F - 1) * (alphabet_size - 1) + (W - 1);
        sinks[idx].add(reinterpret_cast<const T_INT&>(kmer));
        num_parent_kmers++;
    }
    std::for_each(sinks.begin(), sinks.end(), [](auto &f) { f.finish(); });
    logger->trace("Total number of non-dummy k-mers: {}", num_parent_kmers);
    return names;
}

template <typename Decoder, typename KMER>
void skip_same_suffix(const KMER &el, Decoder &decoder, size_t suf) {
    while (!decoder.empty()) {
        KMER kmer(decoder.top());
        if (!KMER::compare_suffix(kmer, el, suf)) {
            break;
        }
        decoder.pop();
    }
}

/**
 * Generates non-redundant dummy-1 source k-mers and dummy sink kmers from #kmers.
 * @return a triplet containing the names of the original k-mer blocks, the dummy-1 source
 * k-mer blocks and the dummy sink k-mers
 */
template <typename T>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::string>
generate_dummy_1_kmers(size_t k,
                       size_t num_threads,
                       const std::filesystem::path &dir,
                       ChunkedWaitQueue<T> *kmers) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer

    // for a DNA alphabet, this will contain 16 chunks, split by kmer[0] and kmer[1]
    std::vector<std::string> original_split_by_F_W = split(k, dir, *kmers);

    const uint8_t alphabet_size = KmerExtractorBOSS::alphabet.size();

    std::vector<Encoder<KMER_INT>> dummy_l1_chunks;
    std::vector<Encoder<KMER_INT>> dummy_sink_chunks;
    std::vector<std::string> dummy_l1_names(alphabet_size);
    std::vector<std::string> dummy_sink_names(alphabet_size);
    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        dummy_l1_names[i] = dir/("dummy_source_1_" + std::to_string(i));
        dummy_sink_names[i] = dir/("dummy_sink_" + std::to_string(i));
        dummy_l1_chunks.emplace_back(dummy_l1_names[i], ENCODER_BUFFER_SIZE);
        dummy_sink_chunks.emplace_back(dummy_sink_names[i], ENCODER_BUFFER_SIZE);
    }

    logger->trace("Generating dummy-1 source k-mers and dummy sink k-mers...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
    for (TAlphabet F = 1; F < alphabet_size; ++F) {  // skip $$..$
        std::vector<std::string> F_chunks(  // chunks with k-mers ***F*
                original_split_by_F_W.begin() + (F - 1) * (alphabet_size - 1),
                original_split_by_F_W.begin() + F * (alphabet_size - 1));
        common::MergeDecoder<KMER_INT> it(F_chunks, false);

        std::vector<std::string> W_chunks;  // chunks with k-mers of the form ****F
        for (TAlphabet c = 1; c < alphabet_size; ++c) {
            W_chunks.push_back(original_split_by_F_W[(c - 1) * (alphabet_size - 1) + (F - 1)]);
        }
        common::ConcatDecoder<KMER_INT> sink_gen_it(W_chunks);
        while (!it.empty()) {
            KMER dummy_source(it.pop());
            // skip k-mers that would generate identical source dummy k-mers
            skip_same_suffix(dummy_source, it, 0);
            dummy_source.to_prev(k + 1, BOSS::kSentinelCode);
            // generate dummy sink k-mers from all non-dummy kmers smaller than |dummy_source|
            while (!sink_gen_it.empty()
                    && sink_gen_it.top() < dummy_source.data()) {
                KMER v(sink_gen_it.pop());
                // skip k-mers with the same suffix as v, as they generate identical dummy
                // sink k-mers
                skip_same_suffix(v, sink_gen_it, 1);
                v.to_next(k + 1, BOSS::kSentinelCode);
                dummy_sink_chunks[F].add(v.data());
            }
            if (!sink_gen_it.empty()) {
                KMER top(sink_gen_it.top());
                if (KMER::compare_suffix(top, dummy_source, 1)) {
                    // The source dummy k-mer #dummy_source generated from #it is
                    // redundant iff it shares its suffix with another real k-mer (#top).
                    // In this case, #top generates a dummy sink k-mer redundant with #it.
                    // So if #dummy_source is redundant, the sink generated from #top is
                    // also redundant - so it's being skipped
                    skip_same_suffix(top, sink_gen_it, 1);
                    continue;
                }
            }
            dummy_l1_chunks[F].add(dummy_source.data());
        }
        // handle leftover sink_gen_it
        while (!sink_gen_it.empty()) {
            KMER v(sink_gen_it.pop());
            skip_same_suffix(v, sink_gen_it, 1);
            v.to_next(k + 1, BOSS::kSentinelCode);
            dummy_sink_chunks[F].add(v.data());
        }
    }

    for (TAlphabet i = 0; i < alphabet_size; ++i) {
        dummy_sink_chunks[i].finish();
        dummy_l1_chunks[i].finish();
    }

    // dummy sink k-mers are partitioned into blocks by F (kmer[1]), so simply
    // concatenating the blocks will result in a single ordered block
    std::string dummy_sink_name = dir/"dummy_sink";
    common::concat(dummy_sink_names, dummy_sink_name);

    // similarly, the 16 blocks of the original k-mers can be concatenated in groups of
    // 4 without destroying the order
    std::vector<std::string> real_split_by_W;
    for (TAlphabet W = 1; W < alphabet_size; ++W) {
        std::vector<std::string> blocks;
        for (TAlphabet F = 1; F < alphabet_size; ++F) {
            blocks.push_back(original_split_by_F_W[(F - 1) * (alphabet_size - 1) + W - 1]);
        }
        real_split_by_W.push_back(dir/("real_split_by_W_" + std::to_string(W)));
        common::concat(blocks, real_split_by_W.back());
    }
    return { real_split_by_W, dummy_l1_names, dummy_sink_name };
}

/** Merges #original_kmers with #reverse_complements and places the result into #kmers */
template <typename T, typename T_INT>
static void merge(common::EliasFanoDecoder<T_INT> &original_kmers,
                  ChunkedWaitQueue<T_INT> &reverse_complements,
                  ChunkedWaitQueue<T> *kmers) {
    auto &kmers_int = reinterpret_cast<ChunkedWaitQueue<T_INT> &>(*kmers);
    auto &it = reverse_complements.begin();
    std::optional<T_INT> orig = original_kmers.next();
    while (it != reverse_complements.end() && orig.has_value()) {
        if (get_first(orig.value()) < get_first(*it)) {
            kmers_int.push(orig.value());
            orig = original_kmers.next();
        } else {
            kmers_int.push(*it);
            ++it;
        }
    }
    while (it != reverse_complements.end()) {
        kmers_int.push(*it);
        ++it;
    }
    while (orig.has_value()) {
        kmers_int.push(orig.value());
        orig = original_kmers.next();
    }
    kmers->shutdown();
}

/**
 * Adds reverse complements
 */
template <typename T, typename T_INT>
void add_reverse_complements(size_t k,
                             const std::filesystem::path &dir,
                             ThreadPool& async_worker,
                             common::SortedSetDisk<T_INT> *rc_set,
                             ChunkedWaitQueue<T> *kmers) {
    logger->trace("Adding reverse complements...");
    common::EliasFanoEncoderBuffered<T_INT> original(dir/"original", ENCODER_BUFFER_SIZE);
    for (auto &it = ++kmers->begin(); it != kmers->end(); ++it) {
        const T &kmer = *it;
        const T &reverse = reverse_complement(k + 1, *it, KmerExtractorBOSS::kComplementCode);
        if (get_first(kmer) != get_first(reverse)) {
            rc_set->add(reinterpret_cast<const T_INT &>(reverse));
            original.add(reinterpret_cast<const T_INT &>(kmer));
        } else {
            if constexpr (utils::is_pair_v<T>) {
                original.add({kmer.first.data(), 2 * kmer.second });
            } else {
                original.add(reinterpret_cast<const T_INT &>(kmer));
            }
        }
    }
    original.finish();
    // start merging #original with #reverse_complements into #kmers
    kmers->reset();
    async_worker.enqueue([rc_set, &dir, kmers](){
      ChunkedWaitQueue<T_INT> &reverse_complements = rc_set->data(true);
      common::EliasFanoDecoder<T_INT> original_kmers(dir/"original");
      merge(original_kmers, reverse_complements, kmers);
    });
}

/**
 * Specialization of recover_dummy_nodes for a disk-based container, such as
 * #SortedSetDisk and #SortedMultisetDisk.
 * The method first traverses the original kmers and generates dummy source k-mers of
 * prefix length 1 (dummy-1). At the next step, the non-redundant dummy-1 kmers are merged
 * with the original k-mers.
 * The method then gradually constructs dummy-i k-mers for i=2..k and writes them into
 * separate files, de-duped and sorted.
 * The final result is obtained by merging the original #kmers and dummy-1 kmers with
 * the dummy-k kmers, for k=2..k
 */
template <class KmerCollector, typename T>
void recover_dummy_nodes_disk(const KmerCollector &kmer_collector,
                              ChunkedWaitQueue<T> *kmers,
                              ThreadPool &async_worker) {
    using KMER = get_first_type_t<T>; // 64/128/256-bit KmerBOSS
    using T_INT = get_int_t<T>; // either KMER_INT or <KMER_INT, count>
    using KMER_INT = typename KMER::WordType; // 64/128/256-bit integer

    const uint8_t alphabet_size = KmerExtractorBOSS::alphabet.size();


    size_t k = kmer_collector.get_k() - 1;
    const std::filesystem::path dir = kmer_collector.tmp_dir();
    size_t num_threads = kmer_collector.num_threads();

    // compute the reverse complements of #kmers and place them in #rc_set
    std::string rc_dir = dir/"rc";
    std::filesystem::create_directory(rc_dir);
    common::SortedSetDisk<T_INT> rc_set(num_threads, kmer_collector.buffer_size(), rc_dir,
                                        std::numeric_limits<size_t>::max());
    if (kmer_collector.is_both_strands_mode()) {
        add_reverse_complements(k, dir, async_worker, &rc_set, kmers);
    }

    std::string dummy_sink_name;
    std::vector<std::string> real_split_by_W;
    std::vector<std::string> dummy_names;
    std::tie(real_split_by_W, dummy_names, dummy_sink_name)
            = generate_dummy_1_kmers(k, num_threads, dir, kmers);

    // stores the sorted original kmers and dummy-1 k-mers
    std::vector<std::string> dummy_chunks = { dummy_sink_name };
    std::vector<std::string> dummy_next_names(alphabet_size);
    // generate dummy k-mers of prefix length 1..k
    logger->trace("Starting generating dummy-1..k source k-mers...");
    for (size_t dummy_pref_len = 1; dummy_pref_len <= k; ++dummy_pref_len) {
        // this will compress all sorted dummy k-mers of given prefix length
        for (const std::string &f : dummy_names) {
            dummy_chunks.push_back(f);
        }

        std::vector<Encoder<KMER_INT>> dummy_next_chunks;
        for (TAlphabet i = 0; i < alphabet_size; ++i) {
            dummy_next_names[i] = dir/("dummy_source_"
                    + std::to_string(dummy_pref_len + 1) + "_" + std::to_string(i));
            dummy_next_chunks.emplace_back(dummy_next_names[i], ENCODER_BUFFER_SIZE);
        }

        KMER prev_kmer(0);
        uint64_t num_kmers = 0;
        const std::function<void(const KMER_INT &)> &write_dummy = [&](const KMER_INT &v) {
            KMER kmer(v);
            kmer.to_prev(k + 1, BOSS::kSentinelCode);
            if (prev_kmer != kmer) {
                dummy_next_chunks[kmer[0]].add(kmer.data());
                prev_kmer = std::move(kmer);
            }
            num_kmers++;
        };
        common::merge_files(dummy_names, write_dummy, false);

        std::for_each(dummy_next_chunks.begin(), dummy_next_chunks.end(),
                      [](auto &v) { v.finish(); });
        std::swap(dummy_names, dummy_next_names);
        logger->trace("Number of dummy k-mers with dummy prefix of length {}: {}",
                      dummy_pref_len, num_kmers);
    }
    // remove the last chunks with .up and .count
    const std::function<void(const KMER_INT &)> on_merge = [](const KMER_INT &) {};
    common::merge_files(dummy_names, on_merge);

    // at this point, we have the original k-mers and dummy-1 k-mers in original_and_dummy_l1,
    // the dummy-x k-mers in dummy_source_{x}, and we merge them all into a single stream
    kmers->reset();
    // add the main dummy source k-mer
    if constexpr (utils::is_pair_v<T>) {
        kmers->push({KMER(0), 0});
    } else {
        kmers->push(KMER(0));
    }
    // push all other dummy and non-dummy k-mers to |kmers|
    async_worker.enqueue([kmers, real_split_by_W, dummy_chunks]() {
        std::function<void(const T_INT &)> on_new_item
                = [kmers](const T_INT &v) { kmers->push(reinterpret_cast<const T &>(v)); };
        common::merge_dummy(real_split_by_W, dummy_chunks, on_new_item);
        kmers->shutdown();
    });
}

inline std::vector<TAlphabet>
encode_filter_suffix_boss(const std::string &filter_suffix) {
    KmerExtractorBOSS kmer_extractor;
    std::vector<TAlphabet> filter_suffix_encoded;
    for (char c : filter_suffix) {
        filter_suffix_encoded.push_back(c == BOSS::kSentinel
                                        ? BOSS::kSentinelCode
                                        : kmer_extractor.encode(c));
    }
    return filter_suffix_encoded;
}

template <typename T_TO, typename T_FROM>
ChunkedWaitQueue<T_TO>& reinterpret_container(ChunkedWaitQueue<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    return reinterpret_cast<ChunkedWaitQueue<T_TO>&>(container);
}

template <typename T_TO, typename T_FROM>
Vector<T_TO>& reinterpret_container(Vector<T_FROM> &container) {
    static_assert(sizeof(T_TO) == sizeof(T_FROM));
    return reinterpret_cast<Vector<T_TO>&>(container);
}


template <typename KmerCollector>
class BOSSChunkConstructor : public IBOSSChunkConstructor {
    friend IBOSSChunkConstructor;

    template <template <typename KMER> class KmerContainer, typename... Args>
    friend std::unique_ptr<IBOSSChunkConstructor>
    initialize_boss_chunk_constructor(size_t k, const Args& ...args);

  private:
    BOSSChunkConstructor(size_t k,
                         bool both_strands_mode,
                         uint8_t bits_per_count,
                         const std::string &filter_suffix,
                         size_t num_threads,
                         double memory_preallocated,
                         const std::filesystem::path &tmp_dir,
                         size_t max_disk_space)
        : kmer_collector_(k + 1,
                          both_strands_mode,
                          encode_filter_suffix_boss(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          tmp_dir,
                          max_disk_space,
                          filter_suffix.empty() /* collect only canonical */),
          bits_per_count_(bits_per_count) {
        if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
            kmer_collector_.add_kmer(std::vector<TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_sequence(std::string_view sequence, uint64_t count) {
        kmer_collector_.add_sequence(sequence, count);
    }

    void add_sequences(const std::function<void(CallString)> &generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk() {
        typename KmerCollector::Data &kmer_ints = kmer_collector_.data();

        using T_INT = typename KmerCollector::Data::value_type;
        using T = get_kmer_t<typename KmerCollector::Kmer, T_INT>;
        auto &kmers = reinterpret_container<T>(kmer_ints);

        if (!kmer_collector_.suffix_length()) {
            Timer timer;
            if constexpr ((std::is_same_v<typename KmerCollector::Data,
                                          ChunkedWaitQueue<T_INT>>)) {
                recover_dummy_nodes_disk(kmer_collector_, &kmers, async_worker_);
            } else {
                // kmer_collector stores (BOSS::k_ + 1)-mers
                static_assert(std::is_same_v<typename KmerCollector::Data, Vector<T_INT>>);
                recover_dummy_nodes(kmer_collector_.get_k() - 1,
                                    kmer_collector_.is_both_strands_mode(),
                                    kmer_collector_.num_threads(), &kmers);
            }
            logger->trace("Dummy source k-mers were reconstructed in {} sec",
                          timer.elapsed());
        }

        BOSS::Chunk *result;

        // kmer_collector stores (BOSS::k_ + 1)-mers
        result = new BOSS::Chunk(kmer_collector_.alphabet_size(), kmer_collector_.get_k() - 1,
                                 kmer_collector_.is_both_strands_mode(), kmers,
                                 bits_per_count_);

        kmer_collector_.clear();

        return result;
    }

    uint64_t get_k() const { return kmer_collector_.get_k() - 1; }

    KmerCollector kmer_collector_;
    uint8_t bits_per_count_;
    /** Used as an async executor for merging chunks from disk */
    ThreadPool async_worker_ = ThreadPool(1, 1);
};

template <template <typename KMER> class KmerContainer, typename... Args>
static std::unique_ptr<IBOSSChunkConstructor>
initialize_boss_chunk_constructor(size_t k, const Args& ...args) {
    if (k < 1 || k > 256 / KmerExtractorBOSS::bits_per_char - 1) {
        logger->error("For succinct graph, k must be between 2 and {}",
                      256 / KmerExtractorBOSS::bits_per_char - 1);
        exit(1);
    }

    if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 64) {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer64>>(k, args...)
        );
    } else if ((k + 1) * KmerExtractorBOSS::bits_per_char <= 128) {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer128>>(k, args...)
        );
    } else {
        return std::unique_ptr<IBOSSChunkConstructor>(
            new BOSSChunkConstructor<KmerContainer<KmerExtractorBOSS::Kmer256>>(k, args...)
        );
    }
}

template <typename KMER>
using KmerSetVector
    = KmerCollector<KMER, KmerExtractorBOSS, common::SortedSet<typename KMER::WordType>>;

template <typename KMER>
using KmerMultsetVector8
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint8_t>>;

template <typename KMER>
using KmerMultsetVector16
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint16_t>>;

template <typename KMER>
using KmerMultsetVector32
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultiset<typename KMER::WordType, uint32_t>>;

template <typename KMER>
using KmerSetDisk
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedSetDisk<typename KMER::WordType>>;

template <typename KMER>
using KmerMultsetDiskVector8
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint8_t>>;

template <typename KMER>
using KmerMultsetDiskVector16
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint16_t>>;

template <typename KMER>
using KmerMultsetDiskVector32
    = KmerCollector<KMER,
                    KmerExtractorBOSS,
                    common::SortedMultisetDisk<typename KMER::WordType, uint32_t>>;

std::unique_ptr<IBOSSChunkConstructor>
IBOSSChunkConstructor::initialize(size_t k,
                                  bool canonical_mode,
                                  uint8_t bits_per_count,
                                  const std::string &filter_suffix,
                                  size_t num_threads,
                                  double memory_preallocated,
                                  kmer::ContainerType container_type,
                                  const std::filesystem::path &tmp_dir,
                                  size_t max_disk_space_bytes) {
#define OTHER_ARGS k, canonical_mode, bits_per_count, filter_suffix, \
                   num_threads, memory_preallocated, tmp_dir, max_disk_space_bytes

    switch (container_type) {
        case kmer::ContainerType::VECTOR:
            if (!bits_per_count) {
                return initialize_boss_chunk_constructor<KmerSetVector>(OTHER_ARGS);
            } else if (bits_per_count <= 8) {
                return initialize_boss_chunk_constructor<KmerMultsetVector8>(OTHER_ARGS);
            } else if (bits_per_count <= 16) {
                return initialize_boss_chunk_constructor<KmerMultsetVector16>(OTHER_ARGS);
            } else if (bits_per_count <= 32) {
                return initialize_boss_chunk_constructor<KmerMultsetVector32>(OTHER_ARGS);
            } else {
                throw std::runtime_error(
                        "Error: trying to allocate too many bits per k-mer count");
            }
        case kmer::ContainerType::VECTOR_DISK:
            if (!bits_per_count) {
                return initialize_boss_chunk_constructor<KmerSetDisk>(OTHER_ARGS);
            } else if (bits_per_count <= 8) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector8>(OTHER_ARGS);
            } else if (bits_per_count <= 16) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector16>(OTHER_ARGS);
            } else if (bits_per_count <= 32) {
                return initialize_boss_chunk_constructor<KmerMultsetDiskVector32>(OTHER_ARGS);
            } else {
                throw std::runtime_error(
                        "Error: trying to allocate too many bits per k-mer count");
            }
        default:
            logger->error("Invalid container type {}", container_type);
            std::exit(1);
    }
}

} // namespace succinct
} // namespace mtg
