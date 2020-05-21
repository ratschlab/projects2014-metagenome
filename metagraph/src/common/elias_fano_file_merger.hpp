#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "common/elias_fano.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"

namespace mg {
namespace common {

/**
 * Heap implemented as a sorted vector.
 *
 * The heap uses a vector as the underlying data structure, thus making it efficient
 * only for a small (<1000) number of elements.
 * @tparam T the actual heap element type; this is the type tested for equality
 */
// Note: profiling shows that using a sorted vector instead of a std::priority queue is
// faster up to ~1000 elements.  Using an unsorted vector (faster insert,
// slower pop()) is ~40% slower. Preventing duplicates in the heap so that we don't
// need to test for dupes at pop time is ~60% slower.
// Note: Profiling indicates that merging elements within the heap is ~30% slower than
// inserting duplicates and merging them when popping from the heap, as we do now
template <typename T, class Compare = std::greater<T>>
class MergeHeap {
    /** The heap stores pairs <Element, SourceIndex> */
    using value_type = std::pair<T, uint32_t>;

  public:
    void emplace(T el, uint32_t idx) {
        auto it = std::lower_bound(els.begin(), els.end(), el,
                                   [this](const value_type &p, const T &v) {
                                       return compare_(p.first, v);
                                   });
        els.emplace(it, el, idx);
    }

    const value_type& top() const { return els.back(); }

    value_type pop() {
        value_type result = els.back();
        els.pop_back();
        return result;
    }

    bool empty() { return els.empty(); }

  private:
    // elements stored in decreasing order of the first tuple member
    std::vector<value_type> els;
    Compare compare_ = Compare();
};

/**
 * Decoder that reads data from several sorted files and merges it into a single sorted
 * stream.
 * @tparam T the type of data being stored
 */
template <typename T>
class MergeDecoder {
  public:
    MergeDecoder(const std::vector<std::string> &source_names) {
        T data_item;
        sources_.resize(source_names.size());
        for (uint32_t i = 0; i < sources_.size(); ++i) {
            sources_[i] = std::ifstream(source_names[i], std::ios::binary);
            if (sources_[i].read(reinterpret_cast<char *>(&data_item), sizeof(data_item))) {
                heap_.emplace(data_item, i);
            }
        }
    }
    std::optional<T> top() {
        if (heap_.empty()) {
            return std::nullopt;
        }
        return heap_.top().first;
    }
    std::optional<T> next() {
        if (heap_.empty()) {
            return std::nullopt;
        }
        auto [result, chunk_index] = heap_.pop();
        T data_item;
        if (sources_[chunk_index].read(reinterpret_cast<char *>(&data_item), sizeof(T))) {
            heap_.emplace(data_item, chunk_index);
        }
        return result;
    }

  private:
    std::vector<std::ifstream> sources_;
    common::MergeHeap<T> heap_;
};

namespace internal {
/**
 * Given a list of n Decoders, containing ordered elements of type T, merge the n
 * sources into a single (ordered) list of type T and invoke #on_new_item for each element
 * @tparam T the type of the  elements to be merged (typically a 64/128 or 256-bit k-mer)
 * @param decoders sources containing sorted lists of type T
 * @param on_new_item callback to invoke when a new element was merged
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T, typename Decoder>
uint64_t merge_files(std::vector<Decoder> &decoders,
                     const std::function<void(const T &)> &on_new_item) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 0;

    MergeHeap<T> merge_heap;
    std::optional<T> data_item;

    for (uint32_t i = 0; i < decoders.size(); ++i) {
        data_item = decoders[i].next();
        if (data_item.has_value()) {
            merge_heap.emplace(data_item.value(), i);
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
        return num_elements_read;
    }

    T last_written = merge_heap.top().first;
    on_new_item(last_written);

    while (!merge_heap.empty()) {
        auto [smallest, chunk_index] = merge_heap.pop();

        if (smallest != last_written) {
            on_new_item(smallest);
            last_written = smallest;
        }

        if ((data_item = decoders[chunk_index].next()).has_value()) {
            merge_heap.emplace(data_item.value(), chunk_index);
            num_elements_read++;
        }
    }

    return num_elements_read;
}
} // internal

/**
 * Merges Elias-Fano sorted compressed files into a single stream.
 */
template <typename T>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    std::vector<EliasFanoDecoder<T>> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.emplace_back(sources[i], remove_sources);
    }

    return internal::merge_files(decoders, on_new_item);
}

template <typename T>
class Decoder {
  public:
    Decoder(const std::string &name) : source_(name, std::ios::binary) {}
    std::optional<T> next() {
        T result;
        if (source_.read(reinterpret_cast<char *>(&result), sizeof(T))) {
            return result;
        }
        return std::nullopt;
    }

  private:
    std::ifstream source_;
};

/** Merges binary files into a single stream */
template <typename T>
uint64_t merge_files_uncompressed(const std::vector<std::string> &sources,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    std::vector<Decoder<T>> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.push_back(Decoder<T>(sources[i]));
    }
    size_t num_elements_read = internal::merge_files(decoders, on_new_item);

    if (remove_sources) {
        std::for_each(sources.begin(), sources.end(),
                      [](const std::string &name) { std::filesystem::remove(name); });
    }

    return num_elements_read;
}

// TODO: these two `merge_files` are almost identical. Merge them into one.
//       Implement the merging mechanism  (remove duplicates, increment
//       counters) in the caller?
/**
 * Given a list of n source files, containing ordered pairs of  <element, count>,
 * merge the n sources (and the corresponding counts) into a single list, ordered by el
 * and delete the original files.
 * If two pairs have the same first element, the counts are added together.
 * @param sources the files containing sorted lists of pairs of type <T, C>
 * @param on_new_item callback to invoke when a new element was merged
 * @param remove_sources if true, remove source files after merging
 *
 * @return the total number of elements read from all files
 *
 * Note: this method blocks until all the data was successfully merged.
 */
template <typename T, typename C>
uint64_t merge_files(const std::vector<std::string> &sources,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 0;

    MergeHeap<std::pair<T, C>, utils::GreaterFirst> merge_heap;
    std::optional<std::pair<T, C>> data_item;
    std::vector<EliasFanoDecoder<std::pair<T, C>>> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.emplace_back(sources[i], remove_sources);
        data_item = decoders.back().next();
        if (data_item.has_value()) {
            merge_heap.emplace(data_item.value(), i);
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
        return num_elements_read;
    }

    // initialize the smallest element
    std::pair<T, C> current = { merge_heap.top().first.first, 0 };

    while (!merge_heap.empty()) {
        auto [smallest, chunk_index] = merge_heap.pop();

        if (smallest.first != current.first) {
            on_new_item(current);
            current = smallest;
        } else {
            if (current.second < std::numeric_limits<C>::max() - smallest.second) {
                current.second += smallest.second;
            } else {
                current.second = std::numeric_limits<C>::max();
            }
        }

        if ((data_item = decoders[chunk_index].next()).has_value()) {
            merge_heap.emplace(data_item.value(), chunk_index);
            num_elements_read++;
        }
    }
    on_new_item(current);

    return num_elements_read;
}

/**
 * Merges the Ts in #source with the Ts in #source_no_count. This is no different than
 * calling merge() for all files.
 * @param source name of a source file containing Elias-Fano encoded INTs
 * @param source_no_count name of a soruce file containing Elias-Fano encoded INTs corresponding to dummy k-mers
 * @param on_new_item callback to invoke for each merged item
 * @param remove_sources if true, the #source and #source_no_count files will be removed
 */

/**
 * Merges the Ts in #source with the Ts in #source_no_count. This is no different than
 * calling merge() for all files.
 * @param source name of a source file containing Elias-Fano encoded INTs
 * @param source_no_count name of a soruce file containing Elias-Fano encoded INTs corresponding to dummy k-mers
 * @param on_new_item callback to invoke for each merged item
 * @param remove_sources if true, the #source and #source_no_count files will be removed
 */
template <typename T>
uint64_t merge_dummy(const std::vector<std::string> &source,
                     std::vector<std::string> source_no_count,
                     const std::function<void(const T &)> &on_new_item,
                     bool remove_sources = true) {
    source_no_count.insert(source_no_count.end(), source.begin(), source.end());
    return merge_files_uncompressed(source_no_count, on_new_item, remove_sources);
}

/**
 * Merges the <T, C> pairs in #sources with the Ts in #source_no_count. The INTs in
 * source_no_count will be assigned a count of 0.
 */
template <typename T, typename C>
uint64_t merge_dummy(const std::vector<std::string> &sources,
                     const std::vector<std::string> &sources_no_count,
                     const std::function<void(const std::pair<T, C> &)> &on_new_item,
                     bool remove_sources = true) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    uint64_t num_elements_read = 0;

    MergeHeap<std::pair<T, C>, utils::GreaterFirst> merge_heap;
    std::pair<T,C> data_item;
    std::vector<std::ifstream> decoders;
    for (uint32_t i = 0; i < sources.size(); ++i) {
        decoders.push_back(std::ifstream(sources[i], std::ios::binary));
        if (decoders.back().read(reinterpret_cast<char *>(&data_item), sizeof(data_item))) {
            merge_heap.emplace(data_item, i);
            num_elements_read++;
        }
    }
    T data_item2;
    for (uint32_t i = 0; i < sources_no_count.size(); ++i) {
        decoders.push_back(std::ifstream(sources_no_count[i], std::ios::binary));
        if (decoders.back().read(reinterpret_cast<char *>(&data_item2), sizeof(data_item2))) {
            merge_heap.emplace({data_item2, 0}, i + sources.size());
            num_elements_read++;
        }
    }

    if (merge_heap.empty()) {
        return 0;
    }

    // initialize the smallest element
    std::pair<T, C> current = { merge_heap.top().first.first, 0 };

    while (!merge_heap.empty()) {
        auto [smallest, chunk_index] = merge_heap.pop();

        if (smallest.first != current.first) {
            on_new_item(current);
            current = smallest;
        } else {
            if (current.second < std::numeric_limits<C>::max() - smallest.second) {
                current.second += smallest.second;
            } else {
                current.second = std::numeric_limits<C>::max();
            }
        }
        if (chunk_index < sources.size()) {
            if (decoders[chunk_index].read(reinterpret_cast<char *>(&data_item), sizeof(data_item))) {
                merge_heap.emplace(data_item, chunk_index);
                num_elements_read++;
            }
        } else {
            if (decoders[chunk_index].read(reinterpret_cast<char *>(&data_item2), sizeof(data_item2))) {
                merge_heap.emplace({ data_item2, 0 }, chunk_index);
                num_elements_read++;
            }
        }
    }
    on_new_item(current);

    if (remove_sources) {
        std::for_each(sources.begin(), sources.end(),
                      [](const std::string &name) { std::filesystem::remove(name); });
        std::for_each(sources_no_count.begin(), sources_no_count.end(),
                      [](const std::string &name) { std::filesystem::remove(name); });
    }

    return num_elements_read;
}


} // namespace common
} // namespace mg
