#ifndef __DBG_BITMAP_CONSTRUCT_HPP__
#define __DBG_BITMAP_CONSTRUCT_HPP__

#include "dbg_bitmap.hpp"
#include "graph/representation/base/dbg_construct.hpp"

namespace mtg {
namespace bitmap_graph {


class IBitmapChunkConstructor : public IGraphChunkConstructor<DBGBitmap::Chunk> {
  public:
    virtual ~IBitmapChunkConstructor() {}

    static IBitmapChunkConstructor* initialize(size_t k,
                                               bool canonical_mode = false,
                                               uint8_t bits_per_count = 0,
                                               const std::string &filter_suffix = "",
                                               size_t num_threads = 1,
                                               double memory_preallocated = 0);

    virtual void add_sequence(std::string_view sequence, uint64_t count = 1) = 0;
    virtual void add_sequences(const std::function<void(CallString)> &generator) = 0;
    virtual void add_sequences(const std::function<void(CallStringCount)> &generator) = 0;

    virtual DBGBitmap::Chunk* build_chunk() = 0;

    virtual size_t get_k() const = 0;
    virtual bool is_canonical_mode() const = 0;

    virtual sdsl::int_vector<> get_weights(uint8_t bits_per_count = 8) = 0;
};


class DBGBitmapConstructor : public IGraphConstructor<DBGBitmap> {
  public:
    // Don't count k-mers if |bits_per_count| is zero.
    DBGBitmapConstructor(size_t k,
                         bool canonical_mode = false,
                         uint8_t bits_per_count = 0,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0);

    void add_sequence(std::string_view sequence, uint64_t count = 1) {
        constructor_->add_sequence(sequence, count);
    }

    void add_sequences(std::vector<std::string>&& sequences) {
        auto seqs = std::make_shared<std::vector<std::string>>(std::move(sequences));
        constructor_->add_sequences(
            [seqs](const CallString &callback) {
                std::for_each(seqs->begin(), seqs->end(), callback);
            }
        );
    }

    void add_sequences(const std::function<void(CallString)> &generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void build_graph(DBGBitmap *graph);
    DBGBitmap::Chunk* build_chunk() { return constructor_->build_chunk(); }

    uint64_t get_k() const { return constructor_->get_k(); }

    static DBGBitmap* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                              bool canonical_mode = false,
                                              bool verbose = false);

    static DBGBitmap* build_graph_from_chunks(uint64_t size,
                                              uint64_t num_kmers,
                                              const std::function<DBGBitmap::Chunk(void)> &next_chunk,
                                              bool canonical_mode = false);

  private:
    std::unique_ptr<IBitmapChunkConstructor> constructor_;
    uint8_t bits_per_count_;
};

} // namespace bitmap_graph
} // namespace mtg

#endif // __DBG_BITMAP_CONSTRUCT_HPP__
