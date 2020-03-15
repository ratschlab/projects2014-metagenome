#ifndef __BIT_VECTOR_SDSL_HPP__
#define __BIT_VECTOR_SDSL_HPP__

#include <cmath>
#include <cstdint>

#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/bit_vector_il.hpp>

#include "vector_algorithm.hpp"
#include "bit_vector.hpp"


template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
class bit_vector_sdsl : public bit_vector {
    template<typename>
    struct is_rrr : std::false_type {};
    template<uint16_t t_bs, class t_rac, uint16_t t_k>
    struct is_rrr<sdsl::rrr_vector<t_bs, t_rac, t_k>> : std::true_type {};

    template<typename>
    struct bv_traits {};

    template<uint32_t k_sblock_rate>
    struct bv_traits<sdsl::hyb_vector<k_sblock_rate>> {
        // hyb_vector doesn't support select
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;
    };
    template<uint32_t t_bs>
    struct bv_traits<sdsl::bit_vector_il<t_bs>> {
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;
    };
    template<uint16_t t_bs, class t_rac, uint16_t t_k>
    struct bv_traits<sdsl::rrr_vector<t_bs, t_rac, t_k>> {
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1;
        // sequential word access for bit_vector_rrr per bit is 21 times faster than select
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 21;
    };

  public:
    explicit bit_vector_sdsl(uint64_t size = 0, bool value = false)
      : bit_vector_sdsl(sdsl::bit_vector(size, value)) {}
    explicit bit_vector_sdsl(const sdsl::bit_vector &vector)
      : vector_(vector), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_) {}
    explicit bit_vector_sdsl(const bit_vector_sdsl &other) { *this = other; }

    bit_vector_sdsl(sdsl::bit_vector&& other) noexcept
      : vector_(std::move(other)), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_) {}
    bit_vector_sdsl(bit_vector_sdsl&& other) noexcept { *this = std::move(other); }
    bit_vector_sdsl(std::initializer_list<bool> init)
      : bit_vector_sdsl(sdsl::bit_vector(init)) {}

    inline bit_vector_sdsl& operator=(const bit_vector_sdsl &other);
    inline bit_vector_sdsl& operator=(bit_vector_sdsl&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override {
        return std::make_unique<bit_vector_sdsl>(*this);
    }

    inline uint64_t rank1(uint64_t id) const override;
    inline uint64_t select0(uint64_t id) const;
    inline uint64_t select1(uint64_t id) const override;
    inline std::pair<bool, uint64_t> inverse_select(uint64_t id) const override;
    inline uint64_t conditional_rank1(uint64_t id) const override;

    inline uint64_t next1(uint64_t id) const override;
    inline uint64_t prev1(uint64_t id) const override;

    inline bool operator[](uint64_t id) const override;
    inline uint64_t get_int(uint64_t id, uint32_t width) const override;

    inline bool load(std::istream &in) override;
    inline void serialize(std::ostream &out) const override;

    inline uint64_t size() const override { return vector_.size(); }

    inline sdsl::bit_vector to_vector() const override;

    inline void call_ones_in_range(uint64_t begin, uint64_t end,
                                   const VoidCall<uint64_t> &callback) const override;

    inline const bv_type& data() const { return vector_; }

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    static inline uint64_t predict_size(uint64_t size, uint64_t num_set_bits);

  private:
    static inline double logbinomial(uint64_t n, uint64_t m);
    static inline double entropy(double q);

    bv_type vector_;
    rank_1_type rk1_;
    select_1_type slct1_;
    select_0_type slct0_;
};


template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>&
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator=(const bit_vector_sdsl &other) {
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>&
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator=(bit_vector_sdsl&& other) noexcept {
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk1_(id >= this->size() ? this->size() : id + 1);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= size() - num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));
    return slct0_(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));
    return slct1_(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
std::pair<bool, uint64_t>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::inverse_select(uint64_t id) const {
    if constexpr(is_rrr<bv_type>{}) {
        if constexpr(bv_type::block_size != 15) {
            std::pair<bool, uint64_t> pair = vector_.inverse_select(id);
            pair.second += pair.first;
            return pair;
        }
    }
    // TODO: implement inverse_select for other bit vectors as well
    return bit_vector::inverse_select(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::conditional_rank1(uint64_t id) const {
    if constexpr(is_rrr<bv_type>{}) {
        if constexpr(bv_type::block_size != 15) {
            return vector_.conditional_rank(id);
        }
    }
    return bit_vector::conditional_rank1(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::next1(uint64_t id) const {
    assert(id < size());
    return ::next1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::prev1(uint64_t id) const {
    assert(id < size());
    return ::prev1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bool
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::get_int(uint64_t id, uint32_t width) const {
    assert(id + width <= size());
    return vector_.get_int(id, width);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bool
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load "
                  << typeid(bv_type).name() << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
void
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::serialize(std::ostream &out) const {
    vector_.serialize(out);

    if (!out.good())
        throw std::ofstream::failure("Error when dumping bit_vector_rrr");
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
sdsl::bit_vector
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::to_vector() const {
    if (num_set_bits()
            < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
        // sparse
        sdsl::bit_vector vector(size(), 0);
        uint64_t max_rank = size() ? rank1(size() - 1) : 0;
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct1_(i)] = 1;
        }
        return vector;

    } else if ((size() - num_set_bits())
            < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
        // dense
        sdsl::bit_vector vector(size(), 1);
        uint64_t max_rank = size() ? rank0(size() - 1) : 0;
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct0_(i)] = 0;
        }
        return vector;

    } else {
        // moderate density
        sdsl::bit_vector vector(size());

        uint64_t i;
        for (i = 0; i + 64 <= vector.size(); i += 64) {
            vector.set_int(i, vector_.get_int(i));
        }
        if (i < vector.size())
            vector.set_int(i, vector_.get_int(i, vector.size() - i), vector.size() - i);

        return vector;
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
void
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::call_ones_in_range(uint64_t begin, uint64_t end,
                     const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    if (num_set_bits()
            < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
        // sparse
        uint64_t num_ones = end ? rank1(end - 1) : 0;
        for (uint64_t r = begin ? rank1(begin - 1) + 1 : 1; r <= num_ones; ++r) {
            callback(select1(r));
        }
    } else if ((size() - num_set_bits())
                < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
        // dense
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = end ? rank0(end - 1) : 0;
        for (uint64_t r = begin ? rank0(begin - 1) + 1 : 1; r <= num_zeros; ++r) {
            zero_pos = select0(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < end) {
            callback(one_pos++);
        }
    } else {
        // moderate density
        assert(begin <= end);
        assert(end <= size());
        ::call_ones(vector_, begin, end, callback);
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::predict_size(uint64_t size, uint64_t num_set_bits) {
    if constexpr(is_rrr<bv_type>{}) {
        uint16_t block_size = bv_type::block_size;
        // TODO: correct the formula for block_size = 15
        return std::ceil(logbinomial(size, num_set_bits))
            + (size + block_size) / block_size * (sdsl::bits::hi(block_size) + 1);
    } else {
        throw std::runtime_error(std::string("Error: unknown space taken for this bit_vector")
                                    + typeid(bv_type).name());
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
double
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::logbinomial(uint64_t n, uint64_t m) {
    return (lgamma(n + 1)
                - lgamma(m + 1)
                - lgamma(n - m + 1)) / log(2);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
double
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::entropy(double q) {
    assert(q >= 0);
    assert(q <= 1);

    if (q == 0 || q == 1)
        return 0;

    return q * log2(q) + (1 - q) * log2(1 - q);
}


template <uint32_t k_sblock_rate = 16>
using bit_vector_hyb
    = bit_vector_sdsl<sdsl::hyb_vector<k_sblock_rate>,
                      typename sdsl::hyb_vector<k_sblock_rate>::rank_1_type,
                      typename sdsl::hyb_vector<k_sblock_rate>::select_1_type,
                      typename sdsl::hyb_vector<k_sblock_rate>::select_0_type>;

template <uint32_t block_size = 512>
using bit_vector_il
    = bit_vector_sdsl<sdsl::bit_vector_il<block_size>,
                      typename sdsl::bit_vector_il<block_size>::rank_1_type,
                      typename sdsl::bit_vector_il<block_size>::select_1_type,
                      typename sdsl::bit_vector_il<block_size>::select_0_type>;

template <uint32_t block_size = 63>
using bit_vector_rrr
    = bit_vector_sdsl<sdsl::rrr_vector<block_size>,
                      typename sdsl::rrr_vector<block_size>::rank_1_type,
                      typename sdsl::rrr_vector<block_size>::select_1_type,
                      typename sdsl::rrr_vector<block_size>::select_0_type>;

#endif // __BIT_VECTOR_SDSL_HPP__
