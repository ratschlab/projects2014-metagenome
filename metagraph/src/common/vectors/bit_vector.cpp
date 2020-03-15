#include "bit_vector.hpp"

#include <type_traits>
#include <cassert>

#include "vector_algorithm.hpp"
#include "bit_vector_adaptive.hpp"
#include "bit_vector_stat.hpp"
#include "bit_vector_sdsl.hpp"
#include "bit_vector_dyn.hpp"
#include "bit_vector_sd.hpp"


std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    return os << bv.to_vector();
}

template <class Vector>
Vector bit_vector::convert_to() {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<Vector*>(this)) {
        // to the same type, no conversion
        return dynamic_cast<Vector&&>(*this);

    } else if (dynamic_cast<bit_vector_stat*>(this)) {
        // stat -> anything else
        return Vector(std::move(dynamic_cast<bit_vector_stat*>(this)->vector_));

    } else if (dynamic_cast<bit_vector_adaptive*>(this)) {
        // adaptive(x) -> anything
        return dynamic_cast<bit_vector_adaptive*>(this)->vector_->convert_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
    }
}
template bit_vector_dyn bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat bit_vector::convert_to<bit_vector_stat>();
template bit_vector_sd bit_vector::convert_to<bit_vector_sd>();
template bit_vector_il<> bit_vector::convert_to<bit_vector_il<>>();
template bit_vector_hyb<> bit_vector::convert_to<bit_vector_hyb<>>();
template bit_vector_rrr<3> bit_vector::convert_to<bit_vector_rrr<3>>();
template bit_vector_rrr<8> bit_vector::convert_to<bit_vector_rrr<8>>();
template bit_vector_rrr<15> bit_vector::convert_to<bit_vector_rrr<15>>();
template bit_vector_rrr<31> bit_vector::convert_to<bit_vector_rrr<31>>();
template bit_vector_rrr<63> bit_vector::convert_to<bit_vector_rrr<63>>();
template bit_vector_rrr<127> bit_vector::convert_to<bit_vector_rrr<127>>();
template bit_vector_rrr<255> bit_vector::convert_to<bit_vector_rrr<255>>();
template sdsl::bit_vector bit_vector::convert_to<sdsl::bit_vector>();
template<> bit_vector_small bit_vector::convert_to() {
    return bit_vector_small(std::move(*this));
}
template<> bit_vector_smart bit_vector::convert_to() {
    return bit_vector_smart(std::move(*this));
}

template <class Vector>
Vector bit_vector::copy_to() const {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<const Vector*>(this)) {
        // copy to the same type, no conversion
        return Vector(dynamic_cast<const Vector&>(*this));

    } else if (dynamic_cast<const bit_vector_stat*>(this)) {
        // copy stat -> anything else
        auto bv = dynamic_cast<const bit_vector_stat*>(this)->vector_;
        return Vector(std::move(bv));

    } else if (dynamic_cast<const bit_vector_adaptive*>(this)) {
        // copy adaptive(x) -> anything
        return dynamic_cast<const bit_vector_adaptive*>(this)->vector_->copy_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
    }
}
template bit_vector_dyn bit_vector::copy_to<bit_vector_dyn>() const;
template bit_vector_stat bit_vector::copy_to<bit_vector_stat>() const;
template bit_vector_sd bit_vector::copy_to<bit_vector_sd>() const;
template bit_vector_il<> bit_vector::copy_to<bit_vector_il<>>() const;
template bit_vector_hyb<> bit_vector::copy_to<bit_vector_hyb<>>() const;
template bit_vector_rrr<3> bit_vector::copy_to<bit_vector_rrr<3>>() const;
template bit_vector_rrr<8> bit_vector::copy_to<bit_vector_rrr<8>>() const;
template bit_vector_rrr<15> bit_vector::copy_to<bit_vector_rrr<15>>() const;
template bit_vector_rrr<31> bit_vector::copy_to<bit_vector_rrr<31>>() const;
template bit_vector_rrr<63> bit_vector::copy_to<bit_vector_rrr<63>>() const;
template bit_vector_rrr<127> bit_vector::copy_to<bit_vector_rrr<127>>() const;
template bit_vector_rrr<255> bit_vector::copy_to<bit_vector_rrr<255>>() const;
template sdsl::bit_vector bit_vector::copy_to<sdsl::bit_vector>() const;
template<> bit_vector_small bit_vector::copy_to() const {
    return bit_vector_small(*this);
}
template<> bit_vector_smart bit_vector::copy_to() const {
    return bit_vector_smart(*this);
}

void bit_vector::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());

    // TODO: tune the coefficient for each representation
    if (num_set_bits() * 3 < size()) {
        // for sparse vectors
        call_ones([other](auto i) { (*other)[i] = true; });

    } else {
        uint64_t i;
        const uint64_t end = size();
        uint64_t *data = other->data();
        for (i = 0; i + 64 <= end; i += 64, ++data) {
            *data |= get_int(i, 64);
        }
        if (i < size()) {
            *data |= get_int(i, size() - i);
        }
    }
}


// indexes are distinct and sorted
sdsl::bit_vector subvector(const bit_vector &col,
                           const std::vector<uint64_t> &indexes) {
    assert(indexes.size() <= col.size());

    sdsl::bit_vector shrinked(indexes.size(), 0);

    uint64_t max_rank = col.num_set_bits();
    if (!max_rank)
        return shrinked;

    // the case of uncompressed vector
    if (dynamic_cast<const bit_vector_stat *>(&col)) {
        for (size_t j = 0; j < indexes.size(); ++j) {
            if (col[indexes[j]])
                shrinked[j] = true;
        }
        return shrinked;
    }

    uint64_t cur_rank = 1;
    uint64_t next_pos = col.select1(1);

    for (size_t j = 0; j < indexes.size(); ++j) {
        if (indexes[j] < next_pos)
            continue;

        if (indexes[j] == next_pos) {
            shrinked[j] = true;
            continue;
        }

        // indexes[j] > next_pos
        if (col[indexes[j]]) {
            shrinked[j] = true;
            continue;
        }

        // we found a zero, update next_pos
        cur_rank = col.rank1(indexes[j]) + 1;
        if (cur_rank > max_rank)
            return shrinked;

        next_pos = col.select1(cur_rank);
    }

    return shrinked;
}
