#ifndef __INT_VECTOR_HPP__
#define __INT_VECTOR_HPP__

#include <functional>

#include <sdsl/int_vector.hpp>

void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   std::function<void(uint64_t, uint64_t)> callback);

void call_nonzeros(const sdsl::int_vector<> &vector,
                   std::function<void(uint64_t, uint64_t)> callback);

#endif // __INT_VECTOR_HPP__
