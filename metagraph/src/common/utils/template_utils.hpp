#ifndef __TEMPLATE_UTILS_HPP__
#define __TEMPLATE_UTILS_HPP__

#include <type_traits>
#include <utility>


namespace utils {

template <typename T>
struct is_pair : std::false_type {};

template <typename T1, typename T2>
struct is_pair<std::pair<T1,T2>> : std::true_type {};

static_assert(is_pair<std::pair<int,size_t>>::value);
static_assert(is_pair<std::pair<uint64_t,size_t>>::value);

static_assert(!is_pair<std::tuple<int,size_t>>::value);
static_assert(!is_pair<int>::value);

template <typename, template <typename, typename...> class>
struct is_instance : public std::false_type {};

template <typename... Ts, template <typename, typename...> class U>
struct is_instance<U<Ts...>, U> : public std::true_type {};

template <class T> struct dependent_false : std::false_type {};

template <typename T, typename... Us>
inline const T& get_first(const std::tuple<T, Us...> &tuple) { return std::get<0>(tuple); }

template <typename T, typename... Us>
inline T& get_first(std::tuple<T, Us...> &tuple) { return std::get<0>(tuple); }

template <typename T, typename U>
inline const T& get_first(const std::pair<T, U> &pair) { return pair.first; }

template <typename T, typename U>
inline T& get_first(std::pair<T, U> &pair) { return pair.first; }

template <typename T>
inline T& get_first(T &value) { return value; }

template <typename T>
inline const T& get_first(const T &value) { return value; }

template <typename T, typename U, typename... Us>
inline const U& get_second(const std::tuple<T, U, Us...> &tuple) { return std::get<1>(tuple); }

template <typename T, typename U, typename... Us>
inline U& get_second(std::tuple<T, U, Us...> &tuple) { return std::get<1>(tuple); }

template <typename T, typename U>
inline const U& get_second(const std::pair<T, U> &pair) { return pair.second; }

template <typename T, typename U>
inline U& get_second(std::pair<T, U> &pair) { return pair.second; }

// return the the type of the first element if T is a pair, or T otherwise
template <typename T, typename = void>
struct get_first_type {
    using type = T;
};

template <typename T>
struct get_first_type<T, std::void_t<typename T::first_type>> {
    using type = typename T::first_type;
};

class GreaterFirst {
  public:
    template <typename T>
    bool operator()(const T &p1, const T &p2) const {
        return get_first(p1) > get_first(p2);
    }
};

struct LessFirst {
    template <typename T>
    bool operator()(const T &p1, const T &p2) const {
        return get_first(p1) < get_first(p2);
    }
};

struct EqualFirst {
    template <typename T>
    bool operator()(const T &p1, const T &p2) const {
        return get_first(p1) == get_first(p2);
    }
};

struct LessSecond {
    template <typename T>
    bool operator()(const T &p1, const T &p2) const {
        return get_second(p1) < get_second(p2);
    }
};

} // namespace utils

#endif // __TEMPLATE_UTILS_HPP__
