#pragma once

#include <tuple>

// -- Test if type T is in tuple Tuple ----------------------------------------

template <typename T, typename Tuple>
struct Is_in;

// recursion anchor
template <typename T>
struct Is_in<T, std::tuple<>> {
   static constexpr bool value = false;
};

template <typename T, typename T_head, typename... Ts>
struct Is_in<T, std::tuple<T_head, Ts...>> {
   static constexpr bool value =
       std::is_same<T, T_head>::value || Is_in<T, std::tuple<Ts...>>::value;
};

#if __cplusplus >= 201703L  // C++ 17 and over
template <class T, typename TupleTs>
inline constexpr bool is_in = Is_in<T, Tuple>::value;
#endif  // C++ 17 and over

// -- Test if all types in Ts are distinct ------------------------------------

template <typename... Ts>
struct Is_set;

// recursion anchor
template <>
struct Is_set<> {
   static constexpr bool value = true;  // empty set
};

template <typename T, typename... Ts>
struct Is_set<T, Ts...> {
   static constexpr bool value =
       Is_set<Ts...>::value && !Is_in<T, std::tuple<Ts...>>::value;
};

template <typename Tuple>
struct Is_tuple_set;

template <typename... Ts>
struct Is_tuple_set<std::tuple<Ts...>> {
   static constexpr bool value = Is_set<Ts...>::value;
};

#if __cplusplus >= 201703L  // C++ 17 and over
template <typename... Ts>
inline constexpr bool is_set = Is_set<Ts...>::value;
#endif  // C++ 17 and over

// -- Utilities ---------------------------------------------------------------

template <std::size_t... ks>
struct Size_list {
   static constexpr std::size_t length = sizeof...(ks);
};

template <std::size_t k, typename List>
struct Append;

template <std::size_t k, std::size_t... ks>
struct Append<k, Size_list<ks...>> {
   typedef Size_list<ks..., k> list;
};

template <std::size_t... ks>
struct Range;

template <>
struct Range<0> {
   typedef Size_list<> list;
};

template <std::size_t k>
struct Range<k> {
   typedef typename Append<k - 1, typename Range<k - 1>::list>::list list;
};

template <std::size_t k>
using range = typename Range<k>::list;

template <std::size_t... ks>
struct Is_a_range;
template <std::size_t k1>
struct Is_a_range<k1> {
   static constexpr bool value = true;
};
template <std::size_t k1, std::size_t k2, std::size_t... kn>
struct Is_a_range<k1, k2, kn...> {
   static constexpr bool value = (k2 == k1 + 1) && Is_a_range<k2, kn...>::value;
};

template <typename Base, std::size_t k>
struct Indexed : Base {
   static constexpr std::size_t index = k;
};

template <typename... Ts>
struct index {
   template <typename SizeList>
   struct with;
   template <std::size_t... ks>
   struct with<Size_list<ks...>> {
      typedef std::tuple<Indexed<Ts, ks>...> type;
   };
};

template <typename... Ts>
struct Enumerate {
   typedef
       typename index<Ts...>::template with<range<sizeof...(Ts)>>::type type;
};

template <typename... Ts>
using enumerate = typename Enumerate<Ts...>::type;

// -- Find the rank of a given type in a set of types -------------------------

#if __cplusplus < 201703L  // before C++ 17
template <class...>
struct disjunction : std::false_type {};
template <class B1>
struct disjunction<B1> : B1 {};
template <class B1, class... Bn>
struct disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>> {};
template <class...>
struct conjunction : std::true_type {};
template <class B1>
struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};
template <bool B>
using bool_constant = std::integral_constant<bool, B>;
#endif  // before C++ 17

template <typename Tuple>
struct tuple_disjunction;
template <typename... Ts>
struct tuple_disjunction<std::tuple<Ts...>> :
#if __cplusplus >= 201703L  // C++ 17 and over
    std::disjunction<Ts...>
#else  // before C++ 17
    disjunction<Ts...>
#endif
{
};

template <typename Tuple>
struct Is_homogeneous;
template <>
struct Is_homogeneous<std::tuple<>> {
   static constexpr bool value = true;
};
template <typename... Ts>
struct Is_homogeneous<std::tuple<Ts...>> {
   typedef std::tuple_element_t<0, std::tuple<Ts...>> first_type;
   static constexpr bool value =
#if __cplusplus >= 201703L  // C++ 17 and over
       std::conjunction
#else  // before C++ 17
       conjunction
#endif
       <std::is_same<first_type, Ts>..., std::true_type>::value;
};

template <typename T, typename Tuple>
struct Find;

template <typename T, typename... Ts>
struct Find<T, std::tuple<Ts...>> {
   typedef std::tuple<Ts...> Tuple;
   static_assert(Is_tuple_set<Tuple>::value, "not a set");
   static_assert(Is_in<T, Tuple>::value, "not in set");
   static constexpr std::size_t value =
       tuple_disjunction<enumerate<std::is_same<T, Ts>...>>::index;
};

#if __cplusplus >= 201703L  // C++ 17 and over
template <typename T, typename Tuple>
inline constexpr std::size_t find = Find<T, Tuple>::value;
#endif  // C++ 17 and over
