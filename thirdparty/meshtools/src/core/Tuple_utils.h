#pragma once

#include <tuple>

template <typename TupleType>
struct SplitTuple;

template <typename T, typename ... Ts>
struct SplitTuple<std::tuple<T, Ts...>>
{
	typedef T First_element_type;
	typedef std::tuple<Ts...> Sub_tuple_type;
};

template <typename T>
struct SplitTuple<std::tuple<T>>
{
	typedef T First_element_type;
};

template <typename T, typename TupleType>
struct IsInTuple;

template <typename T, typename ... Ts>
struct IsInTuple<T, std::tuple<Ts...>>
{
	typedef std::tuple<Ts...> Tuple_type;
	typedef typename SplitTuple<Tuple_type>::First_element_type First_element_type;
	typedef typename SplitTuple<Tuple_type>::Sub_tuple_type Sub_tuple_type;
	constexpr static bool value = std::is_same<T, First_element_type>::value || IsInTuple<T, Sub_tuple_type>::value;
};

template <typename T1, typename T2>
struct IsInTuple<T1, std::tuple<T2>>
{
	constexpr static bool value = std::is_same<T1, T2>::value;
};

template <typename TupleType>
struct IsSingleTypeTuple;

template <typename ... Ts>
struct IsSingleTypeTuple<std::tuple<Ts...>>
{
	typedef std::tuple<Ts...> Tuple_type;
	typedef typename SplitTuple<Tuple_type>::First_element_type First_element_type;
	typedef typename SplitTuple<Tuple_type>::Sub_tuple_type Sub_tuple_type;
	constexpr static bool value = std::is_same<First_element_type, typename SplitTuple<Sub_tuple_type>::First_element_type>::value && IsSingleTypeTuple<Sub_tuple_type>::value;
};

template <typename T>
struct IsSingleTypeTuple<std::tuple<T>>
{
	typedef T First_element_type;
	constexpr static bool value = true;
};


template <typename TupleType>
inline constexpr bool is_single_type() { return IsSingleTypeTuple<TupleType>::value; }

template <typename TupleType, typename T, bool>
struct ExtendTupleSet;

template <typename T, typename ... Ts>
struct ExtendTupleSet<std::tuple<Ts...>, T, true>
{
	typedef std::tuple<Ts..., T> Extended_tuple_type;
};

template <typename T, typename ... Ts>
struct ExtendTupleSet<std::tuple<Ts...>, T, false>
{
	typedef std::tuple<Ts...> Extended_tuple_type;
};

template <typename TupleType>
struct TupleSet;

template <typename ... Ts>
struct TupleSet<std::tuple<Ts...>>
{
	typedef std::tuple<Ts...> Tuple_type;
	typedef typename SplitTuple<Tuple_type>::First_element_type First_element_type;
	typedef typename SplitTuple<Tuple_type>::Sub_tuple_type Sub_tuple_type;
	typedef typename TupleSet<Sub_tuple_type>::type Sub_tuple_set_type;
	typedef typename ExtendTupleSet<Sub_tuple_set_type, First_element_type, !IsInTuple<First_element_type, Sub_tuple_set_type>::value>::Extended_tuple_type type;
};

template <typename T>
struct TupleSet<std::tuple<T>>
{
	typedef std::tuple<T> Tuple_type;
	typedef T First_element_type;
	typedef Tuple_type type;
};

//template <typename TupleType, typename T>
//struct ExtendTuple;
//
//template <typename T, typename ... Ts>
//struct ExtendTuple<std::tuple<Ts...>, T>
//{
//	typedef std::tuple<Ts..., T> Extended_tuple_type;
//};
//

template <typename ... TupleTypes>
struct ConcatenatedTuple;

template <typename TupleType, typename ... TupleTypes>
struct ConcatenatedTuple<TupleType, TupleTypes...>
{
	typedef typename ConcatenatedTuple<TupleType, typename ConcatenatedTuple<TupleTypes...>::type>::type type;
};

template <typename ... Ts1, typename ... Ts2>
struct ConcatenatedTuple<std::tuple<Ts1...>, std::tuple<Ts2...>>
{
	typedef std::tuple<Ts1..., Ts2...> type;
};

template <typename ... Ts>
struct ConcatenatedTuple<std::tuple<Ts...>>
{
	typedef std::tuple<Ts...> type;
};

template <typename ... TupleTypes>
struct ConcatenatedTupleSet;

template <typename TupleType, typename ... TupleTypes>
struct ConcatenatedTupleSet<TupleType, TupleTypes...>
{
	typedef typename TupleSet<typename ConcatenatedTuple<typename TupleSet<TupleType>::type, typename ConcatenatedTupleSet<TupleTypes...>::type>::type>::type type;
};

template <typename ... Ts1, typename ... Ts2>
struct ConcatenatedTupleSet<std::tuple<Ts1...>, std::tuple<Ts2...>>
{
	typedef typename TupleSet<std::tuple<Ts1..., Ts2...>>::type type;
};

template <typename ... Ts>
struct ConcatenatedTupleSet<std::tuple<Ts...>>
{
	typedef typename TupleSet<std::tuple<Ts...>>::type type;
};

