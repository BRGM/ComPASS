#pragma once

#include <tuple>
#include <array>

template <std::size_t ... ns>
struct Range 
{
	static constexpr auto as_tuple() { return std::make_tuple(ns...); }
	template <typename T>
	static auto expand_array(T *p) {
		return std::array<T, sizeof...(ns)>{ {*(p+ns)... } };
	}
	template <typename T>
	static auto expand_array(const T *p) {
		return std::array<T, sizeof...(ns)>{ {*(p + ns)... } };
	}
};

template <typename RangeType, std::size_t n>
struct ExpandedRange;
	
template <std::size_t n, std::size_t ... ns>
struct ExpandedRange<Range<ns...>, n>
{
	typedef Range<ns...,n> type;
};

template <std::size_t n>
struct RangeBuilder
{
	typedef typename ExpandedRange<typename RangeBuilder<n - 1>::range_type, n>::type  range_type;
};

template <>
struct RangeBuilder<0>
{
	typedef Range<0> range_type;
};

template <typename T, std::size_t n>
inline auto make_array(T *p)
{
	return Range<n>::expand_array(p);
}

template <typename T, std::size_t n>
inline auto make_array(const T *p)
{
	return Range<n>::expand_array(p);
}