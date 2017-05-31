#pragma once

#include <functional>

template <typename T>
struct XArrayWrapper
{
	typedef T wrapped_type;
	T * pointer;
	std::size_t length;
	XArrayWrapper() :
		pointer{ nullptr },
		length{ 0 } {}
	static auto retrieve(std::function <void(XArrayWrapper&)> bind)
	{
		auto wrapper = XArrayWrapper{};
		bind(wrapper);
		return wrapper;
	}
};


