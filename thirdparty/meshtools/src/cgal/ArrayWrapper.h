#pragma once

#include <vector>

class ArrayWrapper
{
protected:
	void * pointer;
	std::size_t length;
	ArrayWrapper() :
		pointer{ nullptr },
		length{ 0 } {}
public:
	static auto make_empty() { return ArrayWrapper(); }
	template <typename T>
	static auto wrap(T * p, std::size_t size) {
		ArrayWrapper wrapper;
		wrapper.pointer = static_cast<void *>(p);
		wrapper.length = size;
		return wrapper;
	}
	template <typename T>
	static auto wrap(std::vector<T>& v) {
		ArrayWrapper wrapper;
		wrapper.pointer = static_cast<void *>(v.data());
		wrapper.length = v.size();
		return wrapper;
	}
	template <typename T, std::size_t ...extend>
};
