#pragma once

#include <vector>

class ArrayWrapper
{
protected:
	const void * pointer;
	size_t length;
	ArrayWrapper() :
		pointer(NULL),
		length(0) {}
public:
	template <typename T>
	static auto wrap(const std::vector<T>& v) {
		ArrayWrapper wrapper;
		wrapper.pointer = static_cast<const void *>(v.data());
		wrapper.length = v.size();
		return wrapper;
	}
};