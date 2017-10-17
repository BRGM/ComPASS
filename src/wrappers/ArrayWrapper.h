//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

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
	static auto wrap(std::vector<T>& v) {
		ArrayWrapper wrapper;
		wrapper.pointer = static_cast<void *>(v.data());
		wrapper.length = v.size();
		return wrapper;
	}
	template <typename T, std::size_t ...extend>
	friend class PyBufferWrapper;
};
