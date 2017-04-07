#pragma once

#include <cstddef>
#include <cassert>
#include <iterator>

struct COC_container
{
protected:
	int * begin_;
	int * end_;
	COC_container() = delete;
	COC_container(int * content, int length) :
		begin_{ content },
		end_{ content + length }
	{
		assert(length >= 0);
	}
public:
	int * begin() const { return begin_; }
	int * end() const { return end_; }
	std::size_t length() const { return std::distance(begin_, end_); }
	friend class COC_iterator;
};

class COC_iterator
{
protected:
	// This is const as COC structure is not supposed to be changed
	const int * container_offset;
	int * container_content;
	COC_iterator() = delete;
	COC_iterator(int * offset, int * content) :
		container_offset{ offset },
		container_content{ content }
	{}
	auto container_length() const {
		return *std::next(container_offset) - *container_offset;
	}
public:
	auto operator*() const {
		return COC_container{ container_content, container_length() };
	}
	auto operator<(const COC_iterator& other) const {
		return container_content < other.container_content;
	}
	auto operator==(const COC_iterator& other) const {
		return container_content == other.container_content;
	}
	auto operator!=(const COC_iterator& other) const {
		return container_content != other.container_content;
	}
	COC_iterator& operator++() {
		// WARNING: container_content must be advanced before offset
		//          otherwise length is unvalidated
		std::advance(container_content, container_length());
		std::advance(container_offset, 1);
	}
	friend class COC;
};

/** Container of containers. */
class COC
{
protected:
	int nb_containers;
	int * container_offset;
	int * container_content;
public:
	auto begin() const {
		return COC_iterator{ container_offset, container_content };
	}
	auto end() const {
		assert(nb_containers >= 0);
		assert(container_offset[nb_containers] >= 0);
		auto offset = container_offset;
		auto content = container_content;
		std::advance(offset, nb_containers);
		std::advance(content, container_offset[nb_containers]);
		return COC_iterator{ offset, content };
	}
	auto operator[](const int i) {
		assert(i >= 0);
		assert(i < nb_containers);
		return *COC_iterator{ &container_offset[i], &container_content[container_offset[i]] };
	}
};