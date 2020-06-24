//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <vector>

template <typename Value_type>
struct GenericCOC_iterator;

template <typename Value_type>
struct GenericCOC;

template <typename Value_type>
struct GenericCOC_container {
   typedef Value_type value_type;

  protected:
   value_type* begin_;
   value_type* end_;
   GenericCOC_container() = delete;
   GenericCOC_container(value_type* content, int length)
       : begin_{content}, end_{content + length} {
      assert(length >= 0);
   }

  public:
   value_type* data() const { return begin_; }
   value_type* begin() const { return begin_; }
   value_type* end() const { return end_; }
   std::size_t length() const { return std::distance(begin_, end_); }
   friend class GenericCOC_iterator<value_type>;
};

template <typename Value_type>
struct GenericCOC_iterator {
   typedef Value_type value_type;

  protected:
   // This is const as COC structure is not supposed to be changed
   const int* container_offset;
   value_type* container_content;
   GenericCOC_iterator() = delete;
   GenericCOC_iterator(int* offset, value_type* content)
       : container_offset{offset}, container_content{content} {}
   auto container_length() const {
      return *std::next(container_offset) - *container_offset;
   }

  public:
   auto operator*() const {
      return GenericCOC_container<value_type>{container_content,
                                              container_length()};
   }
   auto operator<(const GenericCOC_iterator& other) const {
      return container_content < other.container_content;
   }
   auto operator==(const GenericCOC_iterator& other) const {
      return container_content == other.container_content;
   }
   auto operator!=(const GenericCOC_iterator& other) const {
      return container_content != other.container_content;
   }
   GenericCOC_iterator& operator++() {
      // WARNING: container_content must be advanced before offset
      //          otherwise length is unvalidated
      std::advance(container_content, container_length());
      std::advance(container_offset, 1);
      return *this;
   }
   friend class GenericCOC<value_type>;
};

/** Container of containers. */
template <typename Value_type>
struct GenericCOC {
   typedef int offset_type;
   typedef Value_type value_type;

  protected:
   std::size_t nb_containers;
   offset_type* container_offset;
   value_type* container_content;

  public:
   GenericCOC()
       : nb_containers{0},
         container_offset{nullptr},
         container_content{nullptr} {}
   auto begin() const {
      return GenericCOC_iterator<value_type>{container_offset,
                                             container_content};
   }
   auto end() const {
      assert(nb_containers >= 0);
      auto offset = container_offset;
      auto content = container_content;
      if (offset != nullptr) {
         assert(nb_containers > 0);
         assert(container_offset != nullptr);
         assert(container_offset[nb_containers] > 0);
         std::advance(offset, nb_containers);
         std::advance(content, container_offset[nb_containers]);
      } else {
         assert(container_offset == nullptr);
      }
      return GenericCOC_iterator<value_type>{offset, content};
   }
   auto operator[](const int i) {
      assert(i >= 0);
      assert(i < nb_containers);
      return *GenericCOC_iterator<value_type>{
          &container_offset[i], &container_content[container_offset[i]]};
   }
   /** Wrap contiguous vector as COC.
   \todo FIXME This is just a convenience function a more robust C++ side API to
   COC<T> is to be developed/used. The problem here is that COC and vectors
   lifetimes are not bound.
   */
   static auto wrap(std::vector<int>& offsets, std::vector<int>& content) {
      assert(offsets.front() == 0);
      assert(offsets.back() == content.size());
      assert(offsets.size() > 0);
      assert(std::all_of(offsets.begin(), offsets.end(),
                         [](const int& i) { return i >= 0; }));
      GenericCOC<value_type> coc;
      coc.nb_containers = offsets.size() - 1;
      coc.container_offset = offsets.data();
      coc.container_content = content.data();
      return coc;
   }
   auto number_of_containers() const { return nb_containers; }
   const offset_type* offset_data() const { return container_offset; }
   offset_type* mutable_offset_data() { return container_offset; }
   const value_type* content_data() const { return container_content; }
   value_type* mutable_content_data() { return container_content; }
   std::size_t size() const {
      if (nb_containers == 0) return 0;
      assert(container_offset);
      assert(container_content);
      assert(container_offset[nb_containers] >= 0);
      return static_cast<std::size_t>(container_offset[nb_containers]);
   }
};

typedef GenericCOC_container<int> COC_container;
typedef GenericCOC_iterator<int> COC_iterator;
typedef GenericCOC<int> COC;
