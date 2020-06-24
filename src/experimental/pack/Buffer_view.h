#pragma once

#include <cassert>
#include <type_traits>

namespace ComPASS {

template <typename T>
struct Buffer_view {
   typedef T value_type;
   typedef std::add_pointer_t<T> value_pointer;
   typedef std::add_pointer_t<std::add_const_t<std::remove_cv_t<T>>>
       const_value_pointer;
   typedef std::add_lvalue_reference_t<T> value_lreference;
   value_pointer first;
   value_pointer past_last;
   auto begin() { return first; }
   auto end() { return past_last; }
   const_value_pointer begin() const { return first; }
   const_value_pointer end() const { return past_last; }
   std::size_t size() const { return past_last - first; }
   value_lreference operator[](std::size_t k) const {
      assert(first + k < past_last);
      return *(first + k);
   }
};

}  // namespace ComPASS
