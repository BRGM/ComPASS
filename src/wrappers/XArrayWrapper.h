//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <cassert>
#include <functional>

template <typename T>
struct XArrayWrapper {
   typedef T wrapped_type;
   T* pointer;
   std::size_t length;
   XArrayWrapper() : pointer{nullptr}, length{0} {}
   XArrayWrapper(T* p, std::size_t n) : pointer{p}, length{n} { assert(p); }
   static auto retrieve(std::function<void(XArrayWrapper&)> bind) {
      auto wrapper = XArrayWrapper{};
      bind(wrapper);
      return wrapper;
   }
   T& operator[](const std::size_t k) {
      assert(k < length);
      return *(pointer + k);
   }
   void fill(const T& t) {
      auto last = pointer + length;
      for (auto p = pointer; p != last; ++p) {
         (*p) = t;
      }
   }
};
