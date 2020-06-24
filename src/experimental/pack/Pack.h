#pragma once

#include <array>
#include <vector>

#include "Buffer_view.h"
#include "Template_utils.h"

namespace ComPASS {

// Cf. FIXME Pack::subpack_type_footprint below: limitations with constexpr
// array before C++17

template <typename Dof, typename word_type>
constexpr std::size_t dof_footprint() {
   typedef typename Dof::type Dof_type;
   static_assert(sizeof(Dof_type) >= sizeof(word_type),
                 "inconsistent word size");
   static_assert(sizeof(Dof_type) % sizeof(word_type) == 0,
                 "inconsistent word size");
   return sizeof(Dof_type) / sizeof(word_type);
}

template <typename word_type, std::size_t k, typename... Dofs>
struct Pack_offsets {
   static void compute(
       const std::array<std::size_t, sizeof...(Dofs)>& pack_sizes,
       std::array<std::size_t, sizeof...(Dofs) + 1>& offsets) {
      Pack_offsets<word_type, k - 1, Dofs...>::compute(pack_sizes, offsets);
      typedef typename std::tuple_element<k - 1, std::tuple<Dofs...>>::type Dof;
      offsets[k] =
          offsets[k - 1] + pack_sizes[k - 1] * dof_footprint<Dof, word_type>();
   }
};

template <typename word_type, typename... Dofs>
struct Pack_offsets<word_type, 0, Dofs...> {
   static void compute(
       const std::array<std::size_t, sizeof...(Dofs)>& pack_sizes,
       std::array<std::size_t, sizeof...(Dofs) + 1>& offsets) {
      offsets.fill(0);
   }
};

template <typename... Dofs>
struct Pack {
   static_assert(Is_set<Dofs...>::value, "duplicate dof type");
   typedef std::tuple<Dofs...> subpack_value_types;

  private:
   static constexpr std::size_t nb_subpacks = sizeof...(Dofs);
   typedef std::array<std::size_t, nb_subpacks> Size_array;
   typedef unsigned char word_type;
   std::vector<word_type> raw_data;
   static_assert(
#if __cplusplus >= 201703L  // C++ 17 and over
       std::conjunction<
           std::bool_constant<sizeof(typename Dofs::type) % sizeof(word_type) ==
                              0>...,
           std::true_type>,
#else
       conjunction<bool_constant<
                       sizeof(typename Dofs::type) % sizeof(word_type) == 0>...,
                   std::true_type>::value,
#endif
       "inconsistent sizes in memory");
   /* FIXME
   // It is valid in C++17 to have the following
   // It works with MSVC 17 but fails with gcc 8.2 because of the need for
   compile time definitions
   // cf.
   https://stackoverflow.com/questions/40690260/undefined-reference-error-for-static-constexpr-member?noredirect=1&lq=1
   static constexpr Size_array subpack_type_footprint{
       sizeof(typename Dofs::type) / sizeof(word_type)...
   };
   // this is replaced by the preceding tedious meta programming workaround
   // cf. Pack_offsets struct above
   */
   std::array<std::size_t, nb_subpacks + 1>
       subpack_offset;       // !!! offsets in words !!!
   Size_array subpack_size;  // pack sizes (nb of subpack elements)
   void rebuild() {
      /* // C++ 17 implementation cf FIXME above
      std::size_t offset = 0;
      for (std::size_t k = 0; k < nb_subpacks; ++k) {
          subpack_offset[k] = offset;
          offset += subpack_size[k] * subpack_type_footprint[k];
      }
      subpack_offset[nb_subpacks] = offset;
      raw_data.resize(offset);
      */
      Pack_offsets<word_type, nb_subpacks, Dofs...>::compute(subpack_size,
                                                             subpack_offset);
      raw_data.resize(subpack_offset.back());
      raw_data.shrink_to_fit();
   }
   template <typename Dof_type, typename Self>
   static auto access_subpack_data(Self& self) {
      static_assert(std::is_same<std::remove_const_t<Self>, Pack>::value,
                    "inconsistent pack type");
      static_assert(Dof_type::is_ComPASS_dof, "not a valid dof");
      constexpr std::size_t k = Find<Dof_type, subpack_value_types>::value;
      typedef typename std::tuple_element_t<k, subpack_value_types>::type
          raw_value_type;
      typedef std::conditional_t<std::is_const<Self>::value,
                                 std::add_const_t<raw_value_type>,
                                 raw_value_type>
          value_type;
      auto raw_begin = self.raw_data.data() + self.subpack_offset[k];
      auto begin = reinterpret_cast<value_type*>(raw_begin);
      return Buffer_view<value_type>{begin, begin + self.subpack_size[k]};
   }
   template <typename... Dof_selection, typename Self>
   static auto access_data(Self& self) {
      static_assert(std::is_same<std::remove_const_t<Self>, Pack>::value,
                    "inconsistent pack type");
      static_assert(sizeof...(Dof_selection) != 0, "no dofs provided");
      static_assert(
          Is_homogeneous<std::tuple<typename Dof_selection::type...>>::value,
          "inconsistent subpack types");
      static_assert(
          Is_a_range<Find<Dof_selection, subpack_value_types>::value...>::value,
          "packs must be contiguous to make a view");
      typedef std::tuple_element_t<0, std::tuple<Dof_selection...>> first_dof;
      auto first_view = access_subpack_data<first_dof>(self);
      if (sizeof...(Dof_selection) == 1) return first_view;
      typedef std::tuple_element_t<sizeof...(Dof_selection) - 1,
                                   std::tuple<Dof_selection...>>
          last_dof;
      auto last_view = access_subpack_data<last_dof>(self);
      static_assert(
          std::is_same<decltype(first_view), decltype(last_view)>::value,
          "inconsistent subpack types");  // redundant with previous assertion
      return decltype(first_view){first_view.begin(), last_view.end()};
   }

  public:
   Pack() = default;
   Pack(Size_array&& sizes) : subpack_size{std::forward<Size_array>(sizes)} {
      rebuild();
   }
   template <typename... Ts>
   Pack(Ts... ns) : subpack_size{static_cast<std::size_t>(ns)...} {
      rebuild();
   }
   template <typename... Dof_selection>
   auto view() {
      return access_data<Dof_selection...>(*this);
   }
   template <typename... Dof_selection>
   auto view() const {
      return access_data<Dof_selection...>(*this);
   }
};

template <typename... Dof_selection, typename Pack_type>
auto view(Pack_type& pack) {
   return pack.template view<Dof_selection...>();
}

}  // namespace ComPASS
