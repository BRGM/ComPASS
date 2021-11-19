#pragma once

#include <utility>

// CHECKME: will be available in C++-23
template <typename Enum>
inline constexpr std::underlying_type_t<Enum> to_underlying(Enum e) noexcept {
   static_assert(std::is_enum_v<Enum>);
   return static_cast<std::underlying_type_t<Enum>>(e);
}

template <typename Enum>
inline constexpr std::underlying_type_t<Enum> enum_to_rank(Enum e) noexcept {
   static_assert(std::is_enum_v<Enum>);
   assert(to_underlying(e) > 0);
   return to_underlying(e) - 1;  // Fortran -> C indexing
}
