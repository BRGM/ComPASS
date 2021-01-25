template <typename T, typename rank_type = int>
inline constexpr rank_type enum_to_rank(const T enum_value) {
   assert(static_cast<rank_type>(enum_value) > 0);
   return static_cast<rank_type>(enum_value) - 1;  // Fortran -> C indexing
}
