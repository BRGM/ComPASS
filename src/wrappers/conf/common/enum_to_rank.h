template <typename T>
inline auto enum_to_rank(const T enum_value) {
   assert(static_cast<int>(enum_value) > 0);
   return static_cast<int>(enum_value) - 1;  // Fortran -> C indexing
}
