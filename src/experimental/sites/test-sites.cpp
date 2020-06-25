#include <array>
#include <iostream>

#include "sites.h"

int main() {
   std::array<int, ComPASS::nb_sites> a;
   std::cout << a[ComPASS::node_site] << std::endl;
   return 0;
}
