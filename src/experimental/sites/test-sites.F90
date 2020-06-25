program test

   use Sites

   implicit none

   integer :: a(nb_sites)

   write (*, *) "Sites:", node_site, cell_site
   write (*, *) a(node_site)

end program test
