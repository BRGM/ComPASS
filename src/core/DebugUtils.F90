!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DebugUtils

  use CommonMPI, only: commRank

   use MeshSchema, only: &
      NbFracLocal_Ncpus, &
      NodebyFractureLocal, NodebyCellLocal, &
      IdNodeLocal, XNodeLocal, XFaceLocal, &
      NodeFlagsLocal, CellFlagsLocal, FaceFlagsLocal

   implicit none

   public :: &
      DebugUtils_is_own_frac_node, &
      DebugUtils_dump_mesh_info

contains

   function DebugUtils_is_own_frac_node(node) result(found)

      integer, intent(in) :: node
      logical             :: found
      integer             :: pf, k

      found = .false.
      do pf = 1, NbFracLocal_Ncpus(commRank + 1)
         do k = NodebyFractureLocal%Pt(pf) + 1, NodebyFractureLocal%Pt(pf + 1)
            if (NodebyFractureLocal%Num(k) == node) then
               found = .true.
               return
            end if
         end do
      end do

   end function DebugUtils_is_own_frac_node

   subroutine DebugUtils_dump_mesh_info() &
       bind(C, name = "debug_utils_dump_mesh_info")

      character(len=1024) :: filename
      integer :: k, Ierr, errcode

      write (filename, "(A5I0.5)") "nodes", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(XNodeLocal, 2)
         write (24, *) XNodeLocal(:, k)
      end do
      close (24)

      write (filename, "(A5I0.5)") "cells", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, NodebyCellLocal%Nb
         write (24, *) NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)) - 1
      end do
      close (24)

      write (filename, "(A8I0.5)") "nodeinfo", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(IdNodeLocal)
         write (24, *) IdNodeLocal(k)%proc, ' ', IdNodeLocal(k)%frac, ' ', &
                       IdNodeLocal(k)%P,    ' ', IdNodeLocal(k)%T
      end do
      close (24)

      write (filename, "(A9I0.5)") "fractures", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, NodebyFractureLocal%Nb
         write (24, *) NodebyFractureLocal%Num(NodebyFractureLocal%Pt(k) + 1:NodebyFractureLocal%Pt(k + 1)) - 1
      end do
      close (24)

      write (filename, "(A9I0.5)") "nodeflags", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(NodeFlagsLocal, 1)
         write (24, *) NodeFlagsLocal(k)
      end do
      close (24)

      write (filename, "(A9I0.5)") "cellflags", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(CellFlagsLocal, 1)
         write (24, *) CellFlagsLocal(k)
      end do
      close (24)

      write (filename, "(A11I0.5)") "facecenters", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(XFaceLocal, 2)
         write (24, *) XFaceLocal(:, k)
      end do
      close (24)

      write (filename, "(A9I0.5)") "faceflags", commRank + 1
      open (unit=24, file=trim(filename), status='REPLACE')
      do k = 1, size(FaceFlagsLocal, 1)
         write (24, *) FaceFlagsLocal(k)
      end do
      close (24)

   end subroutine DebugUtils_dump_mesh_info

end module DebugUtils
