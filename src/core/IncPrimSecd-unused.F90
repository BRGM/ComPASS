   ! WARNING: this algorithm implementation has NOT been checked

   !> \brief  Choice of primary and secondary unknowns
   !! from the matrix dFsurdX with the glouton algorithm
   !! by minimizing the successives angles
   !! \todo IncPrimSecd_IncSecondGlouton not implemented yet ?
   subroutine IncPrimSecd_IncSecondGlouton(cv_info, inc, dFsurdX, &
                                           NumIncTotalPrimCV, NumIncTotalSecondCV)

      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      double precision, intent(in) :: dFsurdX(NbIncTotalMax, NbEqFermetureMax)
      integer, intent(out) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(out) :: NumIncTotalSecondCV(NbEqFermetureMax)

      double precision :: &
         dFsurdXnrm2(NbIncTotalMax), &
         ctrit(NbIncTotalMax), ctrit_max

      integer :: ctrit_maxidxs(NbIncTotalMax)

      double precision :: rnormProjVinc, ss

      double precision :: &
         BaseOrthonormale(NbEqFermetureMax, NbIncTotalMax)

      integer :: is, i, j, nj, j1

      ! FIXME
      call CommonMPI_abort("j1 is used unitialized")

      ! two steps for NumIncTotalSecondcv
      ! 1. first
      ! 2. others

      ! dFsurdXnrm2(i): norm 2 of line i of dFsurdX
      do j = 1, cv_info%NbIncTotal
      do i = 1, cv_info%NbEqFermeture
         dFsurdXnrm2(j) = dFsurdXnrm2(j) + dFsurdX(j, i)**2
      end do
      dFsurdXnrm2(j) = dsqrt(dFsurdXnrm2(j))
      end do

      ! loop for choosing secd inconnues (size is NbEqFermeture), index: is
      do is = 1, cv_info%NbEqFermeture

         ! compute ctrit
         if (is == 1) then

            ctrit(:) = dFsurdXnrm2(:)
         else !

            rnormProjVinc = 0.d0
            do j = 1, cv_info%NbIncTotal ! i: loop index of P T C S

               ss = 0.d0
               do i = 1, cv_info%NbEqFermeture
                  ss = ss + BaseOrthonormale(i, j)*dFsurdX(j, i)
               end do

               rnormProjVinc = rnormProjVinc + ss**2
            end do

            rnormProjVinc = sqrt(rnormProjVinc)

            ! maximise distance = rnormeVinc - rnormeProjVinc
            !   where rnormVinc = dFsurdXnrm2(j)
            ctrit(is) = dFsurdXnrm2(is) - rnormProjVinc
         end if

         ! max of ctrit
         ctrit_max = -100.d0
         do j = 1, cv_info%NbIncPTC
         if (ctrit(j) > ctrit_max) then
            ctrit_max = ctrit(j)
         end if
         end do

         ! ctrit_maxidxs: all elements that takes the max value
         nj = 0
         do j = 1, cv_info%NbIncTotal
         if (abs(ctrit(j) - ctrit_max) < eps) then
            ctrit_maxidxs(nj) = j
            nj = nj + 1
         end if
         end do

         ! which one is second ? (i1, j1)
         ! FIXME: is j1 initialized?
         NumIncTotalSecondCV(is) = j1

         ! update BaseOrthonomale
         do j = 1, cv_info%NbEqFermeture
            BaseOrthonormale(j, is) = dFsurdX(is, j)
         end do

         do i = 1, is - 1

            ss = 0.d0
            do j = 1, cv_info%NbEqFermeture
               ss = ss + BaseOrthonormale(j, i)*dFsurdX(is, j)
            end do

            BaseOrthonormale(:, is) = &
               BaseOrthonormale(:, is) - ss*BaseOrthonormale(:, i)
         end do

         ! normalisation
         ss = 0.d0
         do j = 1, cv_info%NbEqFermeture
            ss = ss + BaseOrthonormale(j, is)**2
         enddo

         ss = dsqrt(ss)
         do j = 1, cv_info%NbEqFermeture
            BaseOrthonormale(j, is) = BaseOrthonormale(j, is)/ss
         enddo

      end do ! end loop of is for choosing second inconnues

   end subroutine IncPrimSecd_IncSecondGlouton
