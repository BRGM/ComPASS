!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INFO:
!     Model for mswells: diphasic law, gas and liquid in a pipe flow
!     See
!        i) "Multi-segmented non-isothermal compositional liquid gas well model" of Castanon Quiroz & Masson
!        ii) DFM model of Shi et al 2005
!
!
!*This module has been adapted from R. Masson's mswells code.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module DFMHydroMSWells

   use mpi, only: MPI_Abort
   use CommonMPI, only: CommonMPI_abort

   implicit none

   double precision, parameter:: AHydro = 1.2d0
   double precision, parameter:: BHydro = 0.3d0
   double precision, parameter:: Sg1Hydro = 0.2d0
   double precision, parameter:: Sg2Hydro = 0.4d0
   double precision, parameter:: rKuHydro = 1.5d0

   private:: &
      f_Gamma, df_Gamma

   public :: &
      DFMHydroMSWells_MixtureVelocity, &
      DFMHydroMSWells_DPhi, &
      DFMHydroMSWells_FluxVsg, &
      DFMHydroMSWells_f_C0, DFMHydroMSWells_df_C0, &
      DFMHydroMSWells_df_sgC0, DFMHydroMSWells_f_sgC0K, DFMHydroMSWells_df_sgC0K, &
      DFMHydroMSWells_f_G, DFMHydroMSWells_df_G
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Mixture velocity function of variation of potential
!
!       edge = ssp, sp parent node of s
!            Dphi = psp - ps + rhom*g*( zsp -zs )
!
!       Um oriented outward to         s -> sp
!
!       viscom = viscosite dynamique du melange
!            rhom   = densite du melange
!       Redge = rayon equivalent de la section de l'arete
!       SizeEdge = longueur de l'arete
!
!       Umnm0  = vitesse du melange a l'instant precedent: not used here
!
   subroutine DFMHydroMSWells_MixtureVelocity( &
      DPhi, rhom, viscom, Redge, SizeEdge, &
      Umnm1, Um, dPhiUm, dRhoUm, dViscoUm)

      double precision, intent(in)  :: DPhi, rhom, viscom, Redge, SizeEdge, Umnm1
      double precision, intent(out) :: Um, dPhiUm, dRhoUm, dViscoUm
      !tmp
      double precision :: fq, g, b, a, dbv, dar, c, daUm, dbUm, dgUm

      fq = 6.d-2 ! Forchheimer coefficient (turbulent value)
      g = DPhi
      b = SizeEdge*8.d0*viscom/Redge**2
      a = SizeEdge*fq*rhom/(4.d0*Redge)
      dbv = SizeEdge*8.d0/Redge**2
      dar = SizeEdge*fq/(4.d0*Redge)

      if (g .ge. 0.0) then
         c = dsqrt(b**2 + 4.d0*g*a)
         Um = -(c - b)/(2.d0*a)
         daUm = -g/(a*c) - (b - c)/(2.d0*a**2)
         dbUm = (1.d0 - b/c)/(2.d0*a)
         dgUm = -1.d0/c
         dPhiUm = dgUm
         dRhoUm = daUm*dar
         dViscoUm = dbUm*dbv
      else
         c = dsqrt(b**2 - 4.d0*g*a)
         Um = (c - b)/(2.d0*a)
         daUm = -g/(a*c) + (b - c)/(2.d0*a**2)
         dbUm = -(1.d0 - b/c)/(2.d0*a)
         dgUm = -1.d0/c
         dPhiUm = dgUm
         dRhoUm = daUm*dar
         dViscoUm = dbUm*dbv

      endif

   end subroutine DFMHydroMSWells_MixtureVelocity

   subroutine DFMHydroMSWells_DPhi(Um, Rhom, Viscom, Redge, SizeEdge, DPhi)

      double precision, intent(in)  :: Rhom, Viscom, Redge, SizeEdge, Um
      double precision, intent(out) :: DPhi
      double precision :: fq, a

!     Darcy Forchheimer avec linearisation du terme quadratique avec Umnm1

      fq = 6.d-2
      a = 8.d0*Viscom/Redge**2 + fq*Rhom*dabs(Um)/(4.d0*Redge)
      DPhi = -a*Um*SizeEdge

   end subroutine DFMHydroMSWells_DPhi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Flux numerique pour la vitesse superficielle du gaz
!
!
!
!       INPUT: Um, rhol, rhog, sigma, factor
!
!
!       OUTPUT: Flux = vitesse superficielle du gaz
!
!
!       sg1, sg2: valeurs des saturations de gaz aux nodes s et sp (parent)
!
!       Um: vitesse du melange
!
!        factor: facteur d'orientation de l'arete
!
!       rhol, rhog: densites des phases liq et gaz
!
!           WARNING: IL FAUT rhol > rhog
!
!       sigma: tension interfaciale liq gaz
!
!        grav: constante de gravite (9.81)
!
!
!       dsg1: derivee wrt sg1
!
!        dsg2: derivee wrt sg2
!
!       dUm: derivee wrt Um
!
!       drhol: derivee wrt rhol
!
!       drhog: derivee wrt rhog
!
!        dsigma: derivee wrt sigma
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine DFMHydroMSWells_FluxVsg(sg1, sg2, Um, rhol, rhog, sigma, grav, factor, &
                                      Flux, &
                                      dsg1, dsg2, dUm, drhol, drhog, dsigma)

      double precision, intent(in)  :: sg1, sg2, Um, rhol, rhog, sigma, grav, factor
      double precision, intent(inout)  :: Flux, dsg1, dsg2, dUm, drhol, drhog, dsigma
      !tmp
      double precision Uc, ratiod, dratiof, dsf, dsigUc, ss, drholUc, drhogUc

      ratiod = rhog/rhol

      if (ratiod .ge. 1.d0) then
         call CommonMPI_abort(" Error in DFMHydroMSWells_FluxVsg: rhog should be greater than rhol!")
      endif

      Uc = factor*(sigma*grav*(rhol - rhog)/rhol**2)**0.25d0

      !Flux
      if (Um .gt. 0.0) then

         Flux = sg1*DFMHydroMSWells_f_C0(sg1)*Um

      else

         Flux = sg2*DFMHydroMSWells_f_C0(sg2)*Um

      endif

      if (Uc .gt. 0.0) then

         Flux = Flux + DFMHydroMSWells_f_sgC0K(sg1)*DFMHydroMSWells_f_G(sg2, ratiod)*Uc

      else

         Flux = Flux + DFMHydroMSWells_f_sgC0K(sg2)*DFMHydroMSWells_f_G(sg1, ratiod)*Uc

      endif

      !derivee wrt sg1

      dsg1 = 0.d0
      if (Um .gt. 0.0) then

         dsg1 = DFMHydroMSWells_df_sgC0(sg1)*Um

      endif

      if (Uc .gt. 0.0) then

         dsg1 = dsg1 + DFMHydroMSWells_df_sgC0K(sg1)*DFMHydroMSWells_f_G(sg2, ratiod)*Uc

      else

         call DFMHydroMSWells_df_G(sg1, ratiod, dsf, dratiof)
         dsg1 = dsg1 + DFMHydroMSWells_f_sgC0K(sg2)*dsf*Uc

      endif

      !derivee wrt sg2

      dsg2 = 0.d0
      if (Um .le. 0.0) then

         dsg2 = DFMHydroMSWells_df_sgC0(sg2)*Um

      endif

      if (Uc .le. 0.0) then

         dsg2 = dsg2 + DFMHydroMSWells_df_sgC0K(sg2)*DFMHydroMSWells_f_G(sg1, ratiod)*Uc

      else

         call DFMHydroMSWells_df_G(sg2, ratiod, dsf, dratiof)
         dsg2 = dsg2 + DFMHydroMSWells_f_sgC0K(sg1)*dsf*Uc

      endif

      !derivee wrt Um
      if (Um .gt. 0.0) then

         dUm = sg1*DFMHydroMSWells_f_C0(sg1)

      else

         dUm = sg2*DFMHydroMSWells_f_C0(sg2)

      endif

      !derivees wrt rhol, rhog, sigma

      dsigUc = factor*(grav*(rhol - rhog)/rhol**2)**0.25d0
      dsigUc = dsigUc*0.25d0/sigma**0.75d0

      ss = factor*0.25d0*((sigma*grav)**0.25d0) &
           *(1.d0/rhol - rhog/rhol**2)**(-0.75d0)

      drholUc = ss*(2.d0*rhog - rhol)/rhol**3

      drhogUc = -ss/rhol**2

      if (Uc .gt. 0.0) then

         call DFMHydroMSWells_df_G(sg2, ratiod, dsf, dratiof)

         dsigma = DFMHydroMSWells_f_sgC0K(sg1)*DFMHydroMSWells_f_G(sg2, ratiod)*dsigUc

         drhol = DFMHydroMSWells_f_sgC0K(sg1)*DFMHydroMSWells_f_G(sg2, ratiod)*drholUc

         drhol = drhol - DFMHydroMSWells_f_sgC0K(sg1)*Uc*dratiof*rhog/rhol**2

         drhog = DFMHydroMSWells_f_sgC0K(sg1)*DFMHydroMSWells_f_G(sg2, ratiod)*drhogUc

         drhog = drhog + DFMHydroMSWells_f_sgC0K(sg1)*Uc*dratiof/rhol

      else

         call DFMHydroMSWells_df_G(sg1, ratiod, dsf, dratiof)

         dsigma = DFMHydroMSWells_f_sgC0K(sg2)*DFMHydroMSWells_f_G(sg1, ratiod)*dsigUc

         drhol = DFMHydroMSWells_f_sgC0K(sg2)*DFMHydroMSWells_f_G(sg1, ratiod)*drholUc

         drhol = drhol - DFMHydroMSWells_f_sgC0K(sg2)*Uc*dratiof*rhog/rhol**2

         drhog = DFMHydroMSWells_f_sgC0K(sg2)*DFMHydroMSWells_f_G(sg1, ratiod)*drhogUc

         drhog = drhog + DFMHydroMSWells_f_sgC0K(sg2)*Uc*dratiof/rhol

      endif

   end subroutine DFMHydroMSWells_FluxVsg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! function de sg ! on suppose Um/Vsgf assez faible pour ne pas intervenir
!
   function f_Gamma(sg) result(s)

      double precision, intent(in) :: sg
      double precision :: s

      s = (sg - BHydro)/(1.d0 - BHydro)
      s = max(0.d0, s)
      s = min(1.d0, s)

   end function f_Gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     derivee de Gamma fonction de sg ! on suppose Um/Vsgf assez faible pour ne pas intervenir
!!
   function df_Gamma(sg) result(ds)

      double precision, intent(in) :: sg
      double precision :: s, ds

      s = (sg - BHydro)/(1.d0 - BHydro)
      ds = 1.d0/(1.d0 - BHydro)

      if (s .lt. 0.0) then
         ds = 0.d0
      endif
      if (s .gt. 1.0) then
         ds = 0.d0
      endif

   end function df_Gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Profile parameter: fonction de sg

   function DFMHydroMSWells_f_C0(sg) result(s)

      double precision, intent(in) :: sg
      double precision ::  gam, s

      gam = f_Gamma(sg)

      s = AHydro/(1.d0 + (AHydro - 1.d0)*gam**2)

   end function DFMHydroMSWells_f_C0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     derivee du Profile parameter fonction de sg
!!
   function DFMHydroMSWells_df_C0(sg) result(df)

      double precision, intent(in) :: sg
      double precision ::  gam, dgam, s, df

      gam = f_Gamma(sg)
      dgam = df_Gamma(sg)

      s = (1.d0 + (AHydro - 1.d0)*gam**2)

      df = -AHydro*(AHydro - 1.d0)*2.d0*gam*dgam/s**2

   end function DFMHydroMSWells_df_C0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     derivee de sg*C0(sg)
!!
   function DFMHydroMSWells_df_sgC0(sg) result(df)

      double precision, intent(in) ::sg
      double precision :: df

      df = sg*DFMHydroMSWells_df_C0(sg) + DFMHydroMSWells_f_C0(sg)

   end function DFMHydroMSWells_df_sgC0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     fonction sg*C0(sg)*K(sg)
!!
!!     Modifiee par rapport a shi 2005, on interpole sg*C0*K au lieu de K
!!!!
   function DFMHydroMSWells_f_sgC0K(sg) result(f)

      double precision, intent(in):: sg
      double precision :: s, s1, s2, f

      if (sg .lt. Sg1Hydro) then

         s = 1.53d0*sg

      else if (sg .gt. Sg2Hydro) then

         s = rKuHydro*DFMHydroMSWells_f_C0(sg)*sg

      else

         s1 = 1.53d0*Sg1Hydro
         s2 = rKuHydro*DFMHydroMSWells_f_C0(Sg2Hydro)*Sg2Hydro
         s = s1*(sg - Sg2Hydro)/(Sg1Hydro - Sg2Hydro)
         s = s + s2*(sg - Sg1Hydro)/(Sg2Hydro - Sg1Hydro)

      endif

      f = s

   end function DFMHydroMSWells_f_sgC0K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     derivee de la fonction sg*C0(sg)*K(sg)
!!
!!     Modifiee par rapport a shi 2005, on interpole sg*C0*K au lieu de K
   function DFMHydroMSWells_df_sgC0K(sg) result(df)

      double precision, intent(in) :: sg
      double precision :: ds, s1, s2, df

      if (sg .lt. Sg1Hydro) then

         ds = 1.53d0

      else if (sg .gt. Sg2Hydro) then

         ds = rKuHydro*DFMHydroMSWells_df_sgC0(sg)

      else

         s1 = 1.53d0*Sg1Hydro
         s2 = rKuHydro*DFMHydroMSWells_f_C0(Sg2Hydro)*Sg2Hydro
         ds = (s1 - s2)/(Sg1Hydro - Sg2Hydro)

      endif

      df = ds

   end function DFMHydroMSWells_df_sgC0K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     G: function de sg  et de ratiodens = rhog/rhol
!!
   function DFMHydroMSWells_f_G(sg, ratiodens) result(f)

      double precision, intent(in):: sg, ratiodens
      double precision :: c, d, s, f

      c = DFMHydroMSWells_f_C0(sg)

      d = dsqrt(ratiodens)

      s = (1.d0 - c*sg)/(sg*c*d + 1.d0 - sg*c)

      f = s

   end function DFMHydroMSWells_f_G
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!     derivee de G wrt sg et ratiodens
!!
   subroutine DFMHydroMSWells_df_G(sg, ratiodens, dsf, dratiof)

      double precision, intent(in) :: sg, ratiodens
      double precision, intent(inout) ::  dsf, dratiof
      double precision c, d, e

      c = DFMHydroMSWells_f_C0(sg)

      d = dsqrt(ratiodens)

      e = 1.d0/(sg*c*d + 1.d0 - sg*c)**2

      dsf = -d*DFMHydroMSWells_df_sgC0(sg)*e

      dratiof = -(1.d0 - c*sg)*e*sg*c/2.d0/d

   end subroutine DFMHydroMSWells_df_G

end module DFMHydroMSWells
