module LeafMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use Thermodynamics
   use DefModel
   use mpi
   use CommonMPI
   use IncCVReservoir
   use MeshSchema
   use VSHydroMSWells
#else
   use iso_c_binding, only: c_double, c_char

   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort
   use Thermodynamics, only: &
#ifdef _THERMIQUE_
      f_MolarEnthalpy, &
#endif
      f_Viscosity, f_MolarDensity
   use DefModel, only: &
      NbPhase, NbComp, &
      NbIncTotalPrimMax, &
      NbPhasePresente_ctx, NumPhasePresente_ctx

#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      GAS_PHASE, LIQUID_PHASE
#endif

   use CommonType, only: CSR

   use IncCVReservoir, only: &
      NumPhasePresente_ctx, NbPhasePresente_ctx

   use MeshSchema, only: &
      XNodeLocal, &
      NodebyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus

   use VSHydroMSWells, only: VSHydroMSWell

#endif
   implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INFO:
!       The purpose of this module is to create the data structures
!       for the  leaves of each mswell, and initialize them.
!       This  only concerns the "only mswell model", i.e., when there is no interaction with the reservoir
!
!TODO: The allocation & initialization process should be done only when we are in the "only mswell model",
!      and  on the python layer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Data to impose BDCs at leaf nodes
   type :: TYPE_LeafDataMSWell

      real(c_double) :: Q(NbIncTotalPrimMax)       !  mass flux
      real(c_double) :: Svel(NbPhase)              !  superficial velocity
      real(c_double) :: Pression, Temperature      !  Input pressure and temperature
      real(c_double) :: Comp(NbComp)
      real(c_double) :: PressionRes                ! Reservoir Pressure
      real(c_double) :: TFluxRes(NbIncTotalPrimMax) ! Flux between reservoir and pressure
      character(c_char) :: &
         BDC !  'f' for flowrate mode without reservoir; 'r' flowrate with reservoir data

   end type TYPE_LeafDataMSWell

   integer, protected :: NbMSWellLeafNodes

   ! MSWell-LeafNodes
   type(CSR), protected :: &
      LeafNodebyMSWellLocal

   type(TYPE_LeafDataMSWell), allocatable, dimension(:), target, protected :: &
      LeafDataMSWell !Values for input at the leaf nodes for mswells

   public :: &
      LeafMSWells_allocate, &
      LeafMSWells_free, &
      LeafMSWells_init

contains

   !> \brief Allocate MSWell-LeafNodes

   subroutine LeafMSWells_allocate()

      integer :: nb_mswells, Nnz, s, num_s, k
      real(c_double) ::  x_s, z_s

      nb_mswells = NodebyMSWellLocal%Nb !number of wells
      LeafNodebyMSWellLocal%Nb = nb_mswells
      if (nb_mswells .eq. 0) then !No MSWells
         allocate (LeafNodebyMSWellLocal%Pt(2))
         allocate (LeafNodebyMSWellLocal%Val(1))
         allocate (LeafNodebyMSWellLocal%Num(1))
         LeafNodebyMSWellLocal%Pt(1) = 0
         LeafNodebyMSWellLocal%Pt(2) = 0
         LeafNodebyMSWellLocal%Val(1) = 0

      else if (nb_mswells .eq. 1) then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! One Vertical well
         !!!number of leafs
         NbMSWellLeafNodes = 1
         Nnz = NbMSWellLeafNodes
         allocate (LeafNodebyMSWellLocal%Pt(nb_mswells + 1))
         allocate (LeafNodebyMSWellLocal%Num(Nnz))
         allocate (LeafNodebyMSWellLocal%Val(Nnz))
         LeafNodebyMSWellLocal%Pt(1) = 0
         LeafNodebyMSWellLocal%Pt(2) = 1
         LeafNodebyMSWellLocal%Val(1) = 1 !Local idx at LeafNodebyMSWellLocal
         k = 1
         LeafNodebyMSWellLocal%Num(1) = NodebyMSWellLocal%Pt(k) + 1 !idx node at mswell
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! One Chair well
         !!!number of leafs
         !!NbMSWellLeafNodes = 2
         !!Nnz = NbMSWellLeafNodes
         !!allocate (LeafNodebyMSWellLocal%Pt(nb_mswells + 1))
         !!allocate (LeafNodebyMSWellLocal%Num(Nnz))
         !!allocate (LeafNodebyMSWellLocal%Val(Nnz))
         !!LeafNodebyMSWellLocal%Pt(1) = 0
         !!LeafNodebyMSWellLocal%Pt(2) = 2

         !!!!!First leaf is the bottom at the origin
         !!k = 1
         !!LeafNodebyMSWellLocal%Num(1) = NodebyMSWellLocal%Pt(k) + 1 !idx node at mswell
         !!LeafNodebyMSWellLocal%Num(1) = 1 !Local idx node at mswell
         !!LeafNodebyMSWellLocal%Val(1) = 1 !Local idx at LeafNodebyMSWellLocal

         !!!!look for the other leaf at the bottom
         !!do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
         !!   num_s = NodebyMSWellLocal%Num(s)

         !!   x_s = XNodeLocal(1, num_s)    ! x-cordinate of node s
         !!   z_s = XNodeLocal(3, num_s)    ! z-cordinate of node s

         !!   if ((abs(z_s) < 1d-8) .and. (abs(x_s) > 1.d0)) then

         !!      LeafNodebyMSWellLocal%Val(2) = 2 !Local idx at LeafNodebyMSWellLocal
         !!      LeafNodebyMSWellLocal%Num(2) = s !Local idx node at mswell
         !!      exit
         !!   endif
         !!end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else if (nb_mswells .eq. 2) then !Two mswells

         NbMSWellLeafNodes = 2
         Nnz = NbMSWellLeafNodes
         allocate (LeafNodebyMSWellLocal%Pt(nb_mswells + 1))
         allocate (LeafNodebyMSWellLocal%Num(Nnz))
         allocate (LeafNodebyMSWellLocal%Val(Nnz))
         !Leaf of the first MSWell
         LeafNodebyMSWellLocal%Pt(1) = 0
         LeafNodebyMSWellLocal%Pt(2) = 1
         LeafNodebyMSWellLocal%Val(1) = 1 !Local idx at LeafNodebyMSWellLocal
         k = 1
         LeafNodebyMSWellLocal%Num(1) = NodebyMSWellLocal%Pt(k) + 1 !idx node at mswell
         !Leaf of the second MSWell
         LeafNodebyMSWellLocal%Pt(3) = 2
         LeafNodebyMSWellLocal%Val(2) = 2 !Local idx at LeafNodebyMSWellLocal
         k = 2
         LeafNodebyMSWellLocal%Num(2) = NodebyMSWellLocal%Pt(k) + 1 !idx node at mswell

      else
         call CommonMPI_abort("Incorrect number of mswells in LeafMSWells_allocate!")
      end if

      allocate (LeafDataMSWell(Nnz))

   end subroutine LeafMSWells_allocate

   !> \brief Deallocate well unknowns vectors
   subroutine LeafMSWells_free()

      deallocate (LeafDataMSWell)
      deallocate (LeafNodebyMSWellLocal%Val)
      deallocate (LeafNodebyMSWellLocal%Num)
      deallocate (LeafNodebyMSWellLocal%Pt)

   end subroutine LeafMSWells_free

   subroutine LeafMSWells_init() &
      bind(C, name="LeafMSWells_init_data")

      integer :: leaf_data_idx_s, nbwells, k, s
      real(c_double) ::  hgas, hliq, zgas, zliquid

      real(c_double) ::  PwIn, TwIn, VslIn, VsgIn, PRes, WIP
      real(c_double) ::  sgas, sliq, viscogas, viscoliq, Tliq, Tgas, Th20, Ten
      character(c_char)    ::  BDCIn
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !User Parameters, TODO: should be done using the Python interface
      !variables at leaf nodes
#if ComPASS_NUMBER_OF_COMPONENTS == 1
      BDCIn = 'r' ! 'f' for flowrate mode without reservoir data; 'r' flowrate with reservoir data
      VsgIn = 0.0d0 !water2phase
#elif ComPASS_NUMBER_OF_COMPONENTS == 2
      BDCIn = 'f' ! 'f' for flowrate mode without reservoir data
      VsgIn = 0.5d0 !immiscible
#endif
      PwIn = 1.d7
      TwIn = 520.d0
      VslIn = 1.0d0 - VsgIn
      PRes = PwIn + 1.d6 !Reservoir pressure at leaf Node
      WIP = 1.d-12 !Well index
      sgas = 0.d0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort("Multi-segmented wells are only implemented for diphasic physics!")

#else

      nbwells = NbMSWellLocal_Ncpus(commRank + 1)

      do k = 1, nbwells

         !Only mswells producers
         if (DataMSWellLocal(k)%Well_type == 'i') cycle

         do s = LeafNodebyMSWellLocal%Pt(k) + 1, LeafNodebyMSWellLocal%Pt(k + 1)

            leaf_data_idx_s = LeafNodebyMSWellLocal%Val(s) ! idx at the LeafData

            LeafDataMSWell(leaf_data_idx_s)%Svel(LIQUID_PHASE) = VslIn
            LeafDataMSWell(leaf_data_idx_s)%Svel(GAS_PHASE) = VsgIn
            LeafDataMSWell(leaf_data_idx_s)%Pression = PwIn
            LeafDataMSWell(leaf_data_idx_s)%Temperature = TwIn
            LeafDataMSWell(leaf_data_idx_s)%Comp(:) = 1 !This works only for one component
            LeafDataMSWell(leaf_data_idx_s)%Pression = PwIn
            LeafDataMSWell(leaf_data_idx_s)%PressionRes = PRes
            LeafDataMSWell(leaf_data_idx_s)%BDC = BDCIn

            ! molar density
            zliquid = f_MolarDensity(LIQUID_PHASE, &
                                     LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                     LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                     LeafDataMSWell(leaf_data_idx_s)%Comp)

            zgas = f_MolarDensity(GAS_PHASE, &
                                  LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                  LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                  LeafDataMSWell(leaf_data_idx_s)%Comp)
            ! viscosite molaire
            viscoliq = f_Viscosity(LIQUID_PHASE, &
                                   LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                   LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                   LeafDataMSWell(leaf_data_idx_s)%Comp)

            viscogas = f_Viscosity(GAS_PHASE, &
                                   LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                   LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                   LeafDataMSWell(leaf_data_idx_s)%Comp)
            ! enthalpy
            hliq = f_MolarEnthalpy(LIQUID_PHASE, &
                                   LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                   LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                   LeafDataMSWell(leaf_data_idx_s)%Comp)

            hgas = f_MolarEnthalpy(GAS_PHASE, &
                                   LeafDataMSWell(leaf_data_idx_s)%Pression, &
                                   LeafDataMSWell(leaf_data_idx_s)%Temperature, &
                                   LeafDataMSWell(leaf_data_idx_s)%Comp)

            sliq = 1.d0 - sgas
            Tgas = zgas/viscogas*sgas
            Tliq = zliquid/viscoliq*sliq
            Th20 = WIP*(Tgas + Tliq)
            Ten = WIP*(Tgas*hgas + Tliq*hliq)

            if (NbComp .eq. 1) then !Water2Phase
               !QIn , without using reservoir data
               LeafDataMSWell(leaf_data_idx_s)%Q(NbComp) = &
                  LeafDataMSWell(leaf_data_idx_s)%Svel(LIQUID_PHASE)*zliquid &
                  + &
                  LeafDataMSWell(leaf_data_idx_s)%Svel(GAS_PHASE)*zgas

               !TFluxRes, using reservoir data
               LeafDataMSWell(leaf_data_idx_s)%TFluxRes(NbComp) = Th20

            else  !Immiscible
               !QIn, without using reservoir data
               LeafDataMSWell(leaf_data_idx_s)%Q(1) = &   !AIR_COMP ==1
                  LeafDataMSWell(leaf_data_idx_s)%Svel(GAS_PHASE)*zgas

               LeafDataMSWell(leaf_data_idx_s)%Q(2) = &   !LIQ_COMP ==2
                  LeafDataMSWell(leaf_data_idx_s)%Svel(LIQUID_PHASE)*zliquid
            endif

#ifdef _THERMIQUE_
            !QIn
            LeafDataMSWell(leaf_data_idx_s)%Q(NbComp + 1) = &
               LeafDataMSWell(leaf_data_idx_s)%Svel(LIQUID_PHASE)*zliquid*hliq &
               + &
               LeafDataMSWell(leaf_data_idx_s)%Svel(GAS_PHASE)*zgas*hgas

            !TFluxRes
            LeafDataMSWell(leaf_data_idx_s)%TFluxRes(NbComp + 1) = Ten
#endif

         end do !LeafNode

         !Initialize Superficial velocity for all node alongs the mswell
         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            VSHydroMSWell(s)%edata%MixVel = VslIn
         end do

      end do !k
#endif
   end subroutine LeafMSWells_init

end module LeafMSWells
