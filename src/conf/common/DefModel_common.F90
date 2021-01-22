   logical(c_bool) :: locked_context(NbContexte)

   logical(c_bool), public :: phase_can_be_present(NbPhase, NbContexte)

   public :: &
      is_context_locked, &
      lock_context, &
      unlock_context, &
      model_number_of_phases, &
      model_number_of_components, &
      model_number_of_contexts, &
      get_model_configuration

   private :: &
      check_context_id

   contains

   subroutine check_context_id(context)
      integer(c_int), intent(in) :: context
      if (context < 1 .or. context > NbContexte) &
         call CommonMPI_abort('Wrong context identifier.')
   end subroutine check_context_id

   logical(c_bool) function is_context_locked(context) result(lock_status) &
      bind(C, name="is_context_locked")
      integer(c_int), value, intent(in) :: context
      call check_context_id(context)
      lock_status = locked_context(context)
   end function is_context_locked

   subroutine lock_context(context) &
      bind(C, name="lock_context")
      integer(c_int), value, intent(in) :: context
      call check_context_id(context)
      locked_context(context) = .true.
   end subroutine lock_context

   subroutine unlock_context(context) &
      bind(C, name="unlock_context")
      integer(c_int), value, intent(in) :: context
      call check_context_id(context)
      locked_context(context) = .false.
   end subroutine unlock_context

   integer(c_int) function model_number_of_phases() result(n) &
      bind(C, name="model_number_of_phases")
      n = NbPhase
   end function model_number_of_phases

   integer(c_int) function model_number_of_components() result(n) &
      bind(C, name="model_number_of_components")
      n = NbComp
   end function model_number_of_components

   integer(c_int) function model_number_of_contexts() result(n) &
      bind(C, name="model_number_of_contexts")
      n = NbContexte
   end function model_number_of_contexts

   pure type(ModelConfiguration) function get_model_configuration() result(configuration)
      configuration%nb_phases = NbPhase
      configuration%nb_components = NbComp
      configuration%nb_contexts = NbContexte
      configuration%IndThermique = IndThermique
      configuration%NbEqEquilibreMax = NbEqEquilibreMax
      configuration%NbIncPTCMax = NbIncPTCMax
      ! CHECKME: normally arrays are deallocated when they go out of scope
      allocate (configuration%NbPhasePresente_ctx(NbContexte))
      allocate (configuration%NumPhasePresente_ctx(NbPhase, NbContexte))
      allocate (configuration%MCP(NbComp, NbPhase))
      allocate (configuration%pssecd(NbEqFermetureMax, NbContexte))
      configuration%NbPhasePresente_ctx = NbPhasePresente_ctx
      configuration%NumPhasePresente_ctx = NumPhasePresente_ctx
      configuration%MCP = MCP
      configuration%pssecd = pssecd
   end function get_model_configuration

   subroutine DefModel_build_phase_table()

      integer :: ic, k

      phase_can_be_present = .false.
      do ic = 1, NbContexte
         do k = 1, NbPhasePresente_ctx(ic)
            phase_can_be_present(NumPhasePresente_ctx(k, ic), ic) = .true.
         end do
      end do

   end subroutine DefModel_build_phase_table
