!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecdTypes

   use iso_c_binding, only: c_int

   use DefModel, only: &
      NbComp, NbPhase, MCP, NbIncPTCMax, &
      NbIncTotalMax, pschoice, psprim, pssecd, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, NbEqEquilibreMax

   use NumbyContext, only: &
      NbEqFermeture_ctx, NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, &
      NumIncComp2NumIncPTC_ctx, NbIncTotalPrim_ctx, NumCompCtilde_ctx, &
      NbEqEquilibre_ctx, NbIncPTC_ctx, NbIncTotal_ctx, NbCompCtilde_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   implicit none

   type ControlVolumeInfo
      integer(c_int) :: &
         NbPhasePresente, NbCompCtilde, &
         NbEqFermeture, NbEqEquilibre, & ! NbEqEquilibre -> Nombre d'Equation d'Equilibre thermodynamique fct du contexte (i.e. égalité des fugacités)
         NbIncPTC, NbIncPTCPrim, NbIncTotal, &
         NbIncTotalPrim, &
         NumPhasePresente(NbPhase), & ! Num -> identifiants de la phase présente
         NumCompCtilde(NbComp), & ! Num -> identifiants des composants absents
         NumCompEqEquilibre(NbEqEquilibreMax), & ! identifiant des composant présents dans au moins 2 phases (donc concernés par égalité fugacités)
         Num2PhasesEqEquilibre(2, NbEqEquilibreMax), & ! phases impliquées dans l'équilibre (cf. NumCompEqEquilibre)
         NumIncPTC2NumIncComp_comp(NbIncPTCMax), & ! Etant donné une ligne du "vecteur inconnu" quel est le composant
         NumIncPTC2NumIncComp_phase(NbIncPTCMax), & ! Etant donné une ligne du "vecteur inconnu" quelle est la phase
         NumIncComp2NumIncPTC(NbComp, NbPhase) ! matrice donnant pour chaque phase et chaque composant la ligne du "vecteur inconnu"
   end type ControlVolumeInfo

   public:: &
      IncPrimSecdTypes_collect_cv_info

contains

   subroutine IncPrimSecdTypes_collect_cv_info(context, cv_info)

      integer, intent(in) :: context
      type(ControlVolumeInfo), intent(out) :: cv_info

      cv_info%NbPhasePresente = NbPhasePresente_ctx(context)
      cv_info%NbCompCtilde = NbCompCtilde_ctx(context)

      cv_info%NbEqFermeture = NbEqFermeture_ctx(context)
      cv_info%NbEqEquilibre = NbEqEquilibre_ctx(context)

      cv_info%NbIncPTC = NbIncPTC_ctx(context)
      cv_info%NbIncPTCPrim = cv_info%NbIncPTC - cv_info%NbEqFermeture

      cv_info%NbIncTotal = NbIncTotal_ctx(context)

      cv_info%NbIncTotalPrim = NbIncTotalPrim_ctx(context)

      cv_info%NumPhasePresente(:) = NumPhasePresente_ctx(:, context)
      cv_info%NumCompCtilde(:) = NumCompCtilde_ctx(:, context)
      cv_info%NumCompEqEquilibre(:) = NumCompEqEquilibre_ctx(:, context)

      cv_info%NumIncPTC2NumIncComp_comp(:) = NumIncPTC2NumIncComp_comp_ctx(:, context)
      cv_info%NumIncPTC2NumIncComp_phase(:) = NumIncPTC2NumIncComp_phase_ctx(:, context)

      cv_info%Num2PhasesEqEquilibre(:, :) = Num2PhasesEqEquilibre_ctx(:, :, context)
      cv_info%NumIncComp2NumIncPTC(:, :) = NumIncComp2NumIncPTC_ctx(:, :, context)

   end subroutine IncPrimSecdTypes_collect_cv_info

end module IncPrimSecdTypes
