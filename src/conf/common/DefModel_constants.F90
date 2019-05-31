  ! ! ****** Constants derived from model (do not edit) ****** ! !

  ! Nombre Max d'eq d'equilibre
  !	       d'eq de fermeture thermodynamique
  !	       d'inc P (T) C
  !	       d'inc P (T) C primaires
  integer, parameter :: &
       NbEqEquilibreMax  = NbComp*(NbPhase-1),           & !< Max number of balance equations
       NbIncPTCMax       = 1 + IndThermique + sum(MCP),  & !< Max number of unknowns P (T) C
       NbIncTotalPrimMax = NbComp + IndThermique,        & !< Max number of primary unknowns
       NbCompThermique   = NbComp + IndThermique           !< Number of balance equations (components) + thermique
