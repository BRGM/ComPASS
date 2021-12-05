set(ComPASSCommonCore_SRCS
    CommonMPI.F90 CommonType.F90 InteroperabilityStructures.F90 MeshInfo.F90
    Physics.F90 SchemeParameters.F90
)

unset(_tmp)
foreach(_src ${ComPASSCommonCore_SRCS})
  list(APPEND _tmp core/${_src})
endforeach()
set(ComPASSCommonCore_SRCS ${_tmp})

set(ComPASSCore_SRCS
    DebugUtils.F90
    GlobalMesh.F90
    IncCV.F90
    IncCVReservoir.F90
    IncCVReservoirTypes.F90
    IncCVWells.F90
    DirichletContribution.F90
    NeumannContribution.F90
    Flux.F90
    Jacobian.F90
    LocalMesh.F90
    IncPrimSecd.F90
    IncPrimSecdTypes.F90
    LoisThermoHydro.cpp
    LoisThermoHydro.F90
    MeshSchema.F90
    Newton.F90
    NumbyContext.F90
    Residu.F90
    Simulation.cpp
    SolvePetsc.F90
    SyncPetsc.F90
    VAGFrac.F90
    NN.F90
    physics/PhysicalConstants.F90
)

set(ComPASSWells_SRCS DefWell.F90 DefFlashWells.F90 WellState.F90)

set(ComPASSMSWells_SRCS DefMSWell.F90)

set(ComPASSFreeFlow_SRCS IncPrimSecdFreeFlow.F90 FreeFlow.F90)

unset(_tmp)
foreach(_src ${ComPASSCore_SRCS})
  list(APPEND _tmp core/${_src})
endforeach()
set(ComPASSCore_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSWells_SRCS})
  list(APPEND _tmp wells/${_src})
endforeach()
set(ComPASSWells_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSMSWells_SRCS})
  list(APPEND _tmp mswells/${_src})
endforeach()
set(ComPASSMSWells_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSFreeFlow_SRCS})
  list(APPEND _tmp freeflow/${_src})
endforeach()
set(ComPASSFreeFlow_SRCS ${_tmp})

set(CONF_SRCS DefFlash.F90 DefModel.F90 Thermodynamics.F90)

set(WRAPPERS_SRCS
    COC.h
    COC_wrappers.cpp
    CommonTypesWrapper.F90
    DebugUtils_wrappers.cpp
    Flux_wrappers.cpp
    Flux_wrappers.F90
    GlobalMeshWrapper.F90
    GlobalMesh_wrappers.cpp
    GlobalVariables_wrappers.F90
    GlobalVariables_wrappers.cpp
    IncCV_wrappers.cpp
    FreeFlow_wrappers.cpp
    IncCV_wrappers.F90
    LinearSystem_wrapper.cpp
    LinearSystemBuilder.h
    LocalMeshWrapper.F90
    MeshUtilities.h
    MeshSchema_wrappers.cpp
    MeshUtilities_wrappers.cpp
    Metis_wrapper.cpp
    Model_wrappers.h
    Model_common_wrappers.cpp
    MeshSchema_wrappers.F90
    NN_wrappers.cpp
    PyBuffer_wrappers.h
    PyBuffer_wrappers.cpp
    Residu_wrappers.cpp
    SolvePetsc_wrappers.F90
    SolvePetsc_wrappers.cpp
    SyncPetsc_wrappers.cpp
    StringWrapper.h
    StringWrapper.F90
    TimeLoop_wrappers.cpp
    VAGFrac_wrappers.cpp
    Well.h
    Well_wrappers.cpp
    Well_wrappers.F90
)

unset(_tmp)
foreach(_src ${WRAPPERS_SRCS})
  list(APPEND _tmp wrappers/${_src})
endforeach()
set(WRAPPERS_SRCS ${_tmp})
