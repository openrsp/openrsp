cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(OpenRSP Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

include(ConfigDoc)

if(NOT "${PARENT_DEFINITIONS}" STREQUAL "")
    foreach (_definition ${PARENT_DEFINITIONS})
        add_definitions(${_definition})
    endforeach()
else()
    add_definitions(-DOPENRSP_STANDALONE)
    add_definitions(-DBUILD_OPENRSP)
    include_directories(${PROJECT_SOURCE_DIR}/include)
endif()
if(NOT "${PARENT_INCLUDE_DIR}" STREQUAL "")
    include_directories(${PARENT_INCLUDE_DIR})
endif()
if(NOT "${PARENT_MODULE_DIR}" STREQUAL "")
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PARENT_MODULE_DIR})
endif()

set(LIBS)

include_directories(${PROJECT_SOURCE_DIR})

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/birefring_old.F90
    ${PROJECT_SOURCE_DIR}/src/openrsp.F90
    ${PROJECT_SOURCE_DIR}/src/openrsp_daldrv_old.F90
    ${PROJECT_SOURCE_DIR}/src/openrsp_old.F90
    ${PROJECT_SOURCE_DIR}/src/prop_contribs_old.F90
    ${PROJECT_SOURCE_DIR}/src/prop_test_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_contribs.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_functions.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_general.F90
    ${PROJECT_SOURCE_DIR}/src/vib_prop_old.F90
    ${PROJECT_SOURCE_DIR}/src/nuc_contributions.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/dalton_ifc.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_molecule.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_xc.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_1el.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_scf.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_f77_memory.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_io.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_pcm.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_rsp_solver.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_basis.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_dirac_gen1int.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/interface_interest.F90
    ${PROJECT_SOURCE_DIR}/src/input/parse_input.F90
    ${PROJECT_SOURCE_DIR}/src/input/openrsp_cfg.F90
    ${PROJECT_SOURCE_DIR}/src/input/input_reader.F90
    ${PROJECT_SOURCE_DIR}/src/external-interest/git_revision_info.f90
    ${PROJECT_SOURCE_DIR}/src/external-interest/module_interest.f90
    ${PROJECT_SOURCE_DIR}/src/external-interest/module_interest_eri.f90
    ${PROJECT_SOURCE_DIR}/src/external-interest/module_interest_eri_rkb.f90
    ${PROJECT_SOURCE_DIR}/src/external-interest/module_interest_hrr.f90
    ${PROJECT_SOURCE_DIR}/src/external-interest/module_interest_osr.f90
    )

include(ConfigCompilerFlags)

add_library(
    openrsp
    ${SOURCES}
    )

install(TARGETS openrsp ARCHIVE DESTINATION lib)
