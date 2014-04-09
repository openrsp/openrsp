cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(OpenRSP Fortran C)

option(ENABLE_INTEREST "Build InteRest code" OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )
include(ConfigParentSettings)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

set(LIBS)

set(SOURCES
    src/input/openrsp_input_reader.F90
    src/input/openrsp_cfg.F90
    src/input/parse_input.F90
    src/interfaces/dalton_ifc.F90
    src/interfaces/interface_1el.F90
    src/interfaces/interface_1el_dirac.F90
    src/interfaces/interface_2el.F90
    src/interfaces/interface_basis.F90
    src/interfaces/interface_dirac_gen1int.F90
    src/interfaces/interface_f77_memory.F90
    src/interfaces/interface_io.F90
    src/interfaces/interface_molecule.F90
    src/interfaces/interface_nuclear.F90
    src/interfaces/interface_rsp_solver.F90
    src/interfaces/interface_scf.F90
    src/interfaces/interface_xc.F90
    src/interfaces/interface_pelib.f90
    src/legacy/legacy_birefringence.f90
    src/legacy/legacy_contributions.F90
    src/legacy/legacy_openrsp.f90
    src/legacy/legacy_openrsp_driver.f90
    src/legacy/legacy_properties.F90
    src/legacy/legacy_vibrational_properties.f90
    src/legacy/rsp_equations_old.F90
    src/openrsp.f90
    src/openrsp_efishgcid.f90
    src/openrsp_general_pv2f.f90
    src/openrsp_general_pv3f.f90
    src/openrsp_general_pv4f.f90
    src/rsp_choose_rule.f90
    src/rsp_contribs.f90
    src/rsp_equations.f90
    src/rsp_field_tuple.f90
    src/rsp_general.f90
    src/rsp_indices_and_addressing.f90
    src/rsp_lof_caching_tmp.f90
    src/rsp_mag_prop.f90
    src/rsp_perturbed_matrices.f90
    src/rsp_perturbed_sdf.f90
    src/rsp_property_caching.f90
    src/rsp_sdf_caching.f90
    src/vib_pv_contribs.f90
    )

if(ENABLE_INTEREST)
    add_definitions(-DENABLE_INTEREST)
    set(SOURCES
        ${SOURCES}
        src/external-interest/git_revision_info.f90
        src/external-interest/module_interest_eri.f90
        src/external-interest/module_interest_hrr.f90
        src/external-interest/module_interest_one.f90
        src/external-interest/module_interest_osr.f90
        src/interfaces/interface_interest.F90
        )
endif()

include(ConfigCompilerFlags)
include(ConfigDocumentation)

add_library(
    openrsp
    ${SOURCES}
    )

install(TARGETS openrsp ARCHIVE DESTINATION lib)
