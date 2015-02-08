# Required version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

# Sets the build mode: if the user specifies -DCMAKE_BUILD_TYPE on the command line,
# take their definition and dump it in the cache along with proper documentation,
# otherwise set CMAKE_BUILD_TYPE prior to calling PROJECT(), from
# http://www.cmake.org/pipermail/cmake/2008-September/023808.html
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ELSE()
    SET(CMAKE_BUILD_TYPE None CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF()

# Project name, enable C, CXX and Fortran
PROJECT(OpenRSP C CXX Fortran)

# Host program could include this file by setting LIB_OPENRSP_PATH
IF("${LIB_OPENRSP_PATH}" STREQUAL "")
  SET(LIB_OPENRSP_PATH ${PROJECT_SOURCE_DIR})
ENDIF()

# Adds cmake directory to CMake module
SET(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers)

#
#OPTION(OPENRSP_AUTO_ERROR_EXIT "Enable automatic exit on error." OFF)
#IF(OPENRSP_AUTO_ERROR_EXIT)
#    ADD_DEFINITIONS(-DOPENRSP_AUTO_ERROR_EXIT)
#ENDIF()
OPTION(OPENRSP_PERTURBATION_FREE "Enable perturbation free." ON)
IF(OPENRSP_PERTURBATION_FREE)
    ADD_DEFINITIONS(-DOPENRSP_PERTURBATION_FREE)
ENDIF()
#
OPTION(OPENRSP_USER_CONTEXT "Enable user context in callback functions." OFF)
IF(OPENRSP_USER_CONTEXT)
    ADD_DEFINITIONS(-DOPENRSP_C_USER_CONTEXT)
ENDIF()

#
SET(LIB_OPENRSP_NAME "openrsp" CACHE STRING "Sets the name of the OpenRSP library.")

# Builds the test suite as excutables, otherwise the test suite is included into the OpenRSP library
OPTION(OPENRSP_TEST_EXECUTABLE "Build the test suite excutables." ON)

# Let the user set different application programming interfaces (API)
OPTION(OPENRSP_CXX_API "Build C++ API." OFF)
OPTION(OPENRSP_F03_API "Build Fortran 2003 API." OFF)

# OpenRSP needs the QcMatrix library as well as the matrix I/O functionality
SET(QCMATRIX_INCLUDE_DIR None CACHE STRING "Sets the include directory of QcMatrix library.")
INCLUDE_DIRECTORIES(${QCMATRIX_INCLUDE_DIR})
IF(OPENRSP_F03_API)
    SET(QCMATRIX_MODULE_DIR None CACHE STRING "Sets the module directory of QcMatrix library.")
    INCLUDE_DIRECTORIES(${QCMATRIX_MODULE_DIR})
    # The Fortran API needs the user context in the C part of OpenRSP
    ADD_DEFINITIONS(-DOPENRSP_C_USER_CONTEXT)
    IF(OPENRSP_USER_CONTEXT)
        ADD_DEFINITIONS(-DOPENRSP_F_USER_CONTEXT)
    ENDIF()
ENDIF()
SET(QCMATRIX_LIB None CACHE STRING "Sets the QcMatrix library (xx/libqcmatrix.a).")

# What is this?
option(ENABLE_INTEREST "Build InteRest code" OFF)

# Header files of OpenRSP library
INCLUDE_DIRECTORIES(${LIB_OPENRSP_PATH}/include)

# Where we put the compiled Fortran modules
SET(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules)

# ConfigParentSettings.cmake configures directories where to place module files
# and where to find header files, this is important if there are module and header
# file dependencies between libraries and where the parent code needs to enforce
# a directory structure to compile correctly
INCLUDE(ConfigParentSettings)

# what is this?
set(LIBS)

# Adds the project binary directory to include directories
#INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Debug mode
IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DOPENRSP_DEBUG")
ENDIF()

# Source codes of OpenRSP library, to be merged with Fortran source codes
INCLUDE(Sources)

# Source codes of C test suite
INCLUDE(TestSuiteC)

# Sets different APIs
IF(OPENRSP_CXX_API)

ENDIF()
IF(OPENRSP_F03_API)
    INCLUDE(FortranAPI)
    # Source codes of Fortran test suite
    INCLUDE(TestSuiteFortran)
ENDIF()

# Fortran source codes, to be merged with C source codes
#SET(OPENRSP_SRCS
#    ${OPENRSP_SRCS}
#    src/interfaces/dalton_ifc.F90
#    src/interfaces/interface_1el.F90
#    src/interfaces/interface_1el_dirac.F90
#    src/interfaces/interface_2el.F90
#    src/interfaces/interface_basis.F90
#    src/interfaces/interface_dirac_gen1int.F90
#    src/interfaces/interface_f77_memory.F90
#    src/interfaces/interface_io.F90
#    src/interfaces/interface_molecule.F90
#    src/interfaces/interface_nuclear.F90
#    src/interfaces/interface_rsp_solver.F90
#    src/interfaces/interface_scf.F90
#    src/interfaces/interface_xc.F90
#    src/interfaces/interface_pelib.f90)
#SET(OPENRSP_SRCS
#    ${OPENRSP_SRCS}
#    src/legacy/legacy_birefringence.f90
#    src/legacy/legacy_contributions.F90
#    src/legacy/legacy_openrsp.f90
#    src/legacy/legacy_openrsp_driver.f90
#    src/legacy/legacy_properties.F90
#    src/legacy/legacy_vibrational_properties.f90
#    src/legacy/rsp_equations_old.F90)
#
#    src/ao_dens/rsp_equations.f90
#
#    src/input/openrsp_input_reader.F90
#    src/input/openrsp_cfg.F90
#    src/input/parse_input.F90
#    src/openrsp.f90
#    src/openrsp_efishgcid.f90
#    src/openrsp_general_pv2f.f90
#    src/openrsp_general_pv3f.f90
#    src/openrsp_general_pv4f.f90
#
#    src/ao_dens/rsp_mag_prop.f90
#
#    src/vib_pv_contribs.f90
#
SET(OPENRSP_SRCS
    ${OPENRSP_SRCS}
    src/ao_dens/rsp_choose_rule.f90
    src/ao_dens/rsp_contribs.f90
    src/ao_dens/rsp_field_tuple.f90
    src/ao_dens/rsp_general.f90
    src/ao_dens/rsp_indices_and_addressing.f90
    src/ao_dens/rsp_lof_caching_tmp.f90
    src/ao_dens/rsp_perturbed_matrices.f90
    src/ao_dens/rsp_perturbed_sdf.f90
    src/ao_dens/rsp_property_caching.f90
    src/ao_dens/rsp_sdf_caching.f90
    src/ao_dens/rsp_pert_table.F90
    src/ao_dens/adapter/openrsp_callback_f.F90
    src/ao_dens/adapter/OpenRSPGetRSPFun_f.F90)

# ??
IF(ENABLE_INTEREST)
    ADD_DEFINITIONS(-DENABLE_INTEREST)
    SET(OPENRSP_SRCS
        ${OPENRSP_SRCS}
        src/external-interest/git_revision_info.f90
        src/external-interest/module_interest_eri.f90
        src/external-interest/module_interest_hrr.f90
        src/external-interest/module_interest_one.f90
        src/external-interest/module_interest_osr.f90
        src/interfaces/interface_interest.F90)
ENDIF()

# do we still need it?
INCLUDE(ConfigCompilerFlags)

# Builds executables for the test suite
IF(OPENRSP_TEST_EXECUTABLE)
    ADD_DEFINITIONS(-DOPENRSP_TEST_EXECUTABLE)
    # Creates the OpenRSP library
    ADD_LIBRARY(${LIB_OPENRSP_NAME} ${OPENRSP_SRCS})
    # Creates C test suite
    ADD_EXECUTABLE(${LIB_OPENRSP_NAME}_c_test ${OPENRSP_C_TEST_SRCS})
    # Linked libraries for the C test suite
    TARGET_LINK_LIBRARIES(${LIB_OPENRSP_NAME}_c_test
                          ${LIB_OPENRSP_NAME}
                          ${QCMATRIX_LIB})
    # Creates Fortran test suite
    IF(OPENRSP_F03_API)
        ADD_EXECUTABLE(${LIB_OPENRSP_NAME}_f_test ${OPENRSP_F_TEST_SRCS})
        # Linked libraries for the Fortran test suite
        TARGET_LINK_LIBRARIES(${LIB_OPENRSP_NAME}_f_test
                              ${LIB_OPENRSP_NAME}
                              ${QCMATRIX_LIB})
    ENDIF()
# Creates the OpenRSP library with the test suite
ELSE()
    ADD_LIBRARY(${LIB_OPENRSP_NAME}
                ${OPENRSP_SRCS}
                ${OPENRSP_C_TEST_SRCS}
                ${OPENRSP_F_TEST_SRCS})
ENDIF()

# Installation
INSTALL(TARGETS ${LIB_OPENRSP_NAME} ARCHIVE DESTINATION lib)
