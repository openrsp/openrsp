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
OPTION(OPENRSP_PERTURBATION_FREE "Enable perturbation free." OFF)
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

# Adds the project binary directory to include directories
#INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Debug mode
IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DOPENRSP_DEBUG")
ENDIF()

# Source codes of OpenRSP library
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

# do we still need it?
#INCLUDE(ConfigCompilerFlags)

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
