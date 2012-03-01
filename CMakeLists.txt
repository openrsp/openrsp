cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(OpenRSP Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
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
if(NOT "${PARENT_INCLUDE_DIRS}" STREQUAL "")
    foreach (_include ${PARENT_INCLUDE_DIRS})
        include_directories(${_include})
    endforeach()
endif()
if(NOT "${PARENT_MODULE_DIRS}" STREQUAL "")
    foreach (_dir ${PARENT_MODULE_DIRS})
        set(CMAKE_Fortran_MODULE_DIRECTORY
            ${CMAKE_Fortran_MODULE_DIRECTORY}
            ${_dir}
            )
    endforeach()
endif()
message("raboof in openrsp: ${CMAKE_Fortran_MODULE_DIRECTORY}")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(LIBS)
include(ConfigExternal)

include_directories(${PROJECT_SOURCE_DIR})

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/dalton_ifc.F90
    ${PROJECT_SOURCE_DIR}/src/basis_set.f90
    ${PROJECT_SOURCE_DIR}/src/birefring_old.F90
    ${PROJECT_SOURCE_DIR}/src/cgto_diff_eri.F90
    ${PROJECT_SOURCE_DIR}/src/eri_basis_loops.f90
    ${PROJECT_SOURCE_DIR}/src/eri_contractions.f90
    ${PROJECT_SOURCE_DIR}/src/openrsp.F90
    ${PROJECT_SOURCE_DIR}/src/openrsp_daldrv_old.F90
    ${PROJECT_SOURCE_DIR}/src/openrsp_old.F90
    ${PROJECT_SOURCE_DIR}/src/prop_contribs_old.F90
    ${PROJECT_SOURCE_DIR}/src/prop_test_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_backend.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_contribs.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations.f90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_functions.F90
    ${PROJECT_SOURCE_DIR}/src/vib_prop_old.F90
    ${PROJECT_SOURCE_DIR}/src/xprecision.F90
    ${PROJECT_SOURCE_DIR}/src/xtrack.F90
    )

add_library(
    openrsp
    ${SOURCES}
    )

add_dependencies(openrsp external)

install(TARGETS openrsp ARCHIVE DESTINATION lib)
