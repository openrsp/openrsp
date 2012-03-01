cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(OpenRSP Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
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
endif()
if(NOT "${PARENT_INCLUDE_DIRS}" STREQUAL "")
    foreach (_include ${PARENT_DEFINITIONS})
        include_directories(${_include})
    endforeach()
endif()

set(LIBS)
include(ConfigExternal)

add_definitions(-DOPENRSP_STANDALONE)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/prop_contribs_old.F90
    ${PROJECT_SOURCE_DIR}/src/prop_test_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_backend.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_contribs.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations.f90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations_old.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_functions.f90
    )

add_library(
    openrsp
    ${SOURCES}
    )
add_dependencies(openrsp external)

#add_executable(
#    xcint.x
#    src/test_driver.F90
#    )
#
#target_link_libraries(
#    xcint.x
#    xcint
#    ${LIBS}
#    )
#
#install(TARGETS xcint ARCHIVE DESTINATION lib)
