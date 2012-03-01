add_custom_target(
    git_update
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

set(ExternalProjectCMakeArgs
    -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    )

include(ExternalProject)
add_custom_target(external)

macro(add_external _project)
    ExternalProject_Add(${_project}
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${PROJECT_SOURCE_DIR}/external
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
        BINARY_DIR ${PROJECT_SOURCE_DIR}/external/${_project}-build
        INSTALL_DIR ${PROJECT_SOURCE_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_SOURCE_DIR}/external/${_project}-build)
    include_directories(${PROJECT_SOURCE_DIR}/external/${_project}-build/modules)
    link_directories(${PROJECT_SOURCE_DIR}/external/lib)
    add_dependencies(external ${_project})
    add_dependencies(${_project} git_update)
    set(LIBS
        ${LIBS}
        ${PROJECT_SOURCE_DIR}/external/lib/lib${_project}.a
        )
endmacro()

add_external(matrix-defop)
