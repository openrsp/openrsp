# Required version of CMake
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Project name
project(OpenRSP C CXX Fortran)

# Host program could include this file by setting LIB_OPENRSP_PATH
if("${LIB_OPENRSP_PATH}" STREQUAL "")
    set(LIB_OPENRSP_PATH ${PROJECT_SOURCE_DIR})
endif()

# Adds cmake directory to CMake module
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${LIB_OPENRSP_PATH}/cmake)

# Sets the build mode: if the user specifies -DCMAKE_BUILD_TYPE on the command
# line, take their definition and dump it in the cache along with proper
# documentation, otherwise set CMAKE_BUILD_TYPE prior to calling project(),
# from http://www.cmake.org/pipermail/cmake/2008-September/023808.html
if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

# Options for making the library
include(cmake/custom/user_context.cmake)
include(Options)

# ConfigParentSettings.cmake configures directories where to place module files
# and where to find header files, this is important if there are module and
# header file dependencies between libraries and where the parent code needs to
# enforce a directory structure to compile correctly
include(ConfigParentSettings)

# Directories containing manuals, header and source code files
include(Directories)

# Chunks of the library
include(Chunks)

# Generates lists of source files from the directories and chunks
include(Sources)

# For developers to build from WEB files
if(OPENRSP_BUILD_WEB)
    # WEB files
    include(WEBFiles)
    # Adds the taret "openrsp_dirs" for making directories for building
    configure_file("${LIB_OPENRSP_PATH}/cmake/cmake_make_dirs.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_make_dirs.cmake"
                   @ONLY)
    add_custom_target(openrsp_dirs
                      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_make_dirs.cmake)
    # Adds the target "latex" to copy LaTeX files to build directory; It has
    # problem to compile pdf in the build directory and with LaTeX files in
    # their original source directory, especially under Cygwin
    configure_file("${LIB_OPENRSP_PATH}/cmake/cmake_copy_latex.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_copy_latex.cmake"
                   @ONLY)
    add_custom_target(latex
                      DEPENDS openrsp_dirs
                      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_copy_latex.cmake)
    # Settings of developer manual
    include(DeveloperManual)
    # Generates main chapter of the developer manual using noweave; For some
    # reason, the generated chapter is not copied by cmake_copy_latex.cmake at
    # first; So we copy it immediately after generating
    find_program(NOWEAVE noweave)
    if(NOT NOWEAVE)
        message(FATAL_ERROR "Could not find noweave.")
    endif()
    add_custom_command(TARGET latex PRE_BUILD
                       COMMAND ${NOWEAVE} -index -n ${OPENRSP_WEB_FILES} >
                               ${OPENRSP_DEVELOPER_WEAVED_TEX}
                       COMMAND ${CMAKE_COMMAND} -E copy
                               ${OPENRSP_DEVELOPER_WEAVED_TEX} ${OPENRSP_LATEX_BUILD_DIR})
    # Generates PDF file of the developer manual using LaTeX
    include(FindLATEX)
    if(NOT (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER))
        message(FATAL_ERROR "Could not find pdflatex, bibtex and/or makeindex.")
    endif()
    # Adds the target "pdf" to make the developer manual, and copies the PDF
    # file to the directory "doc"
    add_custom_target(pdf
                      COMMAND ${PDFLATEX_COMPILER} ${OPENRSP_DEVELOPER_MANUAL}.tex
                      COMMAND ${BIBTEX_COMPILER} ${OPENRSP_DEVELOPER_MANUAL}.aux
                      COMMAND ${MAKEINDEX_COMPILER} ${OPENRSP_DEVELOPER_MANUAL}.idx
                      COMMAND ${PDFLATEX_COMPILER} ${OPENRSP_DEVELOPER_MANUAL}.tex
                      COMMAND ${PDFLATEX_COMPILER} ${OPENRSP_DEVELOPER_MANUAL}.tex
                      COMMAND ${CMAKE_COMMAND} -E copy
                              ${OPENRSP_DEVELOPER_MANUAL}.pdf ${OPENRSP_DOC_DIR}
                      DEPENDS latex openrsp_dirs
                      WORKING_DIRECTORY ${OPENRSP_LATEX_BUILD_DIR})
    # We could use -L'#line %L "%F"%N' for C codes using ${NOTANGLE} for debugging?
    #
    # Generates source codes using notangle
    find_program(NOTANGLE notangle)
    if(NOT NOTANGLE)
        message(FATAL_ERROR "Could not find notangle.")
    endif()
    # Header files
    foreach(CHUNK ${OPENRSP_HEADER_CHUNKS})
        set(TANGLED_SRC ${OPENRSP_HEADER_DIR}/${CHUNK})
        add_custom_command(OUTPUT ${TANGLED_SRC}
                           DEPENDS openrsp_dirs ${OPENRSP_WEB_FILES}
                           COMMAND ${NOTANGLE} -R${CHUNK}
                                   ${OPENRSP_WEB_FILES} > ${TANGLED_SRC})
    endforeach()
    set_source_files_properties(${OPENRSP_HEADER_FILES}
                                PROPERTIES GENERATED TRUE)
    # C source codes
    foreach(CHUNK ${OPENRSP_C_CHUNKS})
        set(TANGLED_SRC ${OPENRSP_C_DIR}/${CHUNK})
        add_custom_command(OUTPUT ${TANGLED_SRC}
                           DEPENDS openrsp_dirs ${OPENRSP_WEB_FILES}
                           COMMAND ${NOTANGLE} -R${CHUNK}
                                   ${OPENRSP_WEB_FILES} > ${TANGLED_SRC})
    endforeach()
    set_source_files_properties(${OPENRSP_C_FILES}
                                PROPERTIES GENERATED TRUE)
    # C test source codes
    foreach(CHUNK ${OPENRSP_C_TEST_CHUNKS})
        set(TANGLED_SRC ${OPENRSP_C_TEST_DIR}/${CHUNK})
        add_custom_command(OUTPUT ${TANGLED_SRC}
                           DEPENDS openrsp_dirs ${OPENRSP_WEB_FILES}
                           COMMAND ${NOTANGLE} -R${CHUNK}
                                   ${OPENRSP_WEB_FILES} > ${TANGLED_SRC})
    endforeach()
    set_source_files_properties(${OPENRSP_C_TEST_FILES}
                                PROPERTIES GENERATED TRUE)
    # Fortran source codes
    foreach(CHUNK ${OPENRSP_FORTRAN_CHUNKS})
        set(TANGLED_SRC ${OPENRSP_FORTRAN_DIR}/${CHUNK})
        add_custom_command(OUTPUT ${TANGLED_SRC}
                           DEPENDS openrsp_dirs ${OPENRSP_WEB_FILES}
                           COMMAND ${NOTANGLE} -R${CHUNK}
                                   ${OPENRSP_WEB_FILES} > ${TANGLED_SRC})
    endforeach()
    set_source_files_properties(${OPENRSP_FORTRAN_FILES}
                                PROPERTIES GENERATED TRUE)
    # Adds intermediate files to be cleaned up
    set(OPENRSP_WEB_CLEAN
        ${OPENRSP_DEVELOPER_WEAVED_TEX}
        ${OPENRSP_MAKE_DIRS})
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${OPENRSP_WEB_CLEAN}")
endif()

# Header files and Fortran module files
include_directories(${OPENRSP_HEADER_DIR})

option(OPENRSP_FORTRAN_API "Build Fortran 2003 APIs." OFF)
if(OPENRSP_FORTRAN_API)
    set(QCMATRIX_MODULE_DIR None
        CACHE STRING "Directory of Fortran modules of QcMatrix library.")
    include_directories(${QCMATRIX_MODULE_DIR})
    # OpenRSP Fortran APIs need to enable user context of OpenRSP C parts
    add_definitions(-DOPENRSP_C_USER_CONTEXT)
    if(OPENRSP_USER_CONTEXT)
        add_definitions(-DOPENRSP_F_USER_CONTEXT)
    endif()
endif()

if(OPENRSP_FORTRAN_API)
    if("${CMAKE_Fortran_MODULE_DIRECTORY}" STREQUAL "")
        set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
    endif()
    include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

# Data for unit testing
include(OpenRSPTestData)

# Builds test suite as executables
option(OPENRSP_TEST_EXECUTABLE "Build test suite as excutables." ON)
if(OPENRSP_TEST_EXECUTABLE)
    add_definitions(-DOPENRSP_TEST_EXECUTABLE)
    if(OPENRSP_FORTRAN_API)
        add_library(openrsp
                    ${OPENRSP_HEADER_FILES}
                    ${OPENRSP_FORTRAN_FILES}
                    ${OPENRSP_C_FILES}
          )
    else()
        add_library(openrsp
                    ${OPENRSP_HEADER_FILES}
                    ${OPENRSP_C_FILES}
          )
    endif()
    add_executable(openrsp_c_test ${OPENRSP_C_TEST_FILES})
    target_link_libraries(openrsp_c_test
                          openrsp
                          ${QCMATRIX_LIB})
    enable_testing()
    add_test(
      NAME openrsp_c_test
      COMMAND $<TARGET_FILE:openrsp_c_test>
      )
# Test suite is included into the library as functions
else()
    if(OPENRSP_FORTRAN_API)
        add_library(openrsp
                    ${OPENRSP_HEADER_FILES}
                    ${OPENRSP_FORTRAN_FILES}
                    ${OPENRSP_C_FILES}
                    ${OPENRSP_C_TEST_FILES}
          )
    else()
        add_library(openrsp
                    ${OPENRSP_HEADER_FILES}
                    ${OPENRSP_C_FILES}
                    ${OPENRSP_C_TEST_FILES}
          )
    endif()
endif()

include(src/ao_dens/CMakeLists.txt)
include(src/ao_dens/adapter/CMakeLists.txt)

# Installs the static library, header files and Fortran module files
install(TARGETS openrsp
        ARCHIVE DESTINATION openrsp/lib)
install(FILES ${OPENRSP_HEADER_FILES}
        DESTINATION openrsp/include)
if(OPENRSP_FORTRAN_API)
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
            DESTINATION openrsp/include)
endif()

# Uninstalls target, modified from (accessing date 2015-09-27)
# http://cmake.org/Wiki/CMake_FAQ#Can_I_do_.22make_uninstall.22_with_CMake.3F
configure_file("${LIB_OPENRSP_PATH}/cmake/cmake_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
               @ONLY)
add_custom_target(uninstall
                  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
