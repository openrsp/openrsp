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

# Project name
PROJECT(OpenRSPWEB C)

# WEB files
SET(OPENRSP_LIB_WEB
    ${CMAKE_SOURCE_DIR}/Introduction.nw
    ${CMAKE_SOURCE_DIR}/OpenRSPPerturbation.nw
    ${CMAKE_SOURCE_DIR}/OpenRSPOverlap.nw
    ${CMAKE_SOURCE_DIR}/OpenRSPFortranAPI.nw)

# Try to find notangle and noweave
FIND_PROGRAM(NOTANGLE notangle)
IF(NOT NOTANGLE)
    MESSAGE(FATAL_ERROR "Could not find notangle.")
ENDIF(NOT NOTANGLE)
FIND_PROGRAM(NOWEAVE noweave)
IF(NOT NOWEAVE)
    MESSAGE(FATAL_ERROR "Could not find noweave.")
ENDIF(NOT NOWEAVE)

# Generate LaTeX files using noweave
SET(OPENRSP_DEVELOPER_LATEX_CHAPTERS ${CMAKE_BINARY_DIR}/OpenRSPDeveloper_chapters.tex)
ADD_CUSTOM_COMMAND(OUTPUT ${OPENRSP_DEVELOPER_LATEX_CHAPTERS}
                   COMMAND ${NOWEAVE} -index -n ${OPENRSP_LIB_WEB} > ${OPENRSP_DEVELOPER_LATEX_CHAPTERS})

# Generate PDF file
SET(OPENRSP_DEVELOPER_LATEX_MAIN ${CMAKE_SOURCE_DIR}/OpenRSPDeveloper.tex)
SET(OPENRSP_DEVELOPER_MANUAL ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.pdf)
ADD_CUSTOM_COMMAND(OUTPUT ${OPENRSP_DEVELOPER_MANUAL}
                   DEPENDS ${OPENRSP_DEVELOPER_LATEX_CHAPTERS}
                   COMMAND cp -r ${CMAKE_SOURCE_DIR}/figures ${CMAKE_BINARY_DIR}
                   COMMAND cp ${CMAKE_SOURCE_DIR}/references.bib ${CMAKE_BINARY_DIR}
                   COMMAND pdflatex ${OPENRSP_DEVELOPER_LATEX_MAIN}
                   COMMAND bibtex OpenRSPDeveloper.aux
                   COMMAND makeindex OpenRSPDeveloper.idx
                   COMMAND pdflatex ${OPENRSP_DEVELOPER_LATEX_MAIN}
                   COMMAND pdflatex ${OPENRSP_DEVELOPER_LATEX_MAIN})
# Add the target "doc"
ADD_CUSTOM_TARGET(doc SOURCES ${OPENRSP_DEVELOPER_MANUAL})

# Intermediate files to be cleaned up
SET(OPENRSP_WEB_CLEAN
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.aux
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.bbl
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.blg
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.idx
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.ilg
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.ind
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.log
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.out
    ${CMAKE_BINARY_DIR}/OpenRSPDeveloper.toc
    ${OPENRSP_DEVELOPER_LATEX_CHAPTERS}
    ${OPENRSP_DEVELOPER_MANUAL}
    ${CMAKE_BINARY_DIR}/figures
    ${CMAKE_BINARY_DIR}/references.bib)

# Add intermediate files to be cleaned up
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${OPENRSP_WEB_CLEAN}")
