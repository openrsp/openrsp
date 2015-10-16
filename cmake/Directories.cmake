# Directories containing manuals, header and source code files
set(OPENRSP_DOC_DIR ${LIB_OPENRSP_PATH}/doc)
set(OPENRSP_HEADER_DIR ${LIB_OPENRSP_PATH}/include)
set(OPENRSP_C_DIR ${LIB_OPENRSP_PATH}/src)
set(OPENRSP_FORTRAN_DIR ${LIB_OPENRSP_PATH}/src/fortran)
set(OPENRSP_C_TEST_DIR ${LIB_OPENRSP_PATH}/tests)
set(OPENRSP_FORTRAN_TEST_DIR ${LIB_OPENRSP_PATH}/tests/fortran)
# For developers to build from WEB files
if(OPENRSP_BUILD_WEB)
    # We process LaTeX files in the build directory
    set(OPENRSP_LATEX_BUILD_DIR ${CMAKE_BINARY_DIR}/latex)
    # We will create these directories before generating and processing LaTeX,
    # header and source code files
    #
    # For the time being, directories "doc", "include" and "src" are kept in
    # the repository
    set(OPENRSP_MAKE_DIRS
        #${OPENRSP_DOC_DIR}
        #${OPENRSP_HEADER_DIR}
        #${OPENRSP_C_DIR}
        #${OPENRSP_FORTRAN_DIR}
        ${OPENRSP_LATEX_BUILD_DIR})
endif()
