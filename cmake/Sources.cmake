# Generates lists of source files from directories and chunks
#
# Header files
if(NOT (OPENRSP_HEADER_DIR AND OPENRSP_HEADER_CHUNKS))
    message(FATAL_ERROR "Directories.cmake and Chunks.cmake not included.")
endif()
set(OPENRSP_HEADER_FILES)
foreach(CHUNK ${OPENRSP_HEADER_CHUNKS})
    set(OPENRSP_HEADER_FILES
        ${OPENRSP_HEADER_FILES}
        ${OPENRSP_HEADER_DIR}/${CHUNK})
endforeach()
# C source codes
if(NOT (OPENRSP_C_DIR AND OPENRSP_C_CHUNKS))
    message(FATAL_ERROR "Directories.cmake and Chunks.cmake not included.")
endif()
set(OPENRSP_C_FILES)
foreach(CHUNK ${OPENRSP_C_CHUNKS})
    set(OPENRSP_C_FILES
        ${OPENRSP_C_FILES}
        ${OPENRSP_C_DIR}/${CHUNK})
endforeach()
# C test source codes
if(NOT (OPENRSP_C_TEST_DIR AND OPENRSP_C_TEST_CHUNKS))
    message(FATAL_ERROR "Directories.cmake and Chunks.cmake not included.")
endif()
set(OPENRSP_C_TEST_FILES)
foreach(CHUNK ${OPENRSP_C_TEST_CHUNKS})
    set(OPENRSP_C_TEST_FILES
        ${OPENRSP_C_TEST_FILES}
        ${OPENRSP_C_TEST_DIR}/${CHUNK})
endforeach()
# Fortran source codes
if(NOT (OPENRSP_FORTRAN_DIR AND OPENRSP_FORTRAN_CHUNKS))
    message(FATAL_ERROR "Directories.cmake and Chunks.cmake not included.")
endif()
set(OPENRSP_FORTRAN_FILES)
foreach(CHUNK ${OPENRSP_FORTRAN_CHUNKS})
    set(OPENRSP_FORTRAN_FILES
        ${OPENRSP_FORTRAN_FILES}
        ${OPENRSP_FORTRAN_DIR}/${CHUNK})
endforeach()
