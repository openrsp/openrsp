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
if(NOT (OPENRSP_SRC_DIR AND OPENRSP_C_CHUNKS))
    message(FATAL_ERROR "Directories.cmake and Chunks.cmake not included.")
endif()
set(OPENRSP_C_FILES)
foreach(CHUNK ${OPENRSP_C_CHUNKS})
    set(OPENRSP_C_FILES
        ${OPENRSP_C_FILES}
        ${OPENRSP_SRC_DIR}/${CHUNK})
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

# Fortran recursive codes and adapters between OpenRSP APIs
set(OPENRSP_AO_DENS_SRCS
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_choose_rule.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_contribs.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_field_tuple.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_general.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_indices_and_addressing.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_lof_caching_tmp.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_perturbed_matrices.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_perturbed_sdf.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_property_caching.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_sdf_caching.f90
    ${LIB_OPENRSP_PATH}/src/ao_dens/rsp_pert_table.F90
    ${LIB_OPENRSP_PATH}/src/ao_dens/adapter/openrsp_callback_f.F90
    ${LIB_OPENRSP_PATH}/src/ao_dens/adapter/OpenRSPGetRSPFun_f.F90)
