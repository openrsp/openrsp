# Source codes of Fortran APIs
SET(OPENRSP_SRCS
    ${OPENRSP_SRCS}
    ${LIB_OPENRSP_PATH}/src/f03/openrsp_api_f.c
    ${LIB_OPENRSP_PATH}/src/f03/rsp_solver_f.F90
    ${LIB_OPENRSP_PATH}/src/f03/rsp_overlap_f.F90
    ${LIB_OPENRSP_PATH}/src/f03/rsp_one_oper_f.F90
    ${LIB_OPENRSP_PATH}/src/f03/rsp_two_oper_f.F90
    ${LIB_OPENRSP_PATH}/src/f03/openrsp_f.F90)
IF(OPENRSP_PERTURBATION_FREE)
    SET(OPENRSP_SRCS
        ${OPENRSP_SRCS}
        ${LIB_OPENRSP_PATH}/src/f03/rsp_pert_f.F90)
ENDIF()
