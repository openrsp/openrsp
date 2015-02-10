# Source codes of Fortran test suite
SET(OPENRSP_F_TEST_SRCS
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_one_oper_exp_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_one_oper_mat_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_two_oper_exp_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_two_oper_mat_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_overlap_exp_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_overlap_mat_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/callback/get_linear_rsp_solution_f.F90
    ${LIB_OPENRSP_PATH}/tests/f90/test_f_OpenRSP_AO.F90
    ${LIB_OPENRSP_PATH}/tests/f90/test_f_OpenRSP.F90)
IF(OPENRSP_PERTURBATION_FREE)
    SET(OPENRSP_F_TEST_SRCS
        ${OPENRSP_F_TEST_SRCS}
        ${LIB_OPENRSP_PATH}/tests/f90/callback/get_pert_comp_f.F90
        ${LIB_OPENRSP_PATH}/tests/f90/callback/get_pert_rank_f.F90)
ENDIF()
