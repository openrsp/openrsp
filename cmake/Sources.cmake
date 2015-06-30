# Source codes of OpenRSP library
SET(OPENRSP_SRCS
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertCreate.c
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertAssemble.c
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertWrite.c
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertGetNumComps.c
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertGetConcatenation.c
    ${LIB_OPENRSP_PATH}/src/perturbation/RSPPertDestroy.c
    ${LIB_OPENRSP_PATH}/src/solver/RSPSolverCreate.c
    ${LIB_OPENRSP_PATH}/src/solver/RSPSolverAssemble.c
    ${LIB_OPENRSP_PATH}/src/solver/RSPSolverWrite.c
    ${LIB_OPENRSP_PATH}/src/solver/RSPSolverGetLinearRSPSolution.c
    ${LIB_OPENRSP_PATH}/src/solver/RSPSolverDestroy.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapCreate.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapAssemble.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapWrite.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapGetMat.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapGetExp.c
    ${LIB_OPENRSP_PATH}/src/overlap/RSPOverlapDestroy.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperCreate.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperAdd.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperAssemble.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperWrite.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperGetMat.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperGetExp.c
    ${LIB_OPENRSP_PATH}/src/one_oper/RSPOneOperDestroy.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperCreate.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperAdd.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperAssemble.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperWrite.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperGetMat.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperGetExp.c
    ${LIB_OPENRSP_PATH}/src/two_oper/RSPTwoOperDestroy.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunCreate.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunAdd.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunAssemble.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunWrite.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunGetMat.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunGetExp.c
    ${LIB_OPENRSP_PATH}/src/xc_fun/RSPXCFunDestroy.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonCreate.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonAssemble.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonWrite.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonGetContributions.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonDestroy.c
    ${LIB_OPENRSP_PATH}/src/nuc_contrib/RSPNucHamiltonGetNumAtoms.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPCreate.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPSetWaveFunction.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPSetPerturbations.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPSetLinearRSPSolver.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPSetPDBS.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPAddOneOper.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPAddTwoOper.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPAddXCFun.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPSetNucContributions.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPAssemble.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPWrite.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPGetRSPFun.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPGetResidue.c
    ${LIB_OPENRSP_PATH}/src/OpenRSPDestroy.c)
# Fortran recursive codes and adapters between OpenRSP APIs
SET(OPENRSP_SRCS
    ${OPENRSP_SRCS}
    src/ao_dens/rsp_choose_rule.f90
    src/ao_dens/rsp_contribs.f90
    src/ao_dens/rsp_field_tuple.f90
    src/ao_dens/rsp_general.f90
    src/ao_dens/rsp_indices_and_addressing.f90
    src/ao_dens/rsp_lof_caching_tmp.f90
    src/ao_dens/rsp_perturbed_matrices.f90
    src/ao_dens/rsp_perturbed_sdf.f90
    src/ao_dens/rsp_property_caching.f90
    src/ao_dens/rsp_sdf_caching.f90
    src/ao_dens/rsp_pert_table.F90
    src/ao_dens/adapter/openrsp_callback_f.F90
    src/ao_dens/adapter/OpenRSPGetRSPFun_f.F90)
