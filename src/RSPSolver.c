/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

*/

#include "RSPSolver.h"

/* <function name='RSPSolverCreate'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Create the context of response equation solver, should be called at first
     <param name='rsp_solver' direction='inout'>
       The context of response equation solver
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_linear_rsp_solution' direction='in'>
       User-specified callback function of linear response equation solver
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPSolverCreate(RSPSolver *rsp_solver,
#if defined(OPENRSP_C_USER_CONTEXT)
                           void *user_ctx,
#endif
                           const GetLinearRSPSolution get_linear_rsp_solution)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_solver->user_ctx = user_ctx;
#endif
    rsp_solver->get_linear_rsp_solution = get_linear_rsp_solution;
    return QSUCCESS;
}

/* <function name='RSPSolverAssemble'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Assembles the context of response equation solver
     <param name='rsp_solver' direction='inout'>
       The context of response equation solver
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPSolverAssemble(RSPSolver *rsp_solver)
{
/*FIXME: to implement? */
    return QSUCCESS;
}

/* <function name='RSPSolverWrite'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Writes the context of response equation solver
     <param name='rsp_solver' direction='in'>
       The context of response equation solver
     </param>
     <param name='fp_solver' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPSolverWrite(const RSPSolver *rsp_solver, FILE *fp_solver)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    if (rsp_solver->user_ctx!=NULL) {
        fprintf(fp_solver, "RSPSolverWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}

/* <function name='RSPSolverGetLinearRSPSolution'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Solve the linear response equation
     <param name='rsp_solver' direction='in'>
       The context of response equation solver
     </param>
     <param name='num_pert' direction='in'>
       Number of different perturbations on the right hand side of the linear
       response equation
     </param>
     <param name='num_comps' direction='in'>
       Number of components of each perturbation, size is <num_pert>
     </param>
     <param name='num_freq_sums' direction='in'>
       For each perturbation, number of complex frequency sums on the left hand
       side of the linear response equation, size is <num_pert>
     </param>
     <param name='freq_sums' direction='in'>
       The complex frequency sums on the left hand side of the linear response
       equation, size is twice of the sum of <num_freq_sums>
     </param>
     <param name='RHS_mat' direction='in'>
       Right-hand-side (RHS) matrices, size is the dot product of <num_comps>
       and <num_freq_sums>, and index of <num_freq_sums> runs faster in memory
     <param name='rsp_param' direction='inout'>
       Solved response parameters, size is the dot product of <num_comps> and
       <num_freq_sums>, and index of <num_freq_sums> runs faster in memory
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver *rsp_solver,
                                         const QInt num_pert,
                                         const QInt *num_comps,
                                         const QInt *num_freq_sums,
                                         const QReal *freq_sums,
                                         QcMat *RHS_mat[],
                                         QcMat *rsp_param[])
{
    rsp_solver->get_linear_rsp_solution(num_pert,
                                        num_comps,
                                        num_freq_sums,
                                        freq_sums,
                                        RHS_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                        rsp_solver->user_ctx,
#endif
                                        rsp_param);
    return QSUCCESS;
}

/* <function name='RSPSolverDestroy'
             attr='private'
             author='Bin Gao'
             date='2014-08-05'>
     Destroys the context of response equation solver, should be called at the end
     <param name='rsp_solver' direction='inout'>
       The context of response equation solver
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPSolverDestroy(RSPSolver *rsp_solver)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_solver->user_ctx = NULL;
#endif
    rsp_solver->get_linear_rsp_solution = NULL;
    return QSUCCESS;
}

