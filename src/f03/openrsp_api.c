/*
   OpenRSP: open-ended library for response theory
   Copyright 2014

   OpenRSP is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OpenRSP is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

   This file implements the functions for Fortran APIs of OpenRSP library.

   2014-07-31, Bin Gao
   * first version
*/

#include "openrsp.h"

QErrorCode f03_api_OpenRSPCreate(QVoid **open_rsp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)malloc(sizeof(OpenRSP));
    if (c_open_rsp==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for c_open_rsp");
    }
    ierr = OpenRSPCreate(c_open_rsp);
    *open_rsp = (QVoid *)(c_open_rsp);
    return ierr;
}

QErrorCode f03_api_OpenRSPSetElecEOM(QVoid **open_rsp,
                                     const QInt elec_eom_type)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPSetElecEOM(c_open_rsp, elec_eom_type);
    return ierr;
}

QErrorCode f03_api_OpenRSPSetSolver(QVoid **open_rsp,
                                    QVoid *user_ctx,
                                    const QVoid *get_rsp_solution)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPSetSolver(c_open_rsp,
                            user_ctx,
                            get_rsp_solution);
    return ierr;
}

#if defined(OPENRSP_PERTURBATION_FREE)
QErrorCode f03_api_OpenRSPSetPerturbations(QVoid **open_rsp,
                                           const QInt num_pert,
                                           const QInt *perturbations,
                                           const QInt *pert_max_orders,
                                           const QInt *pert_sizes,
                                           QVoid *user_ctx,
                                           const QVoid *get_pert_comp,
                                           const QVoid *get_pert_rank)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPSetPerturbations(c_open_rsp,
                                   num_pert,
                                   perturbations,
                                   pert_max_orders,
                                   pert_sizes,
                                   user_ctx,
                                   get_pert_comp,
                                   get_pert_rank);
    return ierr;
}
#endif

QErrorCode f03_api_OpenRSPSetPDBS(QVoid **open_rsp,
                                  const QInt num_pert,
                                  const QInt *perturbations,
                                  const QInt *pert_max_orders,
                                  QVoid *user_ctx,
                                  const QVoid *get_overlap_mat,
                                  const QVoid *get_overlap_exp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPSetPDBS(c_open_rsp,
                          num_pert,
                          perturbations,
                          pert_max_orders,
                          user_ctx,
                          get_overlap_mat,
                          get_overlap_exp);
    return ierr;
}

QErrorCode f03_api_OpenRSPAddOneOper(QVoid **open_rsp,
                                     const QInt num_pert,
                                     const QInt *perturbations,
                                     const QInt *pert_max_orders,
                                     QVoid *user_ctx,
                                     const QVoid *get_one_oper_mat,
                                     const QVoid *get_one_oper_exp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPAddOneOper(c_open_rsp,
                             num_pert,
                             perturbations,
                             pert_max_orders,
                             user_ctx,
                             get_one_oper_mat,
                             get_one_oper_exp);
    return ierr;
}

QErrorCode f03_api_OpenRSPAddTwoOper(QVoid **open_rsp,
                                     const QInt num_pert,
                                     const QInt *perturbations,
                                     const QInt *pert_max_orders,
                                     QVoid *user_ctx,
                                     const QVoid *get_two_oper_mat,
                                     const QVoid *get_two_oper_exp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPAddTwoOper(c_open_rsp,
                             num_pert,
                             perturbations,
                             pert_max_orders,
                             user_ctx,
                             get_two_oper_mat,
                             get_two_oper_exp);
    return ierr;
}

QErrorCode f03_api_OpenRSPAssemble(QVoid **open_rsp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPAssemble(c_open_rsp);
    return ierr;
}

QErrorCode f03_api_OpenRSPWrite(const QVoid **open_rsp, const QChar *file_name)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPWrite(c_open_rsp, file_name);
    return ierr;
}

QErrorCode f03_api_OpenRSPGetRSPFun(QVoid **open_rsp,
                                    const QVoid *ref_ham,
                                    const QVoid *ref_state,
                                    const QVoid *ref_overlap,
                                    const QInt num_pert,
                                    const QInt *perturbations,
                                    const QInt *pert_orders,
                                    const QReal *pert_freqs,
                                    const QInt *kn_rule,
                                    const QInt size_rsp_fun,
                                    QReal *rsp_fun)
{
    OpenRSP *c_open_rsp;
    QMat *c_ref_ham;
    QMat *c_ref_state;
    QMat *c_ref_overlap;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    c_ref_ham = (QMat *)ref_ham;
    c_ref_state = (QMat *)ref_state;
    c_ref_overlap = (QMat *)ref_overlap;
    ierr = OpenRSPGetRSPFun(c_open_rsp,
                            c_ref_ham,
                            c_ref_state,
                            c_ref_overlap,
                            num_pert,
                            perturbations,
                            pert_orders,
                            pert_freqs,
                            kn_rule,
                            size_rsp_fun,
                            rsp_fun);
    return ierr;
}

QErrorCode f03_api_OpenRSPDestroy(QVoid **open_rsp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPDestroy(c_open_rsp);
    *open_rsp = NULL;
    open_rsp = NULL;
    return ierr;
}
