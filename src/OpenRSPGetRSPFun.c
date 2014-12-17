/* OpenRSP: open-ended library for response theory
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

   This file implements the function OpenRSPGetRSPFun().

   2014-07-31, Bin Gao:
   * first version
*/

#include "openrsp.h"

QVoid OpenRSPGetRSPFun_f(const QInt num_pert,
                         const QInt *perturbations,
                         const QInt *pert_orders,
                         const QReal *pert_freqs,
                         const QInt kn_rule[2],
                         const QMat *ref_ham,
                         const QMat *ref_overlap,
                         const QMat *ref_state,
                         RSPSolver *rsp_solver,
                         RSPNucContrib *nuc_contrib,
                         RSPOverlap *overlap,
                         RSPOneOper *one_oper,
                         RSPTwoOper *two_oper,
                         RSPXCFun *xc_fun,
                         const QInt size_rsp_fun,
                         QReal *rsp_fun);

/*@% \brief gets the response function for given perturbations
     \author Bin Gao
     \date 2014-07-31
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QMat:struct]{in} ref_ham Hamiltonian of referenced state
     \param[QMat:struct]{in} ref_state electronic state of referenced state
     \param[QMat:struct]{in} ref_overlap overlap integral matrix of referenced state
     \param[QInt:int]{in} num_pert number of perturbations
     \param[QInt:int]{in} perturbations the perturbations
     \param[QInt:int]{in} pert_orders orders of the perturbations
     \param[QReal:real]{in} pert_freqs frequencies of the perturbations
     \param[QInt:int]{in} kn_rule contains numbers k and n for the kn rule
     \param[QInt:int]{in} size_rsp_fun size of the response function, equals to
         the product of sizes of \var{perturbations}
     \param[QReal:real]{out} rsp_fun the response function
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPGetRSPFun(OpenRSP *open_rsp,
                            const QMat *ref_ham,
                            const QMat *ref_state,
                            const QMat *ref_overlap,
                            const QInt num_pert,
                            const QInt *perturbations,
                            const QInt *pert_orders,
                            const QReal *pert_freqs,
                            const QInt kn_rule[2],
                            const QInt size_rsp_fun,
                            QReal *rsp_fun)
{
    //QErrorCode ierr;  /* error information */
    if (open_rsp->assembled==QFALSE) {
        QErrorExit(FILE_AND_LINE, "OpenRSPAssemble() should be invoked before any calculation");
    }
    switch (open_rsp->elec_EOM_type) {
    /* density matrix-based response theory */
    case ELEC_AO_D_MATRIX:
        OpenRSPGetRSPFun_f(num_pert,
                           perturbations,
                           pert_orders,
                           pert_freqs,
                           kn_rule,
                           ref_ham,
                           ref_overlap,
                           ref_state,
                           open_rsp->rsp_solver,
                           open_rsp->nuc_contrib,
                           open_rsp->overlap,
                           open_rsp->one_oper,
                           open_rsp->two_oper,
                           open_rsp->xc_fun,
                           //id_outp,
                           size_rsp_fun,
                           rsp_fun);
        break;
    /* molecular orbital (MO) coefficient matrix-based response theory */
    case ELEC_MO_C_MATRIX:
        break;
    /* couple cluster-based response theory */
    case ELEC_COUPLED_CLUSTER:
        break;
    default:
        printf("OpenRSPGetRSPFun>> type of EOM of electrons %d\n", open_rsp->elec_EOM_type);
        QErrorExit(FILE_AND_LINE, "invalid type of EOM of electrons");
    }
    return QSUCCESS;
}
