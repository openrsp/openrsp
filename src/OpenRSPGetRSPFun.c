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
     \param[QInt:int]{in} kn_rule contains the perturbation a and numbers k and n
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
                            const QInt kn_rule[],
                            const QInt size_rsp_fun,
                            QReal *rsp_fun)
{
    QInt num_var_dens=1; /* */
    QMat **var_ao_dens; /* */
    QInt num_int;     /* */
    QMat **val_int;    /* */
    QChar mat_label[3];  /* */
    QInt imat;        /* */

    QErrorCode ierr;  /* error information */
    if (open_rsp->assembled==QFALSE) {
        QErrorExit(FILE_AND_LINE, "OpenRSPAssemble() should be invoked before any calculation");
    }

    var_ao_dens = (QMat **)malloc(num_var_dens*sizeof(QMat*));
    if (var_ao_dens==NULL) {
        printf("OpenRSPGetRSPFun>> number of variable AO density matrices %"QINT_FMT"\n",
               num_var_dens);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for var_ao_dens");
    }
    for (imat=0; imat<num_var_dens; imat++) {
        var_ao_dens[imat] = (QMat *)malloc(sizeof(QMat));
        if (var_ao_dens[imat]==NULL) {
            printf("OpenRSPGetRSPFun>> variable AO density matrix %"QINT_FMT"\n",
                   imat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for var_ao_dens[imat]");
        }
        ierr = QMatCreate(var_ao_dens[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatCreate");
        ierr = QMatBlockCreate(var_ao_dens[imat], 1);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatBlockCreate");
        ierr = QMatSetDimMat(var_ao_dens[imat], 2);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatSetDimMat");
    }

    num_int = 4;
    val_int = (QMat **)malloc(num_int*sizeof(QMat*));
    if (val_int==NULL) {
        printf("OpenRSPGetRSPFun>> number of integral matrices %"QINT_FMT"\n", num_int);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for val_int");
    }
    for (imat=0; imat<num_int; imat++) {
        val_int[imat] = (QMat *)malloc(sizeof(QMat));
        if (val_int[imat]==NULL) {
            printf("OpenRSPGetRSPFun>> integral matrix %"QINT_FMT"\n", imat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for val_int[imat]");
        }
        ierr = QMatCreate(val_int[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatCreate");
        ierr = QMatBlockCreate(val_int[imat], 1);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatBlockCreate");
        ierr = QMatSetDimMat(val_int[imat], 2);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatSetDimMat");
    }
    /* gets overlap integrals */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapGetMat(open_rsp->overlap,
                                num_pert,          /* bra */
                                perturbations,
                                pert_orders,
                                pert_freqs,
                                0,                 /* ket */
                                NULL,
                                NULL,
                                NULL,
                                0,                 /* total */
                                NULL,
                                NULL,
                                num_int,
                                val_int);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapGetMat");
        mat_label[0] = 'S';
        for (imat=0; imat<num_int; imat++) {
            sprintf(&mat_label[1], "%"QINT_FMT"", imat);
            ierr = QMatWrite(val_int[imat], mat_label, ASCII_VIEW);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatWrite(S)");
        }
    }
    /* gets integral matrices of the linked list of one-electron operators */
    if (open_rsp->one_oper!=NULL) {
        ierr = RSPOneOperGetMat(open_rsp->one_oper,
                                num_pert,
                                perturbations,
                                pert_orders,
                                num_int,
                                val_int);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperGetMat");
        mat_label[0] = 'h';
        for (imat=0; imat<num_int; imat++) {
            sprintf(&mat_label[1], "%"QINT_FMT"", imat);
            ierr = QMatWrite(val_int[imat], mat_label, ASCII_VIEW);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatWrite(h)");
        }
    }
    /* gets integral matrices of the linked list of two-electron operators */
    if (open_rsp->two_oper!=NULL) {
        ierr = RSPTwoOperGetMat(open_rsp->two_oper,
                                num_pert,
                                perturbations,
                                pert_orders,
                                num_var_dens,
                                var_ao_dens,
                                num_int,
                                val_int);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperGetMat");
        mat_label[0] = 'G';
        for (imat=0; imat<num_int; imat++) {
            sprintf(&mat_label[1], "%"QINT_FMT"", imat);
            ierr = QMatWrite(val_int[imat], mat_label, ASCII_VIEW);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatWrite(G)");
        }
    }

    for (imat=0; imat<num_var_dens; imat++) {
        ierr = QMatDestroy(var_ao_dens[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatDestroy");
        free(var_ao_dens[imat]);
        var_ao_dens[imat] = NULL;
    }
    free(var_ao_dens);
    var_ao_dens = NULL;

    for (imat=0; imat<num_int; imat++) {
        ierr = QMatDestroy(val_int[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QMatDestroy");
        free(val_int[imat]);
        val_int[imat] = NULL;
    }
    free(val_int);
    val_int = NULL;

    return QSUCCESS;
}
