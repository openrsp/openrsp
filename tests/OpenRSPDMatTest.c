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


   This file tests the AO density matrix-based response theory.

   2014-07-31, Bin Gao
   * first version
*/

#include "OpenRSPDMatTest.h"

#if defined(QCMATRIX_SINGLE_PRECISION)
#define QCREAL_FMT "f"
#else
#define QCREAL_FMT "lf"
#endif

QErrorCode OpenRSPDMatTest(OpenRSP *open_rsp, FILE *fp_log)
{
#include "OpenRSPTestPerturbations.h"
    /* overlap operator */
    QInt overlap_num_pert=2;
    QcPertInt overlap_pert_labels[]={PERT_GEOMETRIC,
                                     PERT_MAGNETIC};
    QInt overlap_pert_orders[]={MAX_ORDER_GEOMETRIC,
                                MAX_ORDER_MAGNETIC};
    /* one-electron Hamiltonian */
    QInt oneham_num_pert=2;
    QcPertInt oneham_pert_labels[]={PERT_GEOMETRIC,
                                    PERT_MAGNETIC};
    QInt oneham_pert_orders[]={MAX_ORDER_GEOMETRIC,
                               MAX_ORDER_MAGNETIC};
    /* external field */
    QInt ext_field_num_pert=3;
    QcPertInt ext_field_pert_labels[]={PERT_GEOMETRIC,
                                       PERT_DIPOLE,
                                       PERT_MAGNETIC};
    QInt ext_field_pert_orders[]={MAX_ORDER_GEOMETRIC,
                                  MAX_ORDER_DIPOLE,
                                  MAX_ORDER_MAGNETIC};
    /* two-electron operator */
    QInt two_oper_num_pert=2;
    QcPertInt two_oper_pert_labels[]={PERT_GEOMETRIC,
                                      PERT_MAGNETIC};
    QInt two_oper_pert_orders[]={MAX_ORDER_GEOMETRIC,
                                 MAX_ORDER_MAGNETIC};
    /* referenced state */
    QcMat *F_unpert[1];
    QcMat *D_unpert[1];
    QcMat *S_unpert[1];
    /* polarizability */
    QInt ALPHA_NUM_PROPS=1;
    QInt ALPHA_LEN_TUPLE[]={2};
    QcPertInt ALPHA_PERT_TUPLE[]={PERT_DIPOLE,PERT_DIPOLE};
    QInt ALPHA_NUM_FREQ_CONFIGS[]={1};
    QReal ALPHA_PERT_FREQS[]={-0.072,0.0,0.072,0.0};
    QInt ALPHA_KN_RULES[]={0};
    /* response functions */
    QInt size_rsp_funs;
    QReal rsp_funs[18];
    QInt ipert,jpert,ival;
    /* error information */
    QErrorCode ierr;

    ///* sets the equation of motion of electrons */
    //ierr = OpenRSPSetElecEOM(open_rsp, ELEC_AO_D_MATRIX);
    //QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetElecEOM");
    //fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPSetElecEOM() passed\n");

    /* sets the context of linear response equation solver */
    ierr = OpenRSPSetLinearRSPSolver(open_rsp,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     (void *)SOLVER_CONTEXT,
#endif
                                     &get_linear_rsp_solution);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetLinearRSPSolver()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPSetLinearRSPSolver() passed\n");

    /* sets the context of overlap operator */
    ierr = OpenRSPSetOverlap(open_rsp,
                             overlap_num_pert,
                             overlap_pert_labels,
                             overlap_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (void *)OVERLAP_CONTEXT,
#endif
                             &get_overlap_mat,
                             &get_overlap_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetOverlap()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPSetOverlap() passed\n");

    /* adds one-electron Hamiltonian */
    ierr = OpenRSPAddOneOper(open_rsp,
                             oneham_num_pert,
                             oneham_pert_labels,
                             oneham_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (void *)ONEHAM_CONTEXT,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(h)");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPAddOneOper(h) passed\n");

    /* adds external field */
    ierr = OpenRSPAddOneOper(open_rsp,
                             ext_field_num_pert,
                             ext_field_pert_labels,
                             ext_field_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (void *)EXT_FIELD_CONTEXT,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(V)");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPAddOneOper(V) passed\n");

    /* adds two-electron operator */
    ierr = OpenRSPAddTwoOper(open_rsp,
                             two_oper_num_pert,
                             two_oper_pert_labels,
                             two_oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (void *)TWO_OPER_CONTEXT,
#endif
                             &get_two_oper_mat,
                             &get_two_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddTwoOper()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPAddTwoOper() passed\n");

    /* assembles the context of response theory calculations */
    ierr = OpenRSPAssemble(open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAssemble()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPAssemble() passed\n");

    /* writes the context of response theory calculations */
    ierr = OpenRSPWrite(open_rsp, fp_log);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPWrite()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPWrite() passed\n");

    /* sets the unperturbed AO-based Fock matrix */
    F_unpert[0] = (QcMat *)malloc(sizeof(QcMat));
    if (F_unpert[0]==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for Fock matrix");
    }
    ierr = OpenRSPTestReadMat(DMAT_HF_FOCK, 0, 1, F_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPTestReadMat(F)");
    /* sets the unperturbed AO-based density matrix */
    D_unpert[0] = (QcMat *)malloc(sizeof(QcMat));
    if (D_unpert[0]==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for density matrix");
    }
    ierr = OpenRSPTestReadMat(DMAT_HF_DENSITY, 0, 1, D_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPTestReadMat(D)");
    /* sets the unperturbed AO-based overlap integral matrix */
    S_unpert[0] = (QcMat *)malloc(sizeof(QcMat));
    if (S_unpert[0]==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for overlap integral matrix");
    }
    ierr = OpenRSPTestReadMat(DMAT_OVERLAP, 0, 1, S_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPTestReadMat(S)");

    /* gets the polarizability */
    size_rsp_funs = 9;
    ierr = OpenRSPGetRSPFun(open_rsp,
                            F_unpert[0],
                            D_unpert[0],
                            S_unpert[0],
                            ALPHA_NUM_PROPS,
                            ALPHA_LEN_TUPLE,
                            ALPHA_PERT_TUPLE,
                            ALPHA_NUM_FREQ_CONFIGS,
                            ALPHA_PERT_FREQS,
                            ALPHA_KN_RULES,
                            0,
                            1e-8,
                            size_rsp_funs,
                            rsp_funs);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPGetRSPFun()");
    fprintf(fp_log, "OpenRSPDMatTest>> OpenRSPGetRSPFun() passed\n");
    for (ipert=0,ival=0; ipert<3; ipert++) {
        for (jpert=0; jpert<3; jpert++) {
            fprintf(fp_log,
                    " (%"QCREAL_FMT",%"QCREAL_FMT")",
                    rsp_funs[ival], rsp_funs[ival+1]);
            ival += 2;
        }
        fprintf(fp_log, "\n");
    }

    /* cleans */
    ierr = QcMatDestroy(F_unpert[0]);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(F)");
    free(F_unpert[0]);
    F_unpert[0] = NULL;
    ierr = QcMatDestroy(D_unpert[0]);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(D)");
    free(D_unpert[0]);
    D_unpert[0] = NULL;
    ierr = QcMatDestroy(S_unpert[0]);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(S)");
    free(S_unpert[0]);
    S_unpert[0] = NULL;

    return QSUCCESS;
}

