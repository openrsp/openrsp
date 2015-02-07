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

   This file tests the APIs of OpenRSP library.

   2014-07-31, Bin Gao:
   * first version
*/

/* header file of OpenRSP library */
#include "openrsp.h"
/* parameters for the test suite */
#include "tests/openrsp_c_test.h"
/* callback functions */
#include "tests/openrsp_c_test_callback.h"

#if defined(OPENRSP_TEST_EXECUTABLE)
QErrorCode main()
{
    FILE *fp_log=stdout;  /* file pointer */
#else
QErrorCode test_c_OpenRSP(FILE *fp_log)
{
#endif
    OpenRSP open_rsp;                                 /* context of response theory calculations */
    QChar *solver_lab = "SOLVER";
    const QInt ALL_PERT_LABELS[NUM_ALL_PERT] = {      /* labels of all perturbations */
        PERT_DIPOLE,PERT_MAGNETIC,PERT_GEOMETRIC};
    const QInt ALL_PERT_MAX_ORDERS[NUM_ALL_PERT] = {  /* maximum allowed orders of all perturbations */
        MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC};
    const QInt ALL_PERT_SIZES[] = {                   /* sizes of all perturbations up to their maximum orders */
        3,                                            /* electric dipole */
        3,6,10,15,21,28,36,                           /* magnetic derivatives */
        36,666,8436,82251,465552,1898232,6257232};    /* geometric derivatives (12 atoms) */
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *pert_lab = "PERT";
#endif
    /* overlap integrals with London atomic orbitals */
    QInt overlap_num_pert = 2;
    QInt overlap_pert_labels[2] = {PERT_MAGNETIC,PERT_GEOMETRIC};
    QInt overlap_pert_orders[2] = {MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_lab = "OVERLAP";
#endif
    /* one-electron Hamiltonian */
    QInt oneham_num_pert = 2;
    QInt oneham_pert_labels[2] = {PERT_MAGNETIC,PERT_GEOMETRIC};
    QInt oneham_pert_orders[2] = {MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_lab = "ONEHAM";
#endif
    /* external field */
    QInt ext_field_num_pert = 3;
    QInt ext_field_pert_labels[3] = {PERT_DIPOLE,PERT_MAGNETIC,PERT_GEOMETRIC};
    QInt ext_field_pert_orders[3] = {MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *ext_field_lab = "EXT_FIELD";
#endif
    /* atoms and origins */
    QInt num_atoms = 2;
    QReal atom_coord[] = {0.0,0.0,0.0, 1.0,1.0,1.0};
    QReal atom_charge[] = {1.0, 2.0};
    QReal dipole_origin[3] = {0.1,0.1,0.1};
    QReal gauge_origin[3] = {0.2,0.2,0.2};
    /* referenced state */
    QMat ref_ham;
    QMat ref_state;
    QMat ref_overlap;
/*FIXME: to move to test_... */
    QInt alpha_num_props = 1;
    QInt alpha_num_pert[1] = {2};
    QInt alpha_pert_labels[2] = {PERT_DIPOLE,PERT_DIPOLE};
    QInt alpha_num_freqs[2] = {1,1};
    QReal alpha_pert_freqs[4] = {-0.072,0.0,0.072,0.0};
    QInt alpha_kn_rules[1] = {0};
    QInt size_rsp_funs = 9;
    QReal rsp_funs[9];

    QErrorCode ierr;   /* error information */

    ierr = OpenRSPCreate(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPCreate");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPCreate() passed\n");

    ierr = OpenRSPSetSolver(&open_rsp,
#if defined(OPENRSP_C_USER_CONTEXT)
                            (QVoid *)solver_lab,
#endif
                            &get_rsp_solution);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetSolver");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetSolver() passed\n");

#if defined(OPENRSP_PERTURBATION_FREE)
    ierr = OpenRSPSetPerturbations(&open_rsp,
                                   NUM_ALL_PERT,
                                   ALL_PERT_LABELS,
                                   ALL_PERT_MAX_ORDERS,
                                   ALL_PERT_SIZES,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   (QVoid *)pert_lab,
#endif
                                   &get_pert_comp,
                                   &get_pert_rank);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetPerturbations");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetPerturbations() passed\n");
#endif

    ierr = OpenRSPSetPDBS(&open_rsp,
                          overlap_num_pert,
                          overlap_pert_labels,
                          overlap_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          (QVoid *)overlap_lab,
#endif
                          &get_overlap_mat,
                          &get_overlap_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetPDBS");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetPDBS() passed\n");

    ierr = OpenRSPAddOneOper(&open_rsp,
                             oneham_num_pert,
                             oneham_pert_labels,
                             oneham_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)oneham_lab,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(h)");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPAddOneOper(h) passed\n");

    ierr = OpenRSPAddOneOper(&open_rsp,
                             ext_field_num_pert,
                             ext_field_pert_labels,
                             ext_field_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)ext_field_lab,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(V)");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPAddOneOper(V) passed\n");

    ierr = OpenRSPSetAtoms(&open_rsp,
                           num_atoms,
                           atom_coord,
                           atom_charge);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetAtoms");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetAtoms() passed\n");

    ierr = OpenRSPSetDipoleOrigin(&open_rsp, dipole_origin);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetDipoleOrigin");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetDipoleOrigin() passed\n");

    ierr = OpenRSPSetGaugeOrigin(&open_rsp, gauge_origin);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetGaugeOrigin");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetGaugeOrigin() passed\n");

    ierr = OpenRSPAssemble(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAssemble");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPAssemble() passed\n");

    ierr = OpenRSPWrite(&open_rsp, OPENRSP_C_LOG);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPWrite");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPWrite() passed\n");

    ierr = OpenRSPGetRSPFun(&open_rsp,
                            &ref_ham,
                            &ref_state,
                            &ref_overlap,
                            alpha_num_props,
                            alpha_num_pert,
                            alpha_pert_labels,
                            alpha_num_freqs,
                            alpha_pert_freqs,
                            alpha_kn_rules,
                            size_rsp_funs,
                            rsp_funs);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPGetRSPFun");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPGetRSPFun() passed\n");

    ierr = OpenRSPDestroy(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDestroy");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPDestroy() passed\n");

    return QSUCCESS;
}
