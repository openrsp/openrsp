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

   This file tests the atomic orbital density matrix-based response theory.

   2014-07-31, Bin Gao:
   * first version
*/

/* header file of OpenRSP library */
#include "openrsp.h"
/* parameters for the test suite */
#include "tests/openrsp_test_param.h"
/* defined perturbations and their maximum orders */
#include "tests/openrsp_c_perturbations.h"
/* callback functions */
#include "tests/ao_dens/openrsp_c_callback.h"

QErrorCode test_c_OpenRSP_AO(OpenRSP *open_rsp, FILE *fp_log)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    /* user defined context of linear response equation solver */
    QChar *solver_context = "SOLVER";
#endif
    /* overlap integrals with London atomic orbitals */
    QInt overlap_num_pert = 2;
    QInt overlap_pert_labels[2] = {PERT_GEOMETRIC,PERT_MAGNETIC};
    QInt overlap_pert_orders[2] = {MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_context = "OVERLAP";
#endif
    /* one-electron Hamiltonian */
    QInt oneham_num_pert = 2;
    QInt oneham_pert_labels[2] = {PERT_GEOMETRIC,PERT_MAGNETIC};
    QInt oneham_pert_orders[2] = {MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_context = "ONEHAM";
#endif
    /* external field */
    QInt ext_field_num_pert = 3;
    QInt ext_field_pert_labels[3] = {PERT_GEOMETRIC,PERT_DIPOLE,PERT_MAGNETIC};
    QInt ext_field_pert_orders[3] = {MAX_ORDER_GEOMETRIC,MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *ext_field_context = "EXT_FIELD";
#endif
    /* two-electron Hamiltonian */
    QInt twoel_num_pert = 2;
    QInt twoel_pert_labels[2] = {PERT_GEOMETRIC,PERT_MAGNETIC};
    QInt twoel_pert_orders[2] = {MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *twoel_context = "NONLAO";
#endif
    /* referenced state */
#include "tests/ao_dens/openrsp_c_ao_dims.h"
#include "tests/ao_dens/openrsp_c_ao_fock.h"
#include "tests/ao_dens/openrsp_c_ao_density.h"
#include "tests/ao_dens/openrsp_c_ao_overlap.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QcMat F_unpert;
    QcMat D_unpert;
    QcMat S_unpert;
    /* polarizability */
    QInt ALPHA_NUM_PROPS = 1;
    QInt ALPHA_LEN_TUPLE[1] = {2};
    QInt ALPHA_PERT_TUPLE[2] = {PERT_DIPOLE,PERT_DIPOLE};
    QInt ALPHA_NUM_FREQ_CONFIGS[1] = {1};
    QReal ALPHA_PERT_FREQS[4] = {-0.072,0.0,0.072,0.0};
    QInt ALPHA_KN_RULES[1] = {0};
    /* response functions */
    QInt size_rsp_funs;
    QReal rsp_funs[18];
    QInt ipert,jpert,ival;
    /* error information */
    QErrorCode ierr;

    /* sets the equation of motion of electrons */
    ierr = OpenRSPSetElecEOM(open_rsp, ELEC_AO_D_MATRIX);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetElecEOM");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPSetElecEOM() passed\n");

    /* sets the context of linear response equation solver */
    ierr = OpenRSPSetLinearRSPSolver(open_rsp,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     (QVoid *)solver_context,
#endif
                                     &get_linear_rsp_solution);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetLinearRSPSolver");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPSetLinearRSPSolver() passed\n");

    /* sets the context of perturbation dependent basis sets */
    ierr = OpenRSPSetPDBS(open_rsp,
                          overlap_num_pert,
                          overlap_pert_labels,
                          overlap_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          (QVoid *)overlap_context,
#endif
                          &get_overlap_mat,
                          &get_overlap_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetPDBS");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPSetPDBS() passed\n");

    /* adds one-electron Hamiltonian */
    ierr = OpenRSPAddOneOper(open_rsp,
                             oneham_num_pert,
                             oneham_pert_labels,
                             oneham_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)oneham_context,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(h)");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPAddOneOper(h) passed\n");

    /* adds external field */
    ierr = OpenRSPAddOneOper(open_rsp,
                             ext_field_num_pert,
                             ext_field_pert_labels,
                             ext_field_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)ext_field_context,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(V)");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPAddOneOper(V) passed\n");

    /* adds two-electron Hamiltonian */
    ierr = OpenRSPAddTwoOper(open_rsp,
                             twoel_num_pert,
                             twoel_pert_labels,
                             twoel_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)twoel_context,
#endif
                             &get_two_oper_mat,
                             &get_two_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddTwoOper()");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPAddTwoOper() passed\n");

    /* assembles the context of response theory calculations */
    ierr = OpenRSPAssemble(open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAssemble");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPAssemble() passed\n");

    /* writes the context of response theory calculations */
    ierr = OpenRSPWrite(open_rsp, OPENRSP_C_LOG);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPWrite");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPWrite() passed\n");

    /* sets the unperturbed Fock matrix */
    ierr = QcMatCreate(&F_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatCreate(F)");
    ierr = QcMatBlockCreate(&F_unpert, 1);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate(F)");
    ierr = QcMatSetSymType(&F_unpert, QSYMMAT);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetSymType(F)");
    ierr = QcMatSetDataType(&F_unpert, 1, idx_block_row, idx_block_col, data_type);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType(F)");
    ierr = QcMatSetDimMat(&F_unpert, NUM_AO, NUM_AO);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat(F)");
    ierr = QcMatAssemble(&F_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble(F)");
    ierr = QcMatSetValues(&F_unpert,
                          IDX_BLOCK_ROW,
                          IDX_BLOCK_COL,
                          IDX_FIRST_ROW,
                          NUM_AO,
                          IDX_FIRST_COL,
                          NUM_AO,
                          values_fock,
                          NULL);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetValues(F)");
    /* sets the unperturbed density matrix */
    ierr = QcMatCreate(&D_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatCreate(D)");
    ierr = QcMatBlockCreate(&D_unpert, 1);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate(D)");
    ierr = QcMatSetSymType(&D_unpert, QSYMMAT);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetSymType(D)");
    ierr = QcMatSetDataType(&D_unpert, 1, idx_block_row, idx_block_col, data_type);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType(D)");
    ierr = QcMatSetDimMat(&D_unpert, NUM_AO, NUM_AO);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat(D)");
    ierr = QcMatAssemble(&D_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble(D)");
    ierr = QcMatSetValues(&D_unpert,
                          IDX_BLOCK_ROW,
                          IDX_BLOCK_COL,
                          IDX_FIRST_ROW,
                          NUM_AO,
                          IDX_FIRST_COL,
                          NUM_AO,
                          values_density,
                          NULL);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetValues(D)");
    /* sets the unperturbed overlap integrals */
    ierr = QcMatCreate(&S_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatCreate(S)");
    ierr = QcMatBlockCreate(&S_unpert, 1);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate(S)");
    ierr = QcMatSetSymType(&S_unpert, QSYMMAT);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetSymType(S)");
    ierr = QcMatSetDataType(&S_unpert, 1, idx_block_row, idx_block_col, data_type);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType(S)");
    ierr = QcMatSetDimMat(&S_unpert, NUM_AO, NUM_AO);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat(S)");
    ierr = QcMatAssemble(&S_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble(S)");
    ierr = QcMatSetValues(&S_unpert,
                          IDX_BLOCK_ROW,
                          IDX_BLOCK_COL,
                          IDX_FIRST_ROW,
                          NUM_AO,
                          IDX_FIRST_COL,
                          NUM_AO,
                          values_overlap,
                          NULL);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetValues(S)");

    /* gets the polarizability */
    size_rsp_funs = 9;
    ierr = OpenRSPGetRSPFun(open_rsp,
                            &F_unpert,
                            &D_unpert,
                            &S_unpert,
                            ALPHA_NUM_PROPS,
                            ALPHA_LEN_TUPLE,
                            ALPHA_PERT_TUPLE,
                            ALPHA_NUM_FREQ_CONFIGS,
                            ALPHA_PERT_FREQS,
                            ALPHA_KN_RULES,
                            size_rsp_funs,
                            rsp_funs);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPGetRSPFun");
    fprintf(fp_log, "test_c_OpenRSP_AO>> OpenRSPGetRSPFun() passed\n");
    for (ipert=0,ival=0; ipert<3; ipert++) {
        for (jpert=0; jpert<3; jpert++) {
            fprintf(fp_log, " (%f,%f)", rsp_funs[ival], rsp_funs[ival+1]);
            ival += 2;
        }
        fprintf(fp_log, "\n");
    }

    /* cleans */
    ierr = QcMatDestroy(&F_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(F)");
    ierr = QcMatDestroy(&D_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(D)");
    ierr = QcMatDestroy(&S_unpert);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatDestroy(S)");

    return QSUCCESS;
}
