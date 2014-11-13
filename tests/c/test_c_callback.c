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

   This file implements the callback functions for the C tests.

   2014-07-31, Bin Gao:
   * first version
*/

/* callback functions */
#include "tests/openrsp_c_test_callback.h"

QVoid get_rsp_solution(const QMat *ref_ham,
                       const QMat *ref_state,
                       const QMat *ref_overlap,
                       const QInt num_freq_sums,
                       const QReal *freq_sums,
                       const QInt size_pert,
                       const QInt num_RHS,
                       QMat *RHS_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       QMat *rsp_param[])
{
}

#if defined(OPENRSP_PERTURBATION_FREE)
QVoid get_pert_comp(const QInt perturbation,
                    const QInt pert_order,
                    const QInt pert_rank,
#if defined(OPENRSP_C_USER_CONTEXT)
                    QVoid *user_ctx,
#endif
                    QInt *pert_num_comp,
                    QInt *pert_components,
                    QInt *pert_comp_orders)
{
}

QVoid get_pert_rank(const QInt perturbation,
                    const QInt pert_num_comp,
                    const QInt *pert_components,
                    const QInt *pert_comp_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                    QVoid *user_ctx,
#endif
                    QInt *pert_rank)
{
}
#endif

QVoid get_overlap_mat(const QInt bra_num_pert,
                      const QInt *bra_perturbations,
                      const QInt *bra_pert_orders,
                      const QInt ket_num_pert,
                      const QInt *ket_perturbations,
                      const QInt *ket_pert_orders,
                      const QInt num_pert,
                      const QInt *perturbations,
                      const QInt *pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                      QVoid *user_ctx,
#endif
                      const QInt num_int,
                      QMat *val_int[])
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_lab;
#endif
    QInt ipert;
    QInt imat;
    QErrorCode ierr;
#if defined(ZERO_BASED_NUMBERING)
    const QInt idx_block_row=0;
    const QInt idx_block_col=0;
    const QInt idx_first_row=0;
    const QInt idx_first_col=0;
#else
    const QInt idx_block_row=1;
    const QInt idx_block_col=1;
    const QInt idx_first_row=1;
    const QInt idx_first_col=1;
#endif
    const QInt num_row_set=2;
    const QInt num_col_set=2;
    QReal values_real[4];
    QReal values_imag[4];
    printf("get_overlap_mat>> bra_num_pert %"QINT_FMT"\n", bra_num_pert);
    for (ipert=0; ipert<bra_num_pert; ipert++) {
        printf("get_overlap_mat>> bra_pert. %"QINT_FMT" %"QINT_FMT"\n",
               bra_perturbations[ipert],
               bra_pert_orders[ipert]);
    }
    printf("get_overlap_mat>> ket_num_pert %"QINT_FMT"\n", ket_num_pert);
    for (ipert=0; ipert<ket_num_pert; ipert++) {
        printf("get_overlap_mat>> ket_pert. %"QINT_FMT" %"QINT_FMT"\n",
               ket_perturbations[ipert],
               ket_pert_orders[ipert]);
    }
    printf("get_overlap_mat>> num_pert %"QINT_FMT"\n", num_pert);
    for (ipert=0; ipert<num_pert; ipert++) {
        printf("get_overlap_mat>> pert. %"QINT_FMT" %"QINT_FMT"\n",
               perturbations[ipert],
               pert_orders[ipert]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_lab = (QChar *)user_ctx;
    printf("get_overlap_mat>> label %s\n", overlap_lab);
    if (strcmp(overlap_lab, "OVERLAP")==0) {
        printf("get_overlap_mat>> overlap integrals\n");
        values_real[0]=0.1;values_real[1]=0.2;values_real[2]=0.3;values_real[3]=0.4;
        values_imag[0]=1.1;values_imag[1]=1.2;values_imag[2]=1.3;values_imag[3]=1.4;
    }
    else {
        printf("get_overlap_mat>> unknown one-electron operator\n");
        exit(1);
    }
#endif
    printf("get_overlap_mat>> num_int %"QINT_FMT"\n", num_int);
    for (imat=0; imat<num_int; imat++) {
        ierr = QMatSetValues(val_int[imat],
                             idx_block_row,
                             idx_block_col,
                             idx_first_row,
                             num_row_set,
                             idx_first_col,
                             num_col_set,
                             values_real,
                             values_imag);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QMatSetValues");
            exit(ierr);
        }
    }
}

QVoid get_overlap_exp(const QInt bra_num_pert,
                      const QInt *bra_perturbations,
                      const QInt *bra_pert_orders,
                      const QInt ket_num_pert,
                      const QInt *ket_perturbations,
                      const QInt *ket_pert_orders,
                      const QInt num_pert,
                      const QInt *perturbations,
                      const QInt *pert_orders,
                      const QInt num_dens,
                      QMat *ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                      QVoid *user_ctx,
#endif
                      const QInt num_exp,
                      QReal *val_exp)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_lab;
#endif
    QInt ipert;
    printf("get_overlap_exp>> bra_num_pert %"QINT_FMT"\n", bra_num_pert);
    for (ipert=0; ipert<bra_num_pert; ipert++) {
        printf("get_overlap_exp>> bra_pert. %"QINT_FMT" %"QINT_FMT"\n",
               bra_perturbations[ipert],
               bra_pert_orders[ipert]);
    }
    printf("get_overlap_exp>> ket_num_pert %"QINT_FMT"\n", ket_num_pert);
    for (ipert=0; ipert<ket_num_pert; ipert++) {
        printf("get_overlap_exp>> ket_pert. %"QINT_FMT" %"QINT_FMT"\n",
               ket_perturbations[ipert],
               ket_pert_orders[ipert]);
    }
    printf("get_overlap_exp>> num_pert %"QINT_FMT"\n", num_pert);
    for (ipert=0; ipert<num_pert; ipert++) {
        printf("get_overlap_exp>> pert. %"QINT_FMT" %"QINT_FMT"\n",
               perturbations[ipert],
               pert_orders[ipert]);
    }
    printf("get_overlap_exp>> num_dens %"QINT_FMT"\n", num_dens);
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_lab = (QChar *)user_ctx;
    printf("get_overlap_exp>> label %s\n", overlap_lab);
    if (strcmp(overlap_lab, "OVERLAP")==0) {
        printf("get_overlap_exp>> overlap integrals\n");
    }
    else {
        printf("get_overlap_exp>> unknown one-electron operator\n");
        exit(1);
    }
#endif
}

QVoid get_one_oper_mat(const QInt num_pert,
                       const QInt *perturbations,
                       const QInt *pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_int,
                       QMat *val_int[])
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_lab = "ONEHAM";
    QChar *ext_field_lab = "EXT_FIELD";
    QChar *one_oper_lab;
#endif
    QInt ipert;
    QInt imat;
    QErrorCode ierr;
#if defined(ZERO_BASED_NUMBERING)
    const QInt idx_block_row=0;
    const QInt idx_block_col=0;
    const QInt idx_first_row=0;
    const QInt idx_first_col=0;
#else
    const QInt idx_block_row=1;
    const QInt idx_block_col=1;
    const QInt idx_first_row=1;
    const QInt idx_first_col=1;
#endif
    const QInt num_row_set=2;
    const QInt num_col_set=2;
    QReal values_real[4];
    QReal values_imag[4];
    printf("get_one_oper_mat>> num_pert %"QINT_FMT"\n", num_pert);
    for (ipert=0; ipert<num_pert; ipert++) {
        printf("get_one_oper_mat>> pert. %"QINT_FMT" %"QINT_FMT"\n",
               perturbations[ipert],
               pert_orders[ipert]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_lab = (QChar *)user_ctx;
    printf("get_one_oper_mat>> label %s\n", one_oper_lab);
    if (strcmp(one_oper_lab, oneham_lab)==0) {
        printf("get_one_oper_mat>> one-electron Hamiltonian\n");
        values_real[0]=0.1;values_real[1]=0.2;values_real[2]=0.3;values_real[3]=0.4;
        values_imag[0]=1.1;values_imag[1]=1.2;values_imag[2]=1.3;values_imag[3]=1.4;
    }
    else if (strcmp(one_oper_lab, ext_field_lab)==0) {
        printf("get_one_oper_mat>> external field\n");
        values_real[0]=2.1;values_real[1]=2.2;values_real[2]=2.3;values_real[3]=2.4;
        values_imag[0]=3.1;values_imag[1]=3.2;values_imag[2]=3.3;values_imag[3]=3.4;
    }
    else {
        printf("get_one_oper_mat>> unknown one-electron operator\n");
        exit(1);
    }
#endif
    printf("get_one_oper_mat>> num_int %"QINT_FMT"\n", num_int);
    for (imat=0; imat<num_int; imat++) {
        ierr = QMatSetValues(val_int[imat],
                             idx_block_row,
                             idx_block_col,
                             idx_first_row,
                             num_row_set,
                             idx_first_col,
                             num_col_set,
                             values_real,
                             values_imag);
        if (ierr!=QSUCCESS) {
            printf("get_one_oper_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QMatSetValues");
            exit(ierr);
        }
    }
}

QVoid get_one_oper_exp(const QInt num_pert,
                       const QInt *perturbations,
                       const QInt *pert_orders,
                       const QInt num_dens,
                       QMat *ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_exp,
                       QReal *val_exp)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_lab = "ONEHAM";
    QChar *ext_field_lab = "EXT_FIELD";
    QChar *one_oper_lab;
#endif
    QInt ipert;
    printf("get_one_oper_exp>> num_pert %"QINT_FMT"\n", num_pert);
    for (ipert=0; ipert<num_pert; ipert++) {
        printf("get_one_oper_exp>> pert. %"QINT_FMT" %"QINT_FMT"\n",
               perturbations[ipert],
               pert_orders[ipert]);
    }
    printf("get_one_oper_exp>> num_dens %"QINT_FMT"\n", num_dens);
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_lab = (QChar *)user_ctx;
    printf("get_one_oper_exp>> label %s\n", one_oper_lab);
    if (strcmp(one_oper_lab, oneham_lab)==0) {
        printf("get_one_oper_exp>> one-electron Hamiltonian\n");
    }
    else if (strcmp(one_oper_lab, ext_field_lab)==0) {
        printf("get_one_oper_exp>> external field\n");
    }
    else {
        printf("get_one_oper_exp>> unknown one-electron operator\n");
        exit(1);
    }
#endif
}
