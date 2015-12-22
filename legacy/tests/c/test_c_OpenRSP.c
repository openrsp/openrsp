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
#include "tests/openrsp_test_param.h"
/* defined perturbations and their maximum orders */
#include "tests/openrsp_c_perturbations.h"
/* atoms and origins */
#include "tests/openrsp_c_molecule.h"
#if defined(OPENRSP_PERTURBATION_FREE)
/* callback functions for perturbations */
#include "tests/openrsp_c_pert_callback.h"
#endif

/* declaration of different test routines */
QErrorCode test_c_OpenRSP_AO(OpenRSP *open_rsp, FILE *fp_log);

#if defined(OPENRSP_TEST_EXECUTABLE)
QErrorCode main()
{
    FILE *fp_log=stdout;  /* file pointer */
#else
QErrorCode test_c_OpenRSP(FILE *fp_log)
{
#endif
    OpenRSP open_rsp;                                 /* context of response theory calculations */
#if defined(OPENRSP_PERTURBATION_FREE)
    const QInt ALL_PERT_LABELS[NUM_ALL_PERT] = {      /* labels of all perturbations */
        PERT_GEOMETRIC,PERT_DIPOLE,PERT_MAGNETIC};
    const QInt ALL_PERT_MAX_ORDERS[NUM_ALL_PERT] = {  /* maximum allowed orders of all perturbations */
        MAX_ORDER_GEOMETRIC,MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC};
    const QInt ALL_PERT_SIZES[] = {                   /* sizes of all perturbations up to their maximum orders */
        12,78,364,1365,4368,12376,31824,              /* geometric derivatives (4 atoms) */
        3,                                            /* electric dipole */
        3,6,10,15,21,28,36};                          /* magnetic derivatives */
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *pert_context = "NRNZGEO";                  /* user defined context for perturbations */
#endif
#endif
    QErrorCode ierr;                                  /* error information */

    /* creates the context of response theory calculations */
    ierr = OpenRSPCreate(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPCreate");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPCreate() passed\n");

#if defined(OPENRSP_PERTURBATION_FREE)
    /* sets information of all perturbations */
    ierr = OpenRSPSetPerturbations(&open_rsp,
                                   NUM_ALL_PERT,
                                   ALL_PERT_LABELS,
                                   ALL_PERT_MAX_ORDERS,
                                   ALL_PERT_SIZES,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   (QVoid *)pert_context,
#endif
                                   &get_pert_comp,
                                   &get_pert_rank);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetPerturbations");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetPerturbations() passed\n");
#endif

    /* sets the nuclear contributions */
    ierr = OpenRSPSetNucContributions(&open_rsp,
                                      NUM_ALL_PERT,
                                      ALL_PERT_LABELS,
                                      ALL_PERT_MAX_ORDERS,
#if defined(OPENRSP_C_USER_CONTEXT)   
                                      (QVoid *)pert_context,
#endif
                                      &get_nuc_contrib,
                                      NUM_ATOMS);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetNucContributions");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPSetNucContributions() passed\n");

    /* tests the density matrix-based response theory */
    ierr = test_c_OpenRSP_AO(&open_rsp, fp_log);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling test_c_OpenRSP_AO");
    fprintf(fp_log, "test_c_OpenRSP>> density matrix-based response theory passed\n");

    /* destroys the context of response theory calculations */
    ierr = OpenRSPDestroy(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDestroy");
    fprintf(fp_log, "test_c_OpenRSP>> OpenRSPDestroy() passed\n");

    return QSUCCESS;
}
