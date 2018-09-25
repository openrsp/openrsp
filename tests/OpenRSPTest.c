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


   This file tests the APIs of OpenRSP library.

   2014-07-31, Bin Gao
   * first version
*/

#include "OpenRSPTest.h"

/* <macrodef name='OPENRSP_TEST_EXECUTABLE'>
     Build test suite as excutables.
   </macrodef> */
#if defined(OPENRSP_TEST_EXECUTABLE)
QErrorCode main()
{
    FILE *fp_log=stdout;  /* file pointer */
#else
QErrorCode OpenRSPTest(FILE *fp_log)
{
#endif
#include "OpenRSPTestPerturbations.h"
    OpenRSP open_rsp;  /* context of response theory calculations */
    QErrorCode ierr;   /* error information */

    /* creates the context of response theory calculations */
    ierr = OpenRSPCreate(&open_rsp, NUM_ATOMS);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPCreate()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPCreate() passed\n");

    /* sets information of all perturbations */
    ierr = OpenRSPSetPerturbations(&open_rsp,
                                   NUM_ALL_PERT,
                                   ALL_PERT_LABELS,
                                   ALL_PERT_MAX_ORDERS,
                                   ALL_PERT_SIZES,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   (void *)PERT_CONTEXT,
#endif
                                   &get_pert_concatenation);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetPerturbations()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPSetPerturbations() passed\n");

    /* adds the nuclear Hamiltonian */
    ierr = OpenRSPAddZeroOper(&open_rsp,
                              NUM_ALL_PERT,
                              ALL_PERT_LABELS,
                              ALL_PERT_MAX_ORDERS,
#if defined(OPENRSP_C_USER_CONTEXT)
                              (void *)ZERO_OPER_CONTEXT,
#endif
                              &get_zero_oper_contrib);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddZeroOper()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPAddZeroOper() passed\n");

    /* tests the density matrix-based response theory */
    ierr = OpenRSPDMatTest(&open_rsp, fp_log);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDMatTest()");
    fprintf(fp_log, "OpenRSPTest>> density matrix-based response theory passed\n");

    /* destroys the context of response theory calculations */
    ierr = OpenRSPDestroy(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDestroy()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPDestroy() passed\n");

    return QSUCCESS;
}
