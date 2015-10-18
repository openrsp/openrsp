/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.


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
    ierr = OpenRSPCreate(&open_rsp);
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

    /* sets the nuclear Hamiltonian */
    ierr = OpenRSPSetNucHamilton(&open_rsp,
                                 NUM_ALL_PERT,
                                 ALL_PERT_LABELS,
                                 ALL_PERT_MAX_ORDERS,
#if defined(OPENRSP_C_USER_CONTEXT)   
                                 (void *)NUC_HAMILTON_CONTEXT,
#endif
                                 &get_nuc_contrib,
                                 NUM_ATOMS);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPSetNucHamilton()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPSetNucHamilton() passed\n");

    /* tests the density matrix-based response theory */
    ierr = OpenRSPDensAOTest(&open_rsp, fp_log);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDensAOTest()");
    fprintf(fp_log, "OpenRSPTest>> density matrix-based response theory passed\n");

    /* destroys the context of response theory calculations */
    ierr = OpenRSPDestroy(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDestroy()");
    fprintf(fp_log, "OpenRSPTest>> OpenRSPDestroy() passed\n");

    return QSUCCESS;
}
