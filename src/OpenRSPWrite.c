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

   This file implements the function OpenRSPWrite().

   2014-07-30, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief writes the context of response theory calculations
     \author Bin Gao
     \date 2014-07-30
     \param[OneRSP:struct]{in} open_rsp the context of response theory calculations
     \param[QChar:char]{in} file_name the name of the file
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPWrite(const OpenRSP *open_rsp, const QChar *file_name)
{
    FILE *fp_rsp;      /* file pointer */
    QInt ipert,isize;  /* incremental recorders */
    QErrorCode ierr;   /* error information */
    /* opens the file */
    fp_rsp = fopen(file_name, "a");
    if (fp_rsp==NULL) {
        printf("OpenRSPWrite>> file: %s\n", file_name);
        QErrorExit(FILE_AND_LINE, "failed to open the file in appending mode");
    }
    fprintf(fp_rsp, "\nOpenRSP library compiled at %s, %s\n", __TIME__, __DATE__);
    /* context of the EOM of electrons */
    /*FIXME: ierr = xxWrite(open_rsp->elec_eom); */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverWrite(open_rsp->rsp_solver, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverWrite");
    }
#if defined(OPENRSP_PERTURBATION_FREE)
    fprintf(fp_rsp,
            "OpenRSPWrite>> number of all perturbations involved in calculations %"QINT_FMT"\n",
            open_rsp->num_pert);
    fprintf(fp_rsp,
            "OpenRSPWrite>> perturbation    maximum-order    sizes-up-to-maximum-order\n");
    for (ipert=0; ipert<open_rsp->num_pert; ipert++) {
        fprintf(fp_rsp,
                "OpenRSPWrite>>  %"QINT_FMT"               %"QINT_FMT"               ",
                open_rsp->perturbations[ipert],
                open_rsp->pert_max_orders[ipert]);
        for (isize=open_rsp->size_ptr[ipert]; isize<open_rsp->size_ptr[ipert+1]; isize++) {
            fprintf(fp_rsp, " %"QINT_FMT"", open_rsp->pert_sizes[isize]);
        }
        fprintf(fp_rsp, "\n");
    }
#endif
    if (open_rsp->overlap!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> overlap integrals\n");
        ierr = RSPOverlapWrite(open_rsp->overlap, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapWrite");
    }
    if (open_rsp->one_oper!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> linked list of one-electron operators\n");
        ierr = RSPOneOperWrite(open_rsp->one_oper, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperWrite");
    }
    if (open_rsp->two_oper!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> linked list of two-electron operators\n");
        ierr = RSPTwoOperWrite(open_rsp->two_oper, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperWrite");
    }
    if (open_rsp->nuc_contrib!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> nuclear contributions\n");
        ierr = RSPNucContribWrite(open_rsp->nuc_contrib, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucContribWrite");
    }
    /* closes the file */
    fclose(fp_rsp);
    return QSUCCESS;
}
