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


   2014-07-31, Bin Gao
   * first version
*/

#include "OpenRSP.h"

QErrorCode OpenRSPCreateFortranAdapter(void **open_rsp, const QInt num_atoms)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)malloc(sizeof(OpenRSP));
    if (c_open_rsp==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for c_open_rsp");
    }
    ierr = OpenRSPCreate(c_open_rsp, num_atoms);
    *open_rsp = (void *)(c_open_rsp);
    return ierr;
}

//QErrorCode f_api_OpenRSPSetElecEOM(void **open_rsp,
//                                   const QInt elec_EOM_type)
//{
//    OpenRSP *c_open_rsp;
//    ElecEOMType c_elec_EOM_type;
//    QErrorCode ierr;
//    /* should be consistent with what defined in src/f03/openrsp_f.F90 */
//    switch (elec_EOM_type) {
//    case 0:
//        c_elec_EOM_type = ELEC_AO_D_MATRIX;
//        break;
//    case 1:
//        c_elec_EOM_type = ELEC_MO_C_MATRIX;
//        break;
//    case 2:
//        c_elec_EOM_type = ELEC_COUPLED_CLUSTER;
//        break;
//    default:
//        return QFAILURE;
//    }
//    c_open_rsp = (OpenRSP *)(*open_rsp);
//    ierr = OpenRSPSetElecEOM(c_open_rsp, c_elec_EOM_type);
//    return ierr;
//}

QErrorCode OpenRSPWriteFortranAdapter(const void *open_rsp,
                                      const char *file_name)
{
    FILE *fp_rsp;
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    fp_rsp = fopen(file_name, "a");
    if (fp_rsp==NULL) {
        printf("OpenRSPWriteFortranAdapter>> file name %s\n", file_name);
        QErrorExit(FILE_AND_LINE, "failed to open the file");
    }
    c_open_rsp = (OpenRSP *)open_rsp;
    ierr = OpenRSPWrite(c_open_rsp, fp_rsp);
    fprintf(fp_rsp, "\n");
    fclose(fp_rsp);
    return ierr;
}

QErrorCode OpenRSPWriteStdOutFortranAdapter(const void *open_rsp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)open_rsp;
    ierr = OpenRSPWrite(c_open_rsp, stdout);
    fprintf(stdout, "\n");
    return ierr;
}

QErrorCode OpenRSPDestroyFortranAdapter(void **open_rsp)
{
    OpenRSP *c_open_rsp;
    QErrorCode ierr;
    c_open_rsp = (OpenRSP *)(*open_rsp);
    ierr = OpenRSPDestroy(c_open_rsp);
    *open_rsp = NULL;
    open_rsp = NULL;
    return ierr;
}

