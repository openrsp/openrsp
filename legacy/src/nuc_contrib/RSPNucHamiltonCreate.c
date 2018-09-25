/* OpenRSP: open-ended library for response theory
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

   This file implements the function RSPNucHamiltonCreate().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief creates the context of nuclear Hamiltonian, should be called at first
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{inout} nuc_hamilton the context of nuclear Hamiltonian
    \param[QInt:int]{in} num_pert number of different perturbation labels that can
        act as perturbations on the nuclear Hamiltonian
    \param[QInt:int]{in} pert_labels all the different perturbation labels
    \param[QInt:int]{in} pert_max_orders maximum allowed order of each perturbation (label)
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetNucContrib:void]{in} get_nuc_contrib user specified function for
        getting nuclear contributions
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonCreate(RSPNucHamilton *nuc_hamilton,
                                const QInt num_pert,
                                const QInt *pert_labels,
                                const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid *user_ctx,
#endif
                                const GetNucContrib get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                const QInt num_atoms)
{
    QInt ipert,jpert;  /* incremental recorder over perturbations */
    if (num_pert>0) {
        nuc_hamilton->num_pert = num_pert;
    }
    else {
        printf("RSPNucHamiltonCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbations");
    }
    nuc_hamilton->pert_labels = (QInt *)malloc(num_pert*sizeof(QInt));
    if (nuc_hamilton->pert_labels==NULL) {
        printf("RSPNucHamiltonCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_labels");
    }
    nuc_hamilton->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (nuc_hamilton->pert_max_orders==NULL) {
        printf("RSPNucHamiltonCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_max_orders");
    }
    for (ipert=0; ipert<num_pert; ipert++) {
        /* each element of \var{pert_labels} should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("RSPNucHamiltonCreate>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       jpert,
                       pert_labels[jpert]);
                printf("RSPNucHamiltonCreate>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "same perturbation not allowed");
            }
        }
        nuc_hamilton->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("RSPNucHamiltonCreate>> order of %"QINT_FMT"-th perturbation (%"QINT_FMT") is %"QINT_FMT"\n",
                   ipert,
                   pert_labels[ipert],
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        nuc_hamilton->pert_max_orders[ipert] = pert_max_orders[ipert];
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    nuc_hamilton->user_ctx = user_ctx;
#endif
    nuc_hamilton->get_nuc_contrib = get_nuc_contrib;
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    nuc_hamilton->num_atoms = num_atoms;
    return QSUCCESS;
}
