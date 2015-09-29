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

*/

#include "RSPPerturbation.h"

/*% \brief gets the information of perturbations (orders and number of components)
        from a given perturbation tuple
    \author Bin Gao
    \date 2015-06-29
    \param[RSPPert:struct]{in} rsp_pert context of all perturbations involved
        in calculations
    \param[QInt:int]{in} len_tuple length of the perturbation tuple, in which
        identical perturbation labels should be consecutive
    \param[QInt:int]{in} pert_tuple the perturbation tuple
    \param[QInt:int]{in} num_freq_configs number of different frequency configurations
    \param[QReal:real]{in} pert_freqs complex frequencies of each perturbation label
        over all frequency configurations
    \var[QInt:int]{out} num_pert number of different perturbations from the given
        perturbation tuple
    \var[QInt:int]{out} pert_orders orders of each perturbation
    \var[QInt:int]{out} pert_num_comps numbers of components of each perturbation,
        from the first order up to the order in \var{pert_orders}
    \return[QErrorCode:int] error information
*/
//QErrorCode RSPPertGetFromTuple(const RSPPert *rsp_pert,
//                               const QInt len_tuple,
//                               const QInt *pert_tuple,
//                               const QInt num_freq_configs,
//                               const QReal *pert_freqs,
//                               QInt *num_pert,
//                               QInt *pert_orders,
//                               QInt *pert_num_comps)
//{
//    QInt ipert,jpert;  /* incremental recorders */
//    QInt first_id;     /* first identical pertubation label in the tuple */
//    QInt last_id;      /* last identical pertubation label in the tuple */
//    QBool non_id;      /* indicates if non-identical label found */
//    /* we first get the consecutive identical pertubation labels */
//    first_id = 0;
//    non_id = QFALSE;
//    for (ipert=first_id; ipert<len_tuple-1; ipert++) {
//        if (pert_tupe[ipert]!=pert_tupe[ipert+1]) {
//            last_id = ipert;
//            non_id = QTRUE;
//            break;
//        }
//    }
//    if (non_id=QTRUE) {
//    }
//    else {
//    }
//    /* loops over all known perturbation labels and checks if they are in the tuple */
//    for (ipert=0; ipert<rsp_pert->num_pert; ipert++) {
//        /* loops over perturbation labels in the tuple */
//        for (ipert=first_id; ipert<len_tuple; ipert++) {
//            /* checks if the given label is known */
//            if (pert_tupe[ipert]==rsp_pert->pert_labels[jpert]) {
//                
//            }
//        }
//    }
//    return QSUCCESS;
//}

