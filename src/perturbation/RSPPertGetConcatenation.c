/* RSPPert: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

   RSPPert is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   RSPPert is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with RSPPert. If not, see <http://www.gnu.org/licenses/>.

   This file implements the function RSPPertGetConcatenation().

   2015-06-28, Bin Gao:
   * first version
*/

#include "eom/openrsp_perturbation.h"

/*% \brief gets the ranks of components of sub-perturbation tuples (with
        same perturbation label) for given components of the corresponding
        concatenated perturbation tuple
    \author Bin Gao
    \date 2015-06-28
    \param[RSPPert:struct]{in} rsp_pert context of all perturbations involved
        in calculations
    \param[QInt:int]{in} pert_label the perturbation label
    \param[QInt:int]{in} first_cat_comp rank of the first component of the
        concatenated perturbation tuple
    \param[QInt:int]{in} num_cat_comps number of components of the concatenated
        perturbation tuple
    \param[QInt:int]{in} num_sub_tuples number of sub-perturbation tuples to
        construct the concatenated perturbation tuple
    \param[QInt:int]{in} len_sub_tuples length of each sub-perturbation tuple,
        size is ``num_sub_tuples``; so that the length of the concatenated
        perturbation is ``sum(len_sub_tuples)``
    \var[QInt:int]{out} rank_sub_comps ranks of components of sub-perturbation
        tuples for the corresponding component of the concatenated perturbation
        tuple, i.e. ``num_cat_comps`` components starting from the one with rank
        ``first_cat_comp``, size is therefore ``num_sub_tuples*num_cat_comps``,
        and arranged as ``[num_cat_comps][num_sub_tuples]``
    \return[QErrorCode:int] error information
*/
QErrorCode RSPPertGetConcatenation(const RSPPert *rsp_pert,
                                   const QInt pert_label,
                                   const QInt first_cat_comp,
                                   const QInt num_cat_comps,
                                   const QInt num_sub_tuples,
                                   const QInt *len_sub_tuples,
                                   QInt *rank_sub_comps)
{
/*FIXME: zero-based or one-based numbering*/
    rsp_pert->get_pert_concatenation(pert_label,
                                     first_cat_comp,
                                     num_cat_comps,
                                     num_sub_tuples,
                                     len_sub_tuples,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     rsp_pert->user_ctx,
#endif
                                     rank_sub_comps);
    return QSUCCESS;
}
