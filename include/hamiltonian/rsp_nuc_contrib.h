/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

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

   This is the header file of nuclear contributions, including (derivatives of)
   nuclear repulsion and nuclei-field interaction.

   2014-12-11, Bin Gao:
   * first version
*/

#if !defined(RSP_NUC_CONTRIB_H)
#define RSP_NUC_CONTRIB_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function to get the nuclear contributions */
typedef QVoid (*GetNucContrib)(const QInt,
                               const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QReal*);

/* context of nuclear Hamiltonian */
typedef struct {
    QInt num_pert;                  /* number of perturbations that the nuclear Hamiltonian depend on */
    QInt *pert_labels;              /* labels of the perturbations */
    QInt *pert_max_orders;          /* maximum allowed orders of the perturbations */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                /* user-defined callback function context */
#endif
    GetNucContrib get_nuc_contrib;  /* user specified function for getting nuclear contributions */
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    QInt num_atoms;
} RSPNucHamilton;

/* functions related to the nuclear contributions */
extern QErrorCode RSPNucHamiltonCreate(RSPNucHamilton*,
                                       const QInt,
                                       const QInt*,
                                       const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                       QVoid*,
#endif
                                       const GetNucContrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                       const QInt);
extern QErrorCode RSPNucHamiltonAssemble(RSPNucHamilton*);
extern QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton*,FILE*);
extern QErrorCode RSPNucHamiltonGetContributions(const RSPNucHamilton*,
                                                 const QInt,
                                                 const QInt*,
                                                 const QInt,
                                                 QReal*);
extern QErrorCode RSPNucHamiltonDestroy(RSPNucHamilton*);
/*FIXME: RSPNucHamiltonGetNumAtoms() to be removed after perturbation free scheme implemented*/
extern QErrorCode RSPNucHamiltonGetNumAtoms(const RSPNucHamilton*,QInt*);

#endif
