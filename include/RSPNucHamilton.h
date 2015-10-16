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


  <header name='RSPNucHamilton.h' author='Bin Gao' date='2014-12-11'>
    The header file of nuclear Hamiltonian used inside OpenRSP
  </header>
*/

#if !defined(RSP_NUCHAMILTON_H)
#define RSP_NUCHAMILTON_H

#include "qcmatrix.h"
#include "RSPPerturbation.h"

typedef void (*GetNucContrib)(const QInt,
                              const QcPertInt*,
                              const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                              void*,
#endif
                              const QInt,
                              QReal*);

typedef struct {
    QInt num_pert_lab;              /* number of different perturbation labels
                                       that can act as perturbations on the
                                       nuclear Hamiltonian */
    QInt nuc_num_pert;              /* number of perturbations on the
                                       nuclear Hamiltonian, only used for
                                       callback functions */
    QInt *pert_max_orders;          /* allowed maximal order of a perturbation
                                       described by exactly one of these
                                       different labels */
    QInt *nuc_pert_orders;          /* orders of perturbations on the
                                       nuclear Hamiltonian, only used for
                                       callback functions */
    QcPertInt *pert_labels;         /* all the different perturbation labels */
    QcPertInt *nuc_pert_labels;     /* labels of perturbations on the
                                       nuclear Hamiltonian, only used for
                                       callback functions */
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                 /* user-defined callback-function context */
#endif
    GetNucContrib get_nuc_contrib;  /* user-specified function for calculating
                                       contribution from the nuclear Hamiltonian */
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    QInt num_atoms;
} RSPNucHamilton;

extern QErrorCode RSPNucHamiltonCreate(RSPNucHamilton*,
                                       const QInt,
                                       const QcPertInt*,
                                       const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                       void*,
#endif
                                       const GetNucContrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                       const QInt);
extern QErrorCode RSPNucHamiltonAssemble(RSPNucHamilton*,const RSPPert*);
extern QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton*,FILE*);
extern QErrorCode RSPNucHamiltonGetContributions(RSPNucHamilton*,
                                                 const QInt,
                                                 const QcPertInt*,
                                                 const QInt,
                                                 QReal*);
extern QErrorCode RSPNucHamiltonDestroy(RSPNucHamilton*);
/*FIXME: RSPNucHamiltonGetNumAtoms() to be removed after perturbation free scheme implemented*/
extern QErrorCode RSPNucHamiltonGetNumAtoms(const RSPNucHamilton*,QInt*);

#endif
