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

   This is the header file of nuclear contributions, including (derivatives of)
   nuclear repulsion and nuclei-field interaction.

   2014-12-11, Bin Gao:
   * first version
*/

#if !defined(RSP_NUC_HAMILTONIAN_H)
#define RSP_NUC_HAMILTONIAN_H

/* QcMatrix library */
#include "qcmatrix.h"

/* context of nuclear Breit-Pauli Hamiltonian */
typedef struct {
    QInt num_atoms;        /* number of atoms */
    QReal *atom_coord;     /* coordinates of atoms */
    QReal *atom_charge;    /* charges of atoms */
    QReal *dipole_origin;  /* dipole origin */
    QReal *gauge_origin;   /* gauge origin */
} RSPNucHamilton;

/* functions related to the nuclear Breit-Pauli Hamiltonian */
extern QErrorCode RSPNucHamiltonCreate(RSPNucHamilton*);
extern QErrorCode RSPNucHamiltonSetGeoPerturbations(RSPNucHamilton*,
                                                    const QInt,
                                                    const QReal*,
                                                    const QReal*);
extern QErrorCode RSPNucHamiltonSetScalarPotential(RSPNucHamilton*,
                                                   const QReal[3]);
extern QErrorCode RSPNucHamiltonSetVectorPotential(RSPNucHamilton*,
                                                   const QReal[3]);
extern QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton*,FILE*);
extern QErrorCode RSPNucHamiltonGetContributions(const RSPNucHamilton*,
                                                 const QInt,
                                                 const QInt*,
                                                 const QInt*,
                                                 const QInt,
                                                 QReal*);
extern QErrorCode RSPNucHamiltonDestroy(RSPNucHamilton*);
extern QErrorCode RSPNucHamiltonGetNumAtoms(const RSPNucHamilton*,QInt*);

#endif
