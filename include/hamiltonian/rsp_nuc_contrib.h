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

   This is the header file of (derivatives of) nuclear repulsion and
   nuclei-field interaction.

   2014-12-11, Bin Gao:
   * first version
*/

#if !defined(RSP_NUC_CONTRIB_H)
#define RSP_NUC_CONTRIB_H

/* QcMatrix library */
#include "qcmatrix.h"

/* context of (derivatives of) nuclear repulsion and nuclei-field interaction */
typedef struct {
    QInt num_atoms;          /* number of atoms */
    QReal *atom_coord;       /* coordinates of atoms */
    QReal *atom_charge;      /* charges of atoms */
    QReal dipole_origin[3];  /* dipole origin */
    QReal gauge_origin[3];   /* gauge origin */
} RSPNucContrib;

/* functions related to the (derivatives of) nuclear repulsion and nuclei-field interaction */
extern QErrorCode RSPNucContribCreate(RSPNucContrib*,
                                      const QInt,
                                      const QReal*,
                                      const QReal*);
extern QErrorCode RSPNucContribSetDipoleOrigin(RSPNucContrib*,const QReal[3]);
extern QErrorCode RSPNucContribSetGaugeOrigin(RSPNucContrib*,const QReal[3]);
extern QErrorCode RSPNucContribWrite(const RSPNucContrib*,FILE*);
extern QErrorCode RSPNucContribGet(const RSPNucContrib*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPNucContribGetNumAtoms(const RSPNucContrib*,QInt*);
extern QErrorCode RSPNucContribDestroy(RSPNucContrib*);

#endif
