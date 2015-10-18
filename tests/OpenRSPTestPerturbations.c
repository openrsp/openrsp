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

#include "OpenRSPTestPerturbations.h"

const QcPertInt PERT_GEOMETRIC = 1;
const QcPertInt PERT_DIPOLE = 2;
const QcPertInt PERT_MAGNETIC = 5;
const QInt MAX_ORDER_GEOMETRIC = 7;
const QInt MAX_ORDER_DIPOLE = 1;
const QInt MAX_ORDER_MAGNETIC = 7;
/* labels of all perturbations */
const QcPertInt ALL_PERT_LABELS[]={PERT_GEOMETRIC,
                                   PERT_DIPOLE,
                                   PERT_MAGNETIC};
/* allowed maximal orders of all perturbations */
const QInt ALL_PERT_MAX_ORDERS[]={MAX_ORDER_GEOMETRIC,
                                  MAX_ORDER_DIPOLE,
                                  MAX_ORDER_MAGNETIC};
/* sizes of all perturbations up to their maximal orders */
const QInt ALL_PERT_SIZES[]={
    12,78,364,1365,4368,12376,31824,  /* geometric derivatives (4 atoms) */
    3,                                /* electric dipole */
    3,6,10,15,21,28,36};              /* magnetic derivatives */
const QInt NUM_ATOMS=4;
