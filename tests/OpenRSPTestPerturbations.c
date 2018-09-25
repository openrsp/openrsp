/*
  OpenRSP: open-ended library for response theory
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

*/

#include "OpenRSPTestPerturbations.h"

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
