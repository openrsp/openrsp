/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud

  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


  This is the header file of perturbations.

  2014-07-31, Bin Gao:
  * first version
*/

#include "OpenRSP.h"

#define NUM_ALL_PERT 3
#define PERT_GEOMETRIC 1
#define PERT_DIPOLE 2
#define PERT_MAGNETIC 5
#define MAX_ORDER_GEOMETRIC 7
#define MAX_ORDER_DIPOLE 1
#define MAX_ORDER_MAGNETIC 7
#define NUM_ATOMS 4

//extern const QcPertInt ALL_PERT_LABELS[NUM_ALL_PERT];
//extern const QInt ALL_PERT_MAX_ORDERS[NUM_ALL_PERT];
//extern const QInt ALL_PERT_SIZES[15];
//extern const QInt NUM_ATOMS;
