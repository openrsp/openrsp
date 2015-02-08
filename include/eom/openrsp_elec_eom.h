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

   This is defines the equation of motion (EOM) of electrons.

   2014-08-04, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_ELEC_EOM_H)
#define OPENRSP_ELEC_EOM_H

/* QcMatrix library */
#include "qcmatrix.h"

/* type of EOM of electrons */
typedef enum {ELEC_AO_D_MATRIX=0,
              ELEC_MO_C_MATRIX=1,
              ELEC_COUPLED_CLUSTER=2} ElecEOMType;

/* implementation-specific functions related to the EOM of electrons */
typedef struct {
    QErrorCode (*eleceomxxx)();
} ElecEOMFun;

/* implementation-specific data of the EOM of electrons */
typedef struct {
    ElecEOMFun *elec_eom_fun;
} ElecEOM;

#endif
