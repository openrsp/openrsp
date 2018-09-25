/* OpenRSP: open-ended library for response theory
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

   This is the header file of electronic wave function.

   2015-06-28, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_WAV_FUN_H)
#define OPENRSP_WAV_FUN_H

/* QcMatrix library */
#include "qcmatrix.h"

/* type of electronic wave function */
typedef enum {ELEC_AO_D_MATRIX=0,
              ELEC_MO_C_MATRIX=1,
              ELEC_COUPLED_CLUSTER=2} ElecWavType;

/* implementation-specific functions related to the EOM of electrons */
typedef struct {
    QErrorCode (*eleceomxxx)();
} ElecEOMFun;

/* implementation-specific data of the EOM of electrons */
typedef struct {
    ElecEOMFun *elec_eom_fun;
} ElecEOM;

#endif
