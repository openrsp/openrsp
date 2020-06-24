/* OpenRSP: open-ended library for response theory
   Copyright 2015-2020 Radovan Bast,
                       Daniel H. Friese,
                       Bin Gao,
                       Dan J. Jonsson,
                       Magnus Ringholm,
                       Simen S. Reine
                       Kenneth Ruud,
                       Andreas Thorvaldsen

   This source code form is subject to the terms of the
   GNU Lesser General Public License, version 2.1.
   If a copy of the GNU LGPL v2.1 was not distributed with this
   code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file is the header file of Lagrangian in response theory.

   2020-05-02, Bin Gao:
   * first version
*/

#pragma once

namespace OpenRSP
{
    /* @class:Lagrangian[Abstract Lagrangian class.]
     */
    class Lagrangian
    {
        public:
            explicit Lagrangian() noexcept = default;
            virtual ~Lagrangian() noexcept = default;
    };
}
