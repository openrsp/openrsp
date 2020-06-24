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

   This file is the header file of Lagrangian in atomic-orbital density matrix
   based response theory.

   2020-05-02, Bin Gao:
   * first version
*/

#pragma once

#include "tSymbolic/QuantChem.hpp"

#include "OpenRSP/Lagrangian.hpp"

namespace OpenRSP
{
    /* @class:LagrangianDAO[Lagrangian in atomic-orbital density matrix based response theory.]

       In the framework of atomic-orbital (AO) density matrix based response
       theory, we need to store:

       . overlap operator,
       . zero-, one- and two-electron operators, and exchange-correlation (XC)
         functional derivative in the Hamiltonian,
       . reference state,
       . excited states, and
       . a linear response solver.

       Because response theory These quantities

       Among them, different operators in the Hamiltonian will be stored in a
       corresponding `std::vector` as one may need several such operators in
       the Hamiltonian.
     */
    class LagrangianDAO: virtual public Lagrangian
    {
        public:
            explicit LagrangianDAO() noexcept {}
            virtual ~LagrangianDAO() noexcept = default;
            inline void set_integator(const std::shared_ptr<tSymbolic::OneElecIntegrator> integator) noexcept
            {
                m_one_integrator = integator;
            }
        private:
            std::vector<std::shared_ptr<tSymbolic::ZeroElecOper>> m_zero_opers;
            std::vector<std::shared_ptr<tSymbolic::OneElecOper>> m_one_opers;
            std::vector<std::shared_ptr<tSymbolic::TwoElecOper>> m_two_opers;
            std::vector<std::shared_ptr<tSymbolic::ExchCorrFunctional>> m_xc_func;
            std::shared_ptr<tSymbolic::OneElecIntegrator> m_one_integrator;
    };
}
