%  OpenRSP: open-ended library for response theory
%  Copyright 2015 Radovan Bast,
%                 Daniel H. Friese,
%                 Bin Gao,
%                 Dan J. Jonsson,
%                 Magnus Ringholm,
%                 Kenneth Ruud,
%                 Andreas Thorvaldsen
%  
%  OpenRSP is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  OpenRSP is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU Lesser General Public License for more details.
%  
%  You should have received a copy of the GNU Lesser General Public License
%  along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.
%
%  This file implements the varaibles in symbolic computations.
%
%  2015-06-28, Bin Gao:
%  * first version

@* Variables in Symbolic Computations.

To calculate (derivatives of) nuclear repulsion and nuclei--field interaction,
we choose the molecular Breit-Pauli Hamiltonian\cite{Jaszunski-HCC-361}
\begin{equation}
  \label{eq-BP-Hamiltonian}
  \hat{H}_{\text{nucl}}^{\text{BP}}
  =\sum_{K}\left[Z_{K}\phi_{K} %
    -\boldsymbol{m}_{K}\cdot\boldsymbol{B} %
    +\frac{1}{2}\sum_{L\ne K} %
    \frac{Z_{K}Z_{L}}{\vert\boldsymbol{R}_{K}-\boldsymbol{R}_{L}\vert}\right],
\end{equation}
where we have neglected the nuclear kinetic energy
$\sum_{K}\frac{\pi_{K}^{2}}{2M_{K}}$. The three terms in
Eq.~(\ref{eq-BP-Hamiltonian}) are respectively the interaction of the nuclear
charge with an external electrostatic potential $\phi_{K}$, the nuclear
spin--Zeeman interaction of the nuclear magnetic moment $\boldsymbol{m}_{K}$
with an external magnetic field induction $\boldsymbol{B}$, and the nuclear
repulsion.

The contributions of nuclear Hamiltonian to the molecular properties are
normally evaluated at zero fields, and we will in the following discuss
geometric, electric and magnetic perturbations. When there is no external
perturbations, the nuclear Hamiltonian is simply the nuclear repulsion

@(RSPSymbolicVar.c@>=

#include <stdlib.h>
#include <string.h>
#include <limits.h>

for (nuc_K=0; nuc_K++; nuc_K<num_nuclei) {
@<NucDistance@>;
}

@ The nuclear distance is calculated as

@<NucDistance@>=
nuc_distance = 0;

