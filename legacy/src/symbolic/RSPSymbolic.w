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
%  Main file of the symbolic computations in OpenRSP.
%
%  2015-06-28, Bin Gao:
%  * first version

\documentclass[baseclass=article]{cweb}

\usepackage{amsmath}

\begin{document}

\title{Symbolic Computations in OpenRSP}
\author{Bin Gao}

\maketitle

\tableofcontents

@* Introduction.

The recursive programming techniques described in
Ref.~\cite{Ringholm-JCC-35-622} make it possible to calculate molecular
properties of arbitrary complexity in an analytical manner. But if we are going
to implement MO and CC based response theories, it may mess up or not be easy
to (re)use the already developed codes.

Therefore, I (Bin Gao) will examine an alternative route to the MO and CC based
response theories, by developing a set of symbolic computation functions that
will be used by response theory codes to get molecular properties still in a
recursive and analytical manner.

These symbolic computation functions will be implemented in
\texttt{src/symbolic} using the literate programming approach (using the
\texttt{CWEB} system). In this way, I can describe the idea of implementation
and the real codes at the same time. After becoming familar with the literate
programming approach, it is actually, in my opinion, an efficient way for
methodology development.

@i RSPSymbolicVar.w

@* Bibliography.
\bibliographystyle{plain}
\renewcommand{\refname}{\vspace*{-4ex}}
\bibliography{references}

@
\end{document}
