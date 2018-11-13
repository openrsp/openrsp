

A brief history
===============


The OpenRSP core functionality
------------------------------

Work on the OpenRSP project began in the mid-2000's when the first work started on what
has become the present version of the program. During this time, the theoretical 
foundation of response theory on which OpenRSP was based was developed and used to create
program routines that were connected to the Dalton quantum chemistry program. This version
of the code was used to compute several response properties for which analytic calculation had not
been carried out before.

In 2008, we have generalized OpenRSP for the DIRAC program package which enabled us to access
a wealth of response properties at the 4-component relativistic level.

In 2011, work was started on a version - then still a part of Dalton - where recursion was used to
achieve an open-ended implementation of the theory, so that one set of routines
could be used to manage the calculation of any response property. This version forms the
basis of the present-day core functionality of the OpenRSP, but was since developed
further to include features such as calculation of single residues of response properties (of use
in the calculation of multiphoton strengths), calculation of multiple properties in one invocation with
reuse of common intermediate results, and restructuring of calls to external routines to reduce
recalculation of various contributions such as perturbed one- and two-electron integrals.


OpenRSP as a modular library with an API
----------------------------------------

In order to make OpenRSP into a modular library that was not tied to any one particular
quantum chemistry program - or *host program* - work began in 2013 on developing an
application programming interface (API) for OpenRSP, involving the creation of clearly defined
interfaces between OpenRSP and other codes, the use of callback routines
in order to abstract the way in which the OpenRSP core asks for contributions from external libraries,
and the development of the QcMatrix library to abstract and mediate matrix operations so 
that OpenRSP is agnostic to the underlying implementation of such operations. The first host program
to make use of this modular functionality is the `LSDalton <http://daltonprogram.org/>`_ quantum chemistry program.


Libraries for external contributions
------------------------------------

During the course of its execution, OpenRSP identifies various contributions that it must get
from libraries external to it in order to be able to assemble the response property or
properties to be calculated, such as perturbed one- and two-electron integral contributions,
exchange-correlation contributions if a density-functional theory calculation is requested,
or solution of so-called response equations. Therefore, the development of libraries that 
can provide such functionality at a sufficient level of generality - although not
necessarily driven by the demands of OpenRSP - has nevertheless been an important
concurrent task, and has resulted in the creation of sophisticated software without which
OpenRSP would not be able to do what it does best. Some of the libraries that are
presently used or have been used by OpenRSP are listed below:

* Gen1Int for the calculation of perturbed one-electron integrals
* cgto-diff-eri for the calculation of perturbed two-electron integrals
* HODI for the calculation of perturbed integrals
* XCint and XCfun for the calculation of exchange-correlation contributions
* A linear response equation solver by Sonia Coriani et al.
* FraME for a polarizable embedding description of molecular surroundings
* PCMSolver for a polarizable continuum description of molecular surroundings
