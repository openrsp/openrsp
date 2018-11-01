
OpenRSP: open-ended response theory
===================================

OpenRSP is program library that uses recursive routines to identify and
assemble contributions to molecular properties ("response properties") as they
are formulated in the theory called "response theory" which is used in theoretical
chemistry.

By its recursive structure, OpenRSP makes it possible to calculate
response properties (hereafter just "properties") of arbitrary complexity in
an analytical manner, not resorting to numerical schemes like
finite difference methods for the calculation of these properties. 


What is the advantage?
----------------------

Compared to analytical methods, numerical approaches
may be associated with a greater degree of uncertainty related to
accuracy and practical feasibility of the calculation. 

Today's programs written 
for the calculation of response properties either do not have a recursive
structure, or they use numerical methods to different extents, or both. In the
cases where existing programs use an analytical approach, they are either not
recursive (which typically means that a new program routine must be written for each
new property for which calculation is desired), or they can only be used for a
limited category of properties. 

The structure of OpenRSP solves the task of
identifying and assembling contributions to molecular properties "once and for
all".


A brief history
---------------

The OpenRSP project has been under development since the mid-2000's, but did
not feature a recursive structure from the beginning. The development of the
approach used in today's recursive structure was started in the autumn of 2011,
and the first results (i.e. the first molecular properties calculated and
regarded as suitable for publication) were obtained about one year later (close
to the end of 2012/beginning of 2013). 

Since then, we have worked on improving the functionality and developing a version
of the program ready for publication. This new work has for
example entailed making it possible to calculate other categories of properties
(so-called residues of response properties usable for the calculation of multiphoton
absorption strengths) and making OpenRSP modular - i.e. cleanly separated from and 
more easily able communicate with other programs. We are working towards releasing 
an official version of OpenRSP in 2018.

.. toctree::
   :maxdepth: 1
   :caption: The people behind OpenRSP

   authors.rst
   citations.rst
   change_log.rst

.. toctree::
   :maxdepth: 1
   :caption: Manual

   manual/notations_and_conventions.rst
   manual/installation.rst
   manual/unit_testing.rst
   manual/sphinx.rst
   manual/api_reference.rst
   manual/callback_functions.rst
   manual/limitations.rst
   manual/future_plans.rst

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial/getting_started.rst
   tutorial/openrsp_context.rst
   tutorial/perturbations.rst
   tutorial/overlap_operator.rst
   tutorial/one_elec_oper.rst
   tutorial/two_elec_oper.rst
   tutorial/xc_fun.rst
   tutorial/zero_elec_oper.rst
   tutorial/response_solver.rst
   tutorial/response_functions.rst
   tutorial/residues.rst

.. toctree::
   :maxdepth: 1
   :caption: Documentation for developers

   developer/background_and_rationale.rst
   developer/openrsp_design.rst
   developer/rules.rst
   developer/coding_standards.rst
   developer/files_and_directories.rst
.. developer/pseudocode.rst
