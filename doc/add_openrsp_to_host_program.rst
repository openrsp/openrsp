.. _chapter_add_openrsp_to_host_program:

Add OpenRSP to a quantum chemistry program
==========================================

If you want to add OpenRSP to a quantum chemistry program, then you are free to do so provided
that you do not violate OpenRSP's LGPL v2.1 software license as described on OpenRSP's 
`GitHub repository <https://github.com/openrsp/openrsp>`_.

In order to enable OpenRSP to work as intended, you must provide routines that connect to the
OpenRSP application programming interface (API) to give OpenRSP access to contributions such as
perturbed one- and two electron integrals, exchange-correlation contributions if calculations at
the density-functional theory (DFT) level is desired, or solution routines for response equations.
Please note that OpenRSP is a program library that manages the calculation of response properties,
and it **cannot** carry out actual such calculations without getting contributions like the ones
mentioned here from program routines that are external to OpenRSP.

A description of the API is under development and will be made available on this website.
In the meantime, questions may be directed to the :ref:`chapter_authors`.
