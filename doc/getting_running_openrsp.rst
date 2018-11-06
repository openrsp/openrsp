.. _chapter_getting_running_openrsp:

Get and run OpenRSP
===================

If you want to get OpenRSP and use it for calculations, then please make note of the following:
OpenRSP is a program library that manages the calculation of response properties, and it 
**cannot** calculate these properties without getting contributions like perturbed one- and
two-electron integrals or solutions of response equations from other codes to which it connects
through the application programming interface (API). This means that if you download and build 
OpenRSP from its `GitHub repository <https://github.com/openrsp/openrsp>`_, the compiled product 
will not on its own be able to calculate response properties. A set of API connections to enable 
OpenRSP to manage response property calculations can for example be made in quantum chemistry 
programs the where necessary routines for these contributions are implemented.

Consequently, in order to use OpenRSP for calculations, it is necessary to use it in a 
*host program* into which OpenRSP has been incorporated in this way, and a list of such programs
is kept at the :ref:`chapter_programs_with_openrsp` page. The specific way in which OpenRSP is 
invoked in a host program - i.e. the way that you can make OpenRSP calculate something in that 
program - is a feature of each such program, and you must therefore follow the relevant 
instructions to achieve this, as may for example be shown in the user manual for the host program
that you want to use.
