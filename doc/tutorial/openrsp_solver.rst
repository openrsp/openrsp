===============================
Linear Response Equation Solver
===============================
.. include:: background.rst

OpenRSP API: OpenRSPSetLinearRSPSolver
======================================
.. include:: background.rst

Both in the response function and the residue calculations, OpenRSP needs
to solve the linear response equation, that will be accomplished by the
linear response equation solver of the host program.

The host program can specify the routine of this solver by calling:

.. c:function:: QErrorCode OpenRSPSetLinearRSPSolver(open_rsp, user_ctx, get_linear_rsp_solution)

and the last argument is the callback function of the solver.

.. nextslide::
   :increment:
.. include:: background.rst

The requirement for this callback function ``get_linear_rsp_solution`` has
been discussed in Chapter 4 "**OPENRSP CALLBACK FUNCTIONS**" of the OpenRSP
Manual.

As a user of OpenRSP, one should be aware that OpenRSP will send to the solver,
multiple right hand side (RHS) vectors of the time-dependent self-consistent-field
(TDSCF) equation. Because solving multiple response parameters (corresponding
to multiple RHS vectors) together may help the convergence of the solving
process.

