.. _chapter_response_solver:

Linear Response Equation Solver
===============================

Both in the response function and the residue calculations, OpenRSP needs to
solve the linear response equation, that will be accomplished by the linear
response equation solver from the host program.

The host program can specify the solver by calling:

.. c:function:: QErrorCode OpenRSPSetLinearRSPSolver(open_rsp, user_ctx, get_linear_rsp_solution)

where the last argument is the callback function of the solver.

The requirement for this callback function :c:func:`get_linear_rsp_solution`
has been discussed in :ref:`chapter_callback_functions`.

As a user of OpenRSP, one should be aware that OpenRSP will send to the solver,
multiple right hand side (RHS) vectors of the time-dependent self-consistent-field
(TDSCF) equation. Because solving multiple response parameters (corresponding
to multiple RHS vectors) together may help the convergence of the solving
process.
