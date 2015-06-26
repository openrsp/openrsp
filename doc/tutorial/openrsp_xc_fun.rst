==================================================
Exchange-Correlation Functionals (not implemented)
==================================================
.. include:: background.rst

OpenRSP API: OpenRSPAddXCFun
============================
.. include:: background.rst

Like the :ref:`slide-one-electron` and :ref:`slide-two-electron`, users
can call the following API **several times**:

.. c:function:: QErrorCode OpenRSPAddXCFun(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_xc_fun_mat, get_xc_fun_exp)

to add several exchange-correlation (XC) functionals into a linked list
in OpenRSP.

Each node in the linked list corresponds to an XC functional, which will
be visited by OpenRSP during calculations to make sure all XC functional
contributions are taken into account.

.. nextslide::
   :increment:
.. include:: background.rst

This API and its callback functions have been described in Chapter 3
"**OPENRSP API REFERENCE**" and Chapter 4 "**OPENRSP CALLBACK FUNCTIONS**"
respectively in the OpenRSP Manual. Users can also find examples in the
OpenRSP unit testing (files in ``tests``).

The arguments ``num_pert``, ``pert_labels`` and ``pert_max_orders``, as
usual, will be checked by OpenRSP if a perturbation tuple already results
in zero XC functional integrals, and so that the callback functions will
not be invoked.

.. nextslide::
   :increment:
.. include:: background.rst

However, due to its **nonlinear** dependency of density matrix, the calculations
of integrals and expectation values of XC functionals are much different from
those of the one- and two-electron operators.

Therefore, additional information regarding the (perturbed) density matrices
will be sent to the callback functions ``get_xc_fun_mat`` and ``get_xc_fun_exp``,
to facilitate such calculations.

Users are therefore recommend to read carefully the descriptions of these two
callback functions in Chapter 4 "**OPENRSP CALLBACK FUNCTIONS**" of the OpenRSP
Manual, in order to prepare their own callback functions.

