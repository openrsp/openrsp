.. _chapter_one_elec_oper:

One-Electron Operators
======================

The strategy of treating zero-, one- and two-electron operators, as well as
exchange-correlation functionals is different from that of overlap operator and
linear response equation solver. For the latter, the OpenRSP API will be
usually called once to set up the approprite callback functions.

Taking the one-electron operators as an example, host programs may have
different one-electron functions for different operators. If they do not want
to provide OpenRSP a general callback function, instead they can call the
following API **several times**:

.. c:function:: QErrorCode OpenRSPAddOneOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_one_oper_mat, get_one_oper_exp)

to **add several one-electron operators** to the electronic Hamiltonian.

Inside OpenRSP, these one-electron operators will be saved in a linked list, in
which each node corresponds to a one-electron operator.

The arguments of this API are similar to those of the overlap operator
:c:func:`OpenRSPSetOverlap`, and has been described in
:ref:`chapter_api_reference`.

The callback functions :c:func:`get_one_oper_mat` and
:c:func:`get_one_oper_exp` are presented in :ref:`chapter_callback_functions`.
Users can also find examples in the OpenRSP unit testing (files in the
directory ``tests``).

The arguments ``num_pert_lab``, ``pert_labels`` and ``pert_max_orders`` will be
also used in a similar way as those of the overlap operator, that OpenRSP will
not invoke the callback functions if a perturbation tuple already result in
zero one-electron integrals.
