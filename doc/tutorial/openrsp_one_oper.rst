.. _slide-one-electron:

======================
One-Electron Operators
======================
.. include:: background.rst

OpenRSP API: OpenRSPAddOneOper
==============================
.. include:: background.rst

The strategy of treating one- and two-electron operators, as well as
exchange-correlation functionals is different from that of overlap
integrals, nuclear Hamiltonian and linear response equation solver.
For the latter, the OpenRSP API will be usually called once to set
up the approprite callback functions.

Taking the one-electron operators as an example, host programs may have
different one-electron functions for different operators. If they do not
want to provide OpenRSP a general callback function, instead they can
call the following API **several times**:

.. c:function:: QErrorCode OpenRSPAddOneOper(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_one_oper_mat, get_one_oper_exp)

to **add several one-electron operators** to the electronic Hamiltonian.

Inside OpenRSP, these one-electron operators will be saved in a linked list,
in which each node corresponds to a one-electron operator.

The arguments of this API are similar to those of the overlap integrals
:c:func:`OpenRSPSetPDBS`, and has been described in Chapter 3
"**OPENRSP API REFERENCE**" of the OpenRSP Manual.

The callback functions ``get_one_oper_mat`` and ``get_one_oper_exp`` are
presented in Chapter 4 "**OPENRSP CALLBACK FUNCTIONS**" of the OpenRSP
Manual. Users can also find examples in the OpenRSP unit testing (files
in ``tests``).

The arguments ``num_pert``, ``pert_labels`` and ``pert_max_orders`` will
be also used in a similar way as those of overlap integrals, that OpenRSP
will not invoke the callback functions if a perturbation tuple already
result in zero one-electron integrals.

