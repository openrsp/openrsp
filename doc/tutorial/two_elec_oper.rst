.. _chapter_two_elec_oper:

Two-Electron Operators
======================

Similar to :ref:`chapter_one_elec_oper`, the two-electron operators will be
saved in a linked list in OpenRSP. Users can call the following API:

.. c:function:: QErrorCode OpenRSPAddTwoOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_two_oper_mat, get_two_oper_exp)

**several times** to add several two-electron operators to the electronic
Hamiltonian.

Each node in the linked list corresponds to a two-electron operator. During
calculations, OpenRSP will walk through the linked list to get correct
contributions from the two-electron part.

This API and its callback functions have been described in
:ref:`chapter_api_reference` and :ref:`chapter_callback_functions`.  Users can
also find examples in the OpenRSP unit testing (files in the directory
``tests``).

The arguments ``num_pert_lab``, ``pert_labels`` and ``pert_max_orders`` are
used in a similar way as those of overlap, one-electron integrals, that OpenRSP
will not invoke the callback functions if a perturbation tuple already results
in zero two-electron integrals.

In this tutorial, we will further discuss the callback function
:c:func:`get_two_oper_exp`:

.. c:function:: QVoid get_two_oper_exp(oper_num_pert, oper_pert_labels, oper_pert_orders, dmat_len_tuple, num_LHS_dmat, LHS_dens_mat, num_RHS_dmat, RHS_dens_mat, user_ctx, num_exp, val_exp)

to calculate the expectation values
:math:`\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\boldsymbol{D}_{\text{LHS}}^{a''b''c''\cdots})
\boldsymbol{D}_{\text{RHS}}^{a'b'c'\cdots}]`.

The calculations of two-electron integrals are most time-consuming, and
different host programs may have different strategies to calculate the
two-electron integrals.

OpenRSP therefore sends both the left hand side and right hand side arrays of
density matrices to the callback function, and the host program can decide how
the expectation values will be calculated.

For instance, if OpenRSP sends:

#. ``dmat_len_tuple=2``,
#. ``num_LHS_dmat[2]={2, 2}``,
#. ``LHS_dens_mat[3]={D1,D2, D3,D4}``,
#. ``num_RHS_dmat[2]={2, 3}``,
#. ``RHS_dens_mat[5]={D5,D6, D7,D8,D9}``,

which means OpenRSP wants the following expectation values back:

| ``/* {D1,D2}`` :math:`\otimes` ``{D5,D6}*/``
| ``[`` :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D1})\texttt{D5}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D1})\texttt{D5}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D1})\texttt{D6}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D1})\texttt{D6}])`,
| :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D2})\texttt{D5}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D2})\texttt{D5}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D2})\texttt{D6}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D2})\texttt{D6}])`,
| ``/* {D3,D4}`` :math:`\otimes` ``{D7,D8,D9}*/``
| :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D7}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D7}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D8}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D8}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D9}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D3})\texttt{D9}])`,
| :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D7}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D7}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D8}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D8}])`,
  :math:`\mathrm{Re}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D9}])`,
  :math:`\mathrm{Im}(\mathrm{Tr}[\boldsymbol{G}^{abc\cdots}(\texttt{D4})\texttt{D9}])` ``]``.
