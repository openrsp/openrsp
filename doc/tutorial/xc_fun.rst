.. _chapter_xc_fun:

Exchange-Correlation Functionals
================================

Like the :ref:`chapter_one_elec_oper` and :ref:`chapter_two_elec_oper`, users
can call the API :c:func:`OpenRSPAddXCFun` **several times** to add several
XC functionals into a linked list in OpenRSP.

Each node in the linked list corresponds to an XC functional, which will be
visited by OpenRSP during calculations to make sure all XC functional
contributions are taken into account.

This API and its callback functions have been described in
:ref:`chapter_api_reference` and :ref:`chapter_callback_functions`. Users can
also find examples in the OpenRSP unit testing (files in the directory
``tests``).

Similar to one- and two-electron operators, arguments ``num_pert_lab``,
``pert_labels`` and ``pert_max_orders`` in the API :c:func:`OpenRSPAddXCFun`,
will be used by OpenRSP to check if a perturbation tuple already results in
zero XC functional integrals, and so that the callback functions will not be
invoked.

However, due to its **nonlinear** dependency of density matrix, the
calculations of integrals and expectation values of XC functionals are much
different from those of the one- and two-electron operators.

Therefore, additional information regarding the (perturbed) density matrices
will be sent to the callback functions :c:func:`get_xc_fun_mat` and
:c:func:`get_xc_fun_exp`, to facilitate such calculations.

Users are therefore recommend to read carefully the descriptions of these two
callback functions in :ref:`chapter_callback_functions`, to prepare their own
callback functions.
