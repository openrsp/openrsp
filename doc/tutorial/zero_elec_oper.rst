.. _chapter_zero_elec_oper:

Zero-Electron Operators
=======================

The computation of zero-electron operator contributions (like nuclear
Hamiltonian) to the molecular properties in OpenRSP do not involve the
electronic quantities (i.e. the (perturbed) density matrices), and will be
calculated after all electronic contributions obtained.

Similar to one- and two-electron operators, users can call the following API
**serveral times** to add a few zero-electron operators to OpenRSP (which are
arranged in a linked list):

.. c:function:: QErrorCode OpenRSPAddZeroOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_zero_oper_contrib)

Please see :ref:`chapter_api_reference` and :ref:`chapter_callback_functions`
for the description of this API and the callback function
:c:func:`get_zero_oper_contrib`.
