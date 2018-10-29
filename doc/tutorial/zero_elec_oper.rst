.. _chapter_zero_elec_oper:

Zero-Electron Operators
=======================

The computation of zero-electron operator contributions (like nuclear
Hamiltonian) to the molecular properties in OpenRSP do not involve the
electronic quantities (i.e. the (perturbed) density matrices), and will be
calculated after all electronic contributions obtained.

Similar to one- and two-electron operators, users can call the API
:c:func:`OpenRSPAddZeroOper` **serveral times** to add a few zero-electron
operators to OpenRSP (which are arranged in a linked list).

Please see :ref:`chapter_api_reference` and :ref:`chapter_callback_functions`
for the description of this API and the callback function
:c:func:`get_zero_oper_contrib`.
