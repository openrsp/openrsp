===================
Nuclear Hamiltonian
===================
.. include:: background.rst
 
OpenRSP API: OpenRSPSetNucContributions
=======================================
.. include:: background.rst

The computation of nuclear contributions to the molecular properties 
in OpenRSP are not involve the electronic quantities (i.e. the (perturbed)
density matrices), and will be calculated after all electronic contributions
obtained.

Similar to :ref:`slide-overlap`, users only need to call the following
API once to set up the nuclear Hamiltonian for its contributions:

.. c:function:: QErrorCode OpenRSPSetNucContributions(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_nuc_contrib, num_atoms)

.. nextslide::
   :increment:
.. include:: background.rst

Please see Chapter 3 "**OPENRSP API REFERENCE**" and Chapter 4
"**OPENRSP CALLBACK FUNCTIONS**" in the OpenRSP Manual for the description
of this API and the callback function ``get_nuc_contrib``.

It should be noted that the last argument ``num_atoms`` in this API is the
number of atoms, which **will be removed** after the perturbation free scheme
implemented in OpenRSP, i.e., after the API :c:func:`OpenRSPSetPerturbations`
and its related core parts of OpenRSP are implemented.

