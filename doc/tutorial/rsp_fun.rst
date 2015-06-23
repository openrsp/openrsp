==================
Response Functions
==================
.. include:: background.rst

OpenRSP API: OpenRSPGetRSPFun
=============================
.. include:: background.rst

The response functions can be calculated by calling the
OpenRSP API:

.. c:function:: QErrorCode OpenRSPGetRSPFun(open_rsp, ref_ham, ref_state, ref_overlap, num_props, len_tuple, pert_tuple, num_freq_configs, pert_freqs, kn_rules, size_rsp_funs, rsp_funs)

Suppose we want to calculate the polarizability :math:`\alpha` and
the first hyperpolarizability :math:`\beta` at different frequencies:

:math:`\alpha(-\omega_{1},\omega_{1})`, :math:`\alpha(-\omega_{2},\omega_{2})`,
:math:`\beta(-\omega_{3}-\omega_{4},\omega_{3},\omega_{4})`,
:math:`\beta(-\omega_{5}-\omega_{6},\omega_{5},\omega_{6})`,
:math:`\beta(-\omega_{7}-\omega_{8},\omega_{7},\omega_{8})`.

* That means we have two properties to calculate (:math:`\alpha` and
  :math:`\beta`) so we should set ``num_props=2``.

.. nextslide::
.. include:: background.rst

The perturbation tuples for :math:`\alpha` and :math:`\beta` are
respectively ``{EL,EL}`` and ``{EL,EL,EL}``, so that

* ``len_tuple[2]={2,3}``, and
* ``pert_tuple[5]={EL,EL,EL,EL,EL}``, where the first two integer ``EL``'s
  are the perturbation tuple for :math:`\alpha` while last three for :math:`\beta`.

There are two frequency configurations for :math:`\alpha` (:math:`\omega_{1}`
and :math:`\omega_{2}`) and three for :math:`\beta` (:math:`\{\omega_{3},\omega_{4}\}`,
:math:`\{\omega_{5},\omega_{6}\}` and :math:`\{\omega_{7},\omega_{8}\}`),
which means:

* ``num_freq_configs[2]={2,3}``
* ``pert_freqs``

.. .. include:: home_icon.rst
