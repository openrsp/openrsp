.. _chapter_response_functions:

Response Functions
==================

The response functions can be calculated by calling:

.. c:function:: QErrorCode OpenRSPGetRSPFun(open_rsp, ref_ham, ref_state, ref_overlap, num_props, len_tuple, pert_tuple, num_freq_configs, pert_freqs, kn_rules, r_flag, write_threshold, size_rsp_funs, rsp_funs)

Suppose we want to calculate the polarizability :math:`\alpha` and the first
hyperpolarizability :math:`\beta` at different frequencies:

:math:`\alpha(-\omega_{1},\omega_{1})`, :math:`\alpha(-\omega_{2},\omega_{2})`,
:math:`\beta(-\omega_{3}-\omega_{4},\omega_{3},\omega_{4})`,
:math:`\beta(-\omega_{5}-\omega_{6},\omega_{5},\omega_{6})`,
:math:`\beta(-\omega_{7}-\omega_{8},\omega_{7},\omega_{8})`.

* That means we have two properties to calculate (:math:`\alpha` and
  :math:`\beta`) so we should set ``num_props=2``.

The perturbation tuples for :math:`\alpha` and :math:`\beta` are
presented by integers, let them be ``{EL,EL}`` and ``{EL,EL,EL}``,
so that

* ``len_tuple[2]={2,3}``, and
* ``pert_tuple[5]={EL,EL,EL,EL,EL}``, where the first two are for
  :math:`\alpha` and the last three for :math:`\beta`.

There are two frequency configurations for :math:`\alpha` (:math:`\omega_{1}`
and :math:`\omega_{2}`) and three for :math:`\beta`
(:math:`\{\omega_{3},\omega_{4}\}`, :math:`\{\omega_{5},\omega_{6}\}` and
:math:`\{\omega_{7},\omega_{8}\}`), as such:

* ``num_freq_configs[2]={2,3}``,
* ``pert_freqs[8]=``:math:`\{\omega_{1},\omega_{2},\omega_{3},\omega_{4},\omega_{5},\omega_{6},\omega_{7},\omega_{8}\}`,

and the frequency of perturbation :math:`a` is not needed.

The argument ``kn_rules`` contains the number :math:`k` of :math:`(k,n)` rule
for each property. If we choose :math:`(0,1)` and :math:`(1,1)` rules for
:math:`\alpha` and :math:`\beta` respectively, we have ``kn_rules[2]={0,1}``.

The choice of appropriate :math:`(k,n)` rule usually affect the efficiency of
calculations. Detailed discussion of the :math:`(k,n)` rule can be found, for
instance in:

#. Magnus Ringholm, *et al.*, J. Comput. Chem. 35, 622 (2014).
#. Andreas J. Thorvaldsen, *et al.*, J. Chem. Phys. 129, 214108 (2008).
#. Kasper Kristensen, *et al.*, J. Chem. Phys. 129, 214103 (2008).

``r_flag`` and ``write_threshold`` respectively controls the restarting scheme
and threshold of tensor element writing (by OpenRSP). Please refer to
:ref:`chapter_api_reference` for more detail.

The calculated results are in ``rsp_funs``, whose size is twice of
``size_rsp_funs`` (because OpenRSP represents a complex number by its real and
imaginary parts).

The ``size_rsp_funs`` equals to the sum of the size of each property to
calculate. In this example, if we use non-redundant representation of
perturbations, :math:`\alpha` will have 6 components (:math:`\alpha_{xx}`,
:math:`\alpha_{xy}`, :math:`\alpha_{xz}`, :math:`\alpha_{yy}`,
:math:`\alpha_{yz}`, :math:`\alpha_{zz}`), and :math:`\beta` 10 components
(:math:`\beta_{xxx}`, :math:`\beta_{xxy}`, :math:`\beta_{xxz}`,
:math:`\beta_{xyy}`, :math:`\beta_{xyz}`, :math:`\beta_{xzz}`,
:math:`\beta_{yyy}`, :math:`\beta_{yyz}`, :math:`\beta_{yzz}`,
:math:`\beta_{zzz}`).

Further considering the number of frequency configurations of :math:`\alpha`
(2) and :math:`\beta` (3), we have ``size_rsp_funs`` as

:math:`6\times2` (for :math:`\alpha`) :math:`+` :math:`10\times3` (for
:math:`\beta`) :math:`=42`.

The results ``rsp_funs`` are in a one-dimensional array, and are arranged in
memory as:

``[num_props][num_freq_configs][pert_tuple][2]``,

that is (where the frequency of perturbation :math:`a` is neglected):

| ``[`` :math:`\mathrm{Re}(\alpha_{xx}(\omega_{1}))`,
  :math:`\mathrm{Im}(\alpha_{xx}(\omega_{1}))`,
  :math:`\cdots`,
  :math:`\mathrm{Re}(\alpha_{zz}(\omega_{1}))`,
  :math:`\mathrm{Im}(\alpha_{zz}(\omega_{1}))`,
| :math:`\mathrm{Re}(\alpha_{xx}(\omega_{2}))`,
  :math:`\mathrm{Im}(\alpha_{xx}(\omega_{2}))`,
  :math:`\cdots`,
  :math:`\mathrm{Re}(\alpha_{zz}(\omega_{2}))`,
  :math:`\mathrm{Im}(\alpha_{zz}(\omega_{2}))`,
| :math:`\mathrm{Re}(\beta_{xxx}(\omega_{3},\omega_{4}))`,
  :math:`\mathrm{Im}(\beta_{xxx}(\omega_{3},\omega_{4}))`,
  :math:`\cdots`,
  :math:`\mathrm{Re}(\beta_{zzz}(\omega_{3},\omega_{4}))`, 
  :math:`\mathrm{Im}(\beta_{zzz}(\omega_{3},\omega_{4}))`,
| :math:`\cdots`,
| :math:`\mathrm{Re}(\beta_{xxx}(\omega_{7},\omega_{8}))`, 
  :math:`\mathrm{Im}(\beta_{xxx}(\omega_{7},\omega_{8}))`, 
  :math:`\cdots`,
  :math:`\mathrm{Re}(\beta_{zzz}(\omega_{7},\omega_{8}))`,
  :math:`\mathrm{Im}(\beta_{zzz}(\omega_{7},\omega_{8}))` ``]``.

Last but not least, to perform such calculations, users have to provide OpenRSP
the knowledge of reference state (usually the ground state) by setting the
``ref_ham``, ``ref_state`` and ``ref_overlap``, which are ``struct QcMat*`` in
C and ``type(QcMat)`` in Fortran. They are respectively the Hamiltonian (Fock
matrix), atomic orbital based density matrix and overlap integral matrix.
