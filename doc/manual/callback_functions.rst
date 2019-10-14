.. _subsection_callback_functions:

Callback Function Scheme
------------------------

To use OpenRSP, users should also prepare different callback functions
needed by OpenRSP. These callback functions will be invoked by OpenRSP
during calculations to get integral matrices or expectation values of
different one- and two-electron operators, exchange-correlation functionals
and nuclear contributions, or to solve the linear response equation.
The callback functions are slightly different for C and Fortran users,
which will be described separately in this chapter.

It should be noted that the arguments in the following callback functions are
over complete. For instance, from the knowledge of perturbations
(``oper_num_pert``, ``oper_pert_labels`` and ``oper_pert_orders``), the
dimension of integral matrices ``num_int`` in the callback function
:c:func:`get_one_oper_mat` can be computed.

.. The need of this argument ``num_int`` is kind of technical issue, and we will
.. give detailed explanation in Section :ref:`section_fortran_api_impl`.

**Last but not least, users should be aware that:**

#. OpenRSP always ask for **complex expectation values** for different one-
   and two-electron operators, exchange-correlation functionals and nuclear
   contributions, and these values are presented in memory that the real
   and imaginary parts of each value are consecutive. This affects:

   #. :c:func:`get_overlap_exp`
   #. :c:func:`get_one_oper_exp`
   #. :c:func:`get_two_oper_exp`
   #. :c:func:`get_xc_fun_exp`
   #. :c:func:`get_zero_oper_contrib`

#. In order to reduce the use of temporary matrices and values, OpenRSP
   requires that calculated integral matrices and expectation values
   should **be added to the returned argument**. OpenRSP will zero the
   entries of these matrices and expectation values at first. This
   requirement affects the callback functions of one- and two-electron
   operators, exchange-correlation functionals and nuclear contributions:

   #. :c:func:`get_overlap_mat` and :c:func:`get_overlap_exp`
   #. :c:func:`get_one_oper_mat` and :c:func:`get_one_oper_exp`
   #. :c:func:`get_two_oper_mat` and :c:func:`get_two_oper_exp`
   #. :c:func:`get_xc_fun_mat` and :c:func:`get_xc_fun_exp`
   #. :c:func:`get_zero_oper_contrib`

OpenRSP Callback Functions (C version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Examples of C callback functions can be found in these files
``tests/OpenRSP*Callback.c``. The detailed information of these callback
functions are given as follows.

.. c:function:: void get_pert_concatenation(pert_label, first_cat_comp, num_cat_comps, num_sub_tuples, len_sub_tuples, user_ctx, rank_sub_comps)

   User specified function for getting the ranks of components of
   sub-perturbation tuples (with the same perturbation label) for given
   components of the corresponding concatenated perturbation tuple, the last
   argument for the function :c:func:`OpenRSPSetPerturbations`.

   :param pert_label: the perturbation label
   :type pert_label: const QcPertInt
   :param first_cat_comp: rank of the first component of the concatenated
       perturbation tuple
   :type first_cat_comp: const QInt
   :param num_cat_comps: number of components of the concatenated perturbation
       tuple
   :type num_cat_comps: const QInt
   :param num_sub_tuples: number of sub-perturbation tuples to construct the
       concatenated perturbation tuple
   :type num_sub_tuples: const QInt
   :param len_sub_tuples: length of each sub-perturbation tuple, size is
       ``num_sub_tuples``; so that the length of the concatenated perturbation
       is ``sum(len_sub_tuples)``
   :type len_sub_tuples: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :var rank_sub_comps: ranks of components of sub-perturbation tuples for
       the corresponding component of the concatenated perturbation tuple,
       i.e. ``num_cat_comps`` components starting from the one with rank
       ``first_cat_comp``, size is therefore ``num_sub_tuples`` :math:`\times`
       ``num_cat_comps``, and arranged as ``[num_cat_comps][num_sub_tuples]``
   :vartype rank_sub_comps: QInt\*
   :rtype: void

**NOTE**: :c:func:`get_pert_concatenation` will not be invoked in the current
release so that users can use a "faked" function for it.

.. c:function:: void get_overlap_mat(bra_num_pert, bra_pert_labels, bra_pert_orders, ket_num_pert, ket_pert_labels, ket_pert_orders, oper_num_pert, oper_pert_labels, oper_pert_orders, user_ctx, num_int, val_int)

   User-specified callback function to calculate integral matrices of overlap
   operator as well as its derivatives with respect to different perturbations,
   the second last argument for the function :c:func:`OpenRSPSetOverlap`.

   :param bra_num_pert: number of perturbations on the bra center
   :type bra_num_pert: const QInt
   :param bra_pert_labels: labels of perturbations on the bra center,
       size is ``bra_num_pert``
   :type bra_pert_labels: const QcPertInt\*
   :param bra_pert_orders: orders of perturbations on the bra center,
       size is ``bra_num_pert``
   :type bra_pert_orders: const QInt\*
   :param ket_num_pert: number of perturbations on the ket center
   :type ket_num_pert: const QInt
   :param ket_pert_labels: labels of perturbations on the ket center,
       size is ``ket_num_pert``
   :type ket_pert_labels: const QcPertInt\*
   :param ket_pert_orders: orders of perturbations on the ket center,
       size is ``ket_num_pert``
   :type ket_pert_orders: const QInt\*
   :param oper_num_pert: number of perturbations on the overlap operator [#]_
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the overlap operator,
       size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the overlap operator,
       size is ``oper_num_pert`` [#]_
   :type oper_pert_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_int: number of the integral matrices, as the product of the sizes
       of perturbations on the bra, the ket and the overlap operator
   :type num_int: const QInt
   :var val_int: the integral matrices to be added, size is ``num_int``, and
       arranged as ``[oper_pert][bra_pert][ket_pert]``
   :vartype val_int: QcMat\*[]
   :rtype: void

.. [#] Here perturbations on the overlap operator represent those acting on the
       whole integral of the overlap operator, i.e. they can act on either the
       bra center or the ket center by applying the rule of derivatives of a
       product.
.. [#] Only overlap integrals perturbed on the bra and/or the ket, and those
       perturbed on the whole integral are needed in the calculations. It means
       that, OpenRSP will only ask for overlap integrals either with
       perturbations on the bra and/or ket (``oper_num_pert=0``), or with
       perturbations on the whole overlap integral (``bra_num_pert=0`` and
       ``ket_num_pert=0``).

.. c:function:: void get_overlap_exp(bra_num_pert, bra_pert_labels, bra_pert_orders, ket_num_pert, ket_pert_labels, ket_pert_orders, oper_num_pert, oper_pert_labels, oper_pert_orders, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   User-specified function for calculating expectation values of the overlap
   operator and its derivatives, the last argument for the function
   :c:func:`OpenRSPSetOverlap`.

   :param bra_num_pert: number of perturbations on the bra center
   :type bra_num_pert: const QInt
   :param bra_pert_labels: labels of perturbations on the bra center,
       size is ``bra_num_pert``
   :type bra_pert_labels: const QcPertInt\*
   :param bra_pert_orders: orders of perturbations on the bra center,
       size is ``bra_num_pert``
   :type bra_pert_orders: const QInt\*
   :param ket_num_pert: number of perturbations on the ket center
   :type ket_num_pert: const QInt
   :param ket_pert_labels: labels of perturbations on the ket center,
       size is ``ket_num_pert``
   :type ket_pert_labels: const QcPertInt\*
   :param ket_pert_orders: orders of perturbations on the ket center,
       size is ``ket_num_pert``
   :type ket_pert_orders: const QInt\*
   :param oper_num_pert: number of perturbations on the overlap operator [#]_
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the overlap operator,
       size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the overlap operator,
       size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param num_dmat: number of atomic orbital (AO) based density matrices
   :type num_dmat: const QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_exp: number of the expectation values, as the product of sizes of
       perturbations on the bra, the ket, the overlap operator and the number
       of density matrices (``num_dmat``)
   :type num_exp: const QInt
   :var val_exp: the expectation values to be added, size is ``2``
       :math:`\times` ``num_exp``, and arranged as
       ``[num_dmat][oper_pert][bra_pert][ket_pert][2]``
   :vartype val_exp: QReal\*
   :rtype: void

.. [#] Similar to the callback function :c:func:`get_overlap_mat`, OpenRSP will
       only ask for expectation values either with perturbations on the bra
       and/or ket (``oper_num_pert=0``), or with perturbations on the whole
       overlap integral (``bra_num_pert=0`` and ``ket_num_pert=0``).

.. c:function:: void get_one_oper_mat(oper_num_pert, oper_pert_labels, oper_pert_orders, user_ctx, num_int, val_int)

   User-specified function for calculating integral matrices of the
   one-electron operator and its derivatives, the second last argument for the
   function :c:func:`OpenRSPAddOneOper`.

   :param oper_num_pert: number of perturbations on the one-electron operator
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the one-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the one-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_int: number of the integral matrices, as the size of
       perturbations that are specified by ``oper_num_pert``,
       ``oper_pert_labels`` and ``oper_pert_orders``
   :type num_int: const QInt
   :var val_int: the integral matrices to be added, size is ``num_int``
   :vartype val_int: QcMat\*[]
   :rtype: void

.. c:function:: void get_one_oper_exp(oper_num_pert, oper_pert_labels, oper_pert_orders, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   User-specified callback function to calculate expectation values of
   one-electron operator as well as its derivatives with respect to different
   perturbations, the last argument for the function
   :c:func:`OpenRSPAddOneOper`.

   :param oper_num_pert: number of perturbations on the one-electron operator
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the one-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the one-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param num_dmat: number of AO based density matrices
   :type num_dmat: const QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_exp: number of expectation values, as the product of the size of
       perturbations on the one-electron operator (specified by
       ``oper_num_pert``, ``oper_pert_labels`` and ``oper_pert_orders``) and
       the number of density matrices (``num_dmat``)
   :type num_exp: const QInt
   :var val_exp: the expectation values to be added, size is ``2``
       :math:`\times` ``num_exp``, and arranged as ``[num_dmat][oper_pert][2]``
   :vartype val_exp: QReal\*
   :rtype: void

.. c:function:: void get_two_oper_mat(oper_num_pert, oper_pert_labels, oper_pert_orders, num_dmat, dens_mat, user_ctx, num_int, val_int)

   User-specified function for calculating integral matrices of the
   two-electron operator and its derivatives, the second last argument for the
   function :c:func:`OpenRSPAddTwoOper`.

   :param oper_num_pert: number of perturbations on the two-electron operator
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the two-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the two-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param num_dmat: number of AO based density matrices
   :type num_dmat: const QInt
   :param dens_mat: the AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating
       :math:`\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D})`,
       where :math:`\texttt{perturbations}` are specified by ``oper_num_pert``,
       ``oper_pert_labels`` and ``oper_pert_orders``.
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_int: number of the integral matrices, as the product of the size
       of perturbations on the two-electron operator (specified by
       ``oper_num_pert``, ``oper_pert_labels`` and ``oper_pert_orders``) and
       the number of AO based density matrices (``num_dmat``)
   :type num_int: const QInt
   :var val_int: the integral matrices to be added, size is ``num_int``,
       and arranged as ``[num_dmat][oper_pert]``
   :vartype val_int: QcMat\*[]
   :rtype: void

.. c:function:: void get_two_oper_exp(oper_num_pert, oper_pert_labels, oper_pert_orders, dmat_len_tuple, num_LHS_dmat, LHS_dens_mat, num_RHS_dmat, RHS_dens_mat, user_ctx, num_exp, val_exp)

   User-specified callback function to calculate expectation values of
   two-electron operator as well as its derivatives with respect to different
   perturbations, the last argument for the function
   :c:func:`OpenRSPAddTwoOper`.

   :param oper_num_pert: number of perturbations on the two-electron operator
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the two-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the two-electron
       operator, size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param dmat_len_tuple: length of different perturbation tuples of the
       left-hand-side (LHS) and right-hand-side (RHS) AO based density
       matrices passed; for instance, if the LHS density matrices passed
       are (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{a}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{ab}`), and the
       RHS density matrices passed are (:math:`\boldsymbol{D}^{b}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{bc}`,
       :math:`\boldsymbol{D}^{d}`), then ``dmat_len_tuple`` equals to `4`,
       and that means we want to calculate
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D})\boldsymbol{D}^{b}]`,
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D}^{a})\boldsymbol{D}^{c}]`,
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D}^{b})\boldsymbol{D}^{bc}]`,
       and :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D}^{ab})\boldsymbol{D}^{d}]`,
       where :math:`\texttt{perturbations}` are specified by ``oper_num_pert``,
       ``oper_pert_labels`` and ``oper_pert_orders``.
   :type dmat_len_tuple: const QInt
   :param num_LHS_dmat: number of LHS AO based density matrices passed for
       each LHS density matrix perturbation tuple, size is ``dmat_len_tuple``;
       sticking with the above example, ``num_LHS_dmat`` will be
       ``{1, N_a, N_b, N_ab}`` where ``N_a``, ``N_b`` and ``N_ab`` are
       respectively the numbers of density matrices for the density matrix
       perturbation tuples ``a``, ``b`` and ``ab``
   :type num_LHS_dmat: const QInt\*
   :param LHS_dens_mat: the LHS AO based density matrices (:math:`\boldsymbol{D}_{\text{LHS}}`)
       for calculating
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D}_{\text{LHS}})\boldsymbol{D}_{\text{RHS}}]`,
       size is the sum of ``num_LHS_dmat``
   :type LHS_dens_mat: QcMat\*[]
   :param num_RHS_dmat: number of RHS AO based density matrices passed for
       each RHS density matrix perturbation tuple, size is ``dmat_len_tuple``;
       sticking with the above example, ``num_RHS_dmat`` will be
       ``{N_b, N_c, N_bc, N_d}`` where ``N_b``, ``N_c`` ``N_bc`` and ``N_d``
       are respectively the numbers of density matrices for the density matrix
       perturbation tuples ``b``, ``c``, ``bc`` and ``d``
   :type num_RHS_dmat: const QInt\*
   :param RHS_dens_mat: the RHS AO based density matrices (:math:`\boldsymbol{D}_{\text{RHS}}`)
       for calculating
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{perturbations}}(\boldsymbol{D}_{\text{LHS}})\boldsymbol{D}_{\text{RHS}}]`,
       size is the sum of ``num_RHS_dmat``
   :type RHS_dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_exp: number of expectation values, as the product of the size
       of perturbations on the two-electron operator (specified by
       ``oper_num_pert``, ``oper_pert_labels`` and ``oper_pert_orders``) and
       the number of pairs of LHS and RHS density matrices, and the number of
       pairs of LHS and RHS density matrices can be computed as the dot product
       of ``num_LHS_dmat`` and ``num_RHS_dmat``
   :type num_exp: const QInt
   :var val_exp: the expectation values to be added, size is ``2``
       :math:`\times` ``num_exp``, and arranged as
       ``[dmat_len_tuple][num_LHS_dmat][num_RHS_dmat][oper_pert][2]``
   :vartype val_exp: QReal\*
   :rtype: void

.. c:function:: void get_xc_fun_mat(xc_len_tuple, xc_pert_tuple, num_freq_configs, pert_freq_category, dmat_num_tuple, dmat_idx_tuple, num_dmat, dens_mat, user_ctx, num_int, val_int)

   User-specified function for calculating integral matrices of the XC
   functional and its derivatives, the second last argument for the function
   :c:func:`OpenRSPAddXCFun`.

   :param xc_len_tuple: length of the perturbation tuple on the XC functional
   :type xc_len_tuple: const QInt
   :param xc_pert_tuple: perturbation tuple on the XC functional, size is
       ``xc_len_tuple``
   :type xc_pert_tuple: const QcPertInt\*
   :param num_freq_configs: the number of different frequency configurations to
       be considered for the perturbation tuple specified by ``xc_pert_tuple``
   :type num_freq_configs: const QInt
   :param pert_freq_category: category of perturbation frequencies, size is
       ``[num_freq_configs][xc_len_tuple]``. Take :math:`\mathcal{E}^{gfff}` as an
       example, suppose we have four different frequency configurations:
       "0.0,0.0,0.0,0.0" (:math:`3N\times 10` unique elements),
       "0.0,-0.2,0.1,0.1" (:math:`3N\times 18` unique elements),
       "0.0,-0,3,0.1,0.2" (:math:`3N\times 27` unique elements) and
       "0.0,-0.1,0.1,0.0" (:math:`3N\times 27` unique elements), the
       ``pert_freq_category`` argument would then be ``(1,1,1,1, 1,2,3,3,
       1,2,3,4, 1,2,3,1)``.
   :type pert_freq_category: const QInt\*
   :param dmat_num_tuple: the number of different perturbation tuples of the
       AO based density matrices passed; for instance, the complete density
       matrix perturbation tuples (canonically ordered) for a property
       :math:`\mathcal{E}^{abc}` (i.e. the perturbation tuple ``xc_pert_tuple``
       is ``abc``) are (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{a}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{ab}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{ac}`,
       :math:`\boldsymbol{D}^{bc}`), and with the :math:`(0,2)` rule, the
       relevant density matrix perturbation tuples become (:math:`\boldsymbol{D}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{c}`,
       :math:`\boldsymbol{D}^{bc}`) that gives the ``dmat_num_tuple`` as `4`
   :type dmat_num_tuple: const QInt
   :param dmat_idx_tuple: indices of the density matrix perturbation tuples
       passed (canonically ordered), size is ``dmat_num_tuple``; sticking with
       the above example, the density matrix perturbation tuples passed are
       (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{b}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{bc}`) and their
       associated indices ``dmat_idx_tuple`` is ``{1, 3, 5, 7}`` because these
       numbers correspond to the positions of the ":math:`(k,n)`-surviving"
       perturbation tuples in the canonically ordered complete density matrix
       perturbation tuples
   :type dmat_idx_tuple: const QInt\*
   :param num_dmat: number of collected AO based density matrices for the
       passed density matrix perturbation tuples (specified by
       ``dmat_idx_tuple``) and all frequency configurations, that is
       ``num_freq_configs`` :math:`\times\sum_{\text{i}}N_{\text{i}}`, where
       :math:`N_{\text{i}}` is the number of density matrices for the density
       matrix perturbation tuple ``dmat_idx_tuple[i]`` for a frequency
       configuration
   :type num_dmat: const QInt
   :param dens_mat: the collected AO based density matrices, size is
       ``num_dmat``, and arranged as ``[num_freq_configs][dmat_idx_tuple]``
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_int: number of the integral matrices, equals to the product of
       the size of perturbations on the XC functional (specified by the
       perturbation tuple ``xc_pert_tuple``) and the number of different
       frequency configurations ``num_freq_configs``
   :type num_int: const QInt
   :var val_int: the integral matrices to be added, size is ``num_int``, and
       arranged as ``[num_freq_configs][xc_pert_tuple]``
   :vartype val_int: QcMat\*[]
   :rtype: void

.. c:function:: void get_xc_fun_exp(xc_len_tuple, xc_pert_tuple, num_freq_configs, pert_freq_category, dmat_num_tuple, dmat_idx_tuple, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   User-specified function for calculating expectation values of the XC
   functional and its derivatives, the last argument for the function
   :c:func:`OpenRSPAddXCFun`.

   :param xc_len_tuple: length of the perturbation tuple on the XC functional
   :type xc_len_tuple: const QInt
   :param xc_pert_tuple: perturbation tuple on the XC functional, size is
       ``xc_len_tuple``
   :type xc_pert_tuple: const QcPertInt\*
   :param num_freq_configs: the number of different frequency configurations to
       be considered for the perturbation tuple specified by ``xc_pert_tuple``
   :type num_freq_configs: const QInt
   :param pert_freq_category: category of perturbation frequencies, size is
       ``[num_freq_configs][xc_len_tuple]``.
   :type pert_freq_category: const QInt\*
   :param dmat_num_tuple: the number of different perturbation tuples of the
       AO based density matrices passed
   :type dmat_num_tuple: const QInt
   :param dmat_idx_tuple: indices of the density matrix perturbation tuples
       passed (canonically ordered), size is ``dmat_num_tuple``
   :type dmat_idx_tuple: const QInt\*
   :param num_dmat: number of collected AO based density matrices for the
       passed density matrix perturbation tuples (specified by
       ``dmat_idx_tuple``) and all frequency configurations, that is
       ``num_freq_configs`` :math:`\times\sum_{\text{i}}N_{\text{i}}`, where
       :math:`N_{\text{i}}` is the number of density matrices for the density
       matrix perturbation tuple ``dmat_idx_tuple[i]`` for a frequency
       configuration
   :type num_dmat: const QInt
   :param dens_mat: the collected AO based density matrices, size is
       ``num_dmat``, and arranged as ``[num_freq_configs][dmat_idx_tuple]``
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param num_exp: number of the expectation values, equals to the product of
       the size of perturbations on the XC functional (specified by the
       perturbation tuple ``xc_pert_tuple``) and the number of different
       frequency configurations ``num_freq_configs``
   :type num_exp: const QInt
   :var val_exp: the expectation values to be added, size is ``2``
       :math:`\times` ``num_exp``, and arranged as
       ``[num_freq_configs][xc_pert_tuple][2]``
   :vartype val_exp: QReal\*
   :rtype: void

.. c:function:: void get_zero_oper_contrib(oper_num_pert, oper_pert_labels, oper_pert_orders, user_ctx, size_pert, val_oper)

   User-specified callback function to calculate contributions from the
   zero-electron operator, the last argument for the function
   :c:func:`OpenRSPAddZeroOper`.

   :param oper_num_pert: number of perturbations on the zero-electron operator
   :type oper_num_pert: const QInt
   :param oper_pert_labels: labels of perturbations on the zero-electron operator,
       size is ``oper_num_pert``
   :type oper_pert_labels: const QcPertInt\*
   :param oper_pert_orders: orders of perturbations on the zero-electron operator,
       size is ``oper_num_pert``
   :type oper_pert_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param size_pert: size of the perturbations on the zero-electron operator
   :type size_pert: const QInt
   :var val_oper: contributions from the zero-electron operator to be added,
       arranged as ``[size_pert][2]``
   :vartype val_oper: QReal\*
   :rtype: void

.. c:function:: void get_linear_rsp_solution(num_pert, num_comps, num_freq_sums, freq_sums, RHS_mat, user_ctx, rsp_param)

   User-specified callback function of linear response equation solver, the
   last argument for the function :c:func:`OpenRSPSetLinearRSPSolver`.

   :param num_pert: number of different perturbations on the right hand side of
       the linear response equation
   :type num_pert: const QInt
   :param num_comps: number of components of each perturbation, size is
       ``num_pert``
   :type num_comps: const QInt\*
   :param num_freq_sums: for each perturbation, number of complex frequency
       sums on the left hand side of the linear response equation, size is
       ``num_pert``
   :type num_freq_sums: const QInt\*
   :param freq_sums: the complex frequency sums on the left hand side of the
       linear response equation, size is twice of the sum of ``num_freq_sums``,
       the real and imaginary parts of each frequency sum are consecutive in
       memory
   :type freq_sums: const QReal\*
   :param RHS_mat: RHS matrices, size is the dot product of ``num_comps`` and
       ``num_freq_sums``, and index of ``num_freq_sums`` runs faster in memory
   :type RHS_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :var rsp_param: solved response parameters, size is the dot product of
       ``num_comps`` and ``num_freq_sums``, and index of ``num_freq_sums`` runs
       faster in memory
   :vartype rsp_param: QcMat\*[]
   :rtype: void

OpenRSP Callback Subroutines (Fortran version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The callback subroutines of Fortran codes take almost the exact arguments as
the callback functions of C codes. One difference is the type convention
between C and Fortran, which has been discussed in Secion
:ref:`subsubsection_fortran_convention`.  Moreover, the pointers of basic types
(integer and real numbers) in the C codes should be converted to corresponding
array in Fortran. The array of ``QcMat`` pointers should be converted to an
array of ``type(QcMat)`` in Fortran. Last, the user-defined callback
function/subroutine context should be replaced by ``type(C_PTR)``.

We will develop Fortran unit testing in next release. For the time being,
interested users can refer to LSDalton for examples of Fortran callback
subroutines.
