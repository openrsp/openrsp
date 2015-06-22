.. _chapter-callback-functions:

OpenRSP Callback Functions
==========================

To use OpenRSP, users should also prepare different callback functions
needed by OpenRSP. These callback functions will be invoked by OpenRSP
during calculations to get integral matrices or expectation values of
different one- and two-electron operators, and exchange-correlation
functionals, or to solve the response equation. The callback functions
are slightly different for C and Fortran users, which will be described
separately in this chapter.

It could be noted that the arguments in the following callback functions
are over complete. For instance, from the knowledge of ``len_tuple`` and
``pert_tuple``, the dimension of integral matrices ``num_int`` in the
callback function :c:func:`get_one_oper_mat` can be computed. The need
of this argument ``num_int`` is kind of technical issue, and we will give
detailed explanation in Section :ref:`section-openrsp-Fortran-APIs`.

OpenRSP Callback Functions (C version)
--------------------------------------

Examples of C callback functions can be found in the directory
``tests/c/callback``. The detailed information of these callback
functions are given as follows.

.. get_pert_comp()

.. get_pert_rank()

.. c:function:: QVoid get_overlap_mat(bra_len_tuple, bra_pert_tuple, ket_len_tuple, ket_pert_tuple, len_tuple, pert_tuple, user_ctx, num_int, val_int)

   Callback function for getting integral matrices of overlap integrals,
   the second last argument for the function :c:func:`OpenRSPSetPDBS`.

   :param bra_len_tuple: length of the perturbation tuple on the bra
   :type bra_len_tuple: QInt
   :param bra_pert_tuple: perturbation tuple on the bra,
       size is ``bra_len_tuple``
   :type bra_pert_tuple: QInt\*
   :param ket_len_tuple: length of the perturbation tuple on the ket
   :type ket_len_tuple: QInt
   :param ket_pert_tuple: perturbation tuple on the ket,
       size is ``ket_len_tuple``
   :type ket_pert_tuple: QInt\*
   :param len_tuple: length of perturbation tuple on the overlap integrals
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the overlap integrals, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the product of
       the sizes of perturbations on the bra, the ket and both of them,
       which are specified by the perturbation tuples on the bra
       (``bra_pert_tuple``), the ket (``ket_pert_tuple``) and overlap
       integrals (``pert_tuple``)
   :type num_int: QInt
   :var val_int: the integral matrices, arranged as
       ``(bra_pert_tuple, ket_pert_tuple, pert_tuple)``.
   :vartype val_int: QcMat\*[]
   :rtype: QVoid

.. c:function:: QVoid get_overlap_exp(bra_len_tuple, bra_pert_tuple, ket_len_tuple, ket_pert_tuple, len_tuple, pert_tuple, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   Callback function for getting expectation values of overlap integrals,
   the last argument for the function :c:func:`OpenRSPSetPDBS`.

   :param bra_len_tuple: length of the perturbation tuple on the bra
   :type bra_len_tuple: QInt
   :param bra_pert_tuple: perturbation tuple on the bra,
       size is ``bra_len_tuple``
   :type bra_pert_tuple: QInt\*
   :param ket_len_tuple: length of the perturbation tuple on the ket
   :type ket_len_tuple: QInt
   :param ket_pert_tuple: perturbation tuple on the ket,
       size is ``ket_len_tuple``
   :type ket_pert_tuple: QInt\*
   :param len_tuple: length of perturbation tuple on the overlap integrals
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the overlap integrals, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param num_dmat: number of atomic orbital (AO) based density matrices
   :type num_dmat: QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of number
       of density matrices (``num_dmat``) and the sizes of perturbations
       on the bra, the ket and overlap integrals
   :type num_exp: QInt
   :var val_exp: the expectation values, arranged as
       ``(num_dmat, bra_pert_tuple, ket_pert_tuple, pert_tuple)``.
   :vartype val_exp: QReal\*
   :rtype: QVoid

.. c:function:: QVoid get_one_oper_mat(len_tuple, pert_tuple, user_ctx, num_int, val_int)

   Callback function for getting integral matrices of a one-electron operator,
   the second last argument for the function :c:func:`OpenRSPAddOneOper`.

   :param len_tuple: length of perturbation tuple on the one-electron operator
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the one-electron operator, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the size of perturbations
       (specified by the perturbation tuple ``pert_tuple``)
   :type num_int: QInt
   :var val_int: the integral matrices
   :vartype val_int: QcMat\*[]
   :rtype: QVoid

.. c:function:: QVoid get_one_oper_exp(len_tuple, pert_tuple, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   Callback function for getting expectation values of a one-electron operator,
   the last argument for the function :c:func:`OpenRSPAddOneOper`.

   :param len_tuple: length of perturbation tuple on the one-electron operator
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the one-electron operator, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param num_dmat: number of AO based density matrices
   :type num_dmat: QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of number
       of density matrices (``num_dmat``) and the size of perturbations
       on the one-electron operator (specified by the perturbation tuple
       ``pert_tuple``)
   :type num_exp: QInt
   :var val_exp: the expectation values, arranged as ``(num_dmat, pert_tuple)``
   :vartype val_exp: QReal\*
   :rtype: QVoid

.. c:function:: QVoid get_two_oper_mat(len_tuple, pert_tuple, num_dmat, dens_mat, user_ctx, num_int, val_int)

   Callback function for getting integral matrices of a two-electron operator,
   the second last argument for the function :c:func:`OpenRSPAddTwoOper`.

   :param len_tuple: length of perturbation tuple on the two-electron operator
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the two-electron operator, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param num_dmat: number of AO based density matrices
   :type num_dmat: QInt
   :param dens_mat: the AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D})`
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the product of number
       of AO based density matrices (``num_dmat``) and the size of perturbations
       on the two-electron operator (specified by the perturbation tuple ``pert_tuple``)
   :type num_int: QInt
   :var val_int: the integral matrices, arranged as ``(num_dmat, pert_tuple)``
   :vartype val_int: QcMat\*[]
   :rtype: QVoid

.. c:function:: QVoid get_two_oper_exp(len_tuple, pert_tuple, len_dmat_tuple, num_LHS_dmat, LHS_dens_mat, num_RHS_dmat, RHS_dens_mat, user_ctx, num_exp, val_exp)

   Callback function for getting expectation values of a two-electron operator,
   the last argument for the function :c:func:`OpenRSPAddTwoOper`.

   :param len_tuple: length of perturbation tuple on the two-electron operator
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the two-electron operator, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param len_dmat_tuple: length of different perturbation tuples of the
       left-hand-side (LHS) and right-hand-side (RHS) AO based density
       matrices passed; for instance, if the LHS density matrices passed
       are (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{a}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{ab}`), and the
       RHS density matrices passed are (:math:`\boldsymbol{D}^{b}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{bc}`,
       :math:`\boldsymbol{D}^{d}`), then ``len_dmat_tuple`` equals to 4,
       and that means we want to calculate
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D})\boldsymbol{D}^{b}]`,
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D}^{a})\boldsymbol{D}^{c}]`,
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D}^{b})\boldsymbol{D}^{bc}]`,
       and :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D}^{ab})\boldsymbol{D}^{d}]`
   :type len_dmat_tuple: QInt
   :param num_LHS_dmat: number of LHS AO based density matrices passed for
       each LHS density matrix perturbation tuple, size is ``len_dmat_tuple``;
       sticking with the above example, ``num_LHS_dmat`` will be
       ``(1, N_a, N_b, N_ab)`` where ``N_a``, ``N_b`` and ``N_ab`` are
       respectively the numbers of density matrices for the density matrix
       perturbation tuples ``a``, ``b`` and ``ab``
   :type num_LHS_dmat: QInt\*
   :param LHS_dens_mat: the LHS AO based density matrices (:math:`\boldsymbol{D}_{\text{LHS}}`)
       for calculating
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D}_{\text{LHS}})\boldsymbol{D}_{\text{RHS}}]`,
       size is :math:`\sum_{\texttt{i}=0}^{\texttt{len\_dmat\_tuple}-1}` ``num_LHS_dmat[i]``
   :type LHS_dens_mat: QcMat\*[]
   :param num_RHS_dmat: number of RHS AO based density matrices passed for
       each RHS density matrix perturbation tuple, size is ``len_dmat_tuple``;
       sticking with the above example, ``num_RHS_dmat`` will be
       ``(N_b, N_c, N_bc, N_d)`` where ``N_b``, ``N_c`` ``N_bc`` and ``N_d``
       are respectively the numbers of density matrices for the density matrix
       perturbation tuples ``b``, ``c``, ``bc`` and ``d``
   :type num_RHS_dmat: QInt\*
   :param RHS_dens_mat: the RHS AO based density matrices (:math:`\boldsymbol{D}_{\text{RHS}}`)
       for calculating
       :math:`\mathrm{Tr}[\boldsymbol{G}^{\texttt{pert\_tuple}}(\boldsymbol{D}_{\text{LHS}})\boldsymbol{D}_{\text{RHS}}]`,
       size is :math:`\sum_{\texttt{i}=0}^{\texttt{len\_dmat\_tuple}-1}` ``num_RHS_dmat[i]``
   :type RHS_dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of the number
       of pairs of LHS and RHS density matrices and the size of perturbations
       on the two-electron operator (specified by the perturbation tuple
       ``pert_tuple``), the number of pairs of LHS and RHS density matrices
       can be computed as :math:`\sum_{\texttt{i}=0}^{\texttt{len\_dmat\_tuple}-1}`
       ``num_LHS_dmat[i]`` :math:`\times` ``num_RHS_dmat[i]``
   :type num_exp: QInt
   :var val_exp: the expectation values, arranged as
       ``(num_RHS_dmat, num_LHS_dmat, len_dmat_tuple, pert_tuple)``
   :vartype val_exp: QReal\*
   :rtype: QVoid

.. c:function:: QVoid get_xc_fun_mat(len_tuple, pert_tuple, num_freq_configs, len_dmat_tuple, idx_dmat_tuple, num_dmat, dens_mat, user_ctx, num_int, val_int)

   Callback function for getting integral matrices of XC functional,
   the second last argument for the function :c:func:`OpenRSPAddXCFun`.

   :param len_tuple: length of perturbation tuple on the XC functional
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the XC functional, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param num_freq_configs: the number of different frequency configurations
       to be considered for the perturbation tuple specified by ``pert_tuple``
   :type num_freq_configs: QInt
   :param len_dmat_tuple: the number of different perturbation tuples of the
       AO based density matrices passed; for instance, the complete density
       matrix perturbation tuples (canonically ordered) for a property
       :math:`\mathcal{E}^{abc}` (i.e. the perturbation tuple ``pert_tuple``
       is ``abc``) are (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{a}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{c}`,
       :math:`\boldsymbol{D}^{ab}`, :math:`\boldsymbol{D}^{ac}`,
       :math:`\boldsymbol{D}^{bc}`), and with the :math:`(0,2)` rule, the
       relevant density matrix perturbation tuples become (:math:`\boldsymbol{D}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{c}`,
       :math:`\boldsymbol{D}^{bc}`) that gives the ``len_dmat_tuple`` as 4
   :type len_dmat_tuple: QInt
   :param idx_dmat_tuple: indices of the density matrix perturbation tuples passed
       (canonically ordered), size is ``len_dmat_tuple``; sticking with the above
       example, the density matrix perturbation tuples passed are (:math:`\boldsymbol{D}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{bc}`)
       and their associated indices ``idx_dmat_tuple`` is (1, 3, 4, 7) because these
       numbers correspond to the positions of the ":math:`(k,n)`-surviving" perturbation
       tuples in the canonically ordered complete density matrix perturbation tuples
   :type idx_dmat_tuple: QInt\*
   :param num_dmat: number of collected AO based density matrices for the passed
       density matrix perturbation tuples (specified by ``idx_dmat_tuple``) and
       all frequency configurations, that is ``num_freq_configs``
       :math:`\times\sum_{\text{i}=0}^{\texttt{len\_dmat\_tuple}-1}N_{\text{i}}`,
       where :math:`N_{\text{i}}` is the number of density matrices for the
       density matrix perturbation tuple ``idx_dmat_tuple[i]`` for a frequency
       configuration
   :type num_dmat: QInt
   :param dens_mat: the collected AO based density matrices, size is ``num_dmat``,
       and arranged as ``(idx_dmat_tuple, num_freq_configs``)
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, equals to the product of
       the size of perturbations on the XC functional (specified by the
       perturbation tuple ``pert_tuple``) and the number of different frequency
       configurations ``num_freq_configs``
   :type num_int: QInt
   :var val_int: the integral matrices to be returned, size is ``num_int``,
       and arranged as (``pert_tuple``, ``num_freq_configs``)
   :vartype val_int: QcMat\*[]
   :rtype: QVoid

.. c:function:: QVoid get_xc_fun_exp(len_tuple, pert_tuple, num_freq_configs, len_dmat_tuple, idx_dmat_tuple, num_dmat, dens_mat, user_ctx, num_exp, val_exp)

   Callback function for getting expectation values of XC functional,
   the last argument for the function :c:func:`OpenRSPAddXCFun`.

   :param len_tuple: length of perturbation tuple on the XC functional
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the XC functional, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param num_freq_configs: the number of different frequency configurations
       to be considered for the perturbation tuple specified by ``pert_tuple``
   :type num_freq_configs: QInt
   :param len_dmat_tuple: the number of different perturbation tuples of the
       AO based density matrices passed
   :type len_dmat_tuple: QInt
   :param idx_dmat_tuple: indices of the density matrix perturbation tuples passed
       (canonically ordered), size is ``len_dmat_tuple``
   :type idx_dmat_tuple: QInt\*
   :param num_dmat: number of collected AO based density matrices for the passed
       density matrix perturbation tuples (specified by ``idx_dmat_tuple``) and
       all frequency configurations, that is ``num_freq_configs``
       :math:`\times\sum_{\text{i}=0}^{\texttt{len\_dmat\_tuple}-1}N_{\text{i}}`,
       where :math:`N_{\text{i}}` is the number of density matrices for the
       density matrix perturbation tuple ``idx_dmat_tuple[i]`` for a frequency
       configuration
   :type num_dmat: QInt
   :param dens_mat: the collected AO based density matrices, size is ``num_dmat``,
       and arranged as ``(idx_dmat_tuple, num_freq_configs``)
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of the expectation values, equals to the product of
       the size of perturbations on the XC functional (specified by the
       perturbation tuple ``pert_tuple``) and the number of different frequency
       configurations ``num_freq_configs``
   :type num_exp: QInt
   :var val_exp: the expectation values to be returned, size is ``num_exp``,
       and arranged as (``pert_tuple``, ``num_freq_configs``)
   :vartype val_exp: QReal\*
   :rtype: QVoid

.. c:function:: QVoid get_nuc_contrib(len_tuple, pert_tuple, user_ctx, size_pert, val_nuc)

   Callback function for getting the nuclear contributions, the last argument
   for the function :c:func:`OpenRSPAddNucContributions`.

   :param len_tuple: length of perturbation tuple on the nuclear Hamiltonian
   :type len_tuple: QInt
   :param pert_tuple: perturbation tuple on the nuclear Hamiltonian, size is
       ``len_tuple``
   :type pert_tuple: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param size_pert: size of the perturbations on the nuclear Hamiltonian,
       as specified by ``pert_tuple``
   :type size_pert: QInt
   :var val_nuc: the nuclear contributions, size is ``size_pert``
   :vartype val_nuc: QReal\*
   :rtype: QVoid

.. c:function:: QVoid get_linear_rsp_solution(size_pert, num_freq_sums, freq_sums, RHS_mat, user_ctx, rsp_param)

   Callback function for the linear response equation solver, the last argument
   for the function :c:func:`OpenRSPSetLinearRSPSolver`.

   :param size_pert: size of perturbations acting on the time-dependent
       self-consistent-field (TDSCF) equation
   :type size_pert: QInt
   :param num_freq_sums: number of complex frequency sums on the left hand side
       of the linear response equation
   :type num_freq_sums: QInt
   :param freq_sums: the complex frequency sums on the left hand side, size is
       ``2`` :math:`\times` ``num_freq_sums``, the real and imaginary parts of
       each frequency sum are consecutive in memory
   :type freq_sums: QReal\*
   :param RHS_mat: RHS matrices, size is ``size_pert`` :math:`\times`
       ``num_freq_sums``, and arranged as (``size_pert``, ``num_freq_sums``)
   :type RHS_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :var rsp_param: solved response parameters, size is ``size_pert`` :math:`\times`
       ``num_freq_sums``, and arranged as (``size_pert``, ``num_freq_sums``)
   :vartype rsp_param: QcMat\*[]
   :rtype: QVoid

.. Host programs will call OpenRSP by sending the excited states, so that we
   do not need the callback function get_rsp_eigen_solution
.. .. c:function:: QVoid get_rsp_eigen_solution(num_excit, eigen_val, user_ctx, eigen_vec)
 
    Callback function for the response eigenvalue equation solver, the last argument
    for the function :c:func:`OpenRSPSetRSPEigenSolver`.
 
    :param num_excit: number of excitations to be solved
    :type num_excit: QInt
    :param eigen_val: solved excitation energies, size is ``num_excit``
    :type eigen_val: QReal\*
    :param user_ctx: user-defined callback function context
    :type user_ctx: QVoid\*
    :var eigen_vec: eigenvectors solved from the eigenvalue problem,
        size is ``num_excit``
    :vartype eigen_vec: QcMat\*[]
    :rtype: QVoid

OpenRSP Callback Subroutines (Fortran version)
----------------------------------------------

The callback subroutines of Fortran codes take almost the exact arguments as
the callback functions of C codes. One difference is the type convention
between C and Fortran, which has been discussed in Secion :ref:`section-Fortran-convention`.
Moreover, the pointers of basic types (integer and real numbers) in the C
codes should be converted to corresponding array in Fortran. The array of
``QcMat`` pointers should be converted to an array of ``type(QcMat)`` in Fortran.
Last, the user-defined callback function/subroutine context should be replaced
by::

    integer, intent(in) :: len_ctx
    character(len=1), intent(in) :: user_ctx(len_ctx)

Examples of Fortran callback subroutines can be found in the directory
``tests/f90/callback``.
