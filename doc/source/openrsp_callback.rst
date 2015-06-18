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

OpenRSP Callback Functions (C version)
--------------------------------------

Examples of C callback functions can be found in the directory
``tests/c/callback``. The detailed information of these callback
functions are given as follows.

.. get_pert_comp()

.. get_pert_rank()

.. function:: get_overlap_mat(bra_num_pert,    \
                              bra_pert_labels, \
                              ket_num_pert,    \
                              ket_pert_labels, \
                              num_pert,        \
                              pert_labels,     \
                              user_ctx,        \
                              num_int,         \
                              val_int)

   Callback function for getting integral matrices of overlap integrals,
   the second last argument for the function :py:meth:`OpenRSPSetPDBS`.

   :param bra_num_pert: number of perturbations on the bra
   :type bra_num_pert: QInt
   :param bra_pert_labels: label for each perturbation on the bra,
       size is ``bra_num_pert``
   :type bra_pert_labels: QInt\*
   :param ket_num_pert: number of perturbations on the ket
   :type ket_num_pert: QInt
   :param ket_pert_labels: label for each perturbation on the ket,
       size is ``ket_num_pert``
   :type ket_pert_labels: QInt\*
   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the product of
       the dimension of each perturbation pattern on the bra (``bra_pert_labels``),
       ket (``ket_pert_labels``) and both of them (``pert_labels``)
   :type num_int: QInt
   :param val_int: the integral matrices, ordered as
       (``perturbation pattern on the bra``, ``perturbation pattern on the ket``,
       ``perturbation pattern``), where the perturbations on the bra vary fastest
       in memory
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_overlap_exp(bra_num_pert,    \
                              bra_pert_labels, \
                              ket_num_pert,    \
                              ket_pert_labels, \
                              num_pert,        \
                              pert_labels,     \
                              num_dmat,        \
                              dens_mat,        \
                              user_ctx,        \
                              num_exp,         \
                              val_exp)

   Callback function for getting expectation values of overlap integrals,
   the last argument for the function :py:meth:`OpenRSPSetPDBS`.

   :param bra_num_pert: number of perturbations on the bra
   :type bra_num_pert: QInt
   :param bra_pert_labels: label for each perturbation on the bra,
       size is ``bra_num_pert``
   :type bra_pert_labels: QInt\*
   :param ket_num_pert: number of perturbations on the ket
   :type ket_num_pert: QInt
   :param ket_pert_labels: label for each perturbation on the ket,
       size is ``ket_num_pert``
   :type ket_pert_labels: QInt\*
   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_dmat: number of atomic orbital (AO) based density matrices
   :type num_dmat: QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of number
       of density matrices (``num_dmat``) and the dimension of each perturbation
       pattern on the bra (``bra_pert_labels``), ket (``ket_pert_labels``)
       and both of them (``pert_labels``)
   :type num_exp: QInt
   :param val_exp: the expectation values, ordered as
       (``perturbation pattern on the bra``, ``perturbation pattern on the ket``,
       ``perturbation pattern``, ``num_dmat``)
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_one_oper_mat(num_pert,    \
                               pert_labels, \
                               user_ctx,    \
                               num_int,     \
                               val_int)

   Callback function for getting integral matrices of a one-electron operator,
   the second last argument for the function :py:meth:`OpenRSPAddOneOper`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the dimension of
       the perturbation pattern (``pert_labels``)
   :type num_int: QInt
   :param val_int: the integral matrices
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_one_oper_exp(num_pert,    \
                               pert_labels, \
                               num_dmat,    \
                               dens_mat,    \
                               user_ctx,    \
                               num_exp,     \
                               val_exp)

   Callback function for getting expectation values of a one-electron operator,
   the last argument for the function :py:meth:`OpenRSPAddOneOper`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_dmat: number of AO based density matrices
   :type num_dmat: QInt
   :param dens_mat: the AO based density matrices
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of the
       dimension of the perturbation pattern (``pert_labels``) and the
       number of density matrices (``num_dmat``)
   :type num_exp: QInt
   :param val_exp: the expectation values, ordered as (``perturbation pattern``,
       ``num_dmat``)
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_two_oper_mat(num_pert,     \
                               pert_labels,  \
                               num_var_dmat, \
                               var_dens_mat, \
                               user_ctx,     \
                               num_int,      \
                               val_int)

   Callback function for getting integral matrices of a two-electron operator,
   the second last argument for the function :py:meth:`OpenRSPAddTwoOper`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_var_dmat: number of variable AO based density matrices
   :type num_var_dmat: QInt
   :param var_dens_mat: the variable AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\boldsymbol{G}(\boldsymbol{D})`
   :type var_dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, as the product of
       the dimension of perturbation pattern (``pert_labels``) and the
       number of variable AO based density matrices (``num_var_dmat``)
   :type num_int: QInt
   :param val_int: the integral matrices, ordered as (``perturbation pattern``,
       ``num_var_dmat``)
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_two_oper_exp(num_pert,       \
                               pert_labels,    \
                               num_var_dmat,   \
                               var_dens_mat,   \
                               num_contr_dmat, \
                               contr_dens_mat, \
                               user_ctx,       \
                               num_exp,        \
                               val_exp)

   Callback function for getting expectation values of a two-electron operator,
   the last argument for the function :py:meth:`OpenRSPAddTwoOper`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_var_dmat: number of variable AO based density matrices
   :type num_var_dmat: QInt
   :param var_dens_mat: the variable AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\boldsymbol{G}(\boldsymbol{D})`
   :type var_dens_mat: QcMat\*[]
   :param num_contr_dmat: number of contracted AO based density matrices
   :type num_contr_dmat: QInt
   :param contr_dens_mat: the contracted AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\mathrm{Tr}[\boldsymbol{G}\boldsymbol{D}]`
   :type contr_dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values, as the product of
       the dimension of perturbation pattern (``pert_labels``), the
       number of variable AO based density matrices (``num_var_dmat``)
       and the number of contracted AO based density matrices (``num_contr_dmat``)
   :type num_exp: QInt
   :param val_exp: the expectation values, ordered as (``perturbation pattern``,
       ``num_var_dmat``, ``num_contr_dmat``)
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_xc_fun_mat(num_pert,         \
                             pert_labels,      \
                             num_freq_configs, \
                             len_dmat_tuple,   \
                             dens_mat_tuple,   \
                             num_dmat,         \
                             dens_mat,         \
                             user_ctx,         \
                             num_int,          \
                             val_int)

   Callback function for getting integral matrices of XC functional,
   the second last argument for the function :py:meth:`OpenRSPAddXCFun`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_freq_configs: for the perturbation pattern specified by ``num_pert``
       and ``pert_labels``, the number of different frequency configurations to
       be considered
   :type num_freq_configs: QInt
   :param len_dmat_tuple: for the perturbation pattern specified by ``num_pert``
       and ``pert_labels``, the length (number of different perturbation patterns)
       of the AO based density matrices passed; for instance, the complete density
       matrix tuple (canonically ordered) for a property :math:`\mathcal{E}^{abc}`
       is (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{a}`, :math:`\boldsymbol{D}^{b}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{ab}`, :math:`\boldsymbol{D}^{ac}`,
       :math:`\boldsymbol{D}^{bc}`), and with the (0,2) rule, the relevant density
       matrices are (:math:`\boldsymbol{D}`, :math:`\boldsymbol{D}^{b}`,
       :math:`\boldsymbol{D}^{c}`, :math:`\boldsymbol{D}^{bc}`) and which gives the
       ``len_dmat_tuple`` as 4
   :type len_dmat_tuple: QInt
   :param dens_mat_tuple: the perturbation tuple of the AO based density matrices
       passed, as a canonically ordered list of all relevant perturbation patterns
       of the density matrices, size is ``len_dmat_tuple``; sticking with the example
       above, the density matrices passed are (:math:`\boldsymbol{D}`,
       :math:`\boldsymbol{D}^{b}`, :math:`\boldsymbol{D}^{c}`,
       :math:`\boldsymbol{D}^{bc}`) and the associated perturbation tuple
       ``dens_mat_tuple`` is (1, 3, 4, 7) because these numbers correspond to the
       positions of the ":math:`(k,n)`-surviving" perturbation patterns in the
       canonically ordered list.
   :type dens_mat_tuple: QInt\*
   :param num_dmat: number of collected AO based density matrices for the given
       perturbation tuple ``dens_mat_tuple`` and all frequency configurations,
       that is ``num_freq_configs``
       :math:`\times\prod_{\text{perturbation pattern}}N_{\text{perturbation pattern}}`,
       where :math:`N_{\text{perturbation pattern}}` is the number of density
       matrices per perturbation pattern for a frequency configuration
   :type num_dmat: QInt
   :param dens_mat: the collected AO based density matrices, size is ``num_dmat``,
       and ordered as (``density matrices for freq. config. #1``,
       ``density matrices for freq. config. #2``, ``...``)
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, equals to the product of
       ``num_freq_configs`` and the dimension of perturbation pattern specified
       by ``num_pert`` and ``pert_labels``
   :type num_int: QInt
   :param val_int: the integral matrices to be returned, size is ``num_int``,
       and ordered as (``perturbation pattern``, ``num_freq_configs``)
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_xc_fun_exp(num_pert,         \
                             pert_labels,      \
                             num_freq_configs, \
                             len_dmat_tuple,   \
                             dens_mat_tuple,   \
                             num_dmat,         \
                             dens_mat,         \
                             user_ctx,         \
                             num_exp,          \
                             val_exp)

   Callback function for getting expectation values of XC functional,
   the last argument for the function :py:meth:`OpenRSPAddXCFun`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_freq_configs: for the perturbation pattern specified by ``num_pert``
       and ``pert_labels``, the number of different frequency configurations to
       be considered
   :type num_freq_configs: QInt
   :param len_dmat_tuple: for the perturbation pattern specified by ``num_pert``
       and ``pert_labels``, the length (number of different perturbation patterns)
       of the AO based density matrices passed
   :type len_dmat_tuple: QInt
   :param dens_mat_tuple: the perturbation tuple of the AO based density matrices
       passed, as a canonically ordered list of all relevant perturbation patterns
       of the density matrices, size is ``len_dmat_tuple``
   :type dens_mat_tuple: QInt\*
   :param num_dmat: number of collected AO based density matrices for the given
       perturbation tuple ``dens_mat_tuple`` and all frequency configurations,
       that is ``num_freq_configs``
       :math:`\times\prod_{\text{perturbation pattern}}N_{\text{perturbation pattern}}`,
       where :math:`N_{\text{perturbation pattern}}` is the number of density
       matrices per perturbation pattern for a frequency configuration
   :type num_dmat: QInt
   :param dens_mat: the collected AO based density matrices, size is ``num_dmat``,
       and ordered as (``density matrices for freq. config. #1``,
       ``density matrices for freq. config. #2``, ``...``)
   :type dens_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of the expectation values, equals to the product of
       ``num_freq_configs`` and the dimension of perturbation pattern specified
       by ``num_pert`` and ``pert_labels``
   :type num_exp: QInt
   :param val_exp: the expectation values to be returned, size is ``num_exp``,
       and ordered as (``perturbation pattern``, ``num_freq_configs``)
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_nuc_contrib(num_pert,    \
                              pert_labels, \
                              user_ctx,    \
                              dim_pert,    \
                              val_nuc)

   Callback function for getting the nuclear contributions, the last argument
   for the function :py:meth:`OpenRSPAddNucContributions`.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param dim_pert: dimension of the perturbation pattern specified by
       ``pert_labels``
   :type dim_pert: QInt
   :param val_nuc: the nuclear contributions, size is ``dim_pert``
   :type val_nuc: QReal\*
   :rtype: QVoid

.. function:: get_linear_rsp_solution(size_pert,     \
                                      num_freq_sums, \
                                      freq_sums,     \
                                      RHS_mat,       \
                                      user_ctx,      \
                                      rsp_param)

   Callback function for the linear response equation solver, the last argument
   for the function :py:meth:`OpenRSPSetLinearRSPSolver`.

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
       ``num_freq_sums``, and ordered as (``size_pert``, ``num_freq_sums``)
   :type RHS_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param rsp_param: solved response parameters, size is ``size_pert`` :math:`\times`
       ``num_freq_sums``, and ordered as (``size_pert``, ``num_freq_sums``)
   :type rsp_param: QcMat\*[]
   :rtype: QVoid

.. .. function:: get_rsp_eigen_solution(num_excit, \
                                        eigen_val, \
                                        user_ctx,  \
                                        eigen_vec)
 
    Callback function for the response eigenvalue equation solver, the last argument
    for the function :py:meth:`OpenRSPSetRSPEigenSolver`.
 
    :param num_excit: number of excitations to be solved
    :type num_excit: QInt
    :param eigen_val: solved excitation energies, size is ``num_excit``
    :type eigen_val: QReal\*
    :param user_ctx: user-defined callback function context
    :type user_ctx: QVoid\*
    :param eigen_vec: eigenvectors solved from the eigenvalue problem,
        size is ``num_excit``
    :type eigen_vec: QcMat\*[]
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
