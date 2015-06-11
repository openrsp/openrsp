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

We name a *perturbation tuple* as a list of perturbations specified
by their labels (``pert_labels``), where the same perturbations will
be always consecutive.

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
   the second last argument for function ``OpenRSPSetPDBS``.

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
   :param num_int: number of the integral matrices
   :type num_int: QInt
   :param val_int: the integral matrices
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_overlap_exp(bra_num_pert,    \
                              bra_pert_labels, \
                              ket_num_pert,    \
                              ket_pert_labels, \
                              num_pert,        \
                              pert_labels,     \
                              num_dens,        \
                              ao_dens,         \
                              user_ctx,        \
                              num_exp,         \
                              val_exp)

   Callback function for getting expectation values of overlap integrals,
   the last argument for function ``OpenRSPSetPDBS``.

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
   :param num_dens: number of atomic orbital (AO) based density matrices
   :type num_dens: QInt
   :param ao_dens: the AO based density matrices
   :type ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values
   :type num_exp: QInt
   :param val_exp: the expectation values
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_one_oper_mat(num_pert,    \
                               pert_labels, \
                               user_ctx,    \
                               num_int,     \
                               val_int)

   Callback function for getting integral matrices of a one-electron operator,
   the second last argument for function ``OpenRSPAddOneOper``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices
   :type num_int: QInt
   :param val_int: the integral matrices
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_one_oper_exp(num_pert,    \
                               pert_labels, \
                               num_dens,    \
                               ao_dens,     \
                               user_ctx,    \
                               num_exp,     \
                               val_exp)

   Callback function for getting expectation values of a one-electron operator,
   the last argument for function ``OpenRSPAddOneOper``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_dens: number of atomic orbital (AO) based density matrices
   :type num_dens: QInt
   :param ao_dens: the AO based density matrices
   :type ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values
   :type num_exp: QInt
   :param val_exp: the expectation values
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_two_oper_mat(num_pert,     \
                               pert_labels,  \
                               num_var_dens, \
                               var_ao_dens,  \
                               user_ctx,     \
                               num_int,      \
                               val_int)

   Callback function for getting integral matrices of a two-electron operator,
   the second last argument for function ``OpenRSPAddTwoOper``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_var_dens: number of variable AO based density matrices
   :type num_var_dens: QInt
   :param var_ao_dens: the variable AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\boldsymbol{G}(\boldsymbol{D})`
   :type var_ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices
   :type num_int: QInt
   :param val_int: the integral matrices
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_two_oper_exp(num_pert,       \
                               pert_labels,    \
                               num_var_dens,   \
                               var_ao_dens,    \
                               num_contr_dens, \
                               contr_ao_dens,  \
                               user_ctx,       \
                               num_exp,        \
                               val_exp)

   Callback function for getting expectation values of a two-electron operator,
   the last argument for function ``OpenRSPAddTwoOper``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_var_dens: number of variable AO based density matrices
   :type num_var_dens: QInt
   :param var_ao_dens: the variable AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\boldsymbol{G}(\boldsymbol{D})`
   :type var_ao_dens: QcMat\*[]
   :param num_contr_dens: number of contracted AO based density matrices
   :type num_contr_dens: QInt
   :param contr_ao_dens: the contracted AO based density matrices (:math:`\boldsymbol{D}`)
       for calculating :math:`\mathrm{Tr}[\boldsymbol{G}\boldsymbol{D}]`
   :type contr_ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of expectation values
   :type num_exp: QInt
   :param val_exp: the expectation values
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_xc_fun_mat(num_pert,             \
                             pert_labels,          \
                             num_freq_configs,     \
                             num_dmat_per_tuple,   \
                             dmat_perts_one_tuple, \
                             num_dens,             \
                             ao_dens,              \
                             user_ctx,             \
                             num_int,              \
                             val_int)

   Callback function for getting integral matrices of exchange-correlation (XC)
   functional, the second last argument for function ``OpenRSPAddXCFun``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_freq_configs: for the perturbation tuple specified above, the number of
       different frequency configurations to be considered
   :type num_freq_configs: QInt
   :param num_dmat_per_tuple: for the perturbation tuple specified above, the number of
       different perturbation patterns in the density matrices passed
   :type num_dmat_per_tuple: QInt
   :param dmat_perts_one_tuple: specify the perturbation pattern in a density matrix
       collection for one frequency tuple as each relevant perturbation tuple's number
       in a canonically ordered listing of all perturbation tuple subsets, size is
       ``num_dmat_per_tuple``
   :type dmat_perts_one_tuple: QInt\*
   :param num_dens: number of collected density matrices, equals to
       ``num_freq_configs`` :math:`\times\prod`
       ``number of density matrices for each perturbation pattern for one frequency configuration``
   :type num_dens: QInt
   :param ao_dens: the collected density matrices, ordered as
       ``array for freq. config. #1``, ``array for freq. config. #2``, ...
   :type ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_int: number of the integral matrices, equal to the product of
       ``num_freq_configs`` and the size of perturbation tuple specified by
       ``pert_labels``
   :type num_int: QInt
   :param val_int: the integral matrices to be returned, ordered as
       ``array for freq. config. #1``, ``array for freq. config. #2``, ...
   :type val_int: QcMat\*[]
   :rtype: QVoid

.. function:: get_xc_fun_exp(num_pert,             \
                             pert_labels,          \
                             num_freq_configs,     \
                             num_dmat_per_tuple,   \
                             dmat_perts_one_tuple, \
                             num_dens,             \
                             ao_dens,              \
                             user_ctx,             \
                             num_exp,              \
                             val_exp)

   Callback function for getting expectation values of a two-electron operator,
   the last argument for function ``OpenRSPAddXCFun``.

   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param pert_labels: label for each perturbation, size is ``num_pert``
   :type pert_labels: QInt\*
   :param num_freq_configs: for the perturbation tuple specified above, the number of
       different frequency configurations to be considered
   :type num_freq_configs: QInt
   :param num_dmat_per_tuple: for the perturbation tuple specified above, the number of
       different perturbation patterns in the density matrices passed
   :type num_dmat_per_tuple: QInt
   :param dmat_perts_one_tuple: specify the perturbation pattern in a density matrix
       collection for one frequency tuple as each relevant perturbation tuple's number
       in a canonically ordered listing of all perturbation tuple subsets, size is
       ``num_dmat_per_tuple``
   :type dmat_perts_one_tuple: QInt\*
   :param num_dens: number of collected density matrices, equals to
       ``num_freq_configs`` :math:`\times\prod`
       ``number of density matrices for each perturbation pattern for one frequency configuration``
   :type num_dens: QInt
   :param ao_dens: the collected density matrices, ordered as
       ``array for freq. config. #1``, ``array for freq. config. #2``, ...
   :type ao_dens: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param num_exp: number of the expectation values, equal to the product of
       ``num_freq_configs`` and the size of perturbation tuple specified by
       ``pert_labels``
   :type num_exp: QInt
   :param val_exp: the expectation values to be returned, ordered as
       ``array for freq. config. #1``, ``array for freq. config. #2``, ...
   :type val_exp: QReal\*
   :rtype: QVoid

.. function:: get_linear_rsp_solution(num_freq_sums, \
                                      freq_sums,     \
                                      size_pert,     \
                                      RHS_mat,       \
                                      user_ctx,      \
                                      rsp_param)

   Callback function for the linear response equation solver, the last argument
   for function ``OpenRSPSetLinearRSPSolver``.

   :param num_freq_sums: number of frequency sums on the left hand side
   :type num_freq_sums: QInt
   :param freq_sums: the frequency sums on the left hand side
   :type freq_sums: QReal\*
   :param size_pert: size of perturbaed matrices
   :type size_pert: QInt
   :param RHS_mat: RHS matrices, size is ``num_freq_sums``:math:`\times`
       ``size_pert``
   :type RHS_mat: QcMat\*[]
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param rsp_param: solved response parameters, size is ``num_freq_sums``:math:`\times`
       ``size_pert``
   :type rsp_param: QcMat\*[]
   :rtype: QVoid

.. function:: get_rsp_eigen_solution()

   Callback function for the response eigenvalue equation solver, the last argument
   for function ``OpenRSPSetRSPEigenSolver``.

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
