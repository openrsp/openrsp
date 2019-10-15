.. _subsection_api_reference:

API Reference
-------------

In order to use OpenRSP, C users should first include the header file
of OpenRSP in their codes::

  #inclde "OpenRSP.h"

while Fortran users should use the OpenRSP module::

  use OpenRSP_f

In this chapter, we will describe all the functions defined in OpenRSP
API for users. These functions should be invoked as::

  ierr = OpenRSP...(...)

where ``ierr`` contains the error information. Users should check if
it equals to ``QSUCCESS`` (constant defined in
`QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_). If not, there
was error happened in the invoked function, and the calculations should
stop.

Functions of OpenRSP API (C version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: QErrorCode OpenRSPCreate(open_rsp, num_atoms)

   Creates the context of response theory calculations, should be called at first.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\* (struct\*)
   :param num_atoms: number of atoms (**to be removed after perturbation free scheme implemented**)
   :type num_atoms: const QInt
   :rtype: QErrorCode (error information)

.. .. c:function:: QErrorCode OpenRSPSetWaveFunction(open_rsp, elec_wav_type)
.. 
..    Sets the (electronic) wave function.
.. 
..    :var open_rsp: context of response theory calculations
..    :vartype open_rsp: OpenRSP\*
..    :param elec_wav_type: the type of (electronic) wave function
..    :type elec_wav_type: ElecWavType (enum)
..    :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPSetLinearRSPSolver(open_rsp, user_ctx, get_linear_rsp_solution)

   Sets the context of linear response equation solver.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_linear_rsp_solution: user-specified callback function of linear
       response equation solver, see the callback function
       :c:func:`get_linear_rsp_solution`
   :type get_linear_rsp_solution: const GetLinearRSPSolution (function
       pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPSetPerturbations(open_rsp, num_pert_lab, pert_labels, pert_max_orders, pert_num_comps, user_ctx, get_pert_concatenation)

   Sets all perturbations involved in response theory calculations.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels involved
       in calculations
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param pert_num_comps: number of components of a perturbation described by
       exactly one of the above different labels, up to the allowed maximal
       order, size is therefore the sum of ``pert_max_orders``
   :type pert_num_comps: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_pert_concatenation: user specified function for getting the ranks
       of components of sub-perturbation tuples (with the same perturbation
       label) for given components of the corresponding concatenated
       perturbation tuple
   :type get_pert_concatenation: const GetPertCat (function pointer void (\*)(...))
   :rtype: QErrorCode

**NOTE**: :c:func:`get_pert_concatenation` will not be invoked in the current
release; OpenRSP will use it after the perturbation free scheme implmented.

.. c:function:: QErrorCode OpenRSPSetOverlap(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_overlap_mat, get_overlap_exp)

   Sets the overlap operator.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels that can
       act on the overlap operator
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_overlap_mat: user-specified callback function to calculate
       integral matrices of overlap operator as well as its derivatives with
       respect to different perturbations, see the callback function
       :c:func:`get_overlap_mat`
   :type get_overlap_mat: const GetOverlapMat (function pointer void (\*)(...))
   :param get_overlap_exp: user-specified callback function to calculate
       expectation values of overlap operator as well as its derivatives with
       respect to different perturbations, see the callback function
       :c:func:`get_overlap_exp`
   :type get_overlap_exp: const GetOverlapExp (function pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddOneOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_one_oper_mat, get_one_oper_exp)

   Adds a one-electron operator to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels that can
       act on the one-electron operator
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_one_oper_mat: user-specified callback function to calculate
       integral matrices of one-electron operator as well as its derivatives
       with respect to different perturbations, see the callback function
       :c:func:`get_one_oper_mat`
   :type get_one_oper_mat: const GetOneOperMat (function pointer void (\*)(...))
   :param get_one_oper_exp: user-specified callback function to calculate
       expectation values of one-electron operator as well as its derivatives
       with respect to different perturbations, see the callback function
       :c:func:`get_one_oper_exp`
   :type get_one_oper_exp: const GetOneOperExp (function pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddTwoOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_two_oper_mat, get_two_oper_exp)

   Adds a two-electron operator to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels that can
       act on the two-electron operator
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_two_oper_mat: user-specified callback function to calculate
       integral matrices of two-electron operator as well as its derivatives
       with respect to different perturbations, see the callback function
       :c:func:`get_two_oper_mat`
   :type get_two_oper_mat: const GetTwoOperMat (function pointer void (\*)(...))
   :param get_two_oper_exp: user-specified callback function to calculate
       expectation values of two-electron operator as well as its derivatives
       with respect to different perturbations, see the callback function
       :c:func:`get_two_oper_exp`
   :type get_two_oper_exp: const GetTwoOperExp (function pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddXCFun(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_xc_fun_mat, get_xc_fun_exp)

   Adds an exchange-correlation (XC) functional to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels that can
       act on the XC functional
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_xc_fun_mat: user-specified callback function to calculate
       integral matrices of XC functional as well as its derivatives with
       respect to different perturbations, see the callback function
       :c:func:`get_xc_fun_mat`
   :type get_xc_fun_mat: const GetXCFunMat (function pointer void (\*)(...))
   :param get_xc_fun_exp: user-specified callback function to calculate
       expectation values of XC functional as well as its derivatives with
       respect to different perturbations, see the callback function
       :c:func:`get_xc_fun_exp`
   :type get_xc_fun_exp: const GetXCFunExp (function pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddZeroOper(open_rsp, num_pert_lab, pert_labels, pert_max_orders, user_ctx, get_zero_oper_contrib)

   Adds a zero-electron operator to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert_lab: number of all *different* perturbation labels that can
       act on the zero-electron operator
   :type num_pert_lab: const QInt
   :param pert_labels: all the *different* perturbation labels involved
   :type pert_labels: const QcPertInt\*
   :param pert_max_orders: allowed maximal order of a perturbation described by
       exactly one of the above different labels
   :type pert_max_orders: const QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: void\*
   :param get_zero_oper_contrib: user-specified function to calculate
       contributions from the zero-electron operator, see the callback function
       :c:func:`get_zero_oper_contrib`
   :type get_zero_oper_contrib: const GetZeroOperContrib (function pointer void (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAssemble(open_rsp)

   Assembles the context of response theory calculations and checks its validity,
   should be called before any function ``OpenRSPGet...()``, otherwise the results
   might be incorrect.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPWrite(open_rsp, fp_rsp)

   Writes the context of response theory calculations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: const OpenRSP\*
   :param fp_rsp: file pointer
   :type fp_rsp: FILE\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPGetRSPFun(open_rsp, ref_ham, ref_state, ref_overlap, num_props, len_tuple, pert_tuple, num_freq_configs, pert_freqs, kn_rules, r_flag, write_threshold, size_rsp_funs, rsp_funs)

   Gets the response functions for given perturbations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param ref_ham: Hamiltonian of referenced state
   :type ref_ham: const QcMat\*
   :param ref_state: electronic state of referenced state
   :type ref_state: const QcMat\*
   :param ref_overlap: overlap integral matrix of referenced state
   :type ref_overlap: const QcMat\*
   :param num_props: number of properties to calculate
   :type num_props: const QInt
   :param len_tuple: length of perturbation tuple for each property,
       size is the number of properties (``num_props``)
   :type len_tuple: const QInt\*
   :param pert_tuple: ordered list of perturbation labels (perturbation
       tuple) for each property, size is ``sum(len_tuple)``, the first
       label of each property is the perturbation :math:`a`
   :type pert_tuple: const QcPertInt\*
   :param num_freq_configs: number of different frequency configurations
       for each property, size is ``num_props``
   :type num_freq_configs: const QInt\*
   :param pert_freqs: complex frequencies of each perturbation label (except
       for the perturbation :math:`a`) over all frequency configurations, size is
       ``2`` :math:`\times`
       ``(dot_product(len_tuple,num_freq_configs)-sum(num_freq_configs))``, and
       arranged as ``[num_freq_configs[i]][len_tuple[i]-1][2]`` (``i`` runs from
       ``0`` to ``num_props-1``) and the real and imaginary parts of each frequency
       are consecutive in memory
   :type pert_freqs: const QReal\*
   :param kn_rules: number :math:`k` for the :math:`(k,n)` rule [#]_ for each
       property (OpenRSP will determine the number :math:`n`), size is the
       number of properties (``num_props``)
   :type kn_rules: const QInt\*
   :param r_flag: flag to determine the restarting setup; `0` means "do not
       load/use any existing restarting data and do not save any new restarting
       data", and `3` means "use any existing restarting data and extend existing
       restarting data with all new restarting data"
   :type r_flag: const QInt
   :param write_threshold: tensor elements with absolute value below
       ``write_threshold`` will not be output by OpenRSP
   :type write_threshold: const QReal
   :param size_rsp_funs: size of the response functions, equals to the sum of
       the size of each property to calculate---which is the product of the
       size of added perturbations (specified by the perturbation tuple
       ``pert_tuple``) and the number of frequency configurations
       ``num_freq_configs`` for each property
   :type size_rsp_funs: const QInt
   :var rsp_funs: the response functions, size is ``2`` :math:`\times`
       ``size_rsp_funs`` and arranged as
       ``[num_props][num_freq_configs][pert_tuple][2]``,
       where the real and imaginary parts of the response functions
       are consecutive in memory
   :vartype rsp_funs: QReal\*
   :rtype: QErrorCode

.. [#] The description of the :math:`(k,n)` rule can be found, for instance,
       in [Ringholm2014]_.

.. c:function:: QErrorCode OpenRSPGetResidue(open_rsp, ref_ham, ref_state, ref_overlap, order_residue, num_excit, excit_energy, eigen_vector, num_props, len_tuple, pert_tuple, residue_num_pert, residue_idx_pert, num_freq_configs, pert_freqs, kn_rules, r_flag, write_threshold, size_residues, residues)

   Gets the residues for given perturbations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param ref_ham: Hamiltonian of referenced state
   :type ref_ham: const QcMat\*
   :param ref_state: electronic state of referenced state
   :type ref_state: const QcMat\*
   :param ref_overlap: overlap integral matrix of referenced state
   :type ref_overlap: const QcMat\*
   :param order_residue: order of residues, that is also the length of
       each excitation tuple
   :type order_residue: const QInt
   :param num_excit: number of excitation tuples that will be used for
       residue calculations
   :type num_excit: const QInt
   :param excit_energy: excitation energies of all tuples, size is
       ``order_residue`` :math:`\times` ``num_excit``, and arranged
       as ``[num_excit][order_residue]``; that is, there will be
       ``order_residue`` frequencies of perturbation labels (or sums
       of frequencies of perturbation labels) respectively equal to
       the ``order_residue`` excitation energies per tuple
       ``excit_energy[i][:]`` (``i`` runs from ``0`` to ``num_excit-1``)
   :type excit_energy: const QReal\*
   :param eigen_vector: eigenvectors (obtained from the generalized
       eigenvalue problem) of all excitation tuples, size is ``order_residue``
       :math:`\times` ``num_excit``, and also arranged in memory
       as ``[num_excit][order_residue]`` so that each eigenvector has
       its corresponding excitation energy in ``excit_energy``
   :type eigen_vector: QcMat\*[]
   :param num_props: number of properties to calculate
   :type num_props: const QInt
   :param len_tuple: length of perturbation tuple for each property,
       size is the number of properties (``num_props``)
   :type len_tuple: const QInt\*
   :param pert_tuple: ordered list of perturbation labels (perturbation
       tuple) for each property, size is ``sum(len_tuple)``, the first
       label of each property is the perturbation :math:`a`
   :type pert_tuple: const QcPertInt\*
   :param residue_num_pert: for each property and each excitation energy
       in the tuple, the number of perturbation labels whose sum of
       frequencies equals to that excitation energy, size is ``order_residue``
       :math:`\times` ``num_props``, and arragned as ``[num_props][order_residue]``;
       a negative ``residue_num_pert[i][j]`` (``i`` runs from ``0`` to
       ``num_props-1``) means that the sum of frequencies of perturbation
       labels equals to ``-excit_energy[:][j]``
   :type residue_num_pert: const QInt\*
   :param residue_idx_pert: for each property and each excitation energy
       in the tuple, the indices of perturbation labels whose sum of
       frequencies equals to that excitation energy, size is
       ``sum(residue_num_pert)``, and arranged as ``[residue_num_pert]``
   :type residue_idx_pert: const QInt\*
   :param num_freq_configs: number of different frequency configurations
       for each property, size is ``num_props``
   :type num_freq_configs: const QInt\*
   :param pert_freqs: complex frequencies of each perturbation label (except
       for the perturbation :math:`a`) over all frequency configurations and
       excitation tuples, size is ``2`` :math:`\times`
       ``(dot_product(len_tuple,num_freq_configs)-sum(num_freq_configs))``
       :math:`\times` ``num_excit``, and arranged as
       ``[num_excit][num_freq_configs[i]][len_tuple[i]-1][2]`` (``i`` runs from
       ``0`` to ``num_props-1``) and the real and imaginary parts of each
       frequency are consecutive in memory; notice that the (sums of)
       frequencies of perturbation labels specified by ``residue_idx_pert``
       should equal to the corresponding excitation energies for all frequency
       configurations and excitation tuples
   :type pert_freqs: const QReal\*
   :param kn_rules: number :math:`k` for the :math:`(k,n)` rule for each property
       (OpenRSP will determine the number :math:`n`), size is the number of
       properties (``num_props``)
   :type kn_rules: const QInt\*
   :param r_flag: flag to determine the restarting setup; `0` means "do not
       load/use any existing restarting data and do not save any new restarting
       data", and `3` means "use any existing restarting data and extend existing
       restarting data with all new restarting data"
   :type r_flag: const QInt
   :param write_threshold: tensor elements with absolute value below
       ``write_threshold`` will not be output by OpenRSP
   :type write_threshold: const QReal
   :param size_residues: size of the residues, equals to the sum of the
       size of each property to calculate---which is the product of the
       size of added perturbations (specified by the perturbation tuple
       ``pert_tuple``), the number of excitation tuples (``num_excit``)
       and the number of frequency configurations ``num_freq_configs``
       for each property
   :type size_residues: const QInt
   :var residues: the residues, size is ``2`` :math:`\times`
       ``size_residues`` and arranged as
       ``[num_props][num_excit][num_freq_configs][pert_tuple][2]``, where
       the real and imaginary parts of the residues are consecutive in memory
   :vartype residues: QReal\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPDestroy(open_rsp)

   Destroys the context of response theory calculations, should be called at the end.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. _subsubsection_fortran_convention:

Functions of OpenRSP API (Fortran version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions of OpenRSP API (Fortran) are similar to those of the C version, except
that an extra ``_f`` should be appended to each function. Other differences are
the (ii) argument types and (iii) callback functions (subroutines for Fortran).
The latter will be described in Chapter :ref:`subsection_callback_functions`. The
former relates to the convention of types in Fortran, please refer to the manual
of `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_ and the following
table for the convention:

.. list-table::
   :header-rows: 1

   * - Type in OpenRSP
     - Fortran
   * - ``struct OpenRSP``
     - ``type(OpenRSP)``
   * - ``void* user_ctx``
     - ``type(C_PTR) user_ctx``
   * - callback functions
     - external subroutines

We also want to mention that users can also pass their own defined Fortran type
as the user-defined callback function context to OpenRSP, by encapsulated into
the ``type(C_PTR) user_ctx``.
