.. _chapter-API-reference:

OpenRSP API Reference
=====================

In order to use OpenRSP, C users should first include the header file
of OpenRSP in their codes::

  #inclde "openrsp.h"

while Fortran users should use the OpenRSP module::

  use openrsp_f

In this chapter, we will describe all the functions defined in OpenRSP
API for users. These functions should be invoked as::

  ierr = OpenRSP...(...)

where ``ierr`` contains the error information. Users should check if
it equals to ``QSUCCESS`` (constant defined in
`QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_). If not, there
was error happened in the invoked function, and the calculations should
stop.

Functions of OpenRSP API (C version)
------------------------------------

.. c:function:: QErrorCode OpenRSPCreate(open_rsp)

   Creates the context of response theory calculations, should be called at first.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\* (struct\*)
   :rtype: QErrorCode (error information)

.. c:function:: QErrorCode OpenRSPSetElecEOM(open_rsp, elec_EOM_type)

   Sets the equation of motion (EOM) of electrons.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param elec_EOM_type: the type of EOM of electrons
   :type elec_EOM_type: ElecEOMType (enum)
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPSetLinearRSPSolver(open_rsp, user_ctx, get_linear_rsp_solution)

   Sets the context of linear response equation solver.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_linear_rsp_solution: user specified function of linear
       response equation solver, see the callback function
       :c:func:`get_linear_rsp_solution`
   :type get_linear_rsp_solution: GetLinearRSPSolution (function
       pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. Host programs will call OpenRSP by sending the excited states, so that we
   do not need the function OpenRSPSetRSPEigenSolver
.. .. c:function:: QErrorCode OpenRSPSetRSPEigenSolver(open_rsp, user_ctx, get_rsp_eigen_solution)
 
    Sets the context of response eigenvalue solver.
 
    :var open_rsp: context of response theory calculations
    :vartype open_rsp: OpenRSP\*
    :param user_ctx: user-defined callback function context
    :type user_ctx: QVoid\*
    :param get_rsp_eigen_solution: user specified function of response
        eigenvalue equation solver, see the callback function
        :c:func:`get_rsp_eigen_solution`
    :type get_rsp_eigen_solution: GetRSPEigenSolution (function
        pointer QVoid (\*)(...))
    :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPSetPerturbations(open_rsp, num_pert, pert_labels, pert_max_orders, pert_num_comps, user_ctx, get_pert_concatenation)

   Sets all perturbation labels involved in response theory calculations.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of all *different* perturbation labels involved
       in calculations
   :type num_pert: QInt
   :param pert_labels: all *different* perturbation labels involved
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param pert_num_comps: number of components of each perturbation (label) up to
       its maximum order, size is the sum of ``pert_max_orders``
   :type pert_num_comps: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_pert_concatenation: user specified function for getting the
       rank of concatenation of several perturbation tuples with the same
       perturbation label
   :type get_pert_concatenation: GetPertCat (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

*FIXME: get_pert_comp and get_pert_rank to be discussed and implemented*

.. c:function:: QErrorCode OpenRSPSetPDBS(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_overlap_mat, get_overlap_exp)

   Sets the context of perturbation dependent basis sets.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of *different* perturbation labels that can
       act as perturbations on the basis sets
   :type num_pert: QInt
   :param pert_labels: all the *different* perturbation labels
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_overlap_mat: user specified function for getting overlap
       integrals, see the callback function :c:func:`get_overlap_mat`
   :type get_overlap_mat: GetOverlapMat (function pointer QVoid (\*)(...))
   :param get_overlap_exp: user specified function for getting expectation
       values of overlap integrals, see the callback function
       :c:func:`get_overlap_exp`
   :type get_overlap_exp: GetOverlapExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddOneOper(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_one_oper_mat, get_one_oper_exp)

   Adds a one-electron operator to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of *different* perturbation labels that can
       act as perturbations on the one-electron operator
   :type num_pert: QInt
   :param pert_labels: all the *different* perturbation labels
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_one_oper_mat: user specified function for getting integral matrices,
       see the callback function :c:func:`get_one_oper_mat`
   :type get_one_oper_mat: GetOneOperMat (function pointer QVoid (\*)(...))
   :param get_one_oper_exp: user specified function for getting expectation values,
       see the callback function :c:func:`get_one_oper_exp`
   :type get_one_oper_exp: GetOneOperExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddTwoOper(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_two_oper_mat, get_two_oper_exp)

   Adds a two-electron operator to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of *different* perturbation labels that can
       act as perturbations on the two-electron operator
   :type num_pert: QInt
   :param pert_labels: all the *different* perturbation labels
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_two_oper_mat: user specified function for getting integral matrices,
       see the callback function :c:func:`get_two_oper_mat`
   :type get_two_oper_mat: GetTwoOperMat (function pointer QVoid (\*)(...))
   :param get_two_oper_exp: user specified function for getting expectation values,
       see the callback function :c:func:`get_two_oper_exp`
   :type get_two_oper_exp: GetTwoOperExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAddXCFun(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_xc_fun_mat, get_xc_fun_exp)

   Adds an exchange-correlation (XC) functional to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of *different* perturbation labels that can
       act as perturbations on the XC functional
   :type num_pert: QInt
   :param pert_labels: all the *different* perturbation labels
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_xc_fun_mat: user specified function for getting integral matrices,
       see the callback function :c:func:`get_xc_fun_mat`
   :type get_xc_fun_mat: GetXCFunMat (function pointer QVoid (\*)(...))
   :param get_xc_fun_exp: user specified function for getting expectation values,
       see the callback function :c:func:`get_xc_fun_exp`
   :type get_xc_fun_exp: GetXCFunExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPSetNucContributions(open_rsp, num_pert, pert_labels, pert_max_orders, user_ctx, get_nuc_contrib)

   Sets the nuclear contributions to the Hamiltonian.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :param num_pert: number of *different* perturbation labels that can
       act as perturbations on the nuclear Hamiltonian
   :type num_pert: QInt
   :param pert_labels: all the *different* perturbation labels
   :type pert_labels: QInt\*
   :param pert_max_orders: maximum allowed order of each perturbation (label)
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_nuc_contrib: user specified function for getting the nuclear
       contributions, see the callback function :c:func:`get_nuc_contrib`
   :type get_nuc_contrib: GetNucContrib (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. The following APIs do not need because the nuclear contributions will also
   be obtained through callback function from the host programs.
.. .. c:function:: QErrorCode OpenRSPSetNucGeoPerturbations(open_rsp, num_atoms, atom_coord, atom_charge)
   
      Sets the context of geometric perturbations for nuclear Hamiltonian.
   
      :var open_rsp: context of response theory calculations
      :vartype open_rsp: OpenRSP\*
      :param num_atoms: number of atoms
      :type num_atoms: QInt
      :param atom_coord: coordinates of atoms
      :type atom_coord: QReal\*
      :param atom_charge: charges of atoms
      :type atom_charge: QReal\*
      :rtype: QErrorCode

.. .. c:function:: QErrorCode OpenRSPSetNucScalarPotential(open_rsp, dipole_origin)
   
     Sets the terms in nuclear Hamiltonian due to the scalar potential.
  
     :var open_rsp: context of response theory calculations
     :vartype open_rsp: OpenRSP\*
     :param dipole_origin: coordinates of dipole origin
     :type dipole_origin: QReal[3]
     :rtype: QErrorCode

.. .. c:function:: OpenRSPSetNucVectorPotential(open_rsp, gauge_origin)
   
      Sets the terms in nuclear Hamiltonian due to the vector potential.
   
      :var open_rsp: context of response theory calculations
      :vartype open_rsp: OpenRSP\*
      :param gauge_origin: coordinates of gauge origin
      :type gauge_origin: QReal[3]
      :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPAssemble(open_rsp)

   Assembles the context of response theory calculations and checks its validity,
   should be called before any function ``OpenRSPGet...()``, otherwise the results
   might be incorrect.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPWrite(open_rsp, file_name)

   Writes the context of response theory calculations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param file_name: the name of the file
   :type file_name: QChar\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPGetRSPFun(open_rsp, ref_ham, ref_state, ref_overlap, num_props, len_tuple, pert_tuple, num_freq_configs, pert_freqs, kn_rules, size_rsp_funs, rsp_funs)

   Gets the response functions for given perturbations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param ref_ham: Hamiltonian of referenced state
   :type ref_ham: QcMat\*
   :param ref_state: electronic state of referenced state
   :type ref_state: QcMat\*
   :param ref_overlap: overlap integral matrix of referenced state
   :type ref_overlap: QcMat\*
   :param num_props: number of properties to calculate
   :type num_props: QInt
   :param len_tuple: length of perturbation tuple for each property,
       size is the number of properties (``num_props``)
   :type len_tuple: QInt\*
   :param pert_tuple: ordered list of perturbation labels (perturbation
       tuple) for each property, size is ``sum(len_tuple)``, the first
       label of each property is the perturbation :math:`a`
   :type pert_tuple: QInt\*
   :param num_freq_configs: number of different frequency configurations
       for each property, size is ``num_props``
   :type num_freq_configs: QInt\*
   :param pert_freqs: complex frequencies of each perturbation label (except
       for the perturbation :math:`a`) over all frequency configurations,
       size is ``2*(dot_product(len_tuple,num_freq_configs)-sum(num_freq_configs))``,
       and arranged as ``[num_freq_configs[i]][len_tuple[i]-1][2]`` (``i``
       runs from ``1`` to ``num_props``) and the real and imaginary parts
       of each frequency are consecutive in memory
   :type pert_freqs: QReal\*
   :param kn_rules: number :math:`k` for the :math:`kn` rule for each property
       (OpenRSP will determine the number :math:`n`), size is the number of
       properties (``num_props``)
   :type kn_rules: QInt\*
   :param size_rsp_funs: size of the response functions, equals to the sum of
       the size of each property to calculate---which is the product of the
       size of added perturbations (specified by the perturbation tuple
       ``pert_tuple``) and the number of frequency configurations
       ``num_freq_configs`` for each property
   :type size_rsp_funs: QInt
   :var rsp_funs: the response functions, size is ``2`` :math:`\times`
       ``size_rsp_funs`` and arranged as
       ``[num_props][num_freq_configs][pert_tuple][2]``,
       where the real and imaginary parts of the response functions
       are consecutive in memory
   :vartype rsp_funs: QReal\*
   :rtype: QErrorCode

.. c:function:: QErrorCode OpenRSPGetResidue(open_rsp, ref_ham, ref_state, ref_overlap, num_excit, excit_energy, eigen_vector, num_props, len_tuple, pert_tuple, order_residue, num_freq_configs, pert_freqs, kn_rules, size_residues, residues)

   Gets the residues for given perturbations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param ref_ham: Hamiltonian of referenced state
   :type ref_ham: QcMat\*
   :param ref_state: electronic state of referenced state
   :type ref_state: QcMat\*
   :param ref_overlap: overlap integral matrix of referenced state
   :type ref_overlap: QcMat\*
   :param num_excit: number of excitations
   :type num_excit: QInt
   :param excit_energy: excitation energies, size is ``num_excit``
   :type excit_energy: QReal\*
   :param eigen_vector: eigenvectors obtained from the generalized
       eigenvalue problem, size is ``num_excit``
   :type eigen_vector: QcMat\*[]
   :param num_props: number of properties to calculate
   :type num_props: QInt
   :param len_tuple: length of perturbation tuple for each property,
       size is the number of properties (``num_props``)
   :type len_tuple: QInt\*
   :param pert_tuple: ordered list of perturbation labels (perturbation
       tuple) for each property, size is ``sum(len_tuple)``, the first
       label of each property is the perturbation :math:`a`
   :type pert_tuple: QInt\*
   :param order_residue: order of residues
   :type order_residue: QInt
   :param num_freq_configs: number of different frequency configurations
       for each property, size is ``num_props``
   :type num_freq_configs: QInt\*
   :param pert_freqs: complex frequencies of each perturbation label (except
       for the perturbation :math:`a`) over all frequency configurations,
       size is ``2*(dot_product(len_tuple,num_freq_configs)-sum(num_freq_configs))``,
       and arranged as ``[num_freq_configs[i]][len_tuple[i]-1][2]`` (``i``
       runs from ``1`` to ``num_props``) and the real and imaginary parts
       of each frequency are consecutive in memory
   :type pert_freqs: QReal\*
   :param kn_rules: number :math:`k` for the :math:`kn` rule for each property
       (OpenRSP will determine the number :math:`n`), size is the number of
       properties (``num_props``)
   :type kn_rules: QInt\*
   :param size_residues: size of the residues, equals to the sum of the
       size of each property to calculate---which is the product of the
       size of added perturbations (specified by the perturbation tuple
       ``pert_tuple``) and the number of frequency configurations
       ``num_freq_configs`` for each property
   :type size_residues: QInt
   :var residues: the residues, size is ``2`` :math:`\times`
       ``size_residues`` and arranged as
       ``[num_props][num_freq_configs][pert_tuple][2]``, where the real
       and imaginary parts of the residues are consecutive in memory
   :vartype residues: QReal\*
   :rtype: QErrorCode

*FIXME:*

#. Which perturbations to which excited states, +/-excitation energy?
#. Will calculating several different properties save time?
#. Are the size_residues and residues OK?

.. c:function:: QErrorCode OpenRSPDestroy(open_rsp)

   Destroys the context of response theory calculations, should be called at the end.

   :var open_rsp: context of response theory calculations
   :vartype open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. _section-Fortran-convention:

Functions of OpenRSP API (Fortran version)
------------------------------------------

Functions of OpenRSP API (Fortran) are similar to those of the C version, except
that an extra ``_f`` should be appended to each function. Other differences are
the (ii) argument types and (iii) callback functions (subroutines for Fortran).
The latter will be described in Chapter :ref:`chapter-callback-functions`. The
former relates to the convention of types in Fortran, please refer to the manual
of `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_ and the following
table for the convention:

.. list-table::
   :header-rows: 1

   * - Type in OpenRSP
     - Fortran
   * - ``struct OpenRSP``
     - ``type(OpenRSP)``
   * - ``enum ElecEOMType``
     - ``integer``
   * - ``QVoid* user_ctx``
     - ``character(len=1) user_ctx(:)``
   * - callback functions
     - external subroutines

We also want to mention that users can also pass their own defined Fortran type
as the user-defined callback function context to OpenRSP (thanks to the Fortran
function ``transfer``). For instance, the following code transfers the ``type(QcMat)``
variable ``A`` to a character array ``enc``::

  type(QcMat) A
  character(len=1), allocatable :: enc(:)
  integer len_enc
  len_enc = size(transfer(A, enc))
  allocate(enc(len_enc))
  enc = transfer(A, enc)

Users could then send ``enc`` to OpenRSP, and which will be passed to callback
functions later on, and could be decoded (in the callback functions) as::

  integer, intent(in) :: len_ctx
  character(len=1), intent(in) :: user_ctx(len_ctx)
  ... ...
  type(QcMat) A
  A = transfer(enc, A)
