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
`QMatrix library <http://repo.ctcc.no/projects/qmatrix>`_). If not, there
was error happened in the invoked function, and the calculations should
stop.

Functions of OpenRSP API (C version)
------------------------------------

.. function:: OpenRSPCreate(open_rsp)

   Creates the context of response theory calculations, should be called at first.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\* (struct\*)
   :rtype: QErrorCode (error information)

.. function:: OpenRSPSetElecEOM(open_rsp, elec_EOM_type)

   Sets the equation of motion (EOM) of electrons.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param elec_EOM_type: the type of EOM of electrons
   :type elec_EOM_type: ElecEOMType (enum)
   :rtype: QErrorCode

.. function:: OpenRSPSetSolver(open_rsp, \
                               user_ctx, \
                               get_rsp_solution)

   Sets the context of response equation solver.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_rsp_solution: user specified function of response equation solver
   :type get_rsp_solution: GetRSPSolution (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. function:: OpenRSPSetPerturbations(open_rsp,        \
                                      num_pert,        \
                                      perturbations,   \
                                      pert_max_orders, \
                                      pert_sizes,      \
                                      user_ctx,        \
                                      get_pert_comp,   \
                                      get_pert_rank)

   Sets all perturbations involved in response theory calculations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param num_pert: number of all perturbations involved in calculations
   :type num_pert: QInt
   :param perturbations: all perturbations involved in calculations
   :type perturbations: QInt\*
   :param pert_max_orders: maximum allowed orders of all perturbations
   :type pert_max_orders: QInt\*
   :param pert_sizes: sizes of all perturbations up to their maximum orders,
       whose dimension is the sum of ``pert_max_orders``
   :type pert_sizes: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_pert_comp: user specified function for getting components of a perturbation
   :type get_pert_comp: GetPertComp (function pointer QVoid (\*)(...))
   :param get_pert_rank: user specified function for getting rank of a perturbation
   :type get_pert_rank: GetPertRank (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. function:: OpenRSPSetPDBS(open_rsp,        \
                             num_pert,        \
                             perturbations,   \
                             pert_max_orders, \
                             user_ctx,        \
                             get_overlap_mat, \
                             get_overlap_exp)

   Sets the context of perturbation dependent basis sets.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param num_pert: number of perturbations that the basis sets depend on
   :type num_pert: QInt
   :param perturbations: perturbations that the basis set depend on
   :type perturbations: QInt\*
   :param pert_max_orders: maximum allowed orders of the perturbations
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_overlap_mat: user specified function for getting overlap integrals
   :type get_overlap_mat: GetOverlapMat (function pointer QVoid (\*)(...))
   :param get_overlap_exp: user specified function for getting expectation values of overlap integrals
   :type get_overlap_exp: GetOverlapExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. function:: OpenRSPAddOneOper(open_rsp,         \
                                num_pert,         \
                                perturbations,    \
                                pert_max_orders,  \
                                user_ctx,         \
                                get_one_oper_mat, \
                                get_one_oper_exp)

   Adds a one-electron operator to the Hamiltonian.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param num_pert: number of perturbations that the one-electron operator depends on
   :type num_pert: QInt
   :param perturbations: perturbations that the one-electron operator depends on
   :type perturbations: QInt\*
   :param pert_max_orders: maximum allowed orders of the perturbations
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_one_oper_mat: user specified function for getting integral matrices
   :type get_one_oper_mat: GetOneOperMat (function pointer QVoid (\*)(...))
   :param get_one_oper_exp: user specified function for getting expectation values
   :type get_one_oper_exp: GetOneOperExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. function:: OpenRSPAddTwoOper(open_rsp,         \
                                num_pert,         \
                                perturbations,    \
                                pert_max_orders,  \
                                user_ctx,         \
                                get_two_oper_mat, \
                                get_two_oper_exp)

   Adds a two-electron operator to the Hamiltonian.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param num_pert: number of perturbations that the two-electron operator depends on
   :type num_pert: QInt
   :param perturbations: perturbations that the two-electron operator depends on
   :type perturbations: QInt\*
   :param pert_max_orders: maximum allowed orders of the perturbations
   :type pert_max_orders: QInt\*
   :param user_ctx: user-defined callback function context
   :type user_ctx: QVoid\*
   :param get_two_oper_mat: user specified function for getting integral matrices
   :type get_two_oper_mat: GetTwoOperMat (function pointer QVoid (\*)(...))
   :param get_two_oper_exp: user specified function for getting expectation values
   :type get_two_oper_exp: GetTwoOperExp (function pointer QVoid (\*)(...))
   :rtype: QErrorCode

.. function:: OpenRSPSetAtoms(open_rsp,   \
                              num_atoms,  \
                              atom_coord, \
                              atom_charge)

   Sets the context of atoms for the nuclear contributions.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param num_atoms: number of atoms
   :type num_atoms: QInt
   :param atom_coord: coordinates of atoms
   :type atom_coord: QReal\*
   :param atom_charge: charges of atoms
   :type atom_charge: QReal\*
   :rtype: QErrorCode

.. function:: OpenRSPSetDipoleOrigin(open_rsp, \
                                     dipole_origin)

   Sets the coordinates of dipole origin.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param dipole_origin: coordinates of dipole origin
   :type dipole_origin: QReal[3]
   :rtype: QErrorCode

.. function:: OpenRSPSetGaugeOrigin(open_rsp, \
                                    gauge_origin)

   Sets the coordinates of gauge origin.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param gauge_origin: coordinates of gauge origin
   :type gauge_origin: QReal[3]
   :rtype: QErrorCode

.. function:: OpenRSPAssemble(open_rsp)

   Assembles the context of response theory calculations and checks its validity,
   should be called before any function ``OpenRSPGet...()``, otherwise the results
   might be incorrect.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. function:: OpenRSPWrite(open_rsp, file_name)

   Writes the context of response theory calculations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param file_name: the name of the file
   :type file_name: QChar\*
   :rtype: QErrorCode

.. function:: OpenRSPGetRSPFun(open_rsp,      \
                               ref_ham,       \
                               ref_state,     \
                               ref_overlap,   \
                               num_pert,      \
                               perturbations, \
                               pert_orders,   \
                               pert_freqs,    \
                               kn_rule,       \
                               size_rsp_fun,  \
                               rsp_fun)

   Gets the response function for given perturbations.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :param ref_ham: Hamiltonian of referenced state
   :type ref_ham: QMat\*
   :param ref_state: electronic state of referenced state
   :type ref_state: QMat\*
   :param ref_overlap: overlap integral matrix of referenced state
   :type ref_overlap: QMat\*
   :param num_pert: number of perturbations
   :type num_pert: QInt
   :param perturbations: the perturbations
   :type perturbations: QInt\*
   :param pert_orders: orders of the perturbations
   :type pert_orders: QInt\*
   :param pert_freqs: frequencies of the perturbations
   :type pert_freqs: QReal\*
   :param kn_rule: contains the perturbation :math:`a` and numbers :math:`k` and :math:`n`
   :type kn_rule: QInt[3]
   :param size_rsp_fun: size of the response function, equals to the
       product of sizes of ``perturbations``
   :type size_rsp_fun: QInt
   :param rsp_fun: the response function
   :type rsp_fun: QReal\*
   :rtype: QErrorCode

.. function:: OpenRSPDestroy(open_rsp)

   Destroys the context of response theory calculations, should be called at the end.

   :param open_rsp: context of response theory calculations
   :type open_rsp: OpenRSP\*
   :rtype: QErrorCode

.. _section-Fortran-convention:

Functions of OpenRSP API (Fortran version)
------------------------------------------

Functions of OpenRSP API (Fortran) are similar to those of the C version, except
that an extra ``_f`` should be appended to each function. Other differences are
the (ii) argument types and (iii) callback functions (subroutines for Fortran).
The latter will be described in Chapter :ref:`chapter-callback-functions`. The
former relates to the convention of types in Fortran, please refer to the manual
of `QMatrix library <http://repo.ctcc.no/projects/qmatrix>`_ and the following
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
function ``transfer``). For instance, the following code transfers the ``type(QMat)``
variable ``A`` to a character array ``enc``::

  type(QMat) A
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
  type(QMat) A
  A = transfer(enc, A)
