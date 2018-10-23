

===========================
Daniel's code documentation
===========================

This document is intended as a rough documentation of the fundamental
procedures of the OpenRSP program suite. Its purposes are both to keep
an overview for the experienced programmer and to be an introduction for
the newcomer building some kind of a bridge between the code and the
corresponding articles. It can also serve as some kind of to-do-list for
present of future or present implementations but not yet implemented
parts of the could should be marked properly. Text which is set in
italic within the pseudocode sections gives some explanation of what the
code does at the particular place.

Please note that this documentation is for internal use only.

The subroutine ``openrsp_calc``
===============================

Treatment of several different input cases:

Magnetizability
---------------

VCD
---

Manual specification
--------------------

Manual input of a set of perturbation of frequencies:

-  Specify :math:`k` and :math:`n`

-  allocate ``perturbation_tuple``

-  loop 1, :math:`N` (:math:`N` number of perturbations):

   -  Characterize each perturbation as MAGnetic, ELectric field or
      GEOmetrical distortion

-  loop 1, :math:`M` (:math:`M` number of frequencies)

   -  Initialize frequencies on tuple

   -  if (:math:`k=1`\ ) then

      -  call ``rsp_prop`` with ``F_unpert``=F, ``D_unpert``=D,
         ``S_unpert``=S

      else

      -  ``F_already`` :math:`\rightarrow` F

      -  ``D_already`` :math:`\rightarrow` D

      -  ``S_already`` :math:`\rightarrow` S

      -  call ``rsp_prop`` with ``F_already``=``F_already``,
         ``D_already``=``D_already``, ``S_already``=``S_already``

Excitations
-----------

Calling the interface to the response solver for the determination of
eigenvalues and eigenvectors.

Residues
--------

Modified version of the manual specification algorithm to calculate the
residues. Several differences have been implemented:

-  The only mandatory input for residues is the calculation of
   excitation energies (beforehead) and the residue order (e.g. the
   number of photons in the corresponding transition property).

-  Operator types, frequencies and states for which the residues are to
   be calculated can be given optionally

-  If the optioms are not set, the residues of the corresponding order
   will be calculated for all states using the electric field as
   operator. All frequencies are by default set to the :math:`n`\ th
   part of the excitation energy of the corresponding state if :math:`n`
   is the order of the residue.

-  For calculating the residue of order :math:`n` an input for the
   response function of order :math:`2\cdot n` is formed which is worked
   through with some modifications.

-  Specify :math:`k` and :math:`n`

-  Set default for frequencies

-  allocate ``perturbation_tuple``

-  :math:`N` = :math:`2\cdot n` (:math:`N` number of perturbations,
   :math:`n` order of the residue)

-  loop 1, :math:`N`\ :

   -  Characterize each perturbation as MAGnetic, ELectric field or
      GEOmetrical distortion and set the corresponding dimension of the
      perturbation either to 3 oder :math:`3\cdot \mathcal N`
      (:math:`\mathcal N` number of atoms)

-  loop 1, :math:`M` (:math:`M` number of frequencies)

   -  Initialize frequencies on tuple

   -  Initialize labels on tuple, if needed set to default

   -  if (:math:`k=1`\ ) then

      -  call ``rsp_prop`` with ``F_unpert``=F, ``D_unpert``=D,
         ``S_unpert``=S

      else

      -  ``F_already`` :math:`\rightarrow` F

      -  ``D_already`` :math:`\rightarrow` D

      -  ``S_already`` :math:`\rightarrow` S

      -  call ``rsp_prop`` with ``F_already``=``F_already``,
         ``D_already``=``D_already``, ``S_already``=``S_already``

Second harmonic generation
--------------------------

PV2F
----

PV3F
----

PV4F
----

Hyper Raman
-----------

efishgcid
---------

???

The subroutine ``rsp_prop``
===========================

Calllist:
``  subroutine rsp_prop(pert_unordered, kn, F_unpert, D_unpert, S_unpert,``
``F_already, D_already, S_already, zeromat_already, file_id, Xf_already, Df_already)``,
last 10 objects optional

-  if present ``Xf_already`` or ``Df_already``, ``do_residues``
   :math:`\rightarrow` true

-  if .not. present ``S_already`` then

   -  ``F_already`` :math:`\rightarrow` F

   -  ``D_already`` :math:`\rightarrow` D

   -  ``S_already`` :math:`\rightarrow` S

   -  initialize ``zeromat``

-  else

   -  initialize ``zeromat_already``

-  Determine some arrays concerning the number of perturbations:

   -  ``num_blks`` :math:`\rightarrow` No. of block (’’races’’) of
      perturbations depending on freq. and operator

   -  ``blk_info`` :math:`\rightarrow` Array containing the start index,
      the number of perturbation and the perturbation numbers on the
      block

   -  ``blk_sizes`` :math:`\rightarrow` Number of non-redundant
      perturbations on the block (taking into account also the
      perturbation dimension

   -  ``property_size`` :math:`\rightarrow` Number of non-redundant
      elements of the result tensor

-  if present ``F_already`` then

   -  call ``get_prop`` with ``F_already``,``D_already``, ``S_already``

   else

   -  call ``get_prop`` with F, D, S and also with Xf and Df if needed

The subroutine ``get_prop``
===========================

Calllist:
``pert, kn, nr_ao, num_blks, blk_sizes, blk_info, property_size, prop, F, D, S,``
``do_residues, Xf, Df``

In this routine the different contributions to the seeked property are
compiled.

-  call ``rsp_fds`` for calculating perturbed F, D and S

-  call ``rsp_energy`` for calculating the energy contributions

-  call ``rsp_xcave_interface`` for calculating the exchange-correlation
   contributions

-  call ``rsp_pulay_kn`` for calculating Pulay :math:`k,n`\ -type
   contributions

-  call ``rsp_pulay_lap`` for calculating Pulay Laplace-contributions

-  call ``rsp_idem_lag`` for calculating indempotency Lagrangian
   contributions

-  call ``rsp_scfe_lag`` for calculating SCF Lagrangian contributions

Algorithm 1: The subroutine ``rsp_fds``
=======================================

Calllist: ``zeromat, pert, kn, F, D, S``

This subroutine is recursive and coordinates the calculation of the
perturbed S, F and D intermediates. It corresponds to **Algorithm 1**
from the Paper on the recursive OpenRSP scheme.

-  if ``pert%n_perturbations`` .gt. 1 then

   -  Make ``pert%n_perturbations`` subsets of ``pert``, with
      ``n_perturbations`` reduced by 1 for each element of the subset.

   -  do :math:`i` = 1, ``pert%n_perturbations``rsp\_energy

      -  if ``sdf_already(D,psub)`` .eqv. false then: (see below for a
         description of this condition)

         -  call ``rsp_fds`` with the :math:`i`\ th element of the
            perturbation subset

   -  if ``sdf_already(D,psub)``.eqv.false then: (see below for a
      description of this condition)

      -  if ``kn_skip(...)``.eqv.false then: (see below for a
         description of this condition)

         -  :math:`k` = 1

         -  do :math:`j`\ =1, ``pert%n_perturbations``

            -  ``pert%pid``(:math:`j`\ ) = :math:`k`

            -  :math:`k = k + 1`

         -  call ``get_fds``

``sdf_already(D,psub)``
-----------------------

This logical function checks whether the current perturbed quantities
have already been calculated. The perturbed quantities S, D and F are
linked to each other by pointers showing which of them is the next in
the perturbation row and whether one is the last in the row.

``kn_skip(...)``
----------------

This function checks whether the current perturbed quantity is needed
due to the global indices :math:`k` and :math:`n`\ .

Algorithm 2: The subroutine ``get_fds``
=======================================

This subroutine coordinates the calculation of the perturbed
intermediates S, F and D and therefore it is also the interface to the
Response solver which is needed for the calculation of D. It corresponds
to what is named **Algorithm 2** in the paper on the recursive OpenRSP
scheme.

-  Determine the size of the current perturbation

-  Calculate perturbed S and store it on the appropriate variable

-  Calculate initial part of perturbed F (call ``rsp_fock_lowerorder``)

-  Determine the superstructure of D (what’s that?)

-  loop over size of superstructure:

   -  call ``rsp_get_matrix_z`` to contstruct :math:`\mathbf{D}_p`

   -  :math:`\mathbf{D}_p = \mathbf{D}_p - \mathbf{D} \cdot \mathbf{S} \cdot \mathbf{D}_p - \mathbf{D}_p \cdot \mathbf{S} \cdot \mathbf{D}`

   -  Calculate two-electron contribution to Fp

   -  Calculate xc-contribution to Fp

   -  Calculate pe-contribution to Fp

   -  Initialize RHS and Xf

   -  Call ``rsp_get_matrix_y`` to calculate RHS

   -  if frequency sum of current perturbation not equal to excitation
      energy then

      -  Solve the linear equation system: call ``rsp_solver_exec``

      -  Calculate ``Dh``

      -  Calculate homogeneous contribution to the Fock matrix, added
         onto Fp

      -  Add ``Dh`` onto ``Dp``

   -  else

      -  Read ``Xf``

      -  contract ``Xf`` and ``RHS``

      -  Read ``Df`` and store on ``Dp``

      -  Scale ``Dp`` with contraction, calculate ``Df`` and write on
         linked list instead of full ``D``

      -  Nullify ``Fp``

      -  Calculate homogeneous contribution to the Fock matrix, added
         onto Fp

   -  end if

   -  add ``Fp`` and ``Dp`` on the corresponding linked lists

Algorithm 3-type routines
=========================

These routines follow **Algorithm 3** from the paper on the recursive
OpenRSP scheme and are therefore recursive. In all Algorithm 3-type
subroutines changes due to the calculation of residues will have to be
made.

The subroutine ``rsp_energy``
-----------------------------

This subroutine in general calculates terms which depend on contractions
of :math:`\boldsymbol{\mathcal E}` and its derivatives with perturbed
and non-pertubed densities.

Coordinates the calculation of the energy contributions.

Calllist: ``rsp_energy(pert, total_num_perturbations, kn,``

``num_p_tuples,p_tuples, density_order, D, property_size, cache, prop)``

``pert,p_tuples`` are of the ``p_tuple`` derived type. ``pert`` is the
original perturbation tuple and is used most for keeping track of the
recursion levels. ``p_tuple`` is the ’’working tuple’’ of this
subroutine. It is an array with ``num_p_tuples`` elements in the current
call of ``rsp_energy``. ``density_order`` is a measurement for the
number of derivatives which are found in the densities.
``total_num_perturbations`` is the value ``n_perturbations`` from the
initial ``pert``-tuple. This value is kept unchanged from the first call
through the whole recursion and is handed over to ``get_energy`` at the
end.

``p_tuples`` works as follows: The first element of the array always
lists the perturbations which are set on the integrals. The second
element lists the perturbations on the first density, the third element
lists the perturbations on the second density etc. ``density_order`` is
the sum of the number of elements on ``p_tuples`` minus the number of
elements on the first element on ``p_tuples``.

-  if ``pert%n_perturbations.ge.1`` then

   -  if ``p_tuples(i)%n_perturbations``.eq.0 then

      -  call ``rsp_energy`` with the first perturbation removed from
         ``pert`` to be ``pert`` and ``p_tuples`` being the first
         element of ``pert``, both prepared by special functions

   -  else

      -  call ``rsp_energy`` with the first perturbation removed from
         ``pert`` to be ``pert`` and ``p_tuples`` being t``p_tuples``
         extended by one, both prepared by special functions

   -  removing one element from ``pert`` corresponds to go one level up
      in recursion. The removed element is the ’’current perturbation’’
      of the actual recursion level.

   -  The two self-calls above correspond to a setting the current
      perturbation on the integrals.

   -  do ``i=``2, ``num_p_tuples``

      -  if ``p_tuples(i)%n_perturbations``.eq.0 ``t_new(i)``
         :math:`\rightarrow` first element of ``pert``

      -  else ``t_new(i)`` ``t_new(i)`` :math:`\rightarrow` ``t_new(i)``
         extended by 1

      -  call ``rsp_energy`` with the first perturbation removed from
         ``pert`` to be ``pert`` and ``p_tuples`` being ``t_new(i)``

      -  This self-call corresponds to setting the current perturbation
         on one of the existing densities

   -  if ``num_p_tuples``.le.3 then

      -  call ``rsp_energy`` with the first perturbation removed from
         ``pert`` to be ``pert`` and ``p_tuples`` being a combination of
         ``p_tuples`` and the first element of ``pert``

      -  This self-call corresponds to a complete chain-rule-like
         derivative setting up a new singly-derived density carrying the
         current perturbation.

-  else

   -  This is the final recursion level

   -  Check whether the contribution is needed due to the :math:`k`\ ,
      :math:`n`\ -values and whether it is a relevant contribution to a
      residue (if needed)

   -  Get the data from cache (call ``property_cache_getdata``

   -  call ``get_energy``

The subroutine ``rsp_fock_lowerorder``
--------------------------------------

This subroutine in general calculates terms which depend on perturbed
derivates of :math:`\boldsymbol{\mathcal F}`\ .

Calllist:
``zeromat, pert, total_num_perturbations,num_p_tuples, p_tuples,``
``density_order, D, property_size, Fp,fock_lowerorder_cache``

This subroutine coordinates the calculation of perturbed Fock matrices
which do not depend on the homogeneous part of the perturbed density. It
is recursive and the first part is similar to the first part of
``rsp_energy``. Therefore most elements on the calllist are simliar.
``fock_lowerorder_cache`` is a special derived type variable to keep
track of the calculated Fock matrix intermediates.

-  if ``pert%n_perturbations.ge.1`` then

   -  Do the same as in the ``pert%n_perturbations.ge.1``-part of
      ``rsp_energy`` (recursive setup of perturbation lists)

-  else

   -  This is the final recursion level

   -  Determine which contributions can be skipped due to the density:
      All terms which contain densities that have the same or a higher
      perturbation order than ``total_num_perturbations`` or which are
      not relevant for a residue calculation (if needed) are sorted out.

   -  Function ``f_l_cache_already``: Check whether the corresponding
      element is found on cache

      -  if not found then

         -  call ``get_fock_lowerorder`` :math:`\rightarrow` calculate
            :math:`\mathbf{D}_P`\ , see JCP 128 2008 214108

         -  modify ``fock_lowerorder_cache`` (done in
            ``get_fock_lowerorder``)

      -  else

         -  Get the corresponding contribution from cache (call
            ``f_l_cache_getdata``

The subroutine ``rsp_pulay_kn``
-------------------------------

This subroutine calculates terms of the :math:`(\mathbf{SW})^k_{n_W}`\ -type.
The recursion scheme in this subroutine is similar to the one in
``rsp_energy`` nevertheless it is less complex since the perturbations
can only be set either on :math:`\mathbf{S}` and :math:`\mathbf{W}`\ .
Therefore there are just two self-calls.

Callist: ``pert, kn, p12, S, D, F, property_size, cache, prop``

The perturbation-tuple type variables ``pert`` and ``p12`` play about
the same role as ``pert`` and ``p_tuples`` in the subroutines above.
``p12`` is a one-d-array with two elements on it.

-  if ``pert%n_perturbations`` .gt. 0 then

   -  call ``rsp_pulay_kn`` with first element removed from ``pert``
      (’’current element’’) and the first element of ``p12`` extended by
      the current element. The second element of ``p12`` is left as it
      is.

   -  Corresponds to putting the current perturbation onto
      :math:`\mathbf{S}`\ .

   -  call ``rsp_pulay_kn`` with first element removed from ``pert``
      (’’current element’’) and the first element of ``p12`` left as it
      is. The second element of ``p12`` is extended by the curreny
      element.

   -  Corresponds to putting the current perturbation onto
      :math:`\mathbf{W}`\ .

-  else

   -  Final recursion level

   -  Look up whether the corresponding contribution is needed due to
      the :math:`k`\ , :math:`n` parameters. :math:`\rightarrow`\ **to-do:
      Identify relevant terms for residue calculation here**

   -  if needed then

      -  call ``get_pulay_kn``

The subroutine ``rsp_pulay_lag``
--------------------------------

Very similar to ``rsp_pulay_kn``. Calculates the
:math:`(\mathbf{S}^a \mathbf{W})^{bc...})_{k_s,n^{'}_Y}`\ -type terms

Calllist: ``pert, kn, p12, S, D, F, property_size, cache, prop``

-  if ``pert%n_perturbations`` .gt. 0 then

   -  Do the same recursion scheme as in ``rsp_pulay_kn``

-  else

   -  Final recursion level

   -  Check whether the term is needed due to the :math:`k` and
      :math:`n`\ -parameters :math:`\rightarrow` **to-do: Check here
      whether term is needed for residues or not**

   -  Check whether the contribution is already in cache

   -  if needed and not on cache then

      -  call ``get_pulay_lag``

The subroutine ``rsp_scfe_lag``
-------------------------------

Calculates the terms of the
:math:`(\lambda^a \mathbf{Y})^{bc...}_{k_{\lambda},n^{'}_{Y}}`\ -type

Calllist: ``pert, kn, p12, S, D, F, property_size, cache, prop``, very
similar to the two routines before

-  if ``pert%n_perturbations`` .gt. 0 then

   -  Do the same recursion scheme as in ``rsp_pulay_kn``, first
      self-call: Perturbation on :math:`\lambda`\ , second self-call:
      Perturbation on :math:`\mathbf{Y}`

-  else

   -  Final recursion level

   -  Check whether the term is needed due to the :math:`k` and
      :math:`n`\ -parameters :math:`\rightarrow` **to-do: Check here
      whether term is needed for residues or not**

   -  Check whether the contribution is already in cache

   -  if needed and not on cache then

      -  call ``get_scfe_lag``

The subroutine ``rsp_idem_lag``
-------------------------------

Calculates the terms of the
:math:`(\zeta^a \mathbf{Z})^{bc...}_{k_{\zeta},n^{'}_{Z}}`\ -type

Calllist: ``pert, kn, p12, S, D, F, property_size, cache, prop``, very
similar to the three routines before

-  if ``pert%n_perturbations`` .gt. 0 then

   -  Do the same recursion scheme as in ``rsp_pulay_kn``, first
      self-call: Perturbation on :math:`\zeta`\ , second self-call:
      Perturbation on :math:`\mathbf{Z}`

-  else

   -  Final recursion level

   -  Check whether the term is needed due to the :math:`k` and
      :math:`n`\ -parameters :math:`\rightarrow` **to-do: Check here
      whether term is needed for residues or not**

   -  Check whether the contribution is already in cache

   -  if needed and not on cache then

      -  call ``get_idem_lag``

The subroutines coordinating the calculation of perturbed contributions
=======================================================================

These routines are called at the final recursion level of the Algorith
3-type routines.

The subroutine ``get_energy``
-----------------------------

Calllist:
``num_p_tuples, total_num_perturbations,p_tuples, density_order, D,``

``property_size, cache, prop``

with ``p_tuples`` being of the ``p_tuple`` data type

-  Assemble elements of ``p_tuples`` in blocks

-  Determine the proper indices of ``p_tuples``

-  Read unperturbed and perturbed densities

-  Calculate the contributions from

   -  perturbed one-electron integrals

   -  perturbed overlap matrices

   -  perturbed two-electron integrals

The subroutine ``get_fock_lowerorder``
--------------------------------------

This routine coordinates the calculation of the contributions to the
Fock matrix-type intermediayes.

Calllist: ``zeromat, num_p_tuples, total_num_perturbations, p_tuples,``
``density_order, D, property_size, Fp,fock_lowerorder_cache``

``p_tuple`` is of the perturbation tuple type (array with dimension
``num_p_tuples``) while ``fock_lowerorder_cache`` is a linked list
variable-type which is used for caching the intermediates. ``zeromat``
is an empty matrix-type variable while ``D`` is a linked list variable
containing the perturbed densities.

-  clone 1st element of ``p_tuples`` to ``t_matrix_newpid`` (for use
   only if all perturbations are on the integrals)

-  Characterize the elements of ``p_tuples`` .w.r.t. blocks etc.

-  if ``total_num_perturbations.gt.p_tuples(1)%n_perturbations`` then

   -  Loop over size of outer indeces

      -  Read perturbed densities in a loop over ``num_p_tuples``
         starting with the 2nd element

      -  ``if num_p_tuples .le. 1`` Calculate one-electron integral
         contributions, results added on ``tmp``

      -  ``if num_p_tuples .le. 2`` Calculate the-electron integral
         contributions, results added on ``tmp``

      -  ``if num_p_tuples .le. 2`` Calculate PE-contributions, results
         added on ``tmp``

      -  Calculate xc-contributions, results added on ``tmp``
         (``num_p_tuples`` is the number of derivatives w.r.t. the
         densities. 1el-contributions become zero at derivatives higher
         than 1, 2el-contractions become zero at derivatives higher than
         2, xc-contributions do not necessarily vanish.

      -  ``if num_p_tuples(1)%n_perturbations .gt. 0`` then

         -  Loop over size of inner indices

            -  Determine offsets

            -  Write ``tmp`` on ``lower_order_contribution``

      -  else

         -  Initialize ``lower_order_contribution`` with 0

      -  end if

   -  Merge all elements of ``p_tuples`` to one merged tuple

   -  Put merged tuple in standard order

   -  Collect all perturbation indices on from merged tuple on one
      1d-array

   -  Characterize the merged tuple w.r.t. blocks etc.

   -  Determine offsets and add ``lower_order_contribution`` onto ``Fp``

-  else

   -  ``if num_p_tuples .le. 1`` Calculate one-electron integral
      contributions, results added on ``Fp``

   -  ``if num_p_tuples .le. 2`` Calculate the-electron integral
      contributions, results added on ``Fp``

   -  ``if num_p_tuples .le. 2`` Calculate PE-contributions, results
      added on ``Fp``

   -  Calculate xc-contributions, results added on ``Fp``
      (``num_p_tuples`` is the number of derivatives w.r.t. the
      densities. 1el-contributions become zero at derivatives higher
      than 1, 2el-contractions become zero at derivatives higher than 2,
      xc-contributions do not necessarily vanish.

-  end if

-  Nullify a large amount of variables.

The subroutine ``rsp_get_matrix_zeta``
--------------------------------------

This subroutine is responsible for the calculation of the :math:`\zeta`\ -Lagrangian
multipliers.

Calllist: ``zeromat, p_tuple_a, kn, superstructure_size, deriv_struct,``
``total_num_perturbations, which_index_is_pid, indices_len,``
``ind, F, D, S, Zeta``

The variables ``p_tuple_a`` and ``deriv_struct`` are of the tuple-type.
The result is being returned on ``Zeta`` (matrix type).

This subroutine reads several sets of D, F and S intermediates and
composes them to the perturbed :math:`\zeta`\ .

The following steps are made:

-  loop over size of the superstructure

   -  For residues: Determine whether :math:`\mathbf{D}` is appropriate.
      Otherwise: Skip contributions.

      -  :math:`\zeta + \mathcal F \mathbf{D} \mathbf{S}`

         -  Reading F with a merge of ``p_tuple_a`` and 1st element of
            ith column of ``current_derivative_term``

         -  Reading D with the 2nd element of ith column of
            ``current_derivative_term``

         -  Reading S with the 3rd element of ith column of
            ``current_derivative_term``

      -  :math:`\zeta- \mathcal F \mathbf{D} \mathbf{S}`

         -  Reading F with the 1st element of ith column of
            ``current_derivative_term``

         -  Keep D

         -  Reading S with a merge of ``p_tuple_a`` and 3rd element of
            ith column of ``current_derivative_term``

      -  :math:`\zeta + \tfrac{1}{2} \omega \mathcal S \mathbf{D} \mathbf{S}`

         -  Reading 1st S with the 1st element of ith column of
            ``current_derivative_term``

         -  Keep D

         -  Keep 2nd S

      -  :math:`\zeta +  \mathbf{S} \mathbf{D} \mathcal F`

         -  Reading S with 1st element of ith column of
            ``current_derivative_term``

         -  Reading D with the 2nd element of ith column of
            ``current_derivative_term``

         -  Reading F with a merge of ``p_tuple_a`` and 3rd element of
            ith column of ``current_derivative_term``

      -  :math:`\zeta- \mathbf{S} \mathbf{D} \mathcal F`

         -  Reading S with a merge of ``p_tuple_a`` and 1st element of
            ith column of ``current_derivative_term``

         -  Keep D

         -  Reading F with the 3rd element of ith column of
            ``current_derivative_term``

      -  :math:`\zeta - \tfrac{1}{2} \omega \mathcal S \mathbf{D} \mathbf{S}`

         -  Keep 1st S

         -  Keep D

         -  Reading 2nd S with the 3rd element of the ith column of
            ``current_derivative_term``

-  Add contribution of a fully perturbed :math:`\boldsymbol{\mathcal F}`
   if allowed due to :math:`k` and :math:`n`\ .

The subroutine ``get_scfe_lag``
-------------------------------

This subroutine coordinates the calculation of the Lagrangian
multipliers :math:`\lambda^a` and and the matrix :math:`\mathbf{Y}`\ .

Calllist: ``p12, kn, F, D, S, property_size, cache, prop``, with ``p12``
being a 1d-array of the perturbation-tuple type containing the
perturbations for :math:`\lambda^a` and :math:`\mathbf{Y}` separately.
``property_size`` is an integer, ``cache`` is of the cache-data type and
``prop`` is a real 1d-array of ``property_size`` containing the
contractions made in this routine.

-  Setup derivative superstructures for :math:`\lambda` and
   :math:`\mathbf{Y}`\ .

-  Characterize the elements of ``p12`` w.r.t. blocks etc.

-  loop over number of outer indices

   -  Check whether current perturbation combination is relevant for
      residues. Two possibilities: 1.) Residue relevant perturbations do
      not contribute to :math:`\lambda\quad \rightarrow` whole lambda is
      calculated, term selection in :math:`\mathbf{Y}`\ . 2.) Residue
      relevant perturbations contribute to
      :math:`\lambda\quad \rightarrow` term selection in :math:`\lambda`\ ,
      whole :math:`\mathbf{Y}` is calculated.

   -  call ``rsp_get_matrix_lambda``

      -  call ``rsp_get_matrix_y``

      -  store contraction of result on ``prop_forcache``

-  merge elements of ``p12`` together, put to standardorder and
   characterize

-  determine offset

-  add elements of ``prop_forcache`` onto ``prop``

The subroutine ``rsp_get_matrix_lambda``
----------------------------------------

Calculates transformations of the
:math:`\mathbf{D}_1 \cdot \mathbf{S} \cdot \mathbf{D}_2 - \mathbf{D}_2 \cdot \mathbf{S} \cdot \mathbf{D}_1`\ -type

Calllist: ``zeromat, p_tuple_a, superstructure_size, deriv_struct,``
``total_num_perturbations, which_index_is_pid, indices_len, ind,``
``select_terms,D, S, L``

``p_tuple_a`` and ``deriv_struct`` are of the tuple-type. ``p_tuple_a``
corresponds to the 1st element of ``p12`` in the calling routine.
``deriv_struct`` is an outcome of a superstructure determination.
``select_terms`` triggers the selection of terms in residue
calculations. ``L`` is a matrix-type variable for :math:`\lambda`\ .

-  loop over superstructure size

   -  merge ``p_tuple_a`` and ``deriv_struct(i,1)`` :math:`\rightarrow`
      ``merged_a``

   -  merge ``p_tuple_a`` and ``deriv_struct(i,3)`` :math:`\rightarrow`
      ``merged_b``

   -  read :math:`\mathbf{D_1}` with ``merged_a``, store on A

   -  read :math:`\mathbf{S}` with ``deriv_struct(i,2)``, store on B

   -  read :math:`\mathbf{D_2}` with ``deriv_struct(i,3)``, store on C

   -  :math:`L = L - A \cdot B \cdot C`\ , check before whether one of
      the densities fits with the residue condition

   -  read :math:`\mathbf{D_1}` with ``deriv_struct(i,1)``, store on A

   -  read :math:`\mathbf{D_2}` with ``merged_b``, store on C

   -  :math:`L = L - A \cdot B \cdot C`\ , check before whether one of
      the densities fits with the residue conditions

The subroutine ``rsp_get_matrix_y``
-----------------------------------

Is used for the calculation of :math:`\mathbf{Y}` intermediates both for
the TD-SCF-part and for the calculation of the right-hand-side for the
response equations.

Callist: ``zeromat, superstructure_size, deriv_struct,``
``total_num_perturbations, which_index_is_pid, indices_len,``
``ind, select_terms, F, D, S, Y``

``deriv_struct`` is of the tuple-type and an outcome of a superstructure
determination. ``select_terms`` triggers the selection of terms in
residue calculations. ``Y`` is a matrix-type variable for
:math:`\mathbf{Y}`\ .

-  Loop over superstructure size

   -  read F with 1st element of ``deriv_struct``, store on A

   -  read D with 2nd element of ``deriv_struct``, store on B

   -  read S with 3rd element of ``deriv_struct``, store on C

   -  :math:`Y = Y + A \cdot B \cdot C`\ , check whether D or F fit with
      residues if needed

   -  read S with 1st element of ``deriv_struct``, store on A

   -  :math:`Y = Y + \omega A \cdot B \cdot C`\ , :math:`\omega` is a
      frequency sum, check whether D fits with residues if needed

   -  read S with 1st element of ``deriv_struct``, store on A

   -  read F with 3rd element of ``deriv_struct``, store on C

   -  :math:`Y = Y - A \cdot B \cdot C`\ , check whether D or F fit with
      residues if needed

   -  read S with 1st element of ``deriv_struct``, store on A

   -  read S with 3rd element of ``deriv_struct``, store on C

   -  :math:`Y = Y - \tfrac{1}{2} \cdot \omega A \cdot B \cdot C`\ ,
      slightly different shape of :math:`\omega` due to different
      combinations of :math:`\omega` on intermediates, check whether D
      fits with residues if needed

The exchange-correlation part
=============================

The subroutine ``rsp_xcave_interface``
--------------------------------------

Coordinates the calculation of the xc-contributions to the energy
derivatives. Handles energy orders up to order 5 which can be composed
of different compositions of electrical field and geometrical
distortions. Pure geometrical distortions need no treatment of perturbed
cx-contributions. The main point for this subroutine is to read to
proper densities and to call the subroutine ``xc_integrate`` which does
the integral-density contractions. ``xc_integrate`` is called in a
system of :math:`n+m` loops where :math:`n` is the number of geometrical
distortions and :math:`m` is the number of electrical field
perturbations which are involved. Every loop runs over :math:`3\cdot N`
steps where :math:`N` is the number of atoms or over three steps where
every step represents one spatial component of the electrical field.

subroutine ``xc_integrate``
---------------------------

-  Store density matrix on file ``dmat``

-  do some initializations

-  Determine level of derivative w.r.t. geometrical distortion

-  read functional which is to be used

-  Read number of batches

-  do loop over number of batches

   -  Read number of grid points

   -  do some allocations and read grid

   -  redefine grid points (call ``xcint_mpi_distribute_points``)

   -  call ``xc_integrate_batch``

   -  symmetrize result

This subroutine sets up a set of batches and calls
``xc_integrate_batch`` in a loop over the number of batches.

subroutine ``xc_integrate_batch``
---------------------------------

This subroutine is up to now only documented for the case that the
``get_ave``-flag is true

-  Parse the functional and determine number of varibles according to
   the functional

-  loop points=1, nr\_points in steps of max\_block\_length

   -  nr\_points has to do with the grid

   -  max\_block\_length is set to 100 by parameter

   -  read AOs (call ``karaoke_get_ao``) and compress them (call
      ``karaoke_compress_ao``)

   -  evaluate the density (call ``evaluate_density``), i.e. do the
      contraction AO:math:`_{bk}\cdot` AO:math:`_{bl}\cdot`
      D:math:`_{kl}` assuming a symmetric density

   -  update number of evaluated electrons

   -  if ``get_ave`` then

      -  Determine the level of derivative w.r.t. geometrical distorions

      -  Determine the start indices of the perturbed density matrices
         which are needed (corresponding to the level of perturbation
         and the choice of :math:`k` and :math:`n` (call
         ``up_dmat_index`` in a set of if-conditions depending on the
         perturbation level).

      -  call ``evaluate_density`` once more, assuming a non-symmetric
         density this time

      -  Construct ``xcin`` from the output of ``evaluate_density``

      -  call ``xc_eval_star`` using ``xcin``, yielding ``xcout`` as
         output

      -  add ``xcout`` to ``energy`` termiwise

The subroutine ``xc_eval_star`` acts merely as some kind of wrapper
which calls ``xceval`` which is obiously written in C.

The handling of the perturbed matrices
======================================

The handling of the perturbed intermediates S, D and F in the recursive
scheme requires some effort and is done using two derived types and a
large amount of functions and subroutines managing them. These derived
types are:

The ``p_tuple`` derived type
----------------------------

This derived type variable contains the informations about the
perturbations that are handled. It contains informations on the type of
the perturbation which is handeled (electric or magnetic field encoded
by ``EL`` and ``MAG``, respectively; geometrical distortion denoted by
``GEO``), its level (first, second, third etc. derivative), its
dimension (number of components e.g. 3 for the first derivative
w.r.t. the electric field) and the corresponding frequencies.

**In order to implement residues** also the handling of excitation
energies and eigenvectors should be managed by ``p_tuple`` with ``EXCI``
as the ’’perturbation’’-type, the corresponding excitation energies as
frequencies and the number of excitations as dimension.

The ``matrix`` derived type
---------------------------

It contains of several quantities describing the size and the properties
(symmetry, hermiticity etc.) of the matrix, the elements and a
``self_pointer`` which is used for copying the matrix.

The ``SDF`` derived type
------------------------

This derived type contains of an allocatable array of
``matrix``-variables called ``data`` as well as a variable of the
``p_tuple``-type, a logical ``last`` and a pointer ``next`` which is of
the ``SDF``-type himself. So it can be imagined as some kind of a shell
for a matrix array which contains additional information concerning its
management.

In the corresponding initialization routine the ``matrix``-array is
allocated by 1 and the ``p_tuple``-variable is completely nullified.
``last`` is set to true indicating that this is the last intermediate on
the corresponding ``SDF``-type variable that has been formed. These
``last`` and ``next``-variables are used in the subroutine ``get_fds``
in the postprocessing of the calculated perturbed variables to keep
track of the calculated perturbed intermediates.

All through the program there are three major variables of this type:
``F``, ``D`` and ``S`` corresponding to the Fock matrix, the density
matrix and the overlap matrix and their pertubed derivatives. In general
the variables of the ``SDF`` derived type work like a linked list. This
means that the unperturbed matrix and its perturbed derivatives are
connected by pointers (``next`` in the ``SDF``-data type). Regarding
e.g. the calculation of a second hyperpolarizability (:math:`4^{th}`
derivative w.r.t. the electric field) the storage of all perturbed Fock
matrices ``F`` looks as follows:

ccccc

Element & Number of matrices & Content of ``p_tuple`` & Value of &
Element ``next``
No. & & & ``last`` & points to
1 & 1 & 0 (unperturbed) & false & 2
2 & 3 (drv. w.r.t. el. field comp. :math:`x,y,z`\ ) & 1:math:`^{st}`
drv. & false & 3
3 & 6 (drv. w.r.t. el. field comp. :math:`xx,xy,xz,yy,yz,xx`\ ) &
2:math:`^{nd}` drv. & false & 4
4 & 10 (drv. w.r.t. el. field comp. :math:`xxx,xxy,xxz,xyy,` &
3:math:`^{rd}` drv. & false & 5
& :math:`xyz,xzz,yyy,yyz,yzz,zzz`\ ) & &
5 & 15 (drv. w.r.t. el. field comp. :math:`xxxx,xxxy,` & 4:math:`^{th}`
drv. & true & 1
& :math:`xxxz,xxyy,xxyz,` & &
& :math:`xxzz,xyyy,xyyz,xyzz,xzzz,` & &
& :math:`yyyy,yyyz,yyzz,yzzz,zzzz` & &

The order of the elements in the table also corresponds to the order in
which these elements are calculated. The number of element in every
perturbation level is called the dimensions of the perturbation.

By following the pointer ``next`` from one perturbation level to the
next every perturbation level can be read by the corresponding routines
although there seems to be only one ``SDF``-type variable for e.g. the
Fock matrix.

**In order to implement the residues** the ``SDF`` derived type is to be
used to manage both the excitation eigenvectors :math:`\mathbf{X}_f` and
the corresponding excitation densities :math:`\mathbf{D}^f` with the two
new variables ``Xf`` and ``Df``.

Subroutines for handling the perturbed intermediates
----------------------------------------------------

All these subroutines and functions are collected in the
``rsp_sdf_caching.f90``.

The subroutine ``sdf_setup_datatype``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subroutine is called at the beginning of the ’’life’’ of a
``SDF``-type variable. It does the following:

-  Make ``next`` show on the present variable

-  Set ``last`` to ``true``

-  allocate ``p_tuple``-type array with zero

-  allocate the ``data``-array with 1 and initialize it with the matrix
   from the callist

This corresponds to setting up the variable which at this point only
contains an unperturbed intermediate.

The subroutine ``sdf_init``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subroutine starts the linked list by adding an input
perturbation-tuple to what is on the input ``SDF``-variable, allocating
the ``matrix``-part of the derived type and to set up the matrix-part of
the variable.

subroutine ``sdf_getdata_s``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Determines the offset of the seeked perturbation

-  Searches for the seeked level of perturbation by following the
   ``next``-pointer between the variables and comparing the inherent
   ``p_tuple``-type variables with the input one.

-  Writes the seeked data on the input variable

function ``sdf_next_element``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Moves the ``next``-pointer from the input array itself to the next one.

function ``sdf_getdata``
~~~~~~~~~~~~~~~~~~~~~~~~

Does about the same as the subroutinte ``sdf_getdata_s``

subroutine ``sdf_add``
~~~~~~~~~~~~~~~~~~~~~~

-  Checks whether the corresponding perturbed quantity has already been
   calculated

-  If not:

   -  Set up a variable ``new_element`` of the ``SDF`` data type

   -  Follow the pointer ``next`` until the last element (highest level
      in perturbation) has been found

   -  if found then

      -  Set ``last`` to ``false``

      -  Make ``next`` point to ``new_element``
      
      
      
Recognition of elements which are relevant for the residues
===========================================================

subroutine ``recognize_contribution``
-------------------------------------

In order to determine which term contributes to a residue and which 
not the recursive subroutine ``recognize_contribution`` examines the
frequencies of all relevant perturbation or all their possible sums
(depending on the input) on a match with the excitation energy or
two excitation energies, respectively. This subroutine works as follows:


  - Arguments Perturbation tuple (:math:`b_N`), number of elements in the sum (``n``), perturbation frequencies (:math:`omega_N`), number of actual sum elements (``j``), actual sum value (:math:`omega`)
  - Start values: number of actual sum elements=1,  actual sum value=0  
  - Result returned on logical ``recognized``
  - ``recognized`` set to false
  - for ``i`` in ``j``, ``n-1`` do
  
    - if ``n`` =1 then
    
       - if :math:`omega+omega_N(i)=` excitation energy then
       
          - ``recognized`` :math:`\leftarrow` true
          
          - return to previous invocation
          
       - end if
       
    - else
    
       - ``recognized`` :math:`\leftarrow` call self( :math:`b_N,n-1,omega+omega_N(i),i+1` )
       
    - end if
    
    - if ``recognized`` exit loops
    
  - end for

      
      

Blocking of perturbations and handling of the perturbation tuple
================================================================

subroutine ``derivative_superstructure``
----------------------------------------

This subroutine is recursive and Fill in the purpose - not yet properly
understood.

Calllist:
``pert, kn, primed,current_derivative_term, superstructure_size,``
``new_element_position, derivative_structure``

The variables ``pert,current_derivative_term`` and
``derivative_structure`` are of the ``p_tuple``-type. ``pert`` is a
scalar, ``current_derivative_term`` has the dimension 3 while
``derivative_structure`` is 2-dimensional with the dimensions
``superstructure_size`` and 3. ``primed`` is a logical.

-  if ``pert%n_perturbations.gt.0`` then

   -  1st self-call. Callist unchanged apart due to the original one
      apart from the first two tuple-type variables:

      #. ``pert`` reduced by the first element.

      #. 1st element of ``current_derivative_term`` extended by 1st
         element of ``pert``; elements 2 and 3 from
         ``current_derivative_term``

   -  2nd self-call. Callist unchanged apart due to the original one
      apart from the first two tuple-type variables:

      #. ``pert`` reduced by the first element.

      #. 1st element of ``current_derivative_term``; 2nd element of
         ``current_derivative_term`` extended by 1st element of
         ``pert``; 3rd element of ``current_derivative_term``

   -  3rd self-call. Callist unchanged apart due to the original one
      apart from the first two tuple-type variables:

      #. ``pert`` reduced by the first element.

      #. 1st and 2nd element of ``current_derivative_term``; 3rd element
         of ``current_derivative_term`` extended by 1st element of
         ``pert``

-  else final recursion level

-  ``new_element_position = new_element_position + 1``

-  ``derivative_structure(new_element_position, :) = current_derivative_term(:)``

Thereby this routine distributes perturbations on a new supertuple. In
every self-call the current perturbation of the recursion level (the
first one of pert) is put to another of the three elements of
``current_derivative_term``.

function ``get_num_blks``
-------------------------

Compares the different elements of the input perturbation tuple with
each other. Elements who are equal concerning perturbation operator and
frequency are put in one block

function ``get_blk_info``
-------------------------

Handles the ``blk_info``-array which has two dimensions: ``num_blks``
and 3. For every block it contains the index of the first perturbation
on the block, the number of perturbations on the block and the dimension
of the perturbation on the block.

function ``get_triangular_sizes``
---------------------------------

This function triggers a cascade of functions which at least determine
how many non-redundant perturbations there are on each block. In
contrast to the determination of the block size is takes into account
the dimension of the perturbation.

function ``get_triangulated_size``
----------------------------------

This function determines the number of non-redundant elements of the
result tensor.

General points for adopting the recursive code to handle residues
=================================================================

Things to be done for adopting the recursive scheme to the calculation
of residues:

#. Put calculation of excitation energies to work :math:`\rightarrow`
   **done**

#. Rewrite input parse such that an input for residues is recognized.

#. Use the manual specification scheme for response functions to also
   calculate the residues.

#. Modify the calculation of :math:`\mathbf{M}`\ , set up a
   functionality to contract it with :math:`\mathbf{X}_f^{*}` and to
   keep track of this result which is the right transition matrix
   element.

#. Introduce the calculation of :math:`\mathbf{D}^f` as a replacment for
   :math:`\mathbf{D}^{k\rightarrow f}_H`

#. Implement a functionality that identifies whether a term is needed
   for the residues of not (must check dependence on a density which
   depends on the excitation energy)

#. Put this function to work in the code and replace
   :math:`\mathbf{D}^{k\rightarrow f}_H` by :math:`\mathbf{D}^f` in the
   following subroutines:

   #. ``rsp_energy``

   #. ``rps_fock_lowerorder``

   #. ``rsp_idem_lag``

   #. ``rps_scfe_lag``

#. This will result in a functionality for calculation of residues
   w.r.t. electric field perturbation

#. Then go on with modifying the following routines

   #. ``rsp_pulay_kn``

   #. ``rsp_scfe_lag``

#. This should result in a fully applicable code

The remaining questions are the following

#. How to keep track of the excitation eigenvectors?

#. How to shape the output of the program?

#. How to keep track of double residues?

The selection of the relevant terms
-----------------------------------

For the calculation of residues several terms vanish from the expression
compared to the calculation of a response function. Namely all terms do
vanish which do not depend on the frequency or frequency sum which tend
towards the excitation energy.

Single residues as transition properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If transition properties are seeked, single residues have to be
calculated. For these residues the :math:`2n+1`\ -rule always has to be
kept since otherwise no appropriate decomposition in a left and a right
transition matrix element can be obtained. Concerning the involved
response function we then can assume to have a :math:`2n+2`\ -rule
concerning the linear equation systems that have to be solved since for
the perturbation combination with its frequency tending towards the
excitation energy no linear equation system has to be solved. Keeping
the :math:`2n+1`\ -rule for setting up the response function fundamental
for the formation of the seeked residue and setting all vanishing terms
to zero we find that there is only one term remaing which is not in
accordance with a :math:`2n+2`\ -rule. Nevertheless for this term the
response equations do not have to be solved.

We can summarize the rules for single residues in the following way:

-  :math:`2n+1`\ -rule for setting up the fundamental response function

-  :math:`2n+1`\ -rule for the solution of the response equations


