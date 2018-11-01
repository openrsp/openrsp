.. _chapter_residues:

Residues
========

The residues can be calculated by calling :c:func:`OpenRSPGetResidue`. It can
be noted that this API takes similar arguments as those of
:c:func:`OpenRSPGetRSPFun`. Users can therefore refer to
:ref:`chapter_response_function` for the explanation of the arguments with the
same names, or :ref:`chapter_api_reference` for the description of
:c:func:`OpenRSPGetResidue`.

In the following, we will focus on the arguments that do not present in
response function calculations, as grouped into 3 categories:

#. ``order_residue``,
#. ``num_excit``, ``excit_energy`` and ``eigen_vector``,
#. ``residue_num_pert`` and ``residue_idx_pert``.

The arguments ``size_residues`` and ``residues`` are merely the size of the
calculated residues and the residues, whose description can be found in
:ref:`chapter_api_reference`.

First, ``order_residue`` is simply the order of residues. That is, there will
be ``order_residue`` frequencies (in ``pert_freqs``) of perturbation labels (or
sums of frequencies of perturbation labels) respectively equal to the
``order_residue`` excitation energies.

The excitation energies will be in ``excit_energy``, corresponding to
``order_residue`` excitations (we say an "excitation tuple"). Considering the
possibility that calculations of the same order residues several "excitation
tuples" may save time, OpenRSP will therefore accept for multiple "excitation
tuples" for residue calculations.

Let us make an example: we want to calculate double residuces, that is,
``order_residue=2``, and we have two excited states ``r1`` and ``s1`` (or an
excitation tupe ``(r1,s1)``) that two (sums of) frequencies of perturbation
labels will equal to their excitation energies.

Now if we have other two excited states ``(r2,s2)`` and ``(r3,s3)`` that will
be used for double residue calculations, we could set ``num_excit=3`` for the
number of excitation tuples used in calculations.

The ``excit_energy`` and ``eigen_vector`` respetively contain excitation
energies and eigenvectors of all the excited states:

#. ``excit_energy[order_residue*num_excit]={energy_r1,energy_s1, energy_r2,energy_s2, energy_r3,energy_s3}``,
#. ``eigen_vector[2*3]={eigvec_r1,eigvec_s1, eigvec_r2,eigvec_s2, eigvec_r3,eigvec_s3}``.

The last group of arguments ``residue_num_pert`` and ``residue_idx_pert``
contain the information, **per property**, of perturbation labels whose (sums
of) frequencies equal to the excitation energies.

For instance, we have two properties to calculate ``num_props=2``, and the
lengths of perturbation tuples for them are ``len_tuple[2]={4,5}``.

For the first property, we want (for every excitation tuple):

* frequency of the first perturbation label :math:`=` the first excitation energy,
* frequency of the third perturbation label :math:`=` the second excitation energy.

For the second property, we want (for every excitation tuple):

* frequency of the first perturbation label :math:`=` the first excitation energy,
* sums of frequencies of the third and fourth perturbation labels :math:`=` the
  second excitation energy.

So, we will have:

#. ``residue_num_pert[order_residue*num_props]={1,1, 1,2}``,
#. ``residue_idx_pert[sum(residue_num_pert)]={1,3, 1,3,4}``,

or in zero-based numbering:

#. ``residue_num_pert[order_residue*num_props]={1,1, 1,2}``,
#. ``residue_idx_pert[sum(residue_num_pert)]={0,2, 0,2,3}``.

Last but not least, it is apparent that the (sums of) frequencies
(``pert_freqs``) of perturbation labels specified by ``residue_idx_pert``
should equal to the corresponding excitation energies for all frequency
configurations and excitation tuples.
