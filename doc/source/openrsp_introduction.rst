.. _chapter-introduction:

Introduction
============

What is OpenRSP?
----------------

OpenRSP is a computer program that uses recursive routines to identify
and assemble contributions to molecular properties ("response properties")
as they are formulated in the theory "response theory" which is used in
theoretical chemistry.

By its recursive structure, OpenRSP should make it possible to calculate
molecular properties (hereafter just "properties") of arbitrary complexity
in an analytical manner. "An analytical manner" means that the program
does not resort to so-called numerical mathematical methods for the
calculation of these properties. Such numerical methods are regarded by
the inventors as less suited for this type of calculation because they,
in comparison with analytical methods, are associated with a greater
degree of uncertaintly related to accuracy and practical feasibility of
the calculation.

Today's programs written for this purpose either do not have a recursive
structure, or they use numerical methods to different extents, or both.
In the cases where existing programs use an analytical approach, they are
either not recursive (which means that a new program subroutine must be
written for each new property for which calculation is desired), or they
can only be used for a limited category of properties. The structure of
OpenRSP solves the task of identifying and assemble contributions to
molecular properties "once and for all".

The OpenRSP project has been under development since the mid-2000's, but
did not feature a recursive structure from the beginning. The development
of the approach used in today's recursive structure was started in the
autumn of 2011, and the first results (i.e. the first molecular properties
calculated and regarded as suitable for publication) were obtained about
one year later (close to the end of 2012/beginning of 2013). Since then,
work on developing a version of the program suitable for publication has
been ongoing. This new work has for example entailed making it possible
to calculate other categories of properties (so-called residues of response
properties) and to write code to make OpenRSP more modular (i.e. cleanly
separated from other programs) and able to more easily communicate with
other programs to which it could be feasible to connect it. It is estimated
that a finished version could be completed during the course of 2015.

OpenRSP Notations and Conventions
---------------------------------

The following notations and conventions will be used through the OpenRSP
library and the documentation: 

Perturbation label
  An integer describing what kind of perturbation added to the molecule; all
  *different* perturbations involved in the response theory calculations should
  be given by calling the function :py:meth:`OpenRSPSetPerturbations`; OpenRSP
  will stop if there is any unspecified perturbation label given afterwards when
  calling the function :py:meth:`OpenRSPGetRSPFun` or :py:meth:`OpenRSPGetResidue`.

Perturbation tuple
  An ordered list of perturbation labels that are added to the molecule.
  As a tuple:

  #. Multiple instances of the same labels are allowed so that
     :math:`(a,b,b,c)\ne(a,b,c)`, and
  #. The perturbation labels are ordered so that :math:`(a,b,c)\ne(a,c,b)`
     (because their corresponding response functions or residues are in
     different shapes).

  We will in the following use an abbreviated form of perturbation tuple as,
  for instance :math:`abc\equiv(a,b,c)`.

Perturbation addressing
  #. The perturbation labels in a tuple are always ordered as that of
     the argument ``pert_labels`` given in the function
     :py:meth:`OpenRSPSetPerturbations`. For example, OpenRSP will return
     response functions :math:`\mathcal{E}^{abbc}` if ``pert_labels = (a,b,c)``,
     and :math:`\mathcal{E}^{acbb}` if ``pert_labels = (a,c,b)``, and
     :math:`\mathcal{E}^{bbac}` if ``pert_labels = (b,a,c)``, and so on.
  #. For each perturbation, there may be different number of components.
     For instance, there will usually be :math:`x,y,z` components for
     the electric perturbation in dipole approximation. The number of
     different components for each perturbation and how these components
     are ordered are decided also by the host program. *FIXME: not implemented yet*
  #. Therefore, the shape of the response functions or residues is totally
     decided by the host program. Take :math:`\mathcal{E}^{abbc}` as an
     example, its shape is :math:`(N_{a},N_{bb},N_{c})`, where :math:`N_{a}`
     and :math:`N_{c}` are respectively the numbers of components of
     the first order of perturbations :math:`a` and :math:`c`, and
     :math:`N_{bb}` is the number of components of the second order of
     perturbation :math:`b`. We will use the notation ``r(a,bb,c)`` for
     the results (response functions or residues), where the leftmost
     index (``a``) runs fastest in memory and the rightmost index (``c``)
     runs slowest.
  #. If there two different frequencies for perturbation ``b``, OpenRSP
     will return ``r(a,b1,b2,c)`` that ``b1`` and ``b2`` stand for the
     components of the first order of perturbation :math:`b`.
