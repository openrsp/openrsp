.. _chapter-introduction:

Introduction
============

What is OpenRSP?
----------------

OpenRSP is a computer program that uses recursive routines to identify and
assemble contributions to molecular properties ("response properties") as they
are formulated in the theory "response theory" which is used in theoretical
chemistry.

By its recursive structure, OpenRSP should make it possible to calculate
molecular properties (hereafter just "properties") of arbitrary complexity in
an analytical manner. "An analytical manner" means that the program does not
resort to so-called numerical mathematical methods for the calculation of these
properties. Such numerical methods are regarded by the inventors as less suited
for this type of calculation because they, in comparison with analytical
methods, are associated with a greater degree of uncertaintly related to
accuracy and practical feasibility of the calculation.

Today's programs written for this purpose either do not have a recursive
structure, or they use numerical methods to different extents, or both. In the
cases where existing programs use an analytical approach, they are either not
recursive (which means that a new program subroutine must be written for each
new property for which calculation is desired), or they can only be used for a
limited category of properties. The structure of OpenRSP solves the task of
identifying and assemble contributions to molecular properties "once and for
all".

The OpenRSP project has been under development since the mid-2000's, but did
not feature a recursive structure from the beginning. The development of the
approach used in today's recursive structure was started in the autumn of 2011,
and the first results (i.e. the first molecular properties calculated and
regarded as suitable for publication) were obtained about one year later (close
to the end of 2012/beginning of 2013). Since then, work on developing a version
of the program suitable for publication has been ongoing. This new work has for
example entailed making it possible to calculate other categories of properties
(so-called residues of response properties) and to write code to make OpenRSP
more modular (i.e. cleanly separated from other programs) and able to more
easily communicate with other programs to which it could be feasible to connect
it. It is estimated that a finished version could be completed during the
course of 2015.

OpenRSP Notations and Conventions
---------------------------------

The following notations and conventions will be used through the OpenRSP
library and the documentation:

Perturbation
  is described by a label, a complex frequency and its order. Any two
  perturbations are different if they have different labels, and/or
  frequencies, and/or orders.

Perturbation label
  An integer distinguishing one perturbation from others; all *different*
  perturbation labels involved in the calculations should be given by calling
  the application programming interface (API)
  :c:func:`OpenRSPSetPerturbations`; OpenRSP will stop if there is any
  unspecified perturbation label given afterwards when calling the APIs
  :c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue`.

Perturbation order
  Each perturbation can acting on molecules once or many times, that is the
  order of the perturbation.

Perturbation components and their ranks
  Each perturbation may have different numbers of components for their
  different orders, the position of each component is called its rank.

  For instance, there will usually be :math:`x,y,z` components for the electric
  dipole perturbation, and their ranks are ``{0,1,2}`` in zero-based numbering,
  or ``{1,2,3}`` in one-based numbering.

  The numbers of different components of perturbations and their ranks are
  totally decided by the host program. OpenRSP will get such information from
  callback functions, that is OpenRSP itself is a perturbation free library.
  *FIXME: perturbtion free scheme not implemented yet*

Perturbation tuple
  An ordered list of perturbation labels, and in which we further require that
  *identical perturbation labels should be consecutive*. That means the tuple
  :math:`(a,b,b,c)` is allowed, but :math:`(a,b,c,b)` is illegal because the
  identical labels :math:`b` are not consecutive.

  As a tuple:

  #. Multiple instances of the same labels are allowed so that
     :math:`(a,b,b,c)\ne(a,b,c)`, and
  #. The perturbation labels are ordered so that :math:`(a,b,c)\ne(a,c,b)`
     (because their corresponding response functions or residues are in
     different shapes).

  We will sometimes use an abbreviated form of perturbation tuple as, for
  instance :math:`abc\equiv(a,b,c)`.

  Obviously, a perturbation tuple :math:`+` its corresponding complex
  frequencies for each perturbation label can be viewed as a set of
  perturbations, in which the number of times a label (with the same frequency)
  appears is the order of the corresponding perturbation.

Canonical order
  #. In OpenRSP, all perturbation tuples are canonically orderd according
     to the argument ``pert_tuple`` in the API :c:func:`OpenRSPGetRSPFun`
     or :c:func:`OpenRSPGetResidue`. For instance, when a perturbation
     tuple :math:`(a,b,c)` given as ``pert_tuple`` in the API
     :c:func:`OpenRSPGetRSPFun`, OpenRSP will use such order (:math:`a>b>c`)
     to arrange all perturbation tuples inside and sent to the callback functions.
  #. Moreover, a collection of several perturbation tuples will also follow
     the canonical order. For instance, a collection of all possible perturbation
     tuples of labels :math:`a,b,c` are :math:`(0,a,b,c,ab,ac,bc,abc)`, where
     :math:`0` means unperturbed quantities that is always the first one
     in the collection.

Perturbation :math:`a`
  The first perturbation label in the tuple sent to OpenRSP APIs
  :c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue`, are the
  perturbation :math:`a` [#]_.

.. [#] Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,
       Poul Jørgensen and Sonia Coriani, J. Chem. Phys., 129,
       214108 (2008).

Perturbation addressing
  #. The addressing of perturbation labels in a tuple is decided by
     (i) the argument ``pert_tuple`` sent to the API :c:func:`OpenRSPGetRSPFun`
     or :c:func:`OpenRSPGetResidue`, and (ii) the canonical order that
     OpenRSP uses.
  #. The addressing of components per perturbation (several consecutive
     identical labels with the same complex frequency) are decided by
     the host program. *FIXME: perturbtion free scheme not implemented yet*
  #. The addressing of a collection of perturbation tuples follows the
     canonical order as aforementioned.

  Therefore, the shape of response functions or residues is mostly
  decided by the host program. Take :math:`\mathcal{E}^{abbc}` as an 
  example, its shape is :math:`(N_{a},N_{bb},N_{c})`, where :math:`N_{a}`
  and :math:`N_{c}` are respectively the numbers of components of 
  the first order of the perturbations :math:`a` and :math:`c`, and
  :math:`N_{bb}` is the number of components of the second order of 
  the perturbation :math:`b`, and

  #. In OpenRSP, we will use notation ``[a][bb][c]`` for :math:`\mathcal{E}^{abbc}`,
     where the leftmost index (``a``) runs slowest in memory and the
     rightmost index (``c``) runs fastest. However, one should be
     aware that the results are still in a one-dimensional array.
  #. If there two different frequencies for the perturbation :math:`b`,
     OpenRSP will return ``[a][b1][b2][c]``, where ``b1`` and ``b2``
     stand for the components of the first order of the perturbation
     :math:`b`.
  #. The notation for a collection of perturbation tuples (still in a
     one-dimensional array) is ``{1,[a],[b],[c],[a][b],[a][c],[b][c],[a][b][c]}``
     for :math:`(0,a,b,c,ab,ac,bc,abc)`, where as aforementioned the
     first one is the unperturbed quantities.

