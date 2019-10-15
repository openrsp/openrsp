.. _subsection_notations_and_conventions:

OpenRSP Notations and Conventions
---------------------------------

The following notations and conventions will be used through the OpenRSP
program and the documentation:

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

  **NOTE**: the above perturbtion free scheme is however not implemented for
  the current release so that OpenRSP will use its own internal representations
  for different perturbations.

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

Category of perturbation frequencies
  We use different integers for distinguishing different values of frequencies
  within a frequency configuration. The category arrary is determined by:

  #. For each frequency configuration, we start at the first perturbation and
     let its frequency value be designated number 1, then
  #. For the next perturbation,

     #. If its frequency value corresponds to a frequency value encountered
        previously in this frequency, then use the same designation as for that
        previously encountered frequency value, or
     #. If its frequency value has not been encountered before, then let that
        frequency value be designated with the first unused number;
  #. Continue like this until the end of the perturbation tuple;
  #. Start the numbering over again at the next frequency configuration.

Canonical order
  #. In OpenRSP, all perturbation tuples are canonically orderd according
     to the argument ``pert_tuple`` in the API :c:func:`OpenRSPGetRSPFun`
     or :c:func:`OpenRSPGetResidue`. For instance, when a perturbation
     tuple :math:`(a,b,c)` given as ``pert_tuple`` in the API
     :c:func:`OpenRSPGetRSPFun`, OpenRSP will use such order (:math:`a>b>c`)
     to arrange all perturbation tuples inside and sent to the callback functions.
  #. Moreover, a collection of several perturbation tuples will also follow
     the canonical order. For instance, a collection of all possible perturbation
     tuples of labels :math:`a,b,c,d` are
     :math:`(0,a,b,ab,c,ac,bc,abc,d,ad,bd,abd,cd,acd,bcd,abcd)`, where
     :math:`0` means unperturbed quantities that is always the first one
     in the collection.

     The rules for generating the above collection are:

     #. When taking a new perturbation into consideration, always do so in
        alphabetical order (and begin with the empty set);
     #. When taking a new perturbation into consideration, the new subsets are
        created by making the union of all previous subsets (including the
        empty set) and the new perturbation (putting the new perturbation
        at the end).

Perturbation :math:`a`
  The first perturbation label in the tuple sent to OpenRSP APIs
  :c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue`, are the
  perturbation :math:`a` [Thorvaldsen2008]_.

Perturbation addressing
  #. The addressing of perturbation labels in a tuple is decided by
     (i) the argument ``pert_tuple`` sent to the API :c:func:`OpenRSPGetRSPFun`
     or :c:func:`OpenRSPGetResidue`, and (ii) the canonical order that
     OpenRSP uses.
  #. The addressing of components per perturbation (several consecutive
     identical labels with the same complex frequency) are decided by
     the host program (**NOTE**: as mentioned above, the perturbtion free
     scheme is not implemented for the current release so that OpenRSP will use
     its own internal subroutines to compute the address of components per
     perturbation).
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
     one-dimensional array) is ``{1,[a],[b],[a][b],[c],[a][c],[b][c],[a][b][c]}``
     for :math:`(0,a,b,ab,c,ac,bc,abc)`, where as aforementioned the
     first one is the unperturbed quantities.
