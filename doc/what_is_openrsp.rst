.. _section_what_is_openrsp:

What is OpenRSP?
================

OpenRSP is a program library that uses recursive routines to identify and
assemble contributions to response properties - that is, molecular properties
as they are expressed in the theory called "response theory" from theoretical
chemistry.

The name of OpenRSP reflects the following features:

* It is a library for the **Open**-ended calculation of **R**\ e\ **SP**\ onse
  properties: It can be used for the calculation of reponse properties to
  arbitrary order.
* It is **Open**-source and `is publicly available
  <https://github.com/openrsp/openrsp>`_  under the LGPL v2.1 license.
* It has an application programming interface that **Open**\ s it to connection
  with other programs that wish to make use of its functionality.

What are response properties?
-----------------------------

Response properties describe how fundamental properties of a molecular system
*respond* to external influences like subjection to an electromagnetic field or
displacement of the atomic nuclei.  They and related properties are essential
for the description of spectroscopic processes and molecular characteristics
like infrared spectroscopy, Raman scattering, multiphoton absorption and
vibrational energy levels. If you have ever done computational work on the
molecular level for phenomena in this category, chances are that response
properties were involved at some stage of the calculation.

Response properties can be categorized by their *order*, that is, the "number
of influences" that were taken into consideration for a given property. The
first such order is called *linear response* and contains much-used properties
like the electric dipole polarizability - i.e.  the first-order change to the
molecular dipole moment in the presence of an electric field - or the Hessian
matrix of nuclear geometric displacements - i.e. the change in the molecular
gradient that would result from displacing each coordinate of the molecular
geometry.

Higher orders of response properties describe the changes that the fundamental
molecular property would undergo upon subjection to more than one external
influence, or upon higher-order interactions with the same influence. Examples
of such properties are the geometric gradient of the electric dipole
polarizability - essential for the description of vibrational Raman spectra -
or the cubic and quartic force constants, i.e. the third- and fourth-order
derivatives of the molecular energy with respect to geometrical displacements -
which may be used to calculate corrections to a description of the vibrational
energy levels stemming from the geometric Hessian.

Why use OpenRSP?
----------------

By its recursive structure, OpenRSP makes it possible to calculate response
properties of arbitrary complexity in an analytical manner, not resorting to
numerical schemes like finite difference methods in the calculation. Compared
to analytical methods, numerical approaches may be associated with a greater
degree of uncertainty related to accuracy and practical feasibility of the
calculation, and we therefore think that analytical calculation should be used
whenever it is practical.

Today's programs written for the calculation of response properties may either
not have a recursive structure, or may use numerical methods to different
extents, or both. In the cases where existing programs use an analytical
approach, they may either be not recursive (which typically means that a new
program routine must be written for each new property for which calculation is
desired), or they may only be usable for a limited category of properties.  As
the complexity of the expressions that must be evaluated in an analytical
approach to yield the desired response property increases rapidly with the
order of response, such analytic calculation of high-order response properties
can quickly become a very complicated task and the implementation of *ad hoc*
program routines for their calculation may be intractable at higher orders.

The structure of OpenRSP, using recursion as a core tool, solves the task of
identifying and assembling contributions to response properties "once and for
all".  When combined with program libraries that can provide the contributions
that OpenRSP identifies, any response property can be calculated fully
analytically as long as those libraries can provide the necessary
contributions. We note, however, that the present version of the code is still
awaiting the completion of functionality to handle perturbations that both
change the basis set and have a nonzero frequency associated with them, but
that such extension is within the scope of the present underlying theory.
