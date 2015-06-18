.. _chapter-introduction:

What is OpenRSP?
================

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
