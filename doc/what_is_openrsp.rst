.. _chapter_what_is_openrsp:

What is OpenRSP?
================

OpenRSP is a program library that uses recursive routines to identify and
assemble contributions to molecular properties ("response properties") as they
are formulated in the theory called "response theory" which is used in
theoretical chemistry.

By its recursive structure, OpenRSP makes it possible to calculate response
properties (hereafter just "properties") of arbitrary complexity in an
analytical manner, not resorting to numerical schemes like finite difference
methods for the calculation of these properties.

Why use OpenRSP?
----------------

Compared to analytical methods, numerical approaches may be associated with a
greater degree of uncertainty related to accuracy and practical feasibility of
the calculation.

Today's programs written for the calculation of response properties either do
not have a recursive structure, or they use numerical methods to different
extents, or both. In the cases where existing programs use an analytical
approach, they are either not recursive (which typically means that a new
program routine must be written for each new property for which calculation is
desired), or they can only be used for a limited category of properties.

The structure of OpenRSP solves the task of identifying and assembling
contributions to molecular properties "once and for all". When combined with
program libraries that can provide the contributions that OpenRSP identifies,
any response property can be calculated fully analytically 
as long as those libraries can provide the necessary contributions.

