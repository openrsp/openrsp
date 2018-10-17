TODO List
=========

.. Give future long- and short-term developments

Short-term
----------

To be done before release:

1. Look through RST files in the root directory and "doc" directory
2. Update "ChangeLog" if we still keep it
3. Update "openrsp.bib" if we still keep it
4. If time allowed, make the unit testing easier and optionally add more unit tests

#. Removes the frequency of perturbation :math:`a` in ``OpenRSPGetRSPFun``
#. Adds nuclear contributions in tests
#. Adds DFT contributions in tests
#. Adds the API for residue calculations
#. Adds integrals and expectation values in callback functions?
#. Gradually moves files in ``src/interfaces`` and ``src/input`` to DALTON
   and LSDALTON, and removes/merges files in ``src/legacy`` and ``src/TODO``
#. Prepares unit testing
#. Prepares tutorial and finish documentation
#. Adds PCM and PE contributions (in host programs)

Long-term
---------

#. Changes to perturbation free scheme
#. Implements the symbolic operations
#. Adds response theory based on molecular orbital coefficients? coupled cluster?
