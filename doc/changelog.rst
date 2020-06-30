

Version history and changelog
=============================


Version 1.0.0 (2020-06-30)
--------------------------

Code
~~~~

- Now requires Fortran 2008
- Rewrote linked list functionality for caching to instead use (reallocating) arrays
- "Number of components" marker in rsp_tensor output file now written as 'NUM_COMPONENTS' instead of 'NUM COMPONENTS'
- Significantly decreased usage of array constructors in function/subroutine arguments
- Now compiles and runs with most compilers but still problems with some Intel/2018 and Intel/2019 setups
- Fixed various memory leak/out-of-bounds errors that sometimes would happen
- OpenRSP now looks for available file units before choosing one to use
- Disabled internal memory limit and memory bookkeeping, may be reinstated later
- If a response tensor is large, then if it's printed at the debugging print level, it's broken down into smaller chunks
- Added stops for various currently unsupported residue calculation setups
- Fixed a bug concerning testing of perturbation frequencies against excitation energy for residue calculations
- Removed some unused residue-related routines
- Calculation setup errors encountered in wrapper routines now cause exit, not just warning and then continuing
- The "excitation" perturbation in a single residue calculation is now given
  the label EX1 in the rsp_tensor file; however, its current implementation
  still results in triplication of the calculation result data due to being
  treated as having three components when in fact it has got only one


Project
~~~~~~~

- Added contribution guide and authorship process guide
- Updated pull request template to solicit agreement to contribution terms
- Various changes to documentation


Version 1.0.0-alpha (2018-11-19)
--------------------------------

New
~~~

- Implemented application programming interface and corresponding developer
  manual using literate programming
- Adopted new file format for printing of final results
- Added support for caching of intermediate contributions and restarting an interrupted calculation
- Implemented recurse-calculate-recurse approach for most contributions
- Added support for passing several (pairs of) arguments for contraction with
  perturbed contributions depending to first (second) order on the
  perturbed/unperturbed density matrix, implemented a similar scheme for Pulay
  and Lagrange-type contributions
- Added support for calculation of several properties in one run with reuse of common intermediate results
- Added support for calculation of single residues of electric dipole polarization properties
- Started using callback function scheme for external contributions and added
  application programming interface: Callback functions now fulfill the role
  previously played by interface files (2015-02-09)
- General response code added (2012-03-19)
- Repository initialized (2010-05-23)


Changed
~~~~~~~

- This is the first changelog entry, so no changes to be mentioned here.
