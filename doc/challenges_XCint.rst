

XCint
=====


General considerations
----------------------

The aim of XCint is to provide arbitrary-order derivatives of the
exchange-correlation (XC) energy and the XC potential.

.. math::

   E_{xc} = \int \, d r \, e_{xc} (n_1, \ldots, n_k)

.. math::

   V_{xc;\kappa\lambda} = \int \, d r \, \frac{\partial e_{xc}}{\partial D_{\kappa\lambda}}
                        = \int \, d r \, \sum_k \frac{\partial e_{xc}}{\partial n_k} \frac{\partial n_k}{\partial D_{\kappa\lambda}}

The energy density e_xc depends on a set of local variables, typically 1, 2, 4,
or 8 depending on whether we work with a LDA or GGA type functional, whether
spin polarization is considered or not, whether the reference is closed shell
or not. For meta-GGA functionals, local hybrids, current density functionals,
etc., the number of variable grows.

In real life we evaluate derivatives of E_xc and V_xc on a numerical grid:

.. math::

   E_{xc} = \sum_i w_i e_{xc} (n_1, \ldots, n_i)

.. math::

   V_{xc;\kappa\lambda} = \sum_i w_i \sum_i \frac{\partial e_{xc}}{\partial n_i} \frac{\partial n_i}{\partial D_{\kappa\lambda}}

Here w_i are the integration weights and the sum typically runs over 10000 to
few million points depending on the size of the molecule.

Densities are evaluated using density matrices and the overlap distribution
(product of two AOs), for example the (number) density:

.. math::

   n = \sum_{\kappa\lambda} \Omega_{\kappa\lambda} D_{\kappa\lambda}

The number density is the easiest density but all other densities (gradient,
current, kinetic energy density, spin density, ...) can be calculated using
very similar routines. It boils down to combining the corrent matrix sub-blocks
and using the appropriate AO derivatives.

The challenges in calculating arbitrary-order derivatives of the
exchange-correlation (XC) energy and the XC potential are:

* AO derivatives
* Functional derivatives
* Density derivatives
* Grid weight derivatives
* Assembling the partial derivatives to form directional derivatives
* Assembling the partial derivatives to form geometric derivatives (special case)
* Integration

In the following we will consider these challenges in detail.


AO derivatives
--------------

It would be easy to write a general routine that calculates arbitrary-order
derivatives of the AOs but this would be very inefficient. But we have a code
that can generate very efficient explicit code to arbitrary order in geometric
derivatives and to arbitrary order in angular momentum. Typically we stop
somewhere but one can go as high as needed. In short, this problem is solved.


Functional derivatives
----------------------

This is solved. XCFun provides functional derivatives to arbitrary order.


Density derivatives
-------------------

Pertubrations that do not modify the overlap are easy. They can be done to
arbitrary order. Geometric and magnetic perturbations are tricky.  Presently we
have code to do n^B, n^BF, n^G, n^GG (almost), n^GF, n^GFF.  We have the
"blueprints" for higher order routines but need to program them.  It is
possible to program them to arbitrary order but it is not done yet.


Grid weight derivatives
-----------------------

Presently XCint does not generate the grid but receives it from outside.  This
has advantages but the downside is that the weights are generated without
knowledge of the perturbation. This can be a problem for geometric derivatives.
Presently the contribution from grid weight derivatives is not considered.


Assembling the partial derivatives to form directional derivatives
------------------------------------------------------------------

Once all perturbed densities are known to the necessary order, then they can be
fed to XCFun as taylor series expansions and XCFun returns either the final
"kernel density" which can be easily distributed over the V_xc matrix or
deliveres the final e_xc density derivative (the latter is not tested but in my
opinion should work).


Integration
-----------

Presently we do not have a numerical grid which adapts to the perturbation.
This is probably not a problem for geometric derivatives.
