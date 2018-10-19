


Keyword reference
=================

.CUSTOM
-------

Description...

.FREQ
-------

Description...

The keyword is followed by at least 2 lines::

  .FREQ
  1
  0.02

.SPORDR
-------

Description...

The keyword is followed by 1 line::

  .SPORDR
  2

.SPPLAB
-------

Description...

The keyword is followed by at least 1 line::

  .SPPLAB
  EL
  GEO

.SPRULE
-------

Description...

The keyword is followed by 2 lines::

  .SPRULE
  1
  1

.THRESH
-------

Description...

The keyword is followed by 1 line::

  .THRESH
   1.0d-5

.XCGRID
-------

Use the XCint grid instead of the default Dalton grid.

The radial grid is generated according to Lindh, Malmqvist, and Gagliardi
`[TCA 106(3), 178-187 (2001)] <http://dx.doi.org/10.1007/s002140100263>`_.

The angular grid is generated according to
Lebedev and Laikov
[A quadrature formula for the sphere of the 131st
algebraic order of accuracy,
Russian Academy of Sciences Doklady Mathematics,
Volume 59, Number 3, 1999, pages 477-481].

The keyword is followed by 3 lines::

  .XCGRID
   1.0d-12        # radial precision
   86             # minimum number of angular points per radial shell
   302            # maximum number of angular points per radial shell

The smaller the radial precision, the better.

The higher the values for minimum and maximum number of angular points, the better.

For the minimum and maximum number of angular points the code will use the following
table and select the closest number with at least the desired precision::

     {6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
    194,  230,  266,  302,  350,  434,  590,  770,  974, 1202,
   1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,
   4934, 5294, 5810}

The pruning is a primitive linear interpolation.
The full angular grid is reached at 0.2 times the Bragg radius of the center.
