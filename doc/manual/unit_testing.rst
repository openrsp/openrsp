.. _chapter_unit_testing:

Unit Testing
============

After successfully building OpenRSP (see :ref:`chapter_installation`), we recommend
users perform the unit testing of OpenRSP.

If ``OPENRSP_TEST_EXECUTABLE`` is enabled, you will have an executable
``openrsp_c_test`` after successfully building OpenRSP. Run this executable for
unit testing.

If ``OPENRSP_TEST_EXECUTABLE`` is disabled, you will need to call the function

.. c:function:: QErrorCode OpenRSPTest(FILE *fp_log)

to perform the unit testing.
