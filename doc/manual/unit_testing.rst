.. _chapter_unit_testing:

Unit Testing
============

After successfully building OpenRSP (see :ref:`chapter_installation`), we recommend
users perform the unit testing of OpenRSP.

If ``OPENRSP_TEST_EXECUTABLE`` enabled, you will have an executable
``openrsp_c_test`` after successfully build ing OpenRSP. Run this executable
for unit testing.

If ``OPENRSP_TEST_EXECUTABLE`` disabled, you need to call the function

.. c:function:: QErrorCode OpenRSPTest(FILE *fp_log)

to perform the unit testing.
