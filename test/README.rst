Description of the test procedure
=================================


How to run tests
----------------

We use the runtest script to run tests.
Go to the source root::

  $ cd /path/to/source

Then get a list of options::

  $ ./runtest --help

To run all tests, type::

  $ ./runtest --all

You can also run individual tests. Have a look at the options.

If you are in the subdirectory of an individual test (example: dft_ac), you can
do this instead to run just this one test::

  $ [cd test/dft_ac]
  $ ../../runtest


How to run tests using CTest
----------------------------

Go to the build directory::

  $ cd build

Now run all tests::

  $ make test

You can also run individual tests (example: dft_ac)::

  $ ctest -R dft_ac

The "-R" means regular expression! So you can run all tests
that start with "dft" using::

  $ ctest -R dft

By the way behind the curtains CTest uses runtest.


Some tests are broken, what now?
--------------------------------

You are a user? Then please inform the developers via the mailing list.

You are a developer? Perhaps it is not your fault. Have a look at the nightly
testing dashboard where you can check whether this test is marked as red
(broken) or green (passing). If it is green there but broken on your machine
it can either mean that the test has trouble with your architecture/compiler
(in this case inform other developers) or that you broke it with recent changes
- then you have to go in and fix it before you git push to origin/master.


Suddenly some red tests appear on the nightly testing dashboard, what happens now?
----------------------------------------------------------------------------------

Miro or Radovan or somebody else examining the nightly testing dashboard will
notice it and write to the developer who is likely to have introduced the
conflicting change.


How to add a new test
---------------------

1. Copy an existing test and choose a sensible name, e.g. not
   my_favorite_test or test2, but for instance symmetry_recognition to signal that
   we are testing the symmetry recognition procedure.

2. Put the input files in the directory and modify the
   menu file to work with these filenames.

3. Put (a) correct output(s) in the result directory and modify the filter file
   (see further below) to take out the part that is relevant for the functionality
   that you want to test.  Don't test everything, e.g. the electronic Hartree-Fock
   energy is already tested at several places, no need to do this in every test.

4. Go to the source root directory and run the script runtest with your new test.  It
   is a good practice to test your test, try to introduce deliberately errors in
   your new code and check whether your new test catches the error - if not,
   modify the test. The script will complain if your filter does not filter out anything.

5. If the test works, look at how long it takes to run on a single PC to see where
   it should be quick or regular. A rule of thumb: < 10 seconds = quick, otherwise regular.

6. Add the test to pamadm (search for existing tests to see how).

7. Use git add and git commit to commit your work to the repository.  Make sure
   to include ALL files, don't forget the result directory that should contain the
   reference output and filter files. You can verify this with "git status".


How does the filtering work?
----------------------------

If one would use simply the unix command diff to compare the new output with
the reference output, one would always see differences, simply because dates
and execution times will be different, because insignificant last digits are
different, or because output formatting has changed.

To focus on relevant differences runscript extracts the relevant numbers
according to keywords in result/filter.  Here you may either use strings or
regular expressions, examples of both can be found when looking at the filter
files.  There are two possibilities, you either know in advance how many lines
you want to test, or you don't (keep in mind that others may later introduce
additional blank or text lines). In the first case you can simply specify the
number of lines to be test (using nlines), in the second case you have to
specify a characteristic text that closes the section.

Note that runtest compares relative errors (DIRAC) and absolute errors (Dalton).
The default accepted relative/absolute error (DIRAC/Dalton) is 0.0001.
You can modify this using $ERROR. If you want to ignore sign, use $IGNORESIGN.
Also not that $ERROR and $IGNORESIGN are only valid for the next subtest.
Example::
  $ERROR 0.000000001
  nlines 8:Electronic energy

  $ERROR 1.E-5
  from:Eigenvalues
  to:HOMO - LUMO

In this example we have two subtest. First we test Electronic energy and 7
following lines using tolerance 0.000000001, then we test Eigenvalues using
tolerance 1.E-5.  We need to specify 1.E-5 for the second subtest otherwise
runtest defaults to 0.0001.
Possible keywords are::
  $ERROR          # accepted error
  $IGNORESIGN     # ignore sign
  from:STRING     # start from string STRING
  fromre:STRING   # start from regex STRING
  to:STRING       # stop at string STRING, do not include end line
  toi:STRING      # stop at string STRING, include end line
  tore:STRING     # stop at regex STRING, do not include end line
  toire:STRING    # stop at regex STRING, include end line
  nlines N:STRING # consider line containing STRING and N-1 following lines
