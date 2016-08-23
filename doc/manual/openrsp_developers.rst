.. _chapter-openrsp-developers:

OpenRSP for Developers
======================

*FIXME: this chapter will be removed and its contents will be in the OpenRSP developer manual.*

In this chapter, we will discuss the idea behind the implementation, that
will be useful for further development and maintenance, and be useful for
new developers to understand the library and to start their work on top
of the current development.

Theoretical Background
----------------------

The density matrix-based quasienergy formulation of the Kohn-Sham density
functional response theory using perturbation- and time-dependent basis
sets can be found in:

1. Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen, Poul Jørgensen
   and Sonia Coriani, J. Chem. Phys., 129(21), 214108 (2008).
2. Radovan Bast, Ulf Ekström, Bin Gao, Trygve Helgaker, Kenneth Ruud and
   Andreas J. Thorvaldsena, Phys. Chem. Chem. Phys., 13, 2627-2651 (2011).

A relativistic implementation can be found in:

3. Radovan Bast, Andreas J. Thorvaldsen, Magnus Ringholm and Kenneth Ruud,
   Chem. Phys., 356(1-3), 177-186 (2009).

The recursive programming techniques implemented in OpenRSP can be found in:

.. _Ringholm2014:

4. Magnus Ringholm, Dan Jonsson and Kenneth Ruud, J. Comput. Chem., 35(8),
   622-633 (2014).

The recursive programming techniques used for the first order residues can
be found in:

5. Daniel H. Friese, Maarten T. P. Beerepoot, Magnus Ringholm and Kenneth Ruud,
   J. Chem. Theory Comput., 11(3), 1129-1144 (2015).

Framework of OpenRSP
--------------------

The name OpenRSP stands for "open-ended response theory", that is,
the library is:

#. open-ended for different levels of theory, i.e., one-, two- and
   four-component levels;
#. open-ended for different wave functions, e.g., atomic-orbital (AO)
   based density matrix, molecular orbital (MO) cofficients and
   coupled cluster (CC);
#. open-ended for different kinds of perturbations; and
#. open-ended for different host programs.

For the time being, OpenRSP has implemented:

#. AO based density matrix response theory (source codes [#]_ in
   ``src/ao_dens``),

.. [#] The codes in ``src/ao_dens`` are written in Fortran, but OpenRSP
       APIs are implemented using C language. Therefore, adapter codes
       between them are implemented in ``src/ao_dens/adapter``, for OpenRSP
       APIs calling the codes of AO based density matrix response theory,
       also for the AO based density matrix response theory codes calling
       the callback functions (as function pointers saved by OpenRSP APIs).

and it works for one-, two- and four-component levels by simply setting
the appropriate Hamiltonian. We are now planning to implement the MO and
CC based response theories.

To make OpenRSP work for any perturbation, we are now trying to implement
the so called **perturbation free scheme**, see Section
:ref:`openrsp-c-support-perturbations`.

In order to make it easy for implementing OpenRSP into different host
programs (written in different programming languages), we agree to use
the **callback function scheme** in OpenRSP in the 2015 Skibotn meeting.
The callback functions are specified by host programs by calling the
OpenRSP APIs (both C and Fortran APIs implemented) during run time,
and will be used by OpenRSP during calculations, to get contributions
from electronic and nuclear Hamiltonian, and to get response parameters
from solving the linear response equation.

Another important issue affects the implementation of OpenRSP into different
host programs is the matrix and its different operations that OpenRSP
extensively depends on. Different host programs can have different types
of matrices (dense and sparse, sequential and parallel) and written by
different programming languages (e.g. C and Fortran).

To best utilize the host program's developed matrix routines (if there is),
and also to remove this complexity of matrix problem from OpenRSP, we also
agree to build OpenRSP on top of the `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_
in the 2015 Skibotn meeting. This matrix library works as an adapter
between OpenRSP and different matrix routines (implemented in different
host programs) that can be written in C and Fortran [#]_.

.. [#] If there is no matrix routines implemented in a host program, it
       can fully use the QcMatrix library that will invoke BLAS and LAPACK
       libraries for matrix operations.

Therefore, a full picture of OpenRSP used in a C host program can be
(the use of OpenRSP in a Fortran host program can be found in Secion
:ref:`section-openrsp-Fortran-APIs`):

.. figure:: /_static/OpenRSP_framework.pdf
   :scale: 100 %
   :align: center

As shown in the above picture, the OpenRSP library is divided into three parts:

#. The "OpenRSP C APIs" have been described in Chapter
   :ref:`chapter-API-reference` which work mostly between the host program driver
   routine and other parts of the OpenRSP library;
#. The "OpenRSP response" is the core part in which the performance of response
   theory will be done;
#. The "OpenRSP C support" will be described in Sections
   :ref:`openrsp-c-support-perturbations`, :ref:`openrsp-c-support-Hamiltonian`
   and :ref:`openrsp-c-support-solver`, which saves the information of
   perturbations, electronic and nuclear Hamiltonian and linear response equation
   solver, and will be used by the "OpenRSP response" part during calculating
   response functions and residues.

Symbolic Computations (not implemented)
---------------------------------------

The recursive programming techniques described in Ref. [Ringholm2014_] make it
possible to calculate molecular properties of arbitrary complexity in an
analytical manner. But if we are going to implement MO and CC based response
theories, it may mess up or not be easy to (re)use the already developed codes.

Therefore, I (Bin Gao) will examine an alternative route to the MO and
CC based response theories, by developing a set of symbolic computation
functions that will be used by response theory codes to get molecular
properties still in a recursive and analytical manner.

These symbolic computation functions will be implemented using the literate
programming approach. In this way, I can describe the idea of implementation
and the real codes at the same time. After becoming familar with the literate
programming approach, it is actually, in my opinion, an efficient way for
methodology development.

The documentation of these symbolic computation functions can be generated from
the ``*.nw`` files in the directory ``web``. Actually, the source codes
(``src/*.c``) can also be generated from the ``*.nw`` files, but I think it
could be better that we keep the generated C codes in the OpenRSP repository
and release with OpenRSP. Because releasing only the ``*.nw`` files will impose
further requirements for the users, they will have to install some software for
literate programming.

.. _openrsp-c-support-perturbations:

Perturbation Free Scheme (not implemented)
------------------------------------------

For different perturbations, there could be **different numbers of components**
and **arranged in different ways** in different host programs. For instance,
there are 9 components for the second order magnetic derivatives in a redundant
way :math:`xx,xy,xz,yx,yy,yz,zx,zy,zz`, but 6 components in a non-redundant way
:math:`xx,xy,xz,yy,yz,zz`. There are at most four centers in different
integrals, non-zero high order (:math:`\ge 5`) geometric derivatives are only
those with at most four differentiated centers.

To take all the above information into account in OpenRSP will make it so
complicated and not necessary, because response theory actually does not depend
on the detailed knowledge of different perturbations. In particular, when all
the (perturbed) integrals and expectation values are computed by the host
program's callback functions, the detailed information of perturbations, i.e.
the number of components and how they are arranged in memory can be hidden from
OpenRSP.

The former can be easily solved by sending the number of components of each
perturbation (label) up to its maximum order to the OpenRSP API
:c:func:`OpenRSPSetPerturbations`.

The latter can be important for OpenRSP, for instance, when the higher order
derivatives with respect to **one perturbation** need to be constructed from
several lower order derivatives. For instance, the second order derivatives may
be constructed from the first order ones in the redundant format:

* :math:`x+x\rightarrow xx,\hspace*{2em}0+0\rightarrow 0`,
* :math:`x+y\rightarrow xy,\hspace*{2em}0+1\rightarrow 1`,
* :math:`x+z\rightarrow xz,\hspace*{2em}0+2\rightarrow 2`,
* :math:`y+x\rightarrow yx,\hspace*{2em}1+0\rightarrow 3`,
* :math:`y+y\rightarrow yy,\hspace*{2em}1+1\rightarrow 4`,
* :math:`y+z\rightarrow yz,\hspace*{2em}1+2\rightarrow 5`,
* :math:`z+x\rightarrow zx,\hspace*{2em}2+0\rightarrow 6`,
* :math:`z+y\rightarrow zy,\hspace*{2em}2+1\rightarrow 7`,
* :math:`z+z\rightarrow zz,\hspace*{2em}2+2\rightarrow 8`,

where we have ranked different components in zero-based numbering (numbers on
the right). However, the ranks can be different in different host programs. To
solve this problem, i.e., the mapping relationship of lower and higher order
derivatives with respect to **one perturbation** [#]_, we ask for a callback
function :c:func:`get_pert_concatenation` from host programs, which is the last
argument of the API :c:func:`OpenRSPSetPerturbations`.

.. [#] We emphasize the derivatives of **one perturbation** because
       components of higher order derivatives of different perturbations
       are simply the direct product of components of lower order derivatives.

This callback function is used by OpenRSP to get the ranks of components of
*sub-perturbation tuples with same perturbation label* (lower order derivatives
with respect to one perturbation) for given components of the corresponding
*concatenated perturbation tuple* (higher order derivatives).

.. _openrsp-c-support-Hamiltonian:

Electronic and Nuclear Hamiltonian
----------------------------------

As aforementioned, the ingradients of electronic and nuclear Hamiltonian are
sent to OpenRSP, and (perturbed) integrals and expectation values will be
computed by the callback functions of host programs. These include:

#. overlap integrals (source codes ``web/RSPOverlap.nw``),
#. one-electron operators (source codes ``web/RSPOneOper.nw``),
#. two-electron operators (source codes ``web/RSPTwoOper.nw``),
#. exchange-correlation functionals (source codes ``web/RSPXCFun.nw``),
#. nuclear Hamiltonian (source codes ``web/RSPNucHamilton.nw``),

where the source codes save the callback functions as function pointers in
different C ``struct``, and take care the invoking of these callback functions
during calculations.

Different from the overlap integrals and nuclear Hamiltonian, the one- and
two-electron operators and XC functionals are saved in three different linked
lists in OpenRSP, in which each node corresponds to an operator. This makes it
possible for host programs to add different callback functions for different
operators, if they do not want to or can not provide OpenRSP a general callback
function.

.. _openrsp-c-support-solver:

Response Equation Solver
------------------------

Similar to overlap integrals and nuclear Hamiltonian, the callback function of
a linear response equation solver is saved as a function pointer in a C ``struct``
in OpenRSP. That will be invoked by OpenRSP for obtaining response parameters,
and the source codes related to the solver are in ``web/RSPSolver.nw``.

OpenRSP will send multiple RHS vectors (or matrices) to the solver, for several
frequency sums on the left hand side of the linear response equation and for
several derivatives with respect to (different) perturbations.

Notice that it would be more common and help the convergence to calculate
several frequencies for the same perturbation, than the other way around. So
the RHS matrices and response parameters are arranged as
``[num_comps][num_freq_sums]`` in the callback function
:c:func:`get_linear_rsp_solution`.

.. _section-openrsp-Fortran-APIs:

Implementation of Fortran APIs
------------------------------

OpenRSP APIs that host programs will use to talk to OpenRSP are written in C
language, with Fortran support by using Fortran 2003 language. The source codes
are in ``web/FortranAPIs.nw``, and the framework of OpenRSP used in a Fortran
host program is shown in the following figure:

.. figure:: /_static/OpenRSP_Fortran_API.pdf
   :scale: 100 %
   :align: center

Two new parts are needed for the use of OpenRSP in a Fortran program:

#. "OpenRSP Fortran APIs", and
#. "OpenRSP Fortran support".

Take one-electron integrals as an example, the callback subroutine
:c:func:`get_one_oper_mat` is declared in the ``interface`` of the OpenRSP
Fortran API ``OpenRSPAddOneOper_f()``::

    function OpenRSPAddOneOper_f(...)
        interface
            subroutine get_one_oper_mat(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
                                        num_int,          &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_one_oper_mat
        end interface
    end function OpenRSPAddOneOper_f

But "OpenRSP C support" codes can not call this subroutine
:c:func:`get_one_oper_mat` directly, because the ``type(QcMat)`` can not be
sent from a C function to a Fortran subroutine directly. Instead, another
"OpenRSP Fortran support" subroutine is implemented in OpenRSP that will be
called by the "OpenRSP C support" codes::

    subroutine RSPOneOperGetMat_f(oper_num_pert,    &
                                  oper_pert_labels, &
                                  oper_pert_orders, &
                                  user_ctx,         &
                                  num_int,          &
                                  val_int)          &
        bind(C, name="RSPOneOperGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: oper_num_pert
        integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
        integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(OneOperFun_f), pointer :: one_oper_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_val_int(:)     !integral matrices
        integer(kind=4) ierr                         !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_val_int(num_int), stat=ierr)
        ... ...
        ierr = QcMat_C_F_POINTER(A=f_val_int, c_A=val_int)
        ... ...
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, one_oper_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call one_oper_fun%get_one_oper_mat(oper_num_pert,         &
                                           oper_pert_labels,      &
                                           oper_pert_orders,      &
                                           ... ...,               &
                                           num_int,               &
                                           f_val_int)
        ! cleans up
        nullify(one_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        ... ...
        deallocate(f_val_int)
    end subroutine RSPOneOperGetMat_f

As shown above, the important thing here is to use the QcMatrix function
``QcMat_C_F_POINTER`` converting an array of C pointers ``val_int`` to an array
of Fortran ``type(QcMat)`` variables ``f_val_int``. For sure, these two point
to the same memory so that any manipulation on the latter equals to that on the
former. Another QcMatrix function ``QcMat_C_NULL_PTR`` is used to clean up the
context of Fortran ``type(QcMat)`` variables ``f_val_int`` (but not the C
pointers ``val_int``).

The procedure when doing a callback of Fortran subroutine can be summarized as:

"OpenRSP response" codes (Fortran) :math:`\Rightarrow` "OpenRSP C support"
codes :math:`\Rightarrow` "OpenRSP Fortran support" subroutine
``RSPOneOperGetMat_f()`` :math:`\Rightarrow` :c:func:`get_one_oper_mat`

One can also notice that, the argument ``num_int`` is needed in the
``interface`` of ``OpenRSPAddOneOper_f()`` and the subroutine
``RSPOneOperGetMat_f()``, and "OpenRSP C support" codes also need to pass
``num_int`` to ``RSPOneOperGetMat_f()`` (from C function to Fortran
subroutine). Therefore, these arguments for the dimension of arrays have to be
passed although they are over complete.

Technical Issues in OpenRSP
---------------------------

#. In OpenRSP APIs (C), we choose to represent complex numbers as their real
   and imaginary parts in an array. It might be efficient for host programs'
   integral codes that all real parts of numbers are put together and
   followed by all imaginary parts, but this loss the requirement that OpenRSP
   works with complex numbers, not an array with real and imaginary parts.

   *FIXME: in OpenRSP Fortran APIs, we should choose complex numbers, right? Because Fortran support complex numbers.*

#. We can simply add the following into ``include/OpenRSP.h``, to make OpenRSP
   be called by C++ programs::

     #ifdef __cplusplus
     extern "C" {
     #endif

     ... ...

     #ifdef __cplusplus
     }
     #endif

   But C++ programs can also use OpenRSP by::

     extern "C" {
         #include "openrsp.h"
     }

   Someone also argues that the former solution makes a C code not a plain C
   code, and therefore prefers the latter solution, see
   `<http://stackoverflow.com/questions/16850992/call-a-c-function-from-c-code>`_.

   *FIXME: Therefore, what is the better choice for OpenRSP?*

