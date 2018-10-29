.. _chapter_openrsp_design:

Design OpenRSP
==============

.. _section_literate_programming:

Literate Programming
--------------------

Currently, we use literate programming for OpenRSP APIs. That is, source codes
in ``include`` and ``src`` are genereated from the WEB files in ``web``, as
well as developer documentation generated also from these WEB files.

However, lots of extra work on maintenance and further development is required
for literate programming. We are considering to change or abandon literate
programming in OpenRSP.

.. _section_symbolic_computations:

Symbolic Computations (not implemented)
---------------------------------------

The recursive programming techniques described in Ref. [Ringholm2014_] make it
possible to calculate molecular properties of arbitrary complexity in an
analytical manner. But if we are going to implement MO and CC based response
theories, it may mess up or not be easy to (re)use the already developed codes.

Therefore, I (Bin Gao) propose a clear separation of symbolic and numerical
computations in OpenRSP, and to develop a set of symbolic computation functions
that can be used by response theory codes to get molecular properties still in
a recursive and analytical manner.

.. _section_perturbation_free:

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

* :math:`x+x\Rightarrow xx`, :math:`0+0\Rightarrow 0`,
* :math:`x+y\Rightarrow xy`, :math:`0+1\Rightarrow 1`,
* :math:`x+z\Rightarrow xz`, :math:`0+2\Rightarrow 2`,
* :math:`y+x\Rightarrow yx`, :math:`1+0\Rightarrow 3`,
* :math:`y+y\Rightarrow yy`, :math:`1+1\Rightarrow 4`,
* :math:`y+z\Rightarrow yz`, :math:`1+2\Rightarrow 5`,
* :math:`z+x\Rightarrow zx`, :math:`2+0\Rightarrow 6`,
* :math:`z+y\Rightarrow zy`, :math:`2+1\Rightarrow 7`,
* :math:`z+z\Rightarrow zz`, :math:`2+2\Rightarrow 8`,

where we have ranked different components in zero-based numbering (numbers on
the right). However, the ranks can be different in different host programs. To
solve this problem, i.e., the mapping relationship of lower and higher order
derivatives with respect to **one perturbation**, we ask for a callback
function :c:func:`get_pert_concatenation` from host programs, which is the last
argument of the API :c:func:`OpenRSPSetPerturbations`.

One should note that, we emphasize the derivatives of **one perturbation**
here. Because components of higher order derivatives of different perturbations
are simply the direct product of components of lower order derivatives.

The callback function :c:func:`get_pert_concatenation` is used by OpenRSP to
get the ranks of components of *sub-perturbation tuples with same perturbation
label* (lower order derivatives with respect to one perturbation) for given
components of the corresponding *concatenated perturbation tuple* (higher order
derivatives).

.. _section_openrsp_Hamiltonian:

Hamiltonian in OpenRSP
----------------------

As aforementioned, the ingradients of Hamiltonian are sent to OpenRSP, and
(perturbed) integrals and expectation values will be computed by the callback
functions of host programs. These include:

#. overlap operator (source codes ``web/RSPOverlap.nw``),
#. one-electron operators (source codes ``web/RSPOneOper.nw``),
#. two-electron operators (source codes ``web/RSPTwoOper.nw``),
#. exchange-correlation functionals (source codes ``web/RSPXCFun.nw``),
#. zero-electron operators (source codes ``web/RSPZeroOper.nw``),

where the source codes save the callback functions as function pointers in
different C ``struct``, and take care the invoking of these callback functions
during calculations.

Different from the overlap operator, zero-, one- and two-electron operators and
XC functionals are saved in three different linked lists in OpenRSP, in which
each node corresponds to an operator. This makes it possible for host programs
to add different callback functions for different operators, if they do not
want to or can not provide OpenRSP a general callback function.

.. _section_respose_solver:

Response Equation Solver
------------------------

Similar to the overlap operator, the callback function of a linear response
equation solver is saved as a function pointer in a C ``struct`` in OpenRSP.
That will be invoked by OpenRSP for obtaining response parameters, and the
source codes related to the solver are in ``web/RSPSolver.nw``.

OpenRSP will send multiple RHS vectors (or matrices) to the solver, for several
frequency sums on the left hand side of the linear response equation and for
several derivatives with respect to (different) perturbations.

Notice that it would be more common and help the convergence to calculate
several frequencies for the same perturbation, than the other way around. So
the RHS matrices and response parameters are arranged as
``[num_comps][num_freq_sums]`` in the callback function
:c:func:`get_linear_rsp_solution`.

.. _section_fortran_api_impl:

Implementation of Fortran APIs
------------------------------

OpenRSP APIs that host programs will use to talk to OpenRSP are written in C
language, with Fortran support by using Fortran 2003 language. The source codes
are in ``web/FortranAPIs.nw``, and the framework of OpenRSP used in a Fortran
host program is shown in the following figure:

.. figure:: /_static/openrsp_fortran_api.png
   :scale: 100 %
   :align: center

Two new parts are needed for the use of OpenRSP in a Fortran program:

#. "OpenRSP Fortran APIs", and
#. "OpenRSP Fortran type".

Take one-electron operators as an example, the callback subroutine
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

But "OpenRSP C struct" codes can not call this subroutine
:c:func:`get_one_oper_mat` directly, because the ``type(QcMat)`` can not be
sent from a C function to a Fortran subroutine directly. Instead, another
"OpenRSP Fortran type" subroutine is implemented in OpenRSP that will be
called by the "OpenRSP C struct" codes::

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

"OpenRSP response" codes (Fortran) :math:`\Rightarrow` "OpenRSP C struct" codes
:math:`\Rightarrow` "OpenRSP Fortran type" subroutine ``RSPOneOperGetMat_f()``
:math:`\Rightarrow` :c:func:`get_one_oper_mat`

One can also notice that, the argument ``num_int`` is needed in the
``interface`` of ``OpenRSPAddOneOper_f()`` and the subroutine
``RSPOneOperGetMat_f()``, and "OpenRSP C struct" codes also need to pass
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

..   *FIXME: in OpenRSP Fortran APIs, we should choose complex numbers, right? Because Fortran support complex numbers.*

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
         #include "OpenRSP.h"
     }

   Someone also argues that the former solution makes a C code not a plain C
   code, and therefore prefers the latter solution, see
   `<http://stackoverflow.com/questions/16850992/call-a-c-function-from-c-code>`_.

   C++ users please decide the better choice by themselves.
