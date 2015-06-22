.. _chapter-openrsp-developers:

OpenRSP for Developers
======================

Theoretical Background
----------------------

Framework of OpenRSP
--------------------

One- and Two-Electron Operators
-------------------------------

Exchange-Correlation Functionals
--------------------------------

Response Equation Solver
------------------------

Implementation of Fortran APIs
------------------------------

OpenRSP APIs that host programs will use to talk to OpenRSP are written in C
language, with Fortran support by using Fortran 2003 language. Take
one-electron integrals as an example, I need to, in the OpenRSP API (Fortran
version) OpenRSPAddOneOper_f(), declare the callback subroutine
:py:meth:`get_one_oper_mat` in the interface::

    function OpenRSPAddOneOper_f(...)
        interface
            subroutine get_one_oper_mat(len_tuple,  &
                                        pert_tuple, &
                                        num_int,    &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_one_oper_mat
        end interface
    end function OpenRSPAddOneOper_f

But the C part of OpenRSP can not call this subroutine :py:meth:`get_one_oper_mat`
directly, because the type(QcMat) can not be sent from a C code to a Fortran
code directly. Instead, another subroutine is implemented in OpenRSP that will
be called by the C part of OpenRSP::

    subroutine RSPOneOperGetMat_f(len_tuple,  &
                                  pert_tuple, &
                                  num_int,    &
                                  val_int)    &
        bind(C, name="RSPOneOperGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: len_tuple
        integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(QcMat), allocatable :: f_val_int(:)
        integer(kind=4) ierr
        allocate(f_val_int(num_int), stat=ierr)
        ! converts val_int (type C_PTR) to f_val_int (type QcMat)
        ! call Fortran callback subroutine, which is saved in a Fortran type one_oper_fun as a procedure pointer
        call one_oper_fun%get_one_oper_mat(..., f_val_int)
    end subroutine RSPOneOperGetMat_f

The procedure when doing a callback becomes:

OpenRSP recursive part (Fortran) -> C part of OpenRSP that will make the callback -> RSPOneOperGetMat_f() -> get_one_oper_mat()

As you could notice, the argument num_int is needed in the interface of
OpenRSPAddOneOper_f() and the subroutine RSPOneOperGetMat_f(), and the C part
of OpenRSP also need to pass num_int to RSPOneOperGetMat_f() (from C to
Fortran). Therefore, these arguments for the dimension of arrays have to be
passed.
