! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module matrix_defop

!> radovan: this file is generic for both DALTON and DIRAC
!>          DALTON or DIRAC specific things are in matrix_genop
!>          which is the backend to matrix_defop
!> 
!> Defined operators for type(matrix)
!> (which in linsca comes from module matrix_operations).
!>
!> The operators are as follows:
!
!>     equals                       A = B
!>     plus                         A + B
!>     minus                        A - B
!>     matrix product               A * B
!>     scale by real(8)             r * A
!>     divide by real(8)            A / r
!>     (conjugate) transpose        dag(A)
!>     trace                        tr(A)
!>     scalar product               dot(A,B)
!>     product trace                tr(A,B)
!>     short-hand freeing           A=0, A(:)=0, A(:,:)=0
!>     norm ie. sqrt(dot(A,A))      norm(A)
!>     mat_to_full                  A = B(:,:)
!>     mat_set_from_full            B(:,:) = A, B(:) = A
!>
!> To increase performance and lower memory needs, formulas are not
!> evaluated by executing the operators immediately, but the operators are
!> first collected into "batches" of the forms
!>
!>     axpby:     Z = a * X(^T) + b * Y   (where Z may be Y)
!>     gemm:      Y = a * X(^T) * Z(^T) + b * Y
!>
!> To create matrices, mat_init (from module matrix_operations) is replaced
!> by a more general init_mat. This is to hide context-specific information:
!> closed/open shell, real/complex, 1/2/4-component, full/sparse,
!> and thus make it possible to write general code:
!>
!>     call mat_init(A,nrow,ncol)       original version
!>     call mat_init(A,B)               same shape
!>     call mat_init(A,B,'T')           transpose shape
!>     call mat_init(A,B,'T',C,'N')     rows and columns from diff. mats.
!>     call mat_init(A,B,zero=.true.)   non-allocated zero of same shape
!>     call mat_init(A,reset=.true.)    reset to defaults, nullifying pointers
!>     call mat_init(A,B,alias='FF')    copy pointers, no new allocation
!>
!> Comparison of the shapes of matrices is done with mat_same
!>
!>     if (mat_same(A,B))               if same shape so that A+B possible
!>     if (mat_same(A,B,ta=.true.))     A^T + B is possible
!>     if (mat_same(A,B,mul=.true.))    if A's cols matches B's rows
!>
!> For printing matrices:
!>
!>     call mat_print(A)                default, print without dressing
!>
!>     call mat_print(A,label='A',width=7,unit=6,decor='{,},')
!>
!> The label is printed on the line preceding the matrix. Width specifies the
!> text width of each column (excluding separator). The number of digits printed
!> is adjusted so that no columns become overfilled (*****). Left brace,
!> column separator, right brace, row separator (decor) can be specified, so the
!> matrix can be copy-pasted as input to python, maple, mathematica, etc.
!>
!> To check the state of matrices: Whether defined: isdef(A),
!> whether zero: iszero(A)
!>
!> In order to pass a "formula", such as A+3*dag(B) etc., as argument X to a subroutine,
!> the subroutine must execute mat_fix_arg(X) before using X inside the subroutine.
!> At the end of the subroutine mat_unfix_arg(X) should be called in order to free X.
!> Neither call has any effect on a non-temporary matrix argument.
!> This is neccessary since the Fortran programmer cannot tell whether the temporary
!> is used inside its own subruoutine, where it should be deleted after use,
!> or inside a called subroutine, in which it should be non-temporary.
!>
!> In order to have type(matrix) as the result of a function, mat_fix_result(X)
!> should be called before returning from the function. This sets X's status as
!> temporary, so that it will be deleted once used.
module matrix_defop

  use matrix_genop, mat_print_nontemp => mat_print, &
                    mat_trace_nontemp => mat_trace

  implicit none

  public matrix
  public mat_init
  public mat_print
  public isdef
  public iszero
  public assignment(=)
  public operator(+)
  public operator(-)
  public operator(*)
  public trps
  public dot
  public tr
  public norm
  public matrix_defop_debug
  private

  interface assignment(=)
     module procedure mat_eq_mat
     module procedure mat_eq_zero
     module procedure mat_eq_zero_1D
     module procedure mat_eq_zero_2D
     module procedure mat_eq_zero_3D
     module procedure mat_eq_zero_4D
  end interface

  interface isdef
     module procedure mat_isdef
  end interface

  interface iszero
     module procedure mat_iszero
  end interface

  interface operator(+)
     module procedure plus_mat
     module procedure mat_plus_mat
  end interface

  interface operator(-)
     module procedure minus_mat
     module procedure mat_minus_mat
  end interface

  interface operator(*)
     module procedure mat_times_mat
     module procedure complex_times_mat
     module procedure real_times_mat
  end interface

  interface trps
     module procedure mat_transpose
  end interface

  interface dot
     module procedure mat_dot_prod
  end interface

  interface tr
     module procedure mat_prod_trace
     module procedure mat_trace
  end interface

  interface norm
     module procedure mat_norm
  end interface


  !> for switching debugging on and off
  logical :: matrix_defop_debug = .false.

contains


  function mat_iszero(A)
    type(matrix), intent(in) :: A
    logical                  :: mat_iszero
    if (.not.mat_isdef(A)) &
       call quit('error: matrix isero(A), A undefined')
    mat_iszero = btest(A%flags, matf_zero)
  end function


  function must_eval(A, msg)
    type(matrix), intent(in) :: A
    character(*), intent(in) :: msg
    logical                  :: must_eval
    if (.not.mat_isdef(A)) call quit(msg)
    must_eval = btest(A%flags, matf_temp)
  end function


  !> evaluate the expression Y += fac * X^tx * Z^tz, as specified in
  !> the fields within Y, by executing either an axpy or a gemm operation.
  !> evaluate the matrix C (+)= fac * A^ta * B^tb, where '+' fac A ta B tb
  !> are speficied in C%flags C%temp_fac C%temp_X C%temp_Y
  subroutine eval_temp(C, overwr)
    type(matrix),      intent(inout) :: C
    logical, optional, intent(in)    :: overwr
    type(matrix), target  :: dupC
    logical :: plus, zeroC, useA
    ! plus decides whether to add to or overwrite C
    plus = .true.
    if (present(overwr)) plus = (.not.overwr)
    ! whether C is zero, ie. unallocated
    zeroC = btest(C%flags, matf_zero)
    ! whether result should go into A, because C is zero or an alias
    useA = ((zeroC .or. btest(C%flags, matf_alias)) &
            .and. associated(C%temp_X) &
            .and. .not.associated(C%temp_Y) &
            .and. (zeroC .or. C%temp_fac == 1) &
            .and. .not.btest(C%flags, matf_temp_tx))
    if (useA) useA = btest(C%temp_X%flags, matf_temp)
    ! pre-apply any scale factor on A, if result will go there
    if (useA .and. C%temp_fac /= 1) then
       call mat_axpy(C%temp_fac, C%temp_X, &
                     .false., .false., C%temp_X)
       C%temp_fac = 1
    end if
    ! swap A and C or reallocate C, if needed
    if (useA .and. zeroC) then
       call mat_dup(C%temp_X, C)
    else if (useA) then
       call mat_dup(C, dupC)
       call mat_dup(dupC%temp_X, C)
       C%temp_X => dupC%temp_X
       nullify(C%temp_Y)
       nullify(C%temp_X%temp_X)
       call mat_free(dupC, nodealloc=.true.)
       C%flags = ibset(C%flags, matf_temp)
       C%temp_X%flags = ibclr(C%temp_X%flags, matf_temp)
    else if (btest(C%flags, matf_alias)) then
       call mat_dup(C, dupC)
       call mat_free(C, nodealloc=.true.)
       call mat_init(C, dupC)
       C%temp_fac = dupC%temp_fac
       C%temp_X  => dupC%temp_X
       C%temp_Y  => dupC%temp_Y
       if (btest(dupC%flags, matf_temp_tx)) &
          C%flags = ibset(C%flags, matf_temp_tx)
       if (btest(dupC%flags, matf_temp_ty)) &
          C%flags = ibset(C%flags, matf_temp_ty)
       call mat_axpy((1d0,0d0), dupC, .false., .false., C)
    else if (zeroC) then
       dupC%flags     =  C%flags
       dupC%temp_fac =  C%temp_fac
       dupC%temp_X   => C%temp_X
       dupC%temp_Y   => C%temp_Y
       call mat_init(C, C)
       C%flags     =  dupC%flags
       C%temp_fac =  dupC%temp_fac
       C%temp_X   => dupC%temp_X
       C%temp_Y   => dupC%temp_Y
       C%flags = ibclr(C%flags, matf_zero)
       plus = .false.
    end if
    ! execute gemm or axpy, delete temporary A and B
    if (associated(C%temp_X) .and. associated(C%temp_Y)) then
       call mat_gemm(C%temp_fac, C%temp_X, &
                     btest(C%flags, matf_temp_tx), C%temp_Y, &
                     btest(C%flags, matf_temp_ty), plus, C)
       if (btest(C%temp_X%flags, matf_temp)) C%temp_X = 0
       if (btest(C%temp_Y%flags, matf_temp)) C%temp_Y = 0
    else if (associated(C%temp_X)) then
       call mat_axpy(C%temp_fac, C%temp_X, &
                     btest(C%flags, matf_temp_tx), plus, C)
       if (btest(C%temp_X%flags, matf_temp)) C%temp_X = 0
    end if
    ! clean-up
    C%flags = ibset(C%flags, matf_temp)
    C%flags = ibclr(C%flags, matf_zero)
    C%flags = ibclr(C%flags, matf_temp_tx)
    C%flags = ibclr(C%flags, matf_temp_ty)
    C%temp_fac = 1
    nullify(C%temp_X)
    nullify(C%temp_Y)
  end subroutine



  !> transfer B to A
  subroutine mat_eq_mat(A, B)
    type(matrix), intent(inout) :: A
    type(matrix)                :: B
    logical :: evalB, haveA, useA, overwr
    ! verify that B is defined, and check whether it needs to be evaluated
    evalB = must_eval(B, 'error: matrix A=B, B undefined')
    ! is A currently allocated
    haveA = (isdef(A) .and. .not.btest(A%flags, matf_zero) &
                      .and. .not.btest(A%flags, matf_alias))
    ! if B is an alias of A, make B the original and A the alias
    if (haveA .and. btest(B%flags, matf_alias)) then
       if (mat_isdup(A,B)) then
          A%flags = ibset(A%flags, matf_alias)
          B%flags = ibclr(B%flags, matf_alias)
          haveA = .false.
       end if
    end if
    ! otherwise, if X is A, make it temporary
    if (haveA .and. associated(B%temp_X)) then
       if (mat_isdup(A, B%temp_X)) then
          if (btest(B%temp_X%flags, matf_alias)) &
             call quit('mat_eq_mat error: attempt to overwrite aliased matrix')
          B%temp_X%flags = ibset(B%temp_X%flags, matf_temp)
          haveA = .false.
       end if
    end if
    ! otherwise, if Y is A, make it temporary
    if (haveA .and. associated(B%temp_Y)) then
       if (mat_isdup(A, B%temp_Y)) then
          if (btest(B%temp_Y%flags, matf_alias)) &
             call quit('mat_eq_mat error: attempt to overwrite aliased matrix')
          B%temp_Y%flags = ibset(B%temp_Y%flags, matf_temp)
          haveA = .false.
       end if
    end if
    ! determine whether As contents are needed
    useA = haveA
    if (useA .and. btest(B%flags, matf_temp)) &
       useA = .not.(btest(B%flags, matf_temp) .and. &
                    .not.btest(B%flags, matf_alias))
    if (useA .and. associated(B%temp_X)) &
       useA = .not.(btest(B%temp_X%flags, matf_temp) .and. &
                    .not.btest(B%flags, matf_temp_tx) .and. &
                    .not.associated(B%temp_Y))
    if (useA) useA = mat_same(A, B)
    ! A = B=0
    if (.not.evalB .and. btest(B%flags, matf_zero)) then
       if (haveA) call mat_free(A)
       call mat_dup(B, A)
    ! A = B
    else if (.not.evalB) then
       if (haveA .and. .not.useA) call mat_free(A)
       if (.not.useA) call mat_init(A, B)
       call mat_axpy((1d0,0d0), B, .false., .false., A)
    ! A = B + X*Y, fits A, and B (and X without Y) non-temp
    else if (useA) then
       overwr = btest(B%flags, matf_zero)
       if (.not.overwr) &
          call mat_axpy((1d0,0d0), B, .false., .false., A)
       A%temp_fac = B%temp_fac
       A%temp_X  => B%temp_X
       A%temp_Y  => B%temp_Y
       if (btest(B%flags, matf_temp_tx)) &
          A%flags = ibset(A%flags, matf_temp_tx)
       if (btest(B%flags, matf_temp_ty)) &
          A%flags = ibset(A%flags, matf_temp_ty)
       call mat_free(B, nodealloc=.true.)
       call eval_temp(A, overwr)
    ! A = B + X*Y, either of B or (X without Y and tx) are temporary
    else
       if (haveA) call mat_free(A)
       call eval_temp(B)
       call mat_dup(B, A)
       call mat_free(B, nodealloc=.true.)
    end if
    ! ensure A is non-temp
    A%flags = ibclr(A%flags, matf_temp)
  end subroutine



  subroutine gen_mat_plus_mat(A, plus, B, C)
    type(matrix), target      :: A, B
    type(matrix), intent(out) :: C
    logical,      intent(in)  :: plus
    logical :: evalA, evalB, zeroA, zeroB
    evalA = must_eval(A, 'error: matrix A+B, A undefined')
    evalB = must_eval(B, 'error: matrix A+B, B undefined')
    zeroA = (btest(A%flags, matf_zero) .and. .not.evalA)
    zeroB = (btest(B%flags, matf_zero) .and. .not.evalB)
    if (.not.mat_same(A,B)) &
       call quit('error: matrix A+B, different shapes')
    if (evalA) call eval_temp(A)
    if (evalB) call eval_temp(B)
    call mat_dup(A, C)
    C%flags = ibset(C%flags, matf_temp)
    if (.not.zeroA) C%flags = ibset(C%flags, matf_alias)
    if (.not.zeroB) then
       C%temp_fac = merge(1,-1,plus)
       C%temp_X => B
    end if
  end subroutine


  !> verify that shapes match, check for zeros
  function mat_plus_mat(A, B) result(C)
    type(matrix), target, intent(in) :: A, B
    type(matrix)                     :: C
    call gen_mat_plus_mat(A, .true., B, C)
  end function


  !> verify that shapes match, check for zeros
  function mat_minus_mat(A, B) result(C)
    type(matrix), target, intent(in) :: A, B
    type(matrix)                     :: C
    call gen_mat_plus_mat(A, .false., B, C)
  end function



  !> verify that shapes match, check for zeroes
  function mat_times_mat(A, B) result(C)
    type(matrix), target :: A, B
    type(matrix)         :: C
    logical              :: evalA, evalB, zeroA, zeroB
    evalA = must_eval(A, 'error: matrix A*B, A undefined')
    evalB = must_eval(B, 'error: matrix A*B, B undefined')
    zeroA = (btest(A%flags, matf_zero) .and. .not.btest(A%flags, matf_temp))
    zeroB = (btest(B%flags, matf_zero) .and. .not.btest(B%flags, matf_temp))
    if (.not.mat_same(A, B, mul=.true.)) &
       call quit('error: matrix A*B, but A%ncol and B%nrow do not match')
    call mat_init(C, A, B, noalloc=.true.)
    C%flags = ibset(C%flags, matf_zero)
    C%flags = ibset(C%flags, matf_temp)
    if (evalA) call eval_temp(A)
    if (evalB) call eval_temp(B)
    if (zeroA .or. zeroB) then
       if (btest(A%flags, matf_temp)) call mat_free(A)
       if (btest(B%flags, matf_temp)) call mat_free(B)
       C%flags = ibset(C%flags, matf_zero)
    else
       C%temp_fac =  1
       C%temp_X   => A
       C%temp_Y   => B
    end if
  end function



  subroutine gen_complex_times_mat(r, A, B)
    complex(8),   intent(in) :: r
    type(matrix), target     :: A
    type(matrix)             :: B
    logical :: eval
    eval = must_eval(A, 'error: matrix r*A, A undefined')
    call mat_init(B, A, noalloc=.true.)
    if (eval) call eval_temp(A)
    if (r==0 .or. iszero(A)) then
       B%flags = ibset(B%flags, matf_zero)
    else
       B%flags = ibset(B%flags, matf_temp)
       B%temp_fac =  r
       B%temp_X   => A
    end if
    if (btest(A%flags, matf_temp)) &
       call mat_free(A, nodealloc = (btest(A%flags, matf_temp) &
                                .or. btest(A%flags, matf_alias)))
  end subroutine


  function complex_times_mat(r, A) result(B)
     complex(8),           intent(in) :: r
     type(matrix), target, intent(in) :: A
     type(matrix), target             :: B
     call gen_complex_times_mat(r, A, B)
  end function


  function real_times_mat(r, A) result(B)
    real(8),              intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(matrix), target             :: B
    call gen_complex_times_mat(r*(1d0,0d0), A, B)
  end function


  function plus_mat(A) result(B)
    type(matrix), target, intent(in) :: A
    type(matrix), target             :: B
    call gen_complex_times_mat((1d0,0d0), A, B)
  end function


  function minus_mat(A) result(B)
    type(matrix), target, intent(in) :: A
    type(matrix), target             :: B
    call gen_complex_times_mat((-1d0,0d0), A, B)
  end function


  !> Transpose the matrix A. If A is temporary and A=a*X or A=a*X*Z,
  !> only swap X and Z, and slip-swap tx and tz.
  !> Otherwise, evaluate A (if tmp), and create a new temporary for B
  function mat_transpose(A) result(B)
    type(matrix), target :: A
    type(matrix)         :: B
    logical :: eval
    eval = must_eval(A, 'error: matrix transpose trps(A), A undefined')
    call mat_init(B, A, ta=.true., noalloc=.true.)
    B%flags = ibset(B%flags, matf_zero)
    if (.not.btest(A%flags, matf_zero)) then
       if (eval) call eval_temp(A)
       B%flags = ibset(B%flags, matf_temp)
       B%flags = ibset(B%flags, matf_temp_tx)
       B%temp_X => A
    end if
  end function



  !> back-end for mat_dot_prod and mat_prod_trace (dot(A,B) and tr(A,B))
  function gen_mat_dot(A, B, trps)
    type(matrix) :: A, B
    complex(8)   :: gen_mat_dot
    logical      :: trps, zeroA, zeroB, evalA, evalB
    evalA = must_eval(A, 'error: matrix dot/tr(A,B), A undefined')
    evalB = must_eval(B, 'error: matrix dot/tr(A,B), B undefined')
    zeroA = (.not.evalA .and. btest(A%flags, matf_zero))
    zeroB = (.not.evalB .and. btest(B%flags, matf_zero))
    if (evalA) then
       !ajt FIXME forward A = fac * X^tx
       call eval_temp(A)
    end if
    if (evalB) then
       !ajt FIXME forward B = fac * X^tx
       call eval_temp(B)
    end if
    ! verify that shapes match
    if (.not.mat_same(A, B, ta=trps)) &
       call quit('error: matrix dot/tr(A,B), A and B have different shapes')
    ! don't calculate if
    if (      zeroA .or. zeroB)  gen_mat_dot = 0
    if (.not.(zeroA .or. zeroB)) gen_mat_dot = mat_dot(A, B, trps)
    ! delete any temporaries
    if (evalA) call mat_free(A)
    if (evalB) call mat_free(B)
  end function



  function mat_dot_prod(A, B)
    type(matrix), target, intent(in) :: A, B
    complex(8)                       :: mat_dot_prod
    mat_dot_prod = gen_mat_dot(A, B, (.false.))
  end function


  function mat_prod_trace(A, B)
    type(matrix) :: A, B
    complex(8)   :: mat_prod_trace
    mat_prod_trace = gen_mat_dot(A, B, (.true.))
  end function


  function mat_trace(A)
    type(matrix) :: A
    complex(8)   :: mat_trace
    logical :: temp, zero
    temp = must_eval(A, 'error: tr(A), A undefined')
    zero = (btest(A%flags, matf_zero) .and. .not.temp)
    if (temp) call eval_temp(A)
    if (.not.zero) mat_trace = mat_trace_nontemp(A)
    if (   zero  ) mat_trace = 0
    if (temp) call mat_free(A)
  end function


  function mat_norm(A)
    type(matrix) :: A
    real(8)      :: mat_norm
    complex(8)   :: norm2
    logical      :: eval, temp, zero
    eval = must_eval(A, 'error: norm(A), A undefined')
    temp = btest(A%flags, matf_temp)
    zero = (btest(A%flags, matf_zero) .and. .not.eval)
    if (temp) call eval_temp(A)
    if (temp) A%flags = ibclr(A%flags, matf_temp)
    norm2 = gen_mat_dot(A, A, (.false.))
    if (temp) A = 0
    mat_norm = sqrt(dreal(norm2))
  end function


  !> Short-hand A=0 frees matrix A if it is defined, otherwise nothing
  subroutine mat_eq_zero(A, zero)
    type(matrix), intent(inout) :: A
    integer,      intent(in)    :: zero
    if (zero/=0) call quit('error: matrix A = z: z must be 0')
    if (isdef(A)) call mat_free(A)
  end subroutine


  subroutine mat_eq_zero_1D(A, zero)
    type(matrix), intent(inout) :: A(:)
    integer,      intent(in)    :: zero
    integer                     :: i
    do i = 1, size(A)
       A(i) = 0
    end do
  end subroutine


  subroutine mat_eq_zero_2D(A, zero)
    type(matrix), intent(inout) :: A(:,:)
    integer,      intent(in)    :: zero
    integer                     :: i,j
    do j = 1, size(A,2)
       do i = 1, size(A,1)
          A(i,j) = 0
       end do
    end do
  end subroutine


  subroutine mat_eq_zero_3D(A, zero)
    type(matrix), intent(inout) :: A(:,:,:)
    integer,      intent(in)    :: zero
    integer                     :: i,j,k
    do k = 1, size(A,3)
       do j = 1, size(A,2)
          do i = 1, size(A,1)
             A(i,j,k) = 0
          end do
       end do
    end do
  end subroutine


  subroutine mat_eq_zero_4D(A, zero)
    type(matrix), intent(inout) :: A(:,:,:,:)
    integer,      intent(in)    :: zero
    integer                     :: i,j,k,l
    do l = 1, size(A,4)
       do k = 1, size(A,3)
          do j = 1, size(A,2)
             do i = 1, size(A,1)
                A(i,j,k,l) = 0
             end do
          end do
       end do
    end do
  end subroutine


  !> wrap matrix_defop's mat_print so that temporary and
  !> zero matrices can be printed
  subroutine mat_print(A, label, unit, width, decor)
    type(matrix)                       :: A
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: decor
    logical :: eval, zero, temp
    eval = must_eval(A, 'error: mat_print(A), A undefined')
    temp = btest(A%flags, matf_temp)
    zero = (btest(A%flags, matf_zero) .and. .not.temp)
    if (eval) call eval_temp(A)
    if (zero) call mat_init(A, A, zero=.true.)
    call mat_print_nontemp(A, label, unit, width, decor)
    if (eval .or. temp .or. zero) call mat_free(A)
  end subroutine


end module


subroutine quit(msg)
  character(*), intent(in) :: msg
  print *, msg
  stop
end subroutine



program test

  use matrix_defop
  use matrix_genop, only: matf_clsh !closed-shell
  implicit none

  !------ hard-code some matrix elements for testing -------
  ! This is the converged orbital coefficient matrix for H2O2 in STO-3G basis
  ! with Hartree-Fock (stored transposed)
  real(8), parameter :: cmo_t_elms(9,12) = reshape((/ &
      0.7029125557d0,-0.7034218990d0,-0.1653576986d0,-0.1717881187d0,-0.0428946404d0, &
     -0.0537205739d0, 0.0461043661d0, 0.0481622768d0,-0.0184913160d0, &
      0.0187240325d0,-0.0155704450d0, 0.5921458033d0, 0.6547115887d0, 0.2330299347d0, &
      0.2693265593d0,-0.2152610765d0,-0.2472100367d0, 0.0996244165d0, &
      0.0009619290d0,-0.0009362593d0, 0.0208000696d0, 0.0373913674d0,-0.0632496650d0, &
     -0.5567077025d0,-0.0724179656d0,-0.3505246009d0,-0.6965714321d0, &
     -0.0025013368d0,-0.0008575570d0,-0.0819785515d0, 0.0797872902d0,-0.2052800653d0, &
      0.0539144898d0,-0.5871997373d0,-0.0542563289d0, 0.1127518318d0, &
      0.0022677548d0,-0.0023751718d0, 0.0691252718d0, 0.0884842073d0,-0.4321676646d0, &
     -0.1652464276d0, 0.2578300957d0, 0.5202903321d0, 0.1464654704d0, &
     -0.7029125558d0,-0.7034218990d0,-0.1653576986d0, 0.1717881187d0,-0.0428946404d0, &
      0.0537205739d0, 0.0461043661d0,-0.0481622768d0,-0.0184913160d0, &
     -0.0187240325d0,-0.0155704450d0, 0.5921458033d0,-0.6547115887d0, 0.2330299347d0, &
     -0.2693265593d0,-0.2152610765d0, 0.2472100367d0, 0.0996244165d0, &
      0.0009619290d0, 0.0009362593d0,-0.0208000696d0, 0.0373913674d0, 0.0632496650d0, &
     -0.5567077025d0, 0.0724179656d0,-0.3505246009d0, 0.6965714321d0, &
     -0.0025013368d0, 0.0008575570d0, 0.0819785515d0, 0.0797872902d0, 0.2052800653d0, &
      0.0539144898d0, 0.5871997373d0,-0.0542563289d0,-0.1127518318d0, &
     -0.0022677548d0,-0.0023751718d0, 0.0691252718d0,-0.0884842073d0,-0.4321676646d0, &
      0.1652464276d0, 0.2578300957d0,-0.5202903321d0, 0.1464654704d0, &
     -0.0043805243d0, 0.0044574988d0, 0.1024238890d0, 0.1734654106d0,-0.3425350997d0, &
     -0.2576428881d0, 0.0000223994d0, 0.2601429825d0,-0.0885684546d0, &
      0.0043805243d0, 0.0044574988d0, 0.1024238890d0,-0.1734654106d0,-0.3425350997d0, &
      0.2576428881d0, 0.0000223994d0,-0.2601429825d0,-0.0885684546d0/), (/9,12/))
  real(8), parameter :: ovl_elms(12,12) = reshape((/ &
      1.0000000000d0, 0.2367039365d0, 0.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0082403834d0, 0.0000000000d0, &
      0.0156536771d0, 0.0000000000d0, 0.0523613022d0, 0.0061118010d0, &
      0.2367039365d0, 1.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 0.0082403834d0, 0.1528147927d0, 0.0000000000d0, &
      0.1977725882d0, 0.0000000000d0, 0.4664441507d0, 0.0904222679d0, &
      0.0000000000d0, 0.0000000000d0, 1.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 0.0726585685d0, &
      0.0000000000d0, 0.0000000000d0, 0.1465895707d0,-0.0187884791d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 1.0000000000d0, &
      0.0000000000d0,-0.0156536771d0,-0.1977725882d0, 0.0000000000d0, &
     -0.2247513423d0, 0.0000000000d0, 0.0675432857d0,-0.0852370092d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      1.0000000000d0, 0.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 0.0726585685d0, 0.3538985305d0, 0.0453594011d0, &
      0.0000000000d0, 0.0082403834d0, 0.0000000000d0,-0.0156536771d0, &
      0.0000000000d0, 1.0000000000d0, 0.2367039365d0, 0.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0061118010d0, 0.0523613022d0, &
      0.0082403834d0, 0.1528147927d0, 0.0000000000d0,-0.1977725882d0, &
      0.0000000000d0, 0.2367039365d0, 1.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0904222679d0, 0.4664441507d0, &
      0.0000000000d0, 0.0000000000d0, 0.0726585685d0, 0.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 1.0000000000d0, &
      0.0000000000d0, 0.0000000000d0, 0.0187884791d0,-0.1465895707d0, &
      0.0156536771d0, 0.1977725882d0, 0.0000000000d0,-0.2247513423d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      1.0000000000d0, 0.0000000000d0, 0.0852370092d0,-0.0675432857d0, &
      0.0000000000d0, 0.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      0.0726585685d0, 0.0000000000d0, 0.0000000000d0, 0.0000000000d0, &
      0.0000000000d0, 1.0000000000d0, 0.0453594011d0, 0.3538985305d0, &
      0.0523613022d0, 0.4664441507d0, 0.1465895707d0, 0.0675432857d0, &
      0.3538985305d0, 0.0061118010d0, 0.0904222679d0, 0.0187884791d0, &
      0.0852370092d0, 0.0453594011d0, 1.0000000000d0, 0.1256680400d0, &
      0.0061118010d0, 0.0904222679d0,-0.0187884791d0,-0.0852370092d0, &
      0.0453594011d0, 0.0523613022d0, 0.4664441507d0,-0.1465895707d0, &
     -0.0675432857d0, 0.3538985305d0, 0.1256680400d0, 1.0000000000d0/), (/12,12/))
  ! matrices
  type(matrix) :: C, S, D, A, B

  ! manually initialize orbital matrix C
  C%nrow  = size(cmo_t_elms,2)
  C%ncol  = size(cmo_t_elms,1)
  C%flags = ibset(0, matf_clsh) !closed shell
  call mat_init(C, C)
  C%elms = transpose(cmo_t_elms)
  ! print in 'python/C++' format
  call mat_print(C, label='orbital', decor='[,],')

  ! manually initialize overlap matrix S
  S%nrow  = size(ovl_elms,1)
  S%ncol  = size(ovl_elms,2)
  S%flags = 0 !normal matrix
  call mat_init(S, S)
  S%elms = ovl_elms
  ! print in plain format
  call mat_print(S, label='overlap')

  ! run tests
  call check_orthonormality_of_C
  call calculate_density_D
  call count_electrons_in_D
  call check_idempotency_of_D

  ! free matrices
  C=0; S=0; D=0; A=0; B=0

contains

  subroutine check_orthonormality_of_C
    ! calculate norm of C^T S C minus identity matrix
    type(matrix) :: id
    integer      :: i
    !------- manually initialize identity matrix
    id%nrow  = size(cmo_t_elms,1)
    id%ncol  = id%nrow
    id%flags = 0
    call mat_init(id,id)
    id%elms(:,:) = 0
    do i = 1, id%nrow
       id%elms(i,i) = 1
    end do

    ! ...

    ! free
    id=0
  end subroutine


  subroutine calculate_density_D
    ! D = C C^T

    D = C * trps(C)! ...

  end subroutine


  subroutine count_electrons_in_D
    ! rewrite Tr C^T S C in terms of D

    print *, 'trSD =', dreal(tr(S,D))! ...

  end subroutine


  subroutine check_idempotency_of_D
    ! idempotency relation is D S D = D

    ! ...

  end subroutine


end program
