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
!> Comparison of the shapes (and content) of matrices is done with mat_same
!>
!>     if (mat_same(A,B)) ..            if same shape and content (alias)
!>     if (mat_same(A,B,'N')) ..        if same shapes
!>     if (mat_same(A,B,'T')) ..        if transposed shapes
!>     if (mat_same(A,B,'T','R')) ..    if A's rows matches B's columns
!>     if (mat_same(A,B,'N','C')) ..    if A's columns matches B's columns
!>
!> For printing matrices, mat_print has been replaced:
!>
!>     call print_mat(A,label='A',colwidth=7,unit=6,braces='{')
!>
!> The label is printed on the line preceding the matrix. Width specifies the
!> text width of each column (excluding separator). The number of digits printed
!> is adjusted so that no columns become overfilled (*****). Braces can either
!> be one character {, [, (, | etc., or four '{,};', where the comma is column
!> separator and semicolon row separator. This is convenient for pasting into
!> python, maple, mathematica, etc.
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

  use matrix_genop, mat_print_nontemp => mat_print

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
  !> are speficied in C%flags C%defop_fac C%defop_A C%defop_B
  subroutine eval_temp(C, noplus)
    type(matrix),      intent(inout) :: C
    logical, optional, intent(in)    :: noplus
    type(matrix), target  :: dupC
    logical :: plus, zeroC, useA
print *, 'i eval_temp', associated(C%defop_A), associated(C%defop_B)
    ! plus decides whether to add to or overwrite C
    plus = .true.
    if (present(noplus)) plus = (.not.noplus)
    ! whether C is zero, ie. unallocated
    zeroC = btest(C%flags, matf_zero)
    ! whether result should go into A, because C is zero or an alias
    useA = ((zeroC .or. btest(C%flags, matf_alias)) &
            .and. associated(C%defop_A) &
            .and. .not.associated(C%defop_B) &
            .and. (zeroC .or. C%defop_fac == 1) &
            .and. .not.btest(C%flags, matf_ta))
    if (useA) useA = btest(C%defop_A%flags, matf_temp)
    ! pre-apply any scale factor on A, if result will go there
print *, 'kaller mat_gemm?', associated(C%defop_A), associated(C%defop_B)
    if (useA .and. C%defop_fac /= 1) then
       call mat_axpy(C%defop_fac, C%defop_A, &
                     .false., .false., C%defop_A)
       C%defop_fac = 1
    end if
print *, 'kaller mat_gemm?', associated(C%defop_A), associated(C%defop_B), useA, zeroC
    ! swap A and C or reallocate C, if needed
    if (useA .and. zeroC) then
       call mat_dup(C%defop_A, C)
    else if (useA) then
       call mat_dup(C, dupC)
       call mat_dup(dupC%defop_A, C)
       C%defop_A => dupC%defop_A
       nullify(C%defop_A%defop_A)
       call mat_free(dupC, nodealloc=.true.)
       C%flags = ibset(C%flags, matf_temp)
    else if (btest(C%flags, matf_alias)) then
       call mat_dup(C, dupC)
       call mat_free(C, nodealloc=.true.)
       call mat_init(C, dupC)
       call mat_axpy((1d0,0d0), dupC, .false., .false., C)
    else if (zeroC) then
       dupC%defop_fac =  C%defop_fac
       dupC%defop_A   => C%defop_A
       dupC%defop_B   => C%defop_B
       call mat_init(C, C)
       C%defop_fac =  dupC%defop_fac
       C%defop_A   => dupC%defop_A
       C%defop_B   => dupC%defop_B
       C%flags = ibclr(C%flags, matf_zero)
       plus = .false.
    end if
    ! execute gemm or axpy, delete temporary A and B
print *, 'kaller mat_gemm?', associated(C%defop_A), associated(C%defop_B)
    if (associated(C%defop_A) .and. associated(C%defop_B)) then
print *, 'kaller mat_gemm', C%defop_fac
       call mat_gemm(C%defop_fac, C%defop_A, &
                     btest(C%flags, matf_ta), C%defop_B, &
                     btest(C%flags, matf_tb), plus, C)
       if (btest(C%defop_A%flags, matf_temp)) C%defop_A = 0
       if (btest(C%defop_B%flags, matf_temp)) C%defop_B = 0
print *, 'etter mat_gemm'
    else if (associated(C%defop_A)) then
print *, 'kaller axpy'
       call mat_axpy(C%defop_fac, C%defop_A, &
                     btest(C%flags, matf_ta), plus, C)
print *, 'etter axpy'
       if (btest(C%defop_A%flags, matf_temp)) C%defop_A = 0
    end if
    ! clean-up
    C%flags = ibset(C%flags, matf_temp)
    C%flags = ibclr(C%flags, matf_ta)
    C%flags = ibclr(C%flags, matf_tb)
    C%defop_fac = 1
    nullify(C%defop_A)
    nullify(C%defop_B)
  end subroutine



  !> transfer B to A
  subroutine mat_eq_mat(A, B)
    type(matrix), intent(inout) :: A
    type(matrix)                :: B
    logical :: defA, evalB, ali, fits, noplus
    defA = isdef(A)
    evalB = must_eval(B, 'error: matrix A=B, B undefined')
    ! if A not alias, and either of B B%defop_A and B%defop_B is
    ! an alias of A, make A an alias of that instead
    if (defA) then
       ali = mat_isdup(A, B)
       if (ali) B%flags = ibclr(B%flags, matf_alias)
       if (.not.ali .and.associated(B%defop_A)) then
          ali = mat_isdup(A, B%defop_A)
          if (ali) B%defop_A%flags = ibclr(B%defop_A%flags, matf_alias)
       end if
       if (.not.ali .and.associated(B%defop_B)) then
          ali = mat_isdup(A, B%defop_B)
          if (ali) B%defop_B%flags = ibclr(B%defop_B%flags, matf_alias)
       end if
       if (ali) A%flags = ibset(A%flags, matf_alias)
    end if
    ! determine whether the current allocation of A is to be used
    fits = (.not.btest(B%flags, matf_zero) &
            .and. btest(B%flags, matf_temp) &
            .and. .not.btest(B%flags, matf_alias))
    if (.not.fits .and. associated(B%defop_A)) &
       fits = (.not.associated(B%defop_B) &
               .and. .not.btest(B%flags, matf_ta) &
               .and. btest(B%defop_A%flags, matf_temp) &
               .and.(btest(B%flags, matf_zero) .or. B%defop_fac == 1))
    if (.not.fits .and. defA) then
       fits = (.not.btest(A%flags, matf_alias) .and. &
               .not.btest(A%flags, matf_zero) .and. mat_same(A,B))
    else
       fits = .false.
    end if
    ! either put A into B or B into A
    noplus = .false.
    if (fits) then
       noplus = btest(B%flags, matf_zero)
       if (.not.noplus) &
          call mat_axpy((1d0,0d0), B, .false., .false., A)
       A%defop_fac = B%defop_fac
       A%defop_A  => B%defop_A
       A%defop_B  => B%defop_B
       if (btest(B%flags, matf_ta)) A%flags = ibset(A%flags, matf_ta)
       if (btest(B%flags, matf_tb)) A%flags = ibset(A%flags, matf_tb)
       call mat_free(B, nodealloc = .true.)
    else
       if (defA) call mat_free(A, nodealloc &
                             = (btest(A%flags, matf_alias) &
                           .or. btest(A%flags, matf_zero)))
       call mat_dup(B, A)
       call mat_free(B, nodealloc = .true.)
    end if
    ! evaluate, make non-temp
    call eval_temp(A, noplus)
    A%flags = ibclr(A%flags, matf_temp)
  end subroutine



  subroutine gen_mat_plus_mat(A, plus, B, C)
    type(matrix), target      :: A, B
    type(matrix), intent(out) :: C
    logical,      intent(in)  :: plus
    logical :: evalA, evalB, zeroA, zeroB
    evalA = must_eval(A, 'error: matrix A+B, A undefined')
    evalB = must_eval(B, 'error: matrix A+B, B undefined')
    zeroA = btest(A%flags, matf_zero)
    zeroB = btest(B%flags, matf_zero)
    if (.not.mat_same(A,B)) &
       call quit('error: matrix A+B, different shapes')
    call mat_dup(A, C)
    C%flags = ibset(C%flags, matf_temp)
    if (.not.zeroA) C%flags = ibset(C%flags, matf_alias)
    if (.not.zeroB) then
       C%defop_fac = merge(1,-1,plus)
       C%defop_A => B
    end if
    if (zeroA .and. btest(A%flags, matf_temp)) call mat_free(A)
    if (zeroB .and. btest(B%flags, matf_temp)) call mat_free(B)
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
    C%flags = ibset(C%flags, matf_temp)
    C%flags = ibset(C%flags, matf_zero)
    if (evalA) call eval_temp(A)
    if (evalB) call eval_temp(B)
    if (zeroA .or. zeroB) then
       if (btest(A%flags, matf_temp)) call mat_free(A)
       if (btest(B%flags, matf_temp)) call mat_free(B)
       C%flags = ibset(C%flags, matf_zero)
    else
       C%defop_fac =  1
       C%defop_A   => A
       C%defop_B   => B
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
       B%defop_fac =  r
       B%defop_A   => A
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
    if (btest(A%flags, matf_zero)) then
       B%flags = ibset(B%flags, matf_zero)
    else
       if (eval) call eval_temp(A)
       B%flags = ibset(B%flags, matf_temp)
       B%flags = ibset(B%flags, matf_zero)
       B%flags = ibset(B%flags, matf_ta)
       B%defop_A => A
    end if
  end function



  !> back-end for mat_dot_prod and mat_prod_trace (dot(A,B) and tr(A,B))
  function gen_mat_dot(A, B, trps)
    type(matrix) :: A, B
    complex(8)   :: gen_mat_dot
    logical      :: trps, zeroA, zeroB, evalA, evalB
    evalA = must_eval(A, 'error: matrix dot/tr(A,B), A undefined')
    evalB = must_eval(B, 'error: matrix dot/tr(A,B), B undefined')
    zeroA = btest(A%flags, matf_zero)
    zeroB = btest(B%flags, matf_zero)
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
    if (      zeroA.or.zeroB)  gen_mat_dot = 0
    if (.not.(zeroA.or.zeroB)) gen_mat_dot = mat_dot(A, B, trps)
    ! delete any temporaries
    if (evalA) call mat_free(A)
    if (evalB) call mat_free(B)
  end function



  function mat_dot_prod(A, B)
    type(matrix),target,intent(in) :: A, B
    complex(8)                     :: mat_dot_prod
    mat_dot_prod = gen_mat_dot(A, B, (.false.))
  end function


  function mat_prod_trace(A, B)
    type(matrix) :: A, B
    complex(8)   :: mat_prod_trace
    mat_prod_trace = gen_mat_dot(A, B, (.true.))
  end function


  function mat_norm(A)
    type(matrix) :: A
    real(8)      :: mat_norm
    complex(8)   :: norm2
    norm2 = gen_mat_dot(A, A, (.false.))
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
  subroutine mat_print(A, label, unit, width, braces)
    type(matrix)                       :: A
    character(*), optional, intent(in) :: label, braces
    integer,      optional, intent(in) :: unit, width
    logical :: eval, zero, temp
    eval = must_eval(A, 'error: print_mat(A), A undefined')
    temp = btest(A%flags, matf_temp)
    zero = btest(A%flags, matf_zero)
    if (eval) call eval_temp(A)
    if (zero) call mat_init(A, A, zero=.true.)
    call mat_print_nontemp(A, label, unit, width, braces)
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
  implicit none
  type(matrix) :: A, B
  real(8), target :: cmo_elms(9,12) = reshape((/ &
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
  A%nrow = size(cmo_elms,1)
  A%ncol = size(cmo_elms,2)
  call mat_init(A, A)
  A%elms = cmo_elms
  call mat_print(A)
  B = A * trps(A)
  call mat_print(B)
end program