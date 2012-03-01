! non-copyrighted template

!> @file
!> Contains module matrix_backend

module matrix_backend

  implicit none

  public matrix
  public mat_nullify
  public mat_setup
  public mat_alloc
  public mat_free
  public mat_coshape
  public mat_axpy
  public mat_gemm
  public mat_dot
  public mat_trace
  public mat_print
  public mat_duplicate
  public mat_magic_value
  public mat_magic_setup
  public matrix_backend_debug


  !> matrix structure
  type matrix
     !> number of rows and columns (of each block)
     integer nrow, ncol
     !> pointers to blocks of elements (nrow,ncol). elms: main block,
     !> elms_i: imaginary part, elms_b: beta part, elms_ib: imaginary of beta
     real(8), pointer :: elms(:,:), elms_i(:,:), &
                         elms_b(:,:), elms_ib(:,:)
     !> flags for complex elements elms_i, closed shell (*2 on mat_dot),
     !> open shell (beta elements elms_b)
     logical complex, closed_shell, open_shell
     !> tag used to spot accidental use of uninitialized and memory-corrupted
     !> matrices, and fail with something other than 'segmentation fault'.
     !> Set to mat_magic_setup by mat_setup, changed to mat_magic_value
     !> by mat_alloc, changed back to mat_setup_magic_falue by mat_free, and
     !> finally zeroed by mat_nullify. Can be checked every time a matrix is used
     integer magic_tag
     !> pointer to self, set by init to mark the matrix' correct location.
     !> Used to distinguish init'ed matrices from copy-in'ed matrices like
     !> "call subr((/A,B,C/))". Set by mat_alloc and nullified by mat_free
     !> and mat_nullify
     type(matrix), pointer :: self_pointer
  end type

  !> Large random number >-2^31, <+2^31 (but also works on 64bit),
  !> not likely to be found in memory by accident. Mat_alloc sets
  !> type(matrix)%magic_tag to this value, all matrix operations
  !> will expect it to be there, and mat_free nullifies it
  integer, parameter :: mat_magic_value = -1981879812

  !> Alternate value of magic_tag to be used between mat_setup and mat_alloc, 
  !> and between mat_free and mat_nullify. Set by mat_setup, changed by mat_alloc,
  !> changed back by mat_free, and nullified by mat_nullify
  integer, parameter :: mat_magic_setup = 825169837

  !> for switching debugging on and off
  logical :: matrix_backend_debug = .false.

  private

  contains

   subroutine quit(a)

   character(*) :: a

   end subroutine


  !> nullify all fields in matrix
  subroutine mat_nullify(A)
    type(matrix), intent(out) :: A
    A%nrow = huge(1)
    A%ncol = huge(1)
    nullify(A%elms)
    nullify(A%elms_i)
    nullify(A%elms_b)
    nullify(A%elms_ib)
    A%complex = .false.
    A%closed_shell = .false.
    A%open_shell = .false.
    A%magic_tag = 0
    nullify(A%self_pointer)
  end subroutine


  !> setup fields in matrix C so it inherits shape of:
  !> A, A^T, A*B, A^T*B, A*B^T or A^T*B^T. 
  !> ajt Note: This is where any point group symmetry would be implemented
  subroutine mat_setup(C, A, B, ta, tb)
    type(matrix),           intent(inout) :: C
    type(matrix),           intent(in)    :: A
    type(matrix), optional, intent(in)    :: B
    logical,      optional, intent(in)    :: ta, tb
    logical taa, tbb
    integer nrow
    if (present(tb) .and. .not.present(B)) &
       call quit('error: mat_setup(C,A,B,ta,tb), tb present without B')
    ! process optionals
    taa = .false. !transpose A
    if (present(ta)) taa = ta
    tbb = .false. !transpose B
    if (present(tb)) tbb = tb
    ! set nrow and ncol, keeping in mind that C may be A
    nrow = merge(A%ncol, A%nrow, taa) !don't overwrite C%nrow yet
    if (present(B)) then
       C%ncol = merge(B%nrow, B%ncol, tbb)
       C%closed_shell = (A%closed_shell .or. B%closed_shell)
       C%open_shell = (A%open_shell .or. B%open_shell)
       C%complex = (A%complex .or. B%complex)
    else
       C%ncol = merge(A%nrow, A%ncol, taa)
       C%closed_shell = A%closed_shell
       C%open_shell = A%open_shell
       C%complex = A%complex
    end if
    C%nrow = nrow
    ! cannot be both closed- and open-shell
    if (C%closed_shell .and. C%open_shell) &
       C%closed_shell = .false.
    ! finally, set magic_tag to mat_magic_setup, to indicate C is set up
    C%magic_tag = mat_magic_setup
  end subroutine


  !> allocate elms of previously set-up matrix
  subroutine mat_alloc(A)
    type(matrix), intent(inout), target :: A
if (matrix_backend_debug) print *, 'alloc'
    ! err if matrix has not been mat_setup'ed
    if (A%magic_tag /= mat_magic_setup) &
       call quit('error: mat_alloc(A), but A is not set up')
    ! allocate elms*
    allocate(A%elms(A%nrow,A%ncol))
    if (A%complex) &
       allocate(A%elms_i(A%nrow,A%ncol))
    if (A%open_shell) &
       allocate(A%elms_b(A%nrow,A%ncol))
    if (A%complex .and. A%open_shell) &
       allocate(A%elms_ib(A%nrow,A%ncol))
    ! set magic_tag
    A%magic_tag = mat_magic_value
    ! set self pointer
    A%self_pointer => A
  end subroutine


  !> free matrix's elements and nullify its fields
  subroutine mat_free(A)
    type(matrix), intent(inout), target :: A
if (matrix_backend_debug) print *, 'free'
    ! err if attempting to free un-allocated matrix
    if (A%magic_tag /= mat_magic_value) &
       call quit('error: mat_free(A), A not allocated')
    ! err if attempting to free a relocated (copy-in'ed) matrix
    if (.not.associated(A%self_pointer,A)) &
       call quit('error: mat_free(A), A has been relocated')
    ! deallocate
    deallocate(A%elms)
    if (A%complex) &
       deallocate(A%elms_i)
    if (A%open_shell) &
       deallocate(A%elms_b)
    if (A%complex .and. A%open_shell) &
       deallocate(A%elms_ib)
    ! nullify
    call mat_nullify(A)
  end subroutine


  !> Compare shapes of A and B, and see whether one of the following
  !> are possible: A=B, A=B^T, A+B, A+B^T, A*B, A^T*B, or A*B^T
  function mat_coshape(A, B, ta, tb, add, mul)
    type(matrix),      intent(in) :: A, B
    logical, optional, intent(in) :: ta, tb, add, mul
    logical taa, tbb, aadd, mmul, mat_coshape
    taa  = .false.
    if (present(ta))  taa  = ta
    tbb  = .false.
    if (present(tb))  tbb  = tb
    aadd = .false.
    if (present(add)) aadd = add
    mmul = .false.
    if (present(mul)) mmul = mul
    ! compare shapes
    if (mmul) then
       mat_coshape = (merge(A%nrow, A%ncol, taa) &
                   == merge(B%ncol, B%nrow, tbb))
    else !+ or =
       mat_coshape = (merge(A%ncol, A%nrow, taa) &
                   == merge(B%ncol, B%nrow, tbb) &
                .and. merge(A%nrow, A%ncol, taa) &
                   == merge(B%nrow, B%ncol, tbb))
    end if
    ! compare data types (=)
    if (.not.aadd .and. .not.mmul) &
       mat_coshape = (mat_coshape .and. &
             (A%complex      .eqv. B%complex) .and. &
             (A%closed_shell .eqv. B%closed_shell) .and. &
             (A%open_shell   .eqv. B%open_shell))
  end function


  !> Matrix axpy: Y+=a*X (tx=F zy=F), Y=a*X (tx=F zy=T),
  !> Y+=a*X^T (tx=T zy=F) or Y=a*X^T (tx=T zy=T)
  subroutine mat_axpy(a, X, tx, zy, Y)
    complex(8),   intent(in)    :: a
    type(matrix), intent(in)    :: X
    logical,      intent(in)    :: tx, zy !transpose X, zero Y
    type(matrix), intent(inout) :: Y
if (matrix_backend_debug) &
print '(a,x,f4.1,x,l,x,l)', ' axpy a tx zy =', dreal(a), tx, zy
    if (a==0 .and. zy) then
       Y%elms = 0 !zeroing
    else if (a==0) then
       ! nothing
    else if (associated(X%elms,Y%elms)) then
       if (tx .and. zy) then
          Y%elms = a * transpose(Y%elms)
       else if (tx) then
          Y%elms = Y%elms + a * transpose(Y%elms)
       else if (zy) then
          Y%elms = a * Y%elms
       else
          Y%elms = (1+a) * Y%elms
       end if
    else if (tx .and. zy) then
       Y%elms = a * transpose(X%elms)
    else if (tx) then
       Y%elms = Y%elms + a * transpose(X%elms)
    else if (zy) then
       Y%elms = a * X%elms
    else
       Y%elms = Y%elms + a * X%elms
    end if
  end subroutine



  !> matrix multiply: C (+)= f * A^ta * B^tb
  !> assumed initialied and with corresponding shapes
  subroutine mat_gemm(f, A, ta, B, tb, zc, C)
    complex(8),   intent(in)    :: f
    type(matrix), intent(in)    :: A, B
    logical,      intent(in)    :: ta, tb, zc
    type(matrix), intent(inout) :: C
if (matrix_backend_debug) &
print '(a,x,f4.1,x,l,x,l,x,l)', ' gemm f ta tb zc =', dreal(f), ta, tb, zc
    if (ta .and. tb .and. zc) then
       C%elms = f * matmul(transpose(A%elms), &
                           transpose(B%elms))
    else if (ta .and. tb) then
       C%elms = C%elms + f * matmul(transpose(A%elms), &
                                    transpose(B%elms))
    else if (ta .and. zc) then
       C%elms = f * matmul(transpose(A%elms), B%elms)
    else if (ta) then
       C%elms = C%elms + f * matmul(transpose(A%elms), B%elms)
    else if (tb .and. zc) then
       C%elms = f * matmul(A%elms, transpose(B%elms))
    else if (tb) then
       C%elms = C%elms + f * matmul(A%elms, transpose(B%elms))
    else if (zc) then
       C%elms = f * matmul(A%elms, B%elms)
    else
       C%elms = C%elms + f * matmul(A%elms, B%elms)
    end if
  end subroutine



  !> matrix dot product: sum_ij Aij Bij (t=F)
  !>  or  product trace: sum_ij Aji Bij (t=T)
  !> A and B are assumed initialized and with equal/transpose shapes
  function mat_dot(A, B, t)
    type(matrix), intent(in) :: A, B
    logical,      intent(in) :: t !transpose
    complex(8)               :: mat_dot
    if (   t  ) mat_dot = sum(transpose(A%elms) * B%elms)
    if (.not.t) mat_dot = sum(A%elms * B%elms)
    ! if either closed-shell, then additional factor 2
    if (A%closed_shell .or. B%closed_shell) &
       mat_dot = 2*mat_dot
  end function


  !> matrix trace: sum_i Aii, A assumed initialized and square
  function mat_trace(A)
    type(matrix), intent(in) :: A
    complex(8)               :: mat_trace
    integer :: i
    if (A%nrow /= A%ncol) &
       call quit('error: mat_trace(A) with A non-square')
    mat_trace = sum((/(A%elms(i,i),i=1,A%nrow)/))
    ! if closed-shell, then additional factor 2
    if (A%closed_shell) mat_trace = 2*mat_trace
  end function



  !> print matrix A to optional unit, optionally starting with 'label',
  !> optionally making each column 'width' wide, using optional separators
  !> 'sep' (left brace, column separator, right brace, row separator)
  subroutine mat_print(A, label, unit, width, sep)
    type(matrix),           intent(in) :: A
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: sep
    integer      uni, wid, pre, dec, i, j, siz
    character(8) fmt
    character(4) spr
    ! process optional argument unit, which defaults to stdout
    uni = 6
    if (present(unit)) uni = unit
    ! set pre to the largest number of digits to be printed
    ! before the decimal point (including any minus sign)
    pre = pre_decimals(size(A%elms), A%elms)
    if (A%complex) &
       pre = max(pre, pre_decimals(size(A%elms_i), A%elms_i))
    if (A%open_shell) &
       pre = max(pre, pre_decimals(size(A%elms_b), A%elms_b))
    if (A%complex .and. A%open_shell) &
       pre = max(pre, pre_decimals(size(A%elms_ib), A%elms_ib))
    ! process optional width
    wid = 9
    if (present(width)) wid = max(width, max(pre,2)+2) !max, to avoid *****
    ! set dec to the number of decimals to be printed
    dec = wid - max(pre,2) - 1
    ! argument label is optional. If present, print that
    if (present(label)) write (uni,'(a)') label
    ! process optional argument braces, defaulting to none
    spr = '    '
    if (present(sep)) spr = sep
    ! create the format string to be used for each element
    write (fmt,'(a2,i2,a1,i2,a1)') '(f', wid, '.', dec, ')'
    ! call the printing routine
    if (A%complex .and. A%open_shell) write (uni,'(a)') '(real alpha part)'
    if (.not.A%complex .and. A%open_shell) write (uni,'(a)') '(alpha part)'
    if (A%complex .and. .not.A%open_shell) write (uni,'(a)') '(real part)'
    call subr(A%nrow, A%ncol, A%elms, wid, spr, fmt, uni)
    write (uni,'()') !blank line'
    ! imaginary (alpha) part
    if (A%complex) then
       if (.not.A%open_shell) write (uni,'(a)') '(imaginary part)'
       if (     A%open_shell) write (uni,'(a)') '(imaginary alpha part)'
       call subr(A%nrow, A%ncol, A%elms_i, wid, spr, fmt, uni)
       write (uni,'()')
    end if
    ! beta part
    if (A%open_shell) then
       if (.not.A%complex) write (uni,'(a)') '(beta part)'
       if (     A%complex) write (uni,'(a)') '(real beta part)'
       call subr(A%nrow, A%ncol, A%elms_b, wid, spr, fmt, uni)
       write (uni,'()')
    end if
    ! imaginary beta part
    if (A%complex .and. A%open_shell) then
       write (uni,'(a)') '(imaginary beta part)'
       call subr(A%nrow, A%ncol, A%elms_ib, wid, spr, fmt, uni)
       write (uni,'()')
    end if

  contains

    ! number of pre-decimals in the printing of number
    function pre_decimals(n, r)
      integer, intent(in) :: n
      real(8), intent(in) :: r(n)
      integer             :: pre_decimals
      real(8) maxr, minr
      maxr = maxval(r)
      minr = minval(r)
      pre_decimals = ceiling(log10(max(maxr,-minr)))
      if (-10*minr > maxr) pre_decimals = pre_decimals + 1
    end function

    subroutine subr(nrow, ncol, elms, wid, spr, fmt, unit)
      integer,        intent(in) :: nrow, ncol, wid, unit
      real(8),        intent(in) :: elms(nrow,ncol)
      character(4),   intent(in) :: spr
      character(8),   intent(in) :: fmt
      character(ncol*(wid+1)+3) line
      integer i, j, l
      do i = 1, nrow
         line(1:1) = merge(spr(1:1),' ',i==1)
         line(2:2) = spr(1:1)
         l = 3
         do j = 1, ncol
            write (line(l:l+wid-1), fmt) elms(i,j)
            l = l + wid 
            if (j/=ncol) line(l:l) = spr(2:2)
            if (j/=ncol) l = l+1
         end do
         line(l:l) = spr(3:3)
         line(l+1:l+1) = merge(spr(3:3), spr(4:4), i==nrow)
         l = l + 2
         if (present(sep)) then
            write (unit,'(a)') line
         else
            write (unit,'(a)') line(2:len(line)-2)
         end if
      end do
    end subroutine

  end subroutine


  !> copy all fields of A into B, creating a duplicate
  subroutine mat_duplicate(A, B)
    type(matrix), intent(in)    :: A
    type(matrix), intent(inout) :: B
    B = A
  end subroutine


  !> whether A and B are duplicates
  function mat_isdup(A, B)
    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B
    logical                  :: mat_isdup
    mat_isdup = associated(B%elms, A%elms)
  end function


end module


! Copyright 2009-2011 Andreas J. Thorvaldsen
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
!>    equals                       A = B
!>    plus                         A + B
!>    minus                        A - B
!>    matrix multiply              A * B
!>    scale by integer             n * A
!>    scale by real(8)             r * A
!>    scale by complex(8)          z * A
!>    transpose                    trps(A)
!>    trace                        tr(A)
!>    scalar product               dot(A,B)
!>    product trace                tr(A,B)
!>    short-hand freeing           A=0, A(:)=0, A(:,:)=0
!>    norm ie. sqrt(dot(A,A))      norm(A)
!>
!> To increase performance and lower memory needs, formulas are not
!> evaluated immediately as they stand, but the operators are
!> first collected into 'proxy' matrices of the forms
!>
!>    axpy:      Y (+)= f * X(^T)             (where X may be Y)
!>    gemm:      C (+)= f * A(^T) * B(^T)
!>
!> For printing matrices, mat_print from matrix_genop is provided:
!>
!>    call mat_print(A)         (default print without decor)
!>    call mat_print(A,label='A',width=7,unit=6,sep='{,},')
!>
!> The label is printed on the line preceding the matrix. Width specifies the
!> text width of each column (excluding separator). The number of digits printed
!> is adjusted so that no columns become overfilled (*****). Left brace,
!> column separator, right brace, row separator (sep) can be specified, so the
!> matrix can be copy-pasted as input to python, maple, mathematica, etc.
!>
!> To check the state of matrices: Whether defined: isdef(A),
!> whether zero: iszero(A) (zero matrices do not allocate memory)
module matrix_defop

  use matrix_backend, mat_print_nonzero => mat_print, &
                    mat_trace_nonzero => mat_trace

  implicit none

  public matrix
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
  public mat_print
  public matrix_defop_debug

  !> for switching debugging on or off
  logical :: matrix_defop_debug = .false.

  interface assignment(=)
     module procedure mat_eq_mat
     module procedure mat_eq_prx
     module procedure mat_eq_zero
     module procedure mat_eq_zero_1D
     module procedure mat_eq_zero_2D
     module procedure mat_eq_zero_3D
     module procedure mat_eq_zero_4D
  end interface

  interface operator(+)
     module procedure mat_plus_mat
     module procedure mat_plus_prx
     module procedure prx_plus_mat
     module procedure prx_plus_prx
     module procedure plus_mat
     module procedure plus_prx
  end interface

  interface operator(-)
     module procedure mat_minus_mat
     module procedure mat_minus_prx
     module procedure prx_minus_mat
     module procedure prx_minus_prx
     module procedure minus_mat
     module procedure minus_prx
  end interface

  interface operator(*)
     module procedure complex_times_mat
     module procedure complex_times_prx
     module procedure real_times_mat
     module procedure real_times_prx
     module procedure integer_times_mat
     module procedure integer_times_prx
     module procedure mat_times_mat
     module procedure mat_times_prx
     module procedure prx_times_mat
     module procedure prx_times_prx
  end interface

  interface trps
     module procedure transpose_mat
     module procedure transpose_prx
  end interface

  interface dot
     module procedure mat_dot_mat
     module procedure mat_dot_prx
     !ajt No idea why gfortran complains about this
     ! module procedure prx_dot_mat
     module procedure prx_dot_prx
  end interface

  interface tr
     module procedure mat_prodtr_mat
     module procedure mat_prodtr_prx
     !ajt No idea why gfortran complains about this
     ! module procedure prx_prodtr_mat
     module procedure prx_prodtr_prx
     module procedure mat_trace
     module procedure prx_trace
  end interface

  interface norm
     module procedure mat_norm
     module procedure prx_norm
  end interface

  interface mat_print
     module procedure print_mat
     module procedure print_prx
  end interface

  !> private type used to contain intermediates during matrix
  !> algebra evaluation
  type proxy
     !> scale factor before A. Have A whenever f!=0
     complex(8) f
     !> C= not C+=, A(^T), have B, B(^T)
     logical zc, ta, hb, tb
     !> matrices in C+=f*A*B
     type(matrix) C, A, B
  end type

  ! dummy matrix at which %self_pointer of all intermediate
  ! matrices proxy%C/A/B points. This is needed because some compilers
  ! (eg fortran) move type(proxy) instances around while evaluating
  type(matrix), target :: intermediate

  private

  contains

   subroutine quit(a)

   character(*) :: a

   end subroutine

  function isdef(A)
    type(matrix), intent(in) :: A
    logical                  :: isdef
    isdef = (A%magic_tag == mat_magic_value .or. &
             A%magic_tag == mat_magic_setup)
  end function

  subroutine assert_def(A, errmsg)
    type(matrix), target, intent(in) :: A
    character(*),         intent(in) :: errmsg
    if (.not.isdef(A)) call quit(errmsg)
  end subroutine

  function iszero(A)
    type(matrix), intent(in) :: A
    logical                  :: iszero
    iszero = (A%magic_tag == mat_magic_setup)
    if (.not.iszero .and. A%magic_tag /= mat_magic_value) &
       call quit('iszero(A) error: called with matrix A undefined')
  end function

  ! private inline. Same as .not.iszero, but without insisting 
  ! that A is defined
  function nz(A)
    type(matrix), intent(in) :: A
    logical                  :: nz
    nz = .not.(A%magic_tag == mat_magic_setup)
  end function

  ! private inline. Tells whether A is an intermediate
  function intm(A)
    type(matrix), intent(in) :: A
    logical                  :: intm
    intm = associated(A%self_pointer, intermediate)
  end function

  ! private inline. Tells whether A is an intermediate
  function alias(A, B)
    type(matrix), target, intent(in) :: A, B
    logical                          :: alias
    alias = associated(A%self_pointer, B)
  end function

  ! free A if its an intermediate, nullify it if it's an alias
  subroutine free_or_nullify(A)
    type(matrix), target, intent(inout) :: A
    if (intm(A)) then
       A%self_pointer => A
       call mat_free(A)
    else
       call mat_nullify(A)
    end if
  end subroutine


  !> evaluate proxy P: C (+)= f * A(^T) {* B(^T)}
  subroutine eval_proxy(P)
    type(proxy), target :: P
    if (.not.intm(P%C) .and. intm(P%A) .and. &
        .not.P%ta .and. .not.P%hb) then
       if (P%f/=1) &
          call mat_axpy(P%f, P%A, .false., .true., P%A)
       if (nz(P%C)) &
          call mat_axpy((1d0,0d0), P%C, .false., .false., P%A)
       call mat_duplicate(P%A, P%C)
       call mat_nullify(P%A)
       P%f = 0
    else if (.not.nz(P%C)) then
       call mat_alloc(P%C)
       P%C%self_pointer => intermediate
       P%zc = .true.
    else if (.not.intm(P%C)) then
       call realloc_matrix(.false., P%C)
    end if
    if (P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, P%zc, P%C)
    else if (P%f/=0) then
       call mat_axpy(P%f, P%A, P%ta, P%zc, P%C)
    end if
    ! free A and B
    if (P%f/=0) call free_or_nullify(P%A)
    if (P%hb  ) call free_or_nullify(P%B)
    call clear_proxy(P)
  end subroutine


  ! clear fields in proxy. Note that C,A,B are not touched
  subroutine clear_proxy(P)
    type(proxy), intent(out) :: P
    P%f  = 0
    P%zc = .false.
    P%ta = .false.
    P%hb = .false.
    P%tb = .false.
  end subroutine


  !> reallocate a matrix, and optionally change its sign
  subroutine realloc_matrix(minus, A)
    logical,              intent(in)    :: minus
    type(matrix), target, intent(inout) :: A
    type(matrix) tmp
    call mat_duplicate(A, tmp)
    call mat_nullify(A)
    call mat_setup(A, tmp)
    call mat_alloc(A)
    call mat_axpy(merge(-1,1,minus)*(1d0,0d0), tmp, .false., .true., A)
    call mat_nullify(tmp)
    A%self_pointer => intermediate
  end subroutine


  !> private to set C to C, -C, A+C, A-C, C-A or -A-C utilizing
  !> (if need be) A's current allocation
  subroutine zerplusmin_mat_plusmin_mat(za, ma, A, mc, C)
    logical,      intent(in) :: za, ma, mc
    type(matrix), target     :: C, A
    logical CisA
    integer fc, fa
    CisA = alias(C, A)
    ! integer factors of C and A
    fa = merge(0, merge(-1,1,ma), CisA)
    fc = merge(-1,1,mc) + merge(merge(-1,1,ma), 0, CisA)
    ! first branch: A=-A A+=C freeC C<=A
    if (intm(A) .and. (.not.intm(C) .or. (fa==1 .and. fc==-1))) then
       if (fa==-1) call mat_axpy((-1d0,0d0), A, .false., .true., A)
       if (fc/= 0) call mat_axpy(fc*(1d0,0d0), C, .false., fa==0, A)
       if (intm(C)) call mat_free(C)
       call mat_duplicate(A, C)
       call mat_nullify(A)
    else !other branch: reaC C*=-1 C+=A freeA
       if (intm(A) .and. CisA) then
          call mat_nullify(A)
          C%self_pointer => intermediate
       else if (.not.intm(C)) then
          call realloc_matrix(fc==-1, C)
          if (fc==-1) fc = 1
       end if
       if (fc/=1) call mat_axpy(fc*(1d0,0d0), C, .false., .true.,  C)
       if (fa/=0) call mat_axpy(fa*(1d0,0d0), A, .false., .false., C)
       if (intm(A) .and. .not.CisA) call mat_free(A)
    end if
  end subroutine


  !> copy B to A. If A is already allocated and co-shape B,
  !> copy B's elements into A's elements. Otherwise A becomes
  !> an alias of B
  subroutine mat_eq_mat(A, B)
    type(matrix), target, intent(inout) :: A
    type(matrix)                        :: B
    logical haveA
    call assert_def(B, 'matrix A=B error: B undefined')
    haveA = (isdef(A) .and. alias(A, A))
    if (haveA) then
       haveA = (nz(B) .and. mat_coshape(A,B))
       if (.not.haveA) call mat_free(A)
    end if
    if (haveA) then
       call mat_axpy((1d0,0d0), B, .false., .true., A)
    else
       call mat_duplicate(B,A)
    end if
  end subroutine


  !> evaluate matrix proxy P and place the result in A
  subroutine mat_eq_prx(A, P)
    type(matrix), target, intent(inout) :: A
    type(proxy),  target, intent(in)    :: P
    call mat_eq_prx_subr(A, P)
  end subroutine


  !> evaluate matrix proxy P and place the result in A
  subroutine mat_eq_prx_subr(M, P)
    type(matrix), target, intent(inout) :: M
    type(proxy),  target                :: P
    logical zeroP, haveM, ABisM, needM, useM
    zeroP = (.not.nz(P%C) .and. P%f==0)
    haveM = (isdef(M) .and. alias(M, M))
    ! if M is among the operands in P=C+f*A*B, make the last occurance
    if (haveM) then !intermediate, and forget about having M
       if (P%hb .and. alias(P%B, M)) then
          P%B%self_pointer => intermediate
          if (alias(P%A, M)) P%A%self_pointer => P%B
          if (alias(P%C, M)) P%C%self_pointer => P%B
          haveM = .false.
       else if (P%f/=0 .and. alias(P%A, M)) then
          P%A%self_pointer => intermediate
          if (alias(P%C, M)) P%C%self_pointer => P%A
          haveM = .false.
       else if (alias(P%C, M)) then
          P%C%self_pointer => intermediate
          haveM = .false.
       end if
    end if
    ! whether M could find use in the evaluation of P. No use for M
    ! if C is intermediate, or A is alone, non-transposed and intermediate
    needM = .not.(zeroP .or. intm(P%C) .or. (P%f/=0 .and. &
                          .not.P%ta .and. .not.P%hb .and. intm(P%A)))
    ! lastly, in order to use M, it must be correctly shaped
    useM = .false.
    if (haveM .and. needM) useM = mat_coshape(M, P%C)
    ! either swap M in as C, or free it
    if (useM .and. .not.nz(P%C)) then
       call mat_duplicate(M, P%C)
       P%C%self_pointer => intermediate
       P%zc = .true. !C is scratch
    else if (useM) then
       call zerplusmin_mat_plusmin_mat(.true., .false., M, .false., P%C)
    else if (haveM) then
       call mat_free(M)
    end if
    ! if nonzero, evaluate
    if (.not.zeroP) call eval_proxy(P)
    ! transfer result back from C to M
    call mat_duplicate(P%C, M)
    call mat_nullify(P%C)
    if (.not.zeroP) M%self_pointer => M !make M owner
  end subroutine



  !------------------------------------------
  !   add and subtract operators A+B, A-B
  !------------------------------------------

  !> A+B: verify that shapes match, check for zeros
  function mat_plus_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    call assert_def(A, 'matrix A+B error: A undefined')
    call assert_def(B, 'matrix A+B error: B undefined')
    if (.not.mat_coshape(A, B, add=.true.)) &
       call quit('matrix A+B error: different shapes')
    ! initialize P
    call clear_proxy(P)
    call mat_nullify(P%C)
    if (.not.nz(A) .and. .not.nz(B)) then
       call mat_setup(P%C, A) !no allocate for zero
    else if (.not.nz(A)) then
       call mat_duplicate(B, P%C)
    else
       call mat_duplicate(A, P%C)
       if (nz(B)) P%f = 1
       if (nz(B)) call mat_duplicate(B, P%A)
    end if
  end function

  function mat_plus_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call assert_def(A, 'matrix A+(P) error: A undefined')
    if (.not.mat_coshape(A, P%C, add=.true.)) &
       call quit('matrix A+(P) error: different shapes')
    call plusmin_mat_plusmin_prx(.false., A, .false., P, R)
  end function

  function prx_plus_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call assert_def(A, 'matrix (P)+A error: A undefined')
    if (.not.mat_coshape(P%C, A, add=.true.)) &
       call quit('matrix (P)+A error: different shapes')
    call plusmin_mat_plusmin_prx(.false., A, .false., P, R)
  end function

  function prx_plus_prx(P, Q) result(R)
    type(proxy), target, intent(in) :: P, Q
    type(proxy), target             :: R
    if (.not.mat_coshape(P%C, Q%C, add=.true.)) &
       call quit('matrix (P)+(Q) error: different shapes')
    call prx_plusminus_prx(P, .false., Q, R)
  end function

  function mat_minus_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    call assert_def(A, 'matrix A-B error: A undefined')
    call assert_def(B, 'matrix A-B error: B undefined')
    if (.not.mat_coshape(A, B, add=.true.)) &
       call quit('matrix A-B error: different shapes')
    ! initialize P
    call clear_proxy(P)
    call mat_duplicate(A, P%C)
    if (nz(B)) P%f = -1
    if (nz(B)) call mat_duplicate(B, P%A)
  end function  

  function mat_minus_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call assert_def(A, 'matrix A-(P) error: A undefined')
    if (.not.mat_coshape(A, P%C, add=.true.)) &
       call quit('matrix A-(P) error: different shapes')
    call plusmin_mat_plusmin_prx(.false., A, .true., P, R)
  end function

  function prx_minus_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call assert_def(A, 'matrix (P)-A error: A undefined')
    if (.not.mat_coshape(A, P%C, add=.true.)) &
       call quit('matrix (P)-A error: different shapes')
    call plusmin_mat_plusmin_prx(.true., A, .false., P, R)
  end function

  function prx_minus_prx(P, Q) result(R)
    type(proxy), target, intent(in) :: P, Q
    type(proxy), target             :: R
    logical pq, qp, havePA, haveQA
    if (.not.mat_coshape(P%C, Q%C, add=.true.)) &
       call quit('matrix (P)-(Q) error: different shapes')
    call prx_plusminus_prx(P, .true., Q , R)
  end function


  !> R=M+P, R=M-P, R=-M+P or R=-M-P. M may well be zero
  subroutine plusmin_mat_plusmin_prx(mm, M, mp, P, R)
    logical,              intent(in)  :: mm, mp
    type(matrix), target              :: M
    type(proxy),  target              :: P
    type(proxy),  target, intent(out) :: R
    logical useM, useC, haveA, useA, neg
    ! whether A is intermediate, and can keep the result
    haveA = (P%f/=0 .and. intm(P%A) .and. .not.P%ta .and. .not.P%hb)
    ! decide which matrix to use for the result: P%C, A or P%A
    useC = (intm(P%C) .and. .not.mp)
    useM = (.not.useC .and. intm(M) .and. .not.mm)
    useA = (.not.useC .and..not.useM .and. haveA .and. P%f == merge(-1,1,mp))
    if (.not.(useM.or.useC.or.useA)) useC = intm(P%C)
    if (.not.(useM.or.useC.or.useA)) useM = intm(M)
    if (.not.(useM.or.useC.or.useA)) useA = haveA
    ! decide whether we calculate minus the result
    neg = ((useC .and. mp) .or. (useM .and. mm) &
           .or. (useA .and. P%f == merge(1,-1,mp)))
    ! if result goes into A, apply any scale factor other than +-1
    if (useA .and. P%f/=1 .and. P%f/=-1) then
       call mat_axpy(merge(-1,1,mp)*P%f, P%A, .false., .true., P%A)
       P%f = merge(-1,1,mp)
    end if
    ! if result goes in A, add it to C
    if (useA) then
       call zerplusmin_mat_plusmin_mat(.false., .false., P%A, P%f==-1, P%C)
       P%f  = 0
       call mat_nullify(P%A)
    end if
    ! add/subtract A to/from C
    if (nz(M) .and. nz(P%C)) then
       call zerplusmin_mat_plusmin_mat(.false., mm.neqv.neg, M, mp.neqv.neg, P%C)
    else if (nz(M)) then
       call mat_duplicate(M, P%C)
       if (intm(M)) M%self_pointer => P%C
       if (mm.neqv.neg) call realloc_matrix(.true., P%C)
    end if
    ! either transfer P to R, or set R to -P
    if (neg) then
       if (P%f/=0) call eval_proxy(P)
       call clear_proxy(R)
       call mat_nullify(R%C)
       call mat_setup(R%C, P%C) !zero
       R%f = -1
       call mat_duplicate(P%C, R%A)
    else
       R = P
       if (mp) R%f = -R%f
    end if
    ! clear P
    call mat_nullify(P%C)
    if (P%f/=0) call mat_nullify(P%A)
    if (P%hb  ) call mat_nullify(P%B)
    call clear_proxy(P)
  end subroutine
  

  subroutine prx_plusminus_prx(P, mq, Q, R)
    type(proxy), target              :: P, Q
    logical,             intent(in)  :: mq
    type(proxy), target, intent(out) :: R
    logical pq, qp, havePA, haveQA
    ! whether P%A and Q%A are intermatiates and available
    havePA = (P%f/=0 .and. intm(P%A) .and. .not.P%ta .and. .not.P%hb)
    haveQA = (Q%f/=0 .and. intm(Q%A) .and. .not.Q%ta .and. .not.Q%hb)
    ! determine in which order to combine P and Q to R:
    ! Firstly, pick one which misses A, and thus doesn't need evaluation
    pq = (P%f==0)
    qp = (.not.pq .and. Q%f==0)
    ! then pick P if P%C is intermediate or P%A is with factor 1
    if (.not.(pq.or.qp)) pq = (intm(P%C) .or. (havePA .and. P%f==1))
    ! then pick Q if Q%A is intermediate and with net factor 1 (--)
    if (.not.(pq.or.qp)) qp = (haveQA .and. Q%f == merge(-1,1,mq))
    ! then pick P if P%A is intermediate with factor -1
    if (.not.(pq.or.qp)) pq = (havePA .and. P%f==-1)
    ! then pick Q if Q%C is intermediate, or if Q%A with factor 1
    if (.not.(pq.or.qp)) qp = (intm(Q%C) .or. (haveQA .and. Q%f == merge(1,-1,mq)))
    ! then pick P if P%A is intermediate, whatever the factor
    if (.not.(pq.or.qp)) pq = havePA
    ! then pick Q if Q%A is intermediate, whatever the factor
    if (.not.(pq.or.qp)) qp = haveQA
    ! combine in the decided upon order, defaulting to evaluating P first
    if (pq .or. .not.qp) then
       if (P%f/=0) call eval_proxy(P)
       call plusmin_mat_plusmin_prx(.false., P%C, mq, Q, R)
       call free_or_nullify(P%C)
    else
       if (Q%f/=0) call eval_proxy(Q)
       call plusmin_mat_plusmin_prx(mq, Q%C, .false., P, R)
       call free_or_nullify(Q%C)
    end if
  end subroutine


  !----------------------------------
  !   scale operators f*A, +A, -A
  !----------------------------------

  function complex_times_mat(r, A) result(P)
     complex(8),           intent(in) :: r
     type(matrix), target, intent(in) :: A
     type(proxy),  target             :: P
     call number_times_mat(r, A, P)
  end function

  function real_times_mat(r, A) result(P)
    real(8),              intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call number_times_mat(r*(1d0,0d0), A, P)
  end function

  function integer_times_mat(r, A) result(P)
    integer,              intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call number_times_mat(r*(1d0,0d0), A, P)
  end function

  function plus_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call number_times_mat((1d0,0d0), A, P)
  end function

  function minus_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call number_times_mat((-1d0,0d0), A, P)
  end function

  subroutine number_times_mat(f, A, P)
    complex(8),           intent(in) :: f
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call assert_def(A, 'matrix r*A error: A undefined')
    call clear_proxy(P)
    if (f==1 .or. .not.nz(A)) then
       call mat_duplicate(A, P%C)
    else
       call mat_nullify(P%C)
       call mat_setup(P%C, A) !now C=zero
       if (f/=0) P%f = f
       if (f/=0) call mat_duplicate(A, P%A)
    end if
  end subroutine

  function complex_times_prx(f, P) result(Q)
    complex(8),          intent(in) :: f
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: Q
    call number_times_prx(f, P, Q)
  end function

  function real_times_prx(f, P) result(Q)
    real(8),             intent(in) :: f
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: Q
    call number_times_prx(f*(1d0,0d0), P, Q)
  end function

  function integer_times_prx(f, P) result(Q)
    integer,             intent(in) :: f
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: Q
    call number_times_prx(f*(1d0,0d0), P, Q)
  end function

  function plus_prx(P) result(Q)
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: Q
    call number_times_prx((1d0,0d0), P, Q)
  end function

  function minus_prx(P) result(Q)
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: Q
    call number_times_prx((-1d0,0d0), P, Q)
  end function


  subroutine number_times_prx(f, P, R)
    complex(8),          intent(in)  :: f
    type(proxy), target              :: P
    type(proxy), target, intent(out) :: R
    ! if multiplying a non-zero P by zero f, start by cleaning up that
    if (f==0 .and. P%f/=0) then
       P%f = 0
       call free_or_nullify(P%A)
       if (P%hb) call free_or_nullify(P%B)
       P%hb = .false.
    end if
    if (f==0 .and. nz(P%C)) then
       call mat_duplicate(P%C, P%A)
       call free_or_nullify(P%C)
       call mat_setup(P%C, P%A)
       call mat_nullify(P%A)
    end if
    ! if no actual scaling needed, just copy P to R
    if (f==1 .or. .not.nz(P%C)) then
       R = P
       R%f = R%f * f
    else
       ! evaluate unless A is untransposed and without B, and
       ! will end up with factor +1 after multiplication
       if (.not.(f * P%f == 1 .and. .not.P%ta .and. .not.P%hb)) &
          call eval_proxy(P)
       call clear_proxy(R)
       if (P%f/=0) call mat_duplicate(P%A, R%C)
       if (P%f==0) call mat_nullify(R%C)
       if (P%f==0) call mat_setup(R%C, P%C)
       if (nz(P%C)) R%f = f
       if (nz(P%C)) call mat_duplicate(P%C, R%A)
    end if
    ! clear P
    call mat_nullify(P%C)
    if (P%f/=0) call mat_nullify(P%A)
    if (P%hb  ) call mat_nullify(P%B)
    call clear_proxy(P)
  end subroutine



  !----------------------------------
  !   transpose operator trps(A)
  !----------------------------------

  !> transpose of matrix
  function transpose_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call assert_def(A, 'matrix trps(A) error: A undefined')
    call clear_proxy(P)
    call mat_nullify(P%C)
    call mat_setup(P%C, A, ta=.true.)
    ! if nonzero, set P%A to A^T
    if (nz(A)) then
       P%f  =  1
       P%ta =  .true.
       call mat_duplicate(A, P%A)
    end if
  end function


  !> transpose of matrix proxy
  function transpose_prx(P) result(R)
    type(proxy), target :: P
    type(proxy), target :: R
    logical swap
    call clear_proxy(R)
    ! if P=f*A or P=f*A*B, make R=f*A^T or R=f*B^T*A^T
    if (.not.nz(P%C) .and. .not.(P%f==1 .and. P%ta .and. .not.P%hb)) then
       call mat_duplicate(P%C, R%C) !zero
       R%f = P%f
       if (P%hb) then
          R%hb = .true.
          R%ta = .not.P%tb
          R%tb = .not.P%ta
          call mat_duplicate(P%B, R%A)
          call mat_duplicate(P%A, R%B)
       else
          R%ta = .not.P%ta
          call mat_duplicate(P%A, R%A)
       end if
    else
       ! must evaluate if not P=C+1*A^T
       if (.not.(P%f==1 .and. P%ta .and. .not.P%hb)) &
          call eval_proxy(P)
       if (P%f==0) call mat_nullify(R%C)
       if (P%f==0) call mat_setup(R%C, P%A, ta=.not.P%ta)
       if (P%f/=0) call mat_duplicate(P%A, R%C)
       R%f = merge(0, 1, .not.nz(P%C))
       if (R%f/=0) call mat_duplicate(P%C, R%A)
    end if
    ! clear P
    call mat_nullify(P%C)
    if (P%f/=0) call mat_nullify(P%A)
    if (P%hb  ) call mat_nullify(P%B)
    call clear_proxy(P)
  end function



  !-----------------------------------------
  !   matrix-matrix multiply operator A*B
  !-----------------------------------------

  !> verify that shapes match, check for zeroes
  function mat_times_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    call assert_def(A, 'matrix A*B error: A undefined')
    call assert_def(B, 'matrix A*B error: B undefined')
    if (.not.mat_coshape(A, B, mul=.true.)) &
       call quit('matrix A*B error: incompatible shapes')
    call clear_proxy(P)
    call mat_nullify(P%C)
    call mat_setup(P%C, A, B) !zero
    if (nz(A) .and. nz(B)) then
       P%f  = 1
       P%hb = .true.
       call mat_duplicate(A, P%A)
       call mat_duplicate(B, P%B)
    end if
  end function

  !> verify that shapes match, call routine below
  function mat_times_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call assert_def(A, 'matrix A*(P) error: A undefined')
    if (.not.mat_coshape(A, P%C, mul=.true.)) &
       call quit('matrix A*(P) error: incompatible shapes')
    call mat_revtimes_prx(A, .false., .false., P, R)
  end function

  !> verify that shapes match, call routine below
  function prx_times_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call assert_def(A, 'matrix (P)*A error: A undefined')
    if (.not.mat_coshape(P%C, A, mul=.true.)) &
       call quit('matrix (P)*A error: incompatible shapes')
    call mat_revtimes_prx(A, .false., .true., P, R)
  end function

  !> verify that shapes match, call routine below
  function prx_times_prx(P, Q) result(R)
    type(proxy),  target, intent(in) :: P, Q
    type(proxy),  target             :: R
    logical pnocb, qnocb
    if (.not.mat_coshape(P%C, Q%C, mul=.true.)) &
       call quit('matrix (P)*(Q) error: incompatible shapes')
    pnocb = (.not.nz(P%C) .and. .not.P%hb)
    qnocb = (.not.nz(Q%C) .and. .not.Q%hb)
    if ((pnocb .or. P%f==0) .or. .not.(qnocb .or. Q%f==0)) then
       call prx_revtimes_prx(P, .false., Q, R)
    else
       call prx_revtimes_prx(Q, .true.,  P, R)
    end if
  end function


  !> R = A^ta * P or R = P * A^ta (rev=true)
  subroutine mat_revtimes_prx(A, ta, rev, P, R)
    logical,              intent(in)  :: ta, rev
    type(matrix), target              :: A
    type(proxy),  target              :: P
    type(proxy),  target, intent(out) :: R
    logical zer
    ! initialize R and R%C
    call clear_proxy(R)
    call mat_nullify(R%C)
    if (.not.rev) call mat_setup(R%C, A, P%C, ta=ta) !zero
    if (     rev) call mat_setup(R%C, P%C, A, tb=ta) !zero
    ! if multiplying a non-zero P by zero A, start by cleaning up that
    zer = (.not.nz(A) .or. (.not.nz(P%C) .and. P%f==0))
    if (zer .and. P%f/=0) then
       P%f = 0
       call free_or_nullify(P%A)
       if (P%hb) call free_or_nullify(P%B)
       P%hb = .false.
    end if
    if (zer .and. nz(P%C)) then
       call mat_duplicate(P%C, P%A)
       call free_or_nullify(P%C)
       call mat_setup(P%C, P%A)
       call mat_nullify(P%A)
    end if
    ! if multiplying non-zero A by zero P, clean that up too
    if (zer .and. intm(A)) then
       call mat_duplicate(A, P%A)
       call free_or_nullify(A)
       call mat_setup(A, P%A)
       call mat_nullify(P%A)
    end if
    ! evaluate P if it has both C and A, or both A and B
    if ((nz(P%C) .and. P%f/=0) .or. P%hb) &
       call eval_proxy(P)
    ! make it R = A*P%C, f*A*P%A or reverse
    if (.not.zer) R%f = merge((1d0,0d0), P%f, P%f==0)
    if (.not.zer) R%hb = .true.
    if (.not.zer .and. rev) then
       R%tb = ta
       call mat_duplicate(A, R%B)
       if (nz(P%C)) call mat_duplicate(P%C, R%A)
       if (P%f/=0) R%ta = P%ta
       if (P%f/=0) call mat_duplicate(P%A, R%A)
    else if (.not.zer) then
       R%ta = ta
       call mat_duplicate(A, R%A)
       if (nz(P%C)) call mat_duplicate(P%C, R%B)
       if (P%f/=0) R%tb = P%ta
       if (P%f/=0) call mat_duplicate(P%A, R%B)
    end if
    ! finally clear P
    call mat_nullify(P%C)
    if (P%f/=0) call mat_nullify(P%A)
    if (P%hb  ) call mat_nullify(P%B)
    call clear_proxy(P)
  end subroutine


  !> evaluate P if needed, then do R = P%C * Q
  subroutine prx_revtimes_prx(P, rev, Q, R)
    logical,              intent(in)  :: rev
    type(proxy),  target              :: P, Q
    type(proxy),  target, intent(out) :: R
    ! evaluate P if it has both C and A, or both A and B
    if ((nz(P%C) .and. P%f/=0) .or. P%hb) &
       call eval_proxy(P)
    ! perform multiplication, apply original factor on P%A
    if (P%f==0) call mat_revtimes_prx(P%C, .false., rev, Q, R)
    if (P%f/=0) call mat_revtimes_prx(P%A,  P%ta,   rev, Q, R)
    if (P%f/=0) R%f = R%f * P%f
    ! finally clear P
    call mat_nullify(P%C)
    if (P%f/=0) call mat_nullify(P%A)
    call clear_proxy(P)
  end subroutine



  !--------------------------------------------
  !   dot(A,B), tr(A,B), tr(A) and norm(A)
  !--------------------------------------------

  function mat_dot_mat(A, B)
    type(matrix), target, intent(in) :: A, B
    complex(8)                       :: mat_dot_mat
    call assert_def(A, 'matrix dot(A,B) error: A undefined')
    call assert_def(B, 'matrix dot(A,B) error: B undefined')
    if (.not.mat_coshape(A, B, add=.true.)) &
       call quit('matrix dot(A,B) error: different shapes')
    mat_dot_mat = 0
    if (nz(A) .and. nz(B)) &
       mat_dot_mat = mat_dot(A, B, .false.)
  end function

  function mat_prodtr_mat(A, B)
    type(matrix), target :: A, B
    complex(8)           :: mat_prodtr_mat
    call assert_def(A, 'matrix tr(A*B) error: A undefined')
    call assert_def(B, 'matrix tr(A*B) error: B undefined')
    if (.not.mat_coshape(A, B, ta=.true., add=.true.)) &
       call quit('matrix tr(A*B) error: incompatible shapes')
    mat_prodtr_mat = 0
    if (nz(A) .and. nz(B)) &
       mat_prodtr_mat = mat_dot(A, B, .true.)
  end function

  function mat_dot_prx(A, P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    complex(8)                       :: mat_dot_prx
    call assert_def(A, 'matrix dot(A,(P)) error: A undefined')
    if (.not.mat_coshape(A, P%C, add=.true.)) &
       call quit('matrix dot(A,(P)) error: different shapes')
    call mat_dot_prx_subr(A, P, .false., mat_dot_prx)
  end function

  function mat_prodtr_prx(A, P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    complex(8)                       :: mat_prodtr_prx
    call assert_def(A, 'matrix tr(A*(P)) error: A undefined')
    if (.not.mat_coshape(A, P%C, ta=.true., add=.true.)) &
       call quit('matrix tr(A*(P)) error: different shapes')
    call mat_dot_prx_subr(A, P, .true., mat_prodtr_prx)
  end function

  function prx_dot_mat(P, A)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    complex(8)                       :: prx_dot_mat
    call assert_def(A, 'matrix dot((P),A) error: A undefined')
    if (.not.mat_coshape(A, P%C, add=.true.)) &
       call quit('matrix dot((P),A) error: different shapes')
    call mat_dot_prx_subr(A, P, .false., prx_dot_mat)
  end function

  function prx_prodtr_mat(P, A)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    complex(8)                       :: prx_prodtr_mat
    call assert_def(A, 'matrix tr((P)*A) error: A undefined')
    if (.not.mat_coshape(A, P%C, ta=.true., add=.true.)) &
       call quit('matrix tr((P)*A) error: different shapes')
    call mat_dot_prx_subr(A, P, .true., prx_prodtr_mat)
  end function

  function prx_dot_prx(P, Q)
    type(proxy),  target, intent(in) :: P, Q
    complex(8)                       :: prx_dot_prx
    if (.not.mat_coshape(P%C, Q%C, add=.true.)) &
       call quit('matrix dot((P),(Q)) error: different shapes')
    call prx_dot_prx_subr(Q, P, .false., prx_dot_prx)
  end function

  function prx_prodtr_prx(P, Q)
    type(proxy),  target, intent(in) :: P, Q
    complex(8)                       :: prx_prodtr_prx
    if (.not.mat_coshape(P%C, Q%C, ta=.true., add=.true.)) &
       call quit('matrix tr((P)*(Q)) error: incompatible shapes')
    call prx_dot_prx_subr(P, Q, .true., prx_prodtr_prx)
  end function


  subroutine mat_dot_prx_subr(A, P, t, r)
    type(matrix), target      :: A
    type(proxy),  target      :: P
    logical,      intent(in)  :: t !t=F:dot(A,P), t=T: tr(A*P)
    complex(8),   intent(out) :: r
    ! evaluate if P is of form C+f*A*B
    if (nz(A) .and. P%hb) &
       call eval_proxy(P)
    ! calculate only if non-zero
    r = 0
    if (nz(A) .and. nz(P%C)) &
       r = r + mat_dot(A, P%C, t)
    if (nz(A) .and. P%f/=0) &
       r = r + mat_dot(A, P%A, t.neqv.P%ta) * P%f
    ! clean-up P
    call free_or_nullify(P%C)
    if (P%f/=0) call free_or_nullify(P%A)
    call clear_proxy(P)
  end subroutine


  subroutine prx_dot_prx_subr(P, Q, t, r)
    type(proxy),  target      :: P, Q
    logical,      intent(in)  :: t !t=F:dot(P,Q), t=T: tr(P*Q)
    complex(8),   intent(out) :: r
    logical pq, qp !dot(P,Q) or dot(Q,P)
    ! firstly prefer to have a second factor of the form C+f*A
    pq = (nz(Q%C) .and. .not.Q%hb)
    qp = (.not.pq .and. nz(P%C) .and. .not.P%hb)
    ! secondly, prefer to have a first factor of the form C or f*A
    if (.not.(pq.or.qp)) pq = (P%f==0 .or. (.not.nz(P%C) .and. .not.P%hb))
    if (.not.(pq.or.qp)) qp = (Q%f==0 .or. (.not.nz(Q%C) .and. .not.Q%hb))
    ! calculate dot(P,Q) or dot(Q,P), defaulting to the former
    if (pq .or. .not.qp) then
       call inner(P, Q)
    else
       call inner(Q, P)
    end if
  contains
    subroutine inner(P, Q)
      type(proxy), target :: P, Q
      ! evaluate P if it involves more than one matrix
      if (nz(P%C) .and. P%f/=0 .and. &
          .not.(.not.nz(Q%C) .and. Q%f==0)) call eval_proxy(P)
      ! now, use mat_dot_prx_subr on either of (P%C,Q) or (P%A,Q)
      r = 0
      if (P%f==0) call mat_dot_prx_subr(P%C, Q, t, r)
      if (P%f/=0) call mat_dot_prx_subr(P%A, Q, P%ta.neqv.t, r)
      if (P%f/=0) r = r * P%f
      call free_or_nullify(P%C)
      if (P%f/=0) call free_or_nullify(P%A)
      call clear_proxy(P)
    end subroutine
  end subroutine


  function mat_trace(A)
    type(matrix), target :: A
    complex(8)           :: mat_trace
    call assert_def(A, 'matrix tr(A) error: A undefined')
    if (.not.mat_coshape(A, A, ta=.true., mul=.true.)) &
       call quit('matrix tr(A) error: A not square')
    mat_trace = 0
    if (nz(A)) mat_trace = mat_trace_nonzero(A)
  end function


  function prx_trace(P) result(r)
    type(proxy), target :: P
    complex(8)          :: r
    if (.not.mat_coshape(P%C, P%C, ta=.true., mul=.true.)) &
       call quit('matrix tr(P) error: P not square')
    r = 0
    if (nz(P%C)) &
       r = r + mat_trace_nonzero(P%C)
    if (P%hb) then
       r = r + P%f * mat_dot(P%A, P%B, t = P%ta.neqv.P%tb)
    else if (P%f/=0) then
       r = r + P%f * mat_trace_nonzero(P%A)
    end if
    ! free
    call free_or_nullify(P%C)
    if (P%f/=0) call free_or_nullify(P%A)
    if (P%hb  ) call free_or_nullify(P%B)
    call clear_proxy(P)
  end function


  function mat_norm(A)
    type(matrix), target :: A
    real(8)              :: mat_norm
    call assert_def(A, 'matrix norm(A) error: A undefined')
    mat_norm = 0
    if (nz(A)) mat_norm = sqrt(real(mat_dot(A, A, .true.)))
  end function


  function prx_norm(P) result (r)
    type(proxy), target :: P
    real(8)             :: r
    ! evaluate if more than one matrix in P, ie. P=C+f*A or P=f*A*B
    if ((nz(P%C) .and. P%f/=0) .or. P%hb) &
       call eval_proxy(P)
    ! calculate norm
    r = 0
    if (P%f/=0) then
       r = abs(P%f) * sqrt(real(mat_dot(P%A, P%A, .false.)))
    else if (nz(P%C)) then
       r = sqrt(real(mat_dot(P%C, P%C, .false.)))
    end if
    ! free
    call free_or_nullify(P%C)
    if (P%f/=0) call free_or_nullify(P%A)
    call clear_proxy(P)
  end function




  !---------------------------------------------------------------------------
  !   operators for freeing: A=0, A(:)=0, A(:,:)=0, A(:,:,:)=0, A(:,:,:,:)=0
  !---------------------------------------------------------------------------

  !> Short-hand A=0 frees matrix A if it is self, otherwise nullifies
  subroutine mat_eq_zero(A, zero)
    type(matrix), target, intent(inout) :: A
    integer,              intent(in)    :: zero
    if (zero/=0) call quit('matrix A = z error: z must be 0')
    if (alias(A, A)) then
       call mat_free(A)
    else
       call mat_nullify(A)
    end if
  end subroutine

  subroutine mat_eq_zero_1D(A, z)
    type(matrix), target, intent(inout) :: A(:)
    integer,              intent(in)    :: z
    integer i
    do i=1, size(A)
       A(i) = z
    end do
  end subroutine

  subroutine mat_eq_zero_2D(A, z)
    type(matrix), target, intent(inout) :: A(:,:)
    integer,              intent(in)    :: z
    integer j
    do j=1, size(A,2)
       A(:,j) = z
    end do
  end subroutine

  subroutine mat_eq_zero_3D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:)
    integer,              intent(in)    :: z
    integer k
    do k=1, size(A,3)
       A(:,:,k) = z
    end do
  end subroutine

  subroutine mat_eq_zero_4D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:,:)
    integer,              intent(in)    :: z
    integer l
    do l=1, size(A,4)
       A(:,:,:,l) = z
    end do
  end subroutine



  !------------
  ! printing
  !------------

  !> wrap matrix_genop's mat_print so that also zero matrix is printed correctly
  subroutine print_mat(A, label, unit, width, sep)
    type(matrix), target               :: A
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: sep
    call assert_def(A, 'matrix mat_print(A) error: A undefined')
    if (.not.nz(A)) then
       call print_prx(0*A, label, unit, width, sep)
    else
       call mat_print_nonzero(A, label, unit, width, sep)
    end if
  end subroutine

  !> print proxy
  subroutine print_prx(P, label, unit, width, sep)
    type(proxy),  target               :: P
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: sep
    ! if zero, allocate and fill with zero
    if (.not.nz(P%C) .and. P%f==0) then
       call mat_alloc(P%C)
       P%C%self_pointer => intermediate
       call mat_axpy((0d0,0d0), P%C, .false., .true., P%C)
    else if (P%f/=0) then
       call eval_proxy(P)
    end if
    ! print
    call mat_print_nonzero(P%C, label, unit, width, sep)
    ! free
    call free_or_nullify(P%C)
    if (P%f/=0) call free_or_nullify(P%A)
    call clear_proxy(P)
  end subroutine


end module




module proprty

  ! Use blocks

  use matrix_defop

  implicit none

  ! Define (dummy) rsp_cfg type

  type rsp_cfg

     integer :: i

  end type

!   ! Define (dummy) matrix type
! 
!   type matrix
! 
!      integer :: i
! 
!   end type

  ! Define pert datatype

  type perts

     integer :: npert ! Number of perturbations
     integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
     character(4), allocatable, dimension(:) :: plab ! Perturbation labels
     integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
     complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
     ! Add other perturbation identification info as needed

  end type

  ! Define S, D, F datatype

  type SDF

     type(SDF), pointer :: next
     logical :: last
     ! Should all of the data attributes be pointers too?
     type(perts) :: perturb
     integer :: data
     ! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!    type(matrix), allocatable, dimension(:) :: data ! Tensor data

  end type

  type propcache

     type(propcache), pointer :: next
     logical :: last
     integer :: tlen
     type(perts), allocatable, dimension(:) :: t_orders
     complex(8), allocatable, dimension(:) :: data
     !integer :: s_sep
     ! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!    type(matrix), allocatable, dimension(:) :: data ! Property data    

  end type 



  contains

! function pert_clone(pert)
! 
! implicit none
! 
! type(perts) :: pert, pert_clone
! 
! 
! pert_clone%npert = pert%npert
! allocate(pert_clone%pdim(pert%npert))
! allocate(pert_clone%plab(pert%npert))
! allocate(pert_clone%pid(pert%npert))
! allocate(pert_clone%freq(pert%npert))
! 
! 
! pert_clone%pdim = pert%pdim
! pert_clone%plab= pert%plab
! pert_clone%pid = pert%pid
! pert_clone%freq = pert%freq
! 
! end function


  function pert_ext(pert, ext)

    implicit none

    type(perts) :: pert, ext, pert_ext
    integer :: i

!   write(*,*) 'extending pert', pert%npert

    allocate(pert_ext%pdim(pert%npert + 1))
    allocate(pert_ext%plab(pert%npert + 1))
    allocate(pert_ext%pid(pert%npert + 1))
    allocate(pert_ext%freq(pert%npert + 1))
!   write(*,*) 'allocated extended pert', pert%npert


    if (pert%npert == 0) then

       pert_ext%npert = pert%npert + 1
       pert_ext%pdim = (/ext%pdim(:)/)
       pert_ext%plab = (/ext%plab(:)/)
       pert_ext%pid = (/ext%pid(:)/)
       pert_ext%freq = (/ext%freq(:)/)

       ! Other perturbation information as needed

    else

       pert_ext%npert = pert%npert + 1
       pert_ext%pdim = (/(pert%pdim(i), i = 1, pert%npert), ext%pdim(:)/)
       pert_ext%plab = (/(pert%plab(i), i = 1, pert%npert), ext%plab(:)/)
       pert_ext%pid = (/(pert%pid(i), i = 1, pert%npert), ext%pid(:)/)
       pert_ext%freq = (/(pert%freq(i), i = 1, pert%npert), ext%freq(:)/)
 
       ! Other perturbation information as needed



end if

  end function


  function pert_getone(pert, which)

    implicit none

    type(perts) :: pert, pert_getone
    integer :: which

    allocate(pert_getone%pdim(1))
    allocate(pert_getone%plab(1))
    allocate(pert_getone%pid(1))
    allocate(pert_getone%freq(1))

!     write(*,*) 'Getting a pert', which
!     write(*,*) pert%plab


    pert_getone%npert = 1
    pert_getone%pdim = (/pert%pdim(which)/)
    pert_getone%plab = (/pert%plab(which)/)
    pert_getone%pid = (/pert%pid(which)/)
    pert_getone%freq = (/pert%freq(which)/)
    ! Other perturbation information as needed

! write(*,*) 'got one pert'
  end function


  function pert_rf(pert)

    implicit none

    type(perts) :: pert, pert_rf

!        write(*,*) 'rf', pert%npert

    allocate(pert_rf%pdim(pert%npert - 1))
    allocate(pert_rf%plab(pert%npert - 1))
    allocate(pert_rf%pid(pert%npert - 1))
    allocate(pert_rf%freq(pert%npert - 1))

    if (pert%npert > 1) then

       pert_rf%npert = pert%npert - 1
       pert_rf%pdim = (/pert%pdim(2:pert%npert)/)
       pert_rf%plab = (/pert%plab(2:pert%npert)/)
       pert_rf%pid = (/pert%pid(2:pert%npert)/)
       pert_rf%freq = (/pert%freq(2:pert%npert)/)
       ! Other perturbation information as needed

    else

       pert_rf%npert = 0
       pert_rf%pdim = (/0/)
       pert_rf%plab = (/'NUTN'/)
       pert_rf%pid = (/0/)
       pert_rf%freq = (/0.0/)
       ! Other perturbation information as needed

!        write(*,*) 'finished pert_rf'

    end if



  end function

  function merge_pert(p1, p2)

    implicit none

    type(perts) :: p1, p2, merge_pert

!  write(*,*) 'allocating in merge_pert', p1%npert, p2%npert

    allocate(merge_pert%pdim(p1%npert + p2%npert))
    allocate(merge_pert%plab(p1%npert + p2%npert))
    allocate(merge_pert%pid(p1%npert + p2%npert))
    allocate(merge_pert%freq(p1%npert + p2%npert))

!  write(*,*) 'merge_pert allocation successful'

    ! It could be unnecessary to distinguish between all the following cases

    if (p1%npert > 0 .AND. p2%npert > 0) then

! write(*,*) 'case 1'

       merge_pert%npert = p1%npert + p2%npert
       merge_pert%pdim = (/p1%pdim(:), p2%pdim(:)/)
       merge_pert%plab = (/p1%plab(:), p2%plab(:)/)
       merge_pert%pid = (/p1%pid(:), p2%pid(:)/)
       merge_pert%freq = (/p1%freq(:), p2%freq(:)/)
       ! Other perturbation information as needed

    elseif (p1%npert > 0 .AND. p2%npert == 0) then

! write(*,*) 'case 2'

       merge_pert%npert = p1%npert
       merge_pert%pdim = p1%pdim(:)
       merge_pert%plab = p1%plab(:)
       merge_pert%pid = p1%pid(:)
       merge_pert%freq = p1%freq(:)
       ! Other perturbation information as needed

    elseif (p1%npert == 0 .AND. p2%npert > 0) then

! write(*,*) 'case 3'



       merge_pert%npert = p2%npert
       merge_pert%pdim = p2%pdim(:)
       merge_pert%plab = p2%plab(:)
       merge_pert%pid = p2%pid(:)
       merge_pert%freq = p2%freq(:)
       ! Other perturbation information as needed

    elseif (p1%npert == 0 .AND. p2%npert == 0) then

! write(*,*) 'case 4'

       ! Maybe update the following lines to return a proper emptypert

       merge_pert%npert = 0
       merge_pert%pdim = (/0/)
       merge_pert%plab = (/'NUTN'/)
       merge_pert%pid = (/0/)
       merge_pert%freq = (/0.0/)
       ! Other perturbation information as needed

    else

       write(*,*) 'Error in merge_pert: Unrecognized size of p1 or p2 or both:', &
                   p1%npert, p2%npert

    end if

  end function


function get_emptypert() result(emptypert)

type(perts) :: emptypert

emptypert%npert = 0
allocate(emptypert%pdim(0))    
allocate(emptypert%plab(0))
allocate(emptypert%pid(0))
allocate(emptypert%freq(0))

end function

! ASSUMES STANDARD ORDERED INPUT

  function pert_lt(p1, p2)

    implicit none

    logical :: pert_lt
    integer :: i
    type(perts) :: p1, p2

!      write(*,*) 'comparing perturbations'
!  write(*,*) 'p2'
! write(*,*) p2%npert
! write(*,*) p2%pdim
! write(*,*) p2%plab
! write(*,*) p2%freq
!  write(*,*) 'p1'
! write(*,*) p1%npert
! write(*,*) p1%pdim
! write(*,*) p1%plab
! write(*,*) p1%freq

! Will this give false negatives?

pert_lt = .FALSE.

    ! Compare number of perturbations
    ! REMEMBER: DECREASING ORDER OF DIFFERENTIATION
    if (p1%npert > p2%npert) then

        pert_lt = .TRUE.


    elseif (p1%npert == p2%npert) then


       do i = 1, p1%npert


if (llt(p1%plab(i), p2%plab(i)) .eqv. .TRUE.) then

pert_lt = .TRUE.  
exit

elseif (p1%plab(i) == p2%plab(i)) then

! IS IT OK TO COMPARE ONLY THE REAL PART LIKE THIS?
if (real(p1%freq(i)) < real(p2%freq(i))) then

pert_lt = .TRUE.  
exit

end if

if (pert_lt .eqv. .TRUE.) exit

end if

end do

end if

  end function



  function pert_standardorder(pert) result(pert_st)

implicit none


type(perts) :: pert, pert_st
integer :: i, j, min_which, nmin
integer :: tmp_pdim, tmp_pid, t_first, t_last
character(4) :: tmp_plab, min_curr
complex(8) :: tmp_freq, min_freq_curr

! write (*,*) 'standard ordering', pert%pid


pert_st%npert = pert%npert
allocate(pert_st%pdim(pert%npert))    
allocate(pert_st%plab(pert%npert))
allocate(pert_st%pid(pert%npert))
allocate(pert_st%freq(pert%npert))

pert_st%pdim = pert%pdim
pert_st%plab = pert%plab
pert_st%pid = pert%pid
pert_st%freq = pert%freq


nmin = 1

! if (pert_st%npert > 1) then

do i = nmin, pert_st%npert

min_curr = pert_st%plab(i)
min_freq_curr = pert_st%freq(i)
min_which = i

do j = i + 1, pert_st%npert

if (lle(pert_st%plab(j), min_curr) .EQV. .TRUE.) then

if (pert_st%plab(j) == min_curr) then

! MR: COMPARING ABSOLUTE VALUE (SQUARE MODULUS) - IS IT STILL SUFFICIENTLY GENERAL?
if (abs(pert_st%freq(j)) < abs(min_freq_curr)) then

! write (*,*) 'freq', pert_st%freq(j), 'is less than', min_freq_curr

min_freq_curr = pert_st%freq(j)
min_which = j

else

! write (*,*) 'freq', pert_st%freq(j), 'not less than', min_freq_curr


end if

else

min_curr = pert_st%plab(j)
min_which = j

end if

end if

end do

tmp_pdim = pert_st%pdim(min_which)
tmp_plab = pert_st%plab(min_which)
tmp_pid = pert_st%pid(min_which)
tmp_freq = pert_st%freq(min_which)

pert_st%pdim(min_which) = pert_st%pdim(i)
pert_st%plab(min_which) = pert_st%plab(i)
pert_st%pid(min_which) = pert_st%pid(i)
pert_st%freq(min_which) = pert_st%freq(i)

pert_st%pdim(i) = tmp_pdim
pert_st%plab(i) = tmp_plab
pert_st%pid(i) = tmp_pid
pert_st%freq(i) = tmp_freq


nmin = nmin + 1

end do


! DOES THIS WORK?

! write (*,*) 'standard ordered before freq comparison', pert_st%pid


t_first = 1
t_last = 1


do while (t_last <= pert_st%npert)



if (t_last < pert_st%npert) then

do while ((pert_st%plab(t_last) == pert_st%plab(t_first)))

t_last = t_last + 1

if (t_last > pert_st%npert) exit

end do

t_last = t_last - 1

end if

! write (*,*) 'first and last', t_first, t_last

do i = t_first, t_last, 1

min_freq_curr = pert_st%freq(i)
min_which = i

do j = i + 1, t_last

if (abs(pert_st%freq(j)) < abs(min_freq_curr)) then

! write (*,*) 'freq', pert_st%freq(j), 'is less than', min_freq_curr

min_freq_curr = pert_st%freq(j)
min_which = j

end if

end do

tmp_pdim = pert_st%pdim(min_which)
tmp_plab = pert_st%plab(min_which)
tmp_pid = pert_st%pid(min_which)
tmp_freq = pert_st%freq(min_which)

pert_st%pdim(min_which) = pert_st%pdim(i)
pert_st%plab(min_which) = pert_st%plab(i)
pert_st%pid(min_which) = pert_st%pid(i)
pert_st%freq(min_which) = pert_st%freq(i)

pert_st%pdim(i) = tmp_pdim
pert_st%plab(i) = tmp_plab
pert_st%pid(i) = tmp_pid
pert_st%freq(i) = tmp_freq



end do

t_last = t_last + 1
t_first = t_last

end do



! write (*,*) 'standard ordered after freq comparison', pert_st%pid
! else

! ! Do nothing

! end if

end function


! Begin propcache linked list manipulation/data retrieval routines



  ! Initialization routine
  subroutine propcache_init(new_elem, tlen, perturbs, propsize, data)

    implicit none

    integer :: i, tlen
    type(propcache) :: new_elem
    type(perts), dimension(tlen) :: perturbs
    integer :: propsize
    complex(8), dimension(propsize) :: data
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!    type(matrix), dimension(product(pert%pdim)) :: data








    new_elem%last = .TRUE.

    new_elem%tlen = tlen

allocate(new_elem%t_orders(tlen))

do i = 1, tlen

    new_elem%t_orders(i)%npert = perturbs(i)%npert

    allocate(new_elem%t_orders(i)%pdim(perturbs(i)%npert))
    allocate(new_elem%t_orders(i)%plab(perturbs(i)%npert))
    allocate(new_elem%t_orders(i)%pid(perturbs(i)%npert))
    allocate(new_elem%t_orders(i)%freq(perturbs(i)%npert))
    
    new_elem%t_orders(i)%pdim = perturbs(i)%pdim
new_elem%t_orders(i)%plab = perturbs(i)%plab
new_elem%t_orders(i)%pid = perturbs(i)%pid
new_elem%t_orders(i)%freq = perturbs(i)%freq
! write(*,*) new_elem%perturb%npert

end do

allocate(new_elem%data(propsize))

    new_elem%data = data


  end subroutine


  function perts_standardorder(tlen, t_orders) result(perts_st)

implicit none

integer :: tlen
type(perts), dimension(tlen) :: t_orders, perts_st
type(perts) :: tmp_pert
integer :: i, j, k, m_which, nmin, max_order_curr, len_curr
integer :: tmp_pdim, tmp_pid, t_first, t_last
character(4) :: tmp_plab
character(4), dimension(:), allocatable :: min_plab_curr
complex(8) :: tmp_freq
complex(8), dimension(:), allocatable :: min_freq_curr

! write(*,*) 'starting perts_standardorder'

do i = 1, tlen

! write(*,*) t_orders(i)%pid


perts_st(i)%npert = t_orders(i)%npert
allocate(perts_st(i)%pdim(t_orders(i)%npert))    
allocate(perts_st(i)%plab(t_orders(i)%npert))
allocate(perts_st(i)%pid(t_orders(i)%npert))
allocate(perts_st(i)%freq(t_orders(i)%npert))

perts_st(i)%pdim = t_orders(i)%pdim
perts_st(i)%plab = t_orders(i)%plab
perts_st(i)%pid = t_orders(i)%pid
perts_st(i)%freq = t_orders(i)%freq

end do

do i = 2, tlen

m_which = i

do j = i + 1, tlen


if (pert_lt(pert_standardorder(t_orders(j)), pert_standardorder(t_orders(m_which)))) then

m_which = j

end if

end do


! It might be possible to do this assignment in a different way

tmp_pert%npert = perts_st(m_which)%npert
allocate(tmp_pert%pdim(tmp_pert%npert))
allocate(tmp_pert%plab(tmp_pert%npert))
allocate(tmp_pert%pid(tmp_pert%npert))
allocate(tmp_pert%freq(tmp_pert%npert))

tmp_pert%pdim = perts_st(m_which)%pdim
tmp_pert%plab = perts_st(m_which)%plab
tmp_pert%pid = perts_st(m_which)%pid
tmp_pert%freq = perts_st(m_which)%freq

deallocate(perts_st(m_which)%pdim)
deallocate(perts_st(m_which)%plab)
deallocate(perts_st(m_which)%pid)
deallocate(perts_st(m_which)%freq)

perts_st(m_which)%npert = perts_st(i)%npert
allocate(perts_st(m_which)%pdim(perts_st(i)%npert))
allocate(perts_st(m_which)%plab(perts_st(i)%npert))
allocate(perts_st(m_which)%pid(perts_st(i)%npert))
allocate(perts_st(m_which)%freq(perts_st(i)%npert))

perts_st(m_which)%pdim = perts_st(i)%pdim
perts_st(m_which)%plab = perts_st(i)%plab
perts_st(m_which)%pid = perts_st(i)%pid
perts_st(m_which)%freq = perts_st(i)%freq

deallocate(perts_st(i)%pdim)
deallocate(perts_st(i)%plab)
deallocate(perts_st(i)%pid)
deallocate(perts_st(i)%freq)


perts_st(i)%npert = tmp_pert%npert
allocate(perts_st(i)%pdim(tmp_pert%npert))
allocate(perts_st(i)%plab(tmp_pert%npert))
allocate(perts_st(i)%pid(tmp_pert%npert))
allocate(perts_st(i)%freq(tmp_pert%npert))


perts_st(i)%pdim = tmp_pert%pdim
perts_st(i)%plab = tmp_pert%plab
perts_st(i)%pid = tmp_pert%pid
perts_st(i)%freq = tmp_pert%freq


deallocate(tmp_pert%pdim)
deallocate(tmp_pert%plab)
deallocate(tmp_pert%pid)
deallocate(tmp_pert%freq)

end do













! ! Order perts(2:tlen) in decreasing differentiation order
! 
! do i = 2, tlen
! 
! max_order_curr = t_orders(i)%npert
! m_which = i
! 
! do j = i + 1, tlen
! 
! if (t_orders(j)%npert > max_order_curr) then
! 
! max_order_curr = t_orders(j)%npert
! m_which = j
! 
! end if
! 
! end do
! 
! ! It might be possible to do this assignment in a different way
! 
! tmp_pert%npert = perts_st(m_which)%npert
! allocate(tmp_pert%pdim(tmp_pert%npert))
! allocate(tmp_pert%plab(tmp_pert%npert))
! allocate(tmp_pert%pid(tmp_pert%npert))
! allocate(tmp_pert%freq(tmp_pert%npert))
! 
! tmp_pert%pdim = perts_st(m_which)%pdim
! tmp_pert%plab = perts_st(m_which)%plab
! tmp_pert%pid = perts_st(m_which)%pid
! tmp_pert%freq = perts_st(m_which)%freq
! 
! deallocate(perts_st(m_which)%pdim)
! deallocate(perts_st(m_which)%plab)
! deallocate(perts_st(m_which)%pid)
! deallocate(perts_st(m_which)%freq)
! 
! perts_st(m_which)%npert = perts_st(i)%npert
! allocate(perts_st(m_which)%pdim(perts_st(i)%npert))
! allocate(perts_st(m_which)%plab(perts_st(i)%npert))
! allocate(perts_st(m_which)%pid(perts_st(i)%npert))
! allocate(perts_st(m_which)%freq(perts_st(i)%npert))
! 
! perts_st(m_which)%pdim = perts_st(i)%pdim
! perts_st(m_which)%plab = perts_st(i)%plab
! perts_st(m_which)%pid = perts_st(i)%pid
! perts_st(m_which)%freq = perts_st(i)%freq
! 
! deallocate(perts_st(i)%pdim)
! deallocate(perts_st(i)%plab)
! deallocate(perts_st(i)%pid)
! deallocate(perts_st(i)%freq)
! 
! 
! perts_st(i)%npert = tmp_pert%npert
! allocate(perts_st(i)%pdim(tmp_pert%npert))
! allocate(perts_st(i)%plab(tmp_pert%npert))
! allocate(perts_st(i)%pid(tmp_pert%npert))
! allocate(perts_st(i)%freq(tmp_pert%npert))
! 
! 
! perts_st(i)%pdim = tmp_pert%pdim
! perts_st(i)%plab = tmp_pert%plab
! perts_st(i)%pid = tmp_pert%pid
! perts_st(i)%freq = tmp_pert%freq
! 
! 
! deallocate(tmp_pert%pdim)
! deallocate(tmp_pert%plab)
! deallocate(tmp_pert%pid)
! deallocate(tmp_pert%freq)
! 
! 
! end do
! 
! 
! ! Order first element of perts
! 
! tmp_pert = pert_standardorder(perts_st(1))
! 
! perts_st(1)%pdim = tmp_pert%pdim
! perts_st(1)%plab = tmp_pert%plab
! perts_st(1)%pid = tmp_pert%pid
! perts_st(1)%freq = tmp_pert%freq
! 
! deallocate(tmp_pert%pdim)
! deallocate(tmp_pert%plab)
! deallocate(tmp_pert%pid)
! deallocate(tmp_pert%freq)
! 
! 
! ! Also, for the same order of differentiation between two elements of perts(2:tlen):
! ! Order alphabetically
! 
! t_first = 2
! t_last = 2
! 
! do while (t_last <= tlen)
! 
! len_curr = perts_st(t_first)%npert
! 
! ! DOES THIS WORK PROPERLY?
! 
! 
! if (t_last < tlen) then
! 
! do while ((perts_st(t_last)%npert == len_curr))
! 
! t_last = t_last + 1
! 
! if (t_last > tlen) exit
! 
! end do
! 
! t_last = t_last - 1
! 
! end if
! ! write (*,*) 'outside tlast loop'
! 
! do i = t_first, t_last
! 
! tmp_pert = pert_standardorder(perts_st(i))
! 
! allocate(min_plab_curr(tmp_pert%npert))
! min_plab_curr = t_orders(i)%plab
! m_which = i
! 
! deallocate(tmp_pert%pdim)
! deallocate(tmp_pert%plab)
! deallocate(tmp_pert%pid)
! deallocate(tmp_pert%freq)
! 
! do j = i + 1, t_last
! 
! tmp_pert = pert_standardorder(perts_st(j))
! 
! do k = 1, tmp_pert%npert
! 
! if (lle(perts_st(j)%plab(k), min_plab_curr(k))) then
! 
! min_plab_curr = perts_st(j)%plab
! m_which = j
! 
! end if
! 
! ! deallocate(min_plab_curr)
! 
! end do
! 
! 
! 
! 
! deallocate(tmp_pert%pdim)
! deallocate(tmp_pert%plab)
! deallocate(tmp_pert%pid)
! deallocate(tmp_pert%freq)
! 
! 
! end do
! 
! deallocate(min_plab_curr)
! 
! ! NOTE: FREQUENCY COMPARISON BETWEEN ELEMENTS OF t_orders IS MISSING
! 
! 
! ! Compare alphabetically
! ! Compare frequencies
! 
! 
! 
! ! It might be possible to do this assignment in a different way
! 
! tmp_pert%npert = perts_st(m_which)%npert
! allocate(tmp_pert%pdim(tmp_pert%npert))
! allocate(tmp_pert%plab(tmp_pert%npert))
! allocate(tmp_pert%pid(tmp_pert%npert))
! allocate(tmp_pert%freq(tmp_pert%npert))
! 
! tmp_pert%pdim = perts_st(m_which)%pdim
! tmp_pert%plab = perts_st(m_which)%plab
! tmp_pert%pid = perts_st(m_which)%pid
! tmp_pert%freq = perts_st(m_which)%freq
! 
! deallocate(perts_st(m_which)%pdim)
! deallocate(perts_st(m_which)%plab)
! deallocate(perts_st(m_which)%pid)
! deallocate(perts_st(m_which)%freq)
! 
! perts_st(m_which)%npert = perts_st(i)%npert
! allocate(perts_st(m_which)%pdim(perts_st(i)%npert))
! allocate(perts_st(m_which)%plab(perts_st(i)%npert))
! allocate(perts_st(m_which)%pid(perts_st(i)%npert))
! allocate(perts_st(m_which)%freq(perts_st(i)%npert))
! 
! perts_st(m_which)%pdim = perts_st(i)%pdim
! perts_st(m_which)%plab = perts_st(i)%plab
! perts_st(m_which)%pid = perts_st(i)%pid
! perts_st(m_which)%freq = perts_st(i)%freq
! 
! deallocate(perts_st(i)%pdim)
! deallocate(perts_st(i)%plab)
! deallocate(perts_st(i)%pid)
! deallocate(perts_st(i)%freq)
! 
! 
! perts_st(i)%npert = tmp_pert%npert
! allocate(perts_st(i)%pdim(tmp_pert%npert))
! allocate(perts_st(i)%plab(tmp_pert%npert))
! allocate(perts_st(i)%pid(tmp_pert%npert))
! allocate(perts_st(i)%freq(tmp_pert%npert))
! 
! 
! perts_st(i)%pdim = tmp_pert%pdim
! perts_st(i)%plab = tmp_pert%plab
! perts_st(i)%pid = tmp_pert%pid
! perts_st(i)%freq = tmp_pert%freq
! 
! 
! deallocate(tmp_pert%pdim)
! deallocate(tmp_pert%plab)
! deallocate(tmp_pert%pid)
! deallocate(tmp_pert%freq)
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! end do
! 
! t_last = t_last + 1
! t_first = t_last
! 
! end do
! 
! 
!   write(*,*) 'finished perts_standardorder'
! 
!   do i = 1, tlen
! 
! perts_st(i) = pert_standardorder(perts_st(i))
! 
!   write(*,*) perts_st(i)%pid
! 
!   end do




! nmin = 1
! 
! 
! 
! ! if (pert%npert > 1) then
! 
! do i = nmin, pert%npert
! 
! min_curr = pert%plab(i)
! min_which = i
! 
! do j = i + 1, pert%npert
! 
! if (llt(pert%plab(j), min_curr) .EQV. .TRUE.) then
! 
! min_curr = pert%plab(j)
! min_which = j
! 
! end if
! 
! end do
! 
! tmp_pdim = pert%pdim(min_which)
! tmp_plab = pert%plab(min_which)
! tmp_pid = pert%pid(min_which)
! tmp_freq = pert%freq(min_which)
! 
! pert%pdim(min_which) = pert%pdim(i)
! pert%plab(min_which) = pert%plab(i)
! pert%pid(min_which) = pert%pid(i)
! pert%freq(min_which) = pert%freq(i)
! 
! pert%pdim(i) = tmp_pdim
! pert%plab(i) = tmp_plab
! pert%pid(i) = tmp_pid
! pert%freq(i) = tmp_freq
! 
! 
! nmin = nmin + 1
! 
! end do


end function


  ! Get next element function
  function propcache_next(inst)

    implicit none

    type(propcache) :: inst, propcache_next

    ! Return pointer to next element
    propcache_next = inst%next

  end function

function propcache_getnext(inst) result(nxt)

implicit none

type(propcache), target :: inst
type(propcache), pointer :: nxt

nxt => inst%next


end function



function perts_compare(tlen, t_orders, perts_st_order)

logical :: perts_compare, acc_compare
integer ::  tlen, i
type(perts), dimension(tlen) :: t_orders, perts_st_order
type(perts) :: ptest

perts_compare = .FALSE.
acc_compare = .TRUE.


do i = 1, tlen


ptest = pert_standardorder(t_orders(i))


!  write(*,*) ptest%pid



! write(*,*) perts_st_order(i)%pid

ptest = pert_standardorder(perts_st_order(i))



! write(*,*) ptest%pid

acc_compare = acc_compare .AND. pert_compare(pert_standardorder(t_orders(i)), &
                                pert_standardorder(perts_st_order(i)))

!  write(*,*) pert_compare(pert_standardorder(t_orders(i)), &
!                                 pert_standardorder(perts_st_order(i)))

end do

if (acc_compare .eqv. .TRUE.) then

perts_compare = .TRUE.

end if

end function



  ! Add element routine

  subroutine propcache_add(inst, tlen, perturbs, propsize, data) 

    implicit none

    integer :: tlen, propsize
    type(propcache), target :: inst
    type(propcache), pointer :: new_elem
    type(propcache), pointer :: new_elem_ptr
    type(propcache), pointer :: nxt
    type(perts), dimension(tlen) :: perturbs
    complex(8), dimension(propsize) :: data

nxt => inst

!  write(*,*) 'In sdf_add'

    ! Make new instance with pert and data
    ! That instance%last = .TRUE.
  allocate(new_elem)

    call propcache_init(new_elem, tlen, perturbs, propsize, data)

new_elem_ptr => new_elem

!  write(*,*) 'got back from sdf_init'

    ! Potentially non-terminating
    ! Could this be done in another way?
    ! Traverse list until last element
    do while (nxt%last .eqv. .FALSE.)
!        write(*,*) 'skipped one ahead'
       nxt => propcache_getnext(nxt)
    end do

!  write(*,*) 'nexted to last'

    ! Set that element%last = 0
    nxt%last = .FALSE.

!  write(*,*) 'removed that last'

    ! Point from the new element to the next of the former last element 
    new_elem%next => nxt%next

!  write(*,*) 'Assigned new element next'
    ! Point from the former last element to the new last element
    nxt%next => new_elem

!  write(*,*) 'Attached new element to list'
!  write(*,*) inst%next%perturb%npert
! write(*,*) inst%next%perturb%plab
! write(*,*) ' is it last'
! write(*,*) inst%next%last
! write(*,*) ' is it last done'
! write(*,*) inst%next%next%perturb%npert
! write(*,*) inst%next%next%perturb%plab
! write(*,*) inst%next%next%last
! write(*,*) inst%next%next%next%perturb%npert
! write(*,*) inst%next%next%next%perturb%plab
! write(*,*) inst%next%next%next%last



  end subroutine





  ! Is that element already calculated?
  function propcache_already(inst, tlen, perturbs)

    implicit none

    logical :: propcache_already
    integer :: passedlast, tlen
    type(propcache), target :: inst
    type(propcache), pointer :: nxt
    type(perts), dimension(tlen) :: perturbs, perts_st_order
    

    nxt => inst

    passedlast = 0

! write(*,*) 'doing standard order'
    
    perts_st_order = perts_standardorder(tlen, perturbs)

! write(*,*) 'returned from standard order'

    propcache_already = .FALSE.

!     write(*,*) 'got to propcache_already'

    ! Potentially non-terminating
    ! Could this be done in another way?
    do while ((passedlast < 2) .AND. (propcache_already .eqv. .FALSE.))

!        write(*,*) 'getting next'
! write(*,*) 'prev npert', nxt%perturb%npert

! write(*,*) 'pointing to next'

       nxt => propcache_getnext(nxt)

! write(*,*) 'pointed to next'

! write(*,*) 'new npert', nxt%perturb%npert

if (nxt%tlen == tlen) then

!  write(*,*) 'comparing perts'

       propcache_already = perts_compare(tlen, perts_standardorder(nxt%tlen, &
                          nxt%t_orders), perts_st_order)

!  write(*,*) 'compared perts', propcache_already


end if

!  write(*,*) 'compared perts'

       if (nxt%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (propcache_already .EQV. .TRUE.) then

        write(*,*) 'propcache_already: Found element in cache'

    else
        write(*,*) 'propcache_already: Element not in cache'

    end if

  end function



  ! Find property contribution corresponding to perts and return data or fail
  ! Returns ONE matrix
  function propcache_getdata(inst, tlen, perturbs, propsize)

    implicit none

    logical :: found
    integer :: i, first, last, passedlast, tlen, propsize
    type(propcache), target :: inst
    type(propcache), pointer :: nxt
    type(perts), dimension(tlen) :: perturbs
    complex(8), dimension(propsize) :: propcache_getdata
    

!     ! Skip through elements to find first and last element to return
!     ! Does Fortran store elements in a different way and, if so, does that apply here?
! 
!     nxt => inst
! 
!     first = 1
! 
!     do i = 1, pert%npert
!        first = first + product(pert%pdim(i:pert%npert))*ind(i)
!  
!        if (i == pert%npert) then
!           last = first + ind(i)
!        end if
! 
!     end do
! 
! !     allocate(sdf_getdata(last - first + 1))
! 
!     ! Find the correct element or fail
! 
!     passedlast = 0
! 
!     ! Potentially non-terminating
!     ! Could this be done in another way?
    do while ((passedlast < 2) .OR. (propcache_already(nxt, tlen, perturbs) .eqv. .FALSE.))


         nxt => inst%next
!        inst = inst%next

       found = perts_compare(tlen, nxt%t_orders, perturbs)

       if (nxt%last .eqv. .TRUE.) then
          passedlast = passedlast +1
       end if

    end do

    ! Returns the data in the found instance or quits
! 
    if (found .eqv. .TRUE.) then

     write(*,*) 'Getting propcache data' 
! Uncomment next line when more functionality is added
     !  propcache_getdata = inst%data(last)

    else

       write(*,*) 'Failed to retrieve data in propcache_getdata: Element not found'

    end if

  end function































! End propcache linked list manipulation/data retrieval routines


! Public block



  ! Define SDF linked list manipulation/data retrieval routines

  ! Initialization routine
  ! Adds non-perturbed matrix (collapsed) as element
  subroutine sdf_init(new_elem, pert, data)

    implicit none

    type(SDF) :: new_elem
    type(perts) :: pert

integer :: data
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!    type(matrix), dimension(product(pert%pdim)) :: data




    new_elem%last = .TRUE.
    new_elem%perturb%npert = pert%npert

    allocate(new_elem%perturb%pdim(pert%npert))
    allocate(new_elem%perturb%plab(pert%npert))
    allocate(new_elem%perturb%pid(pert%npert))
    allocate(new_elem%perturb%freq(pert%npert))
    
    new_elem%perturb%pdim = pert%pdim
new_elem%perturb%plab = pert%plab
new_elem%perturb%pid = pert%pid
new_elem%perturb%freq = pert%freq
! write(*,*) new_elem%perturb%npert
    new_elem%data = data

  end subroutine


  subroutine sdf_standardorder(pert, f_data, d_data, s_data)

implicit none


type(perts) :: pert
integer :: f_data, d_data, s_data
integer :: i, j,  min_which, nmin
integer :: tmp_pdim, tmp_pid
character(4) :: tmp_plab, min_curr
complex(8) :: tmp_freq
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!type(matrix), dimension(product(pert%pdim))  :: f_data, d_data, s_data

nmin = 1



! if (pert%npert > 1) then

do i = nmin, pert%npert

min_curr = pert%plab(i)
min_which = i

do j = i + 1, pert%npert

if (llt(pert%plab(j), min_curr) .EQV. .TRUE.) then

min_curr = pert%plab(j)
min_which = j

end if

end do

tmp_pdim = pert%pdim(min_which)
tmp_plab = pert%plab(min_which)
tmp_pid = pert%pid(min_which)
tmp_freq = pert%freq(min_which)

pert%pdim(min_which) = pert%pdim(i)
pert%plab(min_which) = pert%plab(i)
pert%pid(min_which) = pert%pid(i)
pert%freq(min_which) = pert%freq(i)

pert%pdim(i) = tmp_pdim
pert%plab(i) = tmp_plab
pert%pid(i) = tmp_pid
pert%freq(i) = tmp_freq


nmin = nmin + 1

end do


! else

! ! Do nothing

! end if

end subroutine


  ! Get next element function
  function sdf_next(inst)

    implicit none

    type(SDF) :: inst, sdf_next

    ! Return pointer to next element
    sdf_next = inst%next

  end function

function sdf_getnext(inst) result(nxt)

implicit none

type(SDF), target :: inst
type(SDF), pointer :: nxt

nxt => inst%next


end function


  ! Add element routine
  ! This routine assumes that the pert and data is already in standard order

  subroutine sdf_add(inst, pert, data) 

    implicit none

    type(SDF), target :: inst
    type(SDF), pointer :: new_elem
!     type(SDF), target :: new_new_elem
    type(SDF), pointer :: new_elem_ptr
    type(SDF), pointer :: nxt
    type(perts) :: pert
integer :: data
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
 !   type(matrix), dimension(product(pert%pdim)) :: data

nxt => inst

!  write(*,*) 'In sdf_add'

    ! Make new instance with pert and data
    ! That instance%last = .TRUE.
  allocate(new_elem)

    call sdf_init(new_elem, pert, data)

new_elem_ptr => new_elem

!  write(*,*) 'got back from sdf_init'

    ! Potentially non-terminating
    ! Could this be done in another way?
    ! Traverse list until last element
    do while (nxt%last .eqv. .FALSE.)
!        write(*,*) 'skipped one ahead'
       nxt => sdf_getnext(nxt)
    end do

!  write(*,*) 'nexted to last'

    ! Set that element%last = 0
    nxt%last = .FALSE.

!  write(*,*) 'removed that last'

    ! Point from the new element to the next of the former last element 
    new_elem%next => nxt%next

!  write(*,*) 'Assigned new element next'
    ! Point from the former last element to the new last element
    nxt%next => new_elem

!  write(*,*) 'Attached new element to list'
!  write(*,*) inst%next%perturb%npert
! write(*,*) inst%next%perturb%plab
! write(*,*) ' is it last'
! write(*,*) inst%next%last
! write(*,*) ' is it last done'
! write(*,*) inst%next%next%perturb%npert
! write(*,*) inst%next%next%perturb%plab
! write(*,*) inst%next%next%last
! write(*,*) inst%next%next%next%perturb%npert
! write(*,*) inst%next%next%next%perturb%plab
! write(*,*) inst%next%next%next%last



  end subroutine


  function plab_comp(npert, plab1, plab2)

    implicit none

    logical :: plab_comp
    integer :: npert, i, j
    character(4) :: plab1(npert), plab2(npert)

! write(*,*) 'in plab_comp'

    plab_comp = .TRUE.

    do i = 1, npert
! write(*,*) 'i is', i

       if (.NOT.(plab1(i) == plab2(i))) then
          plab_comp = .FALSE.
       end if
    end do

  end function

  function pfreq_comp(n, p1, p2)

    implicit none

    logical :: pfreq_comp
    integer :: i, n
    complex(8), dimension(n) :: p1, p2
    
    pfreq_comp = .TRUE.

    do i = 1, n

       if ((p1(i) == p2(i)) .EQV. .FALSE.) then

          pfreq_comp = .FALSE.

       end if

    end do



  end function














  function pert_compare(p1, p2)

    implicit none

    logical :: pert_compare
    type(perts) :: p1, p2

!      write(*,*) 'comparing perturbations'
!  write(*,*) 'p2'
! write(*,*) p2%npert
! write(*,*) p2%pdim
! write(*,*) p2%plab
! write(*,*) p2%freq
!  write(*,*) 'p1'
! write(*,*) p1%npert
! write(*,*) p1%pdim
! write(*,*) p1%plab
! write(*,*) p1%freq


    ! Compare number of perturbations
    if (p1%npert == p2%npert) then

! write(*,*) 'npert was equal'
       ! Compare perturbation labels
       if (plab_comp(p1%npert, p1%plab, p2%plab) .eqv. .TRUE.) then
!  write(*,*) 'plab was equal'

          ! Compare perturbation frequencies
          ! Must compare element-wise or is this OK?
          if (pfreq_comp(p1%npert, p1%freq, p2%freq) .eqv. .TRUE.) then
! write(*,*) 'pfreq was equal'

             ! Pseudo for comparing other perturbation info
             ! if (p1%other == p2%other) then
             pert_compare = .TRUE.
             ! else
             ! pert_compare = .FALSE.
             ! end if

          else 
!              write(*,*) 'freq false', p1%freq, p2%freq

             pert_compare = .FALSE.
          end if
       else 

!        write(*,*) 'plab false', p1%plab, p2%plab
          pert_compare = .FALSE.
       end if
    else 

!        write(*,*) 'npert false', p1%npert, p2%npert
       pert_compare = .FALSE.
    end if

  end function





  ! Is that element already calculated?
  function sdf_already(inst, pert)

    implicit none

    logical :: sdf_already
    type(SDF), target :: inst
    type(SDF), pointer :: nxt
    type(perts) :: pert, pert_st_order
    integer :: passedlast

    nxt => inst

    passedlast = 0
    
    pert_st_order = pert_standardorder(pert)

    sdf_already = .FALSE.

!     write(*,*) 'got to sdf_already'

    ! Potentially non-terminating
    ! Could this be done in another way?
    do while ((passedlast < 2) .AND. (sdf_already .eqv. .FALSE.))

!        write(*,*) 'getting next'
! write(*,*) 'prev npert', nxt%perturb%npert

       nxt => sdf_getnext(nxt)

! write(*,*) 'new npert', nxt%perturb%npert

       sdf_already = pert_compare(nxt%perturb, pert_st_order)

! write(*,*) 'compared perts'

       if (nxt%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (sdf_already .EQV. .TRUE.) then

!        write(*,*) 'sdf_already: Found element in cache'

    else
!        write(*,*) 'sdf_already: Element not in cache'

    end if

  end function

  ! Find element ind corresponding to pert and return data or fail
  ! Returns ONE matrix
  function sdf_getdata(inst, pert, ind)

    implicit none

    logical :: found
    type(SDF), target :: inst
    type(SDF), pointer :: nxt
    type(perts) :: pert
    type(matrix) :: sdf_getdata
    integer, dimension(pert%npert) :: ind
    integer :: i, first, last, passedlast

    ! Skip through elements to find first and last element to return
    ! Does Fortran store elements in a different way and, if so, does that apply here?

    nxt => inst

    first = 1

    do i = 1, pert%npert
       first = first + product(pert%pdim(i:pert%npert))*ind(i)
 
       if (i == pert%npert) then
          last = first + ind(i)
       end if

    end do

!     allocate(sdf_getdata(last - first + 1))

    ! Find the correct element or fail

    passedlast = 0

    ! Potentially non-terminating
    ! Could this be done in another way?
    do while ((passedlast < 2) .OR. (sdf_already(nxt, pert) .eqv. .FALSE.))


         nxt => inst%next
!        inst = inst%next

       found = pert_compare(nxt%perturb, pert)

       if (nxt%last .eqv. .TRUE.) then
          passedlast = passedlast +1
       end if

    end do

    ! Returns the data in the found instance or quits

    if (found .eqv. .TRUE.) then

     write(*,*) 'Getting SDF data', pert%pid
! Uncomment next line when more functionality is added
     !  sdf_getdata = inst%data(last)

    else

       write(*,*) 'Failed to retrieve data in sdf_getdata: Element not found'

    end if

  end function

  ! Find out if kn rules say that this term should be skipped
  function kn_skip(npert, pertid, kn)

    implicit none

    logical :: kn_skip, pert_hasfirst
    integer :: npert, i
    integer, dimension(npert) :: pertid
    integer, dimension(2) :: kn

! write(*,*) 'inside kn_skip', size(pertid)

    kn_skip = .FALSE.


    pert_hasfirst = .FALSE.

    do i = 1, size(pertid)

! write(*,*) 'inside kn_skip i'
       if (pertid(i) == 1) then
          pert_hasfirst = .TRUE.
       end if
    end do

! write(*,*) 'inside kn_skip'
   
    if (pert_hasfirst .eqv. .TRUE.) then

       if (kn(1) < size(pertid)) then

          kn_skip = .TRUE.
       end if
! MR: I think this should be left out, but it is kept in comments if this turns out to be wrong
!
!        if (kn(1) < (size(pertid) - 1)) then
!           kn_skip = .TRUE.
!        end if

    ! Only n can require skip now
    else

       if (kn(2) < size(pertid)) then

          kn_skip = .TRUE.
       end if

    end if

  end function


  ! Find out if lagrangian kn rules say that this term should be skipped
! NOTE: THIS FUNCTION IS JUST A COPY OF kn_skip RIGHT NOW
! FIX THIS

  function kn_skiplag(npert, pertid, kn)

    implicit none

    logical :: kn_skiplag, pert_hasfirst
    integer :: npert, i
    integer, dimension(npert) :: pertid
    integer, dimension(2) :: kn


! write(*,*) 'inside kn_skiplag'

    kn_skiplag = .FALSE.


    pert_hasfirst = .FALSE.

    do i = 1, size(pertid)
       if (pertid(i) == 1) then
          pert_hasfirst = .TRUE.
       end if
    end do
   
    if (pert_hasfirst .eqv. .TRUE.) then

       if (kn(1) < size(pertid)) then

          kn_skiplag = .TRUE.
       end if
! MR: I think this should be left out, but it is kept in comments if this turns out to be wrong
!
!        if (kn(1) < (size(pertid) - 1)) then
!           kn_skip = .TRUE.
!        end if

    ! Only n can require skip now
    else

       if (kn(2) < size(pertid)) then

          kn_skiplag = .TRUE.
       end if

    end if

  end function







  function get_ncarray(total_order, tlen, t_orders)

    implicit none

    integer :: total_order, tlen, i, j, k
    integer, dimension(total_order) :: get_ncarray
    type(perts), dimension(tlen) :: t_orders

    do i = 1, total_order
       do j = 1, tlen
          do k = 1, t_orders(j)%npert

             if (t_orders(j)%pid(k) == i) then
                get_ncarray(i) = t_orders(j)%pdim(k)
             end if

          end do
       end do

    end do

  end function

  function nc_only(total_order, thisorder, tlen, t_orders, ncarray)

    implicit none

    integer :: i, j, total_order, thisorder, tlen
    integer, dimension(total_order) :: ncarray
    integer, dimension(total_order) :: nc_only
    type(perts), dimension(tlen) :: t_orders

    do i = 1, size(ncarray)
       nc_only(i) = 1
    end do

    do i = 1, tlen
       do j = 1, t_orders(i)%npert
          nc_only(t_orders(i)%pid(j)) = ncarray(t_orders(i)%pid(j))
       end do
    end do

  end function

  function nc_onlysmall(total_order, thisorder, tlen, t_orders, ncarray)

    implicit none

    integer :: i, j, k, total_order, thisorder, tlen
    integer, dimension(total_order) :: ncarray
    integer, dimension(thisorder) :: nc_onlysmall
    type(perts), dimension(tlen) :: t_orders

    k = 1

    do i = 1, tlen
       do j = 1, t_orders(i)%npert

          nc_onlysmall(k) = ncarray(t_orders(i)%pid(j))
          k = k + 1

       end do
    end do

  end function



  recursive subroutine make_outerindices(tot_outer, lvl, ncarray, offset, outer_indices)

    implicit none

    integer :: i, j, k, tot_outer, lvl, offset
    integer, dimension(tot_outer) :: ncarray
    integer, dimension(product(ncarray), tot_outer) :: outer_indices

    k = 1

    if (tot_outer > 0) then
       do i = 1, ncarray(lvl)

          if (lvl < tot_outer) then

             call make_outerindices(tot_outer, lvl + 1, ncarray, &
             k + offset - 1, outer_indices)

          end if

          if (lvl <= tot_outer) then

             do j = 1, product(ncarray(lvl:size(ncarray)))/ncarray(lvl)

                outer_indices(k + offset, lvl) = i
                k = k + 1

             end do

          end if

       end do

    else

    end if

  end subroutine


  function make_outerwhichpert(totpert, tlen, t_orders)

    implicit none

    integer :: i, j, k, totpert, tlen
    type(perts), dimension(tlen) :: t_orders
    integer, dimension(totpert) :: make_outerwhichpert

    do i = 1, totpert

       make_outerwhichpert(i) = 0

    end do

    k = 1

    do i = 2, tlen
       do j = 1, t_orders(i)%npert

          make_outerwhichpert(t_orders(i)%pid(j)) = k
          k = k + 1

       end do
    end do

  end function


  function make_outerwhichpertbig(totpert, tlen, t_orders)

    implicit none

    integer :: i, j, k, totpert, tlen
    type(perts), dimension(tlen) :: t_orders
    integer, dimension(totpert) :: make_outerwhichpertbig

    do i = 1, totpert

       make_outerwhichpertbig(i) = 0

    end do

    k = 1

    do i = 1, tlen
       do j = 1, t_orders(i)%npert

          make_outerwhichpertbig(t_orders(i)%pid(j)) = k
          k = k + 1

       end do
    end do

  end function


  subroutine get_prop_offsets(totpert, tlen, tot_outer, ncarray, ncinner, &
  t_orders, outer_index, o_which, inner_offsets)

    implicit none

    logical :: isinner
    integer :: i, j, k, totpert, tot_outer, outerprod, tlen
    integer, dimension(totpert) :: ncarray, ncprod
    integer, dimension(totpert) :: ncinner, o_which, ncwhich
    logical, dimension(totpert) :: isouter
    integer, dimension(tot_outer) :: outer_index
    type(perts), dimension(tlen) :: t_orders
    integer, dimension(product(ncinner)) :: inner_offsets
    integer, dimension(product(ncinner), totpert - tot_outer) :: inner_indices

    call make_outerindices(t_orders(1)%npert, 1, t_orders(1)%pdim, 0, inner_indices)

    outerprod = 0
    k = 0

    do i = 1, size(ncarray)

       isinner = .FALSE.

       ncprod(i) = product(ncarray(i:size(ncarray)))/ncarray(i)

       do j = 1, t_orders(1)%npert
          if (t_orders(1)%pid(j) == i) then

             isinner = .TRUE.

          end if
       end do

       if (isinner .eqv. .FALSE.) then

          isouter(i) = .TRUE.

       else 

          isouter(i) = .FALSE.

       end if

    end do

    do i = 1, size(isouter)
       if (isouter(i) .eqv. .TRUE.) then

          k = k + 1
          outerprod = outerprod + ncprod(i)*(outer_index(k) - 1)

       end if
    end do

    do i = 1, size(inner_indices, 1)
 
       inner_offsets(i) = outerprod
 
       do j = 1, size(inner_indices, 2)

          inner_offsets(i) = inner_offsets(i) + (inner_indices(i,j) - 1) * &
          ncprod(t_orders(1)%pid(j))
 
       end do

       inner_offsets(i) = inner_offsets(i) + 1

    end do

  end subroutine


  function get_pidoutersmall(totouter, len_outer, o_orders)

    implicit none

    integer :: totouter, len_outer, i, j, k
    integer, dimension(totouter) :: get_pidoutersmall
    type(perts), dimension(len_outer) :: o_orders

    k = 1

    do i = 1, len_outer
       do j = 1, o_orders(i)%npert

          get_pidoutersmall(k) = o_orders(i)%pid(j)
          k = k + 1

       end do
    end do

  end function


  subroutine sortdimbypid(totpert, totouter, pids, dims, dimsouter, whichs)

    implicit none

    ! DIMENSIONS ARE "SMALL" (NOT INCLUDING INNER INDICES)

    integer :: totouter, totpert, s, i, j, whichmax, whatmax
    integer, dimension(totouter) :: b, d, pids, dimsouter
    integer, dimension(totpert) :: whichs, dims

    do i = 1, totpert

       whichs(i) = 0

    end do

    s = totouter
    j = totouter
    d = pids

    do while (j > 0)

       whatmax = 0

       ! At which index is the pid largest?

       do i = 1, s
          if (d(i) > whatmax) then

             ! It is currently largest at index i
             whatmax = d(i)
             whichmax = i

          end if
       end do

       ! Then, put the dimension of that pid at the current end of the array to be returned

       b(j) = dims(whatmax)

       ! j is the (current) highest outer index

       whichs(j) = whatmax
       j = j - 1
       d(whichmax) = 0

    end do

    dimsouter = b

  end subroutine

! Some explanations
! It is assumed that the prop array will be collapsed
! The tmp array is shaped into its proper rank before being sent to the
! various average/nucpot routines
! The outer indices are the indices that are not "fully filled" in the
! average/nucpot calls - The outer indices are the indices that do not
! correspond to perturbations with respect to which the energy is
! differentiated. They are the indices of differentiation of the
! density matrices.
! The inner offsets are the memory positions in the prop tensor where the
! elements of the "recollapsed" tmp tensor (in regular order) should be put
! They correspond to the memory positions of all the inner indices (the 
! indices that correspond to perturbations with respect to which the energy
! is differentiated), but offset by some particular outer indices in the
! outermost loop.

  subroutine get_energy(mol, tlen, totpert, t_orders, d_order, D, propsize, &
                        energy_cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts), dimension(tlen) :: t_orders
    type(SDF) :: D
    type(propcache) :: energy_cache
    integer :: i, j, tlen, totpert, d_order, propsize
    integer, dimension(totpert) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, pidoutersmall
    integer, allocatable, dimension(:,:) :: outer_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(propsize) :: prop
    complex(8), dimension(propsize) :: prop_forcache

    ncarray = get_ncarray(totpert, tlen, t_orders)
    ncouter = nc_only(totpert, totpert - t_orders(1)%npert, tlen - 1, &
              t_orders(2:tlen), ncarray)
    ncinner = nc_only(totpert, t_orders(1)%npert, 1, t_orders(1), ncarray)

    allocate(ncoutersmall(totpert - t_orders(1)%npert))
    allocate(pidoutersmall(totpert - t_orders(1)%npert))

    ncoutersmall = nc_onlysmall(totpert, totpert - t_orders(1)%npert, tlen - 1, &
                   t_orders(2:tlen), ncarray)
    pidoutersmall = get_pidoutersmall(totpert - t_orders(1)%npert, tlen - 1, &
                    t_orders(2:tlen))

    allocate(o_whichpert(totpert))
    allocate(o_whichpertbig(totpert))
    allocate(outer_indices(product(ncoutersmall),size(ncoutersmall)))
    allocate(tmp(product(ncarray)))

    ! Make array where entry %pid is the rank of the outer index under consideration

    o_whichpert = make_outerwhichpert(totpert, tlen, t_orders)
    o_whichpertbig = make_outerwhichpertbig(totpert, tlen, t_orders)

    call sortdimbypid(totpert, totpert - t_orders(1)%npert, pidoutersmall, &
         ncarray, ncoutersmall, o_whichpert)

    ! FIX RESHAPE LATER. MAYBE DO COLLAPSED ARRAYS ALL THE WAY
    ! tmp = reshape(tmp, ncarray)

    if (totpert > t_orders(1)%npert) then

       ! Make outer indices

       call make_outerindices(totpert - t_orders(1)%npert, 1, ncoutersmall, 0, outer_indices)

       ! Make another array with the inverse information? (the rank of the outer index under
       ! consideration would be %pid)

       allocate(inner_offsets(product(ncinner)))

!         write(*,*), size(outer_indices, 1)

       do i = 1, size(outer_indices, 1)

!           write(*,*) i

          tmp = 0.0

          if (d_order == 0) then

             ! write(*,*) 'Calling contributions for d_order = 0'

             ! Is this a good place to put the nucpot contribution?
             ! call rsp_nucpot(mol, t_orders(1)%npert, t_orders(1)%plab, (0d0,0d0)*(/0, j = 1, &
             !                 t_orders(1)*%npert/), (/0, j = 1, t_orders(1)*%npert/), nc, tmp) 
             ! call rsp_oneave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, D, tmp)
             ! call rsp_twoave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/), nc, &
             !                 (/ sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )
             ! call rsp_excave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, &
             !                 (/ sdf_getdata(D, emptypert(1), (//)),
             !                 sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )

          elseif (d_order == 1) then

             ! write(*,*) 'Calling contributions for d_order = 1'

             ! call rsp_oneave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, &
             !                 sdf_getdata(D, t_orders(2), outer_indices(i,:)), tmp)
             ! ! NOT SURE IF THE INDICES IN THIS CALL IS CORRECT
             ! call rsp_twoave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/), nc, &
             !                 (/ sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )
             ! call rsp_excave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, &
             !                 (/ sdf_getdata(D, emptypert(1), (//)),
             !                 sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )

          elseif (d_order == 2) then

             ! write(*,*) 'Calling contributions for d_order = 2'

             ! call rsp_twoave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/), nc, &
             !                 (/ sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )
             ! call rsp_excave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, &
             !                 (/ sdf_getdata(D, emptypert(1), (//)),
             !                 sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )

          else

             ! write(*,*) 'Calling contributions for d_order = 3 or higher: ', d_order

             ! Is it possible to combine explicit and implicit within 
             ! the same array declaration as in the sdf_getdata calls below?
             ! call rsp_excave(mol, t_orders(1)%npert, t_orders(1)%plab, &
             !                 (/1, j = 1, t_orders(1)*%npert/)), nc, &
             !                 (/ sdf_getdata(D, emptypert(1), (//)),
             !                 sdf_getdata(D, t_orders(k), (/ &
             !                 outer_indices(i,o_whichpert(t_orders(k)%pid(j))), &
             !                 j = 1, t_orders(k)%npert /)), k = 2, tlen /) )

          end if

          ! Get inner offsets

          call get_prop_offsets(totpert, tlen, totpert - t_orders(1)%npert, &
               ncarray, ncinner, t_orders,outer_indices(i,:), &
               o_whichpertbig, inner_offsets)

          do j = 1, size(inner_offsets)

             prop(inner_offsets(j)) = prop(inner_offsets(j)) + tmp(j)
             prop_forcache (inner_offsets(j)) = prop_forcache(inner_offsets(j)) + tmp(j)

          end do

       end do

        deallocate(inner_offsets)
    else

       write(*,*) 'all indices inner'

    end if

    call propcache_add(energy_cache, tlen, t_orders, propsize, prop_forcache)    


    deallocate(ncoutersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_whichpertbig)
    deallocate(outer_indices)
    deallocate(tmp)


  end subroutine

  ! Calculate and add all the energy contributions

  recursive subroutine rsp_ener(mol, pert, totpert, kn, tlen, t_orders, &
                       d_order, D, propsize, energy_cache, prop)

    implicit none

    logical :: e_knskip
    type(rsp_cfg) :: mol
    type(perts) :: pert
    integer, dimension(2) :: kn
    integer :: tlen, d_order, i, j, totpert, propsize
    type(perts), dimension(tlen) :: t_orders, t_new
    type(SDF) :: D
    type(propcache) :: energy_cache
    complex(8), dimension(propsize) :: prop

    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'

    if (t_orders(1)%npert == 0) then

       call rsp_ener(mol, pert_rf(pert), totpert, kn, tlen, &
       (/pert_getone(pert,1), t_orders(2:size(t_orders))/), &
       d_order, D, propsize, energy_cache, prop)

    else


       call rsp_ener(mol, pert_rf(pert), totpert, kn, tlen, &
       (/pert_ext(t_orders(1), pert_getone(pert,1)), t_orders(2:size(t_orders))/), &
       d_order, D, propsize, energy_cache, prop)

    end if
    
       ! 2. Differentiate all of the contraction densities in turn

       ! Find the number of terms

       do i = 2, tlen

          t_new = t_orders

          if (t_orders(i)%npert == 0) then

             t_new(i) = pert_getone(pert, 1)

          else

             t_new(i) = pert_ext(t_new(i), pert_getone(pert, 1))

          end if

          call rsp_ener(mol, pert_rf(pert), totpert, kn, tlen, &
          t_new, d_order + 1, D, propsize, energy_cache, prop)

       end do


       ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_ener(mol, pert_rf(pert), totpert, kn, tlen + 1, &
       (/t_orders(:), pert_getone(pert, 1)/), &
       d_order + 1, D, propsize, energy_cache, prop)





! Then (at the lowest recursion level): Call another subroutine 
! eval_eterm (or some other name)
! It evaluates each of the E terms found by the rsp_ener recursion
! It determines whether to call oneave, twoave, excave 
! However, before that, check to see if evaluated before 'up to isomorphism'
! If so: Just add that property instead (but at the right "ways", e.g. (:,i,j,:))
! Should that check happen inside eval_eterm or before? Inside for now

    else


    e_knskip = .FALSE.


       write(*,*) 'getting energy contrib'



       do i = 1, tlen
 
          write(*,*) 'i', i, ' and pid', t_orders(i)%pid

          if (i > 1) then

             if(kn_skip(t_orders(i)%npert, t_orders(i)%pid, kn) .EQV. .TRUE.) then

                e_knskip = .TRUE.

             end if
          
          end if


       end do

       if (e_knskip .EQV. .FALSE.) then

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

write(*,*) 'Evaluating propcache_already'

if (propcache_already(energy_cache, tlen, t_orders) .EQV. .TRUE.) then

open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

write(*,*) 'Getting values from cache'
write(*,*) ' '
       
else

          call get_energy(mol, tlen, totpert, t_orders, d_order, D, propsize, energy_cache, prop)

          write(*,*) 'Calculated energy contribution'
          write(*,*) ' '


end if

       else

          write(*,*) 'energy contribution was k-n skipped'
          write(*,*) ' '

       end if



    end if


  end subroutine




  recursive function d_superstr_dry(mol, pert, kn, primed, d_curr) result(d_size)

    implicit none

    logical :: primed
    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(3) :: d_curr
    integer, dimension(2) :: kn
    integer :: i, d_size

!  write(*,*) 'just started d_superstr_dry'

! write(*,*) 'perts', pert%pid
! write(*,*) 'd curr', d_curr(1)%pid, d_curr(2)%pid, d_curr(3)%pid

d_size = 0

! write(*,*) 'inside superstr_dry'

    if (pert%npert > 0) then

       ! Removed loop because it seems to (redundantly) duplicate recursive calls
       ! do i = 1, 3

!  write(*,*) 'recursion 1, d_size is now', d_size

       d_size = d_size + d_superstr_dry(mol, pert_rf(pert), kn, primed, &
                         (/pert_ext(d_curr(1), pert_getone(pert, 1)) , d_curr(2:3)/))
!  write(*,*) 'recursion 2, d_size is now', d_size

       d_size = d_size + d_superstr_dry(mol, pert_rf(pert), kn, primed, &
                         (/d_curr(1), pert_ext(d_curr(2), pert_getone(pert, 1)), &
                         d_curr(3)/))
!  write(*,*) 'recursion 3, d_size is now', d_size

       d_size = d_size + d_superstr_dry(mol, pert_rf(pert), kn, primed, &
                         (/d_curr(1:2), pert_ext(d_curr(3), pert_getone(pert, 1))/))

! write(*,*) 'after recursion 3, d_size is now', d_size

       ! end do

    else

!        write(*,*) 'stage 4'

       if (primed .EQV. .TRUE.) then

!           write(*,*) 'These are the pids (primed):' 
!           write(*,*) d_curr(1)%pid
!           write(*,*) d_curr(2)%pid
!           write(*,*) d_curr(3)%pid

          if ( ( ( ( d_curr(1)%npert <= kn(2)) .AND.&
              d_curr(2)%npert <= kn(2) ) .AND. &
              d_curr(3)%npert <= kn(2) ) .eqv. .TRUE.) then

!              write(*,*) 'Added to d_size'

             d_size = 1

          else

!              write(*,*) 'Did not add to d_size'

             d_size = 0

          end if


       else


!           write(*,*) 'These are the pids (nonprimed):' 
!           write(*,*) d_curr(1)%pid
!           write(*,*) d_curr(2)%pid
!           write(*,*) d_curr(3)%pid


          if ( ( ( kn_skip(d_curr(1)%npert, d_curr(1)%pid, kn) .OR. &
                   kn_skip(d_curr(2)%npert, d_curr(2)%pid, kn) ) .OR. &
                   kn_skip(d_curr(3)%npert, d_curr(3)%pid, kn) ) .eqv. .FALSE.) then

!              write(*,*) 'Added to d_size'

             d_size = 1

          else

!              write(*,*) 'Did not add to d_size'

             d_size = 0

          end if

       end if

    end if

  end function





  recursive subroutine d_superstr(mol, pert, kn, primed, d_curr, d_size, incr, d_struct)

    implicit none

    logical :: primed
    integer :: i, d_size, incr
    integer, dimension(2) :: kn    
    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(3) :: d_curr
    type(perts), dimension(d_size, 3) :: d_struct
    
    

    if (pert%npert > 0) then

       ! Removed loop because it seems to (redundantly) duplicate recursive calls
       ! do i = 1, 3
! write(*,*) 'calling rec 1'

       call d_superstr(mol, pert_rf(pert), kn, primed, &
            (/pert_ext(d_curr(1), pert_getone(pert, 1)), d_curr(2:3)/), &
             d_size, incr, d_struct)
! write(*,*) 'calling rec 2', incr

       call d_superstr(mol, pert_rf(pert), kn, primed, &
            (/d_curr(1), pert_ext(d_curr(2), pert_getone(pert, 1)), d_curr(3)/), &
             d_size, incr, d_struct)
! write(*,*) 'calling rec 3'

       call d_superstr(mol, pert_rf(pert), kn, primed, &
            (/d_curr(1:2), pert_ext(d_curr(3), pert_getone(pert, 1))/), &
             d_size, incr, d_struct)

       ! end do

    else


       if (primed .EQV. .TRUE.) then

!  write(*,*) 'checking and possibly increasing', incr

          if ( ( ( ( d_curr(1)%npert <= kn(2)) .AND.&
              d_curr(2)%npert <= kn(2) ) .AND. &
              d_curr(3)%npert <= kn(2) ) .eqv. .TRUE.) then

!  write(*,*) 'allowed', incr
! write(*,*) 'npert', d_curr(1)%npert, d_curr(2)%npert, d_curr(3)%npert
!  write(*,*) 'pid', d_curr(1)%pid, d_curr(2)%pid, d_curr(3)%pid
             incr = incr + 1
             d_struct(incr, :) = d_curr
 
! write(*,*) 'increased', incr
              
          end if

 

       else

!  write(*,*) 'checking and increasing nonprimed', incr

          if ( ( ( kn_skip(d_curr(1)%npert, d_curr(1)%pid, kn) .OR. &
                   kn_skip(d_curr(2)%npert, d_curr(2)%pid, kn) ) .OR. &
                   kn_skip(d_curr(3)%npert, d_curr(3)%pid, kn) ) .eqv. .FALSE.) then



incr = incr + 1 

! write(*,*) 'adding to structure', trythis

! trythis = trythis + 1

! write(*,*) 'adding to structure', trythis

             d_struct(incr, :) = d_curr(:)

! write(*,*) 'added to structure'             

! incr = trythis
! write(*,*) 'added to structure 2'       
          end if

       end if

    end if

  end subroutine




function g_i_ind(ord, totpert, ind)

implicit none

type(perts) :: ord
integer :: i, totpert
integer, dimension(ord%npert) :: g_i_ind
integer, dimension(totpert) :: ind

! The dimension of ind was changed from 'tot' to totpert. It is believed that 'tot' was
! just a naming error

do i = 1, ord%npert

! Added index i to rhs below

g_i_ind(i) = ind(ord%pid(i))

end do

end function



function f_ev(d_ord)

implicit none

type(perts) :: d_ord
complex(8) :: f_ev
integer :: i

! write(*,*) 'in f_ev'

f_ev = 0.0

if (d_ord%npert > 0) then

do i = 1, d_ord%npert

f_ev = f_ev + d_ord%freq(i)

end do

end if

! write(*,*) 'f_ev done'

end function

! FOR ALL MATRIX FUNCTIONS BELOW: Made corrections to g_i_ind calls

! Get ONE W at some perturbation indices ind
! Look at again for sign errors

function rsp_getw(mol, d_size, d_ord, totpert, ind, F, D, S) result(W)

implicit none

integer :: i, totpert, d_size
type(rsp_cfg) :: mol
type(perts), dimension(d_size, 3) :: d_ord
integer, dimension(totpert) :: ind
type(sdf) :: F, D, S
type(matrix) :: W




! write(*,*) 'NOTE: GETW: UNCOMMENT MULTIPLICATIONS WHEN MORE FUNCTIONALITY IS ADDED'

do i = 1, d_size

W = W !+ sdf_getdata(D, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(F, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

! write(*,*) 'got the first part'

W = W !- ((1.0)/(2.0)) * f_ev(d_ord(i,1)) * &
!         sdf_getdata(D, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(F, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

! write(*,*) 'got the second part'

W = W !+ ((1.0)/(2.0)) * f_ev(d_ord(i,3)) * &
!         sdf_getdata(D, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(F, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

! write(*,*) 'got the third part'

end do

end function


! Get ONE Y at some perturbation indices ind
! Look at again for sign errors

function rsp_gety(mol, d_size, d_ord, totpert, ind, F, D, S) result(Y)

implicit none

integer :: i, totpert, d_size
type(rsp_cfg) :: mol
type(perts), dimension(d_size, 3) :: d_ord
integer, dimension(totpert) :: ind
type(sdf) :: F, D, S
type(matrix) :: Y

! write(*,*) 'NOTE: GETY: UNCOMMENT MULTIPLICATIONS WHEN MORE FUNCTIONALITY IS ADDED'

do i = 1, d_size

Y = Y !+ sdf_getdata(F, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

Y = Y !- sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(F, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

Y = Y !- ((1.0)/(2.0)) * f_ev(d_ord(i,3)) * &
!         sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

Y = Y !- ((1.0)/(2.0)) * f_ev(d_ord(i,1)) * &
!         sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

end do

end function


! Get ONE Z at some perturbation indices ind
! Look at again for sign errors

function rsp_getz(mol, d_size, d_ord, kn, totpert, ind, F, D, S) result(Z)

implicit none

integer :: i, totpert, d_size
type(rsp_cfg) :: mol
type(perts), dimension(d_size, 3) :: d_ord
type(perts) :: merged_pert
integer, dimension(2) :: kn
integer, dimension(totpert) :: ind
type(sdf) :: F, D, S
type(matrix) :: Z

! write(*,*) 'NOTE: GETZ: UNCOMMENT MULTIPLICATIONS WHEN MORE FUNCTIONALITY IS ADDED'

do i = 1, d_size

Z = Z !+ sdf_getdata(D, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

! write(*,*) 'merging pert', d_ord(i,1)%npert, d_ord(i,2)%npert, d_ord(i,3)%npert
merged_pert = merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3)))

! write(*,*) 'merged pert'

if (kn_skip(totpert, merged_pert%pid, kn) .eqv. .FALSE.) then

Z = Z !- sdf_getdata(D, merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), &
!         g_i_ind(merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), totpert, ind))

end if

end do

end function


! Get ONE lambda at some perturbation indices ind
! Look at again for sign errors

function rsp_getl(mol, pert_a, d_size, d_ord, totpert, ind, D, S) result(L)

implicit none

integer :: i, totpert, d_size
type(rsp_cfg) :: mol
type(perts) :: pert_a
type(perts), dimension(d_size, 3) :: d_ord
integer, dimension(totpert) :: ind
type(sdf) :: D, S
type(matrix) :: L


! write(*,*) 'NOTE: GETL: UNCOMMENT MULTIPLICATIONS WHEN MORE FUNCTIONALITY IS ADDED'

do i = 1, d_size

L = L !+ sdf_getdata(D, merge_pert(pert_a, d_ord(i,1)), &
!                     g_i_ind(merge_pert(pert_a, d_ord(i,1)), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

L = L !- sdf_getdata(D, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!         sdf_getdata(S, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!         sdf_getdata(D, merge_pert(pert_a, d_ord(i,3)), &
!                     g_i_ind(merge_pert(pert_a, d_ord(i,3)), totpert, ind))

end do

end function


! Get ONE zeta at some perturbation indices ind
! Look at again for sign errors

function rsp_getc(mol, pert_a, kn, d_size, d_ord, totpert, ind, F, D, S) result(Cmat)

implicit none

integer :: i, totpert, d_size
type(rsp_cfg) :: mol
type(perts) :: pert_a, merged_pert
type(perts), dimension(d_size, 3) :: d_ord
integer, dimension(2) :: kn
integer, dimension(totpert) :: ind
type(sdf) :: F, D, S
type(matrix) :: Cmat


! write(*,*) 'NOTE: GETC: UNCOMMENT MULTIPLICATIONS WHEN MORE FUNCTIONALITY IS ADDED'

do i = 1, d_size

! NOTE: The arguments of the f_ev calls were written after the rest of the
! code. If the results are incorrect, this would be a relevant thing to
! look at again.
! Also made: corrections to g_i_ind calls

!  Cmat = sdf_getdata(F, merge_pert(pert_a, d_ord(i,1)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,1)), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

merged_pert = merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3)))


if (kn_skip(totpert, merged_pert%pid, kn) .eqv. .FALSE.) then

! SHOULD THIS BE UNCOMMENTED? LOOK AT AGAIN
!  Cmat = Cmat - ((1.0)/(2.0)) * &
!         sdf_getdata(F, merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), &
!         g_i_ind(merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), totpert, ind))

end if

 Cmat = Cmat !- sdf_getdata(F, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, merge_pert(pert_a, d_ord(i,3)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,3)), totpert, ind))

 Cmat = Cmat !+ ((1.0)/(2.0))*f_ev(d_ord(i,1)) * &
!          sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, merge_pert(pert_a, d_ord(i,3)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,3)), totpert, ind))

 Cmat = Cmat !+ f_ev(d_ord(i,2)) * &
!          sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, merge_pert(pert_a, d_ord(i,3)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,3)), totpert, ind))

 Cmat = Cmat !+ sdf_getdata(S, d_ord(i,1), g_i_ind(d_ord(i,1), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(F, merge_pert(pert_a, d_ord(i,3)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,3)), totpert, ind))

if (kn_skip(totpert, merged_pert%pid, kn) .eqv. .FALSE.) then
! SHOULD THIS BE UNCOMMENTED? LOOK AT AGAIN 
!  Cmat = Cmat - ((1.0)/(2.0))* &
!         sdf_getdata(F, merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), &
!         g_i_ind(merge_pert(d_ord(i,1), merge_pert(d_ord(i,2), d_ord(i,3))), totpert, ind))

end if

 Cmat = Cmat !- sdf_getdata(S, merge_pert(pert_a, d_ord(i,1)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,1)), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(F, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

 Cmat = Cmat !+ ((1.0)/(2.0))*f_ev(d_ord(i,3)) * &
!          sdf_getdata(S, merge_pert(pert_a, d_ord(i,1)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,1)), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

 Cmat = Cmat !+ f_ev(d_ord(i,2)) * &
!          sdf_getdata(S, merge_pert(pert_a, d_ord(i,1)), &
!                      g_i_ind(merge_pert(pert_a, d_ord(i,1)), totpert, ind)) * &
!          sdf_getdata(D, d_ord(i,2), g_i_ind(d_ord(i,2), totpert, ind)) * &
!          sdf_getdata(S, d_ord(i,3), g_i_ind(d_ord(i,3), totpert, ind))

end do

end function












  subroutine get_pulay_kn(mol, p12, kn, F, D, S, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert, emptypert
    type(perts), dimension(2) :: p12
    type(perts), dimension(:,:), allocatable :: d_ordb
    type(SDF) :: S, D, F
    type(propcache) :: cache
    type(matrix) :: W
    integer :: i, j, sstr_incr
    integer :: propsize, d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets
    integer, allocatable, dimension(:,:) :: outer_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(propsize) :: prop
    complex(8), dimension(propsize) :: prop_forcache

! write(*,*) 'getting supsize'

! write(*,*) 'p12(1)', p12(1)%pid
! write(*,*) p12(2)%plab
! write(*,*) p12(2)%freq


! write(*,*) 'p12(2)', p12(2)%pid
! write(*,*) p12(2)%plab
! write(*,*) p12(2)%freq

    d_supsize = d_superstr_dry(mol, p12(2), kn, .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

!    write(*,*) 'supsize', d_supsize

    allocate(d_ordb(d_supsize, 3))

! write(*,*) 'getting superstr'
    ! Get superstructure for W

sstr_incr = 0

    call d_superstr(mol, p12(2), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, sstr_incr, d_ordb)

! IT MIGHT BE NECESSARY TO ZERO ANY TMP ARRAY BETWEEN OUTER INDEX TUPLES
! THIS COULD ALSO BE TRUE FOR OTHER TYPES OF CONTRIBUTIONS

! 1. Make indices for W (outer indices)
! 2. Allocate some tmp array
! 3. For each outer index tuple:
! - Make the corresponding W
! - Make the appropriate call to rsp_ovlave and store the result in tmp
! - Make offsets for S (inner indices)
! - Put the elements of tmp in the appropriate place in prop


! write(*,*) 'allocating'

allocate(ncarray(p12(1)%npert + p12(2)%npert))
allocate(ncinner(p12(1)%npert + p12(2)%npert))
! Changed product argument on next line from p12(1) to p12(1)%pdim
allocate(tmp(product(p12(1)%pdim)))
allocate(inner_offsets(product(p12(1)%pdim)))
allocate(outer_indices(product(p12(2)%pdim), p12(2)%npert))

ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
          p12(1)%npert, 1, p12(1), ncarray)

! write(*,*) 'made outer indices'

! Changed 'tot_outer', to p12(2)%npert in call on next line
call make_outerindices(p12(2)%npert, 1, p12(2)%pdim, 0, outer_indices)

do i = 1, size(outer_indices, 1)

tmp = 0.0

! W = mol%zeromat
! UNCOMMENT NEXT OR PREVIOUS LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! W = 0

do j = 1, size(d_ordb, 1)

! NOTE THAT THE RSP_GETW FUNCTIONS SHOULD MAYBE RETURN MATRIX TYPES OR ELSE
! THE ASSIGNMENT TO W COULD BE WRT %ELMS

! SWITCH TO NEXT LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! W = W + rsp_getw(mol, d_ordb(j,:), outer_indices(i,:), F, D, S)

! write(*,*) 'getting W', i, j

W = rsp_getw(mol, d_supsize, d_ordb, p12(1)%npert + p12(2)%npert, outer_indices(i,:), F, D, S)

! write(*,*) 'got W'
            
end do

! MAKE SURE THAT THE SIGN OF THE BELOW EXPRESSION IS OK

! UNCOMMENT NEXT TWO LINES WHEN INTRODUCING THIS CODE TO OPENRSP
! call rsp_ovlave(mol, p12(1)%npert, p12(1)%plab, (/ (j/j, j = 1, p12(1)%npert) /), &
!      ncinner, W, tmp)

! THE SECOND NCARRAY ARGUMENT (THE SECOND LAST ARGUMENT IN TOTAL) IS DUMMY
call get_prop_offsets(size(ncarray), 2, p12(2)%npert, &
     ncarray, ncinner, p12, outer_indices(i,:), &
     ncarray, inner_offsets)

do j = 1, size(inner_offsets)

prop(inner_offsets(j)) = prop(inner_offsets(j)) + tmp(j)
prop_forcache(inner_offsets(j)) = prop_forcache(inner_offsets(j)) + tmp(j)

end do

end do


call propcache_add(cache, 2, p12, propsize, prop_forcache)    


deallocate(d_ordb)
deallocate(ncarray)
deallocate(ncinner)
deallocate(tmp)
deallocate(inner_offsets)
deallocate(outer_indices)



end subroutine




  ! Calculate and add all the Pulay-type (-SW) contributions
  ! Should maybe have a recursive W routine give all W before calling this routine
  ! No, it can be handled with other existing functions

  recursive subroutine rsp_pulay_kn(mol, pert, kn, p12, S, D, F, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert, pdb
    type(perts), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(propcache) :: cache
    integer :: propsize, i
    integer, dimension(2) :: kn
    complex(8), dimension(propsize) :: prop
    
    if (pert%npert > 0) then

!        write(*,*) 'Making recursion 1'
pdb = pert_getone(pert, 1)

! write(*,*) 'p12(1)'
! write(*,*) p12(1)%npert
! write(*,*) p12(1)%pid
! write(*,*) p12(1)%plab
! write(*,*) p12(1)%freq
! 
! write(*,*) 'p12(2)'
! write(*,*) p12(2)%npert
! write(*,*) p12(2)%pid
! write(*,*) p12(2)%plab
! write(*,*) p12(2)%freq

       call rsp_pulay_kn(mol, pert_rf(pert), kn, &
       (/pert_ext(p12(1), pert_getone(pert, 1)), p12(2)/), S, D, F, propsize, &
       cache, prop)

! write(*,*) 'Making recursion 2'

       call rsp_pulay_kn(mol, pert_rf(pert), kn, &
       (/p12(1), pert_ext(p12(2), pert_getone(pert, 1))/), S, D, F, propsize, &
       cache, prop)

! write(*,*) 'Made both recursions'

    else

       if (kn_skip(p12(2)%npert, p12(2)%pid, kn) .EQV. .FALSE.) then


        write(*,*) 'Getting pulay_kn'
        write(*,*) 'p1', p12(1)%pid
        write(*,*) 'p2', p12(2)%pid

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)


if (propcache_already(cache, 2, p12) .EQV. .TRUE.) then

write(*,*) 'Getting values from cache'
write(*,*) ' '

open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

       
else


       ! At lowest level:
        call get_pulay_kn(mol, p12, kn, F, D, S, propsize, cache, prop)

        write(*,*) 'Calculated pulay_kn'
        write(*,*) ' '

end if

       else

       write(*,*) 'pulay kn contribution was k-n skipped:'
         write(*,*) 'p1 ', p12(1)%pid 
write(*,*) 'p2 ', p12(2)%pid 
write(*,*) ' '

       end if 


    end if

  end subroutine







  subroutine get_pulaylag(mol, p12, kn, F, D, S, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert, emptypert
    type(perts), dimension(2) :: p12
    type(perts), dimension(:,:), allocatable :: d_ordb
    type(SDF) :: S, D, F
    type(propcache) :: cache
    type(matrix) :: W
    integer :: i, j, k ,m, incr
    integer :: propsize, d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets
    integer, allocatable, dimension(:) :: outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(propsize) :: prop
    complex(8), dimension(propsize) :: prop_forcache

!     write(*,*) 'inside pulay lag'
!     write(*,*) 'p 1', p12(1)%pid
!     write(*,*) 'p 2', p12(2)%pid


    d_supsize = d_superstr_dry(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/))
   


! write(*,*) 'got supsize', d_supsize

    allocate(d_ordb(d_supsize, 3))
! write(*,*) 'allocated d_ordb'

    ! Get superstructure for W

incr = 0

    call d_superstr(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, incr, d_ordb)

! write(*,*) 'got superstructure'

allocate(ncarray(p12(1)%npert + p12(2)%npert))
allocate(ncinner(p12(1)%npert + p12(2)%npert))
allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))
! Changed product argument on next line from p12(1) to p12(1)%pdim
allocate(tmp(product(p12(1)%pdim)))
allocate(inner_offsets(product(p12(1)%pdim)))
allocate(outer_indices(product(p12(2)%pdim), p12(2)%npert))

! write(*,*) 'allocation was successful'

ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
          p12(1)%npert, 1, p12(1), ncarray)

! Changed 'tot_outer', to p12(2)%npert in call on next line
call make_outerindices(p12(2)%npert, 1, p12(2)%pdim, 0, outer_indices)

do i = 1, size(outer_indices, 1)

tmp = 0.0

! W = mol%zeromat
! UNCOMMENT NEXT OR PREVIOUS LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! W = 0

! do j = 1, size(d_ordb, 1)

! NOTE THAT THE RSP_GETW FUNCTIONS SHOULD MAYBE RETURN MATRIX TYPES OR ELSE
! THE ASSIGNMENT TO W COULD BE WRT %ELMS

! MAYBE LOOK AGAIN AT THE RHS INDICES IN THE outer_ind_b_large ASSIGNMENTS


m = 1

do k = 1, p12(2)%npert

outer_ind_b_large(p12(2)%pid(k)) = outer_indices(i,m)

m = m + 1

end do


! SWITCH TO NEXT TWO LINES WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! W = W + rsp_getw(mol, d_ordb(j,:), p12(1)%npert + p12(2)%npert, &
!                  outer_ind_b_large, F, D, S)
! write(*,*) 'getting w'

W = rsp_getw(mol, d_supsize, d_ordb, p12(1)%npert + p12(2)%npert, &
                 outer_ind_b_large, F, D, S)
            
! end do

! MAKE SURE THAT THE SIGN OF THE BELOW EXPRESSION IS OK

! UNCOMMENT THE NEXT TWO LINES WHEN INTRODUCING THIS CODE TO OPENRSP
! call rsp_ovlave(mol, p12(1)%npert, p12(1)%plab, (/ (j/j, j = 1, p12(1)%npert) /), &
!      ncinner, W, tmp)

! THE SECOND NCARRAY ARGUMENT (THE SECOND LAST ARGUMENT IN TOTAL) IS DUMMY
call get_prop_offsets(size(ncarray), 2, p12(2)%npert, &
     ncarray, ncinner, p12, outer_indices(i,:), &
     ncarray, inner_offsets)

do j = 1, size(inner_offsets)

prop(inner_offsets(j)) = prop(inner_offsets(j)) + tmp(j)
prop_forcache(inner_offsets(j)) = prop_forcache(inner_offsets(j)) + tmp(j)

end do

end do

call propcache_add(cache, 2, p12, propsize, prop_forcache)


deallocate(d_ordb)
deallocate(ncarray)
deallocate(ncinner)
deallocate(outer_ind_b_large)
deallocate(tmp)
deallocate(inner_offsets)
deallocate(outer_indices)





end subroutine




  ! Calculate and add all Pulay-type contributions arising from Lagrange multipliers
  ! Should maybe have a recursive W routine give all W before calling this routine
  ! No, it can be handled with other existing functions

  recursive subroutine rsp_pulay_lag(mol, pert, kn, p12, S, D, F, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(propcache) :: cache
    integer :: propsize, i
    integer, dimension(2) :: kn
    complex(8), dimension(propsize) :: prop
    
    if (pert%npert > 0) then

       call rsp_pulay_lag(mol, pert_rf(pert), kn, &
       (/pert_ext(p12(1), pert_getone(pert, 1)), p12(2)/), S, D, F, propsize, &
       cache, prop)
       call rsp_pulay_lag(mol, pert_rf(pert), kn, &
       (/p12(1), pert_ext(p12(2), pert_getone(pert, 1))/), S, D, F, propsize, &
       cache, prop)

    else

       ! At lowest level:
       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then



       write(*,*) 'Getting pulay lagrange contribution'
       write(*,*) 'p1', p12(1)%pid
       write(*,*) 'p2', p12(2)%pid

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

if (propcache_already(cache, 2, p12) .EQV. .TRUE.) then

write(*,*) 'Getting values from cache'
write(*,*) ' '
       
open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

else


       call get_pulaylag(mol, p12, kn, F, D, S, propsize, cache, prop)

       write(*,*) 'got pulay lagrange contribution'
       write(*,*) ' '

end if

       else

write(*,*) 'pulay lagrange contribution was k-n skipped:'
         write(*,*) 'p1 ', p12(1)%pid 
write(*,*) 'p2 ', p12(2)%pid 
write(*,*) ' '

end if
    end if

  end subroutine





















  subroutine get_idem_lag(mol, p12, kn, F, D, S, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert, emptypert
    type(perts), dimension(2) :: p12
    type(perts), dimension(:,:), allocatable :: d_orda, d_ordb
    type(SDF) :: S, D, F
    type(propcache) :: cache
    type(matrix) :: C, Z
    integer :: i, j, k, m, n, p, incr1, incr2
    integer :: propsize, offset
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, ncprod
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
!     complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(propsize) :: prop
    complex(8), dimension(propsize) :: prop_forcache

!     write(*,*) 'inside idem lag'
!     write(*,*) 'p 1', p12(1)%pid
!     write(*,*) 'p 2', p12(2)%pid

d_supsize = 0

    d_supsize(1) = d_superstr_dry(mol, pert_rf(p12(1)), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = d_superstr_dry(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/))

! write(*,*) 'got supsize', d_supsize
    allocate(d_orda(d_supsize(1), 3))
    allocate(d_ordb(d_supsize(2), 3))

! write(*,*) 'allocated d_ord'
    ! Get superstructure for matrix A, matrix B

incr1 = 0
incr2 = 0

    call d_superstr(mol, pert_rf(p12(1)), kn, .FALSE., & 
         (/get_emptypert(), get_emptypert(), get_emptypert()/), d_supsize(1), incr1, d_orda)
    call d_superstr(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), d_supsize(2), incr2, d_ordb)

! write(*,*) 'got superstructures'


allocate(ncarray(p12(1)%npert + p12(2)%npert))
allocate(ncinner(p12(1)%npert + p12(2)%npert))
allocate(ncprod(p12(1)%npert + p12(2)%npert))
allocate(outer_ind_a_large(p12(1)%npert + p12(2)%npert))
allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))
! allocate(tmp(product(p12(1))))
! allocate(inner_offsets(product(p12(1)%pdim)))
allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%npert))
allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%npert))

! write(*,*) 'allocation was successful'

ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
          p12(1)%npert, 1, p12(1), ncarray)


do i = 1, size(ncarray)

ncprod(i) = product(ncarray(i:size(ncarray)))/ncarray(i)

end do



call make_outerindices(p12(1)%npert, 1, p12(1)%pdim, 0, outer_indices_a)
call make_outerindices(p12(2)%npert, 1, p12(2)%pdim, 0, outer_indices_b)






offset = 0

do i = 1, size(outer_indices_a, 1)

! FOR BOTH C AND Z: SHOULD THEY BE ZEROED AT EACH i, j?
! MAKE SURE THAT THE FUNCTIONS RETURN THE CORRECT INFORMATION FOR C AND Z
! FOR EXAMPLE, SHOULD THE INFORMATION FROM rsp_getc BE PUT IN C%ELMS OR JUST C?
! MAYBE MAKE THE FUNCTIONS RETURN A MATRIX TYPE

! Get matrix A
! STARTING THE LINE WITH C MIGHT BE INTERPRETED AS A COMMENT
! TAKE ANOTHER LOOK AT THIS IF THE RESULTS ARE INCORRECT

! MAYBE LOOK AGAIN AT THE RHS INDICES IN THE outer_ind_a_large AND
! outer_ind_b_large ASSIGNMENTS

! write(*,*) 'oind a at i', i, outer_indices_a(i,:)

p = 1

do n = 1, p12(1)%npert

outer_ind_a_large(p12(1)%pid(n)) = outer_indices_a(i,p)

p = p + 1

end do





 C = rsp_getc(mol, pert_getone(p12(1), 1), kn, d_supsize(1), d_orda, &
              p12(1)%npert + p12(2)%npert, outer_ind_a_large, F, D, S)


do j = 1, size(outer_indices_b, 1)

m = 1

do k = 1, p12(2)%npert

outer_ind_b_large(p12(2)%pid(k)) = outer_indices_b(j,m)

m = m + 1

end do

! Get matrix B

! write(*,*) 'oind', outer_indices_b(j,:)
! write(*,*) 'getting z', j, size(outer_indices_b, 1)

Z = rsp_getz(mol, d_supsize(2), d_ordb, kn, p12(1)%npert + p12(2)%npert, &
             outer_ind_b_large, F, D, S)

! write(*,*) 'got z'

! Calculate offset

do k = 1, p12(1)%npert

offset = offset + ncprod(k)*(outer_indices_a(i, k) - 1)

end do

do k = 1, p12(2)%npert

offset = offset + ncprod(k)*(outer_indices_b(j, k) - 1)

end do

offset = offset + 1


! Add the trace to property at offset
! MAKE SURE THE SIGN BELOW IS CORRECT
! FIND OUT IF THE ARGUMENTS SHOULD BE SEPARATED BY COMMAS
! THE REASON COULD MAYBE BE THAT TIME IS SAVED OVER DOING THE ENTIRE MATRIX MULTIPLICATION

! UNCOMMENT NEXT LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! prop(offset) = -tr(C, Z)
! USE PROPER VALUES ON NEXT LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! prop_forcache(offset) = prop_forcache(offset) -tr(C, Z)



offset = 0

end do

end do

call propcache_add(cache, 2, p12, propsize, prop_forcache) 


deallocate(d_orda)
deallocate(d_ordb)
deallocate(ncarray)
deallocate(ncinner)
deallocate(ncprod)
deallocate(outer_ind_a_large)
deallocate(outer_ind_b_large)
deallocate(outer_indices_a)
deallocate(outer_indices_b)





end subroutine




  ! Calculate and add all contributions arising from idempotency Lagrange multipliers
  ! Should maybe have a recursive W routine give all W before calling this routine
  ! No, it can be handled with other existing functions

  recursive subroutine rsp_idem_lag(mol, pert, kn, p12, S, D, F, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(propcache) :: cache
    integer :: propsize, i
    integer, dimension(2) :: kn
    complex(8), dimension(propsize) :: prop
    
    if (pert%npert > 0) then

       call rsp_idem_lag(mol, pert_rf(pert), kn, &
       (/pert_ext(p12(1), pert_getone(pert, 1)), p12(2)/), S, D, F, propsize, &
       cache, prop)
       call rsp_idem_lag(mol, pert_rf(pert), kn, &
       (/p12(1), pert_ext(p12(2), pert_getone(pert, 1))/), S, D, F, propsize, &
       cache, prop)

    else

       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then

       write(*,*) 'Getting idempotency lagrange contribution'
       write(*,*) 'p1', p12(1)%pid
       write(*,*) 'p2', p12(2)%pid

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

if (propcache_already(cache, 2, p12) .EQV. .TRUE.) then

write(*,*) 'Getting values from cache'
write(*,*) ' '

open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)
      
else



       ! At lowest level:
       call get_idem_lag(mol, p12, kn, F, D, S, propsize, cache, prop)

        write(*,*) 'Calculated idempotency lagrange contribution'
        write(*,*) ' '

end if


       else

write(*,*) 'idempotency lagrange contribution was k-n skipped:'
         write(*,*) 'p1 ', p12(1)%pid 
write(*,*) 'p2 ', p12(2)%pid 
write(*,*) ' '

end if



    end if

  end subroutine


































  subroutine get_scfe_lag(mol, p12, kn, F, D, S, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert, emptypert
    type(perts), dimension(2) :: p12
    type(perts), dimension(:,:), allocatable :: d_orda, d_ordb
    type(SDF) :: S, D, F
    type(propcache) :: cache
    type(matrix) :: L, Y
    integer :: i, j, k, m, n, p, incr1, incr2
    integer :: propsize, offset
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, ncprod
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
!     complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(propsize) :: prop
    complex(8), dimension(propsize) :: prop_forcache


d_supsize = 0

    d_supsize(1) = d_superstr_dry(mol, pert_rf(p12(1)), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = d_superstr_dry(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/))


    allocate(d_orda(d_supsize(1), 3))
    allocate(d_ordb(d_supsize(2), 3))


    ! Get superstructure for matrix A, matrix B

incr1 = 0
incr2 = 0

    call d_superstr(mol, pert_rf(p12(1)), kn, .FALSE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/),& 
                    d_supsize(1), incr1, d_orda)
    call d_superstr(mol, p12(2), kn, .TRUE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                    d_supsize(2), incr2, d_ordb)

allocate(ncarray(p12(1)%npert + p12(2)%npert))
allocate(ncinner(p12(1)%npert + p12(2)%npert))
allocate(ncprod(p12(1)%npert + p12(2)%npert))
! allocate(tmp(product(p12(1))))
! allocate(inner_offsets(product(p12(1)%pdim)))
allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%npert))
allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%npert))
allocate(outer_ind_a_large(p12(1)%npert + p12(2)%npert))
allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))

ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
          p12(1)%npert, 1, p12(1), ncarray)


do i = 1, size(ncarray)

ncprod(i) = product(ncarray(i:size(ncarray)))/ncarray(i)

end do



call make_outerindices(p12(1)%npert, 1, p12(1)%pdim, 0, outer_indices_a)
call make_outerindices(p12(2)%npert, 1, p12(2)%pdim, 0, outer_indices_b)






offset = 0

do i = 1, size(outer_indices_a, 1)

! FOR BOTH L AND Y: SHOULD THEY BE ZEROED AT EACH i, j?
! MAKE SURE THAT THE FUNCTIONS RETURN THE CORRECT INFORMATION FOR C AND Z
! FOR EXAMPLE, SHOULD THE INFORMATION FROM rsp_getc BE PUT IN C%ELMS OR JUST C?
! MAYBE MAKE THE FUNCTIONS RETURN A MATRIX TYPE

! Get matrix A
! STARTING THE LINE WITH C MIGHT BE INTERPRETED AS A COMMENT
! TAKE ANOTHER LOOK AT THIS IF THE RESULTS ARE INCORRECT

! MAYBE LOOK AGAIN AT THE RHS INDICES IN THE outer_ind_a_large AND
! outer_ind_b_large ASSIGNMENTS


p = 1

do n = 1, p12(1)%npert

outer_ind_a_large(p12(1)%pid(n)) = outer_indices_a(i,p)

p = p + 1

end do




 L = rsp_getl(mol, pert_getone(p12(1), 1), d_supsize(1), d_orda, &
              p12(1)%npert + p12(2)%npert, outer_ind_a_large, D, S)


do j = 1, size(outer_indices_b, 1)

m = 1

do k = 1, p12(2)%npert

outer_ind_b_large(p12(2)%pid(k)) = outer_indices_b(j,m)

m = m + 1

end do

 

! Get matrix B

Y = rsp_gety(mol, d_supsize(2), d_ordb, p12(1)%npert + p12(2)%npert, &
             outer_ind_b_large, F, D, S)

! Calculate offset

do k = 1, p12(1)%npert

offset = offset + ncprod(k)*(outer_indices_a(i, k) - 1)

end do

do k = 1, p12(2)%npert

offset = offset + ncprod(k)*(outer_indices_b(j, k) - 1)

end do

offset = offset + 1


! Add the trace to property at offset
! MAKE SURE THE SIGN BELOW IS CORRECT
! FIND OUT IF THE ARGUMENTS SHOULD BE SEPARATED BY COMMAS
! THE REASON COULD MAYBE BE THAT TIME IS SAVED OVER DOING THE ENTIRE MATRIX MULTIPLICATION

! UNCOMMENT NEXT LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! prop(offset) = prop(offset) - tr(L, Y)
! USE PROPER VALUES ON NEXT LINE WHEN MATRIX FUNCTIONALITY IS INTRODUCED
! prop_forcache(offset) = prop_forcache(offset) - tr(L, Y)


offset = 0

end do

end do

call propcache_add(cache, 2, p12, propsize, prop_forcache)


deallocate(d_orda)
deallocate(d_ordb)
deallocate(ncarray)
deallocate(ncinner)
deallocate(ncprod)
deallocate(outer_indices_a)
deallocate(outer_indices_b)
deallocate(outer_ind_a_large)
deallocate(outer_ind_b_large)




end subroutine




  ! Calculate and add all contributions arising from idempotency Lagrange multipliers
  ! Should maybe have a recursive W routine give all W before calling this routine
  ! No, it can be handled with other existing functions

  recursive subroutine rsp_scfe_lag(mol, pert, kn, p12, S, D, F, propsize, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(propcache) :: cache
    integer :: propsize, i
    integer, dimension(2) :: kn
    complex(8), dimension(propsize) :: prop
    
    if (pert%npert > 0) then

       call rsp_scfe_lag(mol, pert_rf(pert), kn, &
       (/pert_ext(p12(1), pert_getone(pert, 1)), p12(2)/), S, D, F, propsize, &
       cache, prop)
       call rsp_scfe_lag(mol, pert_rf(pert), kn, &
       (/p12(1), pert_ext(p12(2), pert_getone(pert, 1))/), S, D, F, propsize, &
       cache, prop)

    else

       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then

       write(*,*) 'Getting scfe lagrange contribution'
       write(*,*) 'p1', p12(1)%pid
       write(*,*) 'p2', p12(2)%pid

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

if (propcache_already(cache, 2, p12) .EQV. .TRUE.) then


open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'T'

close(257)

write(*,*) 'Getting values from cache'
write(*,*) ' '
       
else


       ! At lowest level:
       call get_scfe_lag(mol, p12, kn, F, D, S, propsize, cache, prop)

        write(*,*) 'Calculated scfe lagrange contribution'
        write(*,*) ' '
end if


       else

write(*,*) 'scfe lagrange contribution was k-n skipped:'
         write(*,*) 'p1 ', p12(1)%pid 
write(*,*) 'p2 ', p12(2)%pid 
write(*,*) ' '



end if


      

    end if

  end subroutine





















  subroutine make_pert_subset(pert, psub)

    implicit none

    type(perts) :: pert
    type(perts), dimension(pert%npert) :: psub
    integer :: i, j, k, m

    do i = 1, pert%npert

       psub(i)%npert = pert%npert - 1

       allocate(psub(i)%pdim(pert%npert - 1))
       allocate(psub(i)%plab(pert%npert - 1))
       allocate(psub(i)%pid(pert%npert - 1))
       allocate(psub(i)%freq(pert%npert - 1))

       k = i

       m = 1

       do j = 1, (pert%npert)

          if (.NOT.(j == k)) then

!              write(*,*) 'i, j, k, m', i, j, k, m

             psub(i)%pdim(m) = pert%pdim(j)
             psub(i)%plab(m) = pert%plab(j)
             psub(i)%pid(m) = pert%pid(j)
             psub(i)%freq(m) = pert%freq(j)

             m = m + 1

          else

!              k = k + 1

          end if

       end do
    end do
  
  end subroutine

subroutine getfds_dummy(pert, F, D, S)

integer :: data_f, data_d, data_s
type(perts) :: pert
type(SDF) :: F, D, S

! Get the appropriate Fock/density/overlap matrix
! Rearrange pert and data in predetermined order

! write(*,*) 'inside getfds_dummy'

call sdf_standardorder(pert, data_f, data_d, data_s)

! write(*,*) 'got sdf in standard order'

data_f = 1
data_d = 1
data_s = 1

! Add to cache
call sdf_add(F, pert, data_f)
!  write(*,*) 'added F to cache', F%perturb%npert, F%next%perturb%npert
call sdf_add(D, pert, data_d)
!  write(*,*) 'added D to cache', D%perturb%npert, D%next%perturb%npert
call sdf_add(S, pert, data_s)
!  write(*,*) 'added S to cache', S%perturb%npert, S%next%perturb%npert


end subroutine






    ! Assume there is a pert_dens routine or equivalent one-line call for all of F, D, S
    ! Could write the logic with pseudo for such a routine/routines for now
    ! F, D, and S are each linked lists of the FDS type

  recursive subroutine rsp_fds(mol, pert, kn, F, D, S)

    implicit none

    type(rsp_cfg) :: mol
    type(perts) :: pert
    type(perts), dimension(pert%npert) :: psub
    integer, dimension(2) :: kn
    integer :: i
    type(SDF) :: F, D, S

    ! Unless at lowest level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then get F, D, S 


!     write(*,*) 'I was called with ', pert%npert, ' ', pert%plab

    if (pert%npert > 1) then



       call make_pert_subset(pert, psub)

       do i = 1, size(psub)

!           write(*,*) 'Recursing with', psub(i)%plab
          call rsp_fds(mol, psub(i), kn, F, D, S)
       end do       

    else

!        write(*,*) 'Not recursing'

    end if

    if (sdf_already(D, pert) .eqv. .FALSE.) then
         
       ! Concatenate relevant S (Make array of pointers)

!        write(*,*) 'Calling ovlint with perts ', pert%plab

       !call rsp_ovlint(mol, pert%npert, pert%plab, first, pert%pdim, S(relevant))

       if (kn_skip(pert%npert, pert%pid, kn) .eqv. .FALSE.) then

          ! Concatenate relevant F, D (Make array of pointers)

                 write(*,*) 'Calling ovlint/fock/density with perts ', pert%plab, &
                            'and pertid ', pert%pid
                 
                 call getfds_dummy(pert, F, D, S)

       else

                 write(*,*) 'Would have called ovlint/fock/density with perts ', &
                            pert%plab, &
                            'and pertid ', pert%pid, ' but it was k-n forbidden'

       end if

    else

       write(*,*) 'FDS for perts ', pert%plab, &
                            'and pertid ', pert%pid, ' was found in cache'

    end if

  end subroutine

function rm_dim(np, dims, skip)

  integer :: skip, i, j, np
  integer, dimension(np) :: dims
  integer, dimension(np - 1) :: rm_dim

  j = 1

  do i = 1, np

     if (i == skip) then

     else

       rm_dim(j) = dims(i)
  
       j = j + 1   

     end if

  end do




end function


recursive function kn_prod(np, lvl, dims, k_or_n, sofar) result (kn_p)

  integer :: np, lvl, k_or_n, sofar, kn_p, i
  integer, dimension(np) :: dims

! write(*,*) 'Called at lvl', lvl

kn_p = 0

if (lvl < k_or_n) then

! write(*,*) 'did recursion'

do i = 1, np

kn_p = kn_p + kn_prod(np - 1, lvl + 1, rm_dim(np, dims, i), k_or_n, dims(i)*sofar)

end do

else


kn_p = sofar

end if

!  write(*,*) 'finally ', kn_p

end function


  ! Calculate the best pair (k,n) for calculation of property
  ! Takes each possible combination (k,n) and estimates the complexity
  ! of the property calculation
  ! Returns the pair (k,n) for which the estimated complexity is lowest

  function get_bestkn(pert)
  
  type(perts) :: pert
  integer, dimension(2) :: get_bestkn
  integer :: min_n, n, csize, minsize



  ! Do the case n = 1 first to establish one value

    minsize = kn_prod(pert%npert - 1, 1, pert%pdim(2:pert%npert), &
              pert%npert - 1 - 1, pert%pdim(1)) 

!     write(*,*) 'The include pert A part of minsize is', minsize

    minsize = minsize + kn_prod(pert%npert - 1, &
              0, pert%pdim(2:pert%npert), 1, 1)

!    write(*,*) 'The size for n = 1 is', minsize

    min_n = 1

  ! Get the products for the pert dimensions as dictated by k and n
  ! Assume that integer division takes care of this

    do n = (pert%npert/2), pert%npert - 1

       csize = kn_prod(pert%npert - 1, 1, pert%pdim(2:pert%npert), &
                 pert%npert - n - 1, pert%pdim(1))

!     write(*,*) 'The include pert A part of minsize is', csize

       csize = csize + kn_prod(pert%npert - 1, &
                 0, pert%pdim(2:pert%npert), n, 1)

! write(*,*) 'The size for n = ', n, ' is', csize


       ! If the products are smaller than the previous min, update 
       ! index identifer and minsize

       if (csize < minsize) then

          min_n = n
          minsize = csize

       end if


    end do

    get_bestkn(2) = min_n
    get_bestkn(1) = pert%npert - min_n - 1


  end function


subroutine propcache_allocate(inst)

type(propcache), pointer :: inst

allocate(inst)

inst%next => inst
inst%last = .TRUE.
inst%tlen = 0

allocate(inst%t_orders(1))
allocate(inst%t_orders(1)%pdim(1))
allocate(inst%t_orders(1)%plab(1))
allocate(inst%t_orders(1)%pid(1))
allocate(inst%t_orders(1)%freq(1))

inst%t_orders(1)%pdim = (/0/)
inst%t_orders(1)%plab = (/'NUTN'/)
inst%t_orders(1)%pid = (/0/)
inst%t_orders(1)%freq = (/0.0/)


allocate(inst%data(1))

end subroutine



  subroutine get_prop(mol, pert, kn, prop, F, D, S)

    type(SDF) :: F, D, S
    type(rsp_cfg) :: mol
    type(perts) :: pert, emptypert
    type(perts), dimension(2) :: emptyperts
    integer, dimension(2) :: kn
    complex(8), dimension(product(pert%pdim)) :: prop
    type(propcache), pointer :: energy_cache, pulay_kn_cache, &
                       pulay_lag_cache, idem_cache, scfe_cache

emptypert%npert = 0
allocate(emptypert%pdim(0))    
allocate(emptypert%plab(0))
allocate(emptypert%pid(0))
allocate(emptypert%freq(0))

emptyperts = (/emptypert, emptypert/)

  ! Get all necessary F, D, S derivatives as dictated by
  ! number of perturbations and kn
  ! Maybe also extend to get W for kn part of rsp_pulay

    call rsp_fds(mol, pert, kn, F, D, S)

  ! Allocate the energy contribution property cache

  

call propcache_allocate(energy_cache)

  ! Calculate and add all the energy contributions
 write(*,*) 'rsp fds call finished'

    call rsp_ener(mol, pert, pert%npert, kn, 1, (/emptypert/), 0, D, &
         product(pert%pdim), energy_cache, prop)

 write(*,*) 'rsp ener call finished'

deallocate(energy_cache)

 
call propcache_allocate(pulay_kn_cache)
 
!   ! Calculate and add all the kn-type Pulay-type ((-SW)_(k,n)_W) contributions
 
    call rsp_pulay_kn(mol, pert, kn, (/emptypert, emptypert/), S, D, F, &
        product(pert%pdim), pulay_kn_cache, prop)


 write(*,*) 'rsp pulay_kn call finished'

deallocate(pulay_kn_cache)

call propcache_allocate(pulay_lag_cache)

!   ! Calculate and add all the Lagrange multiplier Pulay-type contributions

    call rsp_pulay_lag(mol, pert_rf(pert), kn, (/pert_getone(pert,1), emptypert/), &
         S, D, F, product(pert%pdim), pulay_lag_cache, prop)

 write(*,*) 'rsp pulay_lag call finished' 

deallocate(pulay_lag_cache)

call propcache_allocate(idem_cache)

!   ! Calculate and add all the idempotency-type contributions
 
    call rsp_idem_lag(mol, pert_rf(pert), kn, (/pert_getone(pert,1), emptypert/), &
         S, D, F, product(pert%pdim), idem_cache, prop)

 write(*,*) 'rsp idem_lag call finished'

deallocate(idem_cache)

call propcache_allocate(scfe_cache)
 
!   ! Calculate and add all the SCF-type contributions
 
    call rsp_scfe_lag(mol, pert_rf(pert), kn, (/pert_getone(pert,1), emptypert/), &
         S, D, F, product(pert%pdim), scfe_cache, prop)

 write(*,*) 'rsp scfe_lag call finished'

deallocate(scfe_cache)


  end subroutine

end module


program testpert

use proprty

implicit none

type(SDF), pointer :: S, D, F
type(rsp_cfg) :: mol
type(perts) :: pert
integer, dimension(2) :: kn
complex(8), dimension (:), allocatable :: prop


! VARIOUS PERTURBATION TUPLE SETUPS
! COMMENT/UNCOMMENT FOR TESTING

! Molecular Hessian
!
! pert%npert = 2
! allocate(pert%pdim(pert%npert))
! allocate(pert%plab(pert%npert))
! allocate(pert%pid(pert%npert))
! allocate(pert%freq(pert%npert))
! 
! pert%pdim = (/12, 12/)
! pert%plab = (/'GEO ', 'GEO '/)
! pert%pid = (/1, 2/)
! pert%freq = (/0.0, 0.0/)
!
! Cubic force field
! 
! pert%npert = 3
! allocate(pert%pdim(pert%npert))
! allocate(pert%plab(pert%npert))
! allocate(pert%pid(pert%npert))
! allocate(pert%freq(pert%npert))
! 
! pert%pdim = (/12, 12, 12/)
! pert%plab = (/'GEO ', 'GEO ', 'GEO '/)
! pert%pid = (/1, 2, 3/)
! pert%freq = (/0.0, 0.0, 0.0/)
!
! Quartic force field
!
! pert%npert = 4
! allocate(pert%pdim(pert%npert))
! allocate(pert%plab(pert%npert))
! allocate(pert%pid(pert%npert))
! allocate(pert%freq(pert%npert))
! 
! pert%pdim = (/12, 12, 12, 12/)
! pert%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO '/)
! pert%pid = (/1, 2, 3, 4/)
! pert%freq = (/0.0, 0.0, 0.0, 0.0/)
!
! Hessian of first hyperpolarizability
!
pert%npert = 5
allocate(pert%pdim(pert%npert))
allocate(pert%plab(pert%npert))
allocate(pert%pid(pert%npert))
allocate(pert%freq(pert%npert))

pert%pdim = (/12, 12, 3, 3, 3/)
pert%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
pert%pid = (/1, 2, 3, 4, 5/)
pert%freq = (/0.0, 0.0, 0.01, 0.02, 0.01/)
! 
! Cubic force field of first hyperpolarizability 
!
! pert%npert = 6
! allocate(pert%pdim(pert%npert))
! allocate(pert%plab(pert%npert))
! allocate(pert%pid(pert%npert))
! allocate(pert%freq(pert%npert))
! 
! pert%pdim = (/6, 6, 6, 3, 3, 3/)
! pert%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
! pert%pid = (/1, 2, 3, 4, 5, 6/)
! pert%freq = (/0.0, 0.0, 0.0, 0.01, 0.02, 0.01/)

open(unit=257, file='totterms', status='replace', action='write') 

write(257,*) 'START'

close(257)

open(unit=257, file='cachehit', status='replace', action='write') 

write(257,*) 'START'

close(257)

kn = get_bestkn(pert)
write(*,*) 'Best k, n was found to be ', kn(1), ' and ', kn(2)

allocate(S)
S%next => S
S%last = .TRUE.
S%perturb%npert = 0
allocate(S%perturb%pdim(1))
allocate(S%perturb%plab(1))
allocate(S%perturb%pid(1))
allocate(S%perturb%freq(1))
S%perturb%pdim = (/0/)
S%perturb%plab = (/'NUTN'/)
S%perturb%pid = (/0/)
S%perturb%freq = (/0.0/)
S%data = 1
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!allocate(S%data(1))

allocate(D)
D%next => D
D%last = .TRUE.
D%perturb%npert = 0
allocate(D%perturb%pdim(1))
allocate(D%perturb%plab(1))
allocate(D%perturb%pid(1))
allocate(D%perturb%freq(1))
D%perturb%pdim = (/0/)
D%perturb%plab = (/'NUTN'/)
D%perturb%pid = (/0/)
D%perturb%freq = (/0.0/)
D%data = 2
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!allocate(D%data(1))

allocate(F)
F%next => F
F%last = .TRUE.
F%perturb%npert = 0
allocate(F%perturb%pdim(1))
allocate(F%perturb%plab(1))
allocate(F%perturb%pid(1))
allocate(F%perturb%freq(1))
F%perturb%pdim = (/0/)
F%perturb%plab = (/'NUTN'/)
F%perturb%pid = (/0/)
F%perturb%freq = (/0.0/)
F%data = 3
! UNCOMMENT NEXT LINE WHEN MORE FUNCTIONALITY IS ADDED
!allocate(F%data(1))



allocate(prop(product(pert%pdim)))

call get_prop(mol, pert, kn, prop, F, D, S)

open(unit=257, file='totterms', status='old', action='write', position='append') 

write(257,*) 'END'

close(257)

open(unit=257, file='cachehit', status='old', action='write', position='append') 

write(257,*) 'END'

close(257)

end program