! non-copyrighted template

!> @file
!> Contains module matrix_genop

module matrix_genop

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
  public matrix_genop_debug


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
  logical matrix_genop_debug = .false.

  private

contains


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
if (matrix_genop_debug) print *, 'alloc'
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
if (matrix_genop_debug) print *, 'free'
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
if (matrix_genop_debug) &
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
if (matrix_genop_debug) &
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
