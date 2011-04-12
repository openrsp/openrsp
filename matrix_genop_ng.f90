!> @file
!> Contains module matrix_defop

module matrix_genop_ng

  implicit none

  public matrix
  public mat_isdef
  public mat_init
  public mat_free
  public mat_same
  public mat_axpy
  public mat_gemm
  public mat_dot
  public mat_trace
  public mat_print
  public mat_dup
  public mat_isdup
  public matrix_genop_debug
  !-- used by matrix_defop --
  public matf_alias
  public matf_zero
  public matf_temp
  public matf_temp_tx
  public matf_temp_ty
  !-- not used by matrix_defop --
  public matf_clsh

  !----- constants for type(matrix)%flags -----
  !> matrix is closed-shell: dot(A,B) = 2*sum(A%elms*B%elms).
  !> True for density and Fock matrices,
  !> false for overlap and 1e-hamil, etc.
  integer, parameter :: matf_clsh  = 1
  !> unused bits in %flags should contain a "magic" number.
  !> This is used to detect memory corruption and accidental
  !> (ab)use of undefined matrices, to thereafter fail gracefully,
  !> while avoiding use of type field defaults (F95 feature).
  integer, parameter :: matf_magic_mask = not(2**20-1)
  integer, parameter :: matf_magic_num = 1765 * 2**20

  !----- flags used by matrix_defop ------
  !> matrix is a temporary, to be deleted after use
  integer, parameter :: matf_temp  = 11
  !> matrix is an alias of (a part of) another matrix,
  !> and should not be overwritten but re-allocated,
  !> and not free'ed, just nullified
  integer, parameter :: matf_alias = 10
  !> matrix is an unallocated zero matrix, which should
  !> never be passed to numerical routines mat_axpy,
  !> mat_gemm or mat_dot
  integer, parameter :: matf_zero  = 7
  !> A transposed in  : C += fac * A^ta * B^tb
  integer, parameter :: matf_temp_tx = 13
  !> B transposed in  : C += fac * A^ta * B^tb
  integer, parameter :: matf_temp_ty    = 14

  !> matrix structure
  type matrix
     !> number of rows of each block
     integer          :: nrow
     !> number of columns of each block
     integer          :: ncol
     !> flags telling which variant of matrix this is
     integer          :: flags
     !> main block of elements. All non-zero matrices have this
     real(8), pointer :: elms(:,:) !elms(nrow,ncol)
     !> during defined operator evaluation (matrix_defop), temporary
     !> matrices have additional contribution: fac * A^ta * B^tb,
     !> where fac, X and Z are here, while ta and tb are kept inside %flags
     complex(8)            :: temp_fac
     type(matrix), pointer :: temp_X, temp_Y
  end type

  !> for switching debugging on and off
  logical :: matrix_genop_debug = .false.

  private

contains


  !> whether A is defined, according to A%flags
  function mat_isdef(A)
    type(matrix), intent(in) :: A
    logical                  :: mat_isdef
    mat_isdef = (iand(A%flags,matf_magic_mask) &
                           == matf_magic_num)
  end function


  !> free A, and reset it to undefined. If nodealloc=true,
  !> don't deallocate A%elms, just nullify
  subroutine mat_free(A, nodealloc)
    type(matrix),   intent(inout) :: A
    logical, optional, intent(in) :: nodealloc
    logical :: nodeall
    ! err if attempting to free undefined matrix
    if (.not.mat_isdef(A)) &
       call quit('error: mat_free(A), a undefined')
    ! nodealloc defaults to false
    nodeall = .false.
    if (present(nodealloc)) nodeall = nodealloc
    ! deallocate
    if (.not.nodeall .and. .not.associated(A%elms)) then
       call quit('error: matrix mat_free(A), expected A%elms to be allocated')
    else if (.not.nodeall) then
       !print *, 'deallocate', loc(A%elms)
       deallocate(A%elms)
    else
       nullify(A%elms)
    end if
    A%nrow  = -huge(1) !will err if used
    A%ncol  =  huge(1) !will err if used
    A%flags = 0 !undefined, since doesn't have matf_magic_num
    A%temp_fac = 1
    nullify(A%temp_X)
    nullify(A%temp_Y)
  end subroutine


  !> initialize C in the shape of A, A^T, A*B, A^T*B, A*B^T or A^T*B^T.
  !> Optionally zero the elements or skip allocation
  subroutine mat_init(C, A, B, ta, tb, zero, noalloc)
    type(matrix),        intent(inout) :: C
    type(matrix),           intent(in) :: A
    type(matrix), optional, intent(in) :: B
    logical,      optional, intent(in) :: ta, tb, zero, noalloc
    logical :: taa, tbb, zer, alc
    integer :: nrow
    if (present(tb) .and. .not.present(B)) &
       call quit('error: mat_init(C,A,B,ta,tb), tb present without B')
    ! process optionals
    taa = .false.; if (present(ta)) taa = ta
    tbb = .false.; if (present(tb)) tbb = tb
    zer = .false.; if (present(zero)) zer = zero
    alc = .true.;  if (present(noalloc)) alc = (.not.noalloc)
    ! set nrow and ncol, keeping in mind that A may be C
    nrow = merge(A%ncol, A%nrow, taa) !don't overwrite C%nrow yet
    if (.not.present(B)) C%ncol = merge(A%nrow, A%ncol, taa)
    if (     present(B)) C%ncol = merge(B%nrow, B%ncol, tbb)
    C%nrow = nrow
    ! inherit flags from A and B, then insert magic
    C%flags = A%flags
    if (present(B)) C%flags = ior(C%flags, B%flags)
    C%flags = iand(C%flags, not(matf_magic_mask))
    C%flags =  ior(C%flags, matf_magic_num)
    ! defop-related flags are not inherited
    C%flags = ibclr(C%flags, matf_zero)
    C%flags = ibclr(C%flags, matf_temp)
    C%flags = ibclr(C%flags, matf_alias)
    C%flags = ibclr(C%flags, matf_temp_tx)
    C%flags = ibclr(C%flags, matf_temp_ty)
    C%temp_fac = 1
    nullify(C%temp_X)
    nullify(C%temp_Y)
    if (alc) allocate(C%elms(C%nrow,C%ncol))
    !if (alc) print *, '--allocate', loc(C%elms)
    if (alc .and. zer) C%elms = 0
    if (.not.alc) nullify(C%elms)
  end subroutine


  !> whether matrices A and B are of same shape
  function mat_same(A, B, ta, tb, mul)
    type(matrix),      intent(in) :: A, B
    logical, optional, intent(in) :: ta, tb, mul
    logical :: taa, tbb, mmul, mat_same
    taa  = .false.; if (present(ta))  taa  = ta
    tbb  = .false.; if (present(tb))  tbb  = tb
    mmul = .false.; if (present(mul)) mmul = mul
    if (mmul) then
       mat_same = (merge(A%nrow, A%ncol, taa) &
                == merge(B%ncol, B%nrow, tbb))
    else
       mat_same = (merge(A%ncol, A%nrow, taa) &
                == merge(B%ncol, B%nrow, tbb) &
             .and. merge(A%nrow, A%ncol, taa) &
                == merge(B%nrow, B%ncol, tbb))
    end if
  end function


  !> Y+=a*X, Y=a*X, Y+=a*X^T or Y=a*X^T
  subroutine mat_axpy(a, X, tx, py, Y)
    complex(8),   intent(in)    :: a
    type(matrix), intent(in)    :: X
    logical,      intent(in)    :: tx, py
    type(matrix), intent(inout) :: Y
    logical :: XisY
    XisY = associated(X%elms,Y%elms)
    if (tx .and. py .and. XisY) then
       Y%elms = Y%elms + a * transpose(Y%elms)
    else if (tx .and. py) then
       Y%elms = Y%elms + a * transpose(X%elms)
    else if (tx .and. XisY) then
       Y%elms = a * transpose(Y%elms)
    else if (tx) then
       Y%elms = a * transpose(X%elms)
    else if (py .and. XisY) then
       Y%elms = Y%elms + a * Y%elms
    else if (py) then
       Y%elms = Y%elms + a * X%elms
    else if (XisY) then
       Y%elms = a * Y%elms
    else
       Y%elms = a * X%elms
    end if
  end subroutine



  !> matrix multiply: C = rc*C + rab * A^ta * B^tb
  !> assumed initialied and with corresponding shapes
  subroutine mat_gemm(fac, A, ta, B, tb, pc, C)
    complex(8),      intent(in) :: fac
    type(matrix),    intent(in) :: A, B
    logical,         intent(in) :: ta, tb, pc
    type(matrix), intent(inout) :: C
    if (ta .and. tb .and. pc) then
       C%elms = C%elms + fac * matmul(transpose(A%elms), &
                                      transpose(B%elms))
    else if (ta .and. tb) then
       C%elms = fac * matmul(transpose(A%elms), &
                             transpose(B%elms))
    else if (ta .and. pc) then
       C%elms = C%elms + fac * matmul(transpose(A%elms), B%elms)
    else if (ta) then
       C%elms = fac * matmul(transpose(A%elms), B%elms)
    else if (tb .and. pc) then
       C%elms = C%elms + fac * matmul(A%elms, transpose(B%elms))
    else if (tb) then
       C%elms = fac * matmul(A%elms, transpose(B%elms))
    else if (pc) then
       C%elms = C%elms + fac * matmul(A%elms, B%elms)
    else
       C%elms = fac * matmul(A%elms, B%elms)
    end if
  end subroutine



  !> matrix dot product: sum_ij Aij^* Bij (ta=F)
  !>  or  product trace: sum_ij Aji Bij   (ta=F)
  !> A and B are assumed initialized and with equal/transpose shapes
  function mat_dot(A, B, t)
    type(matrix), intent(in) :: A, B
    logical,      intent(in) :: t !transpose
    complex(8)               :: mat_dot
    if (   t  ) mat_dot = sum(transpose(A%elms) * B%elms)
    if (.not.t) mat_dot = sum(A%elms * B%elms)
    ! if closed-shell, then additional factor 2
    if (btest(A%flags, matf_clsh) .or. &
        btest(B%flags, matf_clsh)) mat_dot = 2*mat_dot
  end function


  !> matrix trace: sum_i Aii, A assumed initialized
  function mat_trace(A)
    type(matrix), intent(in) :: A
    complex(8)               :: mat_trace
    integer :: i
    if (A%nrow /= A%ncol) &
       call quit('error: mat_trace(A) with A non-square')
    mat_trace = sum((/(A%elms(i,i),i=1,A%nrow)/))
    ! if closed-shell, then additional factor 2
    if (btest(A%flags, matf_clsh)) mat_trace = 2*mat_trace
  end function



  !> print matrix A to optional unit, optionally starting with label, &
  !> optionally making each column 'width' wide, with optional braces
  !> 'decor' (left brace, column separator, right brace, row separator)
  subroutine mat_print(A, label, unit, width, decor)
    type(matrix),           intent(in) :: A
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: decor
    integer      :: uni, colw, dec, i, j, siz
    character(8) :: fmt
    character(4) :: brac
    ! process optional argument unit, which defaults to stdout
    uni = 6
    if (present(unit)) uni = unit
    ! set d to the largest number of digits to be printed
    ! before the decimal point (including - signs)
    dec = ceiling(max(log10(max(1d0, maxval(A%elms)))+1, &
                      log10(max(1d0,-minval(A%elms)))+2))
    ! process optional width
    colw = 9
    if (present(width)) colw = max(width,dec+2) !max, to avoid stars *****
    ! set dec to the number of decimals to be printed
    dec = colw - dec - 1
    ! argument label is optional. If present, print that
    if (present(label)) write (uni,'(a)') label
    ! process optional argument braces, defaulting to none
    brac = '    '
    if (present(decor)) brac = decor
    ! create the format string to be used for each element
    fmt = '(fww.dd)'
    write (fmt(3:7),'(i2,a1,i2)') colw, '.', dec
    ! call the printing routine
    call subr(A%nrow, A%ncol, A%elms, colw, fmt, brac, uni)
    ! final empty line
    write (uni,'()')

  contains

    subroutine subr(nrow, ncol, elms, colw, fmt, brac, unit)
      integer,        intent(in) :: nrow, ncol, colw, unit
      real(8),        intent(in) :: elms(nrow,ncol)
      character(8),   intent(in) :: fmt
      character(4),   intent(in) :: brac
      character(ncol*(colw+1)+3) :: line
      integer :: i, j, l
      do i = 1, nrow
         line(1:1) = merge(brac(1:1),' ',i==1)
         line(2:2) = brac(1:1)
         l = 3
         do j = 1, ncol
            write (line(l:l+colw-1), fmt) elms(i,j)
            l = l + colw 
            if (j/=ncol) line(l:l) = brac(2:2)
            if (j/=ncol) l = l+1
         end do
         line(l:l) = brac(3:3)
         line(l+1:l+1) = merge(brac(3:3), brac(4:4), i==nrow)
         l = l + 2
         if (brac=='    ') write (unit,'(a)') line(3:len(line)-2)
         if (brac/='    ') write (unit,'(a)') line
      end do
    end subroutine

  end subroutine


  !> copy all fields of A into B, creating a duplicate
  subroutine mat_dup(A, B)
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
