!> @file
!> Contains module matrix_defop

module matrix_genop

  implicit none

  public matrix
  public mat_isdef
  public mat_init
  public mat_free
  public mat_same
  public mat_axtpy
  public mat_gemm
  public mat_dot
  public mat_print
  public mat_dup
  public matrix_genop_debug

  !----- constants for type(matrix)%flags -----
  !> matrix is closed-shell: dot(A,B) = 2*sum(A%elms*B%elms).
  !> True for density and Fock matrices,
  !> false for overlap and 1e-hamil, etc.
  integer, parameter :: matf_clsh  = 1

  !> matrix is an unallocated zero matrix
  integer, parameter :: matf_zero  = 7

  !> matrix is an alias of (a part of) another matrix,
  !> and should not be overwritten but re-allocated,
  !> and not free'ed, just nullified
  integer, parameter :: matf_alias = 10

  !----- flags used by matrix_defop only ------
  !> matrix is a temporary, to be deleted after use
  integer, parameter :: matf_temp  = 11
  !> plus or equals in: Y += fac * X^tx * Z^tz
  integer, parameter :: matf_plus  = 12
  !> X transposed in  : Y += fac * X^tx * Z^tz
  integer, parameter :: matf_tx    = 13
  !> Y transposed in  : Y += fac * X^tx * Z^tz
  integer, parameter :: matf_ty    = 14

  !> unused bits in %flags should contain a "magic" number.
  !> This is used to detect memory corruption and accidental
  !> (ab)use of undefined matrices, to thereafter fail gracefully
  integer, parameter :: matf_magic_mask = not(2**20-1)
  integer, parameter :: matf_magic_num = 1765 * 2**20


  !> matrix structure
  type matrix
     !> number of rows of each block
     integer          :: nrow
     !> number of columns of each block
     integer          :: ncol
     !> flags telling which variant of matrix this is
     integer          :: flags
     !> main block of elements. All non-zero matrices have this
     real(8), pointer :: elms(:,:) ! elms(nrow,ncol)
     !> temporary matrices (during defop evaluation) have additional
     !> contribution: fac * X^tx * Z^tz, where fac, X and Z are here,
     !> while tx and tz are kept inside %flags
     complex(8)            :: defop_fac
     type(matrix), pointer :: defop_X, defop_Z
  end type


  !> for switching debugging on and off
  logical :: matrix_genop_debug = .false.

  private

contains


  !> whether A is defined, according to A%flags
  function mat_isdef(A)
    type(matrix), intent(in) :: A
    logical :: mat_isdef
    mat_isdef = (iand(A%flags,matf_magic_mask) &
                           == matf_magic_num)
  end function


  !> free A, and reset it to undefined state
  subroutine mat_free(A)
    type(matrix), intent(inout) :: A
    if (.not.mat_isdef(A)) then
       call quit('error: mat_free(A), a undefined')
    else if (btest(A%flags, matf_zero)) then
       if (associated(A%elms)) &
          call quit('error: mat_free(A), expected A%elms to be allocated')
    else if (btest(A%flags, matf_alias) then
       if (.not.associated(A%elms)) &
          call quit('error: mat_free(A), expected A%elms to be deallocated')
    else
       deallocate(A%elms)
    end if
    nullify(A%elms)
    A%nrow  = -huge(1) !will cause error if used
    A%ncol  =  huge(1) !will cause error if used
    A%flags = 0 !undefined, since doesn't have matf_magic_num
    A%defop_fac = 1
    nullify(A%defop_X)
    nullify(A%defop_Z)
  end subroutine


  !> allocate matrix, according to its A%flags, A%nrow, A%ncol,
  !> ASSUMING that A's fields are set correctly
  subroutine mat_alloc(A)
     type(matrix), intent(inout) :: A
     if (associated(A%elms))
        call quit('error: mat_alloc(A) expected A%elms to be deallocated')
     allocate(A%elms(A%nrow, A%ncol))
     ! clear any zero flag
     A%flags = ibclr(A%flags, matf_zero)
  end subroutine


  !> initialize C in the shape of A, A^T, A*B, A^T*B, A*B^T or A^T*B^T,
  !> optionally not allocate and optionally zero the allocated elements
  subroutine mat_init(C, A, B, ta, tb, zero, noalloc)
    type(matrix),        intent(inout) :: C
    type(matrix),           intent(in) :: A
    type(matrix), optional, intent(in) :: B
    logical,      optional, intent(in) :: ta, tb, zero, alloc
    logical :: zer, alc, trp, taa, tbb
    integer :: nrow
    if (present(tb) .and. .not.present(B)) &
       call quit('error: mat_init(C,A,B,ta,tb), tb present without B')
    trp = .false.; if (present(trps))  trp = trps
    alc = .true.;  if (present(alloc)) alc = alloc
    zer = .false.; if (present(zero))  zer = zero
    taa = .false.; if (present(taa))   taa = ta
    tbb = .false.; if (present(tbb))   tbb = tb
    if (mat_isdef(C)) call mat_free(C)
    C%nrow = merge(A%ncol, A%nrow, taa)
    if (.not.present(B)) C%ncol = merge(A%nrow, A%ncol, taa)
    if (     present(B)) C%ncol = merge(B%nrow, B%ncol, tbb)
    ! inherit flags from A
    C%flags = ior(iand(A%flags,not(matf_magic_mask)), &
                                   matf_magic_num)
    C%defop_fac = 0
    nullify(C%defop_X)
    nullify(C%defop_Z)
    if (alc) call mat_alloc(C)
    if (alc .and. zer) C%elms = 0
    if (.not.alc) C%flags = ibset(C%flags, mat_unalc)
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


   !> copy all fields of A into B
   subroutine mat_dup(A, B)
     type(matrix), intent(in)    :: A
     type(matrix), intent(inout) :: B
     B = A
   end subroutine


   !> copy all fields of A into B
   function mat_isdup(A, B)
     type(matrix), intent(in) :: A
     type(matrix), intent(in) :: B
     logical :: mat_isdup
     mat_isdup = associated(B%elms, A%elms)
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
      else if (tx .and. XisY)
         Y%elms = a * transpose(Y%elms)
      else if (tx)
         Y%elms = a * transpose(X%elms)
      else if (py .and. XisY)
         Y%elms = Y*elms + a * Y%elms
      else if (py)
         Y%elms = Y*elms + a * X%elms
      else if (XisY)
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
      else if (ta .and. pc)
         C%elms = C%elms + fac * matmul(transpose(A%elms), B%elms)
      else if (ta)
         C%elms = fac * matmul(transpose(A%elms), B%elms)
      else if (tb .and. pc)
         C%elms = C%elms + fac * matmul(A%elms, transpose(B%elms))
      else if (tb)
         C%elms = fac * matmul(A%elms, transpose(B%elms))
      else if (pc)
         C%elms = C%elmx + fac * matmul(A%elms, B%elms)
      else
         C%elms = fac * matmul(A%elms, B%elms)
      end if
   end subroutine



   !> matrix dot product: sum_ij Aij^* Bij (ta=F)
   !>  or  product trace: sum_ij Aji Bij   (ta=F)
   !> A and B are assumed initialied and with equal/transpose shapes
   function mat_dot(A, B, t)
      type(matrix), intent(in) :: A, B
      logical,      intent(in) :: t !transpose
      complex(8)               :: mat_dot
      if (t) then
         mat_dot = sum(A%elms * B%elms)
      else
         mat_dot = sum(transpose(A%elms) * B%elms)
      end if
      ! if closed-shell, then additional factor 2
      if (btest(A%flags, matf_clsh) .or. &
          btest(B%flags, matf_clsh)) mat_dot = 2*mat_dot
   end function



   !> print matrix A to optional unit, optionally starting with label, &
   !> optionally making each column 'width' wide, with optional 'braces'
   subroutine mat_print(A, label, unit, width, braces)
      type(matrix),           intent(in) :: A
      character(*), optional, intent(in) :: label, braces
      integer,      optional, intent(in) :: unit, width
      integer                            :: uni, colw, dec, i, j, siz
      character(8)                       :: fmt
      character(4)                       :: brac
      ! process optional argument unit, which defaults to stdout
      uni = 6
      if (present(unit)) uni = unit
      ! set d to the largest number of digits to be printed
      ! before the decimal point (including - signs)
      dec = max(max(1, ceiling(log10( maxval(A%elms)))), &
                max(2, ceiling(log10(-minval(A%elms)))+1))
      ! process optional width
      colw = 9
      if (present(width)) colw = max(width,dec+2) !max, to avoid stars *****
      ! set dec to the number of decimals to be printed
      dec = colw - dec - 1
      ! argument label is optional. If present, print that
      if (present(label)) write (uni,'(a)') label
      ! process optional argument braces, defaulting to none
      brac = '    '
      if (present(braces)) then
         siz = len(braces)
         if (siz<1 .or. siz>4) &
            call quit('matrix_genop.mat_print: Argument braces has wrong length')
         brac = braces
         if (siz==1) brac(2:2) = ','
         if (siz<=2) &
            brac(3:3) = merge( ')', merge( '}', &
                        merge(']', brac(1:1), brac(1:1)=='['), &
                              brac(1:1)=='{' ), brac(1:1)=='(' )
         if (siz<=3) brac(4:4) = ','
      end if
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
         character(ncol*(colw+1)+4) :: line
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


end module
