! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module matrix_defop

!> module matrix_genop contains slightly generalized
!> versions of the operations in DALTON-linsca-ng's module
!> matrix_operations (Dalton)
!>
!> written by Andreas J. Thorvaldsen
!> radovan.bast@uit.no (2009-02-18): generalized for DIRAC
!>
!> mainly intended as hidden back-end for matrix_defop
!> mat_init is generalized so that another matrix can fill
!> the role of nrow and ncol
module matrix_genop

#ifdef DALTON_AO_RSP
  use xmatrix,                               &
    only: matrix,                            &
          mat_init_full    => xmatrix_set,   &
          mat_free_nonzero => xmatrix_clean, &
          mat_mul          => xmatrix_gemm,  &
          mat_zero         => xmatrix_zero,  &
          mat_daxpy        => xmatrix_axpy,  &
          mat_scal         => xmatrix_scal,  &
          mat_trans        => xmatrix_trans, &
          mat_trab         => xmatrix_trdot, &
          mat_dotproduct   => xmatrix_dot
#else
  use matrix_operations,                  &
      only: matrix,                       &
            mat_init_full    => mat_init, &
            mat_free_nonzero => mat_free, &
            mat_mul,                      &
            mat_zero,                     &
            mat_daxpy,                    &
            mat_scal,                     &
            mat_trans,                    &
            mat_trab,                     &
            mat_dotproduct
#endif

   implicit none

   ! ajt LSDALTON has replaced the (global) quit with lsquit
   !     with unit (lupri) as extra argument, which doesn't
   !     exist in DIRAC. For now, this macro gets around that.
#ifdef LSDALTON_ONLY
#define quit(msg) lsquit(msg,-1)
#endif

   public matrix
   public mat_init
   public mat_free
   public mat_same
   public mat_isdef
   public mat_iszero
   public mat_axtpy
   public mat_gemm
   public mat_dot
   public mat_trace_nontmp
   public mat_print_nontmp
   public mat_real_part
   public mat_imag_part
#ifdef PRG_DIRAC
   public copy_matrix_sym
#endif

   private

   interface mat_init
      module procedure mat_init_full
      module procedure mat_init_from
   end interface

#ifdef PRG_DIRAC
#include "dgroup.h"
#include "priunit.h"
#endif /* ifdef PRG_DIRAC */

contains

   subroutine mat_free(A)
      type(matrix),intent(inout) :: A
      !if (associated(A%elms)) print *,'mat_free',loc(A),' elms=',loc(A%elms)
      if (associated(A%elms)) call mat_free_nonzero(A)
      A%nrow = -1; A%ncol = -1
   end subroutine


   subroutine mat_init_from(C, A, ta, B, tb, complex, reset, zero, alias)
      type(matrix),           intent(inout) :: C
      type(matrix), optional, intent(in)    :: A, B
      character*1,  optional, intent(in)    :: ta, tb
      logical,      optional, intent(in)    :: complex, zero, reset
      character*2,  optional, intent(in)    :: alias
      logical :: taa, tbb, zer, res, cplx
      integer :: nrow, ncol
      if (present(ta) .and. .not.present(A)) &
         call quit('mat_init input error: ta cannot be present without A')
      if (present(B) .and. .not.present(ta)) &
         call quit('mat_init input error: B cannot be present without A')
      if (present(tb) .and. .not.present(B)) &
         call quit('mat_init input error: tb cannot be present without B')
      if (present(reset) .and. present(A)) &
         call quit('mat_init input error: reset cannot be present with A')
      if (present(zero) .and. .not.present(A)) &
         call quit('mat_init input error: zero cannot be present without A')
      if (present(alias) .and. .not.present(A)) &
         call quit('mat_init input error: alias cannot be present without A')
      if (present(alias) .and. present(ta)) &
         call quit('mat_init input error: alias cannot be present with ta')
      taa = .false.
      if (present(ta)) then
         if (ta/='N' .and. ta/='T') &
             call quit("mat_init input error: ta is not 'N' or 'T'")
         taa = (ta=='T')
      end if
      tbb = .false.
      if (present(tb)) then
         if (tb/='N' .and. tb/='T') &
             call quit("mat_init input error: tb is not 'N' or 'T'")
         tbb = (tb=='T')
      end if
      res = .false.
      if (present(reset)) res = reset
      zer = .false.
      if (present(zero)) zer = zero
      cplx = .false.
      if (present(A)) cplx = A%complex
      if (present(B)) cplx = (cplx .or. B%complex)
      if (present(complex)) cplx = complex
      C%complex = cplx
      nrow = -1
      if (present(A)) nrow = merge(A%ncol, A%nrow, taa)
      ncol = -1
      if (present(A)) ncol = merge(A%nrow, A%ncol, taa)
      if (present(B)) ncol = merge(B%nrow, B%ncol, tbb)
      C%nrow = nrow
      C%ncol = ncol
      if (zer .or. res) then
         nullify(C%elms)
#ifdef PRG_DIRAC
         C%ih_sym = 1
         C%tr_sym = 1
         C%irep   = 0
#endif
      else if (present(alias)) then
         C = A !copies all fields
         if (alias=='RF') then !real full (alpha and beta)
            if (associated(A%elms)) C%elms => A%elms(1:nrow*ncol)
            C%complex = .false.
         else if (alias=='IF') then !imaginary full (alpha and beta)
            if (associated(A%elms) .and. A%complex) then
                C%elms => A%elms( nrow*ncol+1 : 2*nrow*ncol ) !ajt fixme
            else
                nullify(C%elms)
            end if
            C%complex = .false.
         else if (alias/='FF') then
             call quit("mat_init input error: alias cannot be '" // alias // "'")
         end if
      else
         !print *,'mat_init',loc(C)
         call mat_init_full(C, nrow, ncol)
#ifdef PRG_DIRAC
         if (present(A) .and. .not. present(B)) then
           call copy_matrix_sym(C, A)
         end if
         if (present(B) .and. .not. present(A)) then
           call copy_matrix_sym(C, B)
         end if
         if (present(A) .and.       present(B)) then
           C%ih_sym = A%ih_sym*B%ih_sym
           C%tr_sym = A%tr_sym*B%tr_sym
           C%irep   = ieor(A%irep, B%irep)
         end if
#endif /* ifdef PRG_DIRAC */
      end if
   end subroutine


   function mat_isdef(A)
      type(matrix),intent(in) :: A
      logical                 :: mat_isdef
      mat_isdef = .not.(A%nrow==-1 .and. A%ncol==-1)
   end function


   function mat_iszero(A)
      type(matrix), intent(in) :: A
      logical                  :: mat_iszero
      if (A%nrow==-1 .and. A%ncol==-1) &
         call quit('mat_iszero(A): matrix A undefined')
      mat_iszero = .not.associated(A%elms)
   end function


   function mat_same(A, B, ta, tb)
      type(matrix),          intent(in) :: A, B
      character*1, optional, intent(in) :: ta, tb
      logical                           :: mat_same
      character*1                       :: taa, tbb
      taa = 'N'
      if (present(ta)) taa = ta
      tbb = 'N'
      if (present(tb)) tbb = tb
      if ( taa=='N' .and. tbb=='N' .or. &
           taa=='T' .and. tbb=='T' ) then
         mat_same = (A%nrow==B%nrow .and. A%ncol==B%ncol)
      else if ( taa=='N' .and. tbb=='T' .or. &
                taa=='T' .and. tbb=='N' ) then
         mat_same = (A%nrow==B%ncol .and. A%ncol==B%nrow)
      else if (taa=='R' .and. tbb=='R') then
         mat_same = (A%nrow==B%nrow)
      else if (taa=='R' .and. tbb=='C') then
         mat_same = (A%nrow==B%ncol)
      else if (taa=='C' .and. tbb=='R') then
         mat_same = (A%ncol==B%nrow)
      else if (taa=='C' .and. tbb=='C') then
         mat_same = (A%ncol==B%ncol)
      else
         call quit('mat_same input error : ta='//taa//' tb='//tbb)
      end if
      if (.not.present(ta) .and. .not.present(tb)) &
         mat_same = mat_same .and. associated(A%elms,B%elms)
   end function


   !> private 'driving' routine for mat_axtpy
   !> add/put block ix from X to block iy in Y
   subroutine blk_axtpy(a, ix, X, tx, py, iy, Y)
      real(8),      intent(in) :: a      !scale factor for X
      integer,      intent(in) :: ix, iy !block indices
      logical,      intent(in) :: tx, py !transpose X or not, add to Y or overwrite Y
      type(matrix), intent(in) :: X
      type(matrix), intent(inout) :: Y
      type(matrix)     :: dupX !since X is intent(in)
      real(8), pointer :: Yelms(:)
      logical          :: Ycomplex, Xis0, XisY
      integer          :: siz
      Xis0 = (ix/=1 .and. .not.X%complex)
      XisY = (associated(X%elms, Y%elms) .and. ix==iy)
      if (Xis0) then
         !nothing
      else if (X%complex) then
         dupX = X !copy all fields/pointers, not data
         dupX%complex = .false.
         siz = size(X%elms)/2
         dupX%elms => X%elms(1 + (ix-1)*siz : ix*siz)
      end if
      Ycomplex = Y%complex
      if (Ycomplex) then
         Y%complex = .false.
         Yelms => Y%elms
         siz = size(Yelms)/2
         Y%elms => Yelms(1 + (iy-1)*siz : iy*siz)
      else if (iy==2) then
         return
      end if
      if (Xis0) then !either Y=0 or Y+=0
         if (.not.py) call mat_zero(Y)
      else if (tx .and. a==1 .and. .not.py .and. .not.XisY) then
         call mat_trans( merge(dupX, X, X%complex), Y )
      else if (XisY .and. .not.py) then !scale operation
         if (a==0) call mat_zero(Y)
         if (a/=0) call mat_scal(a, Y)
      else if (.not.XisY) then !axpy operation
         if (.not.py) call mat_zero(Y) !since no out-of-place scal exists
         if (a/=0) call mat_daxpy(a, merge(dupX, X, X%complex), Y)
      end if
      if (Ycomplex) then
         Y%complex = .true.
         Y%elms => Yelms
      end if
   end subroutine


   subroutine mat_axtpy(a, X, tx, py, Y)
      complex(8),   intent(in)    :: a
      type(matrix), intent(in)    :: X
      logical,      intent(in)    :: tx, py
      type(matrix), intent(inout) :: Y
      type(matrix)       :: Z !temporary
      logical            :: XisY, useZ
      logical, parameter :: trps=.true., norm=.false., &
                            plus=.true., eqls=.false.
      XisY = associated(X%elms,Y%elms)
      useZ = (((XisY .or. (tx .and. X%complex)) &
              .and. dimag(a) /= 0) .or. (tx .and. py))
      if (useZ) call mat_init(Z, Y, complex=.false.)
      if (tx .and. py) then
         call blk_axtpy(      1d0, 1, X, trps, eqls, 1, Z) !Z  = X^T
         call blk_axtpy( dimag(a), 1, Z, norm, plus, 2, Y) !Y´+= a´* Z
         call blk_axtpy( dreal(a), 1, Z, norm, plus, 1, Y) !Y += a * Z
         call blk_axtpy(      1d0, 2, X, trps, eqls, 1, Z) !Z  = X´^T
         call blk_axtpy( dimag(a), 1, Z, norm, plus, 1, Y) !Y += a´* Z
         call blk_axtpy(-dreal(a), 1, Z, norm, plus, 2, Y) !Y´-= a * Z
      else if (tx .or. XisY) then
         if (tx) &
            call blk_axtpy(1d0, 1, X, trps, eqls, 1, Y)
         if (tx .and. X%complex) &
            call blk_axtpy(1d0, 2, X, trps, eqls, 2, Y)
         if (dimag(a) /= 0) &
            call blk_axtpy(1d0, 2, Y, norm, eqls, 1, Z) !Z = Y´
                                                       !Y´ = -a * X´^T (or a * Y´)
         call blk_axtpy(dreal(a)*merge(-1,1,tx), 2, Y, norm, eqls, 2, Y)
         call blk_axtpy(dimag(a), 1, Y, norm, plus, 2, Y) !Y´+= a´* X^T (or a´* Y)
         call blk_axtpy(dreal(a), 1, Y, norm, eqls, 1, Y) !Y += a * X^T (or a * Y)
                                                         !Y += a´* X´T (or -a´* Y)
         call blk_axtpy(dimag(a)*merge(1,-1,tx), 1, Z, norm, plus, 1, Y)
      else
         call blk_axtpy( dreal(a), 2, X, norm, py  , 2, Y) !Y´+= a * X´
         call blk_axtpy( dimag(a), 1, X, norm, plus, 2, Y) !Y´+= a´* X
         call blk_axtpy( dreal(a), 1, X, norm, py  , 1, Y) !Y += a * X
         call blk_axtpy(-dimag(a), 2, X, norm, plus, 1, Y) !Y -= a´* X´
      end if
      if (useZ) call mat_free(Z)
   end subroutine


   subroutine blk_gemm(rab, ia, A, ta, ib, B, tb, pc, ic, C)
      real(8),      intent(in)    :: rab
      integer,      intent(in)    :: ia, ib, ic
      logical,      intent(in)    :: ta, tb, pc
      type(matrix), intent(in)    :: B
      type(matrix), intent(inout) :: A, C
      type(matrix)     :: dupB
      logical          :: Acomplex, Ccomplex, Ais0, Bis0
      real(8), pointer :: Aelms(:), Celms(:)
      integer          :: siz
      Ccomplex = C%complex
      if (Ccomplex) then
         C%complex = .false.
         Celms => C%elms
         siz = size(Celms)/2
         C%elms => Celms( 1 + (ic-1)*siz : ic*siz )
      else if (ic==2) then
         return
      end if
      Bis0 = ( ib/=1 .and. .not.B%complex )
      if (.not.Bis0 .and. B%complex) then
         dupB = B !copy all fields
         dupB%complex = .false.
         siz = size(B%elms)/2
         dupB%elms => B%elms( 1 + (ib-1)*siz : ib*siz )
      end if
      Ais0 = ( ia/=1 .and. .not.A%complex )
      Acomplex = A%complex
      if (.not.Ais0 .and. A%complex) then
         A%complex = .false.
         Aelms => A%elms
         siz = size(Aelms)/2
         A%elms => Aelms( 1 + (ia-1)*siz : ia*siz )
      end if
      if (rab==0 .or. Ais0 .or. Bis0) then !either C=0 or C+=0
         if (.not.pc) call mat_zero(C)
      else
         call mat_mul(A, merge(dupB, B, B%complex), merge('T','N',ta), &
                      merge('T','N',tb), rab, merge(1d0, 0d0, pc), C)
      end if
      if (.not.Ais0 .and. Acomplex) then
         A%complex = .true.
         A%elms => Aelms
      end if
      if (Ccomplex) then
         C%complex = .true.
         C%elms => Celms
      end if
   end subroutine


   !> matrix multiply: C = rc*C + rab * A^ta * B^tb
   !> assumed initialied and with corresponding shapes
   subroutine mat_gemm(fab, A, ta, B, tb, pc, C)
      complex(8),      intent(in) :: fab
      type(matrix),    intent(in) :: A, B
      logical,         intent(in) :: ta, tb, pc
      type(matrix), intent(inout) :: C
      type(matrix)       :: X, Y !temporary matrices
      logical, parameter :: trps=.true., norm=.false., &
                            plus=.true., eqls=.false.
      !print *,'dgemm ta=',ta,' tb=',tb,' ',C%nrow,C%ncol,A%ncol,A%nrow,B%ncol,B%nrow
      !call mat_mul(A, B, ta, tb, dreal(fab), dreal(fc), C)
      call mat_init(X, A, complex=.false.)
      call blk_axtpy( dreal(fab), 1, A, norm, eqls, 1, X)
      call blk_axtpy(-dimag(fab), 2, A, norm, plus, 1, X)
      call mat_init(Y, C, complex=.false.)
      call blk_gemm(1d0, 1, X, ta, 2, B, tb, eqls, 1, Y)
      call blk_axtpy( 1d0, 1, Y, norm, pc, 1, C)
      call blk_axtpy( 1d0, 1, Y, norm, pc, 2, C)
      call blk_axtpy( dreal(fab), 2, A, norm, eqls, 1, X)
      call blk_axtpy( dimag(fab), 1, A, norm, plus, 1, X)
      call blk_gemm(1d0, 1, X, ta, 1, B, tb, eqls, 1, Y)
      call blk_axtpy( 1d0, 1, Y, norm, plus, 2, C)
      call blk_axtpy(-1d0, 1, Y, norm, plus, 1, C)
      call mat_free(Y)
      call blk_axtpy( dreal(fab) + dimag(fab), 1, A, norm, eqls, 1, X)
      call blk_axtpy( dreal(fab) - dimag(fab), 2, A, norm, plus, 1, X)
      call mat_init(Y, B, complex=.false.)
      call blk_axtpy( 1d0, 1, B, norm, eqls, 1, Y)
      call blk_axtpy(-1d0, 2, B, norm, plus, 1, Y)
      call blk_gemm(1d0, 1, X, ta, 1, Y, tb, plus, 1, C)
      call mat_free(Y)
      call mat_free(X)
   end subroutine


   !> matrix dot product: sum_ij Aij^* Bij (ta=F)
   !>  or  product trace: sum_ij Aji Bij   (ta=F)
   !> A and B are assumed initialied and with equal/transpose shapes
   function mat_dot(A, B, ta) result(r)
      type(matrix), intent(in) :: A, B
      logical,      intent(in) :: ta
      complex(8)               :: r
      r = blk_dot(1,1) + (0d0,1d0)*merge(1,-1,ta) * blk_dot(2,1) &
        + (0d0,1d0)*blk_dot(1,2) + merge(-1,1,ta) * blk_dot(2,2)
   contains
      function blk_dot(ia, ib) result(r)
         integer :: ia, ib
         real(8) :: r
         type(matrix) :: dupA, dupB
         logical :: Xis0
         integer :: siz
         r = 0
         if (ia/=1 .and. .not.A%complex) return
         if (ib/=1 .and. .not.B%complex) return
         if (A%complex) then
            dupA = A !copy all fields
            dupA%complex = .false.
            siz = size(A%elms)/2
            dupA%elms => A%elms( 1 + (ia-1)*siz : ia*siz )
         end if
         if (B%complex) then
            dupB = B !copy all fields
            dupB%complex = .false.
            siz = size(B%elms)/2
            dupB%elms => B%elms( 1 + (ib-1)*siz : ib*siz )
         end if
         if (ta) then
            r = mat_trab( merge(dupA, A, A%complex), &
                          merge(dupB, B, B%complex)  )
         else
            r = mat_dotproduct( merge(dupA, A, A%complex), &
                                merge(dupB, B, B%complex)  )
         end if
      end function
   end function


   !> matrix sum Aii, for A square and assumed initialied
   function mat_trace_nontmp(A) result(r)
      type(matrix), intent(in) :: A
      real(8)                  :: r
      integer                  :: i
      if (A%nrow/=A%ncol) call quit('mat_trace_nontmp(A) : A not square')
      r = 0d0
      if (associated(A%elms)) &
         r = sum( (/ ( A%elms(1 + (i-1)*(A%nrow+1) ), i=1,A%nrow ) /) )
   end function


   !> print matrix A to unit, starting with label, making each column colwidth wide,
   !> and using braces braces
   subroutine mat_print_nontmp(A, label, unit, colwidth, braces)
      type(matrix),           intent(in) :: A
      integer,      optional, intent(in) :: unit, colwidth
      character(*), optional, intent(in) :: label, braces
      integer                            :: uni, colw, dec, i, j, siz
      character*8                        :: fmt
      character*4                        :: brac
#ifndef PRG_DIRAC
      uni = 6 !process optional argument unit, which defaults to stdout
      if (present(unit)) uni = unit
      !set d to the largest number of digits to be printed
      dec = 1 !before the decimal point (including - signs)
      if (associated(A%elms))                                    &
        dec = max( max( 1, ceiling(log10(maxval(A%elms))) ),     &
                   max( 2, ceiling(log10(-minval(A%elms))) + 1 ) )
      colw = 9 !process default for width
      if (present(colwidth)) colw = max(colwidth,dec+2) !max, to avoid stars *****
      dec = colw - dec - 1 !set d to the number of decimals to be printed
      !argument label is optional. If present, print that
      if (present(label)) write (uni,'(a)') label
      !process optional argument braces, defaulting to none
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
      !create the format string to be used for each element
      fmt = '(fww.dd)'
      write (fmt(3:7),'(i2,a1,i2)') colw, '.', dec
      !call the printing routine
      if (A%complex) then
         siz = size(A%elms)/2
         write (uni,'(a)') '    (real part)'
         call subr(A%nrow, A%ncol, A%elms(1:siz), colw, fmt, brac, uni)
         write (uni,'(a)') '    (imaginary part)'
         call subr(A%nrow, A%ncol, A%elms(siz+1:siz*2), colw, fmt, brac, uni)
      else
         call subr(A%nrow, A%ncol, A%elms, colw, fmt, brac, uni)
      end if
      write (uni,'()') !final empty line
#else
      write(*, *) 'matrix = ', label
      write(*, *) 'ih_sym = ', A%ih_sym
      write(*, *) 'tr_sym = ', A%tr_sym
      write(*, *) 'irep   = ', A%irep

      if (associated(A%elms)) then
        call prqmat(A%elms,            &
                    A%nrow,            &
                    A%ncol,            &
                    A%nrow,            &
                    A%ncol,            &
                    nz,                &
                    ipqtoq(1, A%irep), &
                    lupri)
      else
        write(*, *) 'all elements zero'
      end if
#endif /* #ifndef PRG_DIRAC */

   contains

      subroutine subr(nrow, ncol, elms, colw, fmt, brac, unit)
         integer,     intent(in) :: nrow, ncol, colw, unit
         real(8),     intent(in) :: elms(nrow,ncol)
         character*8, intent(in) :: fmt
         character*4, intent(in) :: brac
         character*( ncol*(colw+1) + 4 ) :: line
         integer :: i, j, l
         do i = 1, nrow
            line(1:1) = merge(brac(1:1),' ',i==1)
            line(2:2) = brac(1:1)
            l = 3
            do j = 1, ncol
               write (line( l : l+colw-1 ), fmt) elms(i,j)
               l = l + colw 
               if (j/=ncol) line(l:l) = brac(2:2)
               if (j/=ncol) l = l+1
            end do
            line(l:l) = brac(3:3)
            line(l+1:l+1) = merge( brac(3:3), brac(4:4), i==nrow )
            l = l + 2
            if (brac=='    ') write (unit,'(a)') line(3:len(line)-2)
            if (brac/='    ') write (unit,'(a)') line
         end do
      end subroutine

   end subroutine


   !> Reference the real part of another matrix (no copy)
   !> Should eventually be replaced by more general mat_init(..)
   function mat_real_part(A) result(B)
      type(matrix), intent(in) :: A
      type(matrix)             :: B
      integer :: siz
      B = A !copy pointers not data
      if (B%complex) then
         siz = size(B%elms)/2
         B%elms => B%elms(1:siz)
         B%complex = .false.
      end if
   end function


   !> Reference the imaginary part of another matrix (no copy)
   !> Should eventually be replaced by more general mat_init(..)
   function mat_imag_part(A) result(B)
      type(matrix), intent(in) :: A
      type(matrix)             :: B
      integer :: siz
      B = A !copy pointers not data
      if (B%complex) then
         siz = size(B%elms)/2
         B%elms => B%elms(siz+1:2*siz)
         B%complex = .false.
      else
         nullify(B%elms) !mark as zero
      end if
   end function

#ifdef PRG_DIRAC
   subroutine copy_matrix_sym(A, B)
     
     type(matrix) :: A, B
     
     A%ih_sym = B%ih_sym
     A%tr_sym = B%tr_sym
     A%irep   = B%irep
   
   end subroutine
#endif

end module
