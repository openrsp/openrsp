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
!>     product trace                tr(A,B),
!>     for zeroing and freeing      A=0, A(:)=0, A(:,:)=0
!>     norm ie. sqrt(dot(A,A))      norm(A)
!>     mat_to_full                  A = B(:,:)
!>     mat_set_from_full            B(:,:) = A, B(:) = A
!>
!> To increase performance and lower memory needs, formulas are not
!> evaluated by executing the operators immediately, but the operators are
!> first collected into "batches" of the forms
!>
!>     axtpby:     Z = a * X(^T) + b * Y   (where Z may be Y)
!>     gemm:       Y = a * X(^T) * Z(^T) + b * Y
!>
!> To create matrices, mat_init (from module matrix_operations) is replaced
!> by a more general init_mat. This is to hide context-specific information:
!> closed/open shell, real/complex, 1/2/4-component, full/sparse,
!> and thus make it possible to write general code:
!>
!>     call init_mat(A,nrow,ncol)       original version
!>     call init_mat(A,B)               same shape
!>     call init_mat(A,B,'T')           transpose shape
!>     call init_mat(A,B,'T',C,'N')     rows and columns from diff. mats.
!>     call init_mat(A,B,zero=.true.)   non-allocated zero of same shape
!>     call init_mat(A,reset=.true.)    reset to defaults, nullifying pointers
!>     call init_mat(A,B,alias='FF')    copy pointers, no new allocation
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

   !ajt Many linsca modules propagate mat_init (no private statement).
   !         To avoid collision, rename it to init_mat
   use matrix_genop, &
         init_mat => mat_init

   implicit none

   ! ajt LSDALTON has replaced the (global) quit with lsquit
   !     with unit (lupri) as extra argument, which doesn't
   !     exist in DIRAC. For now, this macro gets around that.
#ifdef LSDALTON_ONLY
#define quit(msg) lsquit(msg,-1)
#endif

   public matrix
   public init_mat
   public print_mat
   public isdef
   public iszero
   public assignment(=)
   public operator(+)
   public operator(-)
   public operator(*)
   !public operator(/)
   public dag
   public dot
   public tr
   public norm
   public re
   public im
   public mat_fix_arg
   public mat_unfix_arg
   public mat_fix_result
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

   interface dag
      module procedure mat_dagger
   end interface

   interface dot
      module procedure mat_dot_prod
   end interface

   interface tr
      module procedure mat_trace
      module procedure mat_prod_trace
   end interface

   interface norm
      module procedure mat_norm
   end interface

   interface re
      module procedure mat_real_part
   end interface

   interface im
      module procedure mat_imag_part
   end interface

   !> specification of A = fx*X^tx * Z^tz + Y    gemm
   !>        or        A = fx*X^tx * Z^tz        gemm
   !>        or        A = fx*X^tx        + Y    axtpby
   !>        or        A = fx*X^tx               axt or by
   !> A points to the matrix where result is to be stored, which is undef.
   !> Flags dx,dy,dz decides whether arguments whould be deleted after use
   type tmp_mat
      type(matrix),pointer :: A,  X,  Z,  Y !associated or not- pointers
      complex(8)           ::    fx         !complex scale factor
      logical              ::    tx, tz     !transpose X or Z ?
      logical              ::    dx, dz, dy !delete X, Y or Z after use ?
   end type

   !> tmp(:) contains a 'stack' of temporary matrices
   integer,parameter         :: max_tmp = 64
   type(tmp_mat),target,save :: tmp(max_tmp)
   integer                   :: top_tmp = 0 !0 as "before index 1, thus empty"

   !> for switching debugging on and off
   logical :: matrix_defop_debug = .false.

contains

   !> if A is defined, return null, if A is undefined and A%elms
   !> points to within tmp(1:max_tmp), return A%elms
   function find_tmp(A,errmsg) result(T)
      type(matrix),target,  intent(in) :: A
      character(*),optional,intent(in) :: errmsg
      type(tmp_mat),pointer            :: T
      if (A%nrow==-1 .and. A%ncol/=-1 .and. &
                   .not.associated(A%elms)) then
         T => tmp(A%ncol)
         T%A => A
         T%A%ncol = -1
      else if (isdef(A)) then
         nullify(T) !not found, not temporary
      else if (present(errmsg)) then
         call quit(errmsg)
      end if !otherwise, just return null
   end function

   !> record the temporary, non-zero matrix M in tmp(:)
   function new_tmp(A, fx, X, tx, dx, Z, tz, dz, fy, Y, dy) result(T)
      type(matrix),        target,intent(inout) ::  A
      type(matrix), target, optional,intent(in) ::  X,  Z,  Y
      complex(8),           optional,intent(in) :: fx,     fy
      logical,              optional,intent(in) :: tx, tz
      logical,              optional,intent(in) :: dx, dz, dy
      type(tmp_mat),pointer :: T
      integer               :: i
      if (top_tmp==max_tmp) then
         print '(a)','error: tmp(1:max_tmp) in linears/matrix_defop.f90 is full'
         print '(a,i2,a)','       increase max_tmp=',max_tmp,' and recompile'
         call quit('matrix new_tmp : tmp(1:max_tmp) is full')
      end if
      do i=1,top_tmp+1
         if (i==top_tmp) cycle !this i taken
         if (i==top_tmp+1) top_tmp = top_tmp+1
         T => tmp(i)
         if (.not.associated(T%A)) exit
      end do
      T%A => A
      nullify(T%X);   if (present(X))  T%X  => X
      nullify(T%Z);   if (present(Z))  T%Z  => Z
      nullify(T%Y);   if (present(Y))  T%Y  => Y
      T%fx = 1d0;     if (present(fx)) T%fx = fx
      T%tx = .false.; if (present(tx)) T%tx = tx
      T%tz = .false.; if (present(tz)) T%tz = tz
      T%dx = .false.; if (present(dx)) T%dx = dx
      T%dz = .false.; if (present(dz)) T%dz = dz
      T%dy = .false.; if (present(dy)) T%dy = dy
      !hack: point A%ncol to T,
      if (.not.isdef(A)) then
         A%ncol = i      !so we can find T when we
         nullify(A%elms) !eventually want to use A
      end if
   end function

   !> remove tmp_mat T from within tmp(:), deallocate t%A if needed
   subroutine del_tmp(T)
      type(tmp_mat),intent(inout) :: T
      nullify(T%A)
      !update last occupied index
      do while (top_tmp/=0)
         if (associated(tmp(top_tmp)%A)) exit
         top_tmp = top_tmp - 1
      end do
   end subroutine

   !> initialize T%A to a zero matrix of the correct shape,
   !> then delete (or not) T%X,T%Z,T%Y, then del_tmp(T)
   subroutine zero_tmp(T)
      type(tmp_mat),intent(inout) :: T
      !if (matrix_defop_debug) print *,'zero_tmp',loc(T%A),T%A%nrow,T%A%ncol
      if (T%A%nrow==-1 .and. T%A%ncol/=-1) T%A%ncol = -1
      if (isdef(T%A)) then
         if (.not.iszero(T%A)) call mat_free(T%A)
      end if
      if (associated(T%Z)) then
         call init_mat(T%A,T%X,merge('T','N',T%tx) &
                          ,T%Z,merge('T','N',T%tz),zero=.true.)
      else
         call init_mat(T%A,T%X,merge('T','N',T%tx),zero=.true.)
      end if
      if (T%dx) call mat_free(T%X)
      if (T%dz) call mat_free(T%Z)
      if (T%dy) call mat_free(T%Y)
      call del_tmp(T)
   end subroutine


   !> evaluate the expression A = f*X*Z + Y, as specified in the fields in T,
   !> by executing either an axtpby or a gemm operation. Then remove T from tmp(:)
   subroutine eval_tmp(T)
      type(tmp_mat),intent(inout) :: T
      logical             :: Azero, complex, fits
      type(matrix),target :: dupA
      !if (matrix_defop_debug) print *,'eval_tmp',loc(T%A),' =',loc(T%X) &
      !                                     ,' *',loc(T%Z),' +',loc(T%Y)
      !if (matrix_defop_debug) print *, 'fx=', T%fx, ' tx=', T%tx, ' tz=', T%tz, &
      !                                ' dx=', T%dx, ' dz=', T%dz, ' dy=', T%dy 
      !hack: remove the T-reference from A%ncol
      if (T%A%nrow==-1 .and. T%A%ncol/=-1 &
          .and..not.associated(T%A%elms)) T%A%ncol = -1
      Azero = .not.isdef(T%A)             !two steps, to prevent
      if (.not.Azero) Azero = iszero(T%A) !iszero(undef) error
      !if A is among the operands X, Z, Y, move to dupA and clear A
      if (associated(T%A,T%X) .or. associated(T%A,T%Z) &
                              .or. associated(T%A,T%Y)) then
         call init_mat(dupA, T%A, alias='FF')
         call init_mat(T%A, reset=.true.) !clear A
         if (associated(T%A,T%X)) T%dx = .true.
         if (associated(T%A,T%X)) T%X => dupA
         if (associated(T%A,T%Z)) T%dz = .true.
         if (associated(T%A,T%Z)) T%Z => dupA
         if (associated(T%A,T%Y)) T%dy = .true.
         if (associated(T%A,T%Y)) T%Y => dupA
      end if
      !decide whether the result should be complex or not
      complex = (dimag(T%fx)/=0 .or. T%X%complex)
      if (associated(T%Z)) complex = (complex .or. T%Z%complex)
      if (associated(T%Y)) complex = (complex .or. T%Y%complex)
      !see whether the result can be placed in Y
      fits = .false.
      if (associated(T%Y)) &
         fits = (T%dy .and. (T%Y%complex .eqv. complex) &
                      .and. .not.associated(T%X,T%Y)    &
                      .and. .not.associated(T%Z,T%Y))
      !if it doesn't fit in Y, see whether the result can be placed in X
      if (.not.fits .and. .not.T%tx .and. .not.associated(T%Z)) then
         fits = (T%dx .and. (T%X%complex .eqv. complex) &
                 .and. dimag(T%fx)==0)
         if (fits) then !if result fits in X, swap X and Y
            ! apply any scale factor to X before swapping
            if (T%fx/=1) call mat_axtpy(T%fx, T%X, .false., .false., T%X)
            ! swap X and Y using Z as temporary
            T%Z => T%X;   T%X => T%Y;   T%Y => T%Z
            T%dz = T%dx;  T%dx = T%dy;  T%dy = T%dz
            nullify(T%Z); T%dz = .false.; T%fx = 1
         end if
      end if
      !if needed, reallocate A
      if (fits) then
         if (.not.Azero) call mat_free(T%A)
         call init_mat(T%A, T%Y, alias='FF')
         if (T%dy) call init_mat(T%Y, reset=.true.)
         T%dy = .false.
      else if (associated(T%Y)) then
         if (.not.Azero) fits = ((T%A%complex .eqv. complex)  &
                                 .and. mat_same(T%A, T%Y, 'N'))
         if (.not.Azero .and. .not.fits) call mat_free(T%A)
         if (.not.fits) call init_mat(T%A, T%Y, 'N', complex=complex)
         !copy Y to the newly allocated A
         call mat_axtpy((1d0,0d0), T%Y, .false., .false., T%A)
      else if (associated(T%Z)) then
         if (.not.Azero) fits = ((T%A%complex .eqv. complex)       &
                .and. mat_same(T%A, T%X, 'R', merge('C','R',T%tx)) &
                .and. mat_same(T%A, T%Z, 'C', merge('R','C',T%tz)))
         if (.not.Azero .and. .not.fits) call mat_free(T%A)
         if (.not.fits) call init_mat( T%A, T%X, merge('T','N',T%tx), &
                                            T%Z, merge('T','N',T%tz), &
                                                 complex=complex )
      else !associated(T%X)
         if (.not.Azero) fits = ((T%A%complex .eqv. complex) .and. &
                            mat_same(T%A, T%X, merge('T','N',T%tx)))
         if (.not.Azero .and. .not.fits) call mat_free(T%A)
         if (.not.fits) call init_mat( T%A, T%X, merge('T','N',T%tx), &
                                                 complex=complex )
      end if
      !execute the requested operation, gemm or axtpy, if any
      if (associated(T%Z)) then
         call mat_gemm(T%fx, T%X, T%tx, T%Z, T%tz, associated(T%Y), T%A)
      else if (associated(T%X)) then
         call mat_axtpy(T%fx, T%X, T%tx, associated(T%Y), T%A)
      end if
      !delete X and/or Z, as specified. Y is never deleted
      if (T%dx .and. .not.associated(T%X,T%Z)  &
               .and. .not.associated(T%X,T%Y)) &
         call mat_free(T%X)
      if (T%dz .and. .not.associated(T%Z,T%Y)) &
         call mat_free(T%Z)
      if (T%dy) call mat_free(T%Y)
      call del_tmp(T)
   end subroutine


   subroutine mat_eq_mat(A,B)
   !denne skal altid tømme tmp(:)
      type(matrix),target,intent(inout) :: A
      type(matrix),       intent(in)    :: B
      type(tmp_mat),pointer             :: T
      logical                           :: defB
      !if (matrix_defop_debug) print *, 'mat_eq_mat', loc(A), ' =', &
      !                  loc(B), B%nrow, B%ncol, associated(B%elms)
      T => find_tmp(B, 'matrix A = B : B undefined')
      if (associated(T)) then
         T%A => A
         call eval_tmp(T)
      else if (mat_same(A, B)) then
         !nothing
      else if (iszero(B)) then
         T => new_tmp(A, X=B)
         call zero_tmp(T)
      else
         T => new_tmp(A, X=B)
         call eval_tmp(T)
      end if

#ifdef PRG_DIRAC
      call copy_matrix_sym(A, B)
#endif

   end subroutine


   subroutine gen_mat_plus_mat(A, plus, B, C)
      type(matrix), target, intent(in) :: A, B
      type(matrix), target             :: C
      type(tmp_mat), pointer           :: TA, TB, TC
      logical                          :: zeroA, zeroB, swap, ltmp, plus
      !if (matrix_defop_debug) print *, 'mat_plus_minus_mat', loc(A), &
      !                                 merge(' + ',' - ',plus), loc(B), ' => ', loc(C)
      !if (matrix_defop_debug) print *, 'A%nrow=', A%nrow, ' A%ncol=', A%ncol, &
      !                                 ' A%elms=', associated(A%elms)
      !if (matrix_defop_debug) print *, 'B%nrow=', B%nrow, ' B%ncol=', B%ncol, &
      !                                 ' B%elms=', associated(B%elms)
      TA => find_tmp(A, 'matrix A + B : A undefined')
      TB => find_tmp(B, 'matrix A + B : B undefined')
      zeroA = .false.; zeroB = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      if (.not.associated(TB)) zeroB = iszero(B)
      !below, the matrix C will be created as C = 1*B^N + A.
      !Decide here whether some operation can be saved by swapping A and B first
      if (associated(TA) .and. associated(TB)) then
         swap = (TB%fx==merge(1,-1,plus) .and. .not.TB%tx &
                              .and. .not.associated(TB%Z) &
                              .and. .not.associated(TB%Y) &
                              .and. .not.associated(TA%Y))
      else if (associated(TA)) then
         swap = (zeroB .or. .not.associated(TA%Y))
      else
         swap = (zeroB .and. .not.zeroA)
      end if
      !create the temporary matrix C = 1*B^N + A (or C = 1*A^N + B)
      if (swap) then
         ltmp = zeroA; zeroA = zeroB; zeroB = ltmp
         TC  => TA;    TA   => TB;    TB   => TC
         TC => new_tmp(C, fx=(1d0,0d0), X=A, dx=associated(TB), &
                                        Y=B, dy=associated(TA))
      else
         TC => new_tmp(C, fx=merge(1,-1,plus)*(1d0,0d0), &
                                 X=B, dx=associated(TB), &
                                 Y=A, dy=associated(TA))

      end if
      ! if A = 1*X^N (simple temporary) promote into C%Y, otherwise evaluate A
      if (associated(TA)) then
         if (.not.TA%tx     .and. .not.associated(TA%Z) &
             .and. TA%fx==1 .and. .not.associated(TA%Y)) then
            TC%dy = TA%dx
            TC%Y => TA%X
            call del_tmp(TA)
         else
            call eval_tmp(TA)
         end if
      end if
      !if B = fx*X*Z (temporary without Y) promote into C%X,C%Z,
      !otherwise evaluate B
      if (associated(TB)) then
         if (.not.associated(TB%Y)) then
            TC%fx = TC%fx*TB%fx
            TC%tx = TB%tx; TC%tz = TB%tz
            TC%dx = TB%dx; TC%dz = TB%dz
            TC%X => TB%X;  TC%Z => TB%Z
            call del_tmp(TB)
         else
            call eval_tmp(TB)
         end if
      end if
      !verify that shapes do indeed match
      if (associated(TC%Z)) then
         if (.not.mat_same(TC%X, TC%Y, 'R',merge('C','R',TC%tx)) .or. &
             .not.mat_same(TC%Z, TC%Y, 'C',merge('R','C',TC%tz)))    &
            call quit('matrix A +/- B : different shapes')
      else
         if (.not.mat_same(TC%X, TC%Y, merge('T','N',TC%tx))) &
            call quit('matrix A +/- B : different shapes')
      end if
      !if either of A and/or B are zero, remove from TC
      if (zeroA.and.zeroB) then
         call zero_tmp(TC)
      else if (zeroA) then
         nullify(TC%Y)
         TC%dy = .false.
      else if (zeroB) then
         TC%fx = 1
         TC%X => TC%Y;  nullify(TC%Y)
         TC%tx = .false.
         TC%dx = TC%dy; TC%dy = .false.
      end if

#ifdef PRG_DIRAC
      call copy_matrix_sym(C, B)
#endif

   end subroutine


   !> verify that shapes match, check for zeros
   function mat_plus_mat(A, B) result(C)
      type(matrix), target, intent(in) :: A, B
      type(matrix), target             :: C
      call gen_mat_plus_mat(A, .true., B, C)
   end function


   !> verify that shapes match, check for zeros
   function mat_minus_mat(A, B) result(C)
      type(matrix), target, intent(in) :: A, B
      type(matrix), target             :: C
      call gen_mat_plus_mat(A, .false., B, C)
   end function


   !> verify that shapes match, check for zeroes
   function mat_times_mat(A, B) result(C)
      type(matrix), target, intent(in) :: A, B
      type(matrix), target             :: C
      type(tmp_mat), pointer           :: TA, TB, TC
      logical                          :: zeroA, zeroB
      !if (matrix_defop_debug) print *, 'mat_times_mat', loc(A), &
      !                            ' *', loc(B), ' => ', loc(C)
      TA => find_tmp(A, 'matrix A * B : A undefined')
      zeroA = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      TB => find_tmp(B, 'matrix A * B : B undefined')
      zeroB = .false.
      if (.not.associated(TB)) zeroB = iszero(B)
      TC => new_tmp(C, X=A, dx=associated(TA), Z=B, dz=associated(TB))
      if (associated(TA)) then
         if (.not.associated(TA%Z) .and. .not.associated(TA%Y)) then
            TC%fx = TC%fx*TA%fx
            TC%tx = TA%tx
            TC%dx = TA%dx
            TC%X => TA%X
            call del_tmp(TA)
         else
            call eval_tmp(TA)
         end if
      end if
      if (associated(TB)) then
         if (.not.associated(TB%Z) .and. .not.associated(TB%Y)) then
            TC%fx = TC%fx*TB%fx
            TC%tz = TB%tx
            TC%dz = TB%dx
            TC%Z => TB%X
            call del_tmp(TB)
         else
            call eval_tmp(TB)
         end if
      end if
      !finally, verify that rows
      if (.not.mat_same(TC%X,TC%Z,merge('R','C',TC%tx),merge('C','R',TC%tz))) &
         call quit('matrix A * B : incompatible shapes')
      !if either of X or Z are zero, zero C
      if (zeroA.or.zeroB) call zero_tmp(TC)

#ifdef PRG_DIRAC
      C%ih_sym = A%ih_sym*B%ih_sym
      C%tr_sym = A%tr_sym*B%tr_sym
      C%irep   = ieor(A%irep, B%irep)
#endif

   end function


   subroutine gen_complex_times_mat(r, A, B)
      complex(8),           intent(in) :: r
      type(matrix), target, intent(in) :: A
      type(matrix), target             :: B
      type(tmp_mat), pointer           :: TA, TB

#ifdef PRG_DIRAC
      call copy_matrix_sym(B, A)
#endif

      !if (matrix_defop_debug) print *,'real_times_mat',r,' *',loc(A),' =>',loc(B)
      TA => find_tmp(A, 'matrix r * A : A undefined')
      if (associated(TA)) then
         if (.not.associated(TA%Y)) then
            !KK&AJT: Quick fix for AIX...
            if(associated(TA%Z)) then
               TB => new_tmp(B, fx=r*TA%fx, X=TA%X, tx=TA%tx, dx=TA%dx, &
                     Z=TA%Z, tz=TA%tz, dz=TA%dz)
            else
               TB => new_tmp(B, fx=r*TA%fx, X=TA%X, tx=TA%tx, dx=TA%dx)
            endif
            call del_tmp(TA)
            if (r==0) call zero_tmp(TB)
            return
         end if
         call eval_tmp(TA)
      end if
      TB => new_tmp(B, X=A, fx=r, dx=associated(TA))
      if (r==0 .or. iszero(A)) call zero_tmp(TB)

   end subroutine

   function complex_times_mat(r, A) result(B)
      complex(8),           intent(in) :: r
      type(matrix), target, intent(in) :: A
      type(matrix), target             :: B
      call gen_complex_times_mat(r, A, B)
   end function

   function real_times_mat(r, A) result(B)
      real(8),            intent(in) :: r
      type(matrix),target,intent(in) :: A
      type(matrix),target            :: B
      call gen_complex_times_mat(r*(1d0,0d0), A, B)
   end function

   function plus_mat(A) result(B)
      type(matrix),target,intent(in) :: A
      type(matrix),target            :: B
      call gen_complex_times_mat((1d0,0d0), A, B)
   end function

   function minus_mat(A) result(B)
      type(matrix),target,intent(in) :: A
      type(matrix),target            :: B
      call gen_complex_times_mat((-1d0,0d0), A, B)
   end function


   !> Transpose the matrix A. If A is temporary and A=a*X or A=a*X*Z,
   !> only swap X and Z, and slip-swap tx and tz.
   !> Otherwise, evaluate A (if tmp), and create a new temporary for B
   function mat_dagger(A) result(B)
      type(matrix),target,intent(in) :: A
      type(matrix),target            :: B
      type(tmp_mat),pointer          :: TA, TB
      !if (matrix_defop_debug) print *,'mat_dagger',loc(A),' => ',loc(B)
      TA => find_tmp(A,'matrix dag(A) : A undefined')
      if (.not.associated(TA)) then
         TB => new_tmp(B, X=A, tx=.true.)
         if (iszero(A)) call zero_tmp(TB)
      else if (associated(TA%Y)) then
         call eval_tmp(TA)
         TB => new_tmp(B, X=A, tx=.true., dx=.true.)
      else if (associated(TA%Z)) then
         TB => new_tmp(B, fx=TA%fx, X=TA%Z, tx=.not.TA%tz, dx=TA%dz, &
                                    Z=TA%X, tz=.not.TA%tx, dz=TA%dx)
         call del_tmp(TA)
      else
         TB => new_tmp(B, fx=TA%fx, X=TA%X, tx=.not.TA%tx, dx=TA%dx)
         call del_tmp(TA)
      end if

#ifdef PRG_DIRAC
      call copy_matrix_sym(B, A)
#endif

   end function


   function gen_mat_dot(A, B, trps) result(r)
   !back-end for mat_dot and mat_trab (dot(A,B) and tr(A,B), resp.)
      type(matrix),target,intent(in) :: A, B
      complex(8)                        :: r
      type(tmp_mat),pointer             :: TA, TB
      logical                           :: zeroA, zeroB, trps, delA, delB
      !if (matrix_defop_debug) print *, 'mat_prod_trace/mat_dot_prod' , &
      !                                 loc(A), loc(B), A%complex, B%complex
      TA => find_tmp(A,'matrix dot/tr(A,B) : A undefined')
      zeroA = .false.
      if (.not.associated(TA)) zeroA = iszero(A)
      TB => find_tmp(B,'matrix dot/tr(A,B) : B undefined')
      zeroB = .false.
      if (.not.associated(TB)) zeroB = iszero(B)
      r = merge(0, 1, zeroA.or.zeroB)
      delA = associated(TA)
      delB = associated(TB)
      if (associated(TA)) then
         if (.not.associated(TA%Z) .and. &
             .not.associated(TA%Y)) then
            r = r * TA%fx
            trps = (trps.neqv.TA%tx)
            call init_mat(TA%A, TA%X, alias='FF')
            call init_mat(TA%X, reset=.true.)
            delA = TA%dx
            call del_tmp(TA)
         else
            call eval_tmp(TA)
         end if
      end if
      if (associated(TB)) then
         if (.not.associated(TB%Z) .and. &
             .not.associated(TB%Y)) then
            r = r * TB%fx
            trps = (trps.neqv.TB%tx)
            call init_mat(TB%A, TB%X, alias='FF')
            call init_mat(TB%X, reset=.true.)
            delB = TB%dx
            call del_tmp(TB)
         else
            call eval_tmp(TB)
         end if
      end if
      !verify that shapes match
      if (.not.mat_same(A, B, merge('T','N',trps))) &
         call quit('matrix dot/tr(A,B) : different shapes')
      !calculate only if prefactor nonzero. "* 2" is due to identical
      if (r/=0) r = r * mat_dot(A, B, trps) * 2 !alpha and beta blocks
      !release any temporary operands
      if (delA) call mat_free(TA%A)
      if (delB) call mat_free(TB%A)
   end function


   function mat_dot_prod(A, B) result(r)
      type(matrix),target,intent(in) :: A, B
      complex(8)                        :: r
      !print *,'mat_dot_prod',loc(A),loc(B)
      r = gen_mat_dot(A, B, (.false.))
   end function

   function mat_prod_trace(A, B) result(r)
      type(matrix),target,intent(in) :: A, B
      complex(8)                        :: r
      !print *,'mat_prod_trace',loc(A),loc(B)
      r = gen_mat_dot(A, B, (.true.))
   end function

   function mat_trace(A) result(r)
      type(matrix),target,intent(inout) :: A
      real(8)                           :: r
      type(tmp_mat),pointer             :: T
      T => find_tmp(A,'matrix r * A : A undefined')
      if (associated(T)) then
         if (associated(T%Z)) then
            if ( .not.mat_same(T%X, T%Z, merge('N','T',T%tx.neqv.T%tz), &
                                         merge('C','R',T%tx)) ) &
               call quit('matrix tr(A) : A not square')
            r = T%fx * mat_dot(T%X, T%Z, T%tx.neqv.T%tz)
         else
            if (.not.mat_same(T%X,T%X)) &
               call quit('matrix tr(A) : A not square')
            r = T%fx * mat_trace_nontmp(T%X)
         end if
         if (associated(T%Y)) &
            r = r + mat_trace_nontmp(T%Y)
         call zero_tmp(T) !this leaves A as zero
         call init_mat(A, reset=.true.) !undefine A
      else if (.not.iszero(A)) then
         r = mat_trace_nontmp(A)
      else
         r = 0.0d0
      end if
   end function

   function mat_norm(A) result(r)
      type(matrix),target,intent(inout) :: A
      real(8)                           :: r
      !print *,'mat_norm',loc(A)
      r = sqrt(gen_mat_dot(A,A,(.false.)))
   end function

   subroutine mat_eq_zero(A,zero)
      type(matrix),intent(inout) :: A
      integer,intent(in)         :: zero
      type(matrix)               :: B
      !if (matrix_defop_debug) print *,'mat_eq_zero A=',loc(A),' elms=',loc(A%elms)
      if (zero/=0) call quit('matrix A = z: z must be 0')
      if (.not.isdef(A)) call quit('A = 0: A undefined')
      call init_mat(B, A, alias='FF')
      call init_mat(A, B, zero=.true.)
      call mat_free(B)
   end subroutine

   subroutine mat_eq_zero_1D(A,zero)
      type(matrix),intent(inout) :: A(:)
      integer,intent(in)         :: zero
      integer                    :: i
      do i = 1, size(A)
         A(i) = 0
      end do
   end subroutine

   subroutine mat_eq_zero_2D(A,zero)
      type(matrix),intent(inout) :: A(:,:)
      integer,intent(in)         :: zero
      integer                    :: i,j
      do j = 1, size(A,2)
         do i = 1, size(A,1)
            A(i,j) = 0
         end do
      end do
   end subroutine

   subroutine mat_eq_zero_3D(A,zero)
      type(matrix),intent(inout) :: A(:,:,:)
      integer,intent(in)         :: zero
      integer                    :: i,j,k
      do k = 1, size(A,3)
         do j = 1, size(A,2)
            do i = 1, size(A,1)
               A(i,j,k) = 0
            end do
         end do
      end do
   end subroutine

   subroutine mat_eq_zero_4D(A,zero)
      type(matrix),intent(inout) :: A(:,:,:,:)
      integer,intent(in)         :: zero
      integer                    :: i,j,k,l
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

   subroutine print_mat(A, unit, colwidth, label, braces)
      type(matrix)                       :: A
      integer,      optional, intent(in) :: unit, colwidth
      character(*), optional, intent(in) :: label, braces
      type(tmp_mat), pointer             :: T
      logical                            :: delA!,defA
      T => find_tmp(A, 'print_mat(A) : A undefined')
      delA = associated(T)
      if (associated(T)) call eval_tmp(T)
      call mat_print_nontmp(A, label, unit, colwidth, braces)
      if (delA) call mat_free(A)
   end subroutine

   subroutine mat_fix_arg(A)
      type(matrix), intent(inout) :: A
      type(tmp_mat), pointer      :: T
      !no action on non-temporary arguments
      if (isdef(A)) return
      T => find_tmp(A, 'mat_fix_arg(A) : A undefined')
      call eval_tmp(T)
      T => new_tmp(A) !hide A in tmp(:)
   end subroutine

   subroutine mat_unfix_arg(A)
      type(matrix), intent(inout) :: A
      type(tmp_mat), pointer      :: T
      T => find_tmp(A) !look for A in tmp
      if (.not.associated(T)) return !it wasn't previously fixed
      call zero_tmp(T)
      call init_mat(A, reset=.true.)
   end subroutine

   subroutine mat_fix_result(A)
      type(matrix), intent(inout) :: A
      type(tmp_mat), pointer      :: T
      type(matrix), save, target  :: fixA
      if (.not.isdef(A)) &
         call quit('mat_fix_result(A) : A undefined')
      if (isdef(fixA)) &
         call quit('mat_fix_result(A) : another fixed result already exists')
      call init_mat(fixA, A, alias='FF') !diplicate A in fixA
      call init_mat(A, reset=.true.)      !reset/clear/undefine A
      T => new_tmp(A, X=fixA, dx=.true.)   !temporary A = 1*fixA, fixA to be deleted
   end subroutine

end module
