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
   public dag
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


   !> for switching debugging on and off
   logical :: matrix_defop_debug = .false.

contains


   !> initialize a temporary, non-zero matrix
   !>     A = fac * X^tx * Z^tz + Y
   subroutine new_tmp(A, fac, X, tx, Z, tz, Y)
      type(matrix),        target,intent(inout) :: A
      complex(8),           optional,intent(in) :: fac
      type(matrix), target, optional,intent(in) :: X, Z, Y
      logical,              optional,intent(in) :: tx, tz
      
      nullify(A%defop_X)
         if (present(X))  T%X  => X
      nullify(A%defop_Z);   if (present(Z))  T%Z  => Z
      nullify(T%Y);   if (present(Y))  T%Y  => Y
      T%fx = 1d0;     if (present(fx)) T%fx = fx
      T%tx = .false.; if (present(tx)) T%tx = tx
      T%tz = .false.; if (present(tz)) T%tz = tz
   end function


   !> initialize T%A to a zero matrix of the correct shape,
   !> then delete (or not) A%defop_X and A%defop_Z
   subroutine zero_tmp(A)
      type(matrix),intent(inout) :: A
      if (T%A%nrow==-1 .and. T%A%ncol/=-1) T%A%ncol = -1
      if (mat_isdef(A)) then
         
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



   !> evaluate the expression Y += fac * X^tx * Z^tz, as specified in
   !> the fields in Y, by executing either an axpy or a gemm operation.
   subroutine eval_temp(Y)
      type(matrix), intent(inout) :: Y
      type(matrix), target :: dupY
      logical :: zero, fits, plus, hx, hz, tx, tz, dx, dz
      ! extract flags
      zero = btest(Y%flags, matf_zero)
      plus = btest(Y%flags, matf_defop_plus)
      hx = .true. !associated(Y%defop_X) always true
      hz = associated(Y%defop_Z)
      tx = btest(Y%flags, matf_defop_tx))
      tz = btest(Y%flags, matf_defop_tz))
      dx = btest(Y%defop_X%flags, matf_temp))
      dz = .false.
      if (hz) dz = btest(Y%defop_Z%flags, matf_temp))
      Y%flags = ibclr(Y%flags, matf_zero)
      Y%flags = ibclr(Y%flags, matf_defop_plus)
      Y%flags = ibclr(Y%flags, matf_defop_tx)
      Y%flags = ibclr(Y%flags, matf_defop_tz)
      ! if Y is an alias or among the operands X, Z, move to dupY
      if (btest(Y%flags, matf_alias) &
            .or. associated(Y%defop_X, Y) &
            .or. associated(Y%defop_Z, Y)) then
         call mat_dup(Y, dupY)
         Y%flags = ibclr(Y%flags, matf_alias)
         zero = .true.
         if (associated(Y%defop_X, Y)) Y%defop_X => dupY
         if (associated(Y%defop_Z, Y)) Y%defop_Z => dupY
      ! if it is Y = fac*X, and X is temporary, switch X and Y
      else if (.not.plus .and. .not.tx .and. dx .and. .not. hz) then
         dupY%defop_X   => Y%defop_X
         dupY%defop_fac =  Y%defop_fac
         if (.not.zero) call mat_free(Y)
         call mat_dup(Y%defop_X, Y)
         
      end if
      ! see whether the result can fit in Y's current allocation
      fits = (.not.zero .and. .not.associated(T%X,T%Y)    &
                        .and. .not.associated(T%Z,T%Y))
      !if it doesn't fit in Y, see whether the result can be placed in X
      if (.not.fits .and. .not.btest(Y%flags, matf_tx) &
                    .and. .not.associated(Y%defop_Z)) then
         fits = btest(A%defop_X%flags, matf_temp)
         if (fits) then !if result fits in X, swap X and Y
            ! apply any scale factor to X before swapping
            if (T%fx/=1) call mat_axtpy(T%fx, T%X, .false., .false., T%X)
            ! swap X and Y using Z as temporary
            T%Z => T%X;   T%X => T%Y;   T%Y => T%Z
            T%dz = T%dx;  T%dx = T%dy;  T%dy = T%dz
            nullify(T%Z); T%dz = .false.; T%fx = 1
         end if
      end if
      !if needed, reallocate Y
      if (fits) then
         if (.not.zero) call mat_free(T%A)
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
      if (associated(Y%defop_Z)) then
         call mat_gemm(Y%defop_fac, Y%defop_X, btest(Y%flags, matf_tx), &
                       Y%defop_Y, btest(Y%flags, matf_tz), &
                       btest(Y%flags, matf_plus), Y)
      else if (associated(T%X)) then
         call mat_axtpy(T%fx, T%X, T%tx, associated(T%Y), T%A)
      end if
      ! if temporary, free X and/or Z
      if (associated(Y%defop_X)) then
         if (btest(Y%defop_X%flags, matf_temp)) &
            call mat_free(Y%defop_X)
         nullify(Y%defop_X)
      end if
      if (assiciated(Y%defop_Z)) then
         if (btest(Y%defop_Z%flags, matf_temp)) &
            call mat_free(Y%defop_Z)
         nullify(Y%defop_Z)
      end if
      ! clear remaining flags
      Y%flags = ibclr(Y%flags, matf_temp)
      Y%flags = ibclr(Y%flags, matf_plus)
      Y%flags = ibclr(Y%flags, matf_tx)
      Y%flags = ibclr(Y%flags, matf_tz)
      Y%defop_fac = 1
   end subroutine


   !> 
   subroutine mat_eq_mat(A, B)
      type(matrix), intent(inout) :: A
      type(matrix)                :: B
      logical :: evalB, freeA
      evalB = must_eval(B, 'error: matrix A=B, B undefined')
      if (evalB .and. associated(B%defop_Z)) then
      else if (evalB .and. associated(B%defop_X)) then
      else if (.not.iszero(B)) then
      else
         ! B is zero, so just move it to A
         if (isdef(A)) call mat_free(A)
         call mat_dup(B, A)
         
      end if
         else
         end if
         B ikke temporær og/eller null
det eneste som skal skje her er at A flyttes inn på Y-plassen eller slettes, før B slettes
hvis Y-plassen er et alias til noen annen enn A må den enten kopieres over til A, eller
Y og X-plassene må byttes (tx=F). 

tilfelle A overskrivbar med både Y og X:
      
      if (isdef(A)) then
         freeA = (iszero(A) .and. temp(B))
         if (.not.freeA) freeA = (.not.mat_same(A,B))
         if (.not.freeA) freeA = (.not.mat_isdup(A,B))
         if (.not.freeA) freeA = .not.alias(B)
             !B er ikke alias og plus true
             !ingen Z eller xt og X alias til A
             !inger Z eller xt og X temp) call mat_free(A)
      if (.not.isdef(A) .and. tmp(B) .and. .not.alias(B)) then
         call mat_dup(B, A)
         B%flags = ibset(B%flags, matf_alias)
         call mat_free(B)
         call eval_tmp(A)
         A%flags = ibclr(A%flags, matf_temp)
      else if (isdef(A) .and. tmp(B) .and. .not.alias(B) .and. .not.plus(B)) then
         
      end if
      call eval_temp(A)
hvis B er alias til A kan A gjøres om til alias og frigjøres,
    før B kopieres til A, alias klareres, B frigjøres
hvis B%defop_X er alias til A, og B alias til noe annet, gjør 

inn til eval_tmp skal A = f*BT, A = f*B med B temp og A uallokert, A = f*B med B ikke-temp,
    A = f*B + A, A = f*B*C
      ! if 
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
         C%flags = ibset(C%flags, matf_plus)
         C%defop_fac = merge(1,-1,plus)
         C%defop_X => B
      end if
      if (zeroA .and. btest(A%flags, matf_temp)) call mat_free(A)
      if (zeroB .and. btest(B%flags, matf_temp)) call mat_free(B)
   end subroutine


   !> verify that shapes match, check for zeros
   function mat_plus_mat(A, B) result(C)
      type(matrix), target :: A, B
      type(matrix),        :: C
      call gen_mat_plus_mat(A, .true., B, C)
   end function


   !> verify that shapes match, check for zeros
   function mat_minus_mat(A, B) result(C)
      type(matrix), target :: A, B
      type(matrix)         :: C
      call gen_mat_plus_mat(A, .false., B, C)
   end function



  !> verify that shapes match, check for zeroes
  function mat_times_mat(A, B) result(C)
    type(matrix), target :: A, B
    type(matrix)         :: C
    logical              :: evalA, evalB, zeroA, zeroB
    evalA = must_eval(A, 'error: matrix A*B, A undefined')
    evalB = must_eval(B, 'error: matrix A*B, B undefined')
    zeroA = btest(A%flags, matf_zero)
    zeroB = btest(B%flags, matf_zero)
    if (.not.mat_same(A, B, mul=.true.) &
       call quit('error: matrix A*B, but A%ncol and B%nrow do not match')
    call mat_init(C, A, B, alloc=.false.)
    C%flags = ibset(C%flags, matf_temp)
    if (evalA) call eval_temp(A)
    if (evalB) call eval_temp(B)
    if (zeroA .or. zeroB) then
       if (btest(A%flags, matf_temp)) call mat_free(A)
       if (btest(B%flags, matf_temp)) call mat_free(B)
       C%flags = ibset(C%flags, matf_zero)
    else
       C%defop_fac =  1
       C%defop_X   => A
       C%defop_Y   => B
    end if
  end function



  subroutine gen_complex_times_mat(r, A, B)
    complex(8),   intent(in) :: r
    type(matrix), target     :: A
    type(matrix),            :: B
    logical :: temp
    temp = istemp(A, 'error: matrix r*A, A undefined')
    call mat_init(B, A, alloc=.false.)
    if (temp) call eval_tmp(A)
    if (r==0 .or. iszero(A)) then
       B%flags = ibset(B%flags, matf_zero)
    else
       B%flags = ibset(B%flags, matf_temp)
       B%defop_fac =  r
       B%defop_X   => A
    end if
    if (temp) call mat_free(A)
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
    logical :: temp
    temp = istemp(A, 'error: matrix transpose trps(A), A undefined')
    call mat_init(B, A, ta=.true., alloc=.false.)
    if (iszero(A)) then
       B%flags = ibset(B%flags, matf_zero)
    else
       if (temp) call eval_temp(A)
       B%flags = ibset(B%flags, matf_temp)
       B%flags = ibset(B%flags, matf_defop_tx)
       B%defop_X => A
    end if
  end function



  !> back-end for mat_dot_prod and mat_prod_trace (dot(A,B) and tr(A,B))
  function gen_mat_dot(A, B, trps) result(r)
    type(matrix) :: A, B
    complex(8)   :: r
    logical      :: trps, za, zb, da, db
    da = istemp(A, 'error: matrix dot/tr(A,B), A undefined')
    db = istemp(B, 'error: matrix dot/tr(A,B), B undefined')
    za = btest(A%flags, matf_zero)
    zb = btest(B%flags, matf_zero)
    if (da) then
       !ajt FIXME forward A = fac * X^tx
       call eval_temp(A)
    end if
    if (db) then
       !ajt FIXME forward B = fac * X^tx
       call eval_tmp(TB)
    end if
    ! verify that shapes match
    if (.not.mat_same(A, B, merge('T','N',trps))) &
       call quit('error: matrix dot/tr(A,B), A and B have different shapes')
    ! don't calculate if
    if (      za.or.zb)  gen_mat_dot = 0
    if (.not.(za.or.zb)) gen_mat_dot = mat_dot(A, B, trps)
    ! delete any temporaries
    if (da) call mat_free(A)
    if (db) call mat_free(B)
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
    mat_norm = gen_mat_dot(A, A, (.false.))
    mat_norm = sqrt(dreal(mat_norm))
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
    logical :: temp, zero
    call assert_def(A, 'error: print_mat(A), A undefined')
    temp = btest(A%flags, matf_temp)
    zero = btest(A%flags, matf_zero)
    if (temp) call eval_temp(A)
    if (zero) call mat_init(A, zero=.true.)
    call mat_print_nontmp(A, label, unit, width, braces)
    if (temp .or. zero) call mat_free(A)
  end subroutine


end module
