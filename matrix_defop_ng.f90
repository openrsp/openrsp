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
module matrix_defop_ng

  use matrix_genop_ng, mat_print_nonzero => mat_print, &
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
    do i=1, size(A); A(i)=z; end do
  end subroutine

  subroutine mat_eq_zero_2D(A, z)
    type(matrix), target, intent(inout) :: A(:,:)
    integer,              intent(in)    :: z
    integer j
    do j=1, size(A,2); A(:,j)=z; end do
  end subroutine

  subroutine mat_eq_zero_3D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:)
    integer,              intent(in)    :: z
    integer k
    do k=1, size(A,3); A(:,:,k)=z; end do
  end subroutine

  subroutine mat_eq_zero_4D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:,:)
    integer,              intent(in)    :: z
    integer l
    do l=1, size(A,4); A(:,:,:,l)=z; end do
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


#ifdef UNIT_TEST

subroutine quit(msg)
  character(*), intent(in) :: msg
  print *, msg
  stop
end subroutine



program test

  use matrix_defop_ng  !also contains mat_print
  use matrix_genop_ng, genop_mat_print => mat_print
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
  type(matrix) C, S, D, A, B

  ! manually initialize orbital matrix C
  call mat_nullify(C)
  C%nrow  = size(cmo_t_elms,2)
  C%ncol  = size(cmo_t_elms,1)
  C%closed_shell = .true.
  C%magic_tag = mat_magic_setup !required before alloc
  call mat_alloc(C)
  C%elms = transpose(cmo_t_elms)
  ! print
  !call mat_print(C, label='orbital')

  ! manually initialize overlap matrix S
  call mat_nullify(S)
  S%nrow  = size(ovl_elms,1)
  S%ncol  = size(ovl_elms,2)
  S%magic_tag = mat_magic_setup !required before alloc
  call mat_alloc(S)
  S%elms = ovl_elms
  ! print in plain format
  !call mat_print(S, label='overlap')

  call mat_nullify(D)
  call mat_setup(D, C, C, tb=.true.)
  call mat_alloc(D)
  D%elms = matmul(C%elms,transpose(C%elms))
  !call mat_print(D, label='density')

  ! run tests
  call check_orthonormality_of_C
  call calculate_density_D
  call count_electrons_in_D
  call check_idempotency_of_D
  call test_whether_minus_works

  ! free matrices
  C=0; S=0; D=0; A=0; B=0

contains

  subroutine check_orthonormality_of_C
    ! calculate norm of C^T S C minus identity matrix
    type(matrix) I, CtSCmI
    integer      j
    !------- manually initialize identity matrix
    call mat_nullify(I)
    call mat_setup(I, C, C, ta=.true.)
    call mat_alloc(I)
    I%elms(:,:) = 0
    do j=1,I%nrow; I%elms(j,j)=1; end do
    !call mat_print(I, label='I')
    CtSCmI = trps(C)*S*C - I
    !call mat_print(CtSCmI, label='CtSC-I')
    print *, 'norm(CtSC-I) =', norm(trps(C)*S*C - I) ! ...
    ! free identity matrix
    I=0; CtSCmI=0
  end subroutine


  subroutine calculate_density_D
    ! D = C C^T
    D = C*trps(C) ! ...
    !call mat_print(D, label='density')
  end subroutine


  subroutine count_electrons_in_D
    ! rewrite Tr C^T S C in terms of D
    print *, 'trSD =', dreal(tr(S,D)) ! ...
  end subroutine


  subroutine check_idempotency_of_D
    ! idempotency relation is D S D = D
    print *, 'norm(DSD-D) =', norm(D*S*D-D) ! ...
  end subroutine


  subroutine test_whether_minus_works
    A = -D
    !call mat_print(A, label='minus D')
  end subroutine


end program
#endif
