! Copyright 2012      Magnus Ringholm
!           2012      Dan Jonsson
!           2009-2011 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_equations

! ajt fixme Add detailed description of contents here
!
! ajt jan10 Added place to store (cache) solutions of response
!           equations, solved_eqs(100), and avoid re-solving the same
!           equations later. Solved_eqs can also hold excitation densities,
!           used when solving for higher-order transition densities,
!           perturbed transition densities, as well as projected perturbed
!           densities.
!
! ajt feb10 Added solving for the non-resonant component of resonant
!           equations. Separated some utility routines from solve_scf_eq.
!
!> Delivers perturbed density matrices (solutions of response equations).
!> Property integrals are obtained through module prop_contribs.
!> Right-hand-sides are contracted here, and solutions are obtained
!> by invoking one of several response solvers, configured at compile
!> time by #define statements.

!  version history (response code, not just this file):
!  2010-04-21 identical (shared) files: dalton/branches/linsca-ng revision  7878
!                                       dirac/trunk               revision 10659


module rsp_equations

  use matrix_defop !matrix type and operators
  use rsp_contribs !integrals and integral contractions

  implicit none

  public rsp_fock     !contract perturbed Fock-matrices
  public rsp_eval_eqs !contract perturbed scf equation residuals
  public rsp_dens     !perturbed (or response-) densities and Fock-matrices
  public rsp_equations_debug  !turn debug printing on or off

  ! ajt The following three are public for the moment, but are intended
  !     be accessed through subroutines eventually
  public rsp_sol, rsp_cached_sols

  !> turn on or off debugging in this module
  logical :: rsp_equations_debug = .false.

  !> Maximum order of equation solution which can be cached.
  !> Defines the size of the fields in type cached_sol.
  !> ajt fixme Should be removed once pert_dens can do arbitrary-order
  !> equations.
  integer, parameter :: max_order = 3


  !> Type for saving (caching) solutions of response equations,
  !> and avoid re-solving the same equations later in the program.
  type rsp_sol
     !> fields, components and frequencies for the eq
     type(rsp_field), pointer :: fld(:)
     !> perturbed overlap, density and Fock for this equation
     type(matrix)             :: S, D, F
     !> pointer to next in linked list
     type(rsp_sol),   pointer :: next
  end type


  !> To keep collection of saved response equation solutions.
  !> ajt fixme Currently hard-coded size 100
  type(rsp_sol), pointer :: rsp_cached_sols

  ! ajt For some reason, if I put this 'private' above type rsp_sol,
  !     Doxygen does not document the fields inside the type, dispite
  private !it being declared public at the top. Maybe a bug in Doxygen.

contains


  !> contract *one* full perturbed Fock matrix, and additionally calculate
  !> the correcponding perturbed overlap matrix, for use in right-hand-sides
  subroutine rsp_fock(fld, D, F, S)
    !> field lables, frequencies and components. All fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> density matrix expansion for fld(:)
    type(matrix),    intent(in) :: D(2**size(fld))
    !> output perturbed Fock matrix
    type(matrix), intent(inout) :: F
    !> optional output perturbed overlap matrix
    type(matrix), intent(inout), optional :: S
    !-----------------------------------------
    type(rsp_field) :: UNP(0) !unperturbed
    type(matrix)    :: dumD
    ! insist that each fld has only one component (thus also F)
    if (any(fld%ncomp /= 1)) &
       call quit('error: rsp_fock expected fld%ncomp = 1')
    ! look through rsp_cached_sols for cached F
    if (.not.iszero(D(size(D))) .and. .not.present(S)) then
       if (isdef(F)) F = 0
       call rsp_sol_lookup(fld, dumD, F)
       if (isdef(dumD)) dumD = 0
       if (isdef(F)) return
    end if
!     ! S and highest-order contribution to F
!     call rsp_oneint(mol, fld, F, S)
!     ! add highest-order Coulomb+exchange+Kohn-Sham (unperturbed D)
!     call rsp_twoint(mol, fld, D(1:1), F)
!     ! add unperturbed Coulomb-exchange+Kohn-Sham contribution
!     call rsp_twoint(mol, UNP, D, F)
    ! additional contributions for each basis-dependent subset of fld
    if (size(fld) > 1) &
       call quit('rsp_fock_ovl, not implemented order size(fld) > 1')
  end subroutine



  !> Contracts the derivative of the SCF equation FDS-SDF and
  !> idempotency relation DSD-D, with respect to 'fld', with terms
  !> above order 'keep' omitted. up_to_order=size(fld)-1 gives
  !> (minus) equation right-hand-sides. up_to_order=size(fld) gives
  !> (minus) equation residuals
  subroutine rsp_eval_eqs(fld, S, D, F, DSD, FDS, up_to_order)
    !> fields, frequencies and components. NB: all fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> perturbation series for overlap matrix
    type(matrix),    intent(in) :: S(2**(size(fld)))
    !> perturbation series for density matrix
    type(matrix),    intent(in) :: D(2**(size(fld)))
    !> perturbation series for Fock matrix
    type(matrix),    intent(in) :: F(2**(size(fld)))
    !> resulting perturbed idempotency condition
    type(matrix), intent(inout) :: DSD
    !> resulting perturbed SCF equation (or residual)
    type(matrix), intent(inout) :: FDS
    !> drop terms containing any factor above this order
    integer,         intent(in) :: up_to_order
    !-----------------------------------------
    ! insist that each fld has only one component (thus also F)
    if (any(fld%ncomp /= 1)) &
       call quit('error: rsp_fock expected fld%ncomp = 1')
    if (size(fld) == 1 .and. up_to_order == 0) then
       DSD = 0*S(1)
       FDS = 0*DSD
    else
       call quit('rsp_pert_eqs: not implemented order size(fld) > 1')
    end if
  end subroutine



  !> Solve fld-perturbed response equation for perturbed density D.
  !> If already stored in the cache (rsp_cached_sols), the solution is fetched
  !> from there. Otherwise, the RHSs are computed before calling the solver
  subroutine rsp_dens(fld, D)
    !> fields, frequencies and components. NB: all fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> output density matrix perturbed by fields fld
    type(matrix), intent(inout) :: D
    !-------------------------------
    type(matrix), dimension(2**size(fld)) :: S, DD, F
    type(matrix) :: FDS, DSD
    integer      :: nanti, sym
    ! insist that each fld has only one component (thus also F)
    if (any(fld%ncomp /= 1)) &
       call quit('error: rsp_dens expected fld%ncomp = 1')
    ! determine the symmetry of the right-hand-sides
    ! ajt FIXME Using sym, rsp_pert_eqs can be optimized
    nanti = count(rsp_field_anti(fld%label))
    sym = merge(merge(1, -1, mod(nanti,2) == 1), &
                0, all(fld%freq==0) .or. (size(fld)==1 .and. &
                .not.any(rsp_field_bas(fld%label))))
    ! search through solved_eqs for an already calculated solution
    if (isdef(D)) D = 0
    call rsp_sol_lookup(fld, D)
    if (isdef(D)) return
    ! calculate all lower-order contributions to the
    ! minus-right-hand-sides DSD and FDS
    !ajt FIXME Disabled for now
    ! call rsp_eval_eqs(mol, fld, S, D, F, DSD, FDS, size(fld)-1)
    DSD = 0*S(1)
    call mat_ensure_alloc(DSD)
    FDS = 0*S(1)
    call mat_ensure_alloc(FDS)
    ! calculate all lower-order contributions to the p-perturbed Fock
    ! matrices, as well as the p-perturbed overlap matrices
!     call rsp_fock(mol, fld,
!     call pert_fock(mol, p, dimp, (/D,Dp(:pd)/), Fp(:pd), Sp, &
!                    comp=ccomp, freq=ffreq)
    ! top DSDp off with DSpD and FDSp with (F-w/2S)DSp-SpD(F+w/2S)
    DSD = DSD + DD(1) * S(size(S)) * DD(1)
    FDS = FDS + (F(1) - sum(fld%freq)/2  * S(1)) * DD(1) * S(size(S)) &
              - S(size(S)) * DD(1) * (F(1) + sum(fld%freq)/2 * S(1))
    S = 0 !free
    call solve_scf_eq(S(1), DD(1), F(1), sym, sum(fld%freq), &
                      1, DSD, FDS, D, F(size(F)))
    print* !ajt Blank line after prints from solve_scf_eq
  end subroutine



  !> search through rsp_cached_sols, looking for a kept solution
  subroutine rsp_sol_lookup(fld, D, F, S)
    type(rsp_field), intent(in) :: fld(:)
    type(matrix), intent(inout) :: D
    type(matrix), optional, intent(inout) :: F, S
    !--------------------------------------------
    type(rsp_sol), pointer :: sol
    type(matrix) :: dumS
    integer      :: nanti
    nanti = count(rsp_field_anti(fld%label))
    sol => rsp_cached_sols
    do while (associated(sol))
       if (size(sol%fld) == size(fld)) then
          if (      all(sol%fld%label ==  fld%label) &
              .and. all(sol%fld%comp  ==  fld%comp)  &
              .and. all(sol%fld%ncomp ==  fld%ncomp) &
              .and.(all(sol%fld%freq  ==  fld%freq)  &
               .or. all(sol%fld%freq  == -fld%freq))) then
             if (all(sol%fld%freq == fld%freq)) then
                D = sol%D
                if (present(F) .and. isdef(sol%F)) F = sol%F
             else !if (all(sol%fld%freq  == -fld%freq)) then
                D = (-1d0)**nanti * trans(sol%D)
                if (present(F) .and. isdef(sol%F)) &
                   F = (-1d0)**nanti * trans(sol%F)
             end if
             if (present(S) .and. isdef(sol%S)) S = sol%S
             exit
          end if
       end if
       sol => sol%next
    end do
  end subroutine


  subroutine rsp_sol_insert(fld, D, F, S)
    type(rsp_field),        intent(in) :: fld
    type(matrix),           intent(in) :: D
    type(matrix), optional, intent(in) :: F, S
  end subroutine


  !> Solve an scf response equation by invoking the configured response
  !> solver on the residuals (-RHS'es) in DSDp and FDSp, obtaining
  !> solutions Dp (perturbed density) and Fp (perturbed Fock)
  !> ajt Should at some point also return the final residuals, at least
  !>     FDSp, for use in Sellers-type higher-precision formulas
  !> ajt aug09 Added argument sym, which specifies the symmetry of the minus-
  !>           right-hand-side FDSp: anti-symmetric=-1, symmetric=+1, general=0
  !> ajt sep09 Added argument neq (number of equations), so several can be
  !>           solved in one go
  !> ajt feb10 For linsca, added comparison with excitation energies,
  !>           and projection
  subroutine rsp_solve_eqs(S0, D0, F0, sym, freq, neq, DSD, FDS, D, F)
    !       in:            -----------------------------------------
    !       out:                                                      ----
    !> unperturbed overlap matrix, density and Fock matrices
    type(matrix),    intent(in) :: S0, D0, F0
    !> FDS anti:-1, symm:+1, general:0
    integer,         intent(in) :: sym
    !> frequency of equation
    complex(8),      intent(in) :: freq
    !> number of equations
    integer,         intent(in) :: neq
    !> DSD-D residual (-RHS)
    type(matrix), intent(inout) :: DSD(neq)
    !> FDS-SDF residual (-RHS)
    type(matrix), intent(inout) :: FDS(neq)
    !----------------------------------------------
    !> resulting perturbed density and Fock
    type(matrix), intent(inout) :: D(neq), F(neq)
    !----------------------------------------------
    type(matrix) :: SD0
    character(4) :: UNP(0)
    real(8)      :: freq1(1), norm_rhs
    integer      :: i, j
    ! ajt Since we are called from pert_dens only, this check is not really needed...
    if (sym /= -1 .and. sym /= 0 .and. sym /= +1) &
       call quit("solve_scf_eq error: Argument 'sym' is not -1, 0 or 1")
    ! prepare the right-hand-sides
    SD0 = S0*D0
    call scf_eq_prep_rhs(D0, SD0, sym, neq, FDS, DSD, F, D)
    ! call solver
    do i=1, neq
       norm_rhs = norm(FDS(i))
       print *, 'before response solver: norm(RHS) = ', norm_rhs

       if (norm_rhs < 1d-10) then
          print *, '=> skipping this equation'
          ! Xph(1) = 0d0 * S0
       else
          freq1(1) = dreal(freq)
          ! call init_mat(Xph(1), Dp(i))
          call rsp_mosolver_exec(FDS(i:i), freq1(1:1), D(i:i)) !Xph(1))
       end if
       ! if (anti-)symmetric Dp (static with anti-/symmetric FDSp,DSDp),
       ! make sure Dp it comes out completely symmetric
       !if (sym /= 0 .and. freq==0) then
       !   Dp(i) = 1/2d0 * Dp(i) + 2d0 * DS*Xph(1)
       !   Dp(i) = Dp(i) - 1d0*sym * dag(Dp(i))
       !else
       !   Dp(i) = Dp(i) + 2d0 * (DS*Xph(1) - Xph(1)*dag(DS))
       !end if
       !Xph(1) = 0
       FDS(i) = 0
    end do
    SD0 = 0
    ! add last contribution G(Dp) to Fp
    ! call rsp_twoint(mol, UNP, (/D0, Dp/), (/neq/), Fp)
  end subroutine



  !> Private subroutine for preparing right-hand-sides before calling the solver
  !> Calculates particular component of D from DSD, adds contribution
  !> from that to FDS, and projects out the particular component from FDS
  !> (this is for numerical stability). Precalculated S0*D0 is taken as input.
  subroutine scf_eq_prep_rhs(D0, SD0, sym, neq, FDS, DSD, F, D)
      type(matrix),  intent(in)    :: D0, SD0, F(neq)
      integer,       intent(in)    :: sym, neq
      type(matrix),  intent(inout) :: FDS(neq), DSD(neq), D(neq)
      integer       :: i
      do i = 1,neq
         ! if the minus-right-hand-sides are (anti-)symmetric, one gets away
         ! with half the number of multiplications (though twice as many lines)
         if (sym /= 0) then
            D(i) = 1/2d0 * DSD(i) - DSD(i)*SD0
            D(i) = D(i) - sym*1d0 * trans(D(i)) !symmetrize or anti-symmetrize
            FDS(i) = -SD0*FDS(i)
            FDS(i) = FDS(i) - sym*1d0 * trans(FDS(i)) + F(i)
            ! ajt Ensure that the perturbed Fock is precisely (anti-)symmetric
            !ajt wrong: F(i) = 1/2d0 * F(i) - 1/2d0*sym * dag(F(i)) !(anti-)symmetrize
         else
            D(i) = DSD(i) - DSD(i)*SD0 - trans(SD0)*DSD(i)
            FDS(i) = FDS(i)*trans(SD0) - SD0*FDS(i) + F(i)
         end if
         DSD(i) = 0
      end do
      ! contract Coulomb-exchange+Kohn-Sham with particular components D(:)
      ! call rsp_twoint(mol, UNP, (/D0,D/), (/neq/), FDS)
      ! complete right-hand-sides
      do i = 1,neq
         if (sym /= 0) then
            FDS(i) = -SD0*FDS(i)
            FDS(i) = FDS(i) + 1d0*sym * trans(FDS(i))
         else
            FDS(i) = FDS(i)*trans(SD0) - SD0*FDS(i)
         end if
      end do
   end subroutine


end module
