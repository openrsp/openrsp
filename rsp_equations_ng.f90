! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module rsp_equations

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


module rsp_equations_ng

  use matrix_defop_ng !matrix type and operators
  use rsp_contribs_ng !integrals and integral contractions

  implicit none

  public rsp_fock    !contract perturbed Fock-matrices
  public rsp_scf_eq  !contract perturbed scf equation residuals
  public rsp_dens    !perturbed (or response-) densities and Fock-matrices
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
  subroutine rsp_fock(mol, fld, D, F, S)
    !> mol/basis/decomp/thresh needed by integrals and solver
    type(rsp_cfg),   intent(in) :: mol
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
       call rsp_sol_lookup(mol, fld, dumD, F)
       if (isdef(dumD)) dumD = 0
       if (isdef(F)) return
    end do
    ! S and highest-order contribution to F
    call rsp_oneint(mol, fld, F, S)
    ! add highest-order Coulomb+exchange+Kohn-Sham (unperturbed D)
    call rsp_twoint(mol, fld, D(1:1), F)
    ! add unperturbed Coulomb-exchange+Kohn-Sham contribution
    call rsp_twoint(mol, UNP, D, F)
    ! additional contributions for each basis-dependent subset of fld
    if (size(fld) > 1) &
       call quit('rsp_fock_ovl, not implemented order size(fld) > 1')
  end subroutine



  !> Contracts the derivative of the SCF equation FDS-SDF and
  !> idempotency relation DSD-D, with respect to 'fld', with terms
  !> above order 'keep' omitted. up_to_order=size(fld)-1 gives
  !> (minus) equation right-hand-sides. up_to_order=size(fld) gives
  !> (minus) equation residuals
  subroutine rsp_pert_eqs(mol, fld, S, D, F, DSD, FDS, up_to_order)
    !> mol/basis/thresh's etc. needed by integrals and solver
    type(rsp_cfg),   intent(in) :: mol
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
       DSD = mol%zeromat
       FDS = mol%zeromat
    else
       call quit('rsp_pert_eqs: not implemented order size(fld) > 1')
    end if
  end subroutine



  !> Solve fld-perturbed response equation for perturbed density D.
  !> If already stored in the cache (rsp_cached_sols), the solution is fetched
  !> from there. Otherwise, the RHSs are computed before calling the solver
  subroutine rsp_dens(mol, fld, D)
    !> mol/basis/decomp/thresh needed by integrals and solver
    type(rsp_cfg),   intent(in) :: mol
    !> fields, frequencies and components. NB: all fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> output density matrix perturbed by fields fld
    type(matrix), intent(inout) :: D
    !-------------------------------
    type(matrix) :: S, F, FDS, DSD
    type(rsp_sol), pointer :: sol
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
    call rsp_sol_lookup(mol, fld, D)
    if (isdef(D)) return
    ! calculate all lower-order contributions to the
    ! minus-right-hand-sides DSD and FDS
    call rsp_ctr_eqs(mol, fld, DSD, FDS, size(fld)-1)
    ! calculate all lower-order contributions to the p-perturbed Fock
    ! matrices, as well as the p-perturbed overlap matrices
    call rsp_fock(mol, fld, 
    call pert_fock(mol, p, dimp, (/D,Dp(:pd)/), Fp(:pd), Sp, &
                   comp=ccomp, freq=ffreq)

    ! top DSDp off with DSpD and FDSp with (F-w/2S)DSp-SpD(F+w/2S)
    do i = 1, pd
       !ajt fixme This will not result in completely (anti-)symmetric
       !          DSDp and FDSp. Use +=sym*dag(..) to achieve this
       DSDp(i) = DSDp(i) + D(1)*Sp(i)*D(1)
       FDSp(i) = FDSp(i)  + (F(1) - sum(ffreq)/2 * S0)*D(1)*Sp(i) &
               - Sp(i)*D(1)*(F(1) + sum(ffreq)/2 * S0)
    end do
    Sp(:) = 0 !free
    do i = 1, pd
       call solve_scf_eq(mol, S0, D(1), F(1), sym, sum(ffreq), 1, &
                         DSDp(i:i), FDSp(i:i), Dp(i:i), Fp(i:i))

    end do
    print* !ajt Blank line after prints from solve_scf_eq
  end subroutine



  !> search through rsp_cached_sols, looking for a kept solution
  subroutine rsp_sol_lookup(mol, fld, D, F, S)
    type(rsp_cfg),   intent(in) :: mol
    type(rsp_field), intent(in) :: fld
    type(matrix), intent(inout) :: D
    type(matrix), optional, intent(inout) :: F, S
    !-----------------------------------------------
    type(rsp_sol), pointer :: sol
    type(matrix) :: dumS
    integer      :: nanti
    nanti = count(rsp_field_anti(fld%label))
    sol => rsp_cached_sols
    do while (associated(sol))
       if (size(sol%fld) == size(fld) then
          if (      all(sol%fld%label ==  fld%label) &
              .and. all(sol%fld%comp  ==  fld%comp)  &
              .and. all(sol%fld%scomp ==  fld%scomp) &
              .and.(all(sol%fld%freq  ==  fld%freq)  &
               .or. all(sol%fld%freq  == -fld%freq))) then
             if (all(sol%fld%freq == fld%freq)) then
                D = sol%D
                if (present(F) .and. isdef(sol%F)) F = sol%F
             else !if (all(sol%fld%freq  == -fld%freq)) then
                D = (-1d0)**nanti * dag(sol%D)
                if (present(F) .and. isdef(sol%F)) &
                   F = (-1d0)**nanti * dag(sol%F)
             end if
             if (present(S) .and. isdef(sol%S)) S = sol%S
             exit
          end if
       end if
       sol => sol%next
    end do
  end subroutine


  subroutine rsp_sol_insert(mol, fld, D, F, S)
    type(rsp_cfg),   intent(in) :: mol
    type(rsp_field), intent(in) :: fld
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
  subroutine rsp_solve_eqs(mol, S0, D0, F0, sym, freq, neq, DSD, FDS, D, F)
    !       in:            -----------------------------------------
    !       out:                                                      ----
    !> info needed by integrals and solver
    type(prop_molcfg), intent(in) :: mol
    !> unperturbed overlap matrix, density and Fock matrices
    type(matrix),      intent(in) :: S0, D0, F0
    !> FDS anti:-1, symm:+1, general:0
    integer,           intent(in) :: sym
    !> frequency of equation
    complex(8),        intent(in) :: freq
    !> number of equations
    integer,           intent(in) :: neq
    !> DSD-D residual (-RHS)
    type(matrix),   intent(inout) :: DSD(neq)
    !> FDS-SDF residual (-RHS)
    type(matrix),   intent(inout) :: FDS(neq)
    !----------------------------------------------
    !> resulting perturbed density and Fock
    type(matrix),   intent(inout) :: D(neq), F(neq)
    !----------------------------------------------
    type(matrix)   :: S0D0
    character*4    :: UNP(0)
    integer        :: i, j
    ! ajt Since we are called from pert_dens only, this check is not really needed...
    if (sym /= -1 .and. sym /= 0 .and. sym /= +1) &
       call quit("solve_scf_eq error: Argument 'sym' is not -1, 0 or 1")
    ! prepare the right-hand-sides
    DS = D0*S0
    call scf_eq_prep_rhs(mol, D0, DS, sym, neq, FDSp, DSDp, Fp, Dp)
#ifndef PRG_DIRAC
    ! check whether the frequency is resonant (+-excitation energy)
    call scf_eq_find_resonances(freq, iexci)
    ! if resonant, project away the resonant component form FDSp,
    ! while storing the +freq and -freq transition moments in tsmom
    if (allocated(iexci)) then
       allocate(Vexci(size(iexci)))
       allocate(tsmom(2, neq, size(iexci)))
       do i=1, size(iexci)
          call proj_resonant_rhs(S0, DS, sym, freq, neq,       &
                                 solved_eqs(iexci(i))%freq(1), &
                                 solved_eqs(iexci(i))%D(1),    &
                                 FDSp, Fp, Vexci(i), tsmom(:,:,i))
       end do
    end if
#endif /* not PRG_DIRAC */
    ! call solver
#ifdef PRG_DIRAC
    do i = 1, neq
       real_freq = dreal(freq)
       R = (FDSp(i) - dag(FDSp(i)))
       R = 0.5d0*R
       if (dreal(dot(R,R)) > aor_cfg_threshold_rhs) then
          R%tr_sym = 1
          call response_solver(real_freq, R, Xph(1))
          Dp(i) = Dp(i) + Xph(1)
       end if

       R = (FDSp(i) + dag(FDSp(i)))
       R = 0.5d0*R
       FDSp(i) = 0
       if (dreal(dot(R,R)) > aor_cfg_threshold_rhs) then
          R%tr_sym = -1
          call response_solver(real_freq, R, Xph(1))
          Dp(i) = Dp(i) + Xph(1)
       end if

       R = 0
    end do
#elif defined(DALTON_AO_RSP)
    do i=1, neq
       nrm = norm(FDSp(i))
       print *, 'before response solver: norm(RHS) = ', nrm

       if (nrm < 1d-10) then
          print *, '=> skipping this equation'
          Xph(1) = 0d0 * S0
       else
          freq1(1) = dreal(freq)
          call init_mat(Xph(1), Dp(i))
          call rsp_mosolver_exec( FDSp(i:i), freq1(1:1), Xph(1) )
       end if
       ! if (anti-)symmetric Dp (static with anti-/symmetric FDSp,DSDp),
       ! make sure Dp it comes out completely symmetric
       if (sym /= 0 .and. freq==0) then
          Dp(i) = 1/2d0 * Dp(i) + 2d0 * DS*Xph(1)
          Dp(i) = Dp(i) - 1d0*sym * dag(Dp(i))
       else
          Dp(i) = Dp(i) + 2d0 * (DS*Xph(1) - Xph(1)*dag(DS))
       end if
       Xph(1) = 0
       FDSp(i) = 0
    end do
#else
    do i=1, neq
       rhs_norm = norm(FDSp(i))
       print *, 'before response solver: norm(RHS) = ', rhs_norm
       if (rhs_norm < 1d-10) then
          print *, '=> skipping this equation'
          Xph(1) = 0d0 * S0
       else if (FDSp(i)%complex .or. cfg_rsp_complex .or. dimag(freq) /= 0) then
          !on the first run, initialize the solvers
          if (first_complex) then
             call rsp_init(2, 2, 2, 1, 2)
             call rsp_complex_init(1, 1, 1, 1, 1)
             first = .false.
             first_complex = .false.
          end if
          has_imag = FDSp(i)%complex
          if (.not.has_imag) then
             call init_mat(reFDSp(1), FDSp(i), alias='FF')
             call init_mat(FDSp(i), reset=.true.)
             FDSp(i) = (1d0 + tiny(1d0)*(0d0,1d0)) * reFDSp(1)
             reFDSp(1) = 0 !free
          end if
          !configure the frequency and damping factor
          freq1(1) = dreal(freq)
          gamma_saved = cfg_rsp_gamma !save, for restore after solver
          if (cfg_rsp_complex) then
             cfg_rsp_gamma = cfg_rsp_gamma + dimag(freq)
          else
             cfg_rsp_gamma = dimag(freq)
          end if
          ! the solver currently doesn't seem to work with zero gamma
          if (abs(cfg_rsp_gamma) < 1d-6) &
             cfg_rsp_gamma = 2d-6
          ! allocate the solution matrix Xph (Xp/2)
          call init_mat(Xph(1), FDSp(i))
          !create aliases for the real and imaginary parts of FDSp and Xph
          call init_mat(reFDSp(1), FDSp(i), alias='RF')
          call init_mat(imFDSp(1), FDSp(i), alias='IF')
          call init_mat(reXph(1), Xph(1), alias='RF')
          call init_mat(imXph(1), Xph(1), alias='IF')
          ! solve for the real and imaginary part of the FDSp simlultaneously
          call rsp_complex_solver(decomp, F0, D0, S0, 1, reFDSp(1:1), &
                                  freq1, 1, reXph(1:1), imXph(1:1),       &
                                  .true., gdi=imFDSp(1:1))
          !clear re and im aliases of FDSp and Xph
          call init_mat(reFDSp(1), reset=.true.)
          call init_mat(imFDSp(1), reset=.true.)
          call init_mat(reXph(1), reset=.true.)
          call init_mat(imXph(1), reset=.true.)
          !restore gamma to avoid breaking other code
          cfg_rsp_gamma = gamma_saved
       else
          if (first) then
             call rsp_init(1, 1, 1, 1, 1)
             first = .false.
          end if
          freq1(1) = dreal(freq)
          call init_mat(Xph(1), Dp(i))
          ! if non-resonant, no excitation vectors to pass
          if (.not.allocated(iexci)) then
             call rsp_solver(decomp, D0, S0, F0, .true., 1, FDSp(i:i), &
                             freq1(1:1), Xph(1))
          else !if resonant, pass excitation vectors Vexci
             call rsp_solver(decomp, D0, S0, F0, .true., 1, FDSp(i:i), &
                             freq1(1:1), Xph(1), Xproject=Vexci)
          end if
       end if
       ! if (anti-)symmetric Dp (static with anti-/symmetric FDSp,DSDp),
       ! make sure Dp it comes out completely symmetric
       if (sym /= 0 .and. freq==0 .and. .not.cfg_rsp_complex) then
          Dp(i) = 1/2d0 * Dp(i) + 2d0 * DS*Xph(1)
          Dp(i) = Dp(i) - 1d0*sym * dag(Dp(i))
       else
          Dp(i) = Dp(i) + 2d0 * (DS*Xph(1) - Xph(1)*dag(DS))
       end if
       Xph(1) = 0
       FDSp(i) = 0
    end do
#endif /* PRG_DIRAC */
    DS = 0
    ! add last contribution G(Dp) to Fp
    call prop_twoint(mol, UNP, (/D0, Dp/), (/neq/), Fp)
  end subroutine



  !> Private subroutine for preparing right-hand-sides before calling the solver
  !> Calculates particular component of D from DSD, adds contribution
  !> from that to FDS, and projects out the particular component from FDS
  !> (this is for numerical stability). Precalculated S0*D0 is taken as input.
  subroutine scf_eq_prep_rhs(mol, D0, S0D0, sym, neq, FDS, DSD, F, D)
      type(prop_molcfg), intent(in) :: mol
      type(matrix),      intent(in) :: D0, S0D0, F(neq)
      integer,           intent(in) :: sym, neq
      type(matrix),   intent(inout) :: FDS(neq), DSD(neq), D(neq)
      type(rsp_cfg) :: UNP(0)
      integer       :: i
      do i=1, neq
         ! if the minus-right-hand-sides are (anti-)symmetric, one gets away
         ! with half the number of multiplications (though twice as many lines)
         if (sym /= 0) then
            D(i) = 1/2d0 * DSD(i) - DSD(i)*S0D0
            D(i) = D(i) - sym*1d0 * dag(D(i)) !symmetrize or anti-symmetrize
            FDS(i) = -S0D0*FDS(i)
            FDS(i) = FDS(i) - sym*1d0 * dag(FDS(i)) + F(i)
            ! ajt Ensure that the perturbed Fock is precisely (anti-)symmetric
            !ajt wrong: F(i) = 1/2d0 * F(i) - 1/2d0*sym * dag(F(i)) !(anti-)symmetrize
         else
            Dp(i) = DSDp(i) - DSDp(i)*dag(DS) - DS*DSDp(i)
            FDSp(i) = FDSp(i)*DS - dag(DS)*FDSp(i) + Fp(i)
         end if
         DSDp(i) = 0
      end do
      ! contract Coulomb-exchange+Kohn-Sham with particular components D(:)
      call rsp_twoint(mol, UNP, (/D0,D/), (/neq/), FDS)
      ! complete right-hand-sides
      do i=1, neq
         if (sym /= 0) then
            FDS(i) = -S0D0*FDS(i)
            FDS(i) = FDS(i) + 1d0*sym * trps(FDS(i))
         else
            FDS(i) = FDS(i)*trps(S0D0) - S0D0*FDS(i)
         end if
      end do
   end subroutine


end module
