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


module rsp_equations

  use matrix_defop !matrix type and operators
  use rsp_contribs !integrals and integral contractions

  implicit none

  public rsp_fock    !contract perturbed Fock-matrices
  public rsp_scf_eq  !contract perturbed scf equation residuals
  public rsp_dens    !perturbed (or response-) densities and Fock-matrices
  public rsp_equations_debug  !turn debug printing on or off

  ! ajt The following three are public for the moment, but are intended
  !     be accessed through subroutines eventually
  public rsp_sol, rsp_kept_sols

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
  type(rsp_sol), pointer :: rsp_kept_sols

  ! ajt For some reason, if I put this 'private' above type rsp_sol,
  !     Doxygen does not document the fields inside the type, dispite
  private !it being declared public at the top. Maybe a bug in Doxygen.

contains


  !> contract *one* full perturbed Fock matrix, and additionally calculate
  !> the correcponding perturbed overlap matrix, for use in right-hand-sides
  subroutine rsp_fock_ovl(mol, fld, D, F, S)
    !> mol/basis/decomp/thresh needed by integrals and solver
    type(rsp_cfg),   intent(in) :: mol
    !> field lables, frequencies and components. All fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> density matrix expansion for fld(:) unperturbed and perturbed
    !> density matrices. size(D) = 2**size(fld)
    type(matrix),    intent(in) :: D(:)
    !> output perturbed Fock and overlap matrices
    type(matrix), intent(inout) :: F, S
    !----------------------------------
    type(rsp_field) :: UNP(0) !unperturbed
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
  subroutine rsp_pert_eqs(mol, fld, FDS, DSD, up_to_order)
    !> mol/basis/thresh's etc. needed by integrals and solver
    type(rsp_cfg),   intent(in) :: mol
    !> fields, frequencies and components. NB: all fld%ncomp must be 1
    type(rsp_field), intent(in) :: fld(:)
    !> resulting perturbed SCF equation (or residual)
    type(matrix), intent(inout) :: FDS
    !> resulting perturbed idempotency condition
    type(matrix), intent(inout) :: DSD
    !> drop terms containing any factor above this order
    integer,         intent(in) :: up_to_order
    !-----------------------------------------
    if (size(fld) == 1 .and. up_to_order == 0) then
       FDS = mol%zeromat
       DSD = mol%zeromat
    else
       call quit('rsp_scf_eq: not implemented order size(fld) > 1')
    end if
  end subroutine



  !> Solve fld-perturbed response equation for perturbed density D.
  !> If already stored in the cache (rsp_kept_sols), the solution is fetched
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
    ! zero highest-order D and F for RHS calculation
    do i = 1, pd
       Dp(i) = 0d0*D(1)
       Fp(i) = 0d0*F(1)
    end do
    ! calculate all lower-order contributions to the
    ! minus-right-hand-sides DSDp and FDSp
    call pert_scf_eq(mol, S0, p, dimp, (/D,Dp(:pd)/), (/F,Fp(:pd)/), FDSp, DSDp, &
                     comp=ccomp, freq=ffreq)
    ! calculate all lower-order contributions to the p-perturbed Fock
    ! matrices, as well as the p-perturbed overlap matrices
    call pert_fock(mol, p, dimp, (/D,Dp(:pd)/), Fp(:pd), Sp, &
                   comp=ccomp, freq=ffreq)

    ! top DSDp off with DSpD and FDSp with (F-w/2S)DSp-SpD(F+w/2S)
    do i = 1, pd
#ifdef PRG_DIRAC
       ! radovan: this is necessary otherwise matrices inherit wrong irep from Sp
       !          which breaks point group symmetry in DIRAC
       Sp(i)%irep = DSDp(i)%irep
#endif
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



  !> search through rsp_kept_sols, looking for a kept solution
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
    sol => rsp_kept_sols
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
  subroutine solve_scf_eq(mol, S0, D0, F0, sym, freq, neq, DSDp, FDSp, Dp, Fp)
    !       in:           --------------------------------------------
    !       out:                                                       ------
#ifdef PRG_DIRAC
    use dirac_interface
    use aor_cfg
#endif
#ifdef DALTON_AO_RSP
    use dalton_ifc
#endif
    !> mol/basis/decomp/thresh needed by integrals and solver
    type(prop_molcfg), intent(in) :: mol
    !> unperturbed overlap matrix
    type(matrix),      intent(in) :: S0
    !> unperturbed density matrix
    type(matrix),      intent(in) :: D0
    !> unperturbed Fock matrix
    type(matrix),      intent(in) :: F0
    !> FDSp anti:-1, symm:+1, general:0
    integer,           intent(in) :: sym
    !> frequency of equation
    complex(8),        intent(in) :: freq
    !> number of equations
    integer,           intent(in) :: neq
    !> DSD-D residual (-RHS)
    type(matrix),   intent(inout) :: DSDp(neq)
    !> FDS-SDF residual (-RHS)
    type(matrix),   intent(inout) :: FDSp(neq)
    !-----------------------------------------
    !> resulting perturbed density
    type(matrix),   intent(inout) :: Dp(neq)
    !> resulting perturbed Fock
    type(matrix),   intent(inout) :: Fp(neq)
    !-----------------------------------------
    type(matrix)   :: DS
    character*4    :: UNP(0)
    integer        :: i, j
    real(8)        :: rhs_norm
    type(matrix)   :: Xph(1)
#ifdef PRG_DIRAC
    ! R for sym/anti of FDSp, Dhom for homogeneous part of Dp
    type(matrix)   :: R, Dhom
    real(8)        :: real_freq
#else
    logical, save  :: first = .true.
    real(8)        :: freq1(1) !duplicate but array(1)
    ! for complex equations
    logical, save  :: first_complex = .true.
    real(8)        :: gamma_saved
    logical        :: has_imag
    type(matrix)   :: reFDSp(1), imFDSp(1), reXph(1), imXph(1)
    ! for resonant equations (real freq coincides with +-an excitation energy).
    ! Excitation energies and densities are stored in solved_eqs.
    ! Excitation and de-excitation should be projected out of the
    ! right-hand-sides, and de-/excitation vectors should be passed to the solver
    integer,      allocatable :: iexci(:)     !indices in solved_eqs
    type(matrix), allocatable :: Vexci(:)     !resonant excitation vector(s)
    complex(8),   allocatable :: tsmom(:,:,:) !transition moments
    ! points to mol%decomp, which in intent(inout) in the solvers
    type(decompItem), pointer :: decomp
    decomp => mol%decomp
#endif /* PRG_DIRAC */
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
#ifndef PRG_DIRAC
    ! if this is a projected resonant equation, add contribution from
    ! the opposite (de-)excitation (then deallocate)
    if (allocated(iexci)) then
       Vexci(:) = 0 !free excitation vectors
       deallocate(Vexci)
       do i=1, size(iexci)
          call resonant_add_nonres_comp(freq, neq,         &
                             solved_eqs(iexci(i))%freq(1), &
                             solved_eqs(iexci(i))%D(1),    &
                             tsmom(:,:,i), Dp)
       end do
       deallocate(iexci)
       deallocate(tsmom)
    end if
#endif /* not PRG_DIRAC */
    DS = 0
    ! add last contribution G(Dp) to Fp
    call prop_twoint(mol, UNP, (/D0, Dp/), (/neq/), Fp)
  end subroutine



  !> Private subroutine for preparing right-hand-sides before calling the solver
  !> Calculates particular component of Dp from DSDp, adds contribution
  !> from that to FDSp, and projects out the particular component from FDSp
  !> (this is for numerical stability). Precalculated DS is taken as input.
  subroutine scf_eq_prep_rhs(mol, D0, DS, sym, neq, FDSp, DSDp, Fp, Dp)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(prop_molcfg), intent(in) :: mol
      type(matrix),      intent(in) :: D0, DS, Fp(neq)
      integer,           intent(in) :: sym, neq
      type(matrix),   intent(inout) :: FDSp(neq), DSDp(neq), Dp(neq)
      character*4 :: UNP(0)
      integer     :: i
      do i=1, neq
         ! if the minus-right-hand-sides are (anti-)symmetric, one gets away
         ! with half the number of multiplications (though twice as many lines)
         if (sym /= 0) then
            Dp(i) = 1/2d0 * DSDp(i) - DS*DSDp(i)
            Dp(i) = Dp(i) - sym*1d0 * dag(Dp(i)) !symmetrize or anti-symmetrize
            FDSp(i) = FDSp(i)*DS
            FDSp(i) = FDSp(i) - sym*1d0 * dag(FDSp(i)) + Fp(i)
            ! ajt Ensure that the perturbed Fock is precisely (anti-)symmetric
            !ajt wrong: Fp(i) = 1/2d0 * Fp(i) - 1/2d0*sym * dag(Fp(i)) !(anti-)symmetrize
         else
            Dp(i) = DSDp(i) - DSDp(i)*dag(DS) - DS*DSDp(i)
            FDSp(i) = FDSp(i)*DS - dag(DS)*FDSp(i) + Fp(i)
         end if
         DSDp(i) = 0
      end do
      ! contract Coulomb-exchange+Kohn-Sham with particular components Dp(;)
      call prop_twoint(mol, UNP, (/D0, Dp/), (/neq/), FDSp)
      ! complete right-hand-sides
      do i=1, neq
         if (sym /= 0) then
            FDSp(i) = FDSp(i)*DS
            FDSp(i) = FDSp(i) + 1d0*sym * dag(FDSp(i))
         else
            FDSp(i) = FDSp(i)*DS - dag(DS)*FDSp(i)
         end if
      end do
   end subroutine



#ifndef PRG_DIRAC
  !> Look through solved_eqs for whether freq is an excitation energy,
  !> and return then index of that entry. If not, return -1
  subroutine scf_eq_find_resonances(freq, iexci)
    complex(8),           intent(in)  :: freq
    integer, allocatable, intent(out) :: iexci(:)
    integer :: n, i, j
    do j=1,2 !j=1 counts, then allocates iexci
       n = 0 !j=2 records indices in iexci
       do i=1, num_solved_eqs
          ! ajt FIXME Should get this 1d-5 tolerance from some thresholds
          if (solved_eqs(i)%order   == 1      .and. &
              solved_eqs(i)%pert(1) == 'EXCI' .and. &
              (abs(solved_eqs(i)%freq(1) - freq) < 1d-5 .or. &
               abs(solved_eqs(i)%freq(1) + freq) < 1d-5)) then
             n = n+1
             if (j==2) iexci(n) = i
          end if
       end do
       if (j==1 .and. n==0) exit !no resonant excitations
       if (j==1 .and. n/=0) allocate(iexci(n))
    end do
  end subroutine
#endif /* not PRG_DIRAC */



#ifndef PRG_DIRAC
  !> Project the excitation densities out of the -right-hand-sides FDSp,
  !> and return the corresponding excitation vector Vexci and
  !> transition moments tsmom. tsmom(1,:) are de-excitation moments (-exci),
  !> while tsmom(2,:) are excitation moments (+exci).
  !> ajt may10 Added forgotten projection of Fp
  subroutine proj_resonant_rhs(S0, DS, sym, freq, neq, &
                               exci, Dx, FDSp, Fp, Vx, tsmom)
    type(matrix), intent(in)    :: S0, DS
    integer,      intent(in)    :: sym, neq
    complex(8),   intent(in)    :: freq, exci
    type(matrix), intent(in)    :: Dx
    type(matrix), intent(inout) :: FDSp(neq), Fp(neq), Vx
    complex(8),   intent(out)   :: tsmom(2,neq)
    type(matrix) :: SDxS, SVxS
    integer      :: i, j, k
    Vx = Dx*dag(DS) - DS*Dx !for dot/tr with FDSp
    SDxS = S0*Dx*S0         !for removal from FDSp
    SVxS = S0*Vx*S0         !for removal from Fp
    ! loop over -right-hand-sides (equations)
    do i=1, neq
       ! excitation transition moment (frequency +exci)
       tsmom(2,i) = dot(Vx, FDSp(i))
       ! de-excitation transition moment (frequency -exci)
       if (sym/=0) tsmom(1,i) = sym * tsmom(2,i)
       if (sym==0) tsmom(1,i) = tr(Vx, FDSp(i))
       ! remove both resonant and opposite-resonant component
       ! from FDSp. This preserves symmetry. The solution's
       ! opposite-resonant component will be added after solving.
       ! Not activated in LINSCA, as the solver takes care of this
       FDSp(i) = FDSp(i) - tsmom(1,i) * dag(SDxS) &
                         - tsmom(2,i) * SDxS
       ! remove resonant component from Fp, so Fp/Dp will
       ! solve the response equation FpDS+FDpS...-wSDpS=0
       if (abs(exci+freq) < abs(exci-freq)) then
          Fp(i) = Fp(i) + tsmom(1,i) * dag(SVxS)
       else
          Fp(i) = Fp(i) - tsmom(2,i) * SVxS
       end if
    end do
    SDxS = 0
    SVxS = 0
  end subroutine
#endif /* not PRG_DIRAC */



#ifndef PRG_DIRAC
  !> For an equation which was resonant (excitation: freq=exci,
  !> de-excitation: freq=-exci), and is being solved for the
  !> non-resonant part, add the contribution to the solution Dp
  !> from the opposite (de-)excitation density (Dx or Dx^T).
  !> The factor before Dx or Dx^T is given by the transition moment
  !> tsmom(2 or 1,:), divided by twice the frequency
  subroutine resonant_add_nonres_comp(freq, neq, exci, Dx, tsmom, Dp)
    complex(8),   intent(in)    :: freq, exci
    integer,      intent(in)    :: neq
    type(matrix), intent(in)    :: Dx
    complex(8),   intent(in)    :: tsmom(2,neq)
    type(matrix), intent(inout) :: Dp(neq)
    integer :: i
    do i=1, neq
       ! if excitation resonant, add de-excitation contrib a*Dx^T
       if (abs(exci-freq) <= abs(exci+freq)) then
          ! in LINSCA, the solver takes care of this
          Dp(i) = Dp(i) + (tsmom(1,i) / (freq+exci)) * dag(Dx)
       else !if de-excitation resonant, add excitation contrib a*Dx
          Dp(i) = Dp(i) + (tsmom(2,i) / (freq-exci)) * Dx
       end if
    end do
  end subroutine
#endif /* not PRG_DIRAC */


end module
