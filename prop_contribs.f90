! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module prop_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices. 
module prop_contribs

   use matrix_defop
#ifdef PRG_DIRAC
   use dirac_interface
   use character_processing
   use dft_main
#endif
#ifdef DALTON_AO_RSP
   use dalton_ifc,    &
     di_GET_GbDs => di_get_gmat
   use xmatrix, only: &
     mat_to_full => xmatrix_get_full
   use xcharacter, only: &
     prefix_zeros => xchar_from_num
#endif
#ifdef VAR_LINSCA
   use dal_interface, only:              &
           di_read_operator_int,         &
           di_GET_GbDs,                  &
           di_get_sigma_xc_cont,         &
           di_get_T_xc_cont,             &
           di_get_MagDeriv_FxD_DFT,      &
           di_get_MagDeriv_GxD_DFT,      &
           di_get_geomDeriv_molgrad_DFT, &
           di_get_geomDeriv_FxD_DFT,     &
           di_get_geomDeriv_GxD_DFT,     &
           di_select_wrk,                &
           di_deselect_wrk

   use matrix_operations, only: &
           mat_to_full,         &
           mat_set_from_full
#endif

   implicit none
   public prop_oneave
   public prop_twoave
   public prop_oneint
   public prop_twoint
   public pert_basdep
   public pert_shape
   public pert_antisym
   public prop_auxlab

   private

   type ptb !private: 'perturbation'
      character*4  :: code  !four-letter abbreviation
      character*64 :: name  !long name
      integer      :: ncomp !number of components (when known, 0 otherwise)
      logical      :: anti  !anti-symmetric (1,3,5th ord.) perturbed integrals
      logical      :: bas   !basis dependent (sa. GEO and MAG)
      logical      :: lin   !one-electron operator linear in field strength (EL)
      logical      :: quad  !one-electron operator quadratic in field strength (MAGO)
   end type

   ! to compactify the table
   logical, parameter:: T = .true., F = .false.

   ! ajt note: EXCI is a ZERO (no) perturbation, and is introduced to
   !           allow the same code to contract response functions and
   !           "generalized transition moments".
   ! ajt note: AUX0..AUX9 are 10 configurable basis-independent 1-electron
   !           perturbations, configured by setting the corresponding
   !           HERMIT integral label in prop_auxlab(0:9).
   type(ptb) :: pert_table(12) = &                        !nc an ba ln qu
       (/ptb('EXCI', 'Generalized "excitation" field'     , 1, F, F, T, T), &
         ptb('AUX*', 'Auxiliary integrals on file'        , 1, F, F, T, F), &
         ptb('EL'  , 'Electric field'                     , 3, F, F, T, F), &
         ptb('VEL' , 'Velocity'                           , 3, T, F, T, F), &
         ptb('MAGO', 'Magnetic field w/o. London orbitals', 3, T, F, F, T), &
         ptb('MAG' , 'Magnetic field with London orbitals', 3, T, T, F, F), &
         ptb('ELGR', 'Electric field gradient'            , 6, F, F, T, F), &
         ptb('VIBM', 'Dispacement along vibrational modes', 0, F, T, F, F), &
         ptb('GEO' , 'Nuclear coordinates'                , 0, F, T, F, F), & !0=uninitialized
         ptb('NUCM', 'Nuclear magnetic moment'            , 0, F, T, F, T), &
         ptb('AOCC', 'AO contraction coefficients'        , 0, F, T, F, F), &
         ptb('AOEX', 'AO exponents'                       , 0, F, T, F, F)/)

   character(8) :: prop_auxlab(0:9)

   integer :: natom=0 !unitialized

   character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

#ifndef PRG_DIRAC
   logical, external :: do_dft
#endif /* PRG_DIRAC */

contains


#ifdef VAR_LINSCA
   function prefix_zeros(n, l)
      integer, intent(in) :: n, l !number, length
      character(l)        :: prefix_zeros !resulting n in ascii
      character(1), parameter :: char0to9(0:9) &
            = (/'0','1','2','3','4','5','6','7','8','9'/)
      integer :: i, k
      k = n
      do i = l, 1, -1
         prefix_zeros(i:i) = char0to9(mod(k,10))
         k = k / 10
      end do
      if (k /= 0) call quit('prefix_zeros error: Argument integer does not fit ' &
                            // 'in the specified number of ASCII caracters')
   end function
#endif /* VAR_LINSCA */


   subroutine prop_oneave(S0, p, D, dime, E, perm, comp, freq, DFD)
   !--------------------------------------------------------------------------
   ! Contracts the 1-electron integrals perturbed by the perturbations p(:)
   ! with the perturbed density matrices in D(:) (e.g. D=(/Dxy/) for a 2nd
   ! order density), and ADDS the result to the property (rsp func) array
   ! E(:). Front for the private subroutine 'oneave' below, checking the
   ! arguments' dimensions, and doing permutations
   ! S0 is passed as argument only as reference to nuclei and basis set
   !--------------------------------------------------------------------------
      type(matrix), intent(in) :: S0   !unperturbed overlap matrix
      character(*), intent(in) :: p(:) !p(np) perturbation lables
      type(matrix), intent(in) :: D(:) !(un)perturbed density matrices to
                        !contract perturbed one-electron integrals with.
                        !If perm present, size(D) = product(dime(perm(np+1:np+nd))),
                        !if perm not present, size(D) = product(dime(np+1:np+nd))
      integer,      intent(in) :: dime(:) !dime(np+nd) = shape(E),
                        !dimensions of property
      !------------------------------------------------------------------------
      complex(8), intent(inout) :: E(*) !property contributions,
                        !works incrementally, thus contributions are
                        !ADDED to E(*). size(E) = product(dime)
      !------------------------------------------------------------------------
      integer,      optional, intent(in) :: perm(:) !perm(np+nd), permutation of indices.
                        !For each dimension of p and D, the corresponding dimension in E.
                        !Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: comp(:) !comp(np), starting component index
                        !for each p. Default 1 1 ... 1
      complex(8),   optional, intent(in) :: freq(:) !freq(np), complex frequencies
                        !for each p, default all zero. Multiply the half-derivative
                        !overlap integrals, thus no contribution if basis independent
                        !of p
      type(matrix), optional, intent(in) :: DFD(:) !optional perturbed energy-weighted
                        !density matrices. Contracted against perturbed overlap integrals
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                      stepe(size(dime)), ddime(size(dime)), idxe(size(dime)), &
                      i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: ffreq(size(p)), Etmp(product(dime)), tmpf
      logical      :: zero, bas
      ! (re-)initialize the module's data on nuclei and basis
      call init_this_module(S0)
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.pert_table(idxp(i))%lin .and. &
             .not.pert_table(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (pert_table(idxp(i))%lin .and. &
             pert_table(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (pert_table(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (pert_table(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! if unperturbed density and anti-symmetric integral, also zero
      j = count((/(pert_table(idxp(i))%anti, i=1,size(p))/))
      zero = (zero .or. (size(p)==size(dime) .and. mod(j,2)==1))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < size(p)) call quit('prop_oneave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < size(p)')
      if (any(dime <= 0)) call quit('prop_oneave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0')
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call quit('prop_oneave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)')
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call quit('prop_oneave ' &
                      // 'argument error: Permutation must contain each number exactly once')
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_oneave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)')
         if (any(comp <= 0)) call quit('prop_oneave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0')
         ccomp = comp
      end if
      if (any(ccomp + ddime(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_oneave argument error: Lowest component index plus ' &
                    // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(p)')
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call quit('prop_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)')
         ffreq = freq
      end if
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddime(k) > ddime(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddime(i);  ddime(i) = ddime(j);  ddime(j) = tmpi
         tmpi = stepe(i);  stepe(i) = stepe(j);  stepe(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
         tmpf = ffreq(i);  ffreq(i) = ffreq(j);  ffreq(j) = tmpf
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(ddime(size(p)+1:size(dime)))
      if (size(D) /= nd) call quit('prop_oneave error: Number of' &
               // 'perturbed densities D does not correspond to dime (and perm)')
      ! verify number of DFD, and that all are defined
      bas = all((/(pert_table(idxp(i))%bas, i=1,size(p))/))
      if (present(DFD)) then
         if (size(DFD) /= nd) call quit('prop_oneave error: Number of' &
               // 'perturbed DFD differs from number of perturbed densities D')
         ! if no basis perturbation (or zero perturbed integrals,
         ! DFD will not be used, so verify they are defined
         if (.not.bas .or. zero) then
            if (.not.all((/(isdef(DFD(i)), i=1,nd)/))) &
               call quit('prop_oneave error: Undefined matrix in argument DFD(:)')
         end if
      end if
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_oneave error: Undefined matrix in argument D(:)')
         return
      end if
      ! everything set up, so call core procedure oneave.
      ! Argument nd=0 is used when averaging over unperturbed density,
      ! in which case also perturbed nuclear attraction should be included
      nd = merge(0, nd, size(p)==size(dime))
      call oneave(S0, size(p), pp, ccomp, ddime(:size(p)), ffreq, nd, D, Etmp, DFD)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      idxe = 0
      do j = 1, product(dime)
         i = 1 + sum(idxe * stepe)
         E(i) = E(i) + Etmp(j)
         do k = 1, size(dime)
            idxe(k) = idxe(k) + 1
            if (idxe(k) /= ddime(k)) exit
            idxe(k) = 0
         end do
      end do
   end subroutine


   subroutine oneave(S0, np, p, c, dp, w, nd, D, E, DFD)
      type(matrix),  intent(in)  :: S0      !unperturbed overlap, to know its dimension
      integer,       intent(in)  :: np      !number of perturbations and order of density
      character(*),  intent(in)  :: p(np)   !perturbation lables
      integer,       intent(in)  :: c(np)   !lowest component of each perturbation
      integer,       intent(in)  :: dp(np)  !dimensions of property integrals
      complex(8),    intent(in)  :: w(np)   !frequency of each p
      integer,       intent(in)  :: nd      !dimensions of property integrals
      type(matrix),  intent(in)  :: D(max(1,nd)) !un-/perturbed density matrices,
                                            !size(D) = product(1+de(np+1:np+nd))
      complex(8),    intent(out) :: E(*)    !resulting one-electron property contributions
                                            !size(E) = product(dp) * nd
      !> un-/perturbed energy-weighted density matrices,
      !> Should have size(DFD) = size(D)
      type(matrix), optional, intent(in) :: DFD(:)
      !-------------------------------------------
      real(8), allocatable :: RR(:) !scratch
      real(8)      :: R(6) !scratch
      integer      :: i, j, k, l, ii, jj, kk, ll
      type(matrix) :: A(6) !scratch matrices
#ifdef PRG_DIRAC
      !> temporary matrix
      type(matrix) :: T
#endif
      A(:) = (/(0d0*S0, i=1,size(A))/) !scratch matrices
      if (np==0) then
         call quit('prop_oneave error: unperturbed one-electron contribution requested')
      else if (np==1 .and. p(1)(1:3)=='AUX') then
         read (p(1)(4:),'(i1)') i
         if (i < 0 .or. i > 9) call quit('prop_oneave error: Index in' &
                  // ' auxiliary label out of range 0..9: ' // p(1))
         call load_oneint(prop_auxlab(i), A(1), S0)
         do j = 0, max(0,nd-1)
            E(1+j) = tr(A(1),D(1+j))
         end do
      else if (np==1 .and. p(1)=='EL') then
         ! contract -dipole integrals '?DIPLEN ' with densities
         do i = 0, dp(1)-1
               ! load -dipole integral from file
!> \todo #ifdef GEN1INT_DALTON
!> \todo                call load_gen1int( xyz(c(1)+i) // 'DIPLEN ', A(1) )
!> \todo #else
               call load_oneint(xyz(c(1)+i) // 'DIPLEN ',A(1))
!> \todo #endif
               ! loop over density matrices
               do j = 0, max(0,nd-1)
                  E(1+i+dp(1)*j) = tr(A(1),D(j+1)) !*2 for total dens hidden in tr
               end do
         end do
         ! no densities means D(1) contains unperturbed density matrix.
         ! Then also add nuclear attraction contribution to -dipole moment
         if (nd==0) then
#if defined(DALTON_AO_RSP)
            call DIPNUC_ifc( natom, R(:3) )
#else
            call DIPNUC_ifc(natom, R(:3), (/(0d0, i=1,9*natom)/))
#endif
            E(:dp(1)) = E(:dp(1)) - R(c(1):c(1)+dp(1)-1) !sign change for -dipole
         end if
      ! Velocity operator. Since anti-symmetric, no unperturbed (nd==0)
      else if (np==1 .and. p(1)=='VEL' .and. nd /= 0) then
         do i = 0, dp(1)-1
            call load_oneint(xyz(c(1)+i) // 'DIPVEL ', A(1))
            do j = 0, nd-1
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
            end do
         end do
      ! No-London magnetic
      else if (np==1 .and. p(1)=='MAGO' .and. nd /= 0) then
         do i = 0, dp(1)-1
            call load_oneint(xyz(c(1)+i) // 'ANGMOM ', A(1))
            do j = 0, nd-1
#ifndef PRG_DIRAC
               E(1+i+dp(1)*j) = tr(A(1), D(j+1)) / 2 !factor 1/2
#else
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
#endif
            end do
         end do
      ! London magnetic, Since anti-symmetric, no unperturbed (nd==0)
      else if (np==1 .and. p(1)=='MAG' .and. nd /= 0) then
         !one-electron contribution
         do i = 0, dp(1)-1
            ! Hamiltonian integrals
            call load_oneint('dh/dB' // xyz(c(1)+i) // '  ', A(1))
            ! -i/2 Tb in A(1), overlap integrals Sb in A(2)
            if (w(1)==0) then !no -i/2 Tb contribution
               call load_oneint('dS/dB' // xyz(c(1)+i) // '  ', A(2))
            else
               call load_oneint('d|S>/dB' // xyz(c(1)+i), A(2))
               A(1) = A(1) - w(1)/2 * A(2)
               A(1) = A(1) - w(1)/2 * dag(A(2))
               A(2) = A(2) - dag(A(2))
            end if
            do j = 0, nd-1
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
               if (present(DFD)) &
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) - tr(A(2), DFD(j+1))
            end do
         end do
      ! Electric field gradient
      else if (np==1 .and. p(1)=='ELGR') then
         ! one-electron integrals contracted with perturbed densities
         do i = 0, dp(1)-1
            ! turn upper triangle 1..6 into square 1..3 x 1..3
            ii = c(1)+i + merge(3, merge(2, 0, c(1)+i > 1), c(1)+i > 3) - 1
            call load_oneint(xyz(1+mod(ii,3)) // xyz(1+ii/3) // 'THETA ', A(1))
            do j = 0, max(0,nd-1)
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
            end do
         end do
         ! if unperturbed density, add nuclear contribution
         if (nd==0) then !to -quadrupole moment
            call QDRNUC_ifc(R(:6))
            E(:dp(1)) = E(:dp(1)) - R(c(1):c(1)+dp(1)-1) !sign change for -theta
         end if
      else if (np==1 .and. p(1)=='GEO') then
         !one-electron integral contribution
         do i = 0, dp(1)-1
            ! perturbed one-electron Hamiltonian integrals
#ifdef PRG_DIRAC
            call load_oneint('G1N' // prefix_zeros(c(1)+i, 3), A(1))
            T = 1.0d0*A(1)
            call load_oneint('G1B' // prefix_zeros(c(1)+i, 3), T)
            A(1) = A(1) + T
            call load_oneint('G1KX' // prefix_zeros(c(1)+i, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 4)
            call load_oneint('G1KY' // prefix_zeros(c(1)+i, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 3)
            call load_oneint('G1KZ' // prefix_zeros(c(1)+i, 3), T)
            call daxpy_iz_jz(1.0d0, T, 1, A(1), 2)
#else
            call load_oneint('1DHAM' // prefix_zeros(c(1)+i, 3), A(1))
#endif
            ! (half-) perturbed overlap -i/2 Tg added to A(1), Sg in A(2)
            if (w(1)==0) then !w=0 means no -i/2 Tg contribution
#ifdef PRG_DIRAC
               call load_oneint('G1O' // prefix_zeros(c(1)+i,3), A(2))
#else
               call load_oneint('1DOVL' // prefix_zeros(c(1)+i,3), A(2))
#endif
               A(2) = -A(2) !1DOVL is really -dS/dg
            else
               call load_oneint('SQHDR' // prefix_zeros(c(1)+i,3), A(2))
               A(2) = -A(2) !SQHDR is really -dS>/dg
               A(1) = A(1) - w(1)/2 * A(2)
               A(1) = A(1) + w(1)/2 * dag(A(2))
               A(2) = A(2) + dag(A(2))
            end if
            do j = 0, max(0,nd-1)
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
               if (present(DFD)) &
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) - tr(A(2), DFD(j+1))
            end do
         end do
         ! nuclear repulsion contribution
         if (nd==0) then
            allocate(RR(3*natom))
            call GRADNN_ifc(natom, RR(:3*natom))
            E(:dp(1)) = E(:dp(1)) + RR(c(1):c(1)+dp(1)-1)
            deallocate(RR)
         end if
      else if (np==2 .and. all(p==(/'MAGO','MAGO'/))) then
#ifdef PRG_DIRAC
!        if diamagnetic contribution comes from positive-negative
!        response, then do not add it here
         if (.not. diamagnetic_via_pn) then
#endif
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(xyz(min(c(1)+i,c(2)+j)) &
                        // xyz(max(c(1)+i,c(2)+j)) // 'SUSCGO', A(1))
               do k = 0, max(0,nd-1)
                  ! ajt Since MAGO is imaginary/anti-symmetri, we get a minus sign here
                  E(1+i+dp(1)*(j+dp(2)*k)) = -tr(A(1),D(1+k))
               end do
            end do
         end do
#ifdef PRG_DIRAC
         endif
#endif
      ! EL MAG has no nd==0 because MAG anti-symmetric
      else if (np==2 .and. all(p==(/'MAG','EL '/)) .and. nd /= 0) then
         do j = 0, dp(2)-1 !EL indices
            do i = 0, dp(1)-1 !MAG indices
               call load_oneint(xyz(c(2)+j) // '-CM1 ' &
                             // xyz(c(1)+i) // ' ', A(1))
               do k = 0, nd-1
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(1), D(k+1))
               end do
            end do
         end do
      else if (np==2 .and. all(p==(/'MAG','MAG'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               ! Hamiltonian integrals and -i/2 Tbb in A(1), Sbb in A(2)
               call load_oneint(xyz(min(c(1)+i,c(2)+j)) &
                        // xyz(max(c(1)+i,c(2)+j)) // 'dh/dB2', A(1))
               ! If both fields static, no -i/2 Tbb to A(1)
               if (w(1)==0 .and. w(2)==0) then
                  call load_oneint('dS/dB2' // xyz(min(c(1)+i,c(2)+j)) &
                                            // xyz(max(c(1)+i,c(2)+j)), A(2))
               else
                  ! only load >> XX XY YY XZ YZ ZZ, use dag(..) for <<
                  call load_oneint('>>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                       // xyz(max(c(1)+i,c(2)+j)), A(2))
                  A(1) = A(1) - (w(1)+w(2))/2 * A(2)
                  A(1) = A(1) + (w(1)+w(2))/2 * dag(A(2))
                  A(2) = A(2) + dag(A(2))
                  ! only load <> XX XY YY XZ YZ ZZ, use dag(..) for the others
                  call load_oneint('<>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                       // xyz(max(c(1)+i,c(2)+j)), A(3))
                  !ajt1009 Changed sign on these two: (think it's correct)
                  A(1) = A(1) - merge(-1, 1, c(1)+i < c(2)+j) &
                                       * (w(1)-w(2))/2 * A(3)
                  A(1) = A(1) + merge(-1, 1, c(1)+i < c(2)+j) &
                                  * (w(1)-w(2))/2 * dag(A(3))
                  A(2) = A(2) + A(3)
                  A(2) = A(2) + dag(A(3))
               end if
               ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
               A(1) = -A(1)
               A(2) = -A(2)
               do k = 0, max(0,nd-1)
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(1),D(k+1))
                  if (present(DFD)) &
                     E(1+i+dp(1)*(j+dp(2)*k)) = E(1+i+dp(1)*(j+dp(2)*k)) &
                                              - tr(A(2), DFD(k+1))
               end do
            end do
         end do
      ! ELGR MAG has no nd==0 because MAG anti-symmetric
      else if (np==2 .and. all(p==(/'ELGR','MAG '/)) .and. nd /= 0) then
         do j = 0, dp(2)-1 !MAG indices
            ! ajt The integrals 'GG-QDB B' are not traceless.
            !     Therefore load all integrals and remove trace manually
            call load_oneint('XX-QDB ' // xyz(c(2)+j), A(1))
            call load_oneint('YY-QDB ' // xyz(c(2)+j), A(3))
            call load_oneint('ZZ-QDB ' // xyz(c(2)+j), A(6))
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint('XY-QDB ' // xyz(c(2)+j), A(2))
            call load_oneint('XZ-QDB ' // xyz(c(2)+j), A(4))
            call load_oneint('YZ-QDB ' // xyz(c(2)+j), A(5))
            do i = 0, dp(1)-1 !ELGR indices
               do k = 0, nd-1 !density indices
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(c(1)+i), D(k+1))
               end do
            end do
            ! ajt Please leave. Old code without trace removal
            ! do i = 0, dp(1)-1 !ELGR indices
            !    ii = c(1)+i + merge(3, merge(2, 0, c(1)+i > 1), c(1)+i > 3) - 1
            !    call load_oneint(xyz(1+mod(ii,3)) // xyz(1+ii/3) &
            !                // '-QDB ' // xyz(c(2)+j), A(1))
            !    do k = 0, nd-1
            !       E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(1),   D(k+1)) &
            !                                - tr(A(2), DFD(k+1))
            !    end do
            ! end do
         end do
      else if (np==2 .and. all(p==(/'GEO','EL '/))) then
         !contract perturbed densities with -dipole gradient integrals
         do j = 0, dp(2)-1 !EL indices
            do i = 0, dp(1)-1 !GEO indices
               call load_oneint(prefix_zeros(c(1)+i,2) // ' DPG ' &
                           // xyz(c(2)+j), A(1))
               do k = 0, max(0,nd-1) !density indices
                  E(1+i+dp(1)*(j+dp(2)*k)) = -tr(A(1), D(k+1)) !minus sign
               end do
            end do
         end do
         ! nd==0 has nuclear repulsion contribution to -dipole gradient
         if (nd==0) then
            allocate(RR(3*(3*natom)))
            call DIPNUC_ifc(natom, R(:3), RR(:3*(3*natom)))
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) &
                                 - RR(c(2)+j + 3*(c(1)+i-1)) !-dipole sign
               end do
            end do
            deallocate(RR)
         end if
      ! GEO MAGO has no nd==0 because MAGO anti-symmetric
      else if (np==2 .and. all(p==(/'GEO ','MAGO'/)) .and. nd /= 0) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(prefix_zeros(c(1)+i,3) // 'AMDR' &
                           // xyz(c(2)+j), A(1))
               do k = 0, nd-1
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(1), D(k+1)) / 2 !factor 1/2
               end do
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO ','ELGR'/))) then
         do i = 0, dp(1)-1 !GEO indices
            ! ajt The integrals '01QDG XY' are not traceless. Therefore,
            !     load all components at once and remove trace before contracting
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XX', A(1))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG YY', A(3))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG ZZ', A(6))
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XY', A(2))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XZ', A(4))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG YZ', A(5))
            do j = 0, dp(2)-1 !ELGR indices
               do k = 0, max(0,nd-1) !density indices
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(c(2)+j), D(k+1)) * (-3/2d0) !factor -3/2
               end do
            end do
         end do
         ! ajt Please leave. Old code without trace removal
         ! do j = 0, de(2)-1
         !    jj = c(2)+j + merge(3, merge(2, 0, c(2)+j > 1), c(2)+j > 3) - 1
         !    do i = 0, de(1)-1
         !       call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG ' &
         !                   // xyz(1+mod(jj,3)) // xyz(1+jj/3), A(1))
         !       do k = 0, max(0,nd-1) !density indices
         !          E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(c(1)+i), D(k+1)) * (-3/2d0) !factor -3/2
         !       end do
         !    end do
         ! end do
         ! nd==0, nuclear repulsion contribution to -quadrupole gradient
         if (nd==0) then
             call quit('prop_oneave: ELGR GEO - nuc.rep. conrib. not implemented')
         end if
      else if (np==2 .and. all(p==(/'GEO','GEO'/))) then
         allocate(RR((3*natom)*(3*natom)))
         ! one-electron integral contribution
         do k = 0, max(0,nd-1)
            !ajt fixme oneint_ave(GG..) misses frequency dependent -i/2 Tgg contribution
            if (present(DFD)) &
               call oneint_ave('GG', D(k+1), DFD(k+1), RR)
            if (.not.present(DFD)) &
               call oneint_ave('GG', D(k+1), 0d0*S0, RR)
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*(j+dp(2)*k)) = RR(c(1)+i + 3*natom*(c(2)+j-1))
               end do
            end do
         end do
         ! for zero-order density, add nuclear contribution to Hessian
         if (nd==0) then
            call HESSNN_ifc(natom, RR)
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) + RR(c(1)+i + 3*natom*(c(2)+j-1))
               end do
            end do
         end if
         deallocate(RR)
      else if (np==3 .and. all(p==(/'MAG','MAG','EL '/))) then
         do k = 0, dp(3)-1
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  call load_oneint(xyz(c(3)+k) // '-CM2' // xyz(min(c(1)+i,c(2)+j)) &
                              // xyz(max(c(1)+i,c(2)+j)) // ' ', A(1))
                  do l = 0, max(0,nd-1)
                     ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
                     E(1+i+dp(1)*(j+dp(2)*(k+dp(3)*l))) = -tr(A(1), D(1+l))
                  end do
               end do
            end do
         end do
      else
         print *,'prop_oneave: no integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call quit('prop_oneave: no such integrals')
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine


   subroutine prop_twoave(p, D, dime, E, perm, comp)
   !--------------------------------------------------------------------------
   ! Contracts the 2-electron integrals perturbed by the perturbations p(:)
   ! with the perturbed density matrix expansion in D(:) (e.g. D=(/D,Dx,Dy,Dxy/)
   ! for a 2nd order expansion), and ADDS the result to the property (rsp func)
   ! array E(:). Front for the private subroutine 'twoave' below, checking the
   ! arguments' dimensions, and doing permutations
   ! D(1) serves as reference to nuclei, basis and model/functional
   !--------------------------------------------------------------------------
      character(*), intent(in) :: p(:) !p(np) perturbation lables
      type(matrix), intent(in) :: D(:) !(un)perturbed density matrices to
                        !contract perturbed one-electron integrals with.
                        !If perm present, size(D) = product(dime(perm(np+1:np+nd))),
                        !if perm not present, size(D) = product(dime(np+1:np+nd))
      integer,      intent(in) :: dime(:) !dime(np+nd) = shape(E),
                        !dimensions of property
      !------------------------------------------------------------------------
      complex(8), intent(inout) :: E(*) !property contributions,
                        !works incrementally, thus contributions are
                        !ADDED to E(*). size(E) = product(dime)
      !------------------------------------------------------------------------
      integer,      optional, intent(in) :: perm(:) !perm(np+nd), permutation of indices.
                        !For each dimension of p and D, the corresponding dimension in E.
                        !Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: comp(:) !comp(np), starting component index
                        !for each p. Default 1 1 ... 1
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                      stepe(size(dime)), ddime(size(dime)), idxe(size(dime)), &
                      i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: Etmp(product(dime))
      logical      :: zero
      ! (re-)initialize the module's data on nuclei and basis
      call init_this_module(D(1))
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero.
      zero = any((/(pert_table(idxp(i))%lin  .or. &
                    pert_table(idxp(i))%quad .or. &
               .not.pert_table(idxp(i))%bas,  i=1,size(p))/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < size(p)) call quit('prop_twoave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < size(p)')
      if (any(dime <= 0)) call quit('prop_twoave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0')
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call quit('prop_twoave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)')
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call quit('prop_twoave ' &
                      // 'argument error: Permutation must contain each number exactly once')
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_twoave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)')
         if (any(comp <= 0)) call quit('prop_twoave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0')
         ccomp = comp
      end if
      if (any(ccomp + ddime(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_twoave argument error: Lowest component index plus ' &
                    // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(p)')
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddime(k) > ddime(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddime(i);  ddime(i) = ddime(j);  ddime(j) = tmpi
         tmpi = stepe(i);  stepe(i) = stepe(j);  stepe(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(1+ddime(size(p)+1:size(dime)))
      if (size(D) /= nd) call quit('prop_twoave error: Number of' &
               // 'perturbed densities D does not correspond to dime (and perm)')
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_twoave error: Undefined matrix in argument D(:)')
         return
      end if
      ! everything set up, so call core procedure oneave.
      ! Argument nd=0 is used when averaging over unperturbed density,
      ! in which case also perturbed nuclear attraction should be included
      call twoave(size(p), size(dime)-size(p), pp, ccomp, ddime, D, Etmp)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      idxe = 0
      do j = 1, product(dime)
         i = 1 + sum(idxe * stepe)
         E(i) = E(i) + Etmp(j)
         do k = 1, size(dime)
            idxe(k) = idxe(k) + 1
            if (idxe(k) /= ddime(k)) exit
            idxe(k) = 0
         end do
      end do
   end subroutine


   subroutine twoave(np, nd, p, c, de, D, E)
      integer,       intent(in) :: np, nd    !number of perturbations and order of density
      character(*),  intent(in) :: p(np)     !perturbation lables
      integer,       intent(in) :: c(np)     !lowest component of each perturbation
      integer,       intent(in) :: de(np+nd) !dimensions of property (E)
      type(matrix),  intent(in) :: D(*)      !un-/perturbed density matrices (expension),
                                             !size(D) = product(1+de(np+1:np+nd))
      complex(8),   intent(out) :: E(*)      !where to ADD property contributions
                                             !(works incrementally), size(E) = product(de)
      real(8), allocatable :: RR(:) !scratch
      real(8)      :: R(6) !scratch
      integer      :: i, j, k, l, ii, jj, kk, ll, pd, pd1
      type(matrix) :: A(6) !scratch matrices
      A(:) = (/(0d0*D(1), i=1,size(A))/) !scratch matrices
      pd  = product(de(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+de(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call quit('prop_twoave: unperturbed energy/integrals requested')
         !The highest order Ds multiply the unperturbed F, which is unknown here,
         !so insist that no highest-order Ds are present
         if (.not.all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) &
            call quit('prop_twoave: unperturbed Fock matrix requested')
         !ajt fixme Stopping at nd=2, for now
         if (nd==2) then
            ! contract first density to Fock matrix, then trace with last
            do i = 0, de(1)-1
               ! Coulomb and exchange
               call twofck('  ', D(2+i), A(1:1))
               ! Kohn-Sham exchange-correlation
               call twofck_ks(1, (/D(1),D(2+i)/), A(1))
               ! trace with first density matrix
               do j = 0, de(2)-1
                  E(1+i+de(1)*j) = tr(A(1),D(2+de(1)+j))
               end do
            end do
         else if (nd==3) then
            ! contract two of first-order densities to Fock, then trace with third
            do j = 0, de(2)-1
               do i = 0, de(1)-1
                  call twofck_ks(2, (/D(1),D(2+i),D(2+de(1)+j)/), A(1))
                  do k = 0, de(3)-1
                     E(1+i+de(1)*(j+de(2)*k)) = tr(A(1),D(2+de(1)+de(2)+k))
                  end do
                  A(1) = 0d0*A(1)
               end do
            end do
            ! contract first-order densities to Fock, then trace with second-order
            !ajt Fixme
            if (.not.all((/(iszero(D(i)), i=2+de(1)+de(2)+de(3),pd1)/))) &
               call quit('prop_twoave: nd = 3 not fully implemented')
         else if (nd > 3) then
            call quit('prop_twoave: nd > 3 not implemented')
         end if
      !London magnetic, no nd==0 because because MAG anti
      else if (np==1 .and. p(1)=='MAG' .and. nd /= 0) then
         ! di_get_MagDeriv_(F,G)xD_DFT takes initialized
         if (do_dft()) then !scratch in A(4:6)
            do i = 4, 6
               call init_mat(A(i), A(i))
            end do
         end if
         if (all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) then
            E(:product(de)) = 0
         else
            ! contract unperturbed density to Fock, then trace with highest-order
            ! Coulomb-exchange
            call twofck('M ', D(1:1), A(1:3))
            ! Kohn-Sham contribution
            if (do_dft()) then
               call di_get_MagDeriv_FxD_DFT(A(4:6), D(1))
               do i = 1, 3
                  A(i) = A(i) + A(3+i)
               end do
            end if
            do j = 0, pd-1
               do i = 0, de(1)-1
                  E(1+i+de(1)*j) = tr(A(c(1)+i), D(pd1-pd+1+j))
               end do
            end do
         end if
         if (nd==1) then
            ! nothing more
         else if (nd==2) then
            ! contract first density to Fock matrix, then trace with second
            do j = 0, de(2)-1
               if (iszero(D(2+j))) cycle
               ! Coulomb-exchange
               call twofck('M ', D(2+j:2+j), A(1:3))
               ! Kohn-Sham
               if (do_dft()) then
                  call di_get_MagDeriv_GxD_DFT(D(1), D(2+j), A(4:6))
                  do i = 1, 3
                     A(i) = A(i) - A(3+i) !negative sign on KS contrib
                  end do
               end if
               ! trace with second density matrix
               do k = 0, de(3)-1
                  do i = 0, de(1)-1
                     E(1+i+de(1)*(j+de(2)*k)) = E(1+i+de(1)*(j+de(2)*k)) &
                                              + tr(A(c(1)+i), D(2+de(2)+k))
                  end do
               end do
            end do
         else
            call quit('prop_twoave: MAG, nd > 2 not implemented')
         end if
      else if (np==1 .and. p(1)=='GEO') then
         allocate(RR(3*natom*2))
         ! highest-order contribution
         do j = 0, pd-1
            if (iszero(D(pd1-pd+1+j))) then
               E( 1+de(1)*j : de(1)*(j+1) ) = 0
               cycle
            end if
            ! Coulomb-exchange
            call twoctr('G', D(1), D(pd1-pd+1+j), RR(:3*natom))
            ! for molgra (nd==0), factor 1/2 on these integrals
            if (nd==0) RR(:3*natom) = RR(:3*natom)/2
            if (do_dft()) then
               if (nd==0) call di_get_geomDeriv_molgrad_DFT( &
                                          RR(3*natom+1:3*natom*2), &
                                          natom, D(1))
               if (nd/=0) call di_get_geomDeriv_FxD_DFT(           &
                                          RR(3*natom+1:3*natom*2), &
                                          natom, D(1), D(pd1-pd+1+j))
               RR(:3*natom) = RR(:3*natom) + RR(3*natom+1:3*natom*2)
            end if
            E( 1+de(1)*j : de(1)*(j+1) ) = RR(c(1):c(1)+de(1)-1)
         end do
         if (nd==0 .or. nd==1) then
            ! nothing more
         else if (nd==2) then
            ! integrals over products of first order densities
            do k = 0, de(3)-1
               do j = 0, de(2)-1
                  ! Coulomb-exchange
                  call twoctr('G', D(2+j), D(2+de(2)+k), RR(:3*natom))
                  if (do_dft()) then
                     call di_get_geomDeriv_GxD_DFT(RR(3*natom+1:3*natom*2), &
                                                   natom, D(1), D(2+j), D(2+de(2)+k))
                     RR(:3*natom) = RR(:3*natom) + RR(3*natom+1:3*natom*2)
                  end if
                  E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            = E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            + RR(c(1):c(1)+de(1)-1)
               end do
            end do
         else
            call quit('prop_energy_reort: GEO, nd > 2 not implemented')
         end if
         deallocate(RR)
         if (do_dft()) print* !after all the "...integrated to nn electrons..." prints
      else if (np==2 .and. all(p==(/'MAG','MAG'/))) then
         ! highest-order contribution
         if (nd /= 0 .and. all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) then
            E(:de(1)*de(2)*pd) = 0
         else
            ! Coulomb-exchange
            call twofck('MM', D(1:1), A(1:6))
            ! Kohn-Sham exchange-correlation
            if (do_dft()) then
               call quit('prop_twoint: MAG MAG not implemented for DFT')
            end if
            do k = 0, pd-1
               do j = 0, de(2)-1
                  do i = 0, de(1)-1
                     ii = min(c(1)+i,c(2)+j) &
                        + max(c(1)+i,c(2)+j) * (max(c(1)+i,c(2)+j)-1) / 2
                     ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
                     E(1+i+de(1)*(j+de(2)*k)) = -tr(A(ii), D(pd1-pd+1+k)) &
                                              * merge(1/2d0,1d0,nd==0) !half for magnetizability
                  end do
               end do
            end do
         end if
         if (nd==0 .or. nd==1) then
             !nothing more
         else if (nd==2) then
            ! integrals over products of first order densities
            ! Contract first density to Fock, then trace with second
            do k = 0, de(3)-1
               if (iszero(D(2+k))) cycle
               call twofck('MM', D(2+k), A(:6))
               if (do_dft()) then
                  call quit('prop_energy_reort: MAG MAG not implemented for DFT')
               end if
               do l = 0, de(4)-1
                  do j = 0, de(2)-1
                     do i = 0, de(1)-1
                        ii = min(c(1)+i,c(2)+j) &
                           + max(c(1)+i,c(2)+j) * (max(c(1)+i,c(2)+j)-1) / 2
                        ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
                        E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                    = E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                    - tr(A(ii),D(2+de(3)+l))
                     end do
                  end do
               end do
            end do
         else
            call quit('prop_twoint: MAG MAG and nd > 2 not implemented')
         end if
      else if (np==2 .and. all(p==(/'GEO','GEO'/))) then
         allocate(RR((3*natom)*(3*natom)))
         ! two-electron integral contribution
         do k = 0, pd-1
            if (iszero(D(pd1-pd+1+k))) then
               E( 1+de(1)*de(2)*k : de(1)*de(2)*(k+1) ) = 0
               cycle
            end if
            ! Coulomb-exchange
            call twoctr('GG', D(1), D(pd1-pd+1+k), RR)
            if (nd==0) RR = RR/2 !factor 1/2 for unperturbed Hessian integrals
            ! Kohn-Sham exchange-correlation
            if (do_dft()) then
               call quit('prop_twoint: GEO GEO, DFT not implemented')
            end if
            do j = 0, de(2)-1
               do i = 0, de(1)-1
                  E(1+i+de(1)*(j+de(2)*k)) = RR(c(1)+i + 3*natom*(c(2)+j-1))
               end do
            end do
         end do
         if (nd==0 .or. nd==1) then
            ! nothing more
         else if (nd==2) then
            ! integrals over products of first order densities
            do l = 0, de(4)-1
               if (iszero(D(2+de(3)+l))) cycle
               do k = 0, de(3)-1
                  if (iszero(D(2+k))) cycle
                  call twoctr('GG', D(2+k), D(2+de(3)+l), RR)
                  if (do_dft()) then
                     call quit('prop_energy_reort: GEO GEO, DFT not implemented')
                  end if
                  do j = 0, de(2)-1
                     do i = 0, de(1)-1
                        E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                     = E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                     + RR(c(1)+i + 3*natom*(c(2)+j-1))
                     end do
                  end do
               end do
            end do
         else
            call quit('prop_twoave: GEO GEO, nd > 2 not implemented')
         end if
         deallocate(RR)
      else
         print *,'prop_twoave: no integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call quit('prop_twoave: no such integrals')
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine


   subroutine prop_oneint(S0, p, dimp, F, S, comp, freq)
   !--------------------------------------------------------------------------
   ! Calculates the 1-electron integrals perturbed by the perturbations p(:)
   ! and ADDS the integrals to the perturbed Fock matrices F(:).
   ! Front for the private subroutine 'oneint' below, checking the
   ! arguments' dimensions, and doing permutations
   ! S0 is passed as argument only as reference to nuclei and basis set
   !--------------------------------------------------------------------------
      type(matrix), intent(in) :: S0   !unperturbed overlap matrix
      character(*), intent(in) :: p(:) !perturbation lables
      integer,      intent(in) :: dimp(:) !shape(F), size(dimp) = size(p)
                        !dimensions of perturbed Fock matrices F(:)
      !------------------------------------------------------------------------
      type(matrix), optional, intent(inout) :: F(*) !perturbed Fock matrices
                        !to fill with perturbed integrals, size(F) = product(dimp)
      type(matrix), optional, intent(inout) :: S(*) !optionally return the
                        !corresponding perturbed overlap matrices
      integer,      optional, intent(in)    :: comp(:) !comp(np), starting
                        !component index for each p. Default 1 1 ... 1
      complex(8),   optional, intent(in)    :: freq(:) !freq(np), complex frequencies
                        !for each p, default all zero. These multiply the
                        !half-derivative overlap integrals, thus no contribution
                        !if basis independent of p
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dimp)), ccomp(size(p)), &
                      stepf(size(dimp)), ddimp(size(dimp)), idxf(size(dimp)), &
                      i, j, k, tmpi
      type(matrix) :: Ftmp(product(dimp)), Stmp(product(dimp))
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: ffreq(size(p)), tmpf
      logical      :: zero, bas
      ! (re-)initialize the module's data on nuclei and basis
      call init_this_module(S0)
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.pert_table(idxp(i))%lin .and. &
             .not.pert_table(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (pert_table(idxp(i))%lin .and. &
             pert_table(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (pert_table(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (pert_table(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! check dimensions argument dimp, verify that dimensions are positive
      if (size(dimp) /= size(p)) call quit('prop_oneint argument error: ' &
               // 'Different number of perturbations and dimensions, size(dimp) /= size(p)')
      if (any(dimp <= 0)) call quit('prop_oneint argument error: ' &
               // 'Perturbations have a zero or negative dimension, dimp <= 0')
      ! compute step lengths in F (cumulative product of dimp)
      stepf(1) = 1
      do i = 2, size(dimp)
         stepf(i) = stepf(i-1)*dimp(i-1)
      end do
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_oneint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)')
         if (any(comp <= 0)) call quit('prop_oneint argument error: ' &
                    // 'Lowest component indices must be positive, but comp <= 0')
         ccomp = comp
      end if
      if (any(ccomp + dimp - 1 > pert_shape(p))) &
         call quit('prop_oneave argument error: Lowest component index plus ' &
                    // 'dimension exceeds dimension of perturbation, comp + dimp > pert_shape(p)')
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call quit('prop_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)')
         ffreq = freq
      end if
      ! sort perturbations p so that idxp is descending
      pp = p
      ddimp = dimp
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddimp)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddimp(k) > ddimp(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddimp(i);  ddimp(i) = ddimp(j);  ddimp(j) = tmpi
         tmpi = stepf(i);  stepf(i) = stepf(j);  stepf(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
         tmpf = ffreq(i);  ffreq(i) = ffreq(j);  ffreq(j) = tmpf
      end do
      ! zero (deallocate) F and S before calling oneint
      do i = 1, product(dimp)
         if (present(F)) F(i) = 0d0*S0
         if (present(S)) S(i) = 0d0*S0
      end do
      ! early return if the result is zero
      if (zero) return
      ! everything set up, so call core procedure oneint
      call oneint(S0, size(p), pp, ccomp, ddimp, ffreq, Ftmp, Stmp)
      ! add oneavg property contribution in temporary array Etmp(:) to
      ! resulting property array E(:), while permuting indices according to pperm
      bas = all(pert_basdep(p))
      idxf = 0
      do j = 1, product(dimp)
         i = 1 + sum(idxf * stepf)
         if (.not.present(F)) Ftmp(j) = 0
         if (.not.present(S)) Stmp(j) = 0
         if (present(F)) &
            call init_mat(F(i), Ftmp(j), alias='FF')
         if (present(S) .and. bas) &
            call init_mat(S(i), Stmp(j), alias='FF')
         do k = 1, size(dimp)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimp(k)) exit
            idxf(k) = 0
         end do
      end do
   end subroutine


   subroutine oneint(S0, np, p, c, dp, w, F, S)
      type(matrix), intent(in) :: S0      !unperturbed overlap matrix
      integer,      intent(in) :: np      !number of perturbations
      character(*), intent(in) :: p(:)    !perturbation lables
      integer,      intent(in) :: c(np)   !lowest component of each perturbation
      integer,      intent(in) :: dp(np)  !dimensions of property (F and S)
      complex(8),   intent(in) :: w(np)   !frequency of each p
      type(matrix), intent(inout) :: F(*) !perturbed Fock matrices
      type(matrix), intent(inout) :: S(*) !perturbed overlap matrices
      integer      :: i, j, k, ii, jj, kk
      type(matrix) :: A(6) !scratch matrices
      if (np==0) then
         call quit('prop_oneint error:' &
                // ' Unperturbed one-electron integrals requested')
      else if (np==1 .and. p(1)(1:3)=='AUX') then
         read (p(1)(4:),'(i1)') i
         if (i < 0 .or. i > 9) call quit('prop_oneint error:' &
                  // ' Index in auxiliary label out of range 0..9')
         call load_oneint(prop_auxlab(i), F(1), S0)
      else if (np==1 .and. p(1)=='EL') then
         do i = 0, dp(1)-1
!> \todo #ifdef GEN1INT_DALTON
!> \todo             call load_gen1int( xyz(c(1)+i) // 'DIPLEN ', F(i+1), S0 )
!> \todo #else
            call load_oneint(xyz(c(1)+i) // 'DIPLEN ', F(i+1), S0)
!> \todo #endif
         end do
      else if (np==1 .and. p(1)=='VEL') then
         do i = 0, dp(1)-1
            call load_oneint(xyz(c(1)+i) // 'DIPVEL ', F(i+1), S0)
         end do
      else if (np==1 .and. p(1)=='MAGO') then
         do i = 0, dp(1)-1
            call load_oneint(xyz(c(1)+i) // 'ANGMOM ', F(i+1), S0)
#ifndef PRG_DIRAC
            F(i+1) = (1/2d0) * F(1+i) !factor 1/2 here
#endif
         end do
      else if (np==1 .and. p(1)=='MAG') then
         do i = 0, dp(1)-1
            call load_oneint('dh/dB' // xyz(c(1)+i) // '  ', F(i+1), S0)
            if (w(1)==0) then !no -i/2 Tb contribution
               call load_oneint('dS/dB' // xyz(c(1)+i) // '  ', S(i+1), S0)
            else
               call load_oneint('d|S>/dB' // xyz(c(1)+i), S(i+1), S0)
               F(i+1) = F(i+1) - w(1)/2 * S(i+1)
               F(i+1) = F(i+1) - w(1)/2 * dag(S(i+1))
               S(i+1) = S(i+1) - dag(S(i+1))
            end if
         end do
      else if (np==1 .and. p(1)=='ELGR') then
         do i = 0, dp(1)-1
            ii = c(1)+i + merge(3, merge(2, 0, c(1)+i > 1), c(1)+i > 3) - 1
            call load_oneint(xyz(1+mod(ii,3)) // xyz(1+ii/3) // 'THETA ', &
                             F(1+i), S0)
         end do
      else if (np==1 .and. p(1)=='GEO') then
         do i = 0, dp(1)-1
#ifdef PRG_DIRAC
            call quit('implement 1dham')
#else
            call load_oneint('1DHAM' // prefix_zeros(c(1)+i,3), F(i+1), S0)
#endif
            if (w(1)==0) then !no -i/2 Tb contribution
               call load_oneint('1DOVL' // prefix_zeros(c(1)+i,3), S(i+1), S0)
               S(i+1) = -S(i+1) !1DOVL is really -dS/dg
            else
               call load_oneint('SQHDR' // prefix_zeros(c(1)+i,3), S(i+1), S0)
               S(i+1) = -S(i+1) !SQHDR is really -dS>/dg
               F(i+1) = F(i+1) - w(1)/2 * S(i+1)
               F(i+1) = F(i+1) + w(1)/2 * dag(S(i+1))
               S(i+1) = S(i+1) + dag(S(i+1))
            end if
         end do
      else if (np==2 .and. all(p==(/'MAG','EL '/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(xyz(c(2)+j) // '-CM1 ' &
                             // xyz(c(1)+i) // ' ', F(1+i+dp(1)*j), S0)
            end do
         end do
      else if (np==2 .and. all(p==(/'MAGO','MAGO'/))) then
#ifdef PRG_DIRAC
!        if diamagnetic contribution comes from positive-negative
!        response, then do not add it here
         if (.not. diamagnetic_via_pn) then
#endif
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(xyz(min(c(1)+i,c(2)+j))              &
                             // xyz(max(c(1)+i,c(2)+j)) // 'SUSCGO', &
                                F(1+i+dp(1)*j), S0)
               !ajt Since MAGO is imaginary/anti-symmetric, we get a minus sign here
               F(1+i+dp(1)*j) = -F(1+i+dp(1)*j)
            end do
         end do
#ifdef PRG_DIRAC
         endif
#endif
      else if (np==2 .and. all(p==(/'MAG','MAG'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               ii = 1+i+dp(1)*j
               ! Hamiltonian integrals and -i/2 Tbb in A(1), Sbb in A(2)
               call load_oneint(xyz(min(c(1)+i,c(2)+j)) &
                             // xyz(max(c(1)+i,c(2)+j)) // 'dh/dB2', F(ii), S0)
               ! If both fields static, no -i/2 Tbb to A(1)
               if (w(1)==0 .and. w(2)==0) then
                  call load_oneint('dS/dB2' // xyz(min(c(1)+i,c(2)+j))  &
                                            // xyz(max(c(1)+i,c(2)+j)), S(ii), S0)
               else
                  !only load >> XX XY YY XZ YZ ZZ, use dag(..) for <<
                  call load_oneint('>>S/B2' // xyz(min(c(1)+i,c(2)+j))  &
                                            // xyz(max(c(1)+i,c(2)+j)), S(ii), S0)
                  F(ii) = F(ii) - (w(1)+w(2))/2 * S(ii)
                  F(ii) = F(ii) + (w(1)+w(2))/2 * dag(S(ii))
                  S(ii) = S(ii) + dag(S(ii))
                  !only load <> XX XY YY XZ YZ ZZ, use dag(..) for the others
                  call load_oneint('<>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                            // xyz(max(c(1)+i,c(2)+j)), A(1), S0)
                  !ajt&kk 0410 Changed sign on these two: (think it's correct)
                  F(ii) = F(ii) - merge(-1, 1, c(1)+i < c(2)+j) &
                                         * (w(1)-w(2))/2 * A(1)
                  F(ii) = F(ii) + merge(-1, 1, c(1)+i < c(2)+j) &
                                         * (w(1)-w(2))/2 * dag(A(1))
                  S(ii) = S(ii) + A(1)
                  S(ii) = S(ii) + dag(A(1))
                  A(1) = 0 !free
               end if
               !ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
               F(ii) = -F(ii)
               S(ii) = -S(ii)
            end do
         end do
      else if (np==2 .and. all(p==(/'ELGR','MAG '/))) then
         do j = 0, dp(2)-1 !MAG indices
            ! ajt The integrals 'XY-QDB Z' are not traceless like THETA.
            !     Therefore load all integrals and remove trace manually
            call load_oneint('XX-QDB ' // xyz(c(2)+j), A(1), S0)
            call load_oneint('YY-QDB ' // xyz(c(2)+j), A(3), S0)
            call load_oneint('ZZ-QDB ' // xyz(c(2)+j), A(6), S0)
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint('XY-QDB ' // xyz(c(2)+j), A(2), S0)
            call load_oneint('XZ-QDB ' // xyz(c(2)+j), A(4), S0)
            call load_oneint('YZ-QDB ' // xyz(c(2)+j), A(5), S0)
            do i = 0, dp(1)-1 !ELGR indices
               F(1+i+dp(1)*j) = A(c(1)+i)
            end do
         end do
         if (dp(2) /= 0) A(:) = 0 !free
      else if (np==2 .and. all(p==(/'GEO','EL '/))) then
         do j = 0, dp(2)-1 !EL indices
            do i = 0, dp(1)-1 !GEO indices
               call load_oneint(prefix_zeros(c(1)+i,2) // ' DPG ' &
                                // xyz(c(2)+j), F(1+i+dp(1)*j), S0)
               F(1+i+dp(1)*j) = -F(1+i+dp(1)*j) !minus sign
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO ','MAGO'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(prefix_zeros(c(1)+i,3) // 'AMDR' &
                                // xyz(c(2)+j), F(1+i+dp(1)*j), S0)
               F(1+i+dp(1)*j) = 1/2d0 * F(1+i+dp(1)*j) !factor 1/2
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO ','ELGR'/))) then
         do i = 0, dp(1)-1 !GEO indices
            ! ajt The integrals '01QDG XY' are not traceless.
            !     Therefore load all integrals and remove trace manually
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XX', A(1))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG YY', A(3))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG ZZ', A(6))
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XY', A(2))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG XZ', A(4))
            call load_oneint(prefix_zeros(c(1)+i,2) // 'QDG YZ', A(5))
            do j = 0, dp(2)-1 !ELGR indices
               F(1+i+dp(1)*j) = -3/2d0 * A(c(2)+j) !factor -3/2
            end do
         end do
         if (dp(1) /= 0) A(:) = 0
      else if (np==3 .and. all(p==(/'MAG','MAG','EL '/))) then
         do k = 0, dp(3)-1
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  call load_oneint(xyz(c(3)+k) // '-CM2'              &
                                   // xyz(min(c(1)+i,c(2)+j))         &
                                   // xyz(max(c(1)+i,c(2)+j)) // ' ', &
                                   F(1+i+dp(1)*(j+dp(2)*k)))
                  !ajt Since MAGO is imaginary/anti-symmetric, we get a minus sign here
                  F(1+i+dp(1)*(j+dp(2)*k)) = -F(1+i+dp(1)*(j+dp(2)*k))
               end do
            end do
         end do
      else
         print *,'prop_oneint: No integrals for these perturbations:', &
                    (' ' // p(i),i=1,np)
         call quit('prop_oneint: No such integrals')
      end if
   end subroutine


   subroutine prop_twoint(p, D, dimf, F, perm, comp)
   !--------------------------------------------------------------------------
   ! Contracts the 2-electron and Kohn-Sham integrals perturbed by the
   ! perturbations p(:) with the perturbed density matrix expansion in D(:)
   ! (e.g. D=(/D,Dx,Dy,Dxy/) for a 2nd order expansion), and ADD the resulting
   ! Fock matrix contibution to the array F(:).
   ! Front for the private subroutine 'twoave' below, checking the arguments'
   ! dimensions, and doing permutations
   ! D(1) serves as reference to nuclei, basis and model/functional
   !--------------------------------------------------------------------------
      character(*), intent(in) :: p(:) !p(np) perturbation lables
      type(matrix), intent(in) :: D(:) !(un)perturbed density matrices to
                        !contract perturbed one-electron integrals with.
                        !If perm present, size(D) = product(1+dime(perm(np+1:np+nd))),
                        !if perm not present, size(D) = product(1+dime(np+1:np+nd))
      integer,      intent(in) :: dimf(:) !dime(np+nd) = shape(F),
                        !dimensions of perturbed Fock matrices F(:)
      !------------------------------------------------------------------------
      type(matrix), intent(inout) :: F(*) !Perturbed Fock matrices,
                        !works incrementally, thus contributions are
                        !ADDED to F(*). size(F) = product(dimf)
      !------------------------------------------------------------------------
      integer,      optional, intent(in) :: perm(:) !perm(np+nd), permutation of indices.
                        !For each dimension of p and D, the corresponding dimension of F.
                        !Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: comp(:) !comp(np), starting component index
                        !for each p. Default 1 1 ... 1
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dimf)), ccomp(size(p)), &
                      stepf(size(dimf)), ddimf(size(dimf)), idxf(size(dimf)), &
                      i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      type(matrix) :: Ftmp(product(dimf))
      logical      :: zero
      ! (re-)initialize the module's data on nuclei and basis
      call init_this_module(D(1))
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero.
      zero = any((/(pert_table(idxp(i))%lin  .or. &
                    pert_table(idxp(i))%quad .or. &
               .not.pert_table(idxp(i))%bas, i=1,size(p))/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dimf) < size(p)) call quit('prop_twoint argument error: ' &
               // 'More perturbations than dimensions of F(:), size(dimf) < size(p)')
      if (any(dimf <= 0)) call quit('prop_twoint argument error: ' &
               // 'Perturbed Fock F(:) has a zero or negative dimension, dimf <= 0')
      ! compute step lengths in E (cumulative products of dimensions)
      stepf(1) = 1
      do i = 2, size(dimf)
         stepf(i) = stepf(i-1)*dimf(i-1)
      end do
      ! reorder dimensions and step lengths in F according to permutation argument perm
      ddimf = dimf
      if (present(perm)) then
         if (size(perm) /= size(dimf)) call quit('prop_twoint argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dimf)')
         ! verify that perm is indeed a permutation
         do i = 1, size(dimf)-1
            if (perm(i) <= 0 .or. perm(i) > size(dimf) .or. &
                any(perm(i) == perm(i+1:size(dimf)))) call quit('prop_twoint ' &
                      // 'argument error: Permutation must contain each number exactly once')
         end do
         ddimf = (/( dimf(perm(i)), i=1,size(dimf))/)
         stepf = (/(stepf(perm(i)), i=1,size(dimf))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_twoint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)')
         if (any(comp <= 0)) call quit('prop_twoint argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0')
         ccomp = comp
      end if
      if (any(ccomp + ddimf(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_twoint argument error: Lowest component index plus ' &
                   // 'dimension exceeds dimension of perturbation, comp + dimf > pert_shape(p)')
      ! sort perturbations p so that idxp is descending
      pp = p
      do i = 1, size(p)
         !find perturbation after i, with highest idxp (secondly highest ddime)
         j = i
         do k = j+1, size(p)
            if (idxp(k) > idxp(j) .or. (idxp(k) == idxp(j) .and. &
                ddimf(k) > ddimf(j))) j = k
         end do
         !swap all entries i and j
         tmpp =    pp(i);     pp(i) =    pp(j);     pp(j) = tmpp
         tmpi =  idxp(i);   idxp(i) =  idxp(j);   idxp(j) = tmpi
         tmpi = ddimf(i);  ddimf(i) = ddimf(j);  ddimf(j) = tmpi
         tmpi = stepf(i);  stepf(i) = stepf(j);  stepf(j) = tmpi
         tmpi = ccomp(i);  ccomp(i) = ccomp(j);  ccomp(j) = tmpi
      end do
      ! verify that we have the correct number of perturbed densities
      nd = product(1+ddimf(size(p)+1:size(dimf)))
      if (size(D) /= nd) call quit('prop_twoint error: Number of' &
               // 'perturbed densities D does not correspond to dimf (and perm)')
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_twoint error: Undefined matrix in argument D(:)')
         return
      end if
      ! permute perturbed Fock matrices in F(:) over to Ftmp(:)
      idxf = 0
      do j = 1, product(dimf)
         i = 1 + sum(idxf * stepf)
         Ftmp(j:j) = F(i:i) !pointer copy (no data copy)
         call init_mat(F(i), reset=.true.)
         do k = 1, size(dimf)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimf(k)) exit
            idxf(k) = 0
         end do
      end do
      ! everything set up, so call core procedure twoint
      call twoint(size(p), size(dimf)-size(p), pp, ccomp, ddimf, D, Ftmp)
      ! 'un-permute' perturbed Fock matrices from Ftmp(:) back into F(:)
      idxf = 0
      do j = 1, product(dimf)
         i = 1 + sum(idxf * stepf)
         F(i:i) = Ftmp(j:j) !pointer copy (no data copy)
         call init_mat(Ftmp(j), reset=.true.)
         do k = 1, size(dimf)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimf(k)) exit
            idxf(k) = 0
         end do
      end do
   end subroutine


   subroutine twoint(np, nd, p, c, df, D, F)
      integer,       intent(in) :: np, nd    !number of perturbations and order of density
      character(*),  intent(in) :: p(np)     !perturbation lables
      integer,       intent(in) :: c(np)     !lowest component of each perturbation
      integer,       intent(in) :: df(np+nd) !dimensions of perturbed Fock matrix F
      type(matrix),  intent(in) :: D(*)      !un-/perturbed density matrices (expension),
                                             !size(D) = product(1+df(np+1:np+nd))
      type(matrix),  intent(inout) :: F(*)   !where to ADD property contributions
                                             !(works incrementally), size(F) = product(df)
      integer      :: i, j, k, l, ii, jj, kk, ll, pd, pd1
      type(matrix) :: A(6) !scratch matrices
      A(:) = (/(0d0*D(1), i=1,size(A))/) !scratch matrices
      pd  = product(df(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+df(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call quit('prop_twoint error: Unperturbed ' &
                           // 'two-electron Fock matrix requested')
         ! di_GET_GbDs and di_get_sigma expects A initialized (allocated)
         call init_mat(A(1), A(1))
         do i = 0, pd-1
            if (iszero(D(pd1-pd+1+i))) cycle
            ! Coulomb-exchange
            call twofck('  ', D(pd1-pd+1+i), A(1:1))
            F(i+1) = F(i+1) + A(1)
            ! Kohn-Sham exchange-correlation
            call twofck_ks(1, (/D(1),D(pd1-pd+1+i)/), F(1+i))
         end do
         if (nd==0 .or. nd==1 .or. .not.do_dft()) then
            ! nothing more
         else if (nd==2) then
            do j = 0, df(2)-1
               do i = 0, df(1)-1
                  call twofck_ks(2, (/D(1),D(2+i),D(2+df(1)+j)/), F(1+i+df(1)*j))
               end do
            end do
         else
#ifdef PRG_DIRAC
            call twofck_ks(pd1 - 1, D, F(1))
#else
            call quit('prop_twoint: nd > 2 not implemented with DFT')
#endif
         end if
         A(1) = 0 !free
      else if (np==1 .and. p(1)=='MAG') then
         do j = 0, pd-1
            if (iszero(D(pd1-pd+1+j))) cycle
            ! Coulomb-exchange
            call twofck('M ', D(pd1-pd+1+j), A(1:3))
            do i = 0, df(1)-1
               F(1+i+df(1)*j) = F(1+i+df(1)*j) + A(c(1)+i)
            end do
            ! Kohn-Sham
            if (do_dft()) then
               if (nd==0) call di_get_MagDeriv_FxD_DFT(A(1:3), D(1))
               if (nd/=0) call di_get_MagDeriv_GxD_DFT(D(1), D(pd1-pd+1+j), A(1:3))
               do i = 0, df(1)-1
                  if (nd==0) F(1+i+df(1)*j) = F(1+i+df(1)*j) + A(c(1)+i)
                  if (nd/=0) F(1+i+df(1)*j) = F(1+i+df(1)*j) - A(c(1)+i) !negative sign
               end do
            end if
         end do
         if (nd==0 .or. nd==1 .or. .not.do_dft()) then
            ! nothing more
         else
            call quit('prop_twoint: MAG, nd > 2 not implemented with DFT')
         end if
      else if (np==1 .and. p(1)=='GEO') then
         do j = 0, pd-1 !highest-order D indices
            if (iszero(D(pd1-pd+1+j))) cycle
            do i = 0, df(1)-1 !GEO indices
               ! Coulomb-exchange
               if (i==0 .or. mod(c(1)+i-1,3)==0) then
                  call twofck('G ', D(pd1-pd+1+j:pd1-pd+1+j), A(1:3), &
                              a=1+(c(1)+i-1)/3)
                  if (do_dft()) then
                     call quit('prop_twoint: GEO not implemented with DFT')
                  end if
               end if
               F(1+i+df(1)*j) = F(1+i+df(1)*j) + A(1+mod(c(1)+i-1,3))
            end do
         end do
         if (nd==0 .or. nd==1 .or. .not.do_dft()) then
            ! nothing more
         else
            call quit('prop_twoint: GEO not implemented with DFT')
         end if
         ! if (do_dft()) print* !after all the "...integrated to nn electrons..." prints
      else if (np==2 .and. all(p==(/'MAG','MAG'/))) then
         do k = 0, pd-1
            if (iszero(D(pd1-pd+1+k))) cycle
            ! Coulomb-exchange
            call twofck('MM', D(1:1), A(1:6))
            ! Kohn-Sham
            if (do_dft()) then
               call quit('prop_twoint: MAG MAG not implemented with DFT')
            end if
            do j = 0, df(2)-1
               do i = 0, df(1)-1
                  ii = min(c(1)+i,c(2)+j) &
                     + max(c(1)+i,c(2)+j) * (max(c(1)+i,c(2)+j)-1) / 2
                  ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
                  F(1+i+df(1)*(j+df(2)*k)) = F(1+i+df(1)*(j+df(2)*k)) - A(ii)
               end do
            end do
         end do
         if (nd==0 .or. nd==1 .or. .not.do_dft()) then
             !nothing more
         else
            call quit('prop_twoint: MAG MAG and nd > 2 not implemented with DFT')
         end if
      else
         print *,'prop_twoint: No integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call quit('prop_twoint: No such integrals')
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine


   subroutine init_this_module(SorD)
   ! Re-initialize the module's data on nuclei, basis and functional,
   ! which will be passed on to the integral programs
      type(matrix) :: SorD
      call get_natoms(natom)
      pert_table(idx('GEO'))%ncomp = 3*natom
      !AllPrt(idx('NS'))%rank = depends on which isotopes
   end subroutine

#ifdef GEN1INT_DALTON
   !> \brief reads one-electron integral matrix from GEN1INT
   !> \author Bin Gao
   !> \date 2009-12-09
   !> \param lab contains the 8 characters label of the integral matrix
   !> \param S0 contains the overlap matrix
   !> \param centers contains the indices of differentiated centers (max. 4)
   !> \param orders contains the orders of geometric derivatives on each differentiated center
   !> \return A contains the integral matrix
   subroutine load_gen1int( lab, A, S0, centers, orders )
      use gen1int
      character*(*),          intent(in)    :: lab
      type(matrix),           intent(inout) :: A
      type(matrix), optional, intent(in)    :: S0
      !> \todo integer, optional, intent(in)         :: num_D
      !> \todo type(matrix), optional, intent(in)    :: D(*)
      !> \todo real(8), optional, intent(out)        :: expval(:)
      integer, optional, intent(in) :: centers(:)
      integer, optional, intent(in) :: orders(:,:)
      character*8 lab_g1int
      integer i, n
      ! if S0 present, initialize to that
      if (present(S0) .and. .not.isdef(A)) call init_mat(A, S0)
      ! if A is zero, get it allocated
      if (iszero(A)) call init_mat(A, A)
      !> \todo the following will be changed after we finish the new version of GEN1INT
      n = len_trim(lab)
      do i = 1, n
        if ( lab(i:i) /= ' ' ) then
          lab_g1int(i:i) = lab(i:i)
        else
          lab_g1int(i:i) = '_'
        end if
      end do
      lab_g1int(n+1:8) = '_'
      call gen1int_integral_fread( label_int1 = lab_g1int, &
                                   vals_g1int = A%elms,    &
                                   atoms_gdcent = centers, &
                                   order_gdcent = orders )
   end subroutine load_gen1int
#endif

   subroutine load_oneint(lab, A, S0)
!radovan: please don't use fixed size, size can vary in DIRAC
!     character*8,            intent(in)    :: lab
      character(*),           intent(in)    :: lab
      type(matrix),           intent(inout) :: A
      type(matrix), optional, intent(in)    :: S0
      ! if S0 present, initialize to that
      if (present(S0) .and. .not.isdef(A)) call init_mat(A, S0)
      ! if A is zero, get it allocated
      if (iszero(A)) call init_mat(A, A)
      ! read integrals from file into A
      call di_read_operator_int(lab, A)
   end subroutine


   subroutine twofck(what, D, F, a, b)
   ! Add (un)perturbed 2e-contributions to Fock matrices
   ! what='  ' -> unperturbed
   ! what='G ' -> 3 geometry-perturbed wrt. nucleus a
   ! what='M ' -> 3 magnetic-field-perturbed (London)
   ! what='GG' -> 9 2nd-order geometry-perturbed wrt. nuclei a, b
   ! what='MM' -> 6 2nd-order magnetic-field-perturbed (London)
   ! what='GM' -> 9 2nd-order mixed geometry and magnetic
      character*2,       intent(in)    :: what
      type(matrix),      intent(in)    :: D(1)
      type(matrix),      intent(inout) :: F(:)
      integer, optional, intent(in)    :: a, b
      type(matrix)     :: dupD, dupF(size(F))
      integer          :: nf, n2, lwrk, i, aa=0
      real(8), pointer :: wrk(:)
      nf = size(F)
      n2 = D(1)%nrow**2
      if (what=='  ' .and. nf==1) then
         lwrk = 0
      else if (what=='M ' .and. nf==3) then
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else if (what=='G ' .and. nf==3) then
         if (.not.present(a)) &
            call quit("prop_integrals twofck: Geometry Fock but atom 'a' not specified")
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
         aa = a !may not be present, replicate
      else if (what=='MM' .and. nf==6) then
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else if (what=='GG' .and. nf==3*3) then
         if (.not.present(a) .or. .not.present(b)) &
            call quit("prop_integrals twock: 2nd ord. geom Fock reqires atoms 'a','b' specified")
         !lwrk = (1+2*nf+50)*n2 + 10000*nbas + 5000000
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else
         call quit('prop_integrals twofck: wrong size(F) for what=' // what)
      end if
      !make usre F is properly allocated
      do i = 1, nf
         if (iszero(F(i)) .or. (D(1)%complex .and. .not.F(i)%complex)) then
            if (.not.iszero(F(i))) F(i) = 0 !deallocate
            call init_mat(F(i), F(i), complex=D(1)%complex)
         end if
      end do
      !allocate work
      if (lwrk/=0) call di_select_wrk(wrk, lwrk)
      !process complex D with two calls
      if (D(1)%complex) then
         call init_mat(dupD, D(1), alias='RF')
         do i = 1, nf
            call init_mat(dupF(i), F(i), alias='RF')
         end do
         call subr(dupD, dupF)
         call init_mat(dupD, D(1), alias='IF')
         do i = 1, nf
            call init_mat(dupF(i), F(i), alias='IF')
         end do
         call subr(dupD, dupF)
      else
         call subr(D(1), F)
      end if
      if (lwrk/=0) call di_deselect_wrk(wrk, lwrk)
   contains
      subroutine subr(D, F)
         type(matrix), intent(in)    :: D
         type(matrix), intent(inout) :: F(nf)
         if (what=='  ') then
            call di_GET_GbDs(D, F(1))
            return
         end if
         wrk(1:n2) = D%elms !ajt fixme
         call GRCONT(wrk( 1+n2+n2*nf : lwrk ), (lwrk-n2-n2*nf),         &
                     wrk( 1+n2 : n2+n2*nf ), n2*nf, (what(1:1) == 'G'), &
                     (what(1:1) == 'M'), merge(1,2,what(2:2)==' '),     &
                     aa, .false., .true., wrk(1:n2), 1)
         do i = 1, nf
            F(i)%elms = wrk( 1+n2*i : n2*(1+i) ) !ajt fixme
         end do
      end subroutine
   end subroutine


   subroutine twofck_ks(n, D, F)
      integer,      intent(in)    :: n
#ifdef PRG_DIRAC
      type(matrix)                :: D(n+1)
#else
      type(matrix), intent(in)    :: D(n+1)
#endif
      type(matrix), intent(inout) :: F
      type(matrix) :: A
      integer      :: i
#ifdef PRG_DIRAC
      integer :: ones(200), nr_dmat, nz, order
#endif

      if (.not. do_dft()) then
!        we want no xc
         return
      end if

#ifndef PRG_DIRAC
!radovan: this is not correct for CR and higher
!         there you can have nonzero integral with one or more
!         density matrices zero
      if (any((/(iszero(D(i)), i = 2, n+1)/))) then
!        at least one of the density matrices is zero
         return
      end if
#endif

      call init_mat(A, D(1))

      if (n==0) then
         call quit('prop_contribs/twofck_ks error: Kohn-Sham contribution ' &
                // 'to unperturbed Fock matrix requested, but not implemented')
#ifndef PRG_DIRAC
      else if (n==1) then
         if (D(2)%complex) &
            call di_get_sigma_xc_cont(re(D(2)), D(1), D(1), A)
         if (.not.D(2)%complex) &
            call di_get_sigma_xc_cont(D(2), D(1), D(1), A)
         F = F + A
         if (D(2)%complex) &
            call di_get_sigma_xc_cont(im(D(2)), D(1), D(1), A)
         if (D(2)%complex) F = F + (0d0,1d0)*A
      else if (n==2) then
         if (D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), re(D(2)), re(D(3)), A)
         if (D(2)%complex .and. .not.D(3)%complex) &
            call di_get_T_xc_cont(D(1), re(D(2)), D(3), A)
         if (.not.D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), D(2), re(D(3)), A)
         if (.not.D(2)%complex .and. .not.D(3)%complex) &
            call di_get_T_xc_cont(D(1), D(2), D(3), A)
         F = F + A
         if (D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), re(D(2)), im(D(3)), A)
         if (.not.D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), D(2), im(D(3)), A)
         if (D(3)%complex) F = F + (0d0,1d0)*A
         if (D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), im(D(2)), re(D(3)), A)
         if (D(2)%complex .and. .not.D(3)%complex) &
            call di_get_T_xc_cont(D(1), im(D(2)), D(3), A)
         if (D(2)%complex) F = F + (0d0,1d0)*A
         if (D(2)%complex .and. D(3)%complex) &
            call di_get_T_xc_cont(D(1), im(D(2)), im(D(3)), A)
         if (D(2)%complex .and. D(3)%complex) F = F - A
      else
         call quit('prop_contribs/twofck_ks error: n > 2 not implemented')
      end if
#else
      else

         ones    = 1
         nz      = size(A%elms)/(A%nrow*A%ncol)
         nr_dmat = size(D)
         do i = 2, nr_dmat
            if (iszero(D(i))) then
               D(i) = tiny(0.0d0)*D(1)
            end if
         end do

         A%elms = 0.0d0
         A%irep = 0
         do i = 2, nr_dmat
            A%irep = ieor(A%irep, D(i)%irep)
         end do

!        that's really ugly here - will make it better later
!        for the moment i'm just happy that it works
         order = 0
         if (nr_dmat > (0)) then
            order = 1
         end if
         if (nr_dmat > (1)) then
            order = 2
         end if
         if (nr_dmat > (1+1)) then
            order = 3
         end if
         if (nr_dmat > (1+2+1)) then
            order = 4
         end if
         if (nr_dmat > (1+3+3+1)) then
            order = 5
         end if
         if (nr_dmat > (1+4+6+4+1)) then
            order = 6
         end if
         if (nr_dmat > (1+5+10+10+5+1)) then
            call quit('fix twofck_ks order pascal triangle')
         end if
         
         if (order < 2) then
            return
         end if

!radovan:
!        in many situations one could pass F directly
!        to dftexc but the gga symmetrization can make problems
!        with hybrid functionals and anti-hermitian contributions (nonzero frequencies)
!        so in other words while A is always symmetric, F is non necessarily so
!        and it would get symmetrized inside dftexc
         call dftexc(A%nrow,                              &
                     nz,                                  &
                     1,                                   &
                     A%elms,                              &
                     D(1)%elms,                           &
                     nr_dmat - 1,                         &
                     (/(D(i)%elms, i = 2, nr_dmat)/),     &
                     (/(D(i)%irep + 1, i = 2, nr_dmat)/), &
                     (/A%irep + 1/),                      &
                     (/ones(1:nr_dmat - 1)/),             &
                     0,                                   &
                     .false.,                             &
                     order,                               &
                     .false.,                             &
                     .false.)

         F = F + A

      end if
#endif /* ifndef PRG_DIRAC */

      A = 0

   end subroutine


   subroutine twoctr(what, Da, Db, E)
   !gm=T -> 3*NATOMS geometriske
   !gm=F -> tre magnetiske
      character(*), intent(in)    :: what
      type(matrix), intent(in)    :: Da, Db
      real(8),      intent(inout) :: E(:)
      real(8), pointer :: wrk(:)
      integer          :: lwrk

      if (natom==0) call quit('prop_contribs/twoctr error: Module not initialized')
      !ajt MagSus contraction doesn't work, so I disabled this. Use twofck instead
      if (what=='M' .or. what=='MM') &
         call quit("prop_contribs/twoctr: 'M' or 'MM' not yet implemented")
      if     (what=='G' .and. size(e) == 3*natom) then
         lwrk = 50*Da%nrow**2 + 10000*Da%nrow + 5000000
      else if (what=='M' .and. size(e)==3) then
         lwrk = 50*Da%nrow**2 + 10000*Da%nrow + 5000000
      else if (what=='GG' .and. size(e) == 3*3*natom**2) then
         lwrk = 50*Da%nrow**2 + 10000*Da%nrow + 10000000
      else if (what=='MM' .and. size(e)==3*3) then
         lwrk = 50*Da%nrow**2 + 10000*Da%nrow + 5000000
      else
          call quit('prop_contribs/twoctr: wrong size(E) for what=' // what)
      end if

#ifdef PRG_DIRAC
!radovan: for DIRAC the above lwrk allocation is too much (and more than needed)
!         so rather use the whole work array and hope for the best
      lwrk = len_f77_work
#endif

      call di_select_wrk(wrk, lwrk)

#ifndef PRG_DIRAC
!radovan: is this necessary?
      wrk(             1:  Da%nrow**2) = Da%elms !call mat_to_full(D ,1d0,A(:,:,1))
      wrk(Da%nrow**2 + 1:2*Da%nrow**2) = Db%elms !call mat_to_full(DD,1d0,A(:,:,2))
#endif

      E = 0.0d0

#ifdef PRG_DIRAC
      call grcont(wrk,                      &
                  lwrk,                     &
#else
!     this does not compile in DIRAC
!     besides matrices have larger size in DIRAC: 4*Da%nrow**2
      call grcont(wrk(2*Da%nrow**2 + 1),    &
                  lwrk - 2*Da%nrow**2,      &
#endif
                  E,                        &
                  size(E),                  &
                  what(1:1) == 'G',         &
                  what(1:1) == 'M',         &
                  len(what),                &
                  0,                        &
                  .true.,                   &
                  .false.,                  &
#ifdef PRG_DIRAC
                  (/Da%elms, Db%elms/),     &
#else
                  wrk(1),                   &
#endif
                  2)

#ifdef PRG_DIRAC
!radovan: there is a factor 8 missing in DIRAC
!         will fix it later so that it's taken care of inside grcont
      E = 8.0d0*E
#endif

      call di_deselect_wrk(wrk, lwrk)

   end subroutine


   subroutine oneint_ave(what, D, DFD, R)
      character(*), intent(in)  :: what
      type(matrix), intent(in)  :: D, DFD
      real(8),      intent(out) :: R(:)
#include <mxcent.h>
#include <taymol.h>
      real(8), pointer :: wrk(:)
      integer          :: lwrk
      call save_D_and_DFD_for_ABACUS(.false., D, DFD)
      lwrk = 50*D%nrow**2+10000*D%nrow+50000000
      call di_select_wrk(wrk, lwrk)
      HESMOL(:3*natom,:3*natom) = 0
      ! SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
      ! &                  NODV,DIFDIP,HFONLY,NCLONE)
      call ONEDRV(wrk,lwrk,0,.true.,len(what),.true.,.true., &
                  .true.,.false.,.true.,.false.)
      ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
      R(1:9*natom**2) = reshape(HESMOL(:3*natom,:3*natom), (/9*natom**2/))
      call di_deselect_wrk(wrk, lwrk)
   end subroutine


   function idx(p)
      character(*) :: p
      integer      :: idx, i, j
      do idx = 1, size(pert_table)
          if (pert_table(idx)%code == p) return
      end do
      if (p(1:min(3,len(p))) == 'AUX') then
         read (p(min(4,len(p)):),'(i1)', iostat=j) i
         if (j /= 0 .or. i < 0 .or. i > 9) &
            call quit('prop_contribs error: ' &
                     // 'Auxiliary label index should be within 0..9, ' // p)
         idx = 1
         return
      end if
      call quit('Perturbation not found: ' // p)
   end function


   function pert_antisym(p)
      character(*), intent(in) :: p(:)
      logical :: pert_antisym(size(p))
      integer :: i
      pert_antisym = (/(pert_table(idx(p(i)))%anti, i=1,size(p))/)
   end function


   function pert_shape(p)
      character(*), intent(in) :: p(:)
      integer :: pert_shape(size(p)), i
      pert_shape = (/(pert_table(idx(p(i)))%ncomp, i=1,size(p))/)
   end function


   function pert_basdep(p)
      character(*), intent(in) :: p(:)
      logical :: pert_basdep(size(p))
      integer :: i
      pert_basdep = (/(pert_table(idx(p(i)))%bas, i=1,size(p))/)
   end function


   subroutine save_D_and_DFD_for_ABACUS(anti, D, DFD)
   !ajt aug09 Same as set_dsofso in fock-eval.f90, but pertrubed
   !          D and DFD is input instead of unperturbed D and F
#if defined(DALTON_AO_RSP)
      use xmatrix, only: mat_to_full => xmatrix_get_full
#else
      use matrix_operations, only: mat_to_full
#endif
      logical,      intent(in) :: anti !whether the integrals D and DFD
         ! are to be contracted with are symmetric or anti-symmetric
      type(matrix), intent(in) :: D, DFD
         ! perturbed (or un-) density and
         ! 'generalized Fock matrix' (energy-weighted density matrix)
      real(8) :: Full(D%nrow,D%nrow), &
                 Dtri(D%nrow*(D%nrow+1)/2), &
                 DFDtri(D%nrow*(D%nrow+1)/2)
      if (iszero(D)) then
         Dtri = 0
      else
         call mat_to_full(D, merge(4d0,2d0,anti), Full, 'D')
         if (.not.anti) call DGEFSP(D%nrow, Full, Dtri)
         if (     anti) call DGETAP(D%nrow, Full, Dtri)
      end if
      if (iszero(DFD)) then
         DFDtri = 0
      else
         call mat_to_full(DFD, merge(4d0,2d0,anti), Full, 'DFD')
         if (.not.anti) call DGEFSP(D%nrow, Full, DFDtri)
         if (     anti) call DGETAP(D%nrow, Full, DFDtri)
      end if
      ! write fo files
#ifdef PRG_DIRAC
      call quit('dirac: implement write_dsofso')
#else
      call write_dsofso(Dtri,DFDtri)
#endif
   end subroutine


end module
