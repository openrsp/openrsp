! Copyright 2009 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module prop_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices.
module prop_contribs_old

   use matrix_defop
   use interface_molecule
   use interface_xc
   use interface_1el
   use interface_f77_memory
   use interface_scf
   use dalton_ifc
   use nuc_contributions
   use interface_io

   implicit none

   public prop_oneave
   public prop_twoave
   public prop_oneint
   public prop_twoint
   public pert_basdep
   public pert_shape
   public pert_antisym
   public prop_auxlab
   public prop_field

   private


   !> Type describing a single field in a response function
   !> or response equation. A response equation (or density)
   !> corresponds to an array of prop_field. Similarly
   !> a response function corresponds to an array of prop_field
   !> whose freqs sum to zero.
   type prop_field
      !> 4-char pert label
      character(4) :: label
      !> frequency
      complex(8)   :: freq
      !> first component
      integer      :: comp
      !> number of components
      integer      :: ncomp
   end type

   !> private struct to collect
   type prop_field_info
      !> four-letter abbreviation
      character(4)  :: code
      !> long name
      character(64) :: name
      !> number of components (when known, 0 otherwise)
      integer       :: ncomp
      !> anti-symmetric (1,3,5th ord.) perturbed integrals
      logical       :: anti
      !> basis dependent (sa. GEO and MAG)
      logical       :: bas
      !> one-electron operator linear in field strength (EL)
      logical       :: lin
      !> one-electron operator quadratic in field strength (MAGO)
      logical       :: quad
   end type


   ! to compactify the table
   logical, parameter :: T = .true.
   logical, parameter :: F = .false.


   !> ajt nov09: AUX0..AUX9 are 10 configurable basis-independent 1-electron
   !>            perturbations, configured by setting the corresponding
   !>            HERMIT integral label in prop_auxlab(0:9).
   !> ajt jan10: EXCI is a ZERO (no) perturbation, and is introduced to
   !>            allow the same code to contract response functions and
   !>            "generalized transition moments".
   !> ajt may10: FREQ is also a ZERO (no) perturbation, and is introduced to
   !>            allow the same code to contract response functions and
   !>            frequency-differentiated response functions.
   type(prop_field_info) :: field_list(13) = &                         !nc an ba ln qu
      (/prop_field_info('EXCI', 'Generalized "excitation" field'      , 1, F, F, T, T), &
        prop_field_info('FREQ', 'Generalized "freqency" field'        , 1, F, F, T, T), &
        prop_field_info('AUX*', 'Auxiliary integrals on file'         , 1, F, F, T, F), &
        prop_field_info('EL'  , 'Electric field'                      , 3, F, F, T, F), &
        prop_field_info('VEL' , 'Velocity'                            , 3, T, F, T, F), &
        prop_field_info('MAGO', 'Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
        prop_field_info('MAG' , 'Magnetic field with London orbitals' , 3, T, T, F, F), &
        prop_field_info('ELGR', 'Electric field gradient'             , 6, F, F, T, F), &
        prop_field_info('VIBM', 'Displacement along vibrational modes',-1, F, T, F, F), &
        prop_field_info('GEO' , 'Nuclear coordinates'                 ,-1, F, T, F, F), & !-1=mol-dep
        prop_field_info('NUCM', 'Nuclear magnetic moment'             ,-1, F, T, F, T), & !-1=mol-dep
        prop_field_info('AOCC', 'AO contraction coefficients'         ,-1, F, T, F, F), & !-1=mol-dep
        prop_field_info('AOEX', 'AO exponents'                        ,-1, F, T, F, F)/)  !-1=mol-dep

   character(8) :: prop_auxlab(0:9)

   character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

   logical, external :: do_dft

contains


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
                           // 'in the specified number of ASCII caracters',-1)
   end function


   !> Contracts the 1-electron integrals perturbed by the perturbations p(:)
   !> with the perturbed density matrices in D(:) (e.g. D=(/Dxy/) for a 2nd
   !> order density), and ADDS the result to the property (rsp func) array
   !> E(:). Front for the private subroutine 'oneave' below, checking the
   !> arguments' dimensions, and doing permutations
   !> S0 is passed as argument only as reference to nuclei and basis set
   subroutine prop_oneave(S0, p, D, dime, E, perm, comp, freq, DFD)
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to
      !> contract perturbed one-electron integrals with.
      !> If perm present, size(D) = product(dime(perm(np+1:np+nd))),
      !> if perm not present, size(D) = product(dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(E), dimensions of property
      integer,           intent(in) :: dime(:)
      !-------------------------------------------------------------
      !> property contributions, works incrementally, thus
      !> contributions are ADDED to E(*). size(E) = product(dime)
      complex(8),     intent(inout) :: E(*)
      !-------------------------------------------------------------
      !> perm(np+nd), permutation of indices.
      !> For each dimension of p and D, the corresponding dimension in E.
      !> Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p.
      !> Default 1 1 ... 1
      integer,      optional, intent(in) :: comp(:)
      !> freq(np), complex frequencies
      !> for each p, default all zero. Multiply the half-derivative
      !> overlap integrals, thus no contribution if basis independent of p
      complex(8),   optional, intent(in) :: freq(:)
      !> optional perturbed energy-weighted density matrices.
      !> Contracted against perturbed overlap integrals
      type(matrix), optional, intent(in) :: DFD(:)
      !-------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                      stepe(size(dime)), ddime(size(dime)), idxe(size(dime)), &
                      i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: ffreq(size(p)), Etmp(product(dime)), tmpf
      logical      :: zero, bas
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.field_list(idxp(i))%lin .and. &
             .not.field_list(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (field_list(idxp(i))%lin .and. &
             field_list(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (field_list(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (field_list(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < size(p)) call quit('prop_oneave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < size(p)',get_print_unit())
      if (any(dime <= 0)) call quit('prop_oneave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0',get_print_unit())
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call quit('prop_oneave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)',get_print_unit())
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call quit('prop_oneave ' &
                      // 'argument error: Permutation must contain each number exactly once',get_print_unit())
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_oneave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',get_print_unit())
         if (any(comp <= 0)) call quit('prop_oneave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0',get_print_unit())
         ccomp = comp
      end if
      if (any(ccomp + ddime(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_oneave argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(p)',get_print_unit())
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call quit('prop_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)',get_print_unit())
         ffreq = freq
      end if
      ! if unperturbed density and anti-symmetric integral, also zero
      j = count((/(field_list(idxp(i))%anti, i=1,size(p))/))
      if (size(p)==size(dime) .and. mod(j,2)==1) then
         ! if not all bas-dep, or all bas-dep and zero frequencies
         if (.not.all((/(field_list(idxp(i))%bas,i=1,size(p))/)) &
             .or. all(ffreq==0)) zero = .true.
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
               // 'perturbed densities D does not correspond to dime (and perm)',get_print_unit())
      ! verify number of DFD, and that all are defined
      bas = all((/(field_list(idxp(i))%bas, i=1,size(p))/))
      if (present(DFD)) then
         if (size(DFD) /= nd) call quit('prop_oneave error: Number of' &
               // 'perturbed DFD differs from number of perturbed densities D',get_print_unit())
         ! if no basis perturbation (or zero perturbed integrals,
         ! DFD will not be used, so verify they are defined
         if (.not.bas .or. zero) then
            if (.not.all((/(isdef(DFD(i)), i=1,nd)/))) &
               call quit('prop_oneave error: Undefined matrix in argument DFD(:)',get_print_unit())
         end if
      end if
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_oneave error: Undefined matrix in argument D(:)',get_print_unit())
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
      !> unperturbed overlap, to know its dimension
      type(matrix),      intent(in)  :: S0
      !> number of perturbations and order of density
      integer,           intent(in)  :: np
      !> perturbation lables
      character(*),      intent(in)  :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in)  :: c(np)
      !> dimensions of property integrals
      integer,           intent(in)  :: dp(np)
      !> frequency of each p
      complex(8),        intent(in)  :: w(np)
      !> dimensions of property integrals
      integer,           intent(in)  :: nd
      !> un-/perturbed density matrices,
      !> size(D) = product(1+de(np+1:np+nd))
      type(matrix),      intent(in)  :: D(max(1,nd))
      !-----------------------------------------------------
      !> resulting one-electron property contributions,
      !> size(E) = product(dp) * nd
      complex(8),        intent(out) :: E(*)
      !> un-/perturbed energy-weighted density matrices,
      !> Should have size(DFD) = size(D)
      type(matrix), optional, intent(in) :: DFD(:)
      !-----------------------------------------------------
      real(8), allocatable :: RR(:) !scratch
      real(8)      :: R(6) !scratch
      integer      :: i, j, k, l, ii, jj, kk, ll, na
      type(matrix) :: A(6) !scratch matrices
      na = get_nr_atoms()
      do i = 1, size(A)
         A(i) = 0*S0 !scratch matrices
      end do
      if (np==0) then
         call quit('prop_oneave error: unperturbed one-electron contribution requested',get_print_unit())
      else if (np==1 .and. p(1)(1:3)=='AUX') then
         read (p(1)(4:),'(i1)') i
         if (i < 0 .or. i > 9) call quit('prop_oneave error: Index in' &
                  // ' auxiliary label out of range 0..9: ' // p(1),get_print_unit())
         call load_oneint(prop_auxlab(i), A(1))
         do j=0, max(0,nd-1)
            E(1+j) = tr(A(1),D(1+j))
         end do
      else if (np==1 .and. p(1)=='EL') then
         ! contract -dipole integrals '?DIPLEN ' with densities
         do i=0, dp(1)-1
               ! load -dipole integral from file
               call load_oneint(xyz(c(1)+i) // 'DIPLEN ',A(1))
               ! loop over density matrices
               do j=0, max(0,nd-1)
                  E(1+i+dp(1)*j) = tr(A(1),D(j+1)) !*2 for total dens hidden in tr
               end do
         end do
         ! no densities means D(1) contains unperturbed density matrix.
         ! Then also add nuclear attraction contribution to -dipole moment
         if (nd==0) then
            call DIPNUC_ifc(R(:3))
            E(:dp(1)) = E(:dp(1)) - R(c(1):c(1)+dp(1)-1) !sign since dE/d(EL) = -dipole
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
               E(1+i+dp(1)*j) = tr(A(1), D(j+1)) / 2 !factor 1/2
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
               A(1) = A(1) - w(1)/2 * trps(A(2))
               A(2) = A(2) - trps(A(2))
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
            call load_oneint('1DHAM' // prefix_zeros(c(1)+i, 3), A(1))
            ! (half-) perturbed overlap -i/2 Tg added to A(1), Sg in A(2)
            if (w(1)==0) then !w=0 means no -i/2 Tg contribution
               call load_oneint('1DOVL' // prefix_zeros(c(1)+i,3), A(2))
               A(2) = -A(2) !1DOVL is really -dS/dg
            else
               call load_oneint('SQHDR' // prefix_zeros(c(1)+i,3), A(2))
               A(2) = -A(2) !SQHDR is really -dS>/dg
               A(1) = A(1) - w(1)/2 * A(2)
               A(1) = A(1) + w(1)/2 * trps(A(2))
               A(2) = A(2) + trps(A(2))
            end if
            do j = 0, max(0,nd-1)
               E(1+i+dp(1)*j) = tr(A(1), D(j+1))
               if (present(DFD)) &
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) - tr(A(2), DFD(j+1))
            end do
         end do
         ! nuclear repulsion contribution
         if (nd==0) then
            allocate(RR(3*na))
            call gradient_nuc(na, RR(:3*na))
            E(:dp(1)) = E(:dp(1)) + RR(c(1):c(1)+dp(1)-1)
            deallocate(RR)
         end if
      else if (np==2 .and. all(p==(/'MAGO','MAGO'/))) then
         do j = 0, dp(2) - 1
            do i = 0, dp(1) - 1
               call load_oneint(xyz(min(c(1)+i,c(2)+j)) &
                        // xyz(max(c(1)+i,c(2)+j)) // 'SUSCGO', A(1))
               do k = 0, max(0,nd-1)
                  ! ajt Since MAGO is imaginary/anti-symmetri, we get a minus sign here
                  E(1+i+dp(1)*(j+dp(2)*k)) = -tr(A(1),D(1+k))
               end do
            end do
         end do
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
                  ! only load >> XX XY YY XZ YZ ZZ, use trps(..) for <<
                  call load_oneint('>>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                       // xyz(max(c(1)+i,c(2)+j)), A(2))
                  A(1) = A(1) - (w(1)+w(2))/2 * A(2)
                  A(1) = A(1) + (w(1)+w(2))/2 * trps(A(2))
                  A(2) = A(2) + trps(A(2))
                  ! only load <> XX XY YY XZ YZ ZZ, use trps(..) for the others
                  call load_oneint('<>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                       // xyz(max(c(1)+i,c(2)+j)), A(3))
                  !ajt1009 Changed sign on these two: (think it's correct)
                  A(1) = A(1) - merge(-1, 1, c(1)+i < c(2)+j) &
                                       * (w(1)-w(2))/2 * A(3)
                  A(1) = A(1) + merge(-1, 1, c(1)+i < c(2)+j) &
                                  * (w(1)-w(2))/2 * trps(A(3))
                  A(2) = A(2) + A(3)
                  A(2) = A(2) + trps(A(3))
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
               call load_oneint(prefix_zeros(c(1)+i,3) // 'DPG ' &
                           // xyz(c(2)+j), A(1))
               do k = 0, max(0,nd-1) !density indices
                  E(1+i+dp(1)*(j+dp(2)*k)) = -tr(A(1), D(k+1)) !minus sign
               end do
            end do
         end do
         ! nd==0 has nuclear repulsion contribution to -dipole gradient
         if (nd==0) then
            allocate(RR((3*na)*3))
            call DPGNUC_ifc(na, RR(:(3*na)*3))
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) &
                                 - RR(c(1)+i + (3*na)*(c(2)+j-1)) !-dipole sign
               end do
            end do
            deallocate(RR)
         end if
      ! GEO MAGO has no nd==0 because MAGO anti-symmetric
      else if (np==2 .and. all(p==(/'GEO ','MAGO'/)) .and. nd/=0) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               call load_oneint(prefix_zeros(c(1)+i,3) // 'AMDR' &
                           // xyz(c(2)+j), A(1))
               do k = 0, nd-1
                  E(1+i+dp(1)*(j+dp(2)*k)) = tr(A(1), D(k+1)) / 2 !factor 1/2
               end do
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO','MAG'/))) then
         ! ajt This is a stub, and is only correct when we contract against
         !     the (symmetric) unperturbed density (nd==0, w1==-w2)
         if (nd /= 0) &
            call quit('prop_oneave: GEO MAG is not fully implemented. ' &
                   // 'Currently, the density must be the unperturbed one.')
         if (w(1) /= -w(2)) &
            call quit('prop_oneave: GEO MAG is not fully implemented. ' &
                   // 'Currently, frequencies must be same with opposite sign.')
         if (w(1)-w(2) /= 0) then
            ! electronic contribution
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  call load_oneint(prefix_zeros(c(1)+i,2) // ' HDB ' &
                                // xyz(c(2)+j), A(1))
                  do k = 0, max(0,nd-1)
                     E(1+i+dp(1)*(j+dp(2)*k)) = (-1) * & !factor -1
                                ((w(1)-w(2))/2 *  tr(A(1), D(k+1)) &
                              +  (w(1)-w(2))/2 * dot(A(1), D(k+1)))
                  end do
               end do
            end do
            ! nuclear contribution if unperturbed density (nd==0)
            if (nd == 0) then
               allocate(RR(3*(3*na)))
               call AATNUC_ifc(3*na, RR)
               do j = 0, dp(2)-1
                  do i = 0, dp(1)-1
                     E(1+i+dp(1)*j) = E(1+i+dp(1)*j) + 2 & !factor 2
                                    * (w(1)-w(2))/2 * RR(1+j+3*i)
                  end do
               end do
               deallocate(RR)
            end if
         end if
      else if (np==2 .and. all(p==(/'GEO ','ELGR'/))) then
         do i = 0, dp(1)-1 !GEO indices
            ! ajt The integrals '01QDG XY' are not traceless. Therefore,
            !     load all components at once and remove trace before contracting
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXX', A(1))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGYY', A(3))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGZZ', A(6))
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXY', A(2))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXZ', A(4))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGYZ', A(5))
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
         allocate(RR(3*na*3*na))
         ! one-electron integral contribution
         do k = 0, max(0,nd-1)
            !ajt fixme oneint_ave(GG..) misses frequency dependent -i/2 Tgg contribution
            if (present(DFD)) then
               call oneint_ave(get_nr_atoms(), 'GG', D(k+1), DFD(k+1), RR)
            endif
            if (.not.present(DFD)) then
               A(1) = 0*S0
               call oneint_ave(get_nr_atoms(), 'GG', D(k+1), A(1), RR)
            endif
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*(j+dp(2)*k)) = RR(c(1)+i + 3*na*(c(2)+j-1))
               end do
            end do
         end do
         ! for zero-order density, add nuclear contribution to Hessian
         if (nd==0) then
            call hessian_nuc(na, RR)
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  E(1+i+dp(1)*j) = E(1+i+dp(1)*j) + RR(c(1)+i + 3*na*(c(2)+j-1))
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
         call quit('prop_oneave: no such integrals',get_print_unit())
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine



   !> Contract the 2-electron integrals perturbed by the perturbations p(:)
   !> with the perturbed density matrix expansion in D(:) (e.g. D=(/D,Dx,Dy,Dxy/)
   !> for a 2nd order expansion), and ADDS the result to the property (rsp func)
   !> array E(:). Front for the private subroutine 'twoave' below, checking the
   !> arguments' dimensions, and doing permutations
   !> D(1) serves as reference to nuclei, basis and model/functional
   subroutine prop_twoave(p, D, dime, E, perm, comp)
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to contract perturbed one-electron
      !> integrals with. If perm present, size(D) = product(dime(perm(np+1:np+nd))),
      !> if perm not present, size(D) = product(dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(E), dimensions of property
      integer,           intent(in) :: dime(:)
      !------------------------------------------------------------------------
      !> property contributions, works incrementally, thus contributions are
      !> ADDED to E(*). size(E) = product(dime)
      complex(8),     intent(inout) :: E(*) !
      !------------------------------------------------------------------------
      !> perm(np+nd), permutation of indices. For each dimension of p and D,
      !> the corresponding dimension in E. Default 1 2 ... np+nd (no permutation)
      integer, optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer, optional, intent(in) :: comp(:)
      !------------------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dime)), ccomp(size(p)), &
                      stepe(size(dime)), ddime(size(dime)), &
                      idxe(size(dime)), i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: Etmp(product(dime))
      logical      :: zero
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero.
      zero = any((/(field_list(idxp(i))%lin  .or. &
                    field_list(idxp(i))%quad .or. &
               .not.field_list(idxp(i))%bas,  i=1,size(p))/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dime) < size(p)) call quit('prop_twoave argument error: ' &
               // 'More perturbations than dimensions of property, size(dime) < size(p)',-1)
      if (any(dime <= 0)) call quit('prop_twoave argument error: ' &
               // 'Property has a zero or negative dimension, dime <= 0',-1)
      ! compute step lengths in E (cumulative products of dimensions)
      stepe(1) = 1
      do i = 2, size(dime)
         stepe(i) = stepe(i-1)*dime(i-1)
      end do
      ! reorder dimensions and step lengths in E according to permutation argument perm
      ddime = dime
      if (present(perm)) then
         if (size(perm) /= size(dime)) call quit('prop_twoave argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dime)',-1)
         ! verify that perm is indeed a permutation
         do i = 1, size(dime)-1
            if (perm(i) <= 0 .or. perm(i) > size(dime) .or. &
                any(perm(i) == perm(i+1:size(dime)))) call quit('prop_twoave ' &
                      // 'argument error: Permutation must contain each number exactly once',-1)
         end do
         ddime = (/( dime(perm(i)), i=1,size(dime))/)
         stepe = (/(stepe(perm(i)), i=1,size(dime))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_twoave argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',-1)
         if (any(comp <= 0)) call quit('prop_twoave argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0',-1)
         ccomp = comp
      end if
      if (any(ccomp + ddime(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_twoave argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dime > pert_shape(p)',-1)
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
               // 'perturbed densities D does not correspond to dime (and perm)',-1)
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_twoave error: Undefined matrix in argument D(:)',-1)
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
      !> number of perturbations and order of density
      integer,           intent(in)  :: np, nd
      !> perturbation lables
      character(*),      intent(in)  :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in)  :: c(np)
      !> dimensions of property (E)
      integer,           intent(in)  :: de(np+nd)
      !> un-/perturbed density matrices (expension),
      !> size(D) = product(1+de(np+1:np+nd))
      type(matrix),      intent(in)  :: D(*)
      !> where to ADD property contributions
      !> (works incrementally), size(E) = product(de)
      complex(8),        intent(out) :: E(*)
      !-------------------------------------
      real(8), allocatable :: RR(:) !scratch
      real(8)      :: R(6) !scratch
      integer      :: i, j, k, l, ii, jj, kk, ll, pd, pd1, na
      type(matrix) :: A(6) !scratch matrices
      na = get_nr_atoms()
      do i = 1, size(A)
         A(i) = 0*D(1) !scratch matrices
      end do
      pd  = product(de(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+de(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call quit('prop_twoave: unperturbed energy/integrals requested',-1)
         !The highest order Ds multiply the unperturbed F, which is unknown here,
         !so insist that no highest-order Ds are present
         if (.not.all((/(iszero(D(i)), i=pd1-pd+1,pd1)/))) &
            call quit('prop_twoave: unperturbed Fock matrix requested',-1)
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
                  A(1) = 0*A(1)
               end do
            end do
            ! contract first-order densities to Fock, then trace with second-order
            !ajt Fixme
            if (.not.all((/(iszero(D(i)), i=2+de(1)+de(2)+de(3),pd1)/))) &
               call quit('prop_twoave: nd = 3 not fully implemented',-1)
         else if (nd > 3) then
            call quit('prop_twoave: nd > 3 not implemented',-1)
         end if
      !London magnetic, no nd==0 because because MAG anti
      else if (np==1 .and. p(1)=='MAG' .and. nd /= 0) then
         ! di_get_MagDeriv_(F,G)xD_DFT takes initialized
         if (do_dft()) then !scratch in A(4:6)
            do i = 4, 6
               if (iszero(A(i))) call mat_ensure_alloc(A(i))
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
         allocate(RR(3*na*2))
         ! highest-order contribution
         do j = 0, pd-1
            if (iszero(D(pd1-pd+1+j))) then
               E( 1+de(1)*j : de(1)*(j+1) ) = 0
               cycle
            end if
            ! Coulomb-exchange
            call twoctr('G', D(1), D(pd1-pd+1+j), RR(:3*na))
            ! for molgra (nd==0), factor 1/2 on these integrals
            if (nd==0) RR(:3*na) = RR(:3*na)/2
            E( 1+de(1)*j : de(1)*(j+1) ) = RR(c(1):c(1)+de(1)-1)
         end do
         if (nd==0 .or. nd==1) then
            ! nothing more
         else if (nd==2) then
            ! integrals over products of first order densities
            do k = 0, de(3)-1
               if (iszero(D(2+de(2)+k))) cycle
               do j = 0, de(2)-1
                  if (iszero(D(2+j))) cycle
                  ! Coulomb-exchange
                  call twoctr('G', D(2+j), D(2+de(2)+k), RR(:3*na))
                  E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            = E( 1+de(1)*(j+de(2)*k) : de(1)*(1+j+de(2)*k) ) &
                            + RR(c(1):c(1)+de(1)-1)
               end do
            end do
         else
            call quit('prop_twoave: GEO, nd > 2 not implemented')
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
               call quit('prop_twoave: MAG MAG not implemented for DFT')
            end if
            do k = 0, pd-1
               do j = 0, de(2)-1
                  do i = 0, de(1)-1
                     ii = min(c(1)+i,c(2)+j) &
                        + max(c(1)+i,c(2)+j) * (max(c(1)+i,c(2)+j)-1) / 2
                     ! ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
                     E(1+i+de(1)*(j+de(2)*k)) = -tr(A(ii), D(pd1-pd+1+k)) &
                                              * merge(1/2d0,1d0,nd==0) !half for magnetizability
!                    merge(1/2d0,1d0,nd==0) means multiply by 0.5 if nd == 0
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
                  call quit('prop_twoave: MAG MAG not implemented for DFT')
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
            call quit('prop_twoave: MAG MAG and nd > 2 not implemented')
         end if
      else if (np==2 .and. all(p==(/'GEO','GEO'/))) then
         allocate(RR(3*na*3*na))
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
               call quit('prop_twoave: GEO GEO, DFT not implemented',get_print_unit())
            end if
            do j = 0, de(2)-1
               do i = 0, de(1)-1
                  E(1+i+de(1)*(j+de(2)*k)) = RR(c(1)+i + 3*na*(c(2)+j-1))
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
                     call quit('prop_twoave: GEO GEO, DFT not implemented',get_print_unit())
                  end if
                  do j = 0, de(2)-1
                     do i = 0, de(1)-1
                        E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                     = E(1+i+de(1)*(j+de(2)*(k+de(3)*l))) &
                                     + RR(c(1)+i + 3*na*(c(2)+j-1))
                     end do
                  end do
               end do
            end do
         else
            call quit('prop_twoave: GEO GEO, nd > 2 not implemented',get_print_unit())
         end if
         deallocate(RR)
      else
         print *,'prop_twoave: no integrals for these perturbations: ', &
                 (p(i)//' ', i=1,np)
         call quit('prop_twoave: no such integrals',-1)
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine



   !> Calculates the 1-electron integrals perturbed by the perturbations p(:)
   !> and ADDS the integrals to the perturbed Fock matrices F(:).
   !> Front for the private subroutine 'oneint' below, checking the
   !> arguments' dimensions, and doing permutations
   !> S0 is passed as argument only as reference to nuclei and basis set
   subroutine prop_oneint(S0, p, dimp, F, S, comp, freq)
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> shape(F), size(dimp) = size(p), dimensions of perturbed
      !> Fock matrices F(:)
      integer,           intent(in) :: dimp(:)
      !-----------------------------------------------------------
      !> perturbed Fock matrices to fill with perturbed integrals,
      !> size(F) = product(dimp)
      type(matrix), target, optional, intent(inout) :: F(*)
      !> optionally return the corresponding perturbed overlap matrices
      type(matrix), target, optional, intent(inout) :: S(*)
      !-----------------------------------------------------------
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer,      optional, intent(in)    :: comp(:)
      !> freq(np), complex frequencies
      !> for each p, default all zero. These multiply the
      !> half-derivative overlap integrals, thus no contribution
      !> if basis independent of p
      complex(8),   optional, intent(in)    :: freq(:)
      !-----------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dimp)), ccomp(size(p)), &
                      stepf(size(dimp)), ddimp(size(dimp)), idxf(size(dimp)), &
                      i, j, k, tmpi
      type(matrix) :: Ftmp(product(dimp)), Stmp(product(dimp))
      character(4) :: pp(size(p)), tmpp
      complex(8)   :: ffreq(size(p)), tmpf
      logical      :: zero, bas
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero. If there is more than 1
      ! linear perturbation (like 'EL') or more than two quadratic (like 'MAGO')
      zero = .false.
      j = 0
      k = 0
      do i = 1, size(p)
         if (.not.field_list(idxp(i))%lin .and. &
             .not.field_list(idxp(i))%quad) cycle
         ! if both linear and quadratic, interpret as ZERO (like EXCI)
         if (field_list(idxp(i))%lin .and. &
             field_list(idxp(i))%quad) zero = .true.
         if (j /= 0) then !if second lin or quad perturbation
            ! if j was EL, then i is second EL (or first MAGO)
            if (field_list(idxp(j))%lin) zero = .true.
            ! if j was MAGO, and different from i, also zero
            if (field_list(idxp(j))%quad .and. &
               idxp(i) /= idxp(j)) zero = .true. !two different quadratic
            !third linear or quadratic perturbation
            if (k /= 0) zero = .true.
            k = j
         end if
         j = i
      end do
      ! check dimensions argument dimp, verify that dimensions are positive
      if (size(dimp) /= size(p)) call quit('prop_oneint argument error: ' &
               // 'Different number of perturbations and dimensions, size(dimp) /= size(p)',get_print_unit())
      if (any(dimp <= 0)) call quit('prop_oneint argument error: ' &
               // 'Perturbations have a zero or negative dimension, dimp <= 0',get_print_unit())
      ! compute step lengths in F (cumulative product of dimp)
      stepf(1) = 1
      do i = 2, size(dimp)
         stepf(i) = stepf(i-1)*dimp(i-1)
      end do
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_oneint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',get_print_unit())
         if (any(comp <= 0)) call quit('prop_oneint argument error: ' &
                    // 'Lowest component indices must be positive, but comp <= 0',get_print_unit())
         ccomp = comp
      end if
      if (any(ccomp + dimp - 1 > pert_shape(p))) &
         call quit('prop_oneint argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dimp - 1 > pert_shape(p)',get_print_unit())
      ! check optional argument freq, default to zero
      ffreq = 0
      if (present(freq)) then
         if (size(freq) /= size(p)) call quit('prop_oneave ' &
               // 'argument error: Wrong number of frequencies, size(freq) /= size(p)',get_print_unit())
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
         if (present(F)) F(i) = 0*S0
         if (present(S)) S(i) = 0*S0
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
            call mat_move(Ftmp(j), F(i))
         if (present(S) .and. bas) &
            call mat_move(Stmp(j), S(i))
         do k = 1, size(dimp)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimp(k)) exit
            idxf(k) = 0
         end do
      end do
   end subroutine



   subroutine oneint(S0, np, p, c, dp, w, F, S)
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> number of perturbations
      integer,           intent(in) :: np
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> lowest component of each perturbation
      integer,           intent(in) :: c(np)
      !> dimensions of property (F and S)
      integer,           intent(in) :: dp(np)
      !> frequency of each p
      complex(8),        intent(in) :: w(np)
      !--------------------------------------------
      !> perturbed Fock matrices
      type(matrix),   intent(inout) :: F(*)
      !> perturbed overlap matrices
      type(matrix),   intent(inout) :: S(*)
      !--------------------------------------------
      integer      :: i, j, k, ii, jj, kk
      type(matrix) :: A(6) !scratch matrices
      real(8)      :: R(1) !dummy
      if (np==0) then
         call quit('prop_oneint error:' &
                // ' Unperturbed one-electron integrals requested',get_print_unit())
      else if (np==1 .and. p(1)(1:3)=='AUX') then
         read (p(1)(4:),'(i1)') i
         if (i < 0 .or. i > 9) call quit('prop_oneint error:' &
                  // ' Index in auxiliary label out of range 0..9',get_print_unit())
         call load_oneint(prop_auxlab(i), F(1))
      else if (np==1 .and. p(1)=='EL') then
         do i = 0, dp(1)-1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint(xyz(c(1)+i) // 'DIPLEN ', F(i+1))
         end do
      else if (np==1 .and. p(1)=='VEL') then
         do i = 0, dp(1)-1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint(xyz(c(1)+i) // 'DIPVEL ', F(i+1))
         end do
      else if (np==1 .and. p(1)=='MAGO') then
         do i = 0, dp(1)-1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint(xyz(c(1)+i) // 'ANGMOM ', F(i+1))
            F(i+1) = (1/2d0) * F(1+i) !factor 1/2 here
         end do
      else if (np==1 .and. p(1)=='MAG') then
         do i = 0, dp(1)-1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint('dh/dB' // xyz(c(1)+i) // '  ', F(i+1))
            S(i+1) = mat_alloc_like(S0)
            if (w(1)==0) then !no -i/2 Tb contribution
               call load_oneint('dS/dB' // xyz(c(1)+i) // '  ', S(i+1))
            else
               call load_oneint('d|S>/dB' // xyz(c(1)+i), S(i+1))
               F(i+1) = F(i+1) - w(1)/2 * S(i+1)
               F(i+1) = F(i+1) - w(1)/2 * trps(S(i+1))
               S(i+1) = S(i+1) - trps(S(i+1))
            end if
         end do
      else if (np==1 .and. p(1)=='ELGR') then
         do i = 0, dp(1)-1
            ii = c(1)+i + merge(3, merge(2, 0, c(1)+i > 1), c(1)+i > 3) - 1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint(xyz(1+mod(ii,3)) // xyz(1+ii/3) // 'THETA ', &
                             F(1+i))
         end do
      else if (np==1 .and. p(1)=='GEO') then
         do i = 0, dp(1)-1
            F(i+1) = mat_alloc_like(S0)
            call load_oneint('1DHAM' // prefix_zeros(c(1)+i,3), F(i+1))
            S(i+1) = mat_alloc_like(S0)
            if (w(1)==0) then !no -i/2 Tb contribution
               call load_oneint('1DOVL' // prefix_zeros(c(1)+i,3), S(i+1))
               S(i+1) = -S(i+1) !1DOVL is really -dS/dg
            else
               call load_oneint('SQHDR' // prefix_zeros(c(1)+i,3), S(i+1))
               S(i+1) = -S(i+1) !SQHDR is really -dS>/dg
               F(i+1) = F(i+1) - w(1)/2 * S(i+1)
               F(i+1) = F(i+1) + w(1)/2 * trps(S(i+1))
               S(i+1) = S(i+1) + trps(S(i+1))
            end if
         end do
      else if (np==2 .and. all(p==(/'MAG','EL '/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               F(1+i+dp(1)*j) = mat_alloc_like(S0)
               call load_oneint(xyz(c(2)+j) // '-CM1 ' &
                             // xyz(c(1)+i) // ' ', F(1+i+dp(1)*j))
            end do
         end do
      else if (np==2 .and. all(p==(/'MAGO','MAGO'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               F(1+i+dp(1)*j) = mat_alloc_like(S0)
               call load_oneint(xyz(min(c(1)+i,c(2)+j))              &
                             // xyz(max(c(1)+i,c(2)+j)) // 'SUSCGO', &
                                F(1+i+dp(1)*j))
               !ajt Since MAGO is imaginary/anti-symmetric, we get a minus sign here
               F(1+i+dp(1)*j) = -F(1+i+dp(1)*j)
            end do
         end do
      else if (np==2 .and. all(p==(/'MAG','MAG'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               ii = 1+i+dp(1)*j
               ! Hamiltonian integrals and -i/2 Tbb in A(1), Sbb in A(2)
               F(ii) = mat_alloc_like(S0)
               call load_oneint(xyz(min(c(1)+i,c(2)+j)) &
                             // xyz(max(c(1)+i,c(2)+j)) // 'dh/dB2', F(ii))
               ! If both fields static, no -i/2 Tbb to A(1)
               S(ii) = mat_alloc_like(S0)
               if (w(1)==0 .and. w(2)==0) then
                  call load_oneint('dS/dB2' // xyz(min(c(1)+i,c(2)+j))  &
                                            // xyz(max(c(1)+i,c(2)+j)), S(ii))
               else
                  !only load >> XX XY YY XZ YZ ZZ, use trps(..) for <<
                  call load_oneint('>>S/B2' // xyz(min(c(1)+i,c(2)+j))  &
                                            // xyz(max(c(1)+i,c(2)+j)), S(ii))
                  F(ii) = F(ii) - (w(1)+w(2))/2 * S(ii)
                  F(ii) = F(ii) + (w(1)+w(2))/2 * trps(S(ii))
                  S(ii) = S(ii) + trps(S(ii))
                  !only load <> XX XY YY XZ YZ ZZ, use trps(..) for the others
                  A(1) = mat_alloc_like(S0)
                  call load_oneint('<>S/B2' // xyz(min(c(1)+i,c(2)+j)) &
                                            // xyz(max(c(1)+i,c(2)+j)), A(1))
                  !ajt&kk 0410 Changed sign on these two: (think it's correct)
                  F(ii) = F(ii) - merge(-1, 1, c(1)+i < c(2)+j) &
                                         * (w(1)-w(2))/2 * A(1)
                  F(ii) = F(ii) + merge(-1, 1, c(1)+i < c(2)+j) &
                                         * (w(1)-w(2))/2 * trps(A(1))
                  S(ii) = S(ii) + A(1)
                  S(ii) = S(ii) + trps(A(1))
                  A(1) = 0 !free
               end if
               !ajt Since MAG is imaginary/anti-symmetric, we get a minus sign here
               F(ii) = -F(ii)
               S(ii) = -S(ii)
            end do
         end do
      else if (np==2 .and. all(p==(/'ELGR','MAG '/))) then
         do i = 1, 6
            A(i) = mat_alloc_like(S0)
         end do
         do j = 0, dp(2)-1 !MAG indices
            ! ajt The integrals 'XY-QDB Z' are not traceless like THETA.
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
               F(1+i+dp(1)*j) = A(c(1)+i)
            end do
         end do
         A(:) = 0 !free
      else if (np==2 .and. all(p==(/'GEO','EL '/))) then
         do j = 0, dp(2)-1 !EL indices
            do i = 0, dp(1)-1 !GEO indices
               F(1+i+dp(1)*j) = mat_alloc_like(S0)
               call load_oneint(prefix_zeros(c(1)+i,2) // ' DPG ' &
                                // xyz(c(2)+j), F(1+i+dp(1)*j))
               F(1+i+dp(1)*j) = -F(1+i+dp(1)*j) !minus sign
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO ','MAGO'/))) then
         do j = 0, dp(2)-1
            do i = 0, dp(1)-1
               F(1+i+dp(1)*j) = mat_alloc_like(S0)
               call load_oneint(prefix_zeros(c(1)+i,3) // 'AMDR' &
                                // xyz(c(2)+j), F(1+i+dp(1)*j))
               F(1+i+dp(1)*j) = 1/2d0 * F(1+i+dp(1)*j) !factor 1/2
            end do
         end do
      else if (np==2 .and. all(p==(/'GEO ','ELGR'/))) then
         do i = 1, 6
            A(i) = mat_alloc_like(S0)
         end do
         do i = 0, dp(1)-1 !GEO indices
            ! ajt The integrals '01QDG XY' are not traceless.
            !     Therefore load all integrals and remove trace manually
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXX', A(1))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGYY', A(3))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGZZ', A(6))
            A(2) = A(1) + A(3) + A(6) !trace in A(2)
            A(1) = A(1) - 1/3d0*A(2)
            A(3) = A(3) - 1/3d0*A(2)
            A(6) = A(6) - 1/3d0*A(2)
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXY', A(2))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGXZ', A(4))
            call load_oneint(prefix_zeros(c(1)+i,3) // 'QDGYZ', A(5))
            do j = 0, dp(2)-1 !ELGR indices
               F(1+i+dp(1)*j) = -3/2d0 * A(c(2)+j) !factor -3/2
            end do
         end do
         A(:) = 0
      else if (np==3 .and. all(p==(/'MAG','MAG','EL '/))) then
         do k = 0, dp(3)-1
            do j = 0, dp(2)-1
               do i = 0, dp(1)-1
                  F(1+i+dp(1)*(j+dp(2)*k)) = mat_alloc_like(S0)
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
         call quit('prop_oneint: No such integrals',get_print_unit())
      end if
   end subroutine



   !> Contracts the 2-electron and Kohn-Sham integrals perturbed by the
   !> perturbations p(:) with the perturbed density matrix expansion in D(:)
   !> (e.g. D=(/D,Dx,Dy,Dxy/) for a 2nd order expansion), and ADD the resulting
   !> Fock matrix contibution to the array F(:).
   !> Front for the private subroutine 'twoave' below, checking the arguments'
   !> dimensions, and doing permutations
   !> D(1) serves as reference to nuclei, basis and model/functional
   subroutine prop_twoint(p, D, dimf, F, perm, comp)
      !> p(np) perturbation lables
      character(*),      intent(in) :: p(:)
      !> (un)perturbed density matrices to contract perturbed
      !> one-electron integrals with.
      !> If perm present, size(D) = product(1+dime(perm(np+1:np+nd))),
      !> if not present, size(D) = product(1+dime(np+1:np+nd))
      type(matrix),      intent(in) :: D(:)
      !> dime(np+nd) = shape(F), dimensions of perturbed Fock matrices F(:)
      integer,           intent(in) :: dimf(:)
      !---------------------------------------------------------------
      !> Perturbed Fock matrices, works incrementally,
      !> thus contributions are ADDED to F(*). size(F) = product(dimf)
      type(matrix), intent(inout) :: F(*)
      !---------------------------------------------------------------
      !> perm(np+nd), permutation of indices.
      !> For each dimension of p and D, the corresponding dimension of F.
      !> Default 1 2 ... np+nd (no permutation)
      integer,      optional, intent(in) :: perm(:)
      !> comp(np), starting component index for each p. Default 1 1 ... 1
      integer,      optional, intent(in) :: comp(:)
      !---------------------------------------------------------------
      integer      :: idxp(size(p)), pperm(size(dimf)), ccomp(size(p)), &
                      stepf(size(dimf)), ddimf(size(dimf)), idxf(size(dimf)), &
                      i, j, k, tmpi, nd
      character(4) :: pp(size(p)), tmpp
      type(matrix) :: Ftmp(product(dimf))
      logical      :: zero
      ! verify that all perturbation labels exist, find index of each
      idxp = (/(idx(p(i)), i=1,size(p))/)
      ! determine whether these integrals are zero.
      zero = any((/(field_list(idxp(i))%lin  .or. &
                    field_list(idxp(i))%quad .or. &
               .not.field_list(idxp(i))%bas, i=1,size(p))/))
      ! check dimensions argument dime, verify that dimensions are positive
      if (size(dimf) < size(p)) call quit('prop_twoint argument error: ' &
               // 'More perturbations than dimensions of F(:), size(dimf) < size(p)',-1)
      if (any(dimf <= 0)) call quit('prop_twoint argument error: ' &
               // 'Perturbed Fock F(:) has a zero or negative dimension, dimf <= 0',-1)
      ! compute step lengths in E (cumulative products of dimensions)
      stepf(1) = 1
      do i = 2, size(dimf)
         stepf(i) = stepf(i-1)*dimf(i-1)
      end do
      ! reorder dimensions and step lengths in F according to permutation argument perm
      ddimf = dimf
      if (present(perm)) then
         if (size(perm) /= size(dimf)) call quit('prop_twoint argument ' &
               // 'error: Wrong length of permutation vector, size(perm) /= size(dimf)',-1)
         ! verify that perm is indeed a permutation
         do i = 1, size(dimf)-1
            if (perm(i) <= 0 .or. perm(i) > size(dimf) .or. &
                any(perm(i) == perm(i+1:size(dimf)))) call quit('prop_twoint ' &
                      // 'argument error: Permutation must contain each number exactly once',-1)
         end do
         ddimf = (/( dimf(perm(i)), i=1,size(dimf))/)
         stepf = (/(stepf(perm(i)), i=1,size(dimf))/)
      end if
      ! check optional arg comp, default to 1 1 ... 1
      ccomp = 1
      if (present(comp)) then
         if (size(comp) /= size(p)) call quit('prop_twoint argument error: ' &
                    // 'Wrong number of lowest component indices, size(comp) /= size(p)',-1)
         if (any(comp <= 0)) call quit('prop_twoint argument error: ' &
                    // 'Lowest component indices must be positive, comp <= 0',-1)
         ccomp = comp
      end if
      if (any(ccomp + ddimf(:size(p)) - 1 > pert_shape(p))) &
         call quit('prop_twoint argument error: Lowest component index plus ' &
                // 'dimension exceeds dimension of perturbation, comp + dimf > pert_shape(p)',-1)
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
               // 'perturbed densities D does not correspond to dimf (and perm)',-1)
      ! arguments checked. If perturbed integrals are zero, verify all D defined, return
      if (zero) then
         if (.not.all((/(isdef(D(i)), i=1,nd)/))) &
            call quit('prop_twoint error: Undefined matrix in argument D(:)',-1)
         return
      end if
      ! permute perturbed Fock matrices in F(:) over to Ftmp(:)
      idxf = 0
      do j = 1, product(dimf)
         i = 1 + sum(idxf * stepf)
         call mat_move(F(i), Ftmp(j))
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
         call mat_move(Ftmp(j), F(i))
         do k = 1, size(dimf)
            idxf(k) = idxf(k) + 1
            if (idxf(k) /= ddimf(k)) exit
            idxf(k) = 0
         end do
      end do
   end subroutine



   subroutine twoint(np, nd, p, c, df, D, F)
      !> number of perturbations and order of density
      integer,           intent(in) :: np, nd
      !> perturbation lables
      character(*),      intent(in) :: p(np)
      !> lowest component of each perturbation
      integer,           intent(in) :: c(np)
      !> dimensions of perturbed Fock matrix F
      integer,           intent(in) :: df(np+nd)
      !> un-/perturbed density matrices (expension),
      !> size(D) = product(1+df(np+1:np+nd))
      type(matrix),      intent(in) :: D(*)
      !--------------------------------------------------
      !> where to ADD property contributions
      !> (works incrementally), size(F) = product(df)
      type(matrix),   intent(inout) :: F(*)
      !--------------------------------------------------
      integer      :: i, j, k, l, ii, jj, kk, ll, pd, pd1
      type(matrix) :: A(6) !scratch matrices
      do i = 1, size(A)
         A(i) = 0*D(1) !scratch matrices
      end do
      pd  = product(df(np+1:np+nd))   !product of density dimensions
      pd1 = product(1+df(np+1:np+nd)) !size of D(*)
      if (np==0) then
         if (nd==0) call quit('prop_twoint error: Unperturbed ' &
                           // 'two-electron Fock matrix requested',-1)
         ! di_get_gmat and di_get_sigma expects A initialized (allocated)
         if (iszero(A(1))) call mat_ensure_alloc(A(1))
         do i = 0, pd-1
            if (iszero(D(pd1-pd+1+i))) cycle
            ! Coulomb-exchange
            call twofck('  ', D(pd1-pd+1+i), A(1:1))
            F(i+1) = F(i+1) + A(1)
            ! Kohn-Sham exchange-correlation
          ! call twofck_ks(1, (/D(1),D(pd1-pd+1+i)/), F(1+i))
            call rsp_xcint(D=(/D(1), D(pd1-pd+1+i)/), F=F(1+i))
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
            call quit('prop_twoint: nd > 2 not implemented with DFT',-1)
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
         call quit('prop_twoint: No such integrals',-1)
      end if
      A(:) = 0 !delete scratch matrices
   end subroutine

   subroutine load_oneint(lab, A)
      character(*),           intent(in)    :: lab
      type(matrix),           intent(inout) :: A
      if (.not.isdef(A)) &
         call quit('load_oneint: matrix A undefined - must be defined',-1)
      ! if A is zero, get it allocated
      if (iszero(A)) call mat_ensure_alloc(A)
      ! read integrals from file into A
      call legacy_read_integrals(lab, A)
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
      n2 = size(D(1)%elms_alpha)
      if (what=='  ' .and. nf==1) then
         lwrk = 0
      else if (what=='M ' .and. nf==3) then
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else if (what=='G ' .and. nf==3) then
         if (.not.present(a)) &
            call quit("prop_integrals twofck: Geometry Fock but atom 'a' not specified",-1)
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
         aa = a !may not be present, replicate
      else if (what=='MM' .and. nf==6) then
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else if (what=='GG' .and. nf==3*3) then
         if (.not.present(a) .or. .not.present(b)) &
            call quit("prop_integrals twock: 2nd ord. geom Fock reqires atoms 'a','b' specified",-1)
         !lwrk = (1+2*nf+50)*n2 + 10000*nbas + 5000000
         lwrk = 50*n2 + 10000*D(1)%nrow + 5000000
      else
         call quit('prop_integrals twofck: wrong size(F) for what=' // what,-1)
      end if
      !make usre F is properly allocated
      do i = 1, nf
         if (iszero(F(i))) call mat_ensure_alloc(F(i))
      end do
      !allocate work
      if (lwrk/=0) then
         call f77_memory_select(work_len=lwrk, work=wrk)
      end if
      !call integral program
      if (what=='  ') then
         call interface_scf_get_g(D(1), F(1))
      else
         wrk(1:n2) = reshape(D(1)%elms_alpha, (/n2/)) !ajt fixme
         call GRCONT(wrk( 1+n2+n2*nf : lwrk ), (lwrk-n2-n2*nf),         &
                     wrk( 1+n2 : n2+n2*nf ), n2*nf, (what(1:1) == 'G'), &
                     (what(1:1) == 'M'), merge(1,2,what(2:2)==' '),     &
                     aa, .false., .true., wrk(1:n2), 1)
         do i = 1, nf
            F(i)%elms_alpha = reshape(wrk(1+n2*i:n2*(1+i)), shape(F(i)%elms_alpha)) !ajt fixme
         end do
      end if
      !deallocate work
      if (lwrk /= 0) then
         call f77_memory_deselect(work_len=lwrk, work=wrk)
      end if
   contains
      subroutine subr(D, F)
         type(matrix), intent(in)    :: D
         type(matrix), intent(inout) :: F(nf)
      end subroutine
   end subroutine



   subroutine twofck_ks(n, D, F)

#ifndef PRG_DIRAC
      use xcint_main
#endif

      integer,      intent(in)    :: n
      type(matrix)                :: D(n+1)
      type(matrix), intent(inout) :: F
      type(matrix) :: A
      integer      :: i

      if (.not. do_dft()) then
!        we want no xc
         return
      end if

      !radovan: this is not correct for CR and higher
      !         there you can have nonzero integral with one or more
      !         density matrices zero
      !ajt: twofck_ks is only supposed to return a single multilinear
      !     contribution Exc[D1]^n D2 D3 D4...Dn+1, so zero when either of
      !     D(2:) are zero.
      if (any((/(iszero(D(i)), i = 2, n+1)/))) then
         ! at least one of the density matrices is zero
         return
      end if

      A = tiny(0.0d0)*D(1)

      select case (n)
         case (1)
            ! nothing to be done
         case (2)
!           call xc_integrate(xc_mat_dim=D(1)%nrow, D=D, xc_nr_dmat=3, F=(/A/))
!           F = F + A
         case default
            call quit('prop_contribs_old/twofck_ks error: contrib not implemented',-1)
      end select

      A = 0

   end subroutine



   !> what=G : 3*natoms geometric
   !> what=M : 3 magnetic
   subroutine twoctr(what, Da, Db, E)
      character(*),      intent(in) :: what
      type(matrix),      intent(in) :: Da, Db
      real(8),        intent(inout) :: E(:)
      real(8), pointer :: wrk(:)
      integer          :: lwrk, na, nb, l

      na = get_nr_atoms()
      nb = Da%nrow

      !ajt MagSus contraction doesn't work, so I disabled this. Use twofck instead
      if (what=='M' .or. what=='MM') &
          & call quit("prop_contribs/twoctr: 'M' or 'MM' not yet implemented",-1)
      if (what=='G' .and. size(e) == 3*na) then
         lwrk = 50*nb**2 + 10000*nb + 5000000
      else if (what=='M' .and. size(e)==3) then
         lwrk = 50*nb**2 + 10000*nb + 5000000
      else if (what=='GG' .and. size(e) == 3*3*na**2) then
         lwrk = 50*nb**2 + 10000*nb + 10000000
      else if (what=='MM' .and. size(e) == 3*3) then
         lwrk = 50*nb**2 + 10000*nb + 5000000
      else
          call quit('prop_contribs/twoctr: wrong size(E) for what=' // what,get_print_unit())
      end if

      call f77_memory_select(work_len=lwrk, work=wrk)

!     this is done because grcont presently needs Da and Db to be consecutive
!     and (/Da%elms_alpha, Db%elms_alpha/) can cause stack overflow, depending
!     on ulimit and compiler flags; -auto-scalar:fine: -auto:overflow
!     a future rewrite of grcont could take a list of pointers
!     (or just type(matrix)) instead of consecutive real(8) arrays
      l = size(Da%elms_alpha)
      wrk(    1:  l) = reshape(Da%elms_alpha, (/l/))
      wrk(l + 1:2*l) = reshape(Db%elms_alpha, (/l/))

      e = 0.0d0


      call grcont(wrk(2*l + 1),    &
                  lwrk - 2*l,      &
                  e,                &
                  int(size(e)),     &
                  what(1:1) == 'G', &
                  what(1:1) == 'M', &
                  int(len(what)),   &
                  0,                &
                  .true.,           &
                  .false.,          &
                  wrk(1),             &
                  2)


      call f77_memory_deselect(work_len=lwrk, work=wrk)

   end subroutine

   function idx(p)
      character(*) :: p
      integer      :: idx, i, j
      do idx = 1, size(field_list)
          if (field_list(idx)%code == p) return
      end do
      if (p(1:min(3,len(p))) == 'AUX') then
         read (p(min(4,len(p)):),'(i1)', iostat=j) i
         if (j /= 0 .or. i < 0 .or. i > 9) &
            call quit('prop_contribs error: ' &
                     // 'Auxiliary label index should be within 0..9, ' // p,-1)
         idx = 1
         return
      end if
      call quit('Perturbation not found: ' // p,-1)
   end function


   function pert_antisym(p)
      character(*), intent(in) :: p(:)
      logical :: pert_antisym(size(p))
      integer :: i
      pert_antisym = (/(field_list(idx(p(i)))%anti, i=1,size(p))/)
   end function


   !> shape (dimensions) of property p(:)
   function pert_shape(p)
      !> field lables
      character(*),      intent(in) :: p(:)
      integer :: pert_shape(size(p)), i
      pert_shape = (/(field_list(idx(p(i)))%ncomp, i=1,size(p))/)
      ! loop through mol-dependent
      do i=1, size(p)
         if (pert_shape(i) /= -1) then
            ! cycle
         else if (p(i) == 'GEO') then
            pert_shape(i) = 3 * get_nr_atoms()
         else
            call quit('pert_shape error: Number of comp. unknown for ' // p(i),get_print_unit())
         end if
      end do
   end function


   function pert_basdep(p)
      character(*), intent(in) :: p(:)
      logical :: pert_basdep(size(p))
      integer :: i
      pert_basdep = (/(field_list(idx(p(i)))%bas, i=1,size(p))/)
   end function

end module
