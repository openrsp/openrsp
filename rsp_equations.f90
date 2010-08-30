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

   use matrix_defop  !matrix type and operators
   use prop_contribs !integrals and integral contractions

   ! ajt LSDALTON has replaced the (global) quit with lsquit
   !     with unit (lupri) as extra argument, which doesn't
   !     exist in DIRAC. For now, this macro gets around that.
#ifdef LSDALTON_ONLY
#define quit(msg) lsquit(msg,-1)
#endif

   implicit none

   public pert_fock    !contract perturbed Fock-matrices
   public pert_scf_eq  !contract perturbed scf equation residuals
   public pert_dens    !perturbed (or response-) densities and Fock-matrices
   public rsp_equations_debug  !turn debug printing on or off

   ! ajt Freqency derivative hack
   public solve_scf_eq

   ! ajt The following three are public for the moment, but are intended
   !     be accessed through subroutines eventually
   public cached_sol, solved_eqs, num_solved_eqs

   !> turn on or off debugging in this file
   logical :: rsp_equations_debug = .false.

   !> Maximum order of equation solution which can be cached.
   !> Defines the size of the fields in type cached_sol.
   !> ajt fixme Should be removed once pert_dens can do arbitrary-order
   !> equations.
   integer, parameter :: max_order = 3


   !> Type for saving (caching) solutions of response equations,
   !> and avoid re-solving the same equations later in the program.
   type cached_sol
      !> order of equation solved
      integer      :: order            = 0
      !> perturbations, NONE=empty
      character(4) :: pert(max_order)  = 'NONE'
      !> number of components
      integer      :: ncomp(max_order) = 0
      !> starting component
      integer      :: scomp(max_order) = 1
      !> frequencies (0d0,0d0)
      complex(8)   :: freq(max_order)  = 0
      !> cached perturbed densities, if any. We could potentially cache 
      !> other matrices too: F, DS, SD, FD, DF, DFDp, FpDS or DpSD.
      !> Also, to save memory, files can be used, with filenames stored here.
      type(matrix), pointer :: D(:) => null()
      !ajt integer :: remaining_uses = 0 !free D when this hits zero.
   end type


   !> To keep collection of saved response equation solutions.
   !> ajt fixme Currently hard-coded size 100
   type(cached_sol), target, save :: solved_eqs(100)

   !> number of solutions currently cached
   integer :: num_solved_eqs = 0

   ! ajt For some reason, if I put this 'private' above type cached_sol,
   !     Doxygen does not document the fields inside the type, dispite
   private !it being declared public at the top. Might be a bug in Doxygen.

contains


   !> Contract perturbed integrals with (an expansion of) perturbed density
   !> matrices to get the corresponding perturbed Fock matrix.
   !> Optionally, also return the corresponding perturbed overlap integrals
   subroutine pert_fock(mol, p, dimp, D, F, S, comp, freq)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(prop_molcfg), intent(in) :: mol
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> dimension of each p, and thus also of F
      integer,           intent(in) :: dimp(:)
      !> density matrix expansion for p(:) unperturbed and perturbed
      !> density matrices. size(D) = product(1+dimp)
      type(matrix),      intent(in) :: D(:)
      !> output, perturbed Fock matrices. Deferred shape (*) to permit
      !> re-ranking. size(F)=product(dimp)
      type(matrix),   intent(inout) :: F(*)
      !--------------------------------------------------------------------
      !> optionally also output perturbed overlap matrices.
      !> Deferred shape (*) to permit re-ranking
      type(matrix), optional, intent(inout) :: S(*)
      !> starting component within each p, default 1 1 .. 1
      integer,      optional, intent(in)    :: comp(:)
      !> frequency of each p default all zero
      complex(8),   optional, intent(in)    :: freq(:)
      !-----------------------------------------------
      integer      :: ccomp(size(p)), nf, nd
      complex(8)   :: ffreq(size(p))
      character(4) :: UNP(0)
      logical      :: bas(size(p))
      if (size(dimp) /= size(p)) call quit('pert_fock: Differing numbers ' &
               // 'of perturbations and dimensions, size(dimp) /= size(p)')
      nf = product(dimp)
      nd = product(1+dimp)
      if (size(D) /= nd) call quit('pert_fock: Wrong number of perturbed ' &
               // 'densities, size(D) /= product(dimp+1)')
      ccomp = 1
      if (present(comp)) then
         ccomp = comp !ajt fixme Check number and bounds
      end if
      ffreq = 0
      if (present(freq)) then
         ffreq = freq !ajt fixme Check number
      end if
      ! fill F with perturbed one-electron integrals, and
      ! optionally S with perturbed overlap integrals
      call prop_oneint(mol, D(1), p, dimp, F, S, comp=ccomp, freq=ffreq)
      ! add p-perturbed Coulomb-exchange+Kohn-Sham contracted against
      call prop_twoint(mol, p, D(1:1), dimp, F, comp=ccomp) !unperturbed D
      ! add unperturbed Coulomb-exchange+Kohn-Sham contribution
      call prop_twoint(mol, UNP, D, dimp, F)
      ! additional contributions for each basis-dependent subset of p
      bas = pert_basdep(p)
      if (size(p)==2 .and. bas(size(p))) &
         call prop_twoint(mol, p(2:2), D(1:1+dimp(1)), &
                          dimp, F, perm=(/2,1/), comp=ccomp(2:2))
      if (size(p)==2 .and. bas(1)) &
         call prop_twoint(mol, p(1:1), (/D(1),D(1+dimp(1)+1:1+dimp(1)+dimp(2))/), &
                          dimp, F, comp=ccomp(1:1))
      if (size(p) > 2 .and. any(bas)) &
         call quit('pert_fock: general case not implemented')
   end subroutine


!ajt: old arbitrary-order routine. Kept here so that I can look
!     at it if I at some point attempt to generalize pert_scf_eq.
!
!     subroutine prtfck(n,m,p,c,b,o,D,F)
!         character(*):: p(n); integer:: n,m,c(n),b(n+m-1),o(n+m-1)
!         type(MatPtr):: D(0:sum(b(n+1:))); type(MatPtr),intent(inout):: F(product(b))
!         type(Matrix):: DD(0:size(D)-1),FF(size(F))
!         integer:: i,j,k,id(m),bb(n+m-1),ie(n+m-1),se(n+m-1),np,nd,ndd,idd,nsy,nan,nim
!         logical:: im,z(0:size(D)-1),sym,ant,zero; real*8:: fac
!         call MatUse(D)
!         call MatDel(F)
!         im = (mod(count((/(PrtAnti(p(i)),i=1,n)/)),2)==1)
!         id = (/(1+sum(b(n+1:n+i)),i=0,m-1)/)
!         z = (/(IsZero(D(i)),i=0,size(D)-1)/)
!         if (z(0)) stop 'PrtFock: reference density zero'
!         bb = (/b(1:n),(count(.not.z(id(i):id(i+1)-1)),i=1,m-1)/)
!         k = 0
!         do i=0,size(D)-1
!             if (z(i)) cycle
!             if (D(i)%t) then
!                 call init_mat(DD(k),D(i)%m%m%ncol,D(i)%m%m%nrow)
!                 call mat_trans(D(i)%m%m,DD(k))
!             else
!                 DD(k) = D(i)%m%m
!             endif
!             k = k+1
!         enddo
!         if (any(bb==0)) then
!             do i=1,n; j=PrtRank(p(i)); enddo
!         else
!             call fckcore(n,m,p,c,bb,DD(0:k-1),FF(1:product(bb)))
!         endif
!         k = 0
!         do i=0,size(D)-1
!             if (z(i)) cycle
!             if (D(i)%t) call mat_free(DD(k))
!             k = k+1
!         enddo
!         !calculate index steps for the reordered output
!         se(o(n+m-1)) = 1
!         do i=n+m-2,1,-1
!             se(o(i)) = se(o(i+1))*b(o(i+1))
!         enddo
!         np = product(b(1:n))
!         nd = product(b(n+1:)); ie = 0
!         ndd = product(bb(n+1:)); idd = 0
!         do i=1,nd
!             zero = any((/(z(id(j)+ie(n+j)),j=1,m-1)/))
!             if (.not.zero) then
!                 nsy = count((/(D(id(j)+ie(n+j))%m%s,j=1,m-1)/))
!                 nan = count((/(D(id(j)+ie(n+j))%m%a,j=1,m-1)/))
!                 nim = count((/(D(id(j)+ie(n+j))%i,j=1,m-1)/))+merge(1,0,im)
!                 fac = product((/(D(id(j)+ie(n+j))%f,j=1,m-1)/))
!                 if (mod(nim,4)>=2) fac = -fac
!                 sym = (nsy+nan==m-1.and.(mod(nan,2)==0.neqv.im))
!                 ant = (nsy+nan==m-1.and.(mod(nan,2)==1.neqv.im))
!                 idd = idd+1
!             endif
!             !decide symmetry and scale factor for this D-tuple
!             do j=idd,idd+(np-1)*ndd,max(ndd,1)
!                 k = 1+sum(ie*se)
!                 if (zero.or.undef(FF(j))) then
!                     F(k) = 0*D(0)
!                 else
!                     call cast(FF(j),F(k),sym=sym,anti=ant)
!                     F(k)%f = fac; F(k)%i = (mod(nim,2)==1)
!                 endif
!                 do k=n,1,-1
!                     if (ie(k)/=b(k)-1) exit
!                     ie(k) = 0
!                 enddo
!                 ie(k) = ie(k)+1
!             enddo
!             ie(1) = 0
!             do j=n+m-1,n+1,-1
!                 if (ie(j)/=b(j)-1) exit
!                 ie(j) = 0
!             enddo
!             ie(j) = ie(j)+1
!         enddo
!         call MatRel(D)
!     end subroutine


   !> Computes the p'th derivative of idempotency relation DSD-D
   !> and SCF equation FDS-SDF, ie. the p-perturbed SCF equations.
   subroutine pert_scf_eq(mol, S0, p, dimp, Dp, Fp, FDSp, DSDp, comp, freq)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(prop_molcfg), intent(in) :: mol
      !> unperturbed overlap matrix
      type(matrix),      intent(in) :: S0
      !> perturbation lables
      character(*),      intent(in) :: p(:)
      !> number of components within each p (dimension)
      integer,           intent(in) :: dimp(:)
      !> perturbation expantion of density matrix
      type(matrix),      intent(in) :: Dp(:)
      !> perturbation expantion of Fock matrix
      type(matrix),      intent(in) :: Fp(:)
      !-----------------------------------
      !> resulting perturbed SCF equation (or residual)
      !> Deferred shape (*) to permit any rank. size(FDSp) = product(dimp)
      type(matrix),   intent(inout) :: FDSp(*)
      !> resulting perturbed idempotency condition
      !> Deferred shape (*) to permit any rank, size(DSDp) = product(dimp)
      type(matrix),   intent(inout) :: DSDp(*)
      !-------------------------------------
#ifdef PRG_DIRAC
      type(matrix) :: SDFp(1), X
#endif
      !-------------------------------------
      !> starting indices for each p, default 1
      integer,    optional, intent(in) :: comp(:)
      !> frequency of each p, default 0
      complex(8), optional, intent(in) :: freq(:)
      !------------------------------------------
      type(matrix), dimension(size(Dp)) :: Sp, DSp, SDp, FDp, DFp
      integer      :: np, nd, i, j, i0, j0, ij, c(size(p)) !copy of comp or 1
      complex(8)   :: w(size(p)) !copy of freq or zero
      character(4) :: UNP(0)
      logical      :: bas(size(p)), presnt
      np = product(dimp)
      nd = product(dimp+1)
      bas = pert_basdep(p)
      if (size(Dp) /= nd) call quit('pert_scf_eq: Wrong number of perturbed ' &
                // 'densities, size(Dp) /= product(dimp+1)')
      if (size(Fp) /= nd) call quit('pert_scf_eq: Wrong number of perturbed ' &
                // 'densities, size(Dp) /= product(dimp+1)')
      c = 1
      if (present(comp)) then
         c = comp !ajt fixme Check number and bounds
      end if
      w = 0
      if (present(freq)) then
         w = freq !ajt fixme Check number
      end if
      ! start by zeroing FDSp and DSDp
      do i = 1, np
         FDSp(i) = 0d0*S0
         DSDp(i) = 0d0*S0
      end do
      ! zero perturbed overlap Sp, and FDp/DFp (which multiply perturbed overlap)
      Sp  = (/(0d0*S0, i=1,size(Sp)) /)
      FDp = (/(0d0*S0, i=1,size(FDp))/)
      DFp = (/(0d0*S0, i=1,size(DFp))/)
      ! Contract -RHS's. The contraction formulas are derivatives of the
      ! rearranged first-order formula:
      ! d/da (F - i/2 d/dt S) D S - S D (F + i/2 tb\b S) =
      !     = (Fa - i/2 Sat) D S - S D (Fa + i/2 Sat)
      !     + (F D - i/2 St D - i S Dt) Sa - Sa (D F + i/2 D St + i Dt S)
      !     + (F Da - i/2 St Da - i/2 S Dat) S - S (Da F + i/2 Da St + i/2 Dat S)
      if (size(p)==0) then
         ! trivial case: Unperturbed equations
         FDSp(1) = Fp(1)*Dp(1)*S0 - S0*Dp(1)*Fp(1)
         DSDp(1) = Dp(1)*S0*Dp(1) - Dp(1)
      else if (size(p)==1) then
         presnt = .not.all((/(iszero(Dp(i)), i=2,1+dimp(1))/))
         if (presnt) then
            DSp(1) = Dp(1)*S0
            if (bas(1)) FDp(1) = Fp(1)*Dp(1)
            if (bas(1)) call prop_oneint(mol, S0, p, dimp, S=Sp(2:1+dimp(1)), comp=c)
            do i0 = 1, dimp(1)
               i = 1+i0
#ifdef PRG_DIRAC
!radovan: this is necessary otherwise matrices inherit wrong irep from Sp
!         which breaks point group symmetry in DIRAC
               Sp(i)%irep = Fp(i)%irep
#endif
               FDSp(i0) =            (Fp(i) - w(1)/2 * Sp(i)) * DSp(1)
               FDSp(i0) = FDSp(i0) - dag(DSp(1)) * (Fp(i) + w(1)/2 * Sp(i))
               FDSp(i0) = FDSp(i0) + FDp(1) * Sp(i) - Sp(i) * dag(FDp(1))
               FDSp(i0) = FDSp(i0) + (Fp(1) - w(1)/2 * S0) * Dp(i) * S0
               FDSp(i0) = FDSp(i0) - S0 * Dp(i) * (Fp(1) + w(1)/2 * S0)
               DSDp(i0) =            Dp(1) * Sp(i) * Dp(1) - Dp(i)
               DSDp(i0) = DSDp(i0) + Dp(i)*dag(DSp(1))
               DSDp(i0) = DSDp(i0) + DSp(1)*Dp(i)
               Sp(i)=0
            end do
            DSp(1)=0; FDp(1)=0
         end if
      else if (size(p)==2) then
         presnt = .not.all((/(iszero(Dp(ij)), ij = 2+dimp(1)+dimp(2), &
                                              (dimp(1)+1)*(dimp(2)+1))/))
         if (presnt) then
            DSp(1) = Dp(1) * S0
            if (all(bas)) FDp(1) = Fp(1) * Dp(1)
            if (all(bas))call prop_oneint(mol, S0, p, dimp, S=Sp(2+dimp(1)+dimp(2) : &
                                          (dimp(1)+1)*(dimp(2)+1)), comp=c)
            do i0 = 1, dimp(1)*dimp(2)
               i = 1+dimp(1)+dimp(2)+i0
#ifdef PRG_DIRAC
!radovan: this is necessary otherwise matrices inherit wrong irep from Sp
!         which breaks point group symmetry in DIRAC
               Sp(i)%irep = Fp(i)%irep
#endif
               FDSp(i0) =            (Fp(i) - (w(1)+w(2))/2 * Sp(i)) * DSp(1)
               FDSp(i0) = FDSp(i0) - dag(DSp(1)) * (Fp(i) + (w(1)+w(2))/2 * Sp(i))
               FDSp(i0) = FDSp(i0) + FDp(1) * Sp(i) - Sp(i) * dag(FDp(1))
               FDSp(i0) = FDSp(i0) + (Fp(1) - (w(1)+w(2))/2 * S0) * Dp(i) * S0
               FDSp(i0) = FDSp(i0) - S0 * Dp(i) * (Fp(1) + (w(1)+w(2))/2 * S0)
               DSDp(i0) =            Dp(1) * Sp(i) * Dp(1) - Dp(i)
               DSDp(i0) = DSDp(i0) + Dp(i)*dag(DSp(1))
               DSDp(i0) = DSDp(i0) + DSp(1)*Dp(i)
               Sp(i)=0
            end do
            DSp(1)=0; FDp(1)=0
         end if
         presnt = .not.all((/(iszero(Dp(i)), i=2,1+dimp(1))/)) .and. &
                  .not.all((/(iszero(Dp(j)), j=2+dimp(1),1+dimp(1)+dimp(2))/))
         if (presnt) then
            if (bas(1)) call prop_oneint(mol, S0, p(1:1), dimp(1:1), &
                                         S=Sp(2:1+dimp(1)), comp=c(1:1))
            if (bas(2)) call prop_oneint(mol, S0, p(2:2), dimp(2:2), S=Sp(2+dimp(1): &
                                         1+dimp(1)+dimp(2)), comp=c(2:2))
            do j0 = 0, dimp(2)-1
               j = 2+dimp(1)+j0
#ifdef PRG_DIRAC
!radovan: this is necessary otherwise matrices inherit wrong irep from Sp
!         which breaks point group symmetry in DIRAC
               Sp(j)%irep = Fp(j)%irep
#endif
               DSp(j) = Dp(j)*S0 + Dp(1)*Sp(j)
               SDp(j) = S0*Dp(j) + Sp(j)*Dp(1)
               if (bas(1)) FDp(j) = (Fp(j) - w(2)/2 * Sp(j)) * Dp(1) &
                                  + (Fp(1) - w(2) * S0) * Dp(j)
               if (bas(1)) DFp(j) = Dp(1) * (Fp(j) + w(2)/2 * Sp(j)) &
                                  + Dp(j) * (Fp(1) + w(2) * S0)
               do i0 = 0, dimp(1)-1
                  i = 2+i0
#ifdef PRG_DIRAC
!radovan: this is necessary otherwise matrices inherit wrong irep from Sp
!         which breaks point group symmetry in DIRAC
                  Sp(i)%irep = Fp(i)%irep
#endif
                  ij = 1+i0+dimp(1)*j0
                  !               (Fa - i/2 Sat) D S - S D (Fa + i/2 Sat)
                  ! +     (F D - i/2 St D - i S Dt) Sa - Sa (D F + i/2 D St + i Dt S)
                  ! + (F Da - i/2 St Da - i/2 S Dat) S - S (Da F + i/2 Da St + i/2 Dat S)
                  FDSp(ij) = FDSp(ij) + (Fp(i) - w(1)/2 * Sp(i)) * DSp(j)
                  FDSp(ij) = FDSp(ij) - SDp(j) * (Fp(i) + w(1)/2 * Sp(i))
                  FDSp(ij) = FDSp(ij) + FDp(j)*Sp(i)
                  FDSp(ij) = FDSp(ij) - Sp(i)*DFp(j)
                  FDSp(ij) = FDSp(ij) + (Fp(j) - (w(1)+w(2))/2 * Sp(j)) * Dp(i) * S0
                  FDSp(ij) = FDSp(ij) - S0 * Dp(i) * (Fp(j) + (w(1)+w(2))/2 * Sp(j))
                  FDSp(ij) = FDSp(ij) + (Fp(1) - w(1)/2 * S0) * Dp(i) * Sp(j)
                  FDSp(ij) = FDSp(ij) - Sp(j) * Dp(i) * (Fp(1) + w(1)/2 * S0)
                  DSDp(ij) = DSDp(ij) + Dp(j) * Sp(i) * Dp(1)
                  DSDp(ij) = DSDp(ij) + Dp(1) * Sp(i) * Dp(j)
                  DSDp(ij) = DSDp(ij) + Dp(i)*SDp(j)
                  DSDp(ij) = DSDp(ij) + DSp(j)*Dp(i)
               end do
               DSp(j)=0; SDp(j)=0; FDp(j)=0; DFp(j)=0; Sp(j)=0
            end do
            Sp(2:1+dimp(1)) = 0
         end if
      else
#ifdef PRG_DIRAC
!        this implements arbitrary order pert_scf_eq under two conditions:
!        1. S is perturbation independent
!        2. there is only one component/perturbation
         do i = 1, size(p)
!           check that S is perturbation independent
            if (bas(i)) then
               call quit('pert_scf_eq: size(p) > 2 and perturbation dependent S not implemented')
            end if
!           check that there is only one component/perturbation
            if (dimp(i) > 1) then
               call quit('pert_scf_eq: size(p) > 2 and more than one component/perturbation not implemented')
            end if
         end do

!radovan:
!        this has to be done in this order
!        otherwise gamma is always zero in d2h
!        no idea why
         DSDp(1) = Dp(nd)*S0*Dp(1) + Dp(1)*S0*Dp(nd) - Dp(nd)
         do i = 2, nd - 1
            j = nd - i + 1
            DSDp(1) = DSDp(1) + Dp(j)*S0*Dp(i)
         end do

         X = -sum(w/2.0d0)*S0

         FDSp(1) = (Fp(1) - X)*Dp(nd)
         SDFp(1) =             Dp(nd)*(Fp(1) + X)
         do i = 1, nd
            j = nd - i + 1
            FDSp(1) = FDSp(1) + Fp(i)*Dp(j)
            SDFp(1) = SDFp(1) + Dp(i)*Fp(j)
         end do

         FDSp(1) = FDSp(1)*S0 - S0*SDFp(1)

         SDFp(1) = 0
         X       = 0

#else
         call quit('pert_scf_eq error: Too high order, size(p) > 2 not implemented')
#endif /* ifdef PRG_DIRAC */
      end if
   end subroutine


!radovan: please do not leave commented-out code
!         (without describing the reason why it is commented out)
!         in the repository

!    subroutine pert_scf_eq(p, dp, S0, D, F, DSD, FDS, comp, freq)
!    !p'th derivative of idempotency DSD-D and SCF equation FDS-SDF
!       character(*),    intent(in) :: p(:) !perturbation lables p(n)
!       integer,         intent(in) :: dp(:) !dp(n), dimension for each p
!       type(matrix),    intent(in) :: S0, D(:) !overlap, un- and perturbed density matrices
!       type(matrix),    intent(in) :: F(:) !perturbed Fock matrices
!       type(matrix), intent(inout) :: DSD(*), FDS(*) !(*) to permit re-ranking
!       integer,    optional, intent(in) :: comp(:) !range within each p, c(2*n) (lower,upper indices)
!       complex(8), optional, intent(in) :: freq(:) !signed complex frequency of each p
!       integer      :: n, cc(2*size(p)), b(size(p)), op(size(p)), ne
!       integer      :: i, j, k, i0, j0, k0, ij, ik, jk, ij0, ik0, jk0, ijk, ijk0, ji, jj
!       complex(8)   :: ww(size(p)+1)
!       type(matrix) :: Sp(3), Tp(3)
!       call prop_check('pert_scf_eq', S0, p, c=comp, w=freq, np=n, op=op, cc=cc, ww=ww, ne=ne)
!       ww(n+1) = 0
!       if (any(op/=(/(i,i=1,n)/))) &
!          call quit('pert_scf_eq: General case not implemented')
!       b = cc(n+1:n+n) - cc(1:n) + 1
!       if (size(D) /= product(b+1)) &
!          call quit('pert_scf_eq: size(D) /= product(b+1)' &
!                    // ' wrong number of density matrices')
!       if (size(F) /= product(b+1)) &
!          call quit('pert_scf_eq: size(F)/=product(b+1)' &
!                    // ' wrong number of Fock matrices')
!       !ajt fixme: This should be written optimized an for general orders
!       if (n==1 .and. all( p=='EL' .or. p=='MAGO' )) then
!          do i = 1, b(1)
!             DSD(i) = D(1+i)*S0*D(1) + D(1)*S0*D(1+i) - D(1+i)
!             FDS(i) = F(1+i)*D(1)*S0 - S0*D(1)*F(1+i) &
!                    + F(1)*D(1+i)*S0 - S0*D(1+i)*F(1) - ww(1)*S0*D(1+i)*S0
!          end do
!       else if (n==2 .and. all( p=='EL' .or. p=='MAGO' )) then
!          do i0 = 1, b(1)
!             do j0 = 1, b(2)
!                ij0 = i0+b(1)*(j0-1)
!                i = 1+i0
!                j = 1+b(1)+j0
!                ij = 1+b(1)+b(2)+ij0
!                DSD(ij0) = D(ij)*S0*D(1)+D(1)*S0*D(ij)-D(ij)+D(j)*S0*D(i)+D(i)*S0*D(j)
!                FDS(ij0) = (F(ij)*D(1)+F(j)*D(i)+F(i)*D(j)+(F(1)-sum(ww/2)*S0)*D(ij))*S0 &
!                         - S0*(D(1)*F(ij)+D(i)*F(j)+D(j)*F(i)+D(ij)*(F(1)+sum(ww/2)*S0))
!             end do
!          end do
!       else if (n==3 .and. all( p=='EL' .or. p=='MAGO' )) then
!          do i0 = 1, b(1)
!             do j0 = 1, b(2)
!                do k0 = 1,b(3)
!                   ij0 = i0+b(1)*(j0-1)
!                   ik0 = i0+b(1)*(k0-1)
!                   jk0 = j0+b(2)*(k0-1)
!                   ijk0 = i0+b(1)*(j0-1)+b(1)*b(2)*(k0-1)
!                   i = 1+i0
!                   j = 1+b(1)+j0
!                   k = 1+b(1)+b(2)+k0
!                   ij = 1+b(1)+b(2)+b(3)+ij0
!                   ik = 1+b(1)+b(2)+b(3)+b(1)*b(2)+ik0
!                   jk = 1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+jk0
!                   ijk = 1+b(1)+b(2)+b(3)+b(1)*b(2)+b(1)*b(3)+b(2)*b(3)+ijk0
!                   DSD(ijk0) = D(ijk)*S0*D(1) + D(1)*S0*D(ijk) - D(ijk) &
!                             + D(jk)*S0*D(i) + D(i)*S0*D(jk) &
!                             + D(ik)*S0*D(j) + D(j)*S0*D(ik) &
!                             + D(ij)*S0*D(k) + D(k)*S0*D(ij)
!                   FDS(ijk0) = ( F(ijk)*D(1) + F(jk)*D(i) + F(ik)*D(j)    &
!                                + F(ij)*D(k) + F(k)*D(ij) + F(j)*D(ik)    &
!                                + F(i)*D(jk) + (F(1)-sum(ww/2)*S0)*D(ijk) ) * S0 &
!                             - S0 * ( D(1)*F(ijk) + D(i)*F(jk) + D(j)*F(ik)   &
!                                    + D(k)*F(ij) + D(ij)*F(k) + D(ik)*F(j)    &
!                                    + D(jk)*F(i) + D(ijk)*(F(1)+sum(ww/2)*S0) )
!                end do
!             end do
!          end do
!       elseif (n==1) then
!          do i = 0, b(1)-1
!             if (iszero(D(2+i))) then
!                DSD(1+i) = 0d0*D(1)
!                FDS(1+i) = 0d0*F(1)
!             else
!                call prop_oneint(S0, p, (/1/), S=Sp(1:1), comp=(/cc(1)+i/))
!                DSD(1+i) = D(2+i)*(S0*D(1)) + (D(1)*S0)*D(2+i) &
!                         + D(1)*Sp(1)*D(1)  - D(2+i)
!                FDS(1+i) = F(2+i)*(D(1)*S0) - (S0*D(1))*F(2+i) &
!                         + (F(1)-ww(1)/2*S0) * (D(2+i)*S0 + D(1)*Sp(1)) &
!                         - (Sp(1)*D(1) + S0*D(2+i)) * (F(1) + ww(1)/2*S0)
!                Sp(1) = 0
!             end if
!          end do
!       elseif (n==2) then
!          do j = 0, b(2)-1
!             do i = 1, b(1)
!                jj = 1+b(1)+j+1
!                ij = i+b(1)*j
!                ji = 1+b(1)+b(2)+ij
!                call prop_overlap(p(1:1), S0, Sp(1), c=(/cc(1)+i-1/) )
!                call prop_overlap(p(1:1), S0, Tp(1), c=(/cc(1)+i-1/), w=(/-ww(1)/2,(0d0,0d0)/) )
!                call prop_overlap(p(2:2), S0, Sp(2), c=(/cc(2)+j/) )
!                call prop_overlap(p(2:2), S0, Tp(2), c=(/cc(2)+j/), w=(/-ww(2)/2,(0d0,0d0)/) )
!                call prop_overlap(p(1:2), S0, Sp(3), c=(/cc(1)+i-1,cc(2)+j/) )
!                call prop_overlap(p(1:2), S0, Tp(3), c=(/cc(1)+i-1,cc(2)+j/), &
!                                  w = (/ -ww(1)/2, -ww(2)/2, (0d0,0d0) /) )
!                DSD(ij) = D(1)*Sp(3)*D(1) + D(1+i)*Sp(2)*D(1) + D(1)*Sp(2)*D(1+i) &
!                        + (D(1+i)*S0 + D(1)*Sp(1)) * D(jj) &
!                        + D(jj) * (S0*D(1+i) + Sp(1)*D(1))
!                if (.not.iszero(D(ji))) &
!                   DSD(ij) = DSD(ij) + D(ji)*S0*D(1) + D(1)*S0*D(ji) - D(ji)
!                FDS(ij) = (F(ji) + Tp(3) - (ww(1)+ww(2))/2 * Sp(3)) * D(1)*S0          &
!                        - S0*D(1) * (F(ji) + Tp(3) + (ww(1)+ww(2))/2 * Sp(3))          &
!                        + F(1)*(D(1)*Sp(3)) - (Sp(3)*D(1))*F(1)                        &
!                        + (F(jj) + Tp(2) - ww(2)/2 * Sp(2)) * (D(1+i)*S0 + D(1)*Sp(1)) &
!                        - (S0*D(1+i) + Sp(1)*D(1)) * (F(jj) + Tp(2) + ww(2)/2 * Sp(2)) &
!                        + (F(1+i) + Tp(1) - ww(1)/2 * Sp(1)) * D(jj)*S0                &
!                        - S0*D(jj) * (F(1+i) + Tp(1) + ww(1)/2 * Sp(1))
!                if (.not.iszero(D(ji))) &
!                   FDS(ij) = FDS(ij) + (F(1) - (ww(1)+ww(2))/2 * S0) * D(ji)*S0 &
!                                      - S0*D(ji) * (F(1) + (ww(1)+ww(2))/2 * S0)
!                if (.not.iszero(Sp(2))) &
!                   FDS(ij) = FDS(ij) &
!                        + ( (F(1+i) + Tp(1) - ww(1)/2 * Sp(1)) * D(1)         &
!                          + (F(1) - ww(1) * S0) * D(1+i) ) * Sp(2)            &
!                        - Sp(2) * ( D(1) * (F(1+i) + Tp(1) + ww(1)/2 * Sp(1)) &
!                                  + D(1+i) * (F(1) + ww(1) *S0) )
!                if (.not.iszero(Sp(1))) &
!                   FDS(ij) = FDS(ij) &
!                        + (F(1)-ww(2)*S0) * D(jj)*Sp(1) &
!                        - Sp(1)*D(jj) * (F(1)+ww(2)*S0)
!                Sp = 0
!                Tp = 0
!             end do
!          end do
!       else
!          stop 'PertSCFeq: General case not implemented yet'
!       end if
!    end subroutine


! Old arbitrary-order routine
!
!     subroutine PrtScfEqs(S0,p,w,c,D,F,DSD,FDS)
!     !n'th derivative of HF-equation FDS-SDF and constraint DSD-D
!         character(*):: p(:)
!         integer:: c(:)
!         real*8:: w(:)
!         type(MatPtr):: S0,D(0:),F(0:),DS,SD,S
!         type(MatPtr),intent(inout):: DSD,FDS
!         integer:: n,i,j,k,l
!         logical:: pk(size(p))
!         call MatUse(S0); call MatUse(D); call MatUse(F)
!         n = size(p)
!         if (size(c)/=   n) stop 'FDSDSD: size(c)/=size(p)'
!         if (size(w)/=   n) stop 'FDSDSD: size(w)/=size(p)'
!         if (size(D)/=2**n) stop 'FDSDSD: size(D)/=2**size(w)'
!         if (size(F)/=2**n) stop 'FDSDSD: size(F)/=2**size(w)'
!         do i=0,2**n-1 !i=der of F (and D in DS(D))
!             j = 2**n-1-i !j=der of DS and SD
!             DS = D(j)*S0; SD = S0*D(j)
!             do k = 1,j !k=der of S
!                 if (iand(k,j)/=k) cycle
!                 pk = (/(iand(k,2**l)/=0,l=0,n-1)/)
!                 S = PS(pack(p,pk),pack(c,pk),S0)
!                 DS = DS + D(j-k)*S
!                 SD = SD + S*D(j-k)
!             enddo
!             if (i==0) then
!                 FDS = (F(0)+(sum(w)/2)*S0)*DS - SD*(F(0)-(sum(w)/2)*S0)
!                 DSD = DS*D(0) - D(j)
!                 cycle
!             endif
!             pk = (/(iand(i,2**l)/=0,l=0,n-1)/)
!             S = PS(pack(p,pk),pack(c,pk),S0,w=(/sum(w,.not.pk),pack(w,pk)/))
!             FDS = FDS + (F(i)+S/2)*DS - SD*(F(i)-Trps(S/2))
!             DSD = DSD + DS*D(i)
!             call MatRel(DS,SD,S)
!         enddo
!         if (sum(w)==0) then
!             FDS = AntiPart(FDS)
!             DSD = SymPart(DSD)
!         endif
!         call MatRel(S0); call MatRel(D); call MatRel(F)
!     end subroutine


   !> Solve the p-perturbed response equations for Dp and Fp.
   !> If already stored in the cache (solved_eqs), the solution is fetched
   !> from there. Otherwise, the RHSs are computed before calling the solver
   subroutine pert_dens(mol, S0, p, dimp, D, F, Dp, Fp, comp, freq)
      !> mol/basis/decomp/thresh needed by integrals and solver
      type(prop_molcfg), intent(in) :: mol
      !> unperturbed overlap matrix
      type(matrix), intent(in) :: S0
      !> perturbation lables
      character(*), intent(in) :: p(:)
      !> dimension of property (each perturbation)
      integer,      intent(in) :: dimp(:)
      !> lower-than-p'th order density matrices
      type(matrix), intent(in) :: D(:)
      !> lower-than-p'th order Fock matrices
      type(matrix), intent(in) :: F(:)
      !-------------------------------
      !> resulting p-perturbed density. Deferred shape (*) do permit any rank
      type(matrix), intent(inout) :: Dp(*)
      !> resulting p-perturbed Fock. Deferred shape (*) do permit any rank
      type(matrix), intent(inout) :: Fp(*)
      !-----------------------------------
      !> starting component for each p. Default 1 1 ... 1
      integer,    optional, intent(in) :: comp(:)
      !> frequency of each p. Default all zero (static)
      complex(8), optional, intent(in) :: freq(:)
      !------------------------------------------
      type(matrix), dimension(product(dimp)) :: DSDp, FDSp, Sp
      type(cached_sol), pointer :: sol
      integer    :: ccomp(size(p)), pd, pd1, sym, nanti, i, j
      complex(8) :: ffreq(size(p)), freqdiff

      if (size(dimp) /= size(p)) call quit('pert_dens: Differing numbers ' &
               // 'of perturbations and dimensions, size(dimp) /= size(p)')
      pd  = product(dimp)
      pd1 = product(dimp+1)
      if (size(D) /= pd1-pd) call quit('pert_dens: Wrong number of lower-order ' &
               // 'perturbed densities, size(D) /= product(dimp+1) - product(dimp)')
      if (size(F) /= pd1-pd) call quit('pert_dens: Wrong number of lower-order ' &
               // 'perturbed densities, size(D) /= product(dimp+1) - product(dimp)')
      ccomp = 1
      if (present(comp)) then
         ccomp = comp !ajt fixme Check number and bounds
      end if
      ffreq = 0
      if (present(freq)) then
         ffreq = freq !ajt fixme Check number
      end if
      ! determine the symmetry of the right-hand-sides
      ! ajt fixme: Using sym, pert_scf_eq can be optimized
      nanti = count(pert_antisym(p))
      sym = merge(merge(1, -1, mod(nanti,2)==1), 0, all(ffreq==0) &
                  .or. (size(p)==1 .and. .not.any(pert_basdep(p))))
      ! search through solved_eqs for an already calculated solution
      do i = 1, num_solved_eqs
         sol => solved_eqs(i)
         if (sol%order /= size(p)) cycle !wrong equation order
         if (any(sol%pert(:size(p))  /= p     .or. &  !wrong perturbations
                 sol%ncomp(:size(p)) /= dimp  .or. &  !wrong shape
                 sol%scomp(:size(p)) /= ccomp)) cycle !wrong start comp
         if (.not.all(sol%freq(:size(p)) ==  ffreq) .and. & !wrong
             .not.all(sol%freq(:size(p)) == -ffreq)) cycle  !frequencies
         ! solution exists. If same freq, just copy, 
         do j = 1, pd  !if -freq copy and +-transpose
            if (all(sol%freq(:size(p)) == ffreq)) then
               Dp(j) = sol%D(j)
            else !copy and +-transpose
               Dp(j) = (-1d0)**nanti * dag(sol%D(j))
            end if
         end do
         ! calculate the corresponding Fock matrix
         call pert_fock(mol, p, dimp, (/D,Dp(:pd)/), Fp(:pd), &
                        comp=ccomp, freq=ffreq)
         return
      end do
      ! if an 'EXCI' equation was not caught in the previous
      ! loop, the corresponding excitation density has not
      ! been saved in solved_eqs(:), so quit
      if (size(p)==1 .and. p(1)=='EXCI') &
         call quit('rsp_equations/pert_dens error: EXCI specified but no ' &
                // 'excitation densities saved in rsp_equations')
      ! zero highest-order D and F for RHS calculation
      do i = 1, pd
         Dp(i) = 0d0*D(1)
         Fp(i) = 0d0*F(1)
      end do
      ! if first order FREQ equation, result is zero (so just return)
      if (size(p)==1 .and. p(1)=='FREQ') return
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
!radovan: this is necessary otherwise matrices inherit wrong irep from Sp
!         which breaks point group symmetry in DIRAC
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
         ! KK/AJT quick fix. Avoid solving the same response equation twice for
         ! second-order response equations (i.e. when the two perturbations are identical).
         SecondOrderEq: if(size(p) == 2) then
           freqdiff = ffreq(1) - ffreq(2)
           ! If the two perturbations are identical
           CopySolutionAlreadyDone: if(dimp(1) == 3 .and. dimp(2) == 3 &
                                  & .and. p(1)==p(2) .and. abs(freqdiff) < 1e-9 ) then
              ! Component 2 equals component 4 (similarly for 3,7 and 6,8)
              if(i==4) then
                 Dp(i) = Dp(2)
                 Fp(i) = Fp(2)
                 cycle
              end if
              if(i==7) then
                 Dp(i) = Dp(3)
                 Fp(i) = Fp(3)
                 cycle
              end if
              if(i==8) then
                 Dp(i) = Dp(6)
                 Fp(i) = Fp(6)
                 cycle
              end if
           end if CopySolutionAlreadyDone
        end if SecondOrderEq

         ! If no symmetry can be exploited, calculate response parameter
         call solve_scf_eq(mol, S0, D(1), F(1), sym, sum(ffreq), 1, &
                           DSDp(i:i), FDSp(i:i), Dp(i:i), Fp(i:i))

      end do
      print* !ajt Blank line after prints from solve_scf_eq
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
      !       in:          --------------------------------------------
      !       out:                                                      ------
#ifdef PRG_DIRAC
      use dirac_interface
      use aor_cfg
#endif
#ifdef DALTON_AO_RSP
      use dalton_ifc
#endif
#ifdef VAR_LINSCA
      use RSPsolver,     only: rsp_init, rsp_solver
      use complexsolver, only: rsp_complex_init, rsp_complex_solver
      use scf_config,    only: cfg_rsp_complex, cfg_rsp_gamma
      use decompMod,     only: decompItem
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
   ! Look through solved_eqs for whether freq is an excitation energy,
   ! and return then index of that entry. If not, return -1
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
   ! Project the excitation densities out of the -right-hand-sides FDSp,
   ! and return the corresponding excitation vector Vexci and
   ! transition moments tsmom. tsmom(1,:) are de-excitation moments (-exci),
   ! while tsmom(2,:) are excitation moments (+exci).
   ! ajt may10 Added forgotten projection of Fp
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
#ifndef VAR_LINSCA
         ! remove both resonant and opposite-resonant component
         ! from FDSp. This preserves symmetry. The solution's
         ! opposite-resonant component will be added after solving.
         ! Not activated in LINSCA, as the solver takes care of this
         FDSp(i) = FDSp(i) - tsmom(1,i) * dag(SDxS) &
                           - tsmom(2,i) * SDxS
#endif
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
   ! For an equation which was resonant (excitation: freq=exci,
   ! de-excitation: freq=-exci), and is being solved for the
   ! non-resonant part, add the contribution to the solution Dp
   ! from the opposite (de-)excitation density (Dx or Dx^T).
   ! The factor before Dx or Dx^T is given by the transition moment
   ! tsmom(2 or 1,:), divided by twice the frequency
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
! In LINSCA, the solver takes care of this
#ifndef VAR_LINSCA
            Dp(i) = Dp(i) + (tsmom(1,i) / (freq+exci)) * dag(Dx)
#endif
         else !if de-excitation resonant, add excitation contrib a*Dx
#ifndef VAR_LINSCA
            Dp(i) = Dp(i) + (tsmom(2,i) / (freq-exci)) * Dx
#endif
         end if
      end do
   end subroutine
#endif /* not PRG_DIRAC */


end module
