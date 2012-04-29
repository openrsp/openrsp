!> @file
!> Contains module prop_test

!> ajt/radovan: Response-related testing routines and some calculations 
module prop_test_old

   use matrix_defop
   use prop_contribs_old
   use rsp_equations_old

   implicit none

   public test_vcd
   public test_mcd
   public elec_dipmom
   public elec_polariz
   public elec_hypolar
   public alt_elec_hypol
   public elec_sechyp
   public alt_elec_sechyp
   public alt2_elec_sechyp
   public magnetiz
   public elec_quadrupole
   public prop_test_gradient

   private

   !physical constants for conversion to non-atomic units
   real(8), parameter:: cm1 = 1/219474.631371d0, & !1 centimeter-to-minus-one in au
                        cvl = 137.03599907d0,    & !c-the-velocity-of-light in au
                        nm  = 10/0.52917706d0,   & !1 nanometer in au
                        pi  = 3.14159265358979323846D0 !acos(-1.d0)

   !field component lables for printing
   character*2:: fc(3) = (/'Fx','Fy','Fz'/), &
                 bc(3) = (/'Bx','By','Bz'/)

contains

   subroutine test_vcd(mol, ng, S, D, F)
   ! "Atomic axial tensor" = -1/2 * Egbw
   ! Egbw = << GEO(0.0) MAG(0.0) FREQ(91) >>
      type(prop_molcfg), intent(in) :: mol
      integer,      intent(in) :: ng
      type(matrix), intent(in) :: S, D, F
      complex(8)    :: Egbw(ng,3), Rn3(ng,3), Rm3(ng,3)
      type(matrix)  :: Db(3), Fb(3), Dbw(3), Fbw(3)
      type(matrix)  :: M3(3), Sb(3), FDSbw(3), DFDbw(3)
      integer       :: j
      call pert_dens(mol, S, (/'MAG'/), (/3/), &
                     (/D/), (/F/), Db, Fb, freq=(/(0d0,0d0)/))
      ! hack: Obtain frequency derivative of one-electron integrals
      !       by HUGE finite difference
      call prop_oneint(mol, S, (/'MAG'/), (/3/), F=Fbw, S=Sb, &
                       freq = (huge(1d0)/128) * (/(1d0,0d0)/))
      call prop_oneint(mol, S, (/'MAG'/), (/3/), F=M3, &
                       freq = (/(0d0,0d0)/))
      do j=1,3
         Fbw(j) = (128/huge(1d0)) * (Fbw(j)-M3(j))
         M3(j) = 0
      end do

      ! contstruct -RHS's for the frequency-differentiated equation
      do j=1,3
         FDSbw(j) = -1/2d0 * S*D*Sb(j) - S*Db(j)*S &
                  -  1/2d0 * Sb(j)*D*S
         Sb(j) = 0
      end do
      ! call solver directly
      call solve_scf_eq(mol, S, D, F, -1, (0d0,0d0), 3, &
                        Sb, FDSbw, Dbw, Fbw)
      ! contract the frequency-differentiated response function
      Egbw = 0
      call prop_oneave(mol, S, (/'GEO','MAG'/), (/D/), shape(Egbw), Egbw, &
                       freq = (/(-1d0,0d0), (1d0,0d0)/))
      ! call print_tensor(shape(Egbw), Egbw, '-i/2TgbwD'); Egbw=0
      ! hack: Obtain frequency-differentiated one-electron integrals
      !       by another HUGE finite difference
      Rn3 = 0
      call prop_oneave(mol, S, (/'GEO'/), (/Db/), shape(Rn3), Rn3, &
                       freq = (huge(1d0)/128) * (/(-1d0,0d0)/))
      Rm3 = 0
      call prop_oneave(mol, S, (/'GEO'/), (/Db/), shape(Rm3), Rm3, &
                       freq = (/(0d0,0d0)/))
      Egbw = Egbw + (128/huge(1d0)) * (Rn3 - Rm3)
      ! call print_tensor(shape(Egbw), Egbw, '-i/2TgwDb'); Egbw=0
      do j=1,3
         DFDbw(j) = Dbw(j)*F*D + D*Fbw(j)*D + D*F*Dbw(j) &
                  + 1/2d0 * Db(j)*S*D - 1/2d0 * D*S*Db(j)
      end do
      call prop_oneave(mol, S, (/'GEO'/), (/Dbw/), shape(Egbw), Egbw, &
                       DFD = DFDbw)
      DFDbw = 0
      ! call print_tensor(shape(Egbw), Egbw, 'HgDbw - SgDFDbw'); Egbw=0
      call prop_twoave(mol, (/'GEO'/), (/D,Dbw/), shape(Egbw), Egbw)
      ! call print_tensor(shape(Egbw), Egbw, 'Gg(D)Dbw'); Egbw=0
      ! print, free, return
      call print_tensor(shape(Egbw), Egbw, 'Egbw')
      Db=0; Fb=0; Dbw=0; Fbw=0
   end subroutine



   subroutine test_mcd(mol, S, D, F)
      type(prop_molcfg), intent(in) :: mol
      type(matrix),      intent(in) :: S, D, F
      complex(8)   :: exci(2), Ebfx(3,3,2)
      type(matrix) :: Df(3), Dx(2), Dfx(3,2), DFDfx(3,2)
      type(matrix) :: Ff(3), Fx(2), Ffx(3,2)
      integer      :: i, j, k
      ! find degenerate excitation
      do i = num_solved_eqs-1, 0, -1
         if (i == 0) exit
         if (solved_eqs(i  )%order /= 1) cycle
         if (solved_eqs(i+1)%order /= 1) cycle
         if (solved_eqs(i  )%pert(1) /= 'EXCI') cycle
         if (solved_eqs(i+1)%pert(1) /= 'EXCI') cycle
         if (abs(solved_eqs(i  )%freq(1) &
               - solved_eqs(i+1)%freq(1)) >= 1d-5) cycle
         exci(2) = -solved_eqs(i  )%freq(1)
         exci(1) = -solved_eqs(i+1)%freq(1)
         exit
      end do
      if (i==0) then
         print *, 'no degenerate excited states found'
         return
      end if
      print *, 'found degenerate excitations: ', real(exci(1)), real(exci(2))
      call pert_dens(mol, S, (/'EL'/), (/3/), &
                     (/D/), (/F/), Df, Ff, freq=(/-(exci(1)+exci(2))/2/))
      call pert_dens(mol, S, (/'EXCI'/), (/1/), &
                     (/D/), (/F/), Dx(1:1), Fx(1:1), freq=exci(1:1))
      call pert_dens(mol, S, (/'EXCI'/), (/1/), &
                     (/D/), (/F/), Dx(2:2), Fx(2:2), freq=exci(2:2))
      call pert_dens(mol, S, (/'EL  ','EXCI'/), (/3,1/), &
                     (/D,Df,Dx(1:1)/), (/F,Ff,Fx(1:1)/),    &
                     Dfx(:,1), Ffx(:,1), freq=(/-exci(1),exci(1)/))
      call pert_dens(mol, S, (/'EL  ','EXCI'/), (/3,1/), &
                     (/D,Df,Dx(2:2)/), (/F,Ff,Fx(2:2)/),    &
                     Dfx(:,2), Ffx(:,2), freq=(/-exci(2),exci(2)/))
      ! contract no-London response function
      Ebfx(:,:,:) = 0
      call prop_oneave(mol, S, (/'MAGO'/), (/Dfx/), shape(Ebfx), Ebfx)
      ! call print_tensor(shape(Ebfx), Ebfx, 'E1bDfx - i/2TbDfx - SbDFDfx'); Ebfx=0
      call print_tensor(shape(Ebfx), -Ebfx, 'London d/db alphax = -Ebfx', &
                        (/0*exci(1),-exci(1),exci(1)/))
      ! contract London response function
      Ebfx(:,:,:) = 0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Dx/), shape(Ebfx), Ebfx)
      ! call print_tensor(shape(Ebfx), Ebfx, 'E1fbDx'); Ebfx=0
      call prop_twoave(mol, (/'MAG'/), (/D,Df,Dx,Dfx/), shape(Ebfx), Ebfx)
      ! call print_tensor(shape(Ebfx), Ebfx, 'E2bDfDx+E2bD0Dfx'); Ebfx=0
      do k = 1, 2
         do j = 1, 3
            DFDfx(j,k) = D*Ff(j)*Dx(k) + Df(j)*Fx(k)*D &
                       + Dx(k)*Ff(j)*D + D*Fx(k)*Df(j) &
                       + Df(j)*(F+(-2*exci(k))/2*S)*Dx(k) &
                       + Dx(k)*(F-(-2*exci(k))/2*S)*Df(j) &
                       + Dfx(j,k)*(F+(0*exci(k))/2*S)*D + D*Ffx(j,k)*D &
                       + D*(F-(0*exci(k))/2*S)*Dfx(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'MAG'/), (/Dfx/), shape(Ebfx), Ebfx, &
                       freq=(/(0d0,0d0)/), DFD=(/DFDfx/))
      DFDfx = 0
      ! call print_tensor(shape(Ebfx), Ebfx, 'E1bDfx - i/2TbDfx - SbDFDfx'); Ebfx=0
      call print_tensor(shape(Ebfx), -Ebfx, 'London d/db alphax = -Ebfx', &
                        (/0*exci(1),-exci(1),exci(1)/))
      Df=0; Ff=0; Dx=0; Fx=0; Dfx=0; Ffx=0
   end subroutine



   subroutine prop_test_gradient(mol, S, D, F, ng)
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in) :: S, D, F
      integer,      intent(in) :: ng !number of geometrical coordinates
      complex(8)   :: Eg(ng) !NB: zero first
      type(matrix) :: DFD

      ! contract first-order geometry-differentiated integrals
      ! with unperturbed density D and energy-weighted density DFD
      DFD = D*F*D
      Eg = 0 !zero first, since prop_one/twoave are incremental
      call prop_oneave(mol, S, (/'GEO'/), (/D/), (/ng/), Eg, DFD=(/DFD/))
      call prop_twoave(mol, (/'GEO'/), (/D/), (/ng/), Eg)
      ! print
      call print_tensor((/ng/), Eg, 'gradient = Eg = E0g - Sg DFD')
      ! free DFD
      DFD = 0

   end subroutine

   subroutine prop_test_Gprime_Df(mol, S, D, F, idum, freq) !Gprime_Df(-freq)
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in) :: S, D, F
      integer,      intent(in) :: idum    !not used
      complex(8),   intent(in) :: freq(:) !frequency
      complex(8)   :: Efb(3,3) = 0
      type(matrix) :: Df(3), Ff(3), DFDf(3)
      integer      :: i, j
      ! the following contribution is zero
      ! call prop_oneave(mol, S, (/'EL ','MAG'/), (/D/), shape(Efb), Efb)
      ! solve EL response equations for the negative frequency
      call pert_dens(mol, S, (/'EL'/), shape(Df), (/D/), (/F/), &
                     Df, Ff, freq=-freq(1:1))
      ! compute the EL-perturbed energy-weighted density matrices,
      ! to be contracted against MAG-perturbed overlap integrals
      do i = 1, 3
         DFDf(i) = Df(i)*(F - freq(1)/2*S)*D + D*Ff(i)*D &
                     + D*(F + freq(1)/2*S)*Df(i)
      end do
      ! contract Df and DFDf with MAG-perturbed integrals
      call prop_oneave(mol, S, (/'MAG'/), (/Df/), shape(Efb), Efb, &
                       perm=(/2,1/), freq=freq(1:1), DFD=DFDf)
      call prop_twoave(mol, (/'MAG'/), (/D,Df/), shape(Efb), Efb, perm=(/2,1/))
      ! print
      call print_tensor(shape(Efb), Efb, 'G-prime = Efb = E1b Df - Sb DFDf')
      ! free Df(:), Ff(:), DFDf(:)
      Df=0; Ff=0; DFDf=0
   end subroutine


   subroutine prop_test_Gprime_Db(mol, S, D, F, idum, freq) !Gprime_Db(+freq)
      type(prop_molcfg), intent(in) :: mol
      type(matrix),      intent(in) :: S, D, F
      integer,           intent(in) :: idum    !not used
      complex(8),        intent(in) :: freq(:) !frequency
      complex(8)   :: Efb(3,3) = 0
      type(matrix) :: Db(3), Fb(3)
      integer      :: i
      ! the following contribution is zero
      call prop_oneave(mol, S, (/'EL ','MAG'/), (/D/), shape(Efb), Efb)
      ! solve EL response equations for the negative frequency
      call pert_dens(mol, S, (/'MAG'/), shape(Db), (/D/), (/F/), &
                     Db, Fb, freq=freq(1:1))
      ! contract Db with -dipole integrals
      call prop_oneave(mol, S, (/'EL'/), (/Db/), shape(Efb), Efb)
      ! print
      call print_tensor(shape(Efb), Efb, 'G-prime = Efb = E1f Db')
      ! free
      Db = 0;  Fb = 0
   end subroutine


   subroutine elec_dipmom(mol, S, D, F, Ef)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(out) :: Ef(3)
      Ef = 0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), Ef)
      call print_tensor((/3/), -Ef, 'mu = -Ef', unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3/), -Ef, 'mu = -Ef')
   end subroutine


   subroutine elec_quadrupole(mol, S, D, F, t)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(out) :: t(6)
      t = 0.0d0
      call prop_oneave(mol, S, (/'ELGR'/), (/D/), (/6/), t)
      call print_tensor((/6/), -t, 'quadrupole = -t')
   end subroutine


   subroutine elec_polariz(mol, S, D, F, freq, Eff)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq
      complex(8),        intent(out) :: Eff(3,3)
      type(matrix) :: Df(3), Ff(3)
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq/))
      Eff = 0 !zero, as prop_oneave works incrementally
      call prop_oneave(mol, S, (/'EL'/), (/Df/), (/3,3/), Eff)
      call print_tensor((/3,3/), -Eff, 'alpha = -Eff', unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3/), -Eff, 'alpha = -Eff')
      !delete response matrices
      Df = 0; Ff = 0
      !change sign, since polarizability is minus quasi-energy derivative
      Eff = -Eff
   end subroutine


   subroutine elec_hypolar(mol, S, D, F, freq, Efff)
      type(prop_molcfg), intent(in)   :: mol
      type(matrix),      intent(in)   :: S, D, F
      complex(8),        intent(in)   :: freq(3)
      complex(8),        intent(out)  :: Efff(3,3,3)
      type(matrix) :: Df(3), De(3), Dfe(3,3)
      type(matrix) :: Ff(3), Fe(3), Ffe(3,3)
      integer      :: i, j
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/elec_hypolar: sum(freq) should be zero!',-1)
      !solve equations
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=freq(2:2))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=freq(3:3))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Df,De/), (/F,Ff,Fe/), &
                     Dfe, Ffe, freq=(/freq(2),freq(3)/))
      !just one energy contribution, no reort or multiplier contributions
      Efff=0
      call prop_oneave(mol, S, (/'EL'/), (/Dfe/), (/3,3,3/), Efff)
      !print
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq, unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq)
      !free
      Df=0; Ff=0; De=0; Fe=0; Dfe=0; Ffe=0
      !change sign to get beta
      Efff = -Efff
   end subroutine


   subroutine alt_elec_hypol(mol, S, D, F, freq, Efff)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq(3)
      complex(8),        intent(out) :: Efff(3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), DeSD(3), FDSfg
      type(matrix) :: Fe(3), Ff(3), Fg(3), FeDS(3), DSDfg, zm
      character(4) :: UNP(0)
      integer      :: i, j, k
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_hypol: sum(w) should be zero!',-1)
      !solve equations
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=freq(1:1))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=freq(2:2))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=freq(3:3))
      !no energy contributions for HF, one for KSDFT
      Efff = 0
      zm = 0*D
      call prop_twoave(mol, UNP, (/D,De,Df,Dg,(zm,i=1,54)/), (/3,3,3/), Efff)
      !gradient Lagrange multiplier contribution -tr (DeSD-) (FDS-)fg
      do i = 1, 3
         DeSD(i) = De(i)*S*D - D*S*De(i)
      end do
      do k = 1, 3
         do j = 1, 3
            FDSfg = (Fe(k)*Df(j) + Ff(j)*De(k)) * S &
              - S * (De(k)*Ff(j) + Df(j)*Fe(k))
            do i = 1, 3
               Efff(i,j,k) = Efff(i,j,k) - tr(DeSD(i),FDSfg)
            end do; FDSfg=0
         end do
      end do; DeSD=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+) (DSD)fg
      do i = 1, 3
         FeDS(i) = Fe(i)*D*S + S*D*Fe(i) - Fe(i)
      end do
      do k = 1, 3
         do j = 1, 3
            DSDfg = De(k)*S*Df(j) + Df(j)*S*De(k)
            do i = 1, 3
               Efff(i,j,k) = Efff(i,j,k) - tr(FeDS(i),DSDfg)
            end do; DSDfg=0
         end do
      end do; FeDS=0
      !print
      call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq, unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3,3/), -Efff, 'beta = -Efff', freq)
      !free
      De=0; Fe=0; Df=0; Ff=0; Dg=0; Fg=0
      !change sign for beta
      Efff = -Efff
   end subroutine


   subroutine elec_sechyp(mol, S, D, F, freq, Effff)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq(4)
      complex(8),        intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Def(3,3),   &
                      Fe(3), Ff(3), Fg(3), Fef(3,3),   &
                      Deg(3,3), Dfg(3,3), Defg(3,3,3), &
                      Feg(3,3), Ffg(3,3), Fefg(3,3,3)
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/elec_sechyp: sum(freq) should be zero!',-1)
      !solve equations
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(2)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(3)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(4)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(2),freq(3)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Dg/), (/F,Fe,Fg/), &
                     Deg, Feg, freq=(/freq(2),freq(4)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dg/), (/F,Ff,Fg/), &
                     Dfg, Ffg, freq=(/freq(3),freq(4)/))
      call pert_dens(mol, S, (/'EL','EL','EL'/), (/3,3,3/), &
                     (/D,De,Df,Dg,Def,Deg,Dfg/), &
                     (/F,Fe,Ff,Fg,Fef,Feg,Ffg/), &
                     Defg, Fefg, freq=(/freq(2),freq(3),freq(4)/))
      !just one property contribution
      Effff=0
      call prop_oneave(mol, S, (/'EL'/), (/Defg/), (/3,3,3,3/), Effff)
      !print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      !free
      De=0; Df=0; Dg=0; Def=0; Deg=0; Dfg=0; Defg=0
      Fe=0; Ff=0; Fg=0; Fef=0; Feg=0; Ffg=0; Fefg=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine


   subroutine alt_elec_sechyp(mol, S, D, F, freq, Effff)
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(4)
      complex(8),   intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Dh(3), DeSD(3), &
                      Dfg(3,3), Dfh(3,3), Dgh(3,3), FDSfgh, &
                      Fe(3), Ff(3), Fg(3), Fh(3), FeDS(3), &
                      Ffg(3,3), Ffh(3,3), Fgh(3,3), DSDfgh, zm
      character(4) :: UNP(0)
      integer      :: i, j, k, l
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_sechyp: sum(w) should be zero!',-1)
      !solve equations
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(1)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(2)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(3)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dh, Fh, freq=(/freq(4)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dg/), (/F,Ff,Fg/), &
                     Dfg, Ffg, freq=(/freq(2),freq(3)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Df,Dh/), (/F,Ff,Fh/), &
                     Dfh, Ffh, freq=(/freq(2),freq(4)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Dg,Dh/), (/F,Fg,Fh/), &
                     Dgh, Fgh, freq=(/freq(3),freq(4)/))
      ! one energy contribution
      Effff=0
      zm = 0*D
      call prop_twoave(mol, UNP, (/D,De,Df,Dg,Dh, &
                              (zm,i=1,18),Dfg,(zm,i=1,9), &
                              Dfh,Dgh,(zm,i=1,189)/), &
                       shape(Effff), Effff)
      ! gradient Lagrange multiplier contribution -tr (DeSD-) (FDS-)fgh
      do i = 1, 3
         DeSD(i) = De(i)*S*D - D*S*De(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               FDSfgh = (Ff(j)*Dgh(k,l) + Fg(k)*Dfh(j,l) &
                      +  Fh(l)*Dfg(j,k) + Ffg(j,k)*Dh(l) &
                      +  Ffh(j,l)*Dg(k) + Fgh(k,l)*Df(j)) * S &
                  - S * (Df(j)*Fgh(k,l) + Dg(k)*Ffh(j,l) &
                      +  Dh(l)*Ffg(j,k) + Dfg(j,k)*Fh(l) &
                      +  Dfh(j,l)*Fg(k) + Dgh(k,l)*Ff(j))
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSD(i),FDSfgh)
               end do; FDSfgh=0
            end do
         end do
      end do; DeSD=0
      ! idempotency Lagrange multiplier contribution -tr (FhDS+) (DSD)efg
      do i = 1, 3
         FeDS(i) = Fe(i)*D*S + S*D*Fe(i) - Fe(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DSDfgh = Df(j)*S*Dgh(k,l) + Dg(k)*S*Dfh(j,l) &
                      + Dh(l)*S*Dfg(j,k) + Dfg(j,k)*S*Dh(l) &
                      + Dfh(j,l)*S*Dg(k) + Dgh(k,l)*S*Df(j)
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDS(i),DSDfgh)
               end do; DSDfgh=0
            end do
         end do
      end do; FeDS=0
      ! print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      !free
      De=0; Df=0; Dg=0; Dh=0; Dfg=0; Dfh=0; Dgh=0
      Fe=0; Ff=0; Fg=0; Fh=0; Ffg=0; Ffh=0; Fgh=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine


   subroutine alt2_elec_sechyp(mol, S, D, F, freq, Effff)
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq(4)
      complex(8),        intent(out) :: Effff(3,3,3,3)
      type(matrix) :: De(3), Df(3), Dg(3), Dh(3), DSD,      &
                      Fe(3), Ff(3), Fg(3), Fh(3), FDS,      &
                      Def(3,3), Deg(3,3), Deh(3,3),         &
                      Fef(3,3), Feg(3,3), Feh(3,3),         &
                      DeSDf(3,3), FDSgh, FeDSf(3,3), DSDgh, &
                      DeSDg(3,3), FDSfh, FeDSg(3,3), DSDfh, &
                      DeSDh(3,3), FDSfg, FeDSh(3,3), DSDfg
      integer      :: i, j, k, l
      if (abs(sum(freq)) > 1d-15) &
         call quit('prop_test/alt_elec_sechyp: sum(freq) should be zero!',-1)
      !solve equations
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(1)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(2)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dg, Fg, freq=(/freq(3)/))
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dh, Fh, freq=(/freq(4)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(1),freq(2)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Dg/), (/F,Fe,Fg/), &
                     Deg, Feg, freq=(/freq(1),freq(3)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Dh/), (/F,Fe,Fh/), &
                     Deh, Feh, freq=(/freq(1),freq(4)/))
      !no energy contributions
      Effff=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)f (FDS-)gh
      do j = 1, 3
         do i = 1, 3
            DeSDf(i,j) = Def(i,j)*S*D  - D*S*Def(i,j)  &
                       + De(i)*S*Df(j) - Df(j)*S*De(i)
         end do
      end do
      do l = 1, 3
         do k = 1, 3
            FDSgh = (Fg(k)*Dh(l) + Fh(l)*Dg(k)) * S &
              - S * (Dh(l)*Fg(k) + Dg(k)*Fh(l))
            do j = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDf(i,j),FDSgh)
               end do
            end do; FDSgh=0
         end do
      end do; DeSDf=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+)f (DSD)gh
      do j = 1, 3
         do i = 1, 3
            FeDSf(i,j) = (Fef(i,j)*D + Fe(i)*Df(j)) * S &
                   + S * (D*Fef(i,j) + Df(j)*Fe(i)) - Fef(i,j)
         end do
      end do
      do l = 1, 3
         do k = 1, 3
            DSDgh = Dg(k)*S*Dh(l) + Dh(l)*S*Dg(k)
            do j = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSf(i,j),DSDgh)
               end do
            end do; DSDgh=0
         end do
      end do; FeDSf=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)g (FDS-)fh
      do k = 1, 3
         do i = 1, 3
            DeSDg(i,k) = Deg(i,k)*S*D  - D*S*Deg(i,k)  &
                       + De(i)*S*Dg(k) - Dg(k)*S*De(i)
         end do
      end do
      do l = 1, 3
         do j = 1, 3
            FDSfh = (Ff(j)*Dh(l) + Fh(l)*Df(j)) * S &
              - S * (Df(j)*Fh(l) + Dh(l)*Ff(j))
            do k = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDg(i,k),FDSfh)
               end do
            end do; FDSfh=0
         end do
      end do; DeSDg=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+)g (DSD)fh
      do k = 1, 3
         do i = 1, 3
            FeDSg(i,k) = (Feg(i,k)*D + Fe(i)*Dg(k)) * S &
                   + S * (D*Feg(i,k) + Dg(k)*Fe(i)) - Feg(i,k)
         end do
      end do
      do l = 1, 3
         do j = 1, 3
            DSDfh = Df(j)*S*Dh(l) + Dh(l)*S*Df(j)
            do k = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSg(i,k),DSDfh)
               end do
            end do; DSDfh=0
         end do
      end do; FeDSg=0
      !gradient Lagrange multiplier contribution -tr (DeSD-)h (FDS-)fg
      do l = 1, 3
         do i = 1, 3
            DeSDh(i,l) = Deh(i,l)*S*D  - D*S*Deh(i,l)  &
                       + De(i)*S*Dh(l) - Dh(l)*S*De(i)
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            FDSfg = (Ff(j)*Dg(k) + Fg(k)*Df(j)) * S &
              - S * (Df(j)*Fg(k) + Dg(k)*Ff(j))
            do l = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSDh(i,l),FDSfg)
               end do
            end do; FDSfg=0
         end do
      end do; DeSDh=0
      !idempotency Lagrange multiplier contribution -tr (FeDS+)h (DSD)fg
      do l = 1, 3
         do i = 1, 3
            FeDSh(i,l) = (Feh(i,l)*D + Fe(i)*Dh(l)) * S &
                   + S * (D*Feh(i,l) + Dh(l)*Fe(i)) - Feh(i,l)
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            DSDfg = Df(j)*S*Dg(k) + Dg(k)*S*Df(j)
            do l = 1, 3
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDSh(i,l),DSDfg)
               end do
            end do; DSDfg=0
         end do
      end do; FeDSh=0
      !print
      call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq, unit=mol%lupri)
      if (mol%lupri /= 6) call print_tensor((/3,3,3,3/), -Effff, 'gamma = -Effff', freq)
      !free
      De=0; Df=0; Dg=0; Dh=0; Def=0; Deg=0; Deh=0
      Fe=0; Ff=0; Fg=0; Dh=0; Fef=0; Feg=0; Feh=0
      !change sign to get gamma
      Effff = -Effff
   end subroutine


   subroutine magnetiz(mol, S, D, F, freq, Eoo, Ebb)
   !no-Lon and London magnetizability
      type(prop_molcfg), intent(in)  :: mol
      type(matrix),      intent(in)  :: S, D, F
      complex(8),        intent(in)  :: freq
      complex(8),        intent(out) :: Eoo(3,3), Ebb(3,3)
      type(matrix) :: Db(3), Fb(3), DFD, DFDb(3)
      integer      :: i, j
      !no-London
      call pert_dens(mol, S, (/'MAGO'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq/))
      Eoo = 0
      call prop_oneave(mol, S, (/'MAGO','MAGO'/), (/D/), (/3,3/), Eoo)
      !call print_tensor( (/3,3/), Eoo, 'E0oo'); Eoo=0
      call prop_oneave(mol, S, (/'MAGO'/), (/Db/), (/3,3/), Eoo)
      !call print_tensor( (/3,3/), Eoo, 'E1oDo'); Eoo=0
      call print_tensor((/3,3/), Eoo, 'no-London Magnetizability = Eoo', (/-freq,freq/))
      Db=0; Fb=0
      !London diamag
      call pert_dens(mol, S, (/'MAG'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq/))
      Ebb = 0
      call prop_twoave(mol, (/'MAG','MAG'/), (/D/), (/3,3/), Ebb)
      DFD = D*F*D
      call prop_oneave(mol, S, (/'MAG','MAG'/), (/D/), (/3,3/), Ebb, &
                       DFD=(/DFD/), freq=(/-freq,freq/))
      DFD = 0
      !call print_tensor( (/3,3/), Ebb, 'E0bb-i/2TbbD-SbbW'); Ebb=0
      !London paramag
      call prop_twoave(mol, (/'MAG'/), (/D,Db/), (/3,3/), Ebb)
      do j = 1, 3
         DFDb(j) = Db(j)*(F+(freq/2)*S)*D + D*Fb(j)*D &
                 +     D*(F-(freq/2)*S)*Db(j)
      end do
      call prop_oneave(mol, S, (/'MAG'/), (/Db/), (/3,3/), Ebb, &
                       DFD=DFDb, freq=(/-freq/))
      Db=0; Fb=0; DFDb=0
      !call print_tensor( (/3,3/), Ebb, 'E1bDb-i/2TbDb-SbDFDb'); Ebb=0
      call print_tensor((/3,3/), Ebb, 'London Magnetizability = Ebb', (/-freq,freq/))
   end subroutine


  subroutine test_grcont(n, nr_atoms, test_london, test_geo)

    integer              :: n, nr_atoms
    logical              :: test_london, test_geo
!   ----------------------------------------------------------------------------
    real(8), allocatable :: work(:)
    integer              :: lwork

    integer              :: i, j, k, l, m
    real(8)              :: D(n, n, 2)
    real(8)              :: F(n, n, 3), g(3, nr_atoms)

    lwork = 10000*n + 50*n*n
    allocate(work(lwork))

    if (test_london) then

      print '(a)', 'London 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l', 'Bx', 'By', 'Bz'

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1

          D(:, :, 1) = 0.0d0
          D(k, l, 1) = 1.0d0
          F          = 0.0d0

          call grcont(work, lwork, F, n*n*3, &
                      .false., .true., 1, 0, .false., .true., D(:, :, 1), 1)
          do j = 1, n
            do i = 1, n
              write(0, '(4i3, 3f20.12)') i, j, k, l, F(i, j, :)
            end do
          end do
        end do
      end do
    end if

    if (test_geo) then

      print '(a)', 'Geometry 2e integrals'
      print '(a, 12x, a2, 18x, a2, 18x, a2)', '  i  j  k  l  a', 'Rx', 'Ry', 'Rz'

      D = 0.0d0

!     do l = 1, n
!       do k = 1, n
      do l = 1, 1
        do k = 1, 1
          do j = 1, n
            do i = 1, n

              D(k, l, 1) = 1.0d0
              D(i, j, 2) = 1.0d0
              g          = 0.0d0

              call grcont(work, lwork, g, 3*nr_atoms, &
                          .true., .false., 1, 0, .true., .false., D, 2)
              D(k, l, 1) = 0.0d0
              D(i, j, 2) = 0.0d0

              do m = 1, nr_atoms
                write(0, '(5i3, 3f20.12)') i, j, k, l, m, g(:, m)
              end do
            end do
          end do
        end do
      end do
    end if

    deallocate(work)

  end subroutine


   subroutine print_tensor(dims, tensor, title, freqs, colwidth, unit)
   !ajt This is used for printing response function tensors
   !TODO: Add comp. lables argument. Make space between blocks of rows
      integer,                intent(in) :: dims(:)
      complex(8),             intent(in) :: tensor(*) !(product(dims))
      character(*), optional, intent(in) :: title
      integer,      optional, intent(in) :: colwidth, unit
      complex(8),   optional, intent(in) :: freqs(:)
      integer     :: uni, siz, dec, cwid, i
      character*9 :: fmt
      uni = 6
      if (present(unit)) uni = unit
      if (present(title)) then
         if (present(freqs)) then
            fmt = '(a,nf6.3)'
            write (fmt(4:4),'(i1)') size(freqs)
            write (uni,fmt) title // ' w:', dreal(freqs)
         else
            write (uni,'(a)') title
         end if
      end if
      siz = product(dims)
      dec = max(2, 1+ceiling(log10(maxval(abs(tensor(1:siz))))))
      !dec = max(dec,ceiling(log10( maxval(tensor(1:siz)))))
      cwid = 15
      if (present(colwidth)) cwid = max(colwidth,dec+2)
      dec = cwid - dec - 1
      fmt = '(fww.dd)'
      write (fmt(3:7),'(i2,a1,i2)') cwid, '.', dec
      if (any(dimag(tensor(1:siz)) /= 0)) write (uni,'(a)') '    (real part)'
      call subr(dims, dreal(tensor(1:siz)), uni, cwid, fmt)
      if (any(dimag(tensor(1:siz)) /= 0)) then
         write (uni,'(a)') '    (imaginary part)'
         call subr(dims, dimag(tensor(1:siz)), uni, cwid, fmt)
      end if
      write (uni,'()') !final blank line
   contains
      subroutine subr(dims, tnz, unit, cwid, fmt)
         integer,     intent(in) :: dims(:), unit, cwid
         real(8),     intent(in) :: tnz(*)
         character*8, intent(in) :: fmt
         integer :: i, j, l
         character*(dims(1)*(cwid+1)) :: line
         do i = 1, product(dims), dims(1)
            l = 1
            do j = 1, dims(1)
               line(l:l) = ' '
               write (line(l+1:l+cwid), fmt) tnz(i+j-1)
               l = l + 1 + cwid
            end do
            write (unit,'(a)') line
         end do
      end subroutine
   end subroutine


!ajt Function for generating Gaussian random real numbers,
!    which I occiationally use for debugging
!
!    !> Random real numbers with Gaussian distribution
!    !> The generator is seeded on the first call
!    !> Two reals are generated each time, and returned
!    !> in consecutive calls. Random_Seed and Random_Number
!    !> are provided by Intel and GNU compilers.
!    !> To get different numbers, change seed.
!    function RandomReal() result(r)
!        real(8),parameter :: pi=3.14159265358979323846d00
!        real(8) :: r,a,x,y
!        real(8),save:: rnext
!        logical,save:: first=.true.,flip=.true.
!        if (first) then !reset fortran's random generator
!            call Random_Seed(put=(/628468422,281827387/))
!            first = .false.
!        end if
!        r = rnext
!        flip = .not.flip
!        if (flip) return !return sav
!        call Random_Number(x)
!        a = 2*Pi*x !angle in [0,2pi]
!        call Random_Number(x)
!        y = sqrt(-log(x))
!        rnext = cos(a)*y
!        r = sin(a)*y
!    end function

end module