!> @file
!> Contains module birefring

!> Calculation and outputting of optical birefringences
!> This file is currently common to Dalton and Dirac
module birefring

   use matrix_defop      !matrix type and operators
   use prop_contribs     !property contractions and integrals
   use rsp_equations     !response equation contractor/solver
   implicit none

   ! ajt LSDALTON has replaced the (global) quit with lsquit
   !     with unit (lupri) as extra argument, which doesn't
   !     exist in DIRAC. For now, this macro gets around that.
#ifdef LSDALTON_ONLY
#define quit(msg) lsquit(msg,mol%lupri)
#endif
 
   public efgb_Jpri_Bten_Bcal
   public efgb_output
   public cme_jones_eta_apri
   public cme_output, jones_output
   public verdet_dpol_db

   private

   !wavenumber, velocity of light, and nanometer (in au), pi
   real(8), parameter :: cm1 = 1/219474.631371d0, &
                         cvl = 137.03599907d0,    &
                         nm  = 10/0.52917706d0,   &
                         pi  = 3.1415926535897931d0

   !field component labels for printing
   character(2) :: fc(3) = (/'Fx','Fy','Fz'/)
   character(2) :: bc(3) = (/'Bx','By','Bz'/)

contains

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


   subroutine verdet_dpol_db(mol, S, D, F, freq, Eoff, Ebff)
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(3)
      complex(8),   intent(out) :: Eoff(3,3,3), Ebff(3,3,3)
      type(matrix) :: Db(3), De(3), Df(3), Def(3,3), FbDS(3), &
                      Fb(3), Fe(3), Ff(3), Fef(3,3), DbSD(3), &
                      Sb(3), DFDef(3,3), FDSef, DSDef
      integer      :: i, j, k
      integer, parameter :: eps(3,3,3) =       &
          reshape((/0, 0, 0, 0, 0,-1, 0,+1, 0, &
                    0, 0,+1, 0, 0, 0,-1, 0, 0, &      
                    0,-1, 0,+1, 0, 0, 0, 0, 0/), (/3,3,3/))
      if (abs(sum(freq)) > 1d-15) &
         call quit('verdet_dpol_db: sum(freq) should be zero!')
      !----------------
      ! solve equations
      call pert_dens(mol, S, (/'MAG'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq(1)/))
      call pert_dens(mol, S, (/'EL' /), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(2)/))
      call pert_dens(mol, S, (/'EL' /), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(3)/))
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(2),freq(3)/))
      !----------------------
      ! no-London Ebff tensor
      Eoff = 0
      call prop_oneave(mol, S, (/'MAGO'/), (/Def/), shape(Eoff), Eoff)
      call print_tensor(shape(Eoff), -Eoff, 'no-Lon d/db alpha = -Eoff', freq)
      print *, 'sum_ijk -eijk*Eoff_ijk=', &
               sum((/(((-eps(i,j,k)*Eoff(i,j,k),i=1,3),j=1,3),k=1,3)/))
      print *
      !-------------------------------------------------------------------
      ! London Ebff, 0+2+1 rule, only energy (oneave/twoave) contributions
      Ebff = 0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Df/), shape(Ebff), Ebff)
      ! call print_tensor(shape(Ebff), Ebff, 'E1ebDf'); Ebff=0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/De/), shape(Ebff), Ebff, perm=(/1,3,2/))
      ! call print_tensor(shape(Ebff), Ebff, 'E1fbDe'); Ebff=0
      call prop_twoave(mol, (/'MAG'/), (/D,De,Df,Def/), shape(Ebff), Ebff)
      do k = 1, 3
         do j = 1, 3
            DFDef(j,k) = D*Fe(j)*Df(k) + De(j)*Ff(k)*D &
                       + Df(k)*Fe(j)*D + D*Ff(k)*De(j) &
                       + De(j)*(F+(freq(2)-freq(3))/2*S)*Df(k) &
                       + Df(k)*(F-(freq(2)-freq(3))/2*S)*De(j) &
                       + Def(j,k)*(F+(freq(2)+freq(3))/2*S)*D + D*Fef(j,k)*D &
                       + D*(F-(freq(2)+freq(3))/2*S)*Def(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'MAG'/), (/Def/), shape(Ebff), Ebff, &
                       freq=(/freq(1)/), DFD=(/DFDef/))
      DFDef = 0
      ! call print_tensor(shape(Ebff), Ebff, 'E2bDeDf + E1bDef - i/2TbDef - Sb(DFD)ef'); Ebff=0
      call print_tensor(shape(Ebff), -Ebff, 'London d/db alpha = -Ebff', freq)
      print *, 'sum_ijk -eijk*Ebff_ijk=', &
               sum((/(((-eps(i,j,k)*Ebff(i,j,k),i=1,3),j=1,3),k=1,3)/))
      print *
      !-----------------------------------------------------
      ! London Effb, 1+1+1 rule (2n+1 rule), energy contribs
      Ebff = 0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Df/), shape(Ebff), Ebff)
      ! call print_tensor(shape(Ebff), Ebff, 'E1ebDf'); Ebff=0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/De/), shape(Ebff), Ebff, perm=(/1,3,2/))
      ! call print_tensor(shape(Ebff), Ebff, 'E1fbDe'); Ebff=0
      call prop_twoave(mol, (/'MAG'/), (/D,De,Df,(0d0*D,i=1,3*3)/), shape(Ebff), Ebff)
      ! call print_tensor(shape(Ebff), Ebff, 'E2bDeDf'); Ebff=0
      ! partial reorth contribution
      call prop_oneint(mol, S, (/'MAG'/), (/3/), S=Sb) !load overlap
      do k = 1, 3
         do j = 1, 3
            DFDef(j,k) = D*Fe(j)*Df(k) + De(j)*Ff(k)*D   &
                       + Df(k)*Fe(j)*D + D*Ff(k)*De(j)   &
                       + De(j)*(F+(freq(2)-freq(3))/2*S)*Df(k) &
                       + Df(k)*(F-(freq(2)-freq(3))/2*S)*De(j)
            do i = 1, 3
               Ebff(i,j,k) = Ebff(i,j,k) - tr(Sb(i),DFDef(j,k))
            end do; DFDef(j,k)=0
         end do
      end do
      ! call print_tensor(shape(Ebff), Ebff, '-Sb(DFD)ef'); Ebff=0
      ! gradient multiplier contribution -tr (DbSD-) (FDS-)ef
      do i = 1, 3
         DbSD(i) = Db(i)*S*D - D*S*Db(i)
      end do
      do k = 1, 3
         do j = 1, 3
            FDSef = (Fe(j)*Df(k) + Ff(k)*De(j))*S - S*(De(j)*Ff(k) + Df(k)*Fe(j))
            do i = 1, 3
               Ebff(i,j,k) = Ebff(i,j,k) - tr(DbSD(i),FDSef)
            end do; FDSef=0
         end do
      end do
      DbSD = 0
      ! call print_tensor(shape(Ebff), Ebff, '-(DbSD)(FDS)fe'); Ebff=0
      ! idempotency multiplier contribution -tr (FbDS+) (DSD)fe
      do i = 1, 3
         FbDS(i) = Fb(i)*D*S + S*D*Fb(i) - Fb(i) &
                 - F*D*Sb(i) - Sb(i)*D*F
         Sb(i) = 0
      end do
      do k = 1, 3
         do j = 1, 3
            DSDef = De(j)*S*Df(k) + Df(k)*S*De(j)
            do i = 1, 3
               Ebff(i,j,k) = Ebff(i,j,k) - tr(FbDS(i),DSDef)
            end do; DSDef=0
         end do
      end do
      FbDS = 0
      ! call print_tensor(shape(Ebff), Ebff, '-(FbDS)(DSD)fe'); Ebff=0
      call print_tensor(shape(Ebff), -Ebff, 'London d/db alpha = -Ebff', freq)
      print *, 'sum_ijk -eijk*Ebff_ijk=', &
               sum((/(((-eps(i,j,k)*Ebff(i,j,k),i=1,3),j=1,3),k=1,3)/))
      print *
      Db=0; Fb=0; De=0; Fe=0; Df=0; Ff=0
   end subroutine


   subroutine efgb_Jpri_Bten_Bcal(mol, S, D, F, freq, Ef, Eq, Eff, Evf, Eof, Ebf, &
                                  Eqf, Efff, Evff, Eoff, Ebff, Eqff, Effq)
   ! Electric-field-gradient-induced birefringence (Buckingham birefringence).
   ! For freq=(w1,w2,w3) calculate the dipole moment in Ef, quadrupole moment in Eq,
   ! polarizability alpha(w1,w2) in Eff, velocity-polarizability(w1,w2) in Evf,
   ! no-London G-prime(w1,w2) in Eof, London G-prime(w1,w2) in Ebf, A-tensor(w1,w2)
   ! in Eqf, hyperpolarizability beta(w1,w2,w3) in Efff, no-London 
   ! electric-field-perturbed G-prime in Eoff, London electric-field-perturbed
   ! G-prime in Ebff, electric-field-perturbed A-tensor in Eqff, and
   ! electric-field-gradient-perturbed polarizability in Effq
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(3)
      complex(8),   intent(out) :: Ef(3), Eq(6), Eff(3,3), Evf(3,3), Eof(3,3), &
                                   Ebf(3,3), Eqf(6,3), Efff(3,3,3), Evff(3,3,3), &
                                   Eoff(3,3,3), Ebff(3,3,3), Eqff(6,3,3), Effq(3,3,6)
      type(matrix) :: De(3), Df(3), Def(3,3), Dte(3,3), Dt(3), DFDef(3,3)
      type(matrix) :: Fe(3), Ff(3), Fef(3,3), Fte(3,3), Ft(3), DFDe(3)
      integer      :: i, j ,k
      if (abs(sum(freq)) > 1d-15) &
         call quit('efgb_Jpri_Bten_Bcal_new: sum(freq) should be zero!')
      ! dipole and quadrupole moments
      Ef=0; Eq=0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), Ef)
      call print_tensor((/3/), -Ef, 'dipmom = -Ef')
      call prop_oneave(mol, S, (/'ELGR'/), (/D/), (/6/), Eq)
      call print_tensor((/6/), -Eq, 'quadru = -Eq')
      ! first order equations for last field, and polarizability
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(3)/))
      Eff = 0 !PS: alpha(-freq(3),freq(3)) is not returned from this subroutine
      call prop_oneave(mol, S, (/'EL'/), (/Df/), shape(Eff), Eff)
      call print_tensor(shape(Eff), -Eff, 'alpha = -Eff', (/-freq(3),freq(3)/))
      ! first order equations for middle field, and various linear responses
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(2)/))
      Eff=0; Evf=0; Eof=0; Ebf=0; Eqf=0
      call prop_oneave(mol, S, (/'EL'/), (/De/), shape(Eff), Eff)
      call print_tensor(shape(Eff), -Eff, 'alpha = -Eff', (/-freq(2),freq(2)/))
      call prop_oneave(mol, S, (/'VEL'/), (/De/), shape(Evf), Evf)
      call print_tensor(shape(Evf), -Evf, 'alpha-vel = -Evf', (/-freq(2),freq(2)/))
      call prop_oneave(mol, S, (/'MAGO'/), (/De/), shape(Eof), Eof)
      call print_tensor(shape(Eof),  Eof, 'G-prime-noL = Eof', (/-freq(2),freq(2)/))
      call prop_oneave(mol, S, (/'ELGR'/), (/De/), shape(Eqf), Eqf)
      call print_tensor(shape(Eqf), -Eqf, 'A-tensor = -Eqf', (/-freq(2),freq(2)/))
      call prop_twoave(mol, (/'MAG'/), (/D,De/), shape(Ebf), Ebf)
      do j = 1, 3
         DFDe(j) = D*(F-freq(2)/2*S)*De(j) + D*Fe(j)*D &
                 + De(j)*(F+freq(2)/2*S)*D
      end do
      call prop_oneave(mol, S, (/'MAG'/), (/De/), shape(Ebf), Ebf, &
                       freq=(/-freq(2)/), DFD=(/DFDe/))
      DFDe = 0
      call print_tensor(shape(Ebf), Ebf, 'G-prime-Lon = Ebf', (/-freq(2),freq(2)/))
      ! second order equations for last two fields
      call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                     Def, Fef, freq=(/freq(2),freq(3)/))
      ! beta, beta-mixed, no-Lon J-prime, B-tensor, London J-prime
      Efff=0; Evff=0; Eoff=0; Ebff=0; Eqff=0
      call prop_oneave(mol, S, (/'EL'/), (/Def/), shape(Efff), Efff)
      call print_tensor(shape(Efff), -Efff, 'd/df alpha = -Efff', freq)
      call prop_oneave(mol, S, (/'VEL'/), (/Def/), shape(Evff), Evff)
      call print_tensor(shape(Evff), -Evff, 'd/df alpha-vel = -Evff', freq)
      call prop_oneave(mol, S, (/'MAGO'/), (/Def/), shape(Eoff), Eoff)
      call print_tensor(shape(Eoff), Eoff, 'd/df G-prime-noL = Eoff', freq)
      call prop_oneave(mol, S, (/'ELGR'/), (/Def/), shape(Eqff), Eqff)
      call print_tensor(shape(Eqff), -Eqff, 'd/df A-tensor = -Eqff', freq)
      call prop_oneave(mol, S, (/'MAG','EL '/), (/De/), shape(Ebff), Ebff, perm=(/1,3,2/))
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Df/), shape(Ebff), Ebff)
      call prop_twoave(mol, (/'MAG'/), (/D,De,Df,Def/), shape(Ebff), Ebff)
      do k = 1, 3
         do j = 1, 3
            DFDef(j,k) = D*Fe(j)*Df(k) + De(j)*Ff(k)*D &
                       + Df(k)*Fe(j)*D + D*Ff(k)*De(j) &
                       + De(j)*(F+(freq(2)-freq(3))/2*S)*Df(k) &
                       + Df(k)*(F-(freq(2)-freq(3))/2*S)*De(j) &
                       + D*(F-(freq(2)+freq(3))/2*S)*Def(j,k)  &
                       + Def(j,k)*(F+(freq(2)+freq(3))/2*S)*D + D*Fef(j,k)*D
         end do
      end do
      call prop_oneave(mol, S, (/'MAG'/), (/Def/), shape(Ebff), Ebff, &
                       freq=(/freq(1)/), DFD=(/DFDef/))
      DFDef = 0
      ! finished with f, free matrices
      Df=0; Ff=0; Def=0; Fef=0
      call print_tensor(shape(Ebff), Ebff, 'd/df G-prime-Lon = Ebff', freq)
      ! second order equations for first two fields.
      ! If possible, solve 6 instead of 3x3 eqs.
      if (freq(1)==freq(2) .or. freq(1)==-freq(2)) then
         do j = 1, 3
            if (freq(1)==freq(2)) Dt(j) = De(j)
            if (freq(1)==freq(2)) Ft(j) = Fe(j)
            if (freq(1)/=freq(2)) Dt(j) = dag(De(j))
            if (freq(1)/=freq(2)) Ft(j) = dag(Fe(j))
            call pert_dens(mol, S, (/'EL','EL'/), (/j,1/), (/D,Dt(:j),De(j)/), &
                           (/F,Ft(:j),Fe(j)/), Dte(:j,j), Fte(:j,j), &
                           freq=(/freq(1),freq(2)/), comp=(/1,j/))
            do i = 1, j-1
               if (freq(1)==freq(2)) Dte(j,i) = Dte(i,j)
               if (freq(1)==freq(2)) Fte(j,i) = Fte(i,j)
               if (freq(1)/=freq(2)) Dte(j,i) = dag(Dte(i,j))
               if (freq(1)/=freq(2)) Fte(j,i) = dag(Fte(i,j))
            end do
         end do
      else
         call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Dt, Ft, freq=(/freq(1)/))
         call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,Dt,De/), (/F,Ft,Fe/), &
                        Dte, Fte, freq=(/freq(1),freq(2)/))
      end if
      ! field-gradient-perturbed polarizability B-calligraphic
      Effq = 0
      call prop_oneave(mol, S, (/'ELGR'/), (/Dte/), shape(Effq), Effq, perm=(/3,1,2/))
      call print_tensor(shape(Effq), -Effq, 'd/dq Alpha = -Effq', freq)
      ! free matrices, change sign on quasi-energy derivatives to get
      ! property-derivatives
      Dt=0; Ft=0; De=0; Fe=0; Dte=0; Fte=0
      Ef=-Ef; Eq=-Eq; Eff=-Eff; Evf=-Evf; Eqf=-Eqf
      Efff=-Efff; Evff=-Evff; Eqff=-Eqff; Effq=-Effq
   end subroutine


   subroutine efgb_output(w, dip, qua, Pff, Pfv, Pfo, Pfb, Pfq, &
                          Pfff, Pffv, Pffo, Pffb, Pffq, Pffs, o)
   !ajt jan09 Based on Dmitry Shcherbin's post-calculation script efgb.py
      real(8),intent(in)    :: w
      real(8),intent(in)    :: dip(3),qua(6),Pff(3,3),Pfv(3,3),Pfo(3,3), &
                               Pfb(3,3),Pfq(3,6),Pfff(3,3,3),Pffv(3,3,3), &
                               Pffo(3,3,3),Pffb(3,3,3),Pffq(3,3,6),Pffs(3,3,6)
      integer,intent(in)    :: o !outputting unit (lupri or stdout)
      real(8) :: tr3,trm3(3,3),Jprio,Jpril,Bten,Bcal
      real(8),parameter:: c1=3.75874d-21,c2=63154.9d0,c3=6.18385d-14
      real(8),parameter:: hz=1/48168836676d0,c4=0.164518d0*c1/cm1
      !calculate tensor averages sum J'_ijk eps_ijk
      Jprio = Pffo(1,2,3) + Pffo(3,1,2) + Pffo(2,3,1) & !no-London
            - Pffo(3,2,1) - Pffo(2,1,3) - Pffo(1,3,2)
      Jpril = Pffb(1,2,3) + Pffb(3,1,2) + Pffb(2,3,1) & !London
            - Pffb(3,2,1) - Pffb(2,1,3) - Pffb(1,3,2)
      !average sum B_ij(ij)
      Bten = Pffs(1,1,1) + Pffs(2,1,2) + Pffs(3,1,4) &
           + Pffs(1,2,2) + Pffs(2,2,3) + Pffs(3,2,5) &
           + Pffs(1,3,4) + Pffs(2,3,5) + Pffs(3,3,6)
      Bcal = Pffq(1,1,1) + Pffq(2,1,2) + Pffq(3,1,4) &
           + Pffq(1,2,2) + Pffq(2,2,3) + Pffq(3,2,5) &
           + Pffq(1,3,4) + Pffq(2,3,5) + Pffq(3,3,6)
      !print
      write (o,'(//11x,a/11x,a/)')                                                  &
               'Electric-field-gradient-induced (Buckingham) birefringence (EFGB)', &
               '-----------------------------------------------------------------'
      write (o,'(6x,a,f12.9,a,i7,a,i4,a/)') 'For frequency w = ',w, &
               ' au =',nint(w/cm1),' cm-1 ~ ', nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(4(6x,a,1x,g13.7/)/2(6x,a,1x,g13.7/)/)') &
               'B-tensor           <<mu;mu,th>>w0 = -Effq(-w,w,0) average     :',Bten,  &
               'B-calligraphic     <<mu;mu,th>>0w = -Effq(-w,0,w) average     :',Bcal,  &
               'no-London J-prime i<<mu;mu,mg>>0w = -Effb(-w,0,w) eps-average :',Jprio, &
               '   London J-prime i<<mu;mu,mg>>0w = -Effb(-w,0,w) eps-average :',Jpril, &
               'no-London Buckingham birefringence 2/15(Bten-Bcal)-2/3w*Jpri  :',       &
                                          2/15d0 * (Bten-Bcal) - 2/(3*w) * Jprio,       &
               '   London Buckingham birefringence 2/15(Bten-Bcal)-2/3w*Jpri  :',       &
                                          2/15d0 * (Bten-Bcal) - 2/(3*w) * Jpril
   end subroutine


   subroutine cme_jones_eta_apri(mol, S, D, F, freq, Ef, Eq, Eoo, Ebb, Eof, Ebf, Eqf, Eff, &
                                 Eoof, Ebbf, Eqof, Eqbf, Eooff, Ebbff, Eqoff, Eqbff)
   ! ajt jan09 For magneto-electric Jones- and magneto-magnetic Cotton-Mouton birefringences. 
   ! For Jones, the four frequencies should be (w,0,-w,0), for CME, (0,0,-w,w)
   ! In the code, 'e' denotes the first electric field, 'f' the second electric field,
   ! 'b' the first magnetic field, 'c' the second magnetic field, and 'q' the electric field gradient
   ! gradient. The response eqs solved are: All first order: (b,c,e,f,q: 3+3+3+3+6), and
   ! second order NOT involving the first field (b or q) (ce,cf,ef, 9+9+9), thus 45 in total.
   ! ajt sep09 Reversed field order
      type(prop_molcfg), intent(in) :: mol
      type(matrix), intent(in)  :: S, D, F
      complex(8),   intent(in)  :: freq(4)
      complex(8),   intent(out) :: Ef(3), Eq(6), Eoo(3,3), Ebb(3,3), Eof(3,3), Ebf(3,3), &
                                   Eqf(6,3), Eff(3,3), Eoof(3,3,3), Ebbf(3,3,3), &
                                   Eqof(6,3,3), Eqbf(6,3,3), Eooff(3,3,3,3), &
                                   Ebbff(3,3,3,3), Eqoff(6,3,3,3), Eqbff(6,3,3,3)
      type(matrix) :: Db(3), Dc(3), De(3), Df(3), Dq(6), Dce(3,3), Dcf(3,3), Def(3,3), &
                      Fb(3), Fc(3), Fe(3), Ff(3), Fq(6), Fce(3,3), Fcf(3,3), Fef(3,3), &
                      DFD, Sb(3), DFDc(3), DFDcef, DbSD(3), DqSD(6), FDScef, &
                      FbDS(3), FqDS(6), DSDcef, DFDe(3), DFDce(3,3), DFDef(3,3)
      integer      :: i, j, k, l
      complex(8)   :: tmp(3,3,3,3)
      if (abs(sum(freq)) > 1d-15) &
         call quit('cme_jones_eta_apri: sum(freq) should be zero!')
      ! dipole and quadrupole moments
      Ef = 0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), Ef)
      call print_tensor((/3/), -Ef, 'dipmom = -Ef')
      Eq = 0
      call prop_oneave(mol, S, (/'ELGR'/), (/D/), (/6/), Eq)
      call print_tensor((/6/), -Eq, 'quadru = -Eq')
      ! solve equations common to no-London and London
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq(4)/))
      Eff = 0
      if (freq(3) == -freq(4)) then
         call prop_oneave(mol, S, (/'EL'/), (/Df/), shape(Eff), Eff)
         call print_tensor(shape(Eff), -Eff, 'polari = -Eff', (/freq(3),freq(4)/))
      end if
      Eqf = 0
      if (freq(1) == -freq(4)) then
         call prop_oneave(mol, S, (/'ELGR'/), (/Df/), shape(Eqf), Eqf)
         call print_tensor(shape(Eqf), -Eqf, 'polari = -Eff', (/freq(1),freq(4)/))
      end if
      call pert_dens(mol, S, (/'ELGR'/), (/6/), (/D/), (/F/), Dq, Fq, freq=(/freq(1)/))
      if (freq(3)==freq(4) .or. freq(3)==-freq(4)) then
         do l = 1, 3
            if (freq(3)==freq(4)) De(l) = Df(l)
            if (freq(3)==freq(4)) Fe(l) = Ff(l)
            if (freq(3)/=freq(4)) De(l) = dag(Df(l))
            if (freq(3)/=freq(4)) Fe(l) = dag(Ff(l))
            call pert_dens(mol, S, (/'EL','EL'/), (/l,1/), (/D,De(:l),Df(l)/), &
                           (/F,Fe(:l),Ff(l)/), Def(:l,l), Fef(:l,l), &
                           freq=(/freq(3),freq(4)/), comp=(/1,l/))
            do k = 1, l-1
               if (freq(3)==freq(4)) Def(l,k) = Def(k,l)
               if (freq(3)==freq(4)) Fef(l,k) = Fef(k,l)
               if (freq(3)/=freq(4)) Def(l,k) = dag(Def(k,l))
               if (freq(3)/=freq(4)) Fef(l,k) = dag(Fef(k,l))
            end do
         end do
      else
         call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), De, Fe, freq=(/freq(3)/))
         call pert_dens(mol, S, (/'EL','EL'/), (/3,3/), (/D,De,Df/), (/F,Fe,Ff/), &
                        Def, Fef, freq=(/freq(3),freq(4)/))
      end if
      !-------------------!
      ! no-London section !
      !-------------------!
      ! solve equations specific to no-London
      call pert_dens(mol, S, (/'MAGO'/), (/3/), (/D/), (/F/), Dc, Fc, freq=(/freq(2)/))
      Eoo = 0
      if (freq(1)==-freq(2)) then
         call prop_oneave(mol, S, (/'MAGO','MAGO'/), (/D/), shape(Eoo), Eoo)
         call prop_oneave(mol, S, (/'MAGO'/), (/Dc/), shape(Eoo), Eoo)
         call print_tensor(shape(Eoo), Eoo, 'magneti-noL = Eoo', (/freq(1),freq(2)/))
      end if
      if (freq(1)==freq(2) .or. freq(1)==-freq(2)) then
         do j = 1, 3
            if (freq(1)==freq(2)) Db(j) = Dc(j)
            if (freq(1)==freq(2)) Fb(j) = Fc(j)
            if (freq(1)/=freq(2)) Db(j) = -dag(Dc(j))
            if (freq(1)/=freq(2)) Fb(j) = -dag(Fc(j))
         end do
      else
         call pert_dens(mol, S, (/'MAGO'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq(1)/))
      end if
      Eof = 0
      if (freq(3)==-freq(1)) then
         call prop_oneave(mol, S, (/'EL'/), (/Db/), shape(Eof), Eof, perm=(/2,1/))
         call print_tensor(shape(Eof), Eof, 'G-prime-noL = Eof', (/freq(1),freq(3)/))
      end if
      call pert_dens(mol, S, (/'MAGO','EL  '/), (/3,3/), (/D,Dc,Df/), (/F,Fc,Ff/), &
                     Dcf, Fcf, freq=(/freq(2),freq(4)/))
      if (freq(3)==freq(4) .or. freq(3)==-freq(4)) then
         do l = 1, 3
            do j = 1, 3
               if (freq(3)==freq(4)) Dce(j,l) = Dcf(j,l)
               if (freq(3)==freq(4)) Fce(j,l) = Fcf(j,l)
               if (freq(3)/=freq(4)) Dce(j,l) = -dag(Dcf(j,l))
               if (freq(3)/=freq(4)) Fce(j,l) = -dag(Fcf(j,l))
            end do
         end do
      else
         call pert_dens(mol, S, (/'MAGO','EL  '/), (/3,3/), (/D,Dc,De/), (/F,Fc,Fe/), &
                        Dce, Fce, freq=(/freq(2),freq(3)/))
      end if
      !--------------------------------------------------
      ! contract no-London d/db G-prime and d/db A-tensor
      Eoof = 0; Eqof = 0
      if (freq(4)==0) then
         call prop_oneave(mol, S, (/'MAGO','MAGO'/), (/De/), shape(Eoof), Eoof) !dia
         call prop_oneave(mol, S, (/'MAGO'/), (/Dce/), shape(Eoof), Eoof) !para
         call print_tensor(shape(Eoof), Eoof, 'd/db G-prime-noL = Eoof', freq(:3))
         call prop_oneave(mol, S, (/'ELGR'/), (/Dce/), shape(Eqof), Eqof)
         call print_tensor(shape(Eqof), -Eqof, 'd/db A-tensor-noL = -Eqof', freq(:3))
      end if
      !-------------------------------------------------------
      ! simultaneously contract no-London hypermagnetizability
      ! d/df d/db G-prime and d/df d/db A-tensor
      ! just one energy contribution for G-prime, none for A-tensor
      Eooff = 0
      Eqoff = 0
      call prop_oneave(mol, S, (/'MAGO','MAGO'/), (/Def/), shape(Eooff), Eooff)
      ! gradient multiplier contribution -tr (DbSD-) (FDS)cef (or Dq instead of Db)
      do i = 1,3
         DbSD(i) = Db(i)*S*D - D*S*Db(i)
      end do
      do i = 1, 6
         DqSD(i) = Dq(i)*S*D - D*S*Dq(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               FDScef = (Fc(j)*Def(k,l) + Fe(k)*Dcf(j,l) + Ff(l)*Dce(j,k) &
                      +  Fef(k,l)*Dc(j) + Fcf(j,l)*De(k) + Fce(j,k)*Df(l)) * S &
                  - S * (Dc(j)*Fef(k,l) + De(k)*Fcf(j,l) + Df(l)*Fce(j,k) &
                      +  Def(k,l)*Fc(j) + Dcf(j,l)*Fe(k) + Dce(j,k)*Ff(l))
               do i = 1, 3
                  Eooff(i,j,k,l) = Eooff(i,j,k,l) - tr(DbSD(i),FDScef)
               end do
               do i = 1, 6
                  Eqoff(i,j,k,l) = Eqoff(i,j,k,l) - tr(DqSD(i),FDScef)
               end do; FDScef=0
            end do
         end do
      end do; DbSD=0; DqSD=0
      ! idempotency multiplier contribution -tr (FbDS+) (DSD)cef
      do i = 1, 3
         FbDS(i) = Fb(i)*D*S + S*D*Fb(i) - Fb(i)
      end do
      do i = 1, 6
         FqDS(i) = Fq(i)*D*S + S*D*Fq(i) - Fq(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DSDcef = Dc(j)*S*Def(k,l) + De(k)*S*Dcf(j,l) + Df(l)*S*Dce(j,k) &
                      + Def(k,l)*S*Dc(j) + Dcf(j,l)*S*De(k) + Dce(j,k)*S*Df(l)
               do i = 1, 3
                  Eooff(i,j,k,l) = Eooff(i,j,k,l) - tr(FbDS(i),DSDcef)
               end do
               do i = 1, 6
                  Eqoff(i,j,k,l) = Eqoff(i,j,k,l) - tr(FqDS(i),DSDcef)
               end do; DSDcef=0
            end do
         end do
      end do; FbDS=0; FqDS=0
      ! print no-London response tensors to stdout
      call print_tensor(shape(Eooff), Eooff, 'd/df d/db G-prime-noL = Eooff', freq)
      call print_tensor(shape(Eqoff),-Eqoff, 'd/df d/db A-tensor-noL = -Eqoff', freq)
      !free no-London-specific matrices
      Db=0; Fb=0; Dc=0; Fc=0; Dce=0; Fce=0; Dcf=0; Fcf=0
      !----------------!
      ! London section !
      !----------------!
      ! solve equations specific to London
      call pert_dens(mol, S, (/'MAG'/), (/3/), (/D/), (/F/), Dc, Fc, freq=(/freq(2)/))
      Ebb = 0
      ! contract London magnetizability
      if (freq(1)==-freq(2)) then
         call prop_twoave(mol, (/'MAG','MAG'/), (/D/), shape(Ebb), Ebb)
         DFD = D*F*D
         call prop_oneave(mol, S, (/'MAG','MAG'/), (/D/), shape(Ebb), Ebb, &
                          freq=(/freq(1),freq(2)/), DFD=(/DFD/))
         DFD = 0
         do j = 1, 3
            DFDc(j) = Dc(j)*(F+freq(2)/2*S)*D + D*Fc(j)*D &
                        + D*(F-freq(2)/2*S)*Dc(j)
         end do
         call prop_twoave(mol, (/'MAG'/), (/D,Dc/), shape(Ebb), Ebb)
         call prop_oneave(mol, S, (/'MAG'/), (/Dc/), shape(Ebb), Ebb, &
                          freq=(/freq(1)/), DFD=(/DFDc/))
         DFDc = 0
         call print_tensor(shape(Ebb), Ebb, 'magneti-Lon = Ebb', (/freq(1),freq(2)/))
      end if
      if (freq(1)==freq(2) .or. freq(1)==-freq(2)) then
         do j = 1, 3
            if (freq(1)==freq(2)) Db(j) = Dc(j)
            if (freq(1)==freq(2)) Fb(j) = Fc(j)
            if (freq(1)/=freq(2)) Db(j) = -dag(Dc(j))
            if (freq(1)/=freq(2)) Fb(j) = -dag(Fc(j))
         end do
      else
         call pert_dens(mol, S, (/'MAG'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq(1)/))
      end if
      Ebf = 0
      if (freq(3)==-freq(1)) then
         call prop_oneave(mol, S, (/'EL'/), (/Db/), shape(Ebf), Ebf, perm=(/2,1/))
         call print_tensor(shape(Ebf), Ebf, 'G-prime-Lon = Ebf', (/freq(1),freq(3)/))
      end if
      call pert_dens(mol, S, (/'MAG','EL '/), (/3,3/), (/D,Dc,Df/), (/F,Fc,Ff/), &
                     Dcf, Fcf, freq=(/freq(2),freq(4)/))
      if (freq(3)==freq(4) .or. freq(3)==-freq(4)) then
         do l = 1, 3
            do j = 1, 3
               if (freq(3)==freq(4)) Dce(j,l) = Dcf(j,l)
               if (freq(3)==freq(4)) Fce(j,l) = Fcf(j,l)
               if (freq(3)/=freq(4)) Dce(j,l) = -dag(Dcf(j,l))
               if (freq(3)/=freq(4)) Fce(j,l) = -dag(Fcf(j,l))
            end do
         end do
      else
         call pert_dens(mol, S, (/'MAG','EL '/), (/3,3/), (/D,Dc,De/), (/F,Fc,Fe/), &
                        Dce, Fce, freq=(/freq(2),freq(3)/))
      end if
      ! magnetic overlap is needed for DFDce, and later for DFDcef, FDScef and FbDS
      call prop_oneint(mol, S, (/'MAG'/), (/3/), S=Sb)
      !-----------------------------------------------
      ! contract London d/db G-prime and d/db A-tensor
      Ebbf = 0
      Eqbf = 0
      if (freq(4)==0) then
         call prop_oneave(mol, S, (/'MAG','MAG','EL '/), (/D/), shape(Ebbf), Ebbf)
         ! call print_tensor(shape(Ebbf), Ebbf, 'E0bce'); Ebbf=0
         call prop_twoave(mol, (/'MAG','MAG'/), (/D,De/), shape(Ebbf), Ebbf)
         do k = 1, 3
            DFDe(k) = De(k)*(F+freq(3)/2*S)*D + D*Fe(k)*D &
                        + D*(F-freq(3)/2*S)*De(k)
         end do
         call prop_oneave(mol, S, (/'MAG','MAG'/), (/De/), shape(Ebbf), Ebbf, &
                          freq=(/freq(1),freq(2)/), DFD=(/DFDe/))
         DFDe = 0
         ! call print_tensor(shape(Ebbf), Ebbf, 'E1bcDe-i/2TbcDe-SbcWe'); Ebbf=0
         call prop_oneave(mol, S, (/'MAG','EL '/), (/Dc/), shape(Ebbf), Ebbf, perm=(/1,3,2/))
         ! call print_tensor(shape(Ebbf), Ebbf, 'E1ecDb'); Ebbf=0
         call prop_twoave(mol, (/'MAG'/), (/D,Dc,De,Dce/), shape(Ebbf), Ebbf)
         do k = 1, 3
            do j = 1, 3
               DFDce(j,k) = Dc(j)*Fe(k)*D + D*(Fc(j)-freq(3)/2*Sb(j))*De(k) &
                          + D*Fe(k)*Dc(j) + De(k)*(Fc(j)+freq(3)/2*Sb(j))*D &
                          + Dc(j)*(F+(freq(2)-freq(3))/2*S)*De(k) &
                          + De(k)*(F-(freq(2)-freq(3))/2*S)*Dc(j) &
                          + D*(F-(freq(2)+freq(3))/2*S)*Dce(j,k) &
                          + Dce(j,k)*(F+(freq(2)+freq(3))/2*S)*D + D*Fce(j,k)*D
            end do
         end do
         call prop_oneave(mol, S, (/'MAG'/), (/Dce/), shape(Ebbf), Ebbf, &
                          freq=(/freq(1)/), DFD=(/DFDce/))
         DFDce = 0
         ! call print_tensor(shape(Ebbf), Ebbf, 'E1cDeb+E2cDeDb-i/2TcDeb-ScWeb'); Ebbf=0
         call print_tensor(shape(Ebbf), Ebbf, 'd/db G-prime-Lon = Ebbf', freq(:3))
         call prop_oneave(mol, S, (/'ELGR','MAG '/), (/De/), shape(Eqbf), Eqbf)
         call prop_oneave(mol, S, (/'ELGR'/), (/Dce/), shape(Eqbf), Eqbf)
         call print_tensor(shape(Eqbf), -Eqbf, 'd/db A-tensor-Lon = -Efbq', freq(:3))
      end if
      !----------------------------------------------------
      ! simultaneously contract London hypermagnetizability
      ! d/df d/db G-prime and d/df d/db A-tensor
      ! 7 energy contributions to hypermag, 2 to A-tensor
      Eqbff = 0
      Ebbff = 0
      call prop_oneave(mol, S, (/'ELGR','MAG '/), (/Def/), shape(Eqbff), Eqbff)
      ! call print_tensor(shape(Eqbff), Eqbff, 'E1qbDef', Eqbff); Eqbff=0
      !ajt This should really be a 3rd order perturbed density, which isn't
      !    implemented in prop_contribs. Rather, fake it with a second order density
      !    NB: This is incorrect with DFT
      call prop_twoave(mol, (/'MAG'/), (/D,Dq,Def,(0d0*D,i=1,6*3*3)/), &
                       (/6,3,3*3/), Eqbff, perm=(/2,1,3/))
      ! call print_tensor(shape(Eqbff), Eqbff, 'E2bDqDef', Eqbff); Eqbff=0
      call prop_oneave(mol, S, (/'MAG','MAG','EL '/), (/Df/), shape(Ebbff), Ebbff)
      call prop_oneave(mol, S, (/'MAG','MAG','EL '/), (/De/), shape(Ebbff), Ebbff, &
                       perm=(/1,2,4,3/))
      ! call print_tensor(shape(Ebbff), Ebbff, 'E1bcfDe + E1bceDf'); Ebbff=0
      call prop_twoave(mol, (/'MAG','MAG'/), (/D,De,Df,Def/), shape(Ebbff), Ebbff)
      do l = 1, 3
         do k = 1, 3
            DFDef(k,l) = De(k)*Ff(l)*D + D*Fe(k)*Df(l) &
                       + Df(l)*Fe(k)*D + D*Ff(l)*De(k) &
                       + De(k)*(F+(freq(3)-freq(4))/2*S)*Df(l) &
                       + Df(l)*(F-(freq(3)-freq(4))/2*S)*De(k) &
                       + D*(F-(freq(3)+freq(4))/2*S)*Def(k,l)  &
                       + Def(k,l)*(F+(freq(3)+freq(4))/2*S)*D + D*Fef(k,l)*D
         end do
      end do
      call prop_oneave(mol, S, (/'MAG','MAG'/), (/Def/), shape(Ebbff), Ebbff, &
                       freq=(/freq(1),freq(2)/), DFD=(/DFDef/))
      DFDef = 0
      ! call print_tensor(shape(Ebbff), Ebbff, 'E1bcDef + E2bcDeDf - 1/2TbcDef - SbcdDFDef'); Ebbff=0
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Dcf/), shape(Ebbff), Ebbff, perm=(/1,3,2,4/))
      call prop_oneave(mol, S, (/'MAG','EL '/), (/Dce/), shape(Ebbff), Ebbff, perm=(/1,4,2,3/))
      ! call print_tensor(shape(Ebbff), Ebbff, 'E1beDcf + E1bfDce'); Ebbff=0
      !ajt This should really be the 3rd order perturbed density
      !    (/D,Db,De,Df,(0d0*D,i=1,2*3*3),Def,(0d0*D,i=1,3*3*3)/), but 3rd order
      !    isn't implemented in prop_contribs. Rather, fake it with a second order
      !    density. NB: Incorrect with DFT, as E^3b Db De Df will be missing
      call prop_twoave(mol, (/'MAG'/), (/D,Db,Def,(0d0*D,i=1,3*3*3)/), &
                       (/3,3,3*3/), Ebbff, perm=(/2,1,3/))
      !ajt This should really be the 3rd order perturbed density
      !    (/D,Dc,De,Df,Dce,Dcf,Def,(0d0*D,i=1,3*3*3)/), but 3rd order isn't
      !    implemented in prop_contribs. Rather, fake it with three second order
      !    densities. NB: Incorrect with DFT, as E^3b Dc De Df will be missing
      call prop_twoave(mol, (/'MAG'/), (/D,Dc,Def,(0d0*D,i=1,3*3*3)/), &
                       (/3,3,3*3/), Ebbff)
      tmp = 0 !need tmp since index e is between c and f
      call prop_twoave(mol, (/'MAG'/), (/D,De,Dcf,(0d0*D,i=1,3*3*3)/), &
                       (/3,3,3*3/), tmp)
      Ebbff = Ebbff + reshape((/((((tmp(i,k,j,l),i=1,3),j=1,3),k=1,3),l=1,3)/),(/3,3,3,3/))
      call prop_twoave(mol, (/'MAG'/), (/D,Df,Dce,(0d0*D,i=1,3*3*3)/), &
                       (/3,3*3,3/), Ebbff, perm=(/1,3,2/))
      ! call print_tensor(shape(Ebbff), Ebbff, 'E2cDbDef + E2bDcDef + E2bDeDcf + E2bDfDce'); Ebbff=0
      ! reorthonormalization contribution -tr Sb (DFD)cef
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DFDcef = Dc(j)*Fe(k)*Df(l) + Dc(j)*Ff(l)*De(k) &
                      + Df(l)*Fe(k)*Dc(j) + De(k)*Ff(l)*Dc(j) &
                      + De(k)*(Fc(j)-(-freq(3)+freq(4))/2*Sb(j))*Df(l) &
                      + Df(l)*(Fc(j)-(+freq(3)-freq(4))/2*Sb(j))*De(k) &
                      + Dc(j)*(F-(-freq(2)+freq(3)+freq(4))/2*S)*Def(k,l) &
                      + Def(k,l)*(F-(+freq(2)-freq(3)-freq(4))/2*S)*Dc(j) &
                      + Dc(j)*Fef(k,l)*D + D*(Fc(j)-(+freq(3)+freq(4))/2*Sb(j))*Def(k,l) &
                      + D*Fef(k,l)*Dc(j) + Def(k,l)*(Fc(j)-(-freq(3)-freq(4))/2*Sb(j))*D &
                      + De(k)*(F-(+freq(2)-freq(3)+freq(4))/2*S)*Dcf(j,l) &
                      + Dcf(j,l)*(F-(-freq(2)+freq(3)-freq(4))/2*S)*De(k) &
                      + De(k)*Fcf(j,l)*D + D*Fe(k)*Dcf(j,l) &
                      + D*Fcf(j,l)*De(k) + Dcf(j,l)*Fe(k)*D &
                      + Df(l)*(F-(+freq(2)+freq(3)-freq(4))/2*S)*Dce(j,k) &
                      + Dce(j,k)*(F-(-freq(2)-freq(3)+freq(4))/2*S)*Df(l) &
                      + Df(l)*Fce(j,k)*D + D*Ff(l)*Dce(j,k) &
                      + D*Fce(j,k)*Df(l) + Dce(j,k)*Ff(l)*D
               do i = 1, 3
                  Ebbff(i,j,k,l) = Ebbff(i,j,k,l) - tr(Sb(i),DFDcef)
               end do; DFDcef=0
            end do
         end do
      end do
      ! call print_tensor(shape(Ebbff), Ebbff, '-Sb DFDcef'); Ebbff=0
      ! gradient multiplier contribution -tr (DbSD-)(FDS)cef (and q intead of b)
      do i = 1, 3
         DbSD(i) = Db(i)*S*D - D*S*Db(i)
      end do
      do i = 1, 6
         DqSD(i) = Dq(i)*S*D - D*S*Dq(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               FDScef = ((Fc(j)-freq(2)/2*Sb(j))*Def(k,l) &
                      +  Fe(k)*Dcf(j,l) + Ff(l)*Dce(j,k)  &
                      +  Fef(k,l)*Dc(j) + Fcf(j,l)*De(k) + Fce(j,k)*Df(l)) * S &
                  - S * (Def(k,l)*(Fc(j)+freq(2)/2*Sb(j)) &
                      +  Dcf(j,l)*Fe(k) + Dce(j,k)*Ff(l)  &
                      +  Dc(j)*Fef(k,l) + De(k)*Fcf(j,l) + Df(l)*Fce(j,k))   &
                      +     ((F-(freq(3)+freq(4))*S)*Def(k,l)                &
                          +  Fe(k)*Df(l) + Ff(l)*De(k) + Fef(k,l)*D) * Sb(j) &
                  - Sb(j) * (D*Fef(k,l) + De(k)*Ff(l) + Df(l)*Fe(k)          &
                          +  Def(k,l)*(F+(freq(3)+freq(4))*S))
               do i = 1, 3
                  Ebbff(i,j,k,l) = Ebbff(i,j,k,l) - tr(DbSD(i),FDScef)
               end do
               do i = 1, 6
                  Eqbff(i,j,k,l) = Eqbff(i,j,k,l) - tr(DqSD(i),FDScef)
               end do; FDScef=0
            end do
         end do
      end do; DbSD=0; DqSD=0
      ! call print_tensor(shape(Ebbff), Ebbff, '-DbSD FDScef'); Ebbff=0
      ! call print_tensor(shape(Eqbff), Eqbff, '-DqSD FDScef', Eqbff); Eqbff=0
      ! idempotency multiplier contribution -tr (FbDS+) (DSD)feb
      do i = 1, 3
         FbDS(i) = Fb(i)*D*S + S*D*Fb(i) - Fb(i) &
                 - F*D*Sb(i) - Sb(i)*D*F
      end do
      do i = 1, 6
         FqDS(i) = Fq(i)*D*S + S*D*Fq(i) - Fq(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DSDcef = D*Sb(j)*Def(k,l) + De(k)*Sb(j)*Df(l) &
                      + Def(k,l)*Sb(j)*D + Df(l)*Sb(j)*De(k) &
                      + Dc(j)*S*Def(k,l) + De(k)*S*Dcf(j,l) + Df(l)*S*Dce(j,k) &
                      + Def(k,l)*S*Dc(j) + Dcf(j,l)*S*De(k) + Dce(j,k)*S*Df(l)
               do i = 1, 3
                  Ebbff(i,j,k,l) = Ebbff(i,j,k,l) - tr(FbDS(i),DSDcef)
               end do
               do i = 1, 6
                  Eqbff(i,j,k,l) = Eqbff(i,j,k,l) - tr(FqDS(i),DSDcef)
               end do; DSDcef=0
            end do
         end do
      end do; FbDS=0; FqDS=0
      ! call print_tensor(shape(Ebbff), Ebbff, '-FbDS DSDcef'); Ebbff=0
      ! call print_tensor(shape(Eqbff), Eqbff, '-FqDS DSDcef'); Eqbff=0
      ! print to stdout
      call print_tensor(shape(Ebbff), Ebbff, 'd/df d/db G-prime-Lon = Ebbff', freq)
      call print_tensor(shape(Eqbff),-Eqbff, 'd/df d/db A-tensor-Lon = -Eqbff', freq)
      !free London-specific and common matrices
      Db=0; Dc=0; De=0; Df=0; Dq=0; Dce=0; Dcf=0; Def=0
      Fb=0; Fc=0; Fe=0; Ff=0; Fq=0; Fce=0; Fcf=0; Fef=0; Sb=0
      !invert signs on a-prime to get properties from quasi-energy derivatives
      Ef=-Ef; Eq=-Eq; Eff=-Eff; Eqf=-Eqf !ajt WRONG: ; Eqoff=-Eqoff; Eqbff=-Eqbff
   end subroutine


   subroutine cme_output(w,Al,Xn,Xi,Etn,Eta,o)
      real(8) :: w,Al(3,3),Xn(3,3),Etn(3,3,3,3),Xi(3,3),Eta(3,3,3,3),deta,dalxi
      integer :: o,i,j,k
      character*7 :: pw
      real(8),parameter :: c1=3.75874d-21,c2=63154.9d0,c3=6.18385d-14
      real(8),parameter :: hz=1/48168836676d0,c4=0.164518d0*c1/cm1
      write (o,'(//26x,a/26x,a/)') 'Polarizability (au)' &
                                  ,'-------------------'
      write (o,'(4x,a,f6.3,a,i7,a,i4,a/)') 'For frequency w = ',w &
              ,' au =',nint(w/cm1),' cm-1 ~ ', nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(a,13x,a2,18x,a2,18x,a2/)') ' -w  w',fc
      write (o,'(3(1x,a,3f20.12/))') (fc(i),Al(i,:),i=1,3)
      if (any(Xn/=0)) then
         write (o,'(//16x,a/16x,a/)') 'No-London static Magnetizability (au)' &
                                     ,'-------------------------------------'
         write (o,'(19x,a2,18x,a2,18x,a2/)') bc
         write (o,'(3(1x,a,3f20.12/))') (bc(i),Xn(i,:),i=1,3)
         write (o,'(//18x,a/18x,a/)') 'No-London Hypermagnetizability - Eta (au)' &
                                     ,'-----------------------------------------'
         write (o,'(4x,a,f6.3,a,i7,a,i4,a/)') 'For frequency w = ',w,' au =' &
                 ,nint(w/cm1),' cm-1 ~ ',nint((2*pi*cvl/w)/nm),' nm'
         write (o,'(a,18x,a2,18x,a2,18x,a2)') ' w;-w',bc
         write (o,'(3(3(/3(1x,3a2,3f20.11/))))') &
                 (((fc(i),fc(j),bc(k),Etn(i,j,k,:),k=1,3),j=1,3),i=1,3)
         write (o,'(//16x,a/16x,a/)') &
                  'No-London Cotton-mouton constant: mC, dn, Ccm' &
                 ,'---------------------------------------------'
         write (o,'(4x,a,f7.4,a,i6,a,f6.2,a,i6,a/)') 'For freq.',w,' au,', &
                 int(w/hz/1d6),' MHz,',1/w/nm,' nm,',int(w/cm1), &
                 ' cm-1, with T in kelvin:'
         deta  = sum((/((Etn(i,j,i,j)-Etn(i,i,j,j)/3,i=1,3),j=1,3)/))/5
         dalxi = sum((/((Al(i,j)*Xn(i,j)-Al(i,i)*Xn(j,j)/3,i=1,3),j=1,3)/))
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'mC  = ',c1*deta,' + ', &
                 c1*c2*dalxi,' / T  (cm3 gauss-2 4pie0 mol-1)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'dn  = ',c3*deta,' + ', &
                 c3*c2*dalxi,' / T  (tesla-2 atm-1)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'Ccm = ',c4*w*deta,' / T + ', &
                 c4*c2*w*dalxi,' / T^2  (cm-1 gauss-2)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') '(calculated from deta =',deta, &
                 ' au, and dalpha*dxi =',dalxi,' au)'
      end if
      if (any(Xi/=0)) then
         write (o,'(//17x,a/17x,a/)') 'London static Magnetizability (au)', &
                                      '----------------------------------'
         write (o,'(19x,a2,18x,a2,18x,a2/)') bc
         write (o,'(3(1x,a,3f20.12/))') (bc(i),Xi(i,:),i=1,3)
         write (o,'(//17x,a/17x,a/)') 'London Hypermagnetizability - Eta (au)', &
                                      '--------------------------------------'
         write (o,'(4x,a,f6.3,a,i7,a,i4,a/)') 'For frequency w = ',w &
                 ,' au =',nint(w/cm1),' cm-1 ~ ',nint((2*pi*cvl/w)/nm),' nm'
         write (o,'(a,18x,a2,18x,a2,18x,a2)') ' w;-w',bc
         write (o,'(3(3(/3(1x,3a2,3f20.11/))))') &
                 (((fc(i),fc(j),bc(k),Eta(i,j,k,:),k=1,3),j=1,3),i=1,3)
         write (o,'(//16x,a/16x,a/)') 'London Cotton-mouton constant: mC, dn, Ccm', &
                                      '------------------------------------------'
         write (o,'(1x,a,f7.4,a,i6,a,f6.2,a,i6,a/)') 'For freq.',w,' au,', &
                 int(w/hz/1d6),' MHz,',1/w/nm,' nm,',int(w/cm1), &
                 ' cm-1, with T in kelvin:'
         deta  = sum((/((Eta(i,j,i,j)-Eta(j,j,i,i)/3,i=1,3),j=1,3)/))/5
         dalxi = sum((/((Al(i,j)*Xi(i,j)-Al(i,i)*Xi(j,j)/3,i=1,3),j=1,3)/))
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'mC  = ',c1*deta,' + ', &
                 c1*c2*dalxi,' / T  (cm3 gauss-2 4pie0 mol-1)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'dn  = ',c3*deta,' + ', &
                 c3*c2*dalxi,' / T  (tesla-2 atm-1)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') 'Ccm = ',c4*w*deta,' / T + ', &
                 c4*c2*w*dalxi,' / T^2  (cm-1 gauss-2)'
         write (o,'(5x,a,g13.6,a,g13.6,a/)') '(calculated from deta =',deta, &
                 ' au, and dalpha*dxi =',dalxi,' au)'
      end if
   end subroutine


   subroutine jones_output(freq, islon, dip, GtenFmwB0Bw, ApriFmwB0Qw, &
                           GtenF0FmwB0Bw, ApriF0FmwB0Qw, unit)
      logical,    intent(in) :: islon
      complex(8), intent(in) :: freq, dip(3), GtenFmwB0Bw(3,3,3), ApriFmwB0Qw(3,3,6)
      complex(8), intent(in) :: GtenF0FmwB0Bw(3,3,3,3), ApriF0FmwB0Qw(3,3,3,6)
      integer, optional, intent(in) :: unit
      integer, parameter :: eps(3,3,3) = &
          reshape( (/0, 0, 0, 0, 0,-1, 0,+1, 0,             &
                     0, 0,+1, 0, 0, 0,-1, 0, 0,             &      
                     0,-1, 0,+1, 0, 0, 0, 0, 0/), (/3,3,3/) )
      integer, parameter :: ij(3,3) = &
          reshape( (/1, 2, 4, 2, 3, 5, 4, 5, 6/), (/3,3/) )
      complex(8)  :: Gfb, Afb, Gb(3), Ab(3)
      integer     :: i, j, k, l, uni
      character*6 :: lonnol
      !unit defaults to 6=stdout
      uni = 6
      if (present(unit)) uni = unit
      !label for printining
      lonnol = merge('London', 'no-Lon', islon)
      !print the raw tensors contributing to the temperature-dependent term
      !call print_tensor( (/3,3,3/), GtenFmwB0Bw,          &
      !                   lonnol//' d/db G-tensor = Efbb', &
      !                   (/-1,0,1/)*freq, unit=unit )
      !call print_tensor( (/3,3,6/), ApriFmwB0Qw,          &
      !                   lonnol//' d/db A-prime = -Efbq', &
      !                   (/-1,0,1/)*freq, unit=unit )
      !print the raw tensors contributing to the temperature-independent term
      !call print_tensor( (/3,3,3,3/), GtenF0FmwB0Bw,            &
      !                   lonnol//' d/df d/db G-tensor = Effbb', &
      !                   (/0,-1,0,1/)*freq, unit=unit )
      !call print_tensor( (/3,3,3,6/), ApriF0FmwB0Qw,            &
      !                   lonnol//' d/df d/db A-prime = -Effbq', &
      !                   (/0,-1,0,1/)*freq, unit=unit ) !colwidth=16, 
      !calculate temperature-dependent averages
      Gb = (/ ( sum( (/ ( 3*GtenFmwB0Bw(j,j,i)  &
                        + 3*GtenFmwB0Bw(i,j,j)  &
                        - 2*GtenFmwB0Bw(j,i,j), &
                     j = 1, 3) /) ), i = 1, 3) /)
      !gfortran mis-compiles these lines, so rewrite without sum()
      !Ab = (/ ( sum( (/ ( ( ( eps(k,j,i) * ApriFmwB0Qw(k,l,ij(l,j))  &
      !                      + eps(l,j,k) * ApriFmwB0Qw(l,k,ij(i,j)), &
      !              j = 1, 3), k = 1, 3), l = 1, 3) /) ), i = 1, 3) /)
      do i=1,3; Ab(i)=0
         do l=1,3; do k=1,3; do j=1,3
            Ab(i) = Ab(i) + eps(k,j,i) * ApriFmwB0Qw(k,l,ij(l,j))  &
                          + eps(l,j,k) * ApriFmwB0Qw(l,k,ij(i,j))
         end do; end do; end do
      end do
      !calculate temperature-independent averages
      Gfb = sum( (/ ( ( 3*GtenF0FmwB0Bw(i,j,j,i)  &
                      + 3*GtenF0FmwB0Bw(i,i,j,j)  &
                      - 2*GtenF0FmwB0Bw(i,j,i,j), &
                          i = 1, 3), j = 1, 3) /) )
      !gfortran mis-compiles these lines, so rewrite without sum()
      !Afb = sum( (/ ( ( ( ( eps(k,j,i) * ApriF0FmwB0Qw(i,k,l,ij(l,j))  &
      !                    + eps(l,j,k) * ApriF0FmwB0Qw(i,l,k,ij(i,j)), &
      !                   i = 1, 3), j = 1, 3), k = 1, 3), l = 1, 3) /) )
      Afb=0; do l=1,3; do k=1,3; do j=1,3; do i=1,3
         Afb = Afb + eps(k,j,i) * ApriF0FmwB0Qw(i,k,l,ij(l,j)) &
                   + eps(l,j,k) * ApriF0FmwB0Qw(i,l,k,ij(i,j))
      end do; end do; end do; end do 
      write (uni,'()') !two final newlines
      write (uni,'(6x,a,f14.11)')  'frequency     : ', dreal(freq)
      write (uni,'(6x,a,3g18.11)') 'dipole moment : ', dreal(dip)
      !print temperature-independent averages
      write (uni,'(6x,a/22x,g18.11)') lonnol // &
            ' average Gfb = sum_ij 3 Gijji + 3 Giijj - 2 Gijij', &
            dreal(Gfb)
      write (uni,'(6x,a/22x,g18.11)') lonnol // &
            ' average Afb = sum_ijkl e_kji Aikllj + e_ljk Ailkij', &
            dreal(Afb)
      write (uni,'(6x,a/22x,g18.11)') lonnol // &
            ' temp.-indep. Jones biref. A0 = Gfb - w/2 * 2/3 * Afb', &
            dreal(Gfb - freq/2 * 2/3 * Afb)
      !print temperature-dependent averages and invariant
      write (uni,'(6x,a/22x,3g18.11)') lonnol // &
            ' average Gb_i = sum_j 3 Gjji + 3 Gijj - 2 Gjij', &
            dreal(Gb)
      write (uni,'(6x,a/22x,3g18.11)') lonnol // &
            ' average Af_i = sum_jkl e_kji Aikllj + e_ljk Ailkij', &
            dreal(Ab)
      write (uni,'(6x,a/22x,3g18.11)') lonnol // &
            ' temp.-dep. Jones factors Gb_i + w/2 * 2/3 * Afb_i', &
            dreal(Gb - freq/2 * 2/3 * Ab)
      write (uni,'(6x,a/22x,g18.11)') lonnol // &
            ' temp.-dep. Jones biref. A1 = mu_i * (Gb_i - w/2 * 2/3 * Ab_i)', &
            dreal(sum( dip * (Gb - freq/2 * 2/3 * Ab) ))
      write (uni,'(/)') !two final newlines
   end subroutine

end module
