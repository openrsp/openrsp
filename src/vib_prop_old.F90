module vib_prop_old
!Calculation and outputting of vibrational (optical) properties.
!This file is currently common to Dalton and Dirac

   use matrix_defop
   use rsp_equations_old
   use prop_contribs_old

   ! ajt LSDALTON has replaced the (global) quit with lsquit
   !     with unit (lupri) as extra argument, which doesn't
   !     exist in DIRAC. For now, this macro gets around that.
#ifdef LSDALTON_ONLY
#define quit lsquit
#endif

   implicit none

   public roa_pol_Gpri_Aten_grads
   public cars_pol_shyp_polgra
   public vib_ana_polari
   public print_shypol
   public vibhyp_hyp_dipgra_polgra
   public vibshyp_shyp_dipg_polg_hypg

   private

   !physical constants for conversion to non-atomic units
   real(8), parameter :: cm1 = 1/219474.631371d0, & !1 centimeter-to-minus-one in au
                         cvl = 137.03599907d0,    & !c-the-velocity-of-light in au
                         nm  = 10/0.52917706d0,   & !1 nanometer in au
                         pi  = 3.14159265358979323846D0 !acos(-1d0)

   !field component lables for printing
   character(2) :: fc(3) = (/'Fx','Fy','Fz'/)

contains


   subroutine print_tensor(dims, tensor, title, freqs, colwidth, unit)
   !ajt This is used for printing response function tensors
   !TODO: Add comp. lables argument. Make space between blocks of rows
      integer,                intent(in) :: dims(:)
      complex(8),             intent(in) :: tensor(*) !(product(dims))
      character(*), optional, intent(in) :: title
      integer,      optional, intent(in) :: unit, colwidth
      complex(8),   optional, intent(in) :: freqs(:)
      integer      :: uni, siz, dec, cwid, i
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


   !> Load Hessian form DALTON.HES, mass-weight, diagonalize
   !> and sqrt eigenvalues for frequencies (w<0 means imaginary).
   !> ajt FIXME Does not work in LSDALTON.
   subroutine load_vib_modes(mol, nc, nq, w, Q)
#ifdef LSDALTON_ONLY
      use files
#define GPOPEN  LSOPEN
#define GPCLOSE LSCLOSE
#elif defined(BUILD_OPENRSP)
      use dalton_ifc
#endif
      !> reference to molecule, geometry, etc.
      type(prop_molcfg), intent(in)  :: mol
      !> number of (cartesian) coordinates
      integer, intent(in)            :: nc
      !> number of vibrational normal mode coordinates
      integer, intent(out)           :: nq
      !> resulting vibrational frequencies
      real(8), intent(out)           :: w(nc)
      !> resulting mass-weighted normal modes Q(:,:nq)
      real(8), intent(out)           :: Q(nc,nc)
      !-----------------------------------------
      real(8)       :: M(nc/3),Mtot,Z(nc/3),G(nc),Rc(nc),H(nc,nc),orig(3)=0
      integer       :: IS(nc/3),i,internal=0,negative=0,unit=-1
      integer       :: ic, jc, ierr
      real(8)       :: diff
      logical       :: ex
      character(22) :: fmt
#ifdef LSDALTON_ONLY
      call quit('Cannot call NUCLEI_ifc or GPINQ, only new integral code is compiled')
#else
      call NUCLEI_ifc(nc/3,Z,IS,G)
      call GPINQ('DALTON.HES','EXIST',ex)
#endif
      if (.not.ex) call quit('load_vib_modes: Hessian file DALTON.HES not found')
#ifdef BUILD_OPENRSP
      call GPOPEN(unit,'DALTON.HES','OLD',' ','FORMATTED',0,.false.)
#else
      call GPOPEN(unit,'DALTON.HES','OLD','FORMATTED')
#endif

      rewind(unit)

!     read number of coordinates
      read(unit, *, iostat=ierr) i
      if (ierr /= 0) then
         call quit('load_vib_modes: could not read number of coordinates')
      end if
      if (i /= nc) then
         call quit('load_vib_modes: DALTON.HES inconsistent with MOLECULE.INP')
      end if

!     read hessian
      do ic = 1, nc
         read(unit, *) !blank line
         do jc = 1, nc
            read(unit, *, iostat=ierr) H(jc, ic)
         end do
      end do
      if (ierr /= 0) then
         call quit('load_vib_modes: could not read Hessian on DALTON.HES')
      end if

!     read geometry
      read(unit, *) !blank line
      do ic = 1, nc
         read(unit, *, iostat=ierr) Rc(ic)
      end do
      if (ierr /= 0) then
         write(*, *) 'WARNING: could not verify DALTON.HES geometry - make sure it is correct'
      else
         do ic = 1, nc
            diff = dabs(G(ic) - Rc(ic))
!           do not compare geometry without numeric tolerance
            if (diff > 1.0d-3) then
               call quit('load_vib_modes: Different geometry in DALTON.HES and MOLECULE.INP; '// &
                         'one or more coordinates differ by more than 1.0d-3 bohr')
            end if
         end do
      end if

      call GPCLOSE(unit,'KEEP')
#ifdef LSDALTON_ONLY
      call quit('Cannot call VIBHES, VIBMAS, or VIBNOR - only new integral code is compiled')
#else
      call VIBHES(1,nc,G,nq,(/(0d0,i=1,nc)/),H,(/(0d0,i=1,6*nc)/) &
                 ,nint(Z),(/(0d0,i=1,2*nc*nc)/)); nq=nc-nq
      call VIBMAS(M,Mtot,IS,nint(Z),nc/3,G,orig,1)
      call VIBNOR(H,M,(/(0d0,i=1,nc)/),(/(0d0,i=1,nc*(nc+1)/2)/) &
                 ,w,Q,(/(0d0,i=1,nc*nc)/),(/(0d0,i=1,30*nc*nc+1)/) &
                 ,30*nc*nc+1,nc,internal,negative,.false.,1)
#endif
      w = sign(sqrt(abs(w)),w)
   end subroutine


   subroutine print_vib_table(nq,nc,ne,lab,wc,ec,tab,o)
      integer           :: nq,nc,ne,wc(nc),ec(ne+1),o
      character(*)      :: lab(nc) !labels
      real(8)           :: tab(nc,nq) !values to be printed
      integer           :: i,j,k,l,w(nc),d(nc) !col width, decimals
      character*(32*nc) :: p !for pre-printing into
      character*8       :: fmt
      !count digits before decimal point
      do j=1,nc
         w(j) = 1+floor(log10(maxval(abs(tab(j,:)))))
         if (w(j)<0.and.w(j)<-wc(j)+3) then !zero column
            w(j) = 0
         elseif (all(tab(j,:)==nint(tab(j,:)))) then !integer column
            d(j) = -1
            if (any(tab(j,:)<0)) w(j) = w(j)+1
            w(j) = max(w(j)+1,wc(j))
         else !real column
            w(j) = 3+max(1,w(j))
            d(j) = max(0,wc(j)-w(j))
            w(j) = max(w(j),wc(j))
         end if
      end do
      !header
      l = 0
      do j=1,nc
         if (w(j)==0) cycle !skip column
         write (p(l+1:),'(a)') lab(j) !lab(j)(1:w(j))
         l = l+w(j)
      end do
      write (o,'(a)') trim(p(:l))
      !dashes ------
      write (o,'(1x,a)') repeat('-',l-1)
      !the table, line by line
      do i=1,nq
         l = 0
         do j=1,nc
            if (w(j)==0) cycle !skip column
            if (d(j)==-1) then
               write (fmt(1:5),'(a,i2,a)') '(i',w(j),')'
               write (p(l+1:),fmt(1:5)) nint(tab(j,i))
            else
               write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
               write (p(l+1:),fmt) tab(j,i)
            end if
            l = l+w(j)
         end do
         write (o,'(a)') p(:l)
      end do
      !dashes ------
      write (o,'(1x,a)') repeat('-',l-1)
      !column sums for 'energy' columns
      l = max(4,sum(w(:ec(1)-1)))
      write (p(1:),'(a)') repeat(' ',l-3)//'sum'
      do j=ec(1),ec(ne+1)-1
         if (w(j)==0.or.d(j)==-1) cycle
         write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
         write (p(l+1:),fmt) sum(tab(j,:))
         l = l+w(j)
      end do
      write (o,'(a)') p(:l)
      !total of 'energy' columns
      l = max(6,sum(w(:ec(1)-1)))
      write (p(1:),'(a)') repeat(' ',l-5)//'total'
      k = 1
      do j=ec(1),ec(ne)
         if (w(j)==0.or.d(j)==-1) cycle
         if (j==ec(k)) then
            write (fmt,'(a,i2,a,i2,a)') '(f',w(j),'.',d(j),')'
            write (p(l+1:),fmt) sum(tab(ec(k):ec(k+1)-1,:))
            k = k+1
         else
            write (p(l+1:),'(a)') repeat(' ',w(j))
         end if
         l = l+w(j)
      end do
      write (o,'(a/)') trim(p(:l))
   end subroutine


   subroutine print_polari(dax,w,pol,o)
   !Print the polarizability tensor pol, along with it's isotropic
   !and dipole-dipole components. dax is the NORMALIZED dipol axis.
      real(8) :: dax(3),w,pol(3,3)
      integer :: o,i,j
      write (o,'(//26x,a/26x,a/)') 'Polarizability (au)' &
                                  ,'-------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(2(8x,a,i1,a,f6.3,a,i7,a,i5,a/))') &
               'w',0,' = ', w,' au =', nint(w/cm1),' cm-1 ~ ', &
                                  nint((2*pi*cvl/w)/nm),' nm', &
               'w',1,' = ',-w,' au =',-nint(w/cm1),' cm-1 ~ ', &
                                 -nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(a,13x,a2,18x,a2,18x,a2/)') ' w0 w1',fc
      write (o,'(3(1x,a,3f20.13/))') (fc(i),pol(i,:),i=1,3)
      write (o,'(a,g16.9)')  '@alpha    isotropic: ', &
               (pol(1,1)+pol(2,2)+pol(3,3))/3
      write (o,'(a,g16.9/)') '@alpha     dipole^2: ', &
               sum((/((pol(i,j)*dax(i)*dax(j),i=1,3),j=1,3)/))
   end subroutine


   subroutine print_shypol(w,dip,pol,hyp,shyp,o)
      real(8) :: w(4),dip(3),pol(3,3,4),hyp(3,3,3,6),shyp(3,3,3,3)
      integer :: o,i,j,k,l,m,ij
      real(8) :: a,b,dax(3)
      write (o,'(//28x,a/28x,a/)') 'Dipole moment (au)', &
                                   '------------------'
      write (o,'(1x,3(18x,a)//(3x,3f20.14/))') fc,dip
      if (sqrt(sum(dip*dip))<0.00001) then
         write (o,'(4x,a/)') 'The molecule has zero dipole moment' &
                           //' - using z as dipole axis'
         dax = (/0d0,0d0,1d0/)
      else
         dax = dip/sqrt(sum(dip*dip))
      end if
      do i=1,4; a=w(i)
         if (any(w(1:i-1)==a.or.w(1:i-1)==-a)) cycle
         if (all(pol(:,:,i)==0)) cycle
         call print_polari(dax,w(i),pol(:,:,i),o)
      end do
      do j=2,4; do i=1,j-1; a=w(i); b=w(j); ij=i+(j-1)*(j-2)/2
         !only print one beta for each unique freqency triple
         if (i==1.and.j==3.and.any(w(2)==(/-a-b,b/))) cycle
         if (i==2.and.j==3.and.any(w(1)==(/-a-b,a,b/))) cycle
         if (i==1.and.j==4.and.any(w(2)==(/-a-2*b,-a-b,-a,-b,0d0,b/))) cycle
         if (i==2.and.j==4.and.any(w(1)==(/-a-2*b,-a-b,-a,-b,0d0,b,a/))) cycle
         if (i==3.and.j==4.and.any(w(2)==(/-2*a-b,-a-2*b,-a-b,-a,-b,0d0,b,a/))) cycle
         if (all(hyp(:,:,:,ij)==0)) cycle
         write (o,'(//24x,a/24x,a/)') 'Hyperpolarizability (au)' &
                                     ,'------------------------'
         write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
         write (o,'(3(8x,a,i1,a,f6.3,a,i7,a,i5,a/))') &
                  'w',0,' = ',   a,' au =', nint(     a/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/a     )/nm),' nm',  &
                  'w',1,' = ',   b,' au =', nint(     b/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/b     )/nm),' nm',  &
                  'w',2,' = ',-a-b,' au =', nint((-a-b)/cm1),' cm-1 ~ ', &
                                      nint((2*pi*cvl/(-a-b))/nm),' nm'
         write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1 w2',fc
         write (o,'(3(/3(1x,2a2,3f20.12/)))') &
                 ((fc(k),fc(l),hyp(k,l,:,ij),l=1,3),k=1,3)
         write (o,'(a,g16.9)') '@beta     dip-iso: ',1/5d0*sum((/(( &
                 (hyp(k,k,l,ij)+hyp(k,l,k,ij)+hyp(l,k,k,ij))*dax(l), &
                  k=1,3),l=1,3)/))
         write (o,'(a,g16.9/)') '@beta    dipole^3: ', &
                  sum((/(((hyp(k,l,m,ij)*dax(k)*dax(l)*dax(m),k=1,3),l=1,3),m=1,3)/))
      end do; end do
      if (all(shyp==0)) return
      write (o,'(//22x,a/22x,a/)') 'Second hyperpolarizability (au)', &
                                   '-------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(4(8x,a,i1,a,f6.3,a,i7,a,i5,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,4)
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1w2 w3',fc
      write (o,'(3(3(/3(1x,3a2,3f20.10/))))') &
              (((fc(i),fc(j),fc(k),shyp(i,j,k,:),k=1,3),j=1,3),i=1,3)
      write (o,'(a,g16.9)') '@gamma    isotropic: ',1/15d0*sum((/(( &
              shyp(i,i,j,j)+shyp(i,j,i,j)+shyp(i,j,j,i),i=1,3),j=1,3)/))
      write (o,'(a,g16.9)') '@gamma    dip^2-iso: ',1/10d0*sum((/(((     &
              (shyp(i,i,k,l)+shyp(i,k,i,l)+shyp(k,i,i,l)                 &
              +shyp(i,k,l,i)+shyp(k,i,l,i)+shyp(i,i,k,l))*dax(k)*dax(l), &
               k=1,3),l=1,3),i=1,3)/))
      write (o,'(a,g16.9)') '@gamma        dip^4: ',                  &
               1*sum((/((((shyp(i,j,k,l)*dax(i)*dax(j)*dax(k)*dax(l), &
                                      i=1,3),j=1,3),k=1,3),l=1,3)/))
      write (o,'(4(a,1x,g13.7/)/)')                                      &
               '@gamma    Isotropic (AABB+ABAB+ABBA)/15:',               &
               sum((/((shyp(i,i,j,j)+shyp(i,j,i,j)+shyp(i,j,j,i),        &
                                            i=1,3),j=1,3)/))/15,         &
               '@gamma                   XXXX-component:',shyp(1,1,1,1), &
               '@gamma                   YYYY-component:',shyp(2,2,2,2), &
               '@gamma                   ZZZZ-component:',shyp(3,3,3,3)
   end subroutine


   !> Print vibrational contributions to polarizabilities
   subroutine vib_ana_polari(mol, w, dip, nc, dipg, polg, hypg, o)
      !> reference to molecule, geometry, etc., needed by load_vib_modes
      type(prop_molcfg), intent(in)  :: mol
      !> frequencies of 4 electric fields
      real(8),           intent(in) :: w(4)
      !> dipole moment, used to determine dipole axis
      real(8),           intent(in) :: dip(3)
      !> number of cartesian geometrical coordinates
      integer,           intent(in) :: nc
      !> dipole gradient (up to 4 of them)
      real(8),           intent(in) :: dipg(3,nc,4)
      !> polarizability gradient (up to 6)
      real(8),           intent(in) :: polg(3,3,nc,6)
      !> hyperpolarizability gradient (up to 4)
      real(8),           intent(in) :: hypg(3,3,3,nc,4)
      !> printing unit
      integer,           intent(in) :: o
      !---------------------------------
      integer :: nq,i,j,k,l,m,n,ij,ik,jk
      real(8) :: dax(3),dipq(3,nc*4),polq(3,3,nc*6),hypq(3,3,3,nc*4),wq(nc),Q(nc,nc)
      if (sqrt(sum(dip*dip))< 0.00001) dax = (/0d0,0d0,1d0/)
      if (sqrt(sum(dip*dip))>=0.00001) dax = dip/sqrt(sum(dip*dip))
      call load_vib_modes(mol,nc,nq,wq,Q)
      write (o,'(//20x,a/20x,a/)') 'Mass-weighted normal coordinates (au)', &
                                   '-------------------------------------'
      call OUTPUT(Q(:,1:nq),1,nc,1,nq,nc,nq,1,o)
      do i=1,4; call DGEMM('N','N',3,nq,nc,1d0,dipg(:,:,i), &
              3,Q(:,1:nq),nc,0d0,dipq(:,1+nq*(i-1):nq*i),3); end do
      do i=1,6; call DGEMM('N','N',3*3,nq,nc,1d0,polg(:,:,:,i), &
              3*3,Q(:,1:nq),nc,0d0,polq(:,:,1+nq*(i-1):nq*i),3*3); end do
      do i=1,4; call DGEMM('N','N',3*3*3,nq,nc,1d0,hypg(:,:,:,:,i), &
              3*3*3,Q(:,1:nq),nc,0d0,hypq(:,:,:,1+nq*(i-1):nq*i),3*3*3); end do
      do i=1,4
         if (any(w(1:i-1)==w(i).or.w(1:i-1)==-w(i))) cycle
         if (all(dipq(:,1+nq*(i-1):nq*i)==0)) cycle
         call print_vib_alpha(dax,w(i),nq,wq(1:nq),dipq(:,1+nq*(i-1):nq*i),o)
      end do
      do n=4,1,-1; if (w(n)/=0) cycle
         i=merge(1,2,n>1); j=merge(2,3,n>2); k=merge(3,4,n>3)
         if (all(dipq(:,1+nq*(i-1):nq*i)==0).and. &
             all(dipq(:,1+nq*(j-1):nq*j)==0).and. &
             all(dipq(:,1+nq*(k-1):nq*k)==0)) cycle
         ij=i+(j-1)*(j-2)/2; ik=i+(k-1)*(k-2)/2; jk=j+(k-1)*(k-2)/2
         if (all(polq(:,:,1+nq*(ij-1):nq*ij)==0).and. &
             all(polq(:,:,1+nq*(ik-1):nq*ik)==0).and. &
             all(polq(:,:,1+nq*(jk-1):nq*jk)==0)) cycle
         call print_vib_beta(dax,(/w(i),w(j),w(k)/),nq,wq(1:nq), &
                 -(/dipq(:,1+nq*(i-1):nq*i),dipq(:,1+nq*(j-1):nq*j),dipq(:,1+nq*(k-1):nq*k)/), &
                 -(/polq(:,:,1+nq*(ij-1):nq*ij),polq(:,:,1+nq*(ik-1):nq*ik), &
                   polq(:,:,1+nq*(jk-1):nq*jk)/),o)
         exit !only one
      end do
      if ((any(dipq(:,1:nq*4) /= 0) .and. &
           any(hypq(:,:,:,1:nq*4) /= 0)) .or. &
          (any(polq(:,:,1:nq*3) /= 0) .and. &
           any(polq(:,:,nq*3+1:nq*6) /= 0))) &
         call print_vib_gamma(dax,w,nq,wq(1:nq),dipq(:,1:nq*4),    &
                              polq(:,:,1:nq*6),hypq(:,:,:,1:nq*4),o)
   end subroutine


   subroutine print_vib_alpha(dax,w,nq,wq,dipg,o)
      real(8) :: dax(3),w,wq(nq),dipg(3,nq)
      integer :: o,nq,i,j,k,m,q
      real(8) :: pol(3,3),R(3,3),t(5,nq)
      write (o,'(//21x,a/21x,a/)') 'Vibrational polarizability (au)', &
                                   '-------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(2(8x,a,i1,a,f6.3,a,i7,a,i5,a/))')             &
               'w',0,' = ', w,' au =', nint(w/cm1),' cm-1 ~ ', &
                                  nint((2*pi*cvl/w)/nm),' nm', &
               'w',1,' = ',-w,' au =',-nint(w/cm1),' cm-1 ~ ', &
                                 -nint((2*pi*cvl/w)/nm),' nm'
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      pol = 0
      do q=1,nq
         R(:,:) = reshape((/((dipg(j,q)*dipg(i,q),i=1,3),j=1,3)/),(/3,3/))
         pol = pol + R(:,:)/(wq(q)**2-w**2)
         t(1:3,q) = (/1d0*q,wq(q)/cm1,wq(q)/)
         t(4:5,q) = (/1d0/(wq(q)**2-w**2),1d0/)*1/3d0*sum((/(R(i,i),i=1,3)/))
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2/)') ' w0 w1',fc
      write (o,'(3(1x,a2,3f20.13/))') (fc(i),pol(i,:),i=1,3)
      write (o,'(a,g16.9)') '@alpha_vib    isotropic:',sum(t(4,:))
      write (o,'(a,g16.9/)') '@alpha_vib     dipole^2:', &
               sum((/((pol(i,j)*dax(i)*dax(j),i=1,3),j=1,3)/))
      if (abs(w)>0.03) t(5,:)=0
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'mm   :  dmu/dQ x dmu/dQ contribution (isotropic)'
      if (abs(w)<=0.03) write (o,'(8x,a)') 'm1   :  <<Q;F>>_w1 induced gradient with freq -w1 =w0'
      write (o,*)
      call print_vib_table(nq,5,1,(/'         ','freq/cm-1',   &
                        '  freq/au','    mm   ','    m0m1 '/), &
                        (/3,9,10,11,11/),(/4,5/),t,o)
   end subroutine


   subroutine print_vib_beta(dax,w,nq,wq,dipg,polg,o)
      real(8) :: dax(3),w(3),wq(nq),dipg(3,nq,3),polg(3,3,nq,3)
      real(8) :: hyp(3,3,3),R(3,3,3,4),t(15,nq)
      integer :: o,nq,i,j,k,n,q
      write (o,'(//20x,a/20x,a/)') 'Vibrational hyperpolarizability (au)', &
                                   '------------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(3(8x,a,i1,a,f6.3,a,i7,a,i5,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,3)
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      hyp = 0
      do q=1,nq
         R(:,:,:,2:4) = reshape((/ &
                 (((dipg(k,q,3)*polg(i,j,q,1),i=1,3),j=1,3),k=1,3), &
                 (((dipg(j,q,2)*polg(i,k,q,2),i=1,3),j=1,3),k=1,3), &
                 (((dipg(i,q,1)*polg(j,k,q,3),i=1,3),j=1,3),k=1,3)/),(/3,3,3,3/))
         R(:,:,:,1) = R(:,:,:,2)/(wq(q)**2-w(3)**2) &
                   + R(:,:,:,3)/(wq(q)**2-w(2)**2) &
                   + R(:,:,:,4)/(wq(q)**2-w(1)**2)
         hyp = hyp + R(:,:,:,1)
         !ajt gfortran seems to compiles this wrongly, ending up going
         !    out of range in k. So use explicit loop over k instead
         !t(1: 4,q) = (/(sum((/(((R(i,i,j,k)+R(i,j,i,k)+R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         !t(5: 8,q) = (/(sum((/(((R(i,i,j,k)-R(i,j,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         !t(9:12,q) = (/(sum((/(((R(i,i,j,k)-R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5,k=1,4)/)
         do k = 1, 4
            t(0+k,q) = sum((/(((R(i,i,j,k)+R(i,j,i,k)+R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5
            t(4+k,q) = sum((/(((R(i,i,j,k)-R(i,j,i,k))*dax(j),i=1,3),j=1,3)/))/5
            t(8+k,q) = sum((/(((R(i,i,j,k)-R(j,i,i,k))*dax(j),i=1,3),j=1,3)/))/5
         end do
         t(:,q) = (/1d0*q,wq(q)/cm1,wq(q),t(1,q),t(5,q),t(9,q),t(2:4,q),t(6:8,q),t(10:12,q)/)
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1 w2',Fc
      write (o,'(3(/3(1x,2a2,3f20.12/)))') ((Fc(i),Fc(j),hyp(i,j,:),j=1,3),i=1,3)
      write (o,'(a,g16.9)') '@beta_vib    dip-isotropic: ',sum(t(4,:))
      write (o,'(a,g16.9/)') '@beta_vib         dipole^3: ', &
               sum((/(((hyp(i,j,k)*dax(i)*dax(j)*dax(k),i=1,3),j=1,3),k=1,3)/))
      if (w(3)==0.or.maxval(abs(t(5,:)))<1d-8*maxval(abs(t(4,:)))) then
          t(5,:)=0; t(10:12,:)=0; end if
      if (w(2)==0.or.maxval(abs(t(6,:)))<1d-8*maxval(abs(t(4,:)))) then
          t(6,:)=0; t(13:15,:)=0; end if
      if (w(3)==w(2).or.w(3)==w(1).or.abs(w(3))>0.03) t(7,:)=0
      if (w(2)==w(1).or.abs(w(2))>0.03) t(8,:)=0
      if (abs(w(1))>0.03) t(9,:)=0
      do i=7,9; if (w(3)==0.and.all(t(i,:)==0)) t(i+3,:)=0; end do
      do i=7,9; if (w(2)==0.and.all(t(i,:)==0)) t(i+6,:)=0; end do
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'ma   :  dmu/dQ x dalpha/dQ contribution (dipole-isotropic)'
      if (any(t(5,:)/=0)) write (o,'(8x,a)') "ma'  :  as ma, but anisotropic: sum(iiz-izi)/5"
      if (any(t(6,:)/=0)) write (o,'(8x,a)') 'ma"  :  as ma, but anisotropic: sum(iiz-zii)/5'
      write (o,'(8x,a)') 'm2   :  <<Q;F>>_w2 induced gradient with freq -w2 =w0+w1'
      write (o,'(8x,a/)') 'a01  :  <<Q;F,F>>_w0,w1 ind. grad. with freq -w0-w1 =w2'
      call print_vib_table(nq,15,3,                                             &
               (/'         ','freq/cm-1','  freq/au','    ma   ',"    ma'  ",   &
                 '    ma"  ','   m2a01 ','   m1a02 ','   m0a12 ',"   m2a01'",   &
                 "   m1a02'","   m0a12'",'   m2a01"','   m1a02"','   m0a12"'/), &
               (/3,9,10,(11,i=1,12)/),(/4,5,6,7/),t,o)
   end subroutine


   subroutine print_vib_gamma(dax,w,nq,wq,dipq,polq,hypq,o)
      real(8) :: dax(3),w(4),wq(nq),dipq(3,nq,4),polq(3,3,nq,6),hypq(3,3,3,nq,4)
      real(8) :: shyp(3,3,3,3),R(3,3,3,3,9),t(30,nq)
      integer :: o,nq,i,j,k,l,n,q
      write (o,'(//19x,a/19x,a/)') 'Vibrational second hyperpolarizability (au)', &
                                   '-------------------------------------------'
      write (o,'(4x,a/)') 'Optical frequencies and wavelengths:'
      write (o,'(4(8x,a,i1,a,f6.3,a,i7,a,i5,a/))')                   &
              ('w',i-1,' = ',w(i),' au =',nint(w(i)/cm1),' cm-1 ~ ', &
                             nint((2*pi*cvl/w(i))/nm),' nm',i=1,4)
      write (o,'(4x,a,3f8.4/)') 'Normalized dipole axis: ',dax
      shyp = 0
      do q=1,nq
         R(:,:,:,:,3:9) = reshape((/ &
                 ((((polq(j,k,q,3)*polq(i,l,q,4),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((polq(i,k,q,2)*polq(j,l,q,5),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((polq(i,j,q,1)*polq(k,l,q,6),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(l,q,4)*hypq(i,j,k,q,1),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(k,q,3)*hypq(i,j,l,q,2),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(j,q,2)*hypq(i,k,l,q,3),i=1,3),j=1,3),k=1,3),l=1,3), &
                 ((((dipq(i,q,1)*hypq(j,k,l,q,4),i=1,3),j=1,3),k=1,3),l=1,3)/), &
                 (/3,3,3,3,7/))
         R(:,:,:,:,1) = R(:,:,:,:,3)/(wq(q)**2-(w(2)+w(3))**2) &
                      + R(:,:,:,:,4)/(wq(q)**2-(w(1)+w(3))**2) &
                      + R(:,:,:,:,5)/(wq(q)**2-(w(1)+w(2))**2)
         R(:,:,:,:,2) = R(:,:,:,:,6)/(wq(q)**2-w(4)**2) &
                      + R(:,:,:,:,7)/(wq(q)**2-w(3)**2) &
                      + R(:,:,:,:,8)/(wq(q)**2-w(2)**2) &
                      + R(:,:,:,:,9)/(wq(q)**2-w(1)**2)
         shyp = shyp + R(:,:,:,:,1) + R(:,:,:,:,2)
         t( 1: 9,q) = (/(sum((/((R(i,i,j,j,k)+R(i,j,i,j,k)+R(i,j,j,i,k),i=1,3),j=1,3)/)),k=1,9)/)/15
         t(10:18,q) = (/(sum((/((R(i,i,j,j,k)-R(i,j,j,i,k),i=1,3),j=1,3)/)),k=1,9)/)/10
         t(19:27,q) = (/(sum((/((R(i,i,j,j,k)+R(i,j,j,i,k)-2*R(i,j,i,j,k),i=1,3),j=1,3)/)),k=1,9)/)/10
         t(:,q) = (/1d0*q,wq(q)/cm1,wq(q),t(1:2,q),t(10:11,q), &
                    t(19:20,q),t(3:9,q),t(12:18,q),t(21:27,q)/)
      end do
      write (o,'(a,13x,a2,18x,a2,18x,a2)') ' w0w1w2 w3',Fc
      write (o,'(3(3(/3(1x,3a2,3f20.10/))))') (((Fc(k),Fc(j),Fc(i),shyp(k,j,i,:),i=1,3),j=1,3),k=1,3)
      write (o,'(a,g16.9)') '@gamma_vib    isotropic: ',sum(t(4:5,:))
      write (o,'(a,g16.9)') '@gamma_vib    dip^2-iso: ',                           &
               sum((/((((shyp(i,i,k,l)+shyp(i,k,i,l)+shyp(k,i,i,l)                 &
                        +shyp(i,k,l,i)+shyp(k,i,l,i)+shyp(i,i,k,l))*dax(k)*dax(l), &
                                                     k=1,3),l=1,3),i=1,3)/)) / 10
      write (o,'(a,g16.9)') '@gamma_vib        dip^4: ',            &
               sum((/((((shyp(i,j,k,l)*dax(i)*dax(j)*dax(k)*dax(l), &
                                   i=1,3),j=1,3),k=1,3),l=1,3)/))
      write (o,'(4(a,1x,g13.7/))')                                           &
               '@gamma_vib    Isotropic (AABB+ABAB+ABBA)/15:',sum(t(4:5,:)), &
               '@gamma_vib                   XXXX-component:',shyp(1,1,1,1), &
               '@gamma_vib                   YYYY-component:',shyp(2,2,2,2), &
               '@gamma_vib                   ZZZZ-component:',shyp(3,3,3,3)
      !zero (skip) repeated and high-frequency columns
      if (maxval(abs(t(6:7,:)))<1d-8*maxval(abs(t(4:5,:)))) then
          t(6:7,:)=0; t(17:23,:)=0; end if
      if (maxval(abs(t(8:9,:)))<1d-8*maxval(abs(t(4:5,:)))) then
          t(8:9,:)=0; t(24:30,:)=0; end if
      if (w(4)==w(3).or.w(4)==w(2).or.abs(w(2)+w(3))>0.03) t(10,:)=0
      if (w(3)==w(2).or.abs(w(1)+w(3))>0.03) t(11,:)=0
      if (abs(w(2)+w(3))>0.03) t(12,:)=0
      if (w(4)==w(3).or.w(4)==w(2).or.w(4)==w(1).or.abs(w(4))>0.03) t(13,:)=0
      if (w(3)==w(2).or.w(3)==w(1).or.abs(w(3))>0.03) t(14,:)=0
      if (w(2)==w(1).or.abs(w(2))>0.03) t(15,:)=0
      if (abs(w(1))>0.03) t(16,:)=0
      if (w(2)==w(3).or.w(2)==w(4)) then; do i=17,23
              if (all(t(i- 7,:)==0)) t(i,:)=0; end do; end if
      if (w(2)==w(4).or.w(1)==w(3)) then; do i=24,30
              if (all(t(i-14,:)==0)) t(i,:)=0; end do; end if
      !print the table
      write (o,'(4x,a/)') 'Table of mode contributions and dominant residues (au):'
      write (o,'(8x,a)') 'aa   :  dalpha/dQ x dalpha/dQ contribution (isotropic)'
      if (any(t(6:7,:)/=0)) write (o,'(8x,a)') "aa'  :  as aa, but anisotropic: sum(iijj-ijji)/10"
      if (any(t(8:9,:)/=0)) write (o,'(8x,a)') 'aa"  :  as aa, but anisotropic: sum(iijj-2ijij+ijji)/10'
      write (o,'(8x,a)') 'mb   :  dmu/dQ x dbeta/dQ  contribution (isotropic)'
      write (o,'(8x,a)') 'm3   :  <<Q;F>>_w3 induced gradient with freq -w3 =w0+w1+w2'
      write (o,'(8x,a)') 'a01  :  <<Q;F,F>>_w0,w1 ind. grad. with freq -w0-w1 =w2+w3'
      write (o,'(8x,a/)') 'b123 :  <<Q;F,F,F>>_w1,w2,w3 ind. grad. with freq -w1-w2-w3 =w0'
      call print_vib_table(nq,30,3,                                                         &
               (/'         ','freq/cm-1','  freq/au','    aa   ','    mb   ',"    aa'  ",   &
                 "    mb'  ",'    aa"  ','    mb"  ','  a12a03 ','  a02a13 ','  a01a23 ',   &
                 '  m3b012 ','  m2b013 ','  m1b023 ','  m0b123 ',"  a12a03'","  a02a13'",   &
                 "  a01a23'","  m3b012'","  m2b013'","  m1b023'","  m0b123'",'  a12a03"',   &
                 '  a02a13"','  a01a23"','  m3b012"','  m2b013"','  m1b023"','  m0b123"'/), &
               (/3,9,10,(11,i=1,27)/),(/4,6,8,10/),t,o)
   end subroutine


   !> For frequency w=freq calculate polarizability -Eff, G-prime Ebf(w,-w),
   !> A-tensor -Eqf, pol. gradient Egff, G-prime gradient Egbf(0,w,-w),
   !> A-tensor gradient Egqf
   subroutine roa_pol_Gpri_Aten_grads(mol, S, D, F, ng, freq, &
                                      Eff, Ebf, Eqf, Egff, Egbf, Egqf)
      !> reference to molecule, solver- and integral settings
      type(prop_molcfg), intent(in)  :: mol
      !> overlap matrix
      type(matrix),      intent(in)  :: S
      !> density matrix
      type(matrix),      intent(in)  :: D
      !> Fock matrix
      type(matrix),      intent(in)  :: F
      !> number of geometrical coordinates
      integer,           intent(in)  :: ng
      !> external laser frequency
      complex(8),        intent(in)  :: freq
      !> polarizability is returned in Eff, G' in Ebf, A-tensor in Eqf,
      !> and the respective gradients in Egff, Egbf and Egqf
      !> ajt FIXME: Sign-name confusion Eff should be -alpha
      complex(8),        intent(out) :: Eff(3,3), Ebf(3,3), Eqf(6,3), &
                                        Egff(ng,3,3),Egbf(ng,3,3),   &
                                        Egqf(ng,6,3)
      type(matrix) :: De(3), Df(3), Db(3), Dq(6), Def(3,3), Dbf(3,3), Dqf(6,3), &
                      Fe(3), Ff(3), Fb(3), Fq(6), Fef(3,3), Fbf(3,3), Fqf(6,3), &
                      DFDef(3,3), DFDbf(3,3), DFDqf(6,3)
      integer      :: i, j, k
      ! electric equations, polarizability, G-prime, A-tensor
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/-freq/))
      Eff=0; Ebf=0; Eqf=0
      call prop_oneave(mol, S, (/'EL'/), (/Df/), shape(Eff), Eff)
      call prop_oneave(mol, S, (/'MAGO'/),(/Df/), shape(Ebf), Ebf)
      call prop_oneave(mol, S, (/'ELGR'/), (/Df/), shape(Eqf), Eqf)
      call print_tensor(shape(Eff), -Eff, 'Alpha = -Eff', (/freq,-freq/))
      call print_tensor(shape(Ebf),  Ebf, 'G-prime = Ebf', (/freq,-freq/))
      call print_tensor(shape(Eqf), -Eqf, 'A-tensor = -Eqf', (/freq,-freq/))
      !--------------------------------------------------
      ! electric-electric equations, solve 6 eqs, not 3x3
      do k = 1, 3;
         De(k) = trps(Df(k))
         Fe(k) = trps(Ff(k))
         call pert_dens(mol, S, (/'EL','EL'/), (/k,1/), &
                        (/D,De(1:k),Df(k)/), (/F,Fe(1:k),Ff(k)/), &
                        Def(1:k,k), Fef(1:k,k), freq=(/freq,-freq/), comp=(/1,k/))
         do j = 1, k-1
            Def(k,j) = trps(Def(j,k))
            Fef(k,j) = trps(Fef(j,k))
         end do
      end do
      !------------------------
      ! polarizability gradient
      Egff = 0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/Df/), shape(Egff), Egff)
      ! call print_tensor(Egff, shape(Egff), 'E1geDf'); Egff=0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/De/), shape(Egff), Egff, perm=(/1,3,2/))
      ! call print_tensor(Egff, shape(Egff), 'E1gfDe'); Egff=0
      call prop_twoave(mol, (/'GEO'/), (/D,De,Df,Def/), shape(Egff), Egff)
      do k = 1, 3
         do j = 1, 3
            DFDef(j,k) = De(j)*Ff(k)*D + De(j)*(F+freq*S)*Df(k) + D*Fe(j)*Df(k) &
                       + Df(k)*Fe(j)*D + Df(k)*(F-freq*S)*De(j) + D*Ff(k)*De(j) &
                       + Def(j,k)*F*D + D*Fef(j,k)*D + D*F*Def(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'GEO'/), (/Def/), shape(Egff), Egff, DFD=(/DFDef/))
      DFDef = 0
      ! call print_tensor(Egff, shape(Egff), 'E1gDef+E2gDeDf-SgDFDef'); Egff=0
      call print_tensor(shape(Egff), -Egff, 'd/dg Polari = -Egff', &
                        (/0*freq,freq,-freq/))
      DFDef=0; De=0; Fe=0; Def=0; Fef=0
      !----------------------------------------
      ! magetic and electric-magnetic equations
      call pert_dens(mol, S, (/'MAGO'/), (/3/), (/D/), (/F/), Db, Fb, freq=(/freq/))
      call pert_dens(mol, S, (/'MAGO','EL  '/), (/3,3/), (/D,Db,Df/), (/F,Fb,Ff/), &
                     Dbf, Fbf, freq=(/freq,-freq/))
      !-----------------
      ! G-prime gradient
      Egbf = 0
      call prop_oneave(mol, S, (/'GEO ','MAGO'/), (/Df/), shape(Egbf), Egbf)
      ! call print_tensor(Egbf, shape(Egbf), 'E1gbDf'); Egbf=0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/Db/), shape(Egbf), Egbf, perm=(/1,3,2/))
      ! call print_tensor(Egbf, shape(Egbf), 'E1gfDb'); Egbf=0
      call prop_twoave(mol, (/'GEO'/), (/D,Db,Df,Dbf/), shape(Egbf), Egbf)
      do k = 1, 3
         do j = 1, 3
            DFDbf(j,k) = Db(j)*Ff(k)*D + Db(j)*(F+freq*S)*Df(k) + D*Fb(j)*Df(k) &
                       + Df(k)*Fb(j)*D + Df(k)*(F-freq*S)*Db(j) + D*Ff(k)*Db(j) &
                       + Dbf(j,k)*F*D + D*Fbf(j,k)*D + D*F*Dbf(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'GEO'/), (/Dbf/), shape(Egbf), Egbf, DFD=(/DFDbf/))
      ! call print_tensor(Egbf, shape(Egbf), 'E2gDbDf + E1gDfb - SgDFDbf'); Egbf=0
      call print_tensor(shape(Egbf), Egbf, 'd/dg G-prime = Egbf', &
                        (/0*freq,freq,-freq/))
      DFDbf=0; Db=0; Fb=0; Dbf=0; Fbf=0
      !---------------------------------------------
      ! quadrupole and quadrupole-electric equations
      call pert_dens(mol, S, (/'ELGR'/), (/6/), (/D/), (/F/), Dq, Fq, freq=(/freq/))
      call pert_dens(mol, S, (/'ELGR','EL  '/), (/6,3/), (/D,Dq,Df/), (/F,Fq,Ff/), &
                     Dqf, Fqf, freq=(/freq,-freq/))
      !------------------
      ! A-tensor gradient
      Egqf = 0
      call prop_oneave(mol, S, (/'GEO ','ELGR'/), (/Df/), shape(Egqf), Egqf)
      ! call print_tensor(shape(Egqf), Egqf, 'E1gqDf'); Egqf=0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/Dq/), shape(Egqf), Egqf, perm=(/1,3,2/))
      ! call print_tensor(shape(Egqf), Egqf, 'E1gfDq'); Egqf=0
      call prop_twoave(mol, (/'GEO'/), (/D,Dq,Df,Dqf/), shape(Egqf), Egqf)
      do k = 1, 3
         do j = 1, 6
            DFDqf(j,k) = Dq(j)*Ff(k)*D + Dq(j)*(F+freq*S)*Df(k) + D*Fq(j)*Df(k) &
                       + Df(k)*Fq(j)*D + Df(k)*(F-freq*S)*Dq(j) + D*Ff(k)*Dq(j) &
                       + Dqf(j,k)*F*D + D*Fqf(j,k)*D + D*F*Dqf(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'GEO'/), (/Dqf/), shape(Egqf), Egqf, DFD=(/DFDqf/))
      ! call print_tensor(shape(Egqf), Egqf, 'E2gDqDf + E1gDqf - SgDFDqf'); Egqf=0
      call print_tensor(shape(Egqf), -Egqf, 'd/dg A-tensor = -Egqf', &
                        (/0*freq,freq,-freq/))
      DFDqf=0; Dq=0; Fq=0; Df=0; Ff=0; Dqf=0; Fqf=0
      ! change sign on quasi-energy derivatives to get properties
      Eff=-Eff; Eqf=-Eqf; Egff=-Egff; Egqf=-Egqf
   end subroutine


   !> For frequency freq (w), calculates polarizability the alpha(-w,w) in Eff,
   !> polarizability gradient d/dg alpha(-w,w) in Eff, and 2nd hyperpolarizability
   !> gamma(-w,w,-w,w) in Effff
   subroutine cars_pol_shyp_polgra(mol, S, D, F, ng, freq, Ef, Eff, Effff, Egff)
      !> reference to molecule, solver- and integral settings
      type(prop_molcfg), intent(in)  :: mol
      !> overlap matrix
      type(matrix),      intent(in)  :: S
      !> density matrix
      type(matrix),      intent(in)  :: D
      !> Fock matrix
      type(matrix),      intent(in)  :: F
      !> number of geometrical coordinates
      integer,           intent(in)  :: ng
      !> external laser frequency
      complex(8),        intent(in)  :: freq
      !> dipole moment returned in Ef, polarizability in in Eff
      !> 2nd hyperpolarizability in Effff and polarizability gradient in Egff
      !> ajt FIXME: Sign-name confusion Ef should be -dipmom
      complex(8),        intent(out) :: Ef(3), Eff(3,3), Effff(3,3,3,3), &
                                        Egff(ng,3,3)
      type(matrix) :: De(3), Df(3), Def(3,3), Dff(3,3), FeDS(3), FDSfef, DSDfef
      type(matrix) :: Fe(3), Ff(3), Fef(3,3), Fff(3,3), DeSD(3), DFDef(3,3)
      integer      :: i, j, k, l
      !dipole moment Ef=-mu
      Ef = 0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), Ef)
      call print_tensor((/3/), -Ef, 'dipmom = -Ef')
      !electric equations, Eff=-alpha
      call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), Df, Ff, freq=(/freq/))
      do i = 1, 3
         De(i) = trps(Df(i))
         Fe(i) = trps(Ff(i))
      end do
      Eff = 0
      call prop_oneave(mol, S, (/'EL'/), (/De/), (/3,3/), Eff)
      call print_tensor((/3,3/), -Eff, 'Alpha = -Eff', (/-freq,freq/))
      !---------------------------------------------------------------------
      ! electric-electric equations De(-w)f(+w), solve 6 instead of 3x3 eqs.
      do k = 1, 3;
         De(k) = trps(Df(k))
         Fe(k) = trps(Ff(k))
         call pert_dens(mol, S, (/'EL','EL'/), (/k,1/), &
                        (/D,De(:k),Df(k)/), (/F,Fe(:k),Ff(k)/), &
                        Def(:k,k), Fef(:k,k), freq=(/-freq,freq/), comp=(/1,k/))
         do j = 1, k-1
            Def(k,j) = trps(Def(j,k))
            Fef(k,j) = trps(Fef(j,k))
         end do
      end do
      !------------------------
      ! polarizability gradient
      Egff = 0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/Df/), shape(Egff), Egff)
      ! call print_tensor(Egff, shape(Egff), 'E1geDf'); Egff=0
      call prop_oneave(mol, S, (/'GEO','EL '/), (/De/), shape(Egff), Egff, perm=(/1,3,2/))
      ! call print_tensor(Egff, shape(Egff), 'E1gfDe'); Egff=0
      call prop_twoave(mol, (/'GEO'/), (/D,De,Df,Def/), shape(Egff), Egff)
      do k = 1, 3
         do j = 1, 3
            DFDef(j,k) = De(j)*Ff(k)*D + De(j)*(F-freq*S)*Df(k) + D*Fe(j)*Df(k) &
                       + Df(k)*Fe(j)*D + Df(k)*(F+freq*S)*De(j) + D*Ff(k)*De(j) &
                       + Def(j,k)*F*D + D*Fef(j,k)*D + D*F*Def(j,k)
         end do
      end do
      call prop_oneave(mol, S, (/'GEO'/), (/Def/), shape(Egff), Egff, DFD=(/DFDef/))
      DFDef = 0
      ! call print_tensor(Egff, shape(Egff), 'E1gDef+E2gDeDf-SgDFDef'); Egff=0
      call print_tensor(shape(Egff), -Egff, 'd/dg Alpha = -Egff', &
                        (/0*freq,-freq,freq/))
      DFDef=0
      !------------------------------------------------------------------
      ! electric-electric equations Df(w)f(w), solve 6 instead of 3x3 eqs.
      do l = 1, 3;
         call pert_dens(mol, S, (/'EL','EL'/), (/l,1/), &
                        (/D,Df(:l),Df(l)/), (/F,Ff(:l),Ff(l)/), &
                        Dff(:l,l), Fff(:l,l), freq=(/freq,freq/), comp=(/1,l/))
         do j = 1, l-1
            Dff(l,j) = Dff(j,l)
            Fff(l,j) = Fff(j,l)
         end do
      end do
      !---------------------------
      ! second hyperpolarizability
      Effff = 0 !no energy (oneave/twoave) contributions
      ! gradient multiplier contribution -tr (DeSD-) (FDS-)fef
      do i = 1, 3
         DeSD(i) = De(i)*S*D - D*S*De(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               FDSfef = (Ff(j)*Def(k,l) + Fe(k)*Dff(j,l) + Ff(l)*Def(k,j)      &
                      +  Fef(k,l)*Df(j) + Fff(j,l)*De(k) + Fef(k,j)*Df(l)) * S &
                  - S * (Df(j)*Fef(k,l) + De(k)*Fff(j,l) + Df(l)*Fef(k,j)      &
                      +  Def(k,l)*Ff(j) + Dff(j,l)*Fe(k) + Def(k,j)*Ff(l))
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(DeSD(i),FDSfef)
               end do; FDSfef=0
            end do
         end do
      end do
      DeSD = 0
      ! idempotency multiplier contribution -tr (FeDS+) (DSD)fef
      do i = 1, 3
         FeDS(i) = Fe(i)*D*S + S*D*Fe(i) - Fe(i)
      end do
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               DSDfef = Df(j)*S*Def(k,l) + De(k)*S*Dff(j,l) + Df(l)*S*Def(k,j) &
                      + Def(k,l)*S*Df(j) + Dff(j,l)*S*De(k) + Def(k,j)*S*Df(l)
               do i = 1, 3
                  Effff(i,j,k,l) = Effff(i,j,k,l) - tr(FeDS(i),DSDfef)
               end do; DSDfef=0
            end do
         end do
      end do
      FeDS = 0
      call print_tensor((/3,3,3,3/), -Effff, 'Gamma = -Effff', (/-freq,freq,-freq,freq/))
      ! delete matrices, reverse signs (from quasi-energy- to property-derivatives)
      De=0; Fe=0; Df=0; Ff=0; Def=0; Fef=0; Dff=0; Fff=0
      Ef=-Ef; Eff=-Eff; Effff=-Effff; Egff=-Egff
   end subroutine


   subroutine vibhyp_hyp_dipgra_polgra(mol, S, D, F, ng, freq, &
                                       dip, Eff, Efff, Egf, Egff)
      !> reference to molecule, solver- and integral settings
      type(prop_molcfg), intent(in)  :: mol
      !> overlap matrix
      type(matrix),      intent(in)  :: S
      !> density matrix
      type(matrix),      intent(in)  :: D
      !> Fock matrix
      type(matrix),      intent(in)  :: F
      !> number of geometrical coordinates
      integer,           intent(in)  :: ng
      !> external laser frequency
      complex(8),        intent(in)  :: freq(3)
      !> dipole moment returned in dip, polarizability in in Eff
      !> hyperpolarizability in Efff, dipole gradient in Egf,
      !> polarizability gradient in Egff.
      !> ajt FIXME: Sign-name confusion Ef should be -dipmom
      complex(8),        intent(out) :: dip(3), Eff(3,3,3), Efff(3,3,3), &
                                        Egf(ng,3,3), Egff(ng,3,3,3)
      type(matrix) :: Df(3,3), Dff(3,3,3), DFDf(3), &
                      Ff(3,3), Fff(3,3,3), DFDff(3,3)
      integer      :: i, j, k, m, n, nj, nk, mj, mk
      logical      :: exists
      ! verify that frequencies sum to zero
      if (abs(sum(freq)) > 1d-15) &
         call quit('vibhyp_hyp_dipgra_polgra: sum(freq) should be zero!',-1)
      ! verify that DALTON.HES exists before doing anything
#ifdef LSDALTON_ONLY
      call quit('Cannot call GPINQ, only new integral code is compiled',-1)
#else
      call GPINQ('DALTON.HES', 'EXIST', exists)
#endif
      if (.not.exists) call quit('vibhyp_hyp_dipgra_polgra: Hessian file', &
                                 ' DALTON.HES not found, but will be needed',-1)
      ! dipole moment
      dip = 0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), dip)
      dip = -dip
      call print_tensor((/3/), dip, 'dipmom = -Ef')
      !----------------------
      ! first order equations
      do n = 1, 3
         ! if frequency previously solved for, copy and skip
         do m = 1, n
            if (freq(m)==freq(n) .or. freq(m)==-freq(n)) exit
         end do
         if (m /= n) then
            do i=1,3
               if (freq(m)==freq(n)) Df(i,n) = Df(i,m)
               if (freq(m)==freq(n)) Ff(i,n) = Ff(i,m)
               if (freq(m)/=freq(n)) Df(i,n) = trps(Df(i,m))
               if (freq(m)/=freq(n)) Ff(i,n) = trps(Ff(i,m))
            end do
            cycle
         end if
         ! new frequency, so solve
         call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), &
                        Df(:,n), Ff(:,n), freq=(/freq(n)/))
      end do
      !-------------------
      ! 3 polarizabilities
      do n = 1, 3
         ! if frequency already calculated, copy and skip
         do m = 1, n
            if (freq(m)==freq(n) .or. freq(m)==-freq(n)) exit
         end do
         if (m /= n) then
            Eff(:,:,n) = Eff(:,:,m)
            cycle
         end if
         ! not previously calculated, so contract
         Eff(:,:,n) = 0
         call prop_oneave(mol, S, (/'EL'/), (/Df(:,n)/), (/3,3/), Eff(:,:,n))
         call print_tensor((/3,3/), -Eff(:,:,n), 'Alpha = -Eff', (/-freq(n),freq(n)/) )
      enddo
      !--------------------------
      ! 3 dipole moment gradients
      do n = 1, 3
         ! if frequency already calculated, copy and skip
         do m = 1, n
            if (freq(m)==freq(n) .or. freq(m)==-freq(n)) exit
         end do
         if (m /= n) then
            Egf(:,:,n) = Egf(:,:,m)
            cycle
         end if
         ! not previously calculated, so contract
         Egf(:,:,n) = 0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/D/), (/ng,3/), Egf(:,:,n))
         ! call print_tensor((/ng,3/), Egf(:,:,n), 'hgf + HgfD0'); Egf(:,:,n)=0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,n)/), (/ng,3/), Egf(:,:,n))
         ! call print_tensor((/ng,3/), Egf(:,:,n), 'GgD0Df'); Egf(:,:,n)=0
         do i = 1, 3
            DFDf(i) = Df(i,n)*(F+freq(n)/2*S)*D + D*Ff(i,n)*D &
                          + D*(F-freq(n)/2*S)*Df(i,n)
         end do
         call prop_oneave(mol, S, (/'GEO'/), (/Df(:,n)/), (/ng,3/), Egf(:,:,n), &
                          DFD=(/DFDf/), freq=(/-freq(n)/))
         DFDf = 0
         ! call print_tensor((/ng,3/), Egf(:,:,n), 'HgDf - i/2TgDf - Sg(DFD)f'); Egf(:,:,n)=0
         call print_tensor((/ng,3/), -Egf(:,:,n), 'd/dg dipmom = -Egf', &
                           (/-freq(n),freq(n)/))
      end do
      !-----------------------
      ! second order equations
      do n = 1, 3
         nj = max(1,n-1) !first electric frequency
         nk = min(3,n+1) !second electric frequency
         ! if these equations have been solved previously, copy and skip
         do m = 1, n
            mj = max(1,m-1)
            mk = min(3,m+1)
            if (freq(mj)== freq(nj) .and. freq(mk)== freq(nk) .or. &
                freq(mj)== freq(nk) .and. freq(mk)== freq(nj) .or. &
                freq(mj)==-freq(nj) .and. freq(mk)==-freq(nk) .or. &
                freq(mj)==-freq(nk) .and. freq(mk)==-freq(nj)) exit
         end do
         if (m /= n) then
            do k = 1, 3
               do j = 1, 3
                  if (freq(mj)==freq(nj)) then
                     Dff(j,k,n) = Dff(j,k,m)
                     Fff(j,k,n) = Fff(j,k,m)
                  else if (freq(mj)==freq(nk)) then
                     Dff(j,k,n) = Dff(k,j,m)
                     Fff(j,k,n) = Fff(k,j,m)
                  else if (freq(mj)==-freq(nj)) then
                     Dff(j,k,n) = trps(Dff(j,k,m))
                     Fff(j,k,n) = trps(Fff(j,k,m))
                  else !freq(mj)==-freq(nk)
                     Dff(j,k,n) = trps(Dff(k,j,m))
                     Fff(j,k,n) = trps(Fff(k,j,m))
                  end if
               end do
            end do
            cycle
         end if
         ! new pair of frequencies, so must solve
         if (freq(nj)==freq(nk) .or. freq(nj)==-freq(nk)) then
            do k = 1, 3
               call pert_dens(mol, S, (/'EL','EL'/), (/k,1/),  &
                              (/D,Df(1:k,nj),Df(k,nk)/),  &
                              (/F,Ff(1:k,nj),Ff(k,nk)/),  &
                              Dff(1:k,k,n), Fff(1:k,k,n), &
                              freq=(/freq(nj),freq(nk)/), comp=(/1,k/))
               do j = 1, k-1
                  if (freq(nj)==freq(nk)) Dff(k,j,n) = Dff(j,k,n)
                  if (freq(nj)==freq(nk)) Fff(k,j,n) = Fff(j,k,n)
                  if (freq(nj)/=freq(nk)) Dff(k,j,n) = trps(Dff(j,k,n))
                  if (freq(nj)/=freq(nk)) Fff(k,j,n) = trps(Fff(j,k,n))
               end do
            end do
         else
            call pert_dens(mol, S, (/'EL','EL'/), (/3,3/),                        &
                           (/D,Df(:,nj),Df(:,nk)/), (/F,Ff(:,nj),Ff(:,nk)/), &
                           Dff(:,:,n), Fff(:,:,n), freq=(/freq(nj),freq(nk)/))
         end if
      end do
      !-----------------------------
      ! 1 first hyperpolarizabilitiy
      Efff = 0
      call prop_oneave(mol, S, (/'EL'/), (/Dff(:,:,3)/), (/3,3,3/), Efff)
      call print_tensor((/3,3,3/), -Efff, 'Beta = -Efff', &
                        (/freq(1),freq(2),freq(3)/))
      !---------------------------
      ! 3 polarizability gradients
      do n = 1, 3
         nj = max(1,n-1) !first electric frequency
         nk = min(3,n+1) !second electric frequency
         ! if already calculated for these frequencies, copy and skip
         do m = 1, n
            mj = max(1,m-1)
            mk = min(3,m+1)
            if (freq(mj)== freq(nj) .and. freq(mk)== freq(nk) .or. &
                freq(mj)== freq(nk) .and. freq(mk)== freq(nj) .or. &
                freq(mj)==-freq(nj) .and. freq(mk)==-freq(nk) .or. &
                freq(mj)==-freq(nk) .and. freq(mk)==-freq(nj)) exit
         end do
         if (m /= n) then
            if (freq(mj)==freq(nj) .or. freq(mj)==-freq(nj)) then
               Egff(:,:,:,n) = Egff(:,:,:,m)
            else
               do k = 1, 3
                  do j = 1, 3
                  Egff(:,j,k,n) = Egff(:,k,j,m)
               end do
               end do
            end if
            cycle
         end if
         ! not calculated before, so calculate
         Egff(:,:,:,n) = 0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Df(:,nk)/), (/ng,3,3/), Egff(:,:,:,n))
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), 'E1geDf'); Egff(:,:,:,n)=0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Df(:,nj)/), (/ng,3,3/), Egff(:,:,:,n), perm=(/1,3,2/))
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), 'E1gfDe'); Egff(:,:,:,n)=0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,nj),Df(:,nk),Dff(:,:,n)/), &
                          (/ng,3,3/), Egff(:,:,:,n))
         do k = 1, 3
            do j = 1, 3
               DFDff(j,k) = Df(j,nj)*Ff(k,nk)*D + D*Ff(k,nk)*Df(j,nj)     &
                          + Df(k,nk)*Ff(j,nj)*D + D*Ff(j,nj)*Df(k,nk)     &
                          + Df(k,nk)*(F+(freq(nk)-freq(nj))/2*S)*Df(j,nj) &
                          + Df(j,nj)*(F+(freq(nj)-freq(nk))/2*S)*Df(k,nk) &
                          + D*(F-(freq(nj)+freq(nk))/2*S)*Dff(j,k,n)      &
                          + Dff(j,k,n)*(F+(freq(nj)+freq(nk))/2*S)*D + D*Fff(j,k,n)*D
            end do
         end do
         call prop_oneave(mol, S, (/'GEO'/), (/Dff(:,:,n)/), (/ng,3,3/), Egff(:,:,:,n), &
                          DFD=(/DFDff/), freq=(/-freq(nj)-freq(nk)/))
         DFDff = 0
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), &
         !                   'E1gDef + E2gDeDf - i/2TgDef - Sg(DFD)ef'); Egff(:,:,:,n)=0
         call print_tensor((/ng,3,3/), -Egff(:,:,:,n), 'd/dg Alpha = -Egff', &
                           (/-freq(nj)-freq(nk),freq(nj),freq(nk)/))
      end do
      ! delete/free perturbed matrices
      Df=0; Dff=0; Ff=0; Fff=0
      ! reverse signs on quasi-energy derivatives to get dipole, alpha and beta
      Eff=-Eff; Efff=-Efff; Egf=-Egf; Egff=-Egff
   end subroutine



   subroutine vibshyp_shyp_dipg_polg_hypg(mol, S, D, F, ng, freq, &
                                          dip, Effff, Egf, Egff, Egfff)
      !> reference to molecule, solver- and integral settings
      type(prop_molcfg), intent(in)  :: mol
      !> overlap matrix
      type(matrix),      intent(in)  :: S
      !> density matrix
      type(matrix),      intent(in)  :: D
      !> Fock matrix
      type(matrix),      intent(in)  :: F
      !> number of geometrical coordinates
      integer,           intent(in)  :: ng
      !> external laser frequencies
      complex(8),        intent(in)  :: freq(4)
      !> dipole moment returned in dip, polarizability in in Eff
      !> 2nd hyperpolarizability in Effff, 4 dipole gradients in Egf,
      !> 6 polarizability gradients in Egff. and 4 hyperpolarizability
      !> gradients in Egfff
      !> ajt FIXME: Sign-name confusion Ef should be -dipmom
      complex(8),        intent(out) :: dip(3), Effff(3,3,3,3), Egf(ng,3,4)
      complex(8),        intent(out) :: Egff(ng,3,3,6), Egfff(ng,3,3,3,4)
      type(matrix) :: Df(3,4), Dff(3,3,6), Dfff(3,3,3), DFDf(3), DFDff(3,3)
      type(matrix) :: Ff(3,4), Fff(3,3,6), Ffff(3,3,3), DFDfff(3,3,3), DFD
      type(matrix) :: FD0, FD1, FD2, FD3
      complex(8)   :: Eff(3,3), Efff(3,3,3), Rg333(ng,3,3,3)
      integer      :: i, j, k, l, m, n, ni, nj, nk, nl, mi, mj, mk, ml
      integer      :: nij, nik, njk, njl, nkl, mij, mik, mjk
      logical      :: exists
      ! verify that frequencies sum to zero
      if (abs(sum(freq)) > 1d-15) &
         call quit('vibgam_shyp_dipg_polg_hypg: sum(freq) should be zero!',-1)
      ! verify that DALTON.HES exists before doing anything
#ifdef LSDALTON_ONLY
      call quit('Cannot call GPINQ, only new integral code is compiled',-1)
#else
      call GPINQ('DALTON.HES', 'EXIST', exists)
#endif
      if (.not.exists) call quit('vibgam_shyp_dipg_polg_hypg: Hessian file', &
                                 ' DALTON.HES not found, but will be needed',-1)
      ! dipole moment
      dip = 0
      call prop_oneave(mol, S, (/'EL'/), (/D/), (/3/), dip)
      dip = -dip
      call print_tensor((/3/), dip, 'dipmom = -Ef')
      !----------------------
      ! first order equations
      do n=1,4
         ! if frequency previously solved for, copy and skip
         do m=1,n
            if (freq(m)==freq(n) .or. freq(m)==-freq(n)) exit
         end do
         if (m /= n) then
            do i=1,3
               if (freq(m)==freq(n)) Df(i,n) = Df(i,m)
               if (freq(m)==freq(n)) Ff(i,n) = Ff(i,m)
               if (freq(m)/=freq(n)) Df(i,n) = trps(Df(i,m))
               if (freq(m)/=freq(n)) Ff(i,n) = trps(Ff(i,m))
               Egf(:,:,n) = Egf(:,:,m)
            end do
            cycle
         end if
         ! new frequency, so solve
         call pert_dens(mol, S, (/'EL'/), (/3/), (/D/), (/F/), &
                        Df(:,n), Ff(:,n), freq=(/freq(n)/))
         ! polarizability
         Eff(:,:) = 0
         call prop_oneave(mol, S, (/'EL'/), (/Df(:,n)/), (/3,3/), Eff(:,:))
         call print_tensor((/3,3/), -Eff(:,:), 'Alpha = -Eff', (/-freq(n),freq(n)/))
         ! dipole moment gradient
         Egf(:,:,n) = 0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/D/), (/ng,3/), Egf(:,:,n))
         ! call print_tensor((/ng,3/), Egf(:,:,n), 'E0gf'); Egf(:,:,n)=0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,n)/), (/ng,3/), Egf(:,:,n))
         ! call print_tensor( (/3,ng/), Efg(:,:,n), 'E1gDf'); Efg(:,:,n)=0
         do i=1,3
            DFDf(i) = Df(i,n)*(F+freq(n)/2*S)*D + D*Ff(i,n)*D &
                          + D*(F-freq(n)/2*S)*Df(i,n)
         end do
         call prop_oneave(mol, S, (/'GEO'/), (/Df(:,n)/), (/ng,3/), Egf(:,:,n), &
                          DFD=(/DFDf/), freq=(/-freq(n)/))
         DFDf = 0
         ! call print_tensor((/ng,3/), Egf(:,:,n), 'E1gDf - i/2TgDf - Sg(DFD)f'); Egf(:,:,n)=0
         call print_tensor((/ng,3/), -Egf(:,:,n), 'd/dg dipmom = -Egf', &
                           (/-freq(n),freq(n)/))
      end do
      !-----------------------
      ! second order equations
      do n=1,6
         nj = merge(merge(2,3,n<2),4,n<4)
         ni = n - (nj-2)*(nj-1)/2
         ! if these equations have been solved previously, copy and skip
         do m=1,n
            mj = merge(merge(2,3,m<2),4,m<4)
            mi = m - (mj-2)*(mj-1)/2
            if ((freq(mi)== freq(ni) .and. freq(mj)== freq(nj)) .or. &
                (freq(mi)== freq(nj) .and. freq(mj)== freq(ni)) .or. &
                (freq(mi)==-freq(ni) .and. freq(mj)==-freq(nj)) .or. &
                (freq(mi)==-freq(nj) .and. freq(mj)==-freq(ni))) exit
         end do
         if (m/=n) then
            do j=1,3
               do i=1,3
                  if (freq(mi)==freq(ni) .and. freq(mj)== freq(nj)) then
                     Dff(i,j,n) = Dff(i,j,m)
                     Fff(i,j,n) = Fff(i,j,m)
                     Egff(:,i,j,n) = Egff(:,i,j,m)
                  else if (freq(mi)==freq(nj) .and. freq(mj)== freq(ni)) then
                     Dff(i,j,n) = Dff(j,i,m)
                     Fff(i,j,n) = Fff(j,i,m)
                     Egff(:,i,j,n) = Egff(:,j,i,m)
                  else if (freq(mi)==-freq(ni) .and. freq(mj)==-freq(nj)) then
                     Dff(i,j,n) = trps(Dff(i,j,m))
                     Fff(i,j,n) = trps(Fff(i,j,m))
                     Egff(:,i,j,n) = Egff(:,i,j,m)
                  else !freq(mi)==-freq(nj) .and. freq(mj)==-freq(ni)
                     Dff(i,j,n) = trps(Dff(j,i,m))
                     Fff(i,j,n) = trps(Fff(j,i,m))
                     Egff(:,i,j,n) = Egff(:,j,i,m)
                  end if
               end do
            end do
            cycle
         end if
         ! new pair of frequencies, so must solve
         ! if same or opposite freqs, only 6 eqns
         if (freq(ni)==freq(nj) .or. freq(ni)==-freq(nj)) then
            do j=1,3
               call pert_dens(mol, S, (/'EL','EL'/), (/j,1/), &
                              (/D,Df(1:j,ni),Df(j,nj)/),      &
                              (/F,Ff(1:j,ni),Ff(j,nj)/),      &
                              Dff(1:j,j,n), Fff(1:j,j,n),     &
                              freq=(/freq(ni),freq(nj)/), comp=(/1,j/))
               do i=1,j-1
                  if (freq(ni)==freq(nj)) Dff(j,i,n) = Dff(i,j,n)
                  if (freq(ni)==freq(nj)) Fff(j,i,n) = Fff(i,j,n)
                  if (freq(ni)/=freq(nj)) Dff(j,i,n) = trps(Dff(i,j,n))
                  if (freq(ni)/=freq(nj)) Fff(j,i,n) = trps(Fff(i,j,n))
               end do
            end do
         else
            call pert_dens(mol, S, (/'EL','EL'/), (/3,3/),                   &
                           (/D,Df(:,ni),Df(:,nj)/), (/F,Ff(:,ni),Ff(:,nj)/), &
                           Dff(:,:,n), Fff(:,:,n), freq=(/freq(ni),freq(nj)/))
         end if
         ! 1st hyperpolarizability
         Efff(:,:,:) = 0
         call prop_oneave(mol, S, (/'EL'/), (/Dff(:,:,n)/), (/3,3,3/), Efff)
         call print_tensor((/3,3,3/), -Efff, 'Beta = -Efff', &
                           (/-freq(ni)-freq(nj),freq(ni),freq(nj)/))
         ! polarizability gradient
         Egff(:,:,:,n) = 0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Df(:,nj)/), (/ng,3,3/), Egff(:,:,:,n))
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), 'E1geDf'); Egff(:,:,:,n)=0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Df(:,ni)/), (/ng,3,3/), &
                          Egff(:,:,:,n), perm=(/1,3,2/))
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), 'E1gfDe'); Egff(:,:,:,n)=0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,ni),Df(:,nj),Dff(:,:,n)/), &
                          (/ng,3,3/), Egff(:,:,:,n))
         ! prepare DFDff
         FD0 = (F + (freq(ni)+freq(nj))/2 * S) * D
         do j=1,3
            FD1 = Ff(j,nj)*D + (F + (freq(ni)-freq(nj))/2 * S) * Df(j,nj)
            do i=1,3
               FD2 = Fff(i,j,n)*D + Ff(j,nj)*Df(i,ni) + Ff(i,ni)*Df(j,nj) &
                   + (F - (freq(ni)+freq(nj))/2 * S) * Dff(i,j,n)
               DFDff(i,j) = D*FD2 + Df(i,ni)*FD1 + Dff(i,j,n)*FD0
            end do
         end do
         FD0=0; FD2=0
         do i=1,3
            FD1 = Ff(i,ni)*D + (F - (freq(ni)-freq(nj))/2 * S) * Df(i,ni)
            do j=1,3
               DFDff(i,j) = DFDff(i,j) + Df(j,nj)*FD1
            end do
         end do
         FD1=0; FD2=0
         call prop_oneave(mol, S, (/'GEO'/), (/Dff(:,:,n)/), (/ng,3,3/), Egff(:,:,:,n), &
                          DFD=(/DFDff/), freq=(/-freq(ni)-freq(nj)/))
         DFDff(:,:) = 0
         ! call print_tensor((/ng,3,3/), Egff(:,:,:,n), &
         !                   'E1gDef + E2gDeDf - i/2TgDef - Sg(DFD)ef'); Egff(:,:,:,n)=0
         call print_tensor((/ng,3,3/), -Egff(:,:,:,n), 'd/dg Alpha = -Egff', &
                           (/-freq(ni)-freq(nj),freq(ni),freq(nj)/))
      end do
      !------------------------
      ! 4 third order equations
      do n=1,4
         ni = merge(1,2,n<4)
         nj = merge(2,3,n<3)
         nk = merge(3,4,n<2)
         nij = ni + (nj-2)*(nj-1)/2
         nik = ni + (nk-2)*(nk-1)/2
         njk = nj + (nk-2)*(nk-1)/2
         !ajt This is overkill, but I'm in a hurry
         do m=1,n
            mi = merge(1,2,m<4)
            mj = merge(2,3,m<3)
            mk = merge(3,4,m<2)
            mij = mi + (mj-2)*(mj-1)/2
            mik = mi + (mk-2)*(mk-1)/2
            mjk = mj + (mk-2)*(mk-1)/2
            if ((freq(mi)== freq(ni) .and. freq(mj)== freq(nj)       &
                                     .and. freq(mk)== freq(nk)) .or. &
                (freq(mi)== freq(nj) .and. freq(mj)== freq(ni)       &
                                     .and. freq(mk)== freq(nk)) .or. &
                (freq(mi)== freq(nk) .and. freq(mj)== freq(nj)       &
                                     .and. freq(mk)== freq(ni)) .or. &
                (freq(mi)== freq(ni) .and. freq(mj)== freq(nk)       &
                                     .and. freq(mk)== freq(nj)) .or. &
                (freq(mi)== freq(nj) .and. freq(mj)== freq(nk)       &
                                     .and. freq(mk)== freq(ni)) .or. &
                (freq(mi)== freq(nk) .and. freq(mj)== freq(ni)       &
                                     .and. freq(mk)== freq(nj)) .or. &
                (freq(mi)==-freq(ni) .and. freq(mj)==-freq(nj)       &
                                     .and. freq(mk)==-freq(nk)) .or. &
                (freq(mi)==-freq(nj) .and. freq(mj)==-freq(ni)       &
                                     .and. freq(mk)==-freq(nk)) .or. &
                (freq(mi)==-freq(nk) .and. freq(mj)==-freq(nj)       &
                                     .and. freq(mk)==-freq(ni)) .or. &
                (freq(mi)==-freq(ni) .and. freq(mj)==-freq(nk)       &
                                     .and. freq(mk)==-freq(nj)) .or. &
                (freq(mi)==-freq(nj) .and. freq(mj)==-freq(nk)       &
                                     .and. freq(mk)==-freq(ni)) .or. &
                (freq(mi)==-freq(nk) .and. freq(mj)==-freq(ni)       &
                                     .and. freq(mk)==-freq(nj))) exit
         end do
         if (m/=n) then
            do k=1,3
               do j=1,3
                  do i=1,3
                     if ((freq(mi)== freq(ni) .and. freq(mj)== freq(nj)            &
                              .and. freq(mk)== freq(nk)) .or. (freq(mi)==-freq(ni) &
                              .and. freq(mj)==-freq(nj) .and. freq(mk)==-freq(nk))) then
                        Egfff(:,i,j,k,n) = Egfff(:,i,j,k,m)
                     else if ((freq(mi)== freq(nj) .and. freq(mj)== freq(ni)       &
                              .and. freq(mk)== freq(nk)) .or. (freq(mi)==-freq(nj) &
                              .and. freq(mj)==-freq(ni) .and. freq(mk)==-freq(nk))) then
                        Egfff(:,i,j,k,n) = Egfff(:,j,i,k,m)
                     else if ((freq(mi)== freq(nk) .and. freq(mj)== freq(nj)       &
                              .and. freq(mk)== freq(ni)) .or. (freq(mi)==-freq(nk) &
                              .and. freq(mj)==-freq(nj) .and. freq(mk)==-freq(ni))) then
                        Egfff(:,i,j,k,n) = Egfff(:,k,j,i,m)
                     else if ((freq(mi)== freq(ni) .and. freq(mj)== freq(nk)       &
                              .and. freq(mk)== freq(nj)) .or. (freq(mi)==-freq(ni) &
                              .and. freq(mj)==-freq(nk) .and. freq(mk)==-freq(nj))) then
                        Egfff(:,i,j,k,n) = Egfff(:,i,k,j,m)
                     else if ((freq(mi)== freq(nj) .and. freq(mj)== freq(nk)       &
                              .and. freq(mk)== freq(ni)) .or. (freq(mi)==-freq(nj) &
                              .and. freq(mj)==-freq(nk) .and. freq(mk)==-freq(ni))) then
                        Egfff(:,i,j,k,n) = Egfff(:,k,i,j,m)
                     else !if ((freq(mi)== freq(nk) .and. freq(mj)== freq(ni)       &
                          !    .and. freq(mk)== freq(nj)) .or. (freq(mi)==-freq(nk) &
                          !    .and. freq(mj)==-freq(ni) .and. freq(mk)==-freq(nj))) then
                        Egfff(:,i,j,k,n) = Egfff(:,j,k,i,m)
                     end if
                  end do
               end do
            end do
            cycle
         end if
         ! new triple of frequencies, so must solve
         if (freq(ni)==freq(nj) .and. freq(ni)==freq(nk)) then
            ! all freqs equal, just 10 eqns
            do k=1,3
               do j=1,k
                  call pert_dens(mol, S, (/'EL','EL','EL'/), (/j,1,1/),        &
                                 (/D,Df(:j,ni),Df(j,nj),Df(k,nk),              &
                                   Dff(:j,j,nij),Dff(:j,k,nik),Dff(j,k,njk)/), &
                                 (/F,Ff(:j,ni),Ff(j,nj),Ff(k,nk),              &
                                   Fff(:j,j,nij),Fff(:j,k,nik),Fff(j,k,njk)/), &
                                 Dfff(:j,j,k), Ffff(:j,j,k),                   &
                                 freq=(/freq(ni),freq(nj),freq(nk)/), comp=(/1,j,k/))
                  do i=1,j-1 !i-j symmetry
                     Dfff(j,i,k) = Dfff(i,j,k)
                     Ffff(j,i,k) = Ffff(i,j,k)
                  end do
               end do
               do i=1,k-1 !i-k symmetry
                  do j=1,k-1
                     Dfff(k,j,i) = Dfff(i,j,k)
                     Ffff(k,j,i) = Ffff(i,j,k)
                  end do
               end do
               do j=1,k-1 !j-k symmetry
                  do i=1,k
                     Dfff(i,k,j) = Dfff(i,j,k)
                     Ffff(i,k,j) = Ffff(i,j,k)
                  end do
               end do
            end do
         else if (freq(nj)==freq(nk) .or. (freq(ni)==0 .and. freq(nj)==-freq(nk))) then
            ! j-k symmetry, 20 eqns
            do k=1,3
               call pert_dens(mol, S, (/'EL','EL','EL'/), (/3,k,1/),        &
                              (/D,Df(:,ni),Df(:k,nj),Df(k,nk),              &
                                Dff(:,:k,nij),Dff(:,k,nik),Dff(:k,k,njk)/), &
                              (/F,Ff(:,ni),Ff(:k,nj),Ff(k,nk),              &
                                Fff(:,:k,nij),Fff(:,k,nik),Fff(:k,k,njk)/), &
                              Dfff(:,:k,k), Ffff(:,:k,k),                   &
                              freq=(/freq(ni),freq(nj),freq(nk)/), comp=(/1,1,k/))
               do j=1,k-1
                  do i=1,3
                     if (freq(nj)==freq(nk)) Dfff(i,k,j) = Dfff(i,j,k)
                     if (freq(nj)==freq(nk)) Ffff(i,k,j) = Ffff(i,j,k)
                     if (freq(nj)/=freq(nk)) Dfff(i,k,j) = trps(Dfff(i,j,k))
                     if (freq(nj)/=freq(nk)) Ffff(i,k,j) = trps(Ffff(i,j,k))
                  end do
               end do
            end do
         else
            ! different frequencies, 27 eqns
            call pert_dens(mol, S, (/'EL','EL','EL'/), (/3,3,3/),      &
                           (/D,Df(:,ni),Df(:,nj),Df(:,nk),             &
                             Dff(:,:,nij),Dff(:,:,nik),Dff(:,:,njk)/), &
                           (/F,Ff(:,ni),Ff(:,nj),Ff(:,nk),             &
                             Fff(:,:,nij),Fff(:,:,nik),Fff(:,:,njk)/), &
                           Dfff(:,:,:), Ffff(:,:,:),                   &
                           freq=(/freq(ni),freq(nj),freq(nk)/))
         end if
         ! 2nd hyperpolarizability
         Effff(:,:,:,:) = 0
         call prop_oneave(mol, S, (/'EL'/), (/Dfff/), (/3,3,3,3/), Effff)
         call print_tensor((/3,3,3,3/), -Effff, 'Gamma = -Effff', &
                           (/-freq(ni)-freq(nj)-freq(nk),freq(ni),freq(nj),freq(nk)/))
         ! 1st hyperpolarizability gradient
         Egfff(:,:,:,:,n) = 0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Dff(:,:,njk)/), &
                          (/ng,3,3,3/), Egfff(:,:,:,:,n))
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E1gdDef'); Egfff(:,:,:,:,n)=0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Dff(:,:,nik)/), &
                          (/ng,3,3,3/), Egfff(:,:,:,:,n), perm=(/1,3,2,4/))
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E1geDdf'); Egfff(:,:,:,:,n)=0
         call prop_oneave(mol, S, (/'GEO','EL '/), (/Dff(:,:,nij)/), &
                          (/ng,3,3,3/), Egfff(:,:,:,:,n), perm=(/1,4,2,3/))
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E1gfDde'); Egfff(:,:,:,:,n)=0
         ! --- two-electron contribution
         ! third-order density is currently not supported in prop_twoave,
         ! but can be 'faked' in terms of three second-order densities (*HF ONLY*)
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,ni),Dff(:,:,njk), &
                          Dfff(:,:,:)/), (/ng,3,3*3/), Egfff(:,:,:,:,n))
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E2GD0Ddef + E2gDdDef'); Egfff(:,:,:,:,n)=0
         Rg333 = 0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,nj),Dff(:,:,nik), &
                          (mol%zeromat,i=1,3*3*3)/), (/ng,3,3*3/), Rg333)
         do j=1,3
            do i=1,3
               Egfff(:,i,j,:,n) = Egfff(:,i,j,:,n) + Rg333(:,j,i,:)
            end do
         end do
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E2gDeDdf'); Egfff(:,:,:,:,n)=0
         call prop_twoave(mol, (/'GEO'/), (/D,Df(:,nk),Dff(:,:,nij), &
                          (mol%zeromat,i=1,3*3*3)/), (/ng,3*3,3/),   &
                          Egfff(:,:,:,:,n), perm=(/1,3,2/))
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), 'E2gDfDde'); Egfff(:,:,:,:,n)=0
         ! --- one-electron contribution
         ! reorth matrices DFDfff, to be contracted with perturbed overlap Sg
         ! first loop: k, j, i
         FD0 = (F + (freq(ni)+freq(nj)+freq(nk))/2 * S) * D
         do k=1,3
            FD1 = Ff(k,nk)*D + (F + (freq(ni)+freq(nj)-freq(nk))/2 * S) * Df(k,nk)
            do j=1,3
               FD2 = Fff(j,k,njk)*D + Ff(k,nk)*Df(j,nj) + Ff(j,nj)*Df(k,nk) &
                   + (F + (freq(ni)-freq(nj)-freq(nk))/2 * S) * Dff(j,k,njk)
               do i=1,3
                  FD3 = Ffff(i,j,k)*D + (F - (freq(ni)+freq(nj)+freq(nk))/2 * S) * Dfff(i,j,k) &
                      + Fff(j,k,njk)*Df(i,ni) + Fff(i,k,nik)*Df(j,nj) + Fff(i,j,nij)*Df(k,nk)  &
                      + Ff(i,ni)*Dff(j,k,njk) + Ff(j,nj)*Dff(i,k,nik) + Ff(k,nk)*Dff(i,j,nij)
                  DFDfff(i,j,k) = D*FD3 + Df(i,ni)*FD2 + Dff(i,j,nij)*FD1 + Dfff(i,j,k)*FD0
               end do
            end do
         end do
         FD0=0; FD3=0
         ! second loop; j, i, k
         do j=1,3
            FD1 = Ff(j,nj)*D + (F + (freq(ni)-freq(nj)+freq(nk))/2 * S) * Df(j,nj)
            do i=1,3
               FD2 = Fff(i,j,nij)*D + Ff(j,nj)*Df(i,ni) + Ff(i,ni)*Df(j,nj) &
                   + (F - (freq(ni)+freq(nj)-freq(nk))/2 * S) * Dff(i,j,nij)
               do k=1,3
                  DFDfff(i,j,k) = DFDfff(i,j,k) + Df(k,nk)*FD2
                  DFDfff(i,j,k) = DFDfff(i,j,k) + Dff(i,k,nik)*FD1
               end do
            end do
         end do
         ! last loop: i, k, j
         do i=1,3
            FD1 = Ff(i,ni)*D + (F - (freq(ni)-freq(nj)-freq(nk))/2 * S) * Df(i,ni)
            do k=1,3
               FD2 = Fff(i,k,nik)*D + Ff(k,nk)*Df(i,ni) + Ff(i,ni)*Df(k,nk) &
                   + (F - (freq(ni)-freq(nj)+freq(nk))/2 * S) * Dff(i,k,nik)
               do j=1,3
                  DFDfff(i,j,k) = DFDfff(i,j,k) + Df(j,nj)*FD2
                  DFDfff(i,j,k) = DFDfff(i,j,k) + Dff(j,k,njk)*FD1
               end do
            end do
         end do
         FD1=0; FD2=0
         ! contract with 1DHAM, SQHDOR and DEROVL integrals
         call prop_oneave(mol, S, (/'GEO'/), (/Dfff/), (/ng,3,3,3/), &
                          Egfff(:,:,:,:,n), DFD = (/DFDfff/),        &
                          freq = (/-freq(ni)-freq(nj)-freq(nk)/))
         DFDfff(:,:,:) = 0
         ! call print_tensor((/ng,3,3,3/), Egfff(:,:,:,:,n), &
         !                   'HgDdef - i/2TgDdef - Sg(DFD)def'); Egfff(:,:,:,:,n)=0
         call print_tensor((/ng,3,3,3/), -Egfff(:,:,:,:,n), 'd/dg Beta = -Egfff', &
                           (/-freq(ni)-freq(nj)-freq(nk),freq(ni),freq(nj),freq(nk)/))
      end do
   end subroutine


end module
