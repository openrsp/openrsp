! -------------------------------------------------------------------------
!
! Library:      InteRest
!
! File:         module_interest_osr.f90 
!
! Licensing:    Copyright 2010 Michal Repisky, and Kenneth Ruud
!
!               This file is part of InteRest.
!
!               InteRest is free software: you can redistribute it and/or modify
!               it under the terms of the GNU Lesse General Public License as published by
!               the Free Software Foundation, either version 3 of the License, or
!               (at your option) any later version.
!               
!               InteRest is distributed in the hope that it will be useful,
!               but WITHOUT ANY WARRANTY; without even the implied warranty of
!               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!               GNU Lesser General Public License for more details.
!               
!               You should have received a copy of the GNU Lesser General Public License
!               along with InteRest.  If not, see <http://www.gnu.org/licenses/>.
!
! Description:  Performing integral calculations over normalized cartesian/spherical
!               unbalanced, or kinetically/magnetically balanced gaussian type functions
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revision:     2.0   
!
! -------------------------------------------------------------------------

module module_interest_osr
  implicit none

  public  interest_osr_initialize
  public  interest_osr_class_overlap
  public  interest_osr_class_dipole
  public  interest_osr_class_nuclear

  integer, public, parameter :: lmx1 = 8 
  integer, public, parameter :: lmx2 = 2*lmx1-1 
  integer, public, parameter :: lmx4 = 4*lmx1-3 
  integer, public, parameter :: ncc(lmx1) = (/1,3,6,10,15,21,28,36/)
  integer, public, parameter :: nsc(lmx1) = (/1,3,5, 7, 9,11,13,15/)
  integer, public, parameter :: ncc2 = (lmx2*(lmx2+1)*(lmx2+2))/6 
  real(8), public, parameter :: pi   = 4.0d0*datan(1.0d0)

  logical, public, save :: is_interest_initialized = .false.


  private

  integer, parameter :: lmxx = lmx4+6
  real(8), parameter :: tf   = 92.0d0  !!dfloat(2*(4*(lmx1-1)+8)+36)

  integer :: mo(6,ncc2)                !fixme: dynamic allocation

  real(8) :: f(13)
  real(8) :: fnx(lmx4)
  real(8) :: gama(0:lmx2-1)
  real(8) :: fnxtab(121*lmxx)          !fixme: dynamic allocation
  real(8) :: vaux(ncc2,lmx4)           !fixme: dynamic allocation
  real(8) :: waux(ncc2,ncc2,lmx2)      !fixme: dynamic allocation + size???
  
  real(8) :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14
  real(8) :: c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27
  real(8) :: c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40
contains

! ---------------------------------------------------------------------------------------
! CLASS_OVERLAP
!
!note: na vstupe musi byt pre RKB bazu: (la-1), (lab+2)
!note: li1 (on output) - dlzka intermediates
!note: na vystupe, 1-elektronove integraly: [la...lab|0]
! ---------------------------------------------------------------------------------------
  subroutine interest_osr_class_overlap(li1,fint,cosr,la,lab,pexp,rxpa,rypa,rzpa)
    !-- output --!
    integer, intent(out) :: li1
    real(8), intent(out) :: fint(*)

    !-- input --!
    integer, intent(in) :: la,lab
    real(8), intent(in) :: cosr
    real(8), intent(in) :: pexp,rxpa,rypa,rzpa

    !-- local --!
    integer :: n,l,ll,lm,ln,lx,ly
    real(8) :: sint(3,lmx2)
    real(8) :: o2p,pexp2
    real(8) :: sxone,syone,szone
    real(8) :: sxtwo,sytwo,sztwo

    sxtwo = 0.0d0
    sytwo = 0.0d0
    sztwo = 0.0d0
    sint(1,1) = dsqrt(pi/pexp)
    sint(2,1) = sint(1,1)
    sint(3,1) = sint(2,1)
    pexp2 = 2.0d0*pexp
    do n=2,lab,1
      o2p   = dfloat(n-2)/pexp2
      sxone = sxtwo
      syone = sytwo
      szone = sztwo
      sxtwo = sint(1,n-1)
      sytwo = sint(2,n-1)
      sztwo = sint(3,n-1)
      sint(1,n) = rxpa*sxtwo + o2p*sxone
      sint(2,n) = rypa*sytwo + o2p*syone
      sint(3,n) = rzpa*sztwo + o2p*szone
    enddo

    li1=0
    n=max(1,la)
    do l=n,lab,1
      do lx=1,l
        ll=l-lx+1
        do ly=1,lx,1
          lm =lx-ly+1
          ln =ly
          li1=li1+1
          fint(li1)=cosr*sint(1,ll)*sint(2,lm)*sint(3,ln)
        enddo
      enddo
    enddo
  end subroutine 

! ---------------------------------------------------------------------------------------
! CLASS_DIPOLE
!
!note: na vstupe musi byt pre RKB bazu: (la-1), (lab+2)
!note: li1 (on output) - dlzka intermediates
!note: na vystupe, 1-elektronove integraly: [la...lab|x/y/z|0]
! ---------------------------------------------------------------------------------------
  subroutine interest_osr_class_dipole(li1,fintx,finty,fintz,cosr,la,lab,pexp,&
                                       rxpa,rypa,rzpa,rxpc,rypc,rzpc)
    !-- output --!
    integer, intent(out) :: li1
    real(8), intent(out) :: fintx(*)
    real(8), intent(out) :: finty(*)
    real(8), intent(out) :: fintz(*)

    !-- input --!
    integer, intent(in) :: la,lab
    real(8), intent(in) :: cosr,pexp
    real(8), intent(in) :: rxpa,rypa,rzpa
    real(8), intent(in) :: rxpc,rypc,rzpc

    !-- local --!
    integer :: n,l,ll,lm,ln,lx,ly
    real(8) :: sint(3,lmx2)
    real(8) :: dipi(3,lmx2)
    real(8) :: o2psi,o2pdi,pexp2
    real(8) :: sxone,syone,szone
    real(8) :: sxtwo,sytwo,sztwo


    o2pdi = 0.0d0
    sxtwo = 0.0d0
    sytwo = 0.0d0
    sztwo = 0.0d0
    sint(1,1) =      dsqrt(pi/pexp)
    sint(2,1) =      sint(1,1)
    sint(3,1) =      sint(2,1)
    dipi(1,1) = rxpc*sint(1,1)
    dipi(2,1) = rypc*sint(2,1)
    dipi(3,1) = rzpc*sint(3,1)
    pexp2 = 2.0d0*pexp
    do n=2,lab,1
      o2psi = o2pdi
      o2pdi = dfloat(n-1)/pexp2
      sxone = sxtwo
      syone = sytwo
      szone = sztwo
      sxtwo = sint(1,n-1)
      sytwo = sint(2,n-1)
      sztwo = sint(3,n-1)
      sint(1,n) = rxpa*sxtwo     + o2psi*sxone
      sint(2,n) = rypa*sytwo     + o2psi*syone
      sint(3,n) = rzpa*sztwo     + o2psi*szone
      dipi(1,n) = rxpc*sint(1,n) + o2pdi*sxtwo
      dipi(2,n) = rypc*sint(2,n) + o2pdi*sytwo
      dipi(3,n) = rzpc*sint(3,n) + o2pdi*sztwo
    enddo
  
    li1=0
    n=max(1,la)
    do l=n,lab,1
      do lx=1,l
        ll=l-lx+1
        do ly=1,lx,1
          lm =lx-ly+1
          ln =ly
          li1=li1+1
          fintx(li1)=cosr*dipi(1,ll)*sint(2,lm)*sint(3,ln)
          finty(li1)=cosr*sint(1,ll)*dipi(2,lm)*sint(3,ln)
          fintz(li1)=cosr*sint(1,ll)*sint(2,lm)*dipi(3,ln)
        enddo
      enddo
    enddo
  end subroutine 

! ---------------------------------------------------------------------------------------
! CLASS_NUCLEAR
!
!note: na vstupe musi byt pre RKB bazu: (la-1), (lc-1), (lab+2), (lcd+2)
!note: na vstupe musi byt pre 2-centrove integraly, lcd=1 
!note: lf1 (on input ) - dlzka finalneho batchu (prvy  electron),t.j.  ncc(la)*ncc(lb)(*5)
!note: lf2 (on input ) - dlzka finalneho batchu (druhy electron),t.j.  ncc(lc)*ncc(ld)(*5)
!note: li1 (on output) - dlzka intermediates    (prvy  electron)
!note: li2 (on output) - dlzka intermediates    (druhy electron)
!
!note: na vystupe, 1-elektronove integraly: [la...lab|0]
!note: na vystupe, 2-elektronove integraly: [lc...lcd,0|la....lab,0], left-side first 
!
! ---------------------------------------------------------------------------------------
!        |  NAI-PN    NAI-GCDM       ERI
! fnbra  |   1.0      sigma/p      q/(p+q)
! fdbra  |   1/2p      1/2p         1/2p
! fnket  |   rdum      rdum        p/(p+q)
! fdket  |   rdum      rdum         1/2q
! ---------------------------------------------------------------------------------------
  subroutine interest_osr_class_nuclear(li1,li2,fint,cosr,t,                                 &
                                        la,lab,lf1,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                        lc,lcd,lf2,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
    !-- output --!
    integer, intent(out) :: li1,li2 
    real(8), intent(out) :: fint(*)

    !-- input --!
    real(8), intent(in) :: cosr,t
    integer, intent(in) :: la,lab,lf1
    integer, intent(in) :: lc,lcd,lf2
    real(8), intent(in) :: fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp
    real(8), intent(in) :: fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq 

    !-- local --!
    integer :: iadr
    real(8) :: tdif
    integer :: m,mmax,nmin 


    mmax = lab + lcd - 1 
    nmin = max(1,(la*(la-1)*(la+1))/6+1)

    !-- [ 0 | * ]--!
    if( t.le.12.0d0 )then
      iadr = nint(10.0d0*t)
      tdif = t - 0.1d0*dfloat(iadr)
      iadr = iadr*lmxx+1
      call nuclear_boysxina(mmax,t,tdif,fnx,fnxtab(iadr))
    elseif( t.gt.12.0d0 .and. t.lt.tf) then         
      call nuclear_boysxinb(mmax,t,fnx)
    else
      call nuclear_boysxinc(mmax,t,fnx)
    endif

    if( mmax.eq.1 )then
      li1=1 
      li2=1 
      fint(1) = cosr*fnx(1)
      return
    endif

    do m=1,mmax 
      vaux(1,m) = cosr*fnx(m)
    enddo

    if( lab.eq.1 )then
      gama(0) = 0.0d0
      call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,1,lc,lcd,     &
                       fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      return
    endif

    !-- [ 1 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      vaux(2,m) = rxpa*vaux(1,m) + rxwp*vaux(1,m+1)
      vaux(3,m) = rypa*vaux(1,m) + rywp*vaux(1,m+1)
      vaux(4,m) = rzpa*vaux(1,m) + rzwp*vaux(1,m+1)
    enddo
    if( lab.eq.2 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,4,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,4,lc,lcd,     &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 2 | * ] --!
    f(1)=fdbra
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(1,m) - fnbra*vaux(1,m+1))

      vaux(  5,m) = rxpa*vaux(  2,m) + rxwp*vaux(  2,m+1) + c1   
      vaux(  6,m) = rxpa*vaux(  3,m) + rxwp*vaux(  3,m+1)
      vaux(  7,m) = rxpa*vaux(  4,m) + rxwp*vaux(  4,m+1)
      vaux(  8,m) = rypa*vaux(  3,m) + rywp*vaux(  3,m+1) + c1
      vaux(  9,m) = rzpa*vaux(  3,m) + rzwp*vaux(  3,m+1)
      vaux( 10,m) = rzpa*vaux(  4,m) + rzwp*vaux(  4,m+1) + c1
    enddo
    if( lab.eq.3 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,10,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,10,lc,lcd,    &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 3 | * ] --!
    f(2)=f(1)+f(1)
    mmax=mmax-1
    do m=1,mmax
      c1 = f(2)*(vaux(2,m) - fnbra*vaux(2,m+1))
      c2 = f(2)*(vaux(3,m) - fnbra*vaux(3,m+1))
      c3 = f(2)*(vaux(4,m) - fnbra*vaux(4,m+1))

      vaux( 11,m) = rxpa*vaux(  5,m) + rxwp*vaux(  5,m+1) + c1
      vaux( 12,m) = rypa*vaux(  5,m) + rywp*vaux(  5,m+1)
      vaux( 13,m) = rzpa*vaux(  5,m) + rzwp*vaux(  5,m+1)
      vaux( 14,m) = rxpa*vaux(  8,m) + rxwp*vaux(  8,m+1)
      vaux( 15,m) = rxpa*vaux(  9,m) + rxwp*vaux(  9,m+1)
      vaux( 16,m) = rxpa*vaux( 10,m) + rxwp*vaux( 10,m+1)
      vaux( 17,m) = rypa*vaux(  8,m) + rywp*vaux(  8,m+1) + c2
      vaux( 18,m) = rzpa*vaux(  8,m) + rzwp*vaux(  8,m+1)
      vaux( 19,m) = rypa*vaux( 10,m) + rywp*vaux( 10,m+1)
      vaux( 20,m) = rzpa*vaux( 10,m) + rzwp*vaux( 10,m+1) + c3
    enddo
    if( lab.eq.4 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,20,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,20,lc,lcd,    &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 4 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux( 5,m) - fnbra*vaux( 5,m+1))
      c2 = f(1)*(vaux( 8,m) - fnbra*vaux( 8,m+1))
      c3 = f(1)*(vaux(10,m) - fnbra*vaux(10,m+1))

      vaux( 21,m) = rxpa*vaux( 11,m) + rxwp*vaux( 11,m+1) + c1*3.0d0
      vaux( 22,m) = rypa*vaux( 11,m) + rywp*vaux( 11,m+1)
      vaux( 23,m) = rzpa*vaux( 11,m) + rzwp*vaux( 11,m+1)
      vaux( 24,m) = rypa*vaux( 12,m) + rywp*vaux( 12,m+1) + c1
      vaux( 25,m) = rzpa*vaux( 12,m) + rzwp*vaux( 12,m+1)
      vaux( 26,m) = rzpa*vaux( 13,m) + rzwp*vaux( 13,m+1) + c1
      vaux( 27,m) = rxpa*vaux( 17,m) + rxwp*vaux( 17,m+1)
      vaux( 28,m) = rxpa*vaux( 18,m) + rxwp*vaux( 18,m+1)
      vaux( 29,m) = rxpa*vaux( 19,m) + rxwp*vaux( 19,m+1)
      vaux( 30,m) = rxpa*vaux( 20,m) + rxwp*vaux( 20,m+1)
      vaux( 31,m) = rypa*vaux( 17,m) + rywp*vaux( 17,m+1) + c2*3.0d0
      vaux( 32,m) = rzpa*vaux( 17,m) + rzwp*vaux( 17,m+1)
      vaux( 33,m) = rzpa*vaux( 18,m) + rzwp*vaux( 18,m+1) + c2
      vaux( 34,m) = rypa*vaux( 20,m) + rywp*vaux( 20,m+1)
      vaux( 35,m) = rzpa*vaux( 20,m) + rzwp*vaux( 20,m+1) + c3*3.0d0
    enddo
    if( lab.eq.5 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,35,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,35,lc,lcd,    &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 5 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(11,m) - fnbra*vaux(11,m+1))
      c2 = f(1)*(vaux(17,m) - fnbra*vaux(17,m+1))
      c3 = f(1)*(vaux(20,m) - fnbra*vaux(20,m+1))

      vaux( 36,m) = rxpa*vaux( 21,m) + rxwp*vaux( 21,m+1) + c1*4.0d0
      vaux( 37,m) = rypa*vaux( 21,m) + rywp*vaux( 21,m+1)
      vaux( 38,m) = rzpa*vaux( 21,m) + rzwp*vaux( 21,m+1)
      vaux( 39,m) = rypa*vaux( 22,m) + rywp*vaux( 22,m+1) + c1         
      vaux( 40,m) = rzpa*vaux( 22,m) + rzwp*vaux( 22,m+1)
      vaux( 41,m) = rzpa*vaux( 23,m) + rzwp*vaux( 23,m+1) + c1        
      vaux( 42,m) = rxpa*vaux( 27,m) + rxwp*vaux( 27,m+1) + c2       
      vaux( 43,m) = rzpa*vaux( 24,m) + rzwp*vaux( 24,m+1)
      vaux( 44,m) = rypa*vaux( 26,m) + rywp*vaux( 26,m+1)
      vaux( 45,m) = rxpa*vaux( 30,m) + rxwp*vaux( 30,m+1) + c3        
      vaux( 46,m) = rxpa*vaux( 31,m) + rxwp*vaux( 31,m+1)
      vaux( 47,m) = rxpa*vaux( 32,m) + rxwp*vaux( 32,m+1)
      vaux( 48,m) = rxpa*vaux( 33,m) + rxwp*vaux( 33,m+1)
      vaux( 49,m) = rxpa*vaux( 34,m) + rxwp*vaux( 34,m+1)
      vaux( 50,m) = rxpa*vaux( 35,m) + rxwp*vaux( 35,m+1)
      vaux( 51,m) = rypa*vaux( 31,m) + rywp*vaux( 31,m+1) + c2*4.0d0
      vaux( 52,m) = rzpa*vaux( 31,m) + rzwp*vaux( 31,m+1)
      vaux( 53,m) = rzpa*vaux( 32,m) + rzwp*vaux( 32,m+1) + c2         
      vaux( 54,m) = rypa*vaux( 34,m) + rywp*vaux( 34,m+1) + c3          
      vaux( 55,m) = rypa*vaux( 35,m) + rywp*vaux( 35,m+1)
      vaux( 56,m) = rzpa*vaux( 35,m) + rzwp*vaux( 35,m+1) + c3*4.0d0
    enddo
    if( lab.eq.6 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,56,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        gama(5) = gama(4) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,56,lc,lcd,    &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 6 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(21,m) - fnbra*vaux(21,m+1))
      c4 = f(2)*(vaux(22,m) - fnbra*vaux(22,m+1))
      c5 = f(2)*(vaux(23,m) - fnbra*vaux(23,m+1))
      c6 = f(1)*(vaux(24,m) - fnbra*vaux(24,m+1))
      c2 = f(1)*(vaux(31,m) - fnbra*vaux(31,m+1))
      c7 = f(2)*(vaux(32,m) - fnbra*vaux(32,m+1))
      c3 = f(1)*(vaux(35,m) - fnbra*vaux(35,m+1))

      vaux( 57,m) = rxpa*vaux( 36,m) + rxwp*vaux( 36,m+1) + c1*5.0d0
      vaux( 58,m) = rypa*vaux( 36,m) + rywp*vaux( 36,m+1)
      vaux( 59,m) = rzpa*vaux( 36,m) + rzwp*vaux( 36,m+1)
      vaux( 60,m) = rypa*vaux( 37,m) + rywp*vaux( 37,m+1) + c1 
      vaux( 61,m) = rzpa*vaux( 37,m) + rzwp*vaux( 37,m+1)
      vaux( 62,m) = rzpa*vaux( 38,m) + rzwp*vaux( 38,m+1) + c1 
      vaux( 63,m) = rypa*vaux( 39,m) + rywp*vaux( 39,m+1) + c4 
      vaux( 64,m) = rzpa*vaux( 39,m) + rzwp*vaux( 39,m+1)
      vaux( 65,m) = rypa*vaux( 41,m) + rywp*vaux( 41,m+1)
      vaux( 66,m) = rzpa*vaux( 41,m) + rzwp*vaux( 41,m+1) + c5 
      vaux( 67,m) = rxpa*vaux( 46,m) + rxwp*vaux( 46,m+1) + c2 
      vaux( 68,m) = rzpa*vaux( 42,m) + rzwp*vaux( 42,m+1)
      vaux( 69,m) = rzpa*vaux( 43,m) + rzwp*vaux( 43,m+1) + c6 
      vaux( 70,m) = rypa*vaux( 45,m) + rywp*vaux( 45,m+1)
      vaux( 71,m) = rxpa*vaux( 50,m) + rxwp*vaux( 50,m+1) + c3 
      vaux( 72,m) = rxpa*vaux( 51,m) + rxwp*vaux( 51,m+1)
      vaux( 73,m) = rxpa*vaux( 52,m) + rxwp*vaux( 52,m+1)
      vaux( 74,m) = rxpa*vaux( 53,m) + rxwp*vaux( 53,m+1)
      vaux( 75,m) = rxpa*vaux( 54,m) + rxwp*vaux( 54,m+1)
      vaux( 76,m) = rxpa*vaux( 55,m) + rxwp*vaux( 55,m+1)
      vaux( 77,m) = rxpa*vaux( 56,m) + rxwp*vaux( 56,m+1)
      vaux( 78,m) = rypa*vaux( 51,m) + rywp*vaux( 51,m+1) + c2*5.0d0
      vaux( 79,m) = rzpa*vaux( 51,m) + rzwp*vaux( 51,m+1)
      vaux( 80,m) = rzpa*vaux( 52,m) + rzwp*vaux( 52,m+1) + c2 
      vaux( 81,m) = rzpa*vaux( 53,m) + rzwp*vaux( 53,m+1) + c7 
      vaux( 82,m) = rypa*vaux( 55,m) + rywp*vaux( 55,m+1) + c3 
      vaux( 83,m) = rypa*vaux( 56,m) + rywp*vaux( 56,m+1)
      vaux( 84,m) = rzpa*vaux( 56,m) + rzwp*vaux( 56,m+1) + c3*5.0d0
    enddo
    if( lab.eq.7 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,84,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        gama(5) = gama(4) + gama(1) 
        gama(6) = gama(5) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,84,lc,lcd,    &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 7 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux( 36,m) - fnbra*vaux( 36,m+1))
      c4 = f(2)*(vaux( 37,m) - fnbra*vaux( 37,m+1))
      c5 = f(2)*(vaux( 38,m) - fnbra*vaux( 38,m+1))
      c6 = f(1)*(vaux( 39,m) - fnbra*vaux( 39,m+1))
      c7 = f(1)*(vaux( 42,m) - fnbra*vaux( 42,m+1))
      c8 = f(1)*(vaux( 45,m) - fnbra*vaux( 45,m+1))
      c9 = f(2)*(vaux( 46,m) - fnbra*vaux( 46,m+1))
      c10= f(2)*(vaux( 50,m) - fnbra*vaux( 50,m+1))
      c2 = f(1)*(vaux( 51,m) - fnbra*vaux( 51,m+1))
      c11= f(2)*(vaux( 52,m) - fnbra*vaux( 52,m+1))
      c12= f(2)*(vaux( 55,m) - fnbra*vaux( 55,m+1))
      c3 = f(1)*(vaux( 56,m) - fnbra*vaux( 56,m+1))

      vaux( 85,m) = rxpa*vaux( 57,m) + rxwp*vaux( 57,m+1) + c1*6.0d0 
      vaux( 86,m) = rypa*vaux( 57,m) + rywp*vaux( 57,m+1)
      vaux( 87,m) = rzpa*vaux( 57,m) + rzwp*vaux( 57,m+1)
      vaux( 88,m) = rypa*vaux( 58,m) + rywp*vaux( 58,m+1) + c1 
      vaux( 89,m) = rzpa*vaux( 58,m) + rzwp*vaux( 58,m+1)
      vaux( 90,m) = rzpa*vaux( 59,m) + rzwp*vaux( 59,m+1) + c1 
      vaux( 91,m) = rypa*vaux( 60,m) + rywp*vaux( 60,m+1) + c4 
      vaux( 92,m) = rzpa*vaux( 60,m) + rzwp*vaux( 60,m+1)
      vaux( 93,m) = rypa*vaux( 62,m) + rywp*vaux( 62,m+1)
      vaux( 94,m) = rzpa*vaux( 62,m) + rzwp*vaux( 62,m+1) + c5 
      vaux( 95,m) = rxpa*vaux( 67,m) + rxwp*vaux( 67,m+1) + c9 
      vaux( 96,m) = rzpa*vaux( 63,m) + rzwp*vaux( 63,m+1)
      vaux( 97,m) = rzpa*vaux( 64,m) + rzwp*vaux( 64,m+1) + c6 
      vaux( 98,m) = rypa*vaux( 66,m) + rywp*vaux( 66,m+1)
      vaux( 99,m) = rxpa*vaux( 71,m) + rxwp*vaux( 71,m+1) + c10
      vaux(100,m) = rxpa*vaux( 72,m) + rxwp*vaux( 72,m+1) + c2 
      vaux(101,m) = rzpa*vaux( 67,m) + rzwp*vaux( 67,m+1)
      vaux(102,m) = rzpa*vaux( 68,m) + rzwp*vaux( 68,m+1) + c7 
      vaux(103,m) = rypa*vaux( 70,m) + rywp*vaux( 70,m+1) + c8 
      vaux(104,m) = rypa*vaux( 71,m) + rywp*vaux( 71,m+1)
      vaux(105,m) = rxpa*vaux( 77,m) + rxwp*vaux( 77,m+1) + c3 
      vaux(106,m) = rxpa*vaux( 78,m) + rxwp*vaux( 78,m+1)
      vaux(107,m) = rxpa*vaux( 79,m) + rxwp*vaux( 79,m+1)
      vaux(108,m) = rxpa*vaux( 80,m) + rxwp*vaux( 80,m+1)
      vaux(109,m) = rxpa*vaux( 81,m) + rxwp*vaux( 81,m+1)
      vaux(110,m) = rxpa*vaux( 82,m) + rxwp*vaux( 82,m+1)
      vaux(111,m) = rxpa*vaux( 83,m) + rxwp*vaux( 83,m+1)
      vaux(112,m) = rxpa*vaux( 84,m) + rxwp*vaux( 84,m+1)
      vaux(113,m) = rypa*vaux( 78,m) + rywp*vaux( 78,m+1) + c2*6.0d0 
      vaux(114,m) = rzpa*vaux( 78,m) + rzwp*vaux( 78,m+1)
      vaux(115,m) = rzpa*vaux( 79,m) + rzwp*vaux( 79,m+1) + c2 
      vaux(116,m) = rzpa*vaux( 80,m) + rzwp*vaux( 80,m+1) + c11 
      vaux(117,m) = rypa*vaux( 82,m) + rywp*vaux( 82,m+1) + c12 
      vaux(118,m) = rypa*vaux( 83,m) + rywp*vaux( 83,m+1) + c3 
      vaux(119,m) = rypa*vaux( 84,m) + rywp*vaux( 84,m+1)
      vaux(120,m) = rzpa*vaux( 84,m) + rzwp*vaux( 84,m+1) + c3*6.0d0 
    enddo
    if( lab.eq.8 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,120,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        gama(5) = gama(4) + gama(1) 
        gama(6) = gama(5) + gama(1) 
        gama(7) = gama(6) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,120,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 8 | * ]--!
    f(3)=f(2)+f(1)
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux( 57,m) - fnbra*vaux( 57,m+1))
      c4 = f(2)*(vaux( 58,m) - fnbra*vaux( 58,m+1))
      c5 = f(2)*(vaux( 59,m) - fnbra*vaux( 59,m+1))
      c6 = f(1)*(vaux( 60,m) - fnbra*vaux( 60,m+1))
      c7 = f(3)*(vaux( 62,m) - fnbra*vaux( 62,m+1))
      c8 = f(1)*(vaux( 63,m) - fnbra*vaux( 63,m+1))
      c9 = f(1)*(vaux( 66,m) - fnbra*vaux( 66,m+1))
      c10= f(1)*(vaux( 67,m) - fnbra*vaux( 67,m+1))
      c11= f(1)*(vaux( 71,m) - fnbra*vaux( 71,m+1))
      c12= f(2)*(vaux( 72,m) - fnbra*vaux( 72,m+1))
      c13= f(2)*(vaux( 77,m) - fnbra*vaux( 77,m+1))
      c2 = f(1)*(vaux( 78,m) - fnbra*vaux( 78,m+1))
      c14= f(2)*(vaux( 79,m) - fnbra*vaux( 79,m+1))
      c15= f(3)*(vaux( 80,m) - fnbra*vaux( 80,m+1))
      c16= f(1)*(vaux( 81,m) - fnbra*vaux( 81,m+1))
      c17= f(2)*(vaux( 83,m) - fnbra*vaux( 83,m+1))
      c3 = f(1)*(vaux( 84,m) - fnbra*vaux( 84,m+1))

      vaux(121,m) = rxpa*vaux( 85,m) + rxwp*vaux( 85,m+1) + c1*7.0d0 
      vaux(122,m) = rypa*vaux( 85,m) + rywp*vaux( 85,m+1)
      vaux(123,m) = rzpa*vaux( 85,m) + rzwp*vaux( 85,m+1)
      vaux(124,m) = rypa*vaux( 86,m) + rywp*vaux( 86,m+1) + c1 
      vaux(125,m) = rzpa*vaux( 86,m) + rzwp*vaux( 86,m+1)
      vaux(126,m) = rzpa*vaux( 87,m) + rzwp*vaux( 87,m+1) + c1 
      vaux(127,m) = rypa*vaux( 88,m) + rywp*vaux( 88,m+1) + c4 
      vaux(128,m) = rzpa*vaux( 88,m) + rzwp*vaux( 88,m+1)
      vaux(129,m) = rypa*vaux( 90,m) + rywp*vaux( 90,m+1)
      vaux(130,m) = rzpa*vaux( 90,m) + rzwp*vaux( 90,m+1) + c5 
      vaux(131,m) = rypa*vaux( 91,m) + rywp*vaux( 91,m+1) + c6*3.0d0
      vaux(132,m) = rzpa*vaux( 91,m) + rzwp*vaux( 91,m+1)
      vaux(133,m) = rzpa*vaux( 92,m) + rzwp*vaux( 92,m+1) + c6 
      vaux(134,m) = rypa*vaux( 94,m) + rywp*vaux( 94,m+1)
      vaux(135,m) = rzpa*vaux( 94,m) + rzwp*vaux( 94,m+1) + c7 
      vaux(136,m) = rxpa*vaux(100,m) + rxwp*vaux(100,m+1) + c12 
      vaux(137,m) = rzpa*vaux( 95,m) + rzwp*vaux( 95,m+1)
      vaux(138,m) = rzpa*vaux( 96,m) + rzwp*vaux( 96,m+1) + c8 
      vaux(139,m) = rypa*vaux( 98,m) + rywp*vaux( 98,m+1) + c9 
      vaux(140,m) = rypa*vaux( 99,m) + rywp*vaux( 99,m+1)
      vaux(141,m) = rxpa*vaux(105,m) + rxwp*vaux(105,m+1) + c13 
      vaux(142,m) = rxpa*vaux(106,m) + rxwp*vaux(106,m+1) + c2 
      vaux(143,m) = rzpa*vaux(100,m) + rzwp*vaux(100,m+1)
      vaux(144,m) = rzpa*vaux(101,m) + rzwp*vaux(101,m+1) + c10 
      vaux(145,m) = rxpa*vaux(109,m) + rxwp*vaux(109,m+1) + c16 
      vaux(146,m) = rypa*vaux(104,m) + rywp*vaux(104,m+1) + c11 
      vaux(147,m) = rypa*vaux(105,m) + rywp*vaux(105,m+1)
      vaux(148,m) = rxpa*vaux(112,m) + rxwp*vaux(112,m+1) + c3 
      vaux(149,m) = rxpa*vaux(113,m) + rxwp*vaux(113,m+1)
      vaux(150,m) = rxpa*vaux(114,m) + rxwp*vaux(114,m+1)
      vaux(151,m) = rxpa*vaux(115,m) + rxwp*vaux(115,m+1)
      vaux(152,m) = rxpa*vaux(116,m) + rxwp*vaux(116,m+1)
      vaux(153,m) = rxpa*vaux(117,m) + rxwp*vaux(117,m+1)
      vaux(154,m) = rxpa*vaux(118,m) + rxwp*vaux(118,m+1)
      vaux(155,m) = rxpa*vaux(119,m) + rxwp*vaux(119,m+1)
      vaux(156,m) = rxpa*vaux(120,m) + rxwp*vaux(120,m+1)
      vaux(157,m) = rypa*vaux(113,m) + rywp*vaux(113,m+1) + c2*7.0d0 
      vaux(158,m) = rzpa*vaux(113,m) + rzwp*vaux(113,m+1)
      vaux(159,m) = rzpa*vaux(114,m) + rzwp*vaux(114,m+1) + c2 
      vaux(160,m) = rzpa*vaux(115,m) + rzwp*vaux(115,m+1) + c14 
      vaux(161,m) = rzpa*vaux(116,m) + rzwp*vaux(116,m+1) + c15 
      vaux(162,m) = rypa*vaux(118,m) + rywp*vaux(118,m+1) + c17 
      vaux(163,m) = rypa*vaux(119,m) + rywp*vaux(119,m+1) + c3 
      vaux(164,m) = rypa*vaux(120,m) + rywp*vaux(120,m+1)
      vaux(165,m) = rzpa*vaux(120,m) + rzwp*vaux(120,m+1) + c3*7.0d0 
    enddo
    if( lab.eq.9 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,165,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        gama(5) = gama(4) + gama(1) 
        gama(6) = gama(5) + gama(1) 
        gama(7) = gama(6) + gama(1) 
        gama(8) = gama(7) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,165,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 9 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux( 85,m) - fnbra*vaux( 85,m+1))
      c2 = f(2)*(vaux( 86,m) - fnbra*vaux( 86,m+1))
      c3 = f(2)*(vaux( 87,m) - fnbra*vaux( 87,m+1))
      c4 = f(1)*(vaux( 88,m) - fnbra*vaux( 88,m+1))
      c5 = f(3)*(vaux( 90,m) - fnbra*vaux( 90,m+1))
      c6 = f(1)*(vaux( 91,m) - fnbra*vaux( 91,m+1))
      c7 = f(1)*(vaux( 94,m) - fnbra*vaux( 94,m+1))
      c8 = f(1)*(vaux( 95,m) - fnbra*vaux( 95,m+1))
      c9 = f(2)*(vaux( 96,m) - fnbra*vaux( 96,m+1))
      c10= f(1)*(vaux( 99,m) - fnbra*vaux( 99,m+1))
      c11= f(1)*(vaux(100,m) - fnbra*vaux(100,m+1))
      c12= f(1)*(vaux(105,m) - fnbra*vaux(105,m+1))
      c13= f(2)*(vaux(106,m) - fnbra*vaux(106,m+1))
      c14= f(2)*(vaux(112,m) - fnbra*vaux(112,m+1))
      c15= f(1)*(vaux(113,m) - fnbra*vaux(113,m+1))
      c16= f(2)*(vaux(114,m) - fnbra*vaux(114,m+1))
      c17= f(3)*(vaux(115,m) - fnbra*vaux(115,m+1))
      c18= f(1)*(vaux(116,m) - fnbra*vaux(116,m+1))
      c19= f(1)*(vaux(117,m) - fnbra*vaux(117,m+1))
      c20= f(3)*(vaux(118,m) - fnbra*vaux(118,m+1))
      c21= f(2)*(vaux(119,m) - fnbra*vaux(119,m+1))
      c22= f(1)*(vaux(120,m) - fnbra*vaux(120,m+1))

      vaux(166,m) = rxpa*vaux(121,m) + rxwp*vaux(121,m+1) + c1*8.0d0 
      vaux(167,m) = rypa*vaux(121,m) + rywp*vaux(121,m+1)
      vaux(168,m) = rzpa*vaux(121,m) + rzwp*vaux(121,m+1)
      vaux(169,m) = rypa*vaux(122,m) + rywp*vaux(122,m+1) + c1 
      vaux(170,m) = rzpa*vaux(122,m) + rzwp*vaux(122,m+1)
      vaux(171,m) = rzpa*vaux(123,m) + rzwp*vaux(123,m+1) + c1 
      vaux(172,m) = rypa*vaux(124,m) + rywp*vaux(124,m+1) + c2 
      vaux(173,m) = rzpa*vaux(124,m) + rzwp*vaux(124,m+1)
      vaux(174,m) = rypa*vaux(126,m) + rywp*vaux(126,m+1)
      vaux(175,m) = rzpa*vaux(126,m) + rzwp*vaux(126,m+1) + c3 
      vaux(176,m) = rypa*vaux(127,m) + rywp*vaux(127,m+1) + c4*3.0d0
      vaux(177,m) = rzpa*vaux(127,m) + rzwp*vaux(127,m+1)
      vaux(178,m) = rzpa*vaux(128,m) + rzwp*vaux(128,m+1) + c4 
      vaux(179,m) = rypa*vaux(130,m) + rywp*vaux(130,m+1)
      vaux(180,m) = rzpa*vaux(130,m) + rzwp*vaux(130,m+1) + c5 
      vaux(181,m) = rxpa*vaux(136,m) + rxwp*vaux(136,m+1) + c11*3.0d0 
      vaux(182,m) = rzpa*vaux(131,m) + rzwp*vaux(131,m+1)
      vaux(183,m) = rzpa*vaux(132,m) + rzwp*vaux(132,m+1) + c6 
      vaux(184,m) = rypa*vaux(134,m) + rywp*vaux(134,m+1) + c7 
      vaux(185,m) = rypa*vaux(135,m) + rywp*vaux(135,m+1)
      vaux(186,m) = rxpa*vaux(141,m) + rxwp*vaux(141,m+1) + c12*3.0d0 
      vaux(187,m) = rxpa*vaux(142,m) + rxwp*vaux(142,m+1) + c13 
      vaux(188,m) = rzpa*vaux(136,m) + rzwp*vaux(136,m+1)
      vaux(189,m) = rzpa*vaux(137,m) + rzwp*vaux(137,m+1) + c8 
      vaux(190,m) = rzpa*vaux(138,m) + rzwp*vaux(138,m+1) + c9 
      vaux(191,m) = rypa*vaux(140,m) + rywp*vaux(140,m+1) + c10 
      vaux(192,m) = rypa*vaux(141,m) + rywp*vaux(141,m+1)
      vaux(193,m) = rxpa*vaux(148,m) + rxwp*vaux(148,m+1) + c14 
      vaux(194,m) = rxpa*vaux(149,m) + rxwp*vaux(149,m+1) + c15 
      vaux(195,m) = rzpa*vaux(142,m) + rzwp*vaux(142,m+1)
      vaux(196,m) = rzpa*vaux(143,m) + rzwp*vaux(143,m+1) + c11 
      vaux(197,m) = rxpa*vaux(152,m) + rxwp*vaux(152,m+1) + c18 
      vaux(198,m) = rxpa*vaux(153,m) + rxwp*vaux(153,m+1) + c19 
      vaux(199,m) = rypa*vaux(147,m) + rywp*vaux(147,m+1) + c12 
      vaux(200,m) = rypa*vaux(148,m) + rywp*vaux(148,m+1)
      vaux(201,m) = rxpa*vaux(156,m) + rxwp*vaux(156,m+1) + c22 
      vaux(202,m) = rxpa*vaux(157,m) + rxwp*vaux(157,m+1)
      vaux(203,m) = rxpa*vaux(158,m) + rxwp*vaux(158,m+1)
      vaux(204,m) = rxpa*vaux(159,m) + rxwp*vaux(159,m+1)
      vaux(205,m) = rxpa*vaux(160,m) + rxwp*vaux(160,m+1)
      vaux(206,m) = rxpa*vaux(161,m) + rxwp*vaux(161,m+1)
      vaux(207,m) = rxpa*vaux(162,m) + rxwp*vaux(162,m+1)
      vaux(208,m) = rxpa*vaux(163,m) + rxwp*vaux(163,m+1)
      vaux(209,m) = rxpa*vaux(164,m) + rxwp*vaux(164,m+1)
      vaux(210,m) = rxpa*vaux(165,m) + rxwp*vaux(165,m+1)
      vaux(211,m) = rypa*vaux(157,m) + rywp*vaux(157,m+1) + c15*8.0d0 
      vaux(212,m) = rzpa*vaux(157,m) + rzwp*vaux(157,m+1)
      vaux(213,m) = rzpa*vaux(158,m) + rzwp*vaux(158,m+1) + c15 
      vaux(214,m) = rzpa*vaux(159,m) + rzwp*vaux(159,m+1) + c16 
      vaux(215,m) = rzpa*vaux(160,m) + rzwp*vaux(160,m+1) + c17 
      vaux(216,m) = rypa*vaux(162,m) + rywp*vaux(162,m+1) + c20 
      vaux(217,m) = rypa*vaux(163,m) + rywp*vaux(163,m+1) + c21 
      vaux(218,m) = rypa*vaux(164,m) + rywp*vaux(164,m+1) + c22 
      vaux(219,m) = rypa*vaux(165,m) + rywp*vaux(165,m+1)
      vaux(220,m) = rzpa*vaux(165,m) + rzwp*vaux(165,m+1) + c22*8.0d0 
    enddo
    if( lab.eq.10 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,220,vaux,fint,li1,li2)
      else
        gama(0) = 0.0d0
        gama(1) = fnket*fdbra 
        gama(2) = gama(1) + gama(1) 
        gama(3) = gama(2) + gama(1) 
        gama(4) = gama(3) + gama(1) 
        gama(5) = gama(4) + gama(1) 
        gama(6) = gama(5) + gama(1) 
        gama(7) = gama(6) + gama(1) 
        gama(8) = gama(7) + gama(1) 
        gama(9) = gama(8) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,220,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [ 10 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(121,m) - fnbra*vaux(121,m+1))
      c2 = f(2)*(vaux(122,m) - fnbra*vaux(122,m+1))
      c3 = f(2)*(vaux(123,m) - fnbra*vaux(123,m+1))
      c4 = f(1)*(vaux(124,m) - fnbra*vaux(124,m+1))
      c5 = f(3)*(vaux(126,m) - fnbra*vaux(126,m+1))
      c6 = f(1)*(vaux(127,m) - fnbra*vaux(127,m+1))
      c7 = f(1)*(vaux(130,m) - fnbra*vaux(130,m+1))
      c8 = f(1)*(vaux(131,m) - fnbra*vaux(131,m+1))
      c9 = f(2)*(vaux(132,m) - fnbra*vaux(132,m+1))
      c10= f(1)*(vaux(135,m) - fnbra*vaux(135,m+1))
      c11= f(1)*(vaux(136,m) - fnbra*vaux(136,m+1))
      c12= f(2)*(vaux(137,m) - fnbra*vaux(137,m+1))
      c13= f(2)*(vaux(140,m) - fnbra*vaux(140,m+1))
      c14= f(1)*(vaux(141,m) - fnbra*vaux(141,m+1))
      c15= f(1)*(vaux(142,m) - fnbra*vaux(142,m+1))
      c16= f(1)*(vaux(148,m) - fnbra*vaux(148,m+1))
      c17= f(2)*(vaux(149,m) - fnbra*vaux(149,m+1))
      c18= f(2)*(vaux(156,m) - fnbra*vaux(156,m+1))
      c19= f(1)*(vaux(157,m) - fnbra*vaux(157,m+1))
      c20= f(2)*(vaux(158,m) - fnbra*vaux(158,m+1))
      c21= f(3)*(vaux(159,m) - fnbra*vaux(159,m+1))
      c22= f(1)*(vaux(160,m) - fnbra*vaux(160,m+1))
      c23= f(1)*(vaux(161,m) - fnbra*vaux(161,m+1))
      c24= f(1)*(vaux(162,m) - fnbra*vaux(162,m+1))
      c25= f(3)*(vaux(163,m) - fnbra*vaux(163,m+1))
      c26= f(2)*(vaux(164,m) - fnbra*vaux(164,m+1))
      c27= f(1)*(vaux(165,m) - fnbra*vaux(165,m+1))

      vaux(221,m) = rxpa*vaux(166,m) + rxwp*vaux(166,m+1) + c1*9.0d0 
      vaux(222,m) = rypa*vaux(166,m) + rywp*vaux(166,m+1)
      vaux(223,m) = rzpa*vaux(166,m) + rzwp*vaux(166,m+1)
      vaux(224,m) = rypa*vaux(167,m) + rywp*vaux(167,m+1) + c1 
      vaux(225,m) = rzpa*vaux(167,m) + rzwp*vaux(167,m+1)
      vaux(226,m) = rzpa*vaux(168,m) + rzwp*vaux(168,m+1) + c1 
      vaux(227,m) = rypa*vaux(169,m) + rywp*vaux(169,m+1) + c2 
      vaux(228,m) = rzpa*vaux(169,m) + rzwp*vaux(169,m+1)
      vaux(229,m) = rypa*vaux(171,m) + rywp*vaux(171,m+1)
      vaux(230,m) = rzpa*vaux(171,m) + rzwp*vaux(171,m+1) + c3 
      vaux(231,m) = rypa*vaux(172,m) + rywp*vaux(172,m+1) + c4*3.0d0 
      vaux(232,m) = rzpa*vaux(172,m) + rzwp*vaux(172,m+1)
      vaux(233,m) = rzpa*vaux(173,m) + rzwp*vaux(173,m+1) + c4 
      vaux(234,m) = rypa*vaux(175,m) + rywp*vaux(175,m+1)
      vaux(235,m) = rzpa*vaux(175,m) + rzwp*vaux(175,m+1) + c5 
      vaux(236,m) = rypa*vaux(176,m) + rywp*vaux(176,m+1) + c6*4.0d0 
      vaux(237,m) = rzpa*vaux(176,m) + rzwp*vaux(176,m+1)
      vaux(238,m) = rzpa*vaux(177,m) + rzwp*vaux(177,m+1) + c6 
      vaux(239,m) = rypa*vaux(179,m) + rywp*vaux(179,m+1) + c7 
      vaux(240,m) = rypa*vaux(180,m) + rywp*vaux(180,m+1)
      vaux(241,m) = rzpa*vaux(180,m) + rzwp*vaux(180,m+1) + c7*4.0d0 
      vaux(242,m) = rxpa*vaux(187,m) + rxwp*vaux(187,m+1) + c15*3.0d0 
      vaux(243,m) = rzpa*vaux(181,m) + rzwp*vaux(181,m+1)
      vaux(244,m) = rzpa*vaux(182,m) + rzwp*vaux(182,m+1) + c8 
      vaux(245,m) = rzpa*vaux(183,m) + rzwp*vaux(183,m+1) + c9 
      vaux(246,m) = rypa*vaux(185,m) + rywp*vaux(185,m+1) + c10 
      vaux(247,m) = rypa*vaux(186,m) + rywp*vaux(186,m+1)
      vaux(248,m) = rxpa*vaux(193,m) + rxwp*vaux(193,m+1) + c16*3.0d0 
      vaux(249,m) = rxpa*vaux(194,m) + rxwp*vaux(194,m+1) + c17 
      vaux(250,m) = rzpa*vaux(187,m) + rzwp*vaux(187,m+1)
      vaux(251,m) = rzpa*vaux(188,m) + rzwp*vaux(188,m+1) + c11 
      vaux(252,m) = rzpa*vaux(189,m) + rzwp*vaux(189,m+1) + c12 
      vaux(253,m) = rypa*vaux(191,m) + rywp*vaux(191,m+1) + c13 
      vaux(254,m) = rypa*vaux(192,m) + rywp*vaux(192,m+1) + c14 
      vaux(255,m) = rypa*vaux(193,m) + rywp*vaux(193,m+1)
      vaux(256,m) = rxpa*vaux(201,m) + rxwp*vaux(201,m+1) + c18 
      vaux(257,m) = rxpa*vaux(202,m) + rxwp*vaux(202,m+1) + c19 
      vaux(258,m) = rzpa*vaux(194,m) + rzwp*vaux(194,m+1)
      vaux(259,m) = rzpa*vaux(195,m) + rzwp*vaux(195,m+1) + c15 
      vaux(260,m) = rxpa*vaux(205,m) + rxwp*vaux(205,m+1) + c22 
      vaux(261,m) = rxpa*vaux(206,m) + rxwp*vaux(206,m+1) + c23 
      vaux(262,m) = rxpa*vaux(207,m) + rxwp*vaux(207,m+1) + c24 
      vaux(263,m) = rypa*vaux(200,m) + rywp*vaux(200,m+1) + c16 
      vaux(264,m) = rypa*vaux(201,m) + rywp*vaux(201,m+1)
      vaux(265,m) = rxpa*vaux(210,m) + rxwp*vaux(210,m+1) + c27 
      vaux(266,m) = rxpa*vaux(211,m) + rxwp*vaux(211,m+1)
      vaux(267,m) = rxpa*vaux(212,m) + rxwp*vaux(212,m+1)
      vaux(268,m) = rxpa*vaux(213,m) + rxwp*vaux(213,m+1)
      vaux(269,m) = rxpa*vaux(214,m) + rxwp*vaux(214,m+1)
      vaux(270,m) = rxpa*vaux(215,m) + rxwp*vaux(215,m+1)
      vaux(271,m) = rxpa*vaux(216,m) + rxwp*vaux(216,m+1)
      vaux(272,m) = rxpa*vaux(217,m) + rxwp*vaux(217,m+1)
      vaux(273,m) = rxpa*vaux(218,m) + rxwp*vaux(218,m+1)
      vaux(274,m) = rxpa*vaux(219,m) + rxwp*vaux(219,m+1)
      vaux(275,m) = rxpa*vaux(220,m) + rxwp*vaux(220,m+1)
      vaux(276,m) = rypa*vaux(211,m) + rywp*vaux(211,m+1) + c19*9.0d0 
      vaux(277,m) = rzpa*vaux(211,m) + rzwp*vaux(211,m+1)
      vaux(278,m) = rzpa*vaux(212,m) + rzwp*vaux(212,m+1) + c19 
      vaux(279,m) = rzpa*vaux(213,m) + rzwp*vaux(213,m+1) + c20 
      vaux(280,m) = rzpa*vaux(214,m) + rzwp*vaux(214,m+1) + c21 
      vaux(281,m) = rzpa*vaux(215,m) + rzwp*vaux(215,m+1) + c22*4.0d0 
      vaux(282,m) = rypa*vaux(217,m) + rywp*vaux(217,m+1) + c25 
      vaux(283,m) = rypa*vaux(218,m) + rywp*vaux(218,m+1) + c26 
      vaux(284,m) = rypa*vaux(219,m) + rywp*vaux(219,m+1) + c27 
      vaux(285,m) = rypa*vaux(220,m) + rywp*vaux(220,m+1)
      vaux(286,m) = rzpa*vaux(220,m) + rzwp*vaux(220,m+1) + c27*9.0d0 
    enddo
    if( lab.eq.11 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,286,vaux,fint,li1,li2)
      else
        gama( 0) = 0.0d0
        gama( 1) = fnket*fdbra 
        gama( 2) = gama(1) + gama(1) 
        gama( 3) = gama(2) + gama(1) 
        gama( 4) = gama(3) + gama(1) 
        gama( 5) = gama(4) + gama(1) 
        gama( 6) = gama(5) + gama(1) 
        gama( 7) = gama(6) + gama(1) 
        gama( 8) = gama(7) + gama(1) 
        gama( 9) = gama(8) + gama(1) 
        gama(10) = gama(9) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,286,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [11 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(166,m) - fnbra*vaux(166,m+1))
      c2 = f(2)*(vaux(167,m) - fnbra*vaux(167,m+1))
      c3 = f(2)*(vaux(168,m) - fnbra*vaux(168,m+1))
      c4 = f(1)*(vaux(169,m) - fnbra*vaux(169,m+1))
      c5 = f(3)*(vaux(171,m) - fnbra*vaux(171,m+1))
      c6 = f(1)*(vaux(172,m) - fnbra*vaux(172,m+1))
      c7 = f(1)*(vaux(175,m) - fnbra*vaux(175,m+1))
      c8 = f(1)*(vaux(176,m) - fnbra*vaux(176,m+1))
      c9 = f(2)*(vaux(177,m) - fnbra*vaux(177,m+1))
      c10= f(1)*(vaux(180,m) - fnbra*vaux(180,m+1))
      c11= f(1)*(vaux(181,m) - fnbra*vaux(181,m+1))
      c12= f(2)*(vaux(182,m) - fnbra*vaux(182,m+1))
      c13= f(2)*(vaux(185,m) - fnbra*vaux(185,m+1))
      c14= f(1)*(vaux(186,m) - fnbra*vaux(186,m+1))
      c15= f(1)*(vaux(187,m) - fnbra*vaux(187,m+1))
      c16= f(2)*(vaux(188,m) - fnbra*vaux(188,m+1))
      c17= f(2)*(vaux(192,m) - fnbra*vaux(192,m+1))
      c18= f(1)*(vaux(193,m) - fnbra*vaux(193,m+1))
      c19= f(1)*(vaux(194,m) - fnbra*vaux(194,m+1))
      c20= f(1)*(vaux(201,m) - fnbra*vaux(201,m+1))
      c21= f(2)*(vaux(202,m) - fnbra*vaux(202,m+1))
      c22= f(2)*(vaux(206,m) - fnbra*vaux(206,m+1))
      c23= f(2)*(vaux(210,m) - fnbra*vaux(210,m+1))
      c24= f(1)*(vaux(211,m) - fnbra*vaux(211,m+1))
      c25= f(2)*(vaux(212,m) - fnbra*vaux(212,m+1))
      c26= f(3)*(vaux(213,m) - fnbra*vaux(213,m+1))
      c27= f(1)*(vaux(214,m) - fnbra*vaux(214,m+1))
      c28= f(1)*(vaux(215,m) - fnbra*vaux(215,m+1))
      c29= f(1)*(vaux(216,m) - fnbra*vaux(216,m+1))
      c30= f(1)*(vaux(217,m) - fnbra*vaux(217,m+1))
      c31= f(3)*(vaux(218,m) - fnbra*vaux(218,m+1))
      c32= f(2)*(vaux(219,m) - fnbra*vaux(219,m+1))
      c33= f(1)*(vaux(220,m) - fnbra*vaux(220,m+1))

      vaux(287,m) = rxpa*vaux(221,m) + rxwp*vaux(221,m+1) + c1*10.0d0
      vaux(288,m) = rypa*vaux(221,m) + rywp*vaux(221,m+1)
      vaux(289,m) = rzpa*vaux(221,m) + rzwp*vaux(221,m+1)
      vaux(290,m) = rypa*vaux(222,m) + rywp*vaux(222,m+1) + c1 
      vaux(291,m) = rzpa*vaux(222,m) + rzwp*vaux(222,m+1)
      vaux(292,m) = rzpa*vaux(223,m) + rzwp*vaux(223,m+1) + c1
      vaux(293,m) = rypa*vaux(224,m) + rywp*vaux(224,m+1) + c2
      vaux(294,m) = rzpa*vaux(224,m) + rzwp*vaux(224,m+1)
      vaux(295,m) = rypa*vaux(226,m) + rywp*vaux(226,m+1)
      vaux(296,m) = rzpa*vaux(226,m) + rzwp*vaux(226,m+1) + c3 
      vaux(297,m) = rypa*vaux(227,m) + rywp*vaux(227,m+1) + c4*3.0d0
      vaux(298,m) = rzpa*vaux(227,m) + rzwp*vaux(227,m+1)
      vaux(299,m) = rzpa*vaux(228,m) + rzwp*vaux(228,m+1) + c4 
      vaux(300,m) = rypa*vaux(230,m) + rywp*vaux(230,m+1)
      vaux(301,m) = rzpa*vaux(230,m) + rzwp*vaux(230,m+1) + c5 
      vaux(302,m) = rypa*vaux(231,m) + rywp*vaux(231,m+1) + c6*4.0d0 
      vaux(303,m) = rzpa*vaux(231,m) + rzwp*vaux(231,m+1)
      vaux(304,m) = rzpa*vaux(232,m) + rzwp*vaux(232,m+1) + c6 
      vaux(305,m) = rypa*vaux(234,m) + rywp*vaux(234,m+1) + c7 
      vaux(306,m) = rypa*vaux(235,m) + rywp*vaux(235,m+1)
      vaux(307,m) = rzpa*vaux(235,m) + rzwp*vaux(235,m+1) + c7*4.0d0 
      vaux(308,m) = rxpa*vaux(242,m) + rxwp*vaux(242,m+1) + c15*4.0d0 
      vaux(309,m) = rzpa*vaux(236,m) + rzwp*vaux(236,m+1)
      vaux(310,m) = rzpa*vaux(237,m) + rzwp*vaux(237,m+1) + c8 
      vaux(311,m) = rzpa*vaux(238,m) + rzwp*vaux(238,m+1) + c9 
      vaux(312,m) = rypa*vaux(240,m) + rywp*vaux(240,m+1) + c10 
      vaux(313,m) = rypa*vaux(241,m) + rywp*vaux(241,m+1)
      vaux(314,m) = rxpa*vaux(248,m) + rxwp*vaux(248,m+1) + c18*4.0d0 
      vaux(315,m) = rxpa*vaux(249,m) + rxwp*vaux(249,m+1) + c19*3.0d0 
      vaux(316,m) = rzpa*vaux(242,m) + rzwp*vaux(242,m+1)
      vaux(317,m) = rzpa*vaux(243,m) + rzwp*vaux(243,m+1) + c11 
      vaux(318,m) = rzpa*vaux(244,m) + rzwp*vaux(244,m+1) + c12 
      vaux(319,m) = rypa*vaux(246,m) + rywp*vaux(246,m+1) + c13 
      vaux(320,m) = rypa*vaux(247,m) + rywp*vaux(247,m+1) + c14 
      vaux(321,m) = rypa*vaux(248,m) + rywp*vaux(248,m+1)
      vaux(322,m) = rxpa*vaux(256,m) + rxwp*vaux(256,m+1) + c20*3.0d0 
      vaux(323,m) = rxpa*vaux(257,m) + rxwp*vaux(257,m+1) + c21 
      vaux(324,m) = rzpa*vaux(249,m) + rzwp*vaux(249,m+1)
      vaux(325,m) = rzpa*vaux(250,m) + rzwp*vaux(250,m+1) + c15 
      vaux(326,m) = rzpa*vaux(251,m) + rzwp*vaux(251,m+1) + c16 
      vaux(327,m) = rxpa*vaux(261,m) + rxwp*vaux(261,m+1) + c22 
      vaux(328,m) = rypa*vaux(254,m) + rywp*vaux(254,m+1) + c17 
      vaux(329,m) = rypa*vaux(255,m) + rywp*vaux(255,m+1) + c18 
      vaux(330,m) = rypa*vaux(256,m) + rywp*vaux(256,m+1)
      vaux(331,m) = rxpa*vaux(265,m) + rxwp*vaux(265,m+1) + c23 
      vaux(332,m) = rxpa*vaux(266,m) + rxwp*vaux(266,m+1) + c24 
      vaux(333,m) = rzpa*vaux(257,m) + rzwp*vaux(257,m+1)
      vaux(334,m) = rzpa*vaux(258,m) + rzwp*vaux(258,m+1) + c19 
      vaux(335,m) = rxpa*vaux(269,m) + rxwp*vaux(269,m+1) + c27 
      vaux(336,m) = rxpa*vaux(270,m) + rxwp*vaux(270,m+1) + c28 
      vaux(337,m) = rxpa*vaux(271,m) + rxwp*vaux(271,m+1) + c29 
      vaux(338,m) = rxpa*vaux(272,m) + rxwp*vaux(272,m+1) + c30 
      vaux(339,m) = rypa*vaux(264,m) + rywp*vaux(264,m+1) + c20 
      vaux(340,m) = rypa*vaux(265,m) + rywp*vaux(265,m+1)
      vaux(341,m) = rxpa*vaux(275,m) + rxwp*vaux(275,m+1) + c33 
      vaux(342,m) = rxpa*vaux(276,m) + rxwp*vaux(276,m+1)
      vaux(343,m) = rxpa*vaux(277,m) + rxwp*vaux(277,m+1)
      vaux(344,m) = rxpa*vaux(278,m) + rxwp*vaux(278,m+1)
      vaux(345,m) = rxpa*vaux(279,m) + rxwp*vaux(279,m+1)
      vaux(346,m) = rxpa*vaux(280,m) + rxwp*vaux(280,m+1)
      vaux(347,m) = rxpa*vaux(281,m) + rxwp*vaux(281,m+1)
      vaux(348,m) = rxpa*vaux(282,m) + rxwp*vaux(282,m+1)
      vaux(349,m) = rxpa*vaux(283,m) + rxwp*vaux(283,m+1)
      vaux(350,m) = rxpa*vaux(284,m) + rxwp*vaux(284,m+1)
      vaux(351,m) = rxpa*vaux(285,m) + rxwp*vaux(285,m+1)
      vaux(352,m) = rxpa*vaux(286,m) + rxwp*vaux(286,m+1)
      vaux(353,m) = rypa*vaux(276,m) + rywp*vaux(276,m+1) + c24*10.0d0
      vaux(354,m) = rzpa*vaux(276,m) + rzwp*vaux(276,m+1)
      vaux(355,m) = rzpa*vaux(277,m) + rzwp*vaux(277,m+1) + c24 
      vaux(356,m) = rzpa*vaux(278,m) + rzwp*vaux(278,m+1) + c25 
      vaux(357,m) = rzpa*vaux(279,m) + rzwp*vaux(279,m+1) + c26 
      vaux(358,m) = rzpa*vaux(280,m) + rzwp*vaux(280,m+1) + c27*4.0d0 
      vaux(359,m) = rypa*vaux(282,m) + rywp*vaux(282,m+1) + c30*4.0d0 
      vaux(360,m) = rypa*vaux(283,m) + rywp*vaux(283,m+1) + c31 
      vaux(361,m) = rypa*vaux(284,m) + rywp*vaux(284,m+1) + c32 
      vaux(362,m) = rypa*vaux(285,m) + rywp*vaux(285,m+1) + c33 
      vaux(363,m) = rypa*vaux(286,m) + rywp*vaux(286,m+1)
      vaux(364,m) = rzpa*vaux(286,m) + rzwp*vaux(286,m+1) + c33*10.0d0
    enddo
    if( lab.eq.12 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,364,vaux,fint,li1,li2)
      else
        gama( 0) = 0.0d0
        gama( 1) = fnket*fdbra 
        gama( 2) = gama( 1) + gama(1) 
        gama( 3) = gama( 2) + gama(1) 
        gama( 4) = gama( 3) + gama(1) 
        gama( 5) = gama( 4) + gama(1) 
        gama( 6) = gama( 5) + gama(1) 
        gama( 7) = gama( 6) + gama(1) 
        gama( 8) = gama( 7) + gama(1) 
        gama( 9) = gama( 8) + gama(1) 
        gama(10) = gama( 9) + gama(1) 
        gama(11) = gama(10) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,364,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [12 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      c1 = f(1)*(vaux(221,m) - fnbra*vaux(221,m+1))
      c2 = f(2)*(vaux(222,m) - fnbra*vaux(222,m+1))
      c3 = f(2)*(vaux(223,m) - fnbra*vaux(223,m+1))
      c4 = f(1)*(vaux(224,m) - fnbra*vaux(224,m+1))
      c5 = f(3)*(vaux(226,m) - fnbra*vaux(226,m+1))
      c6 = f(1)*(vaux(227,m) - fnbra*vaux(227,m+1))
      c7 = f(1)*(vaux(230,m) - fnbra*vaux(230,m+1))
      c8 = f(1)*(vaux(231,m) - fnbra*vaux(231,m+1))
      c9 = f(2)*(vaux(232,m) - fnbra*vaux(232,m+1))
      c10= f(1)*(vaux(235,m) - fnbra*vaux(235,m+1))
      c11= f(1)*(vaux(236,m) - fnbra*vaux(236,m+1))
      c12= f(2)*(vaux(237,m) - fnbra*vaux(237,m+1))
      c13= f(2)*(vaux(240,m) - fnbra*vaux(240,m+1))
      c14= f(1)*(vaux(241,m) - fnbra*vaux(241,m+1))
      c15= f(1)*(vaux(242,m) - fnbra*vaux(242,m+1))
      c16= f(2)*(vaux(243,m) - fnbra*vaux(243,m+1))
      c17= f(3)*(vaux(244,m) - fnbra*vaux(244,m+1))
      c18= f(2)*(vaux(247,m) - fnbra*vaux(247,m+1))
      c19= f(1)*(vaux(248,m) - fnbra*vaux(248,m+1))
      c20= f(1)*(vaux(249,m) - fnbra*vaux(249,m+1))
      c21= f(2)*(vaux(250,m) - fnbra*vaux(250,m+1))
      c22= f(2)*(vaux(255,m) - fnbra*vaux(255,m+1))
      c23= f(1)*(vaux(256,m) - fnbra*vaux(256,m+1))
      c24= f(1)*(vaux(257,m) - fnbra*vaux(257,m+1))
      c25= f(1)*(vaux(265,m) - fnbra*vaux(265,m+1))
      c26= f(2)*(vaux(266,m) - fnbra*vaux(266,m+1))
      c27= f(2)*(vaux(270,m) - fnbra*vaux(270,m+1))
      c28= f(2)*(vaux(271,m) - fnbra*vaux(271,m+1))
      c29= f(2)*(vaux(275,m) - fnbra*vaux(275,m+1))
      c30= f(1)*(vaux(276,m) - fnbra*vaux(276,m+1))
      c31= f(2)*(vaux(277,m) - fnbra*vaux(277,m+1))
      c32= f(3)*(vaux(278,m) - fnbra*vaux(278,m+1))
      c33= f(1)*(vaux(279,m) - fnbra*vaux(279,m+1))
      c34= f(1)*(vaux(280,m) - fnbra*vaux(280,m+1))
      c35= f(1)*(vaux(281,m) - fnbra*vaux(281,m+1))
      c36= f(1)*(vaux(282,m) - fnbra*vaux(282,m+1))
      c37= f(1)*(vaux(283,m) - fnbra*vaux(283,m+1))
      c38= f(3)*(vaux(284,m) - fnbra*vaux(284,m+1))
      c39= f(2)*(vaux(285,m) - fnbra*vaux(285,m+1))
      c40= f(1)*(vaux(286,m) - fnbra*vaux(286,m+1))

      vaux(365,m) = rxpa*vaux(287,m) + rxwp*vaux(287,m+1) + c1*11.0d0
      vaux(366,m) = rypa*vaux(287,m) + rywp*vaux(287,m+1)
      vaux(367,m) = rzpa*vaux(287,m) + rzwp*vaux(287,m+1)
      vaux(368,m) = rypa*vaux(288,m) + rywp*vaux(288,m+1) + c1 
      vaux(369,m) = rzpa*vaux(288,m) + rzwp*vaux(288,m+1)
      vaux(370,m) = rzpa*vaux(289,m) + rzwp*vaux(289,m+1) + c1
      vaux(371,m) = rypa*vaux(290,m) + rywp*vaux(290,m+1) + c2 
      vaux(372,m) = rzpa*vaux(290,m) + rzwp*vaux(290,m+1)
      vaux(373,m) = rypa*vaux(292,m) + rywp*vaux(292,m+1)
      vaux(374,m) = rzpa*vaux(292,m) + rzwp*vaux(292,m+1) + c3 
      vaux(375,m) = rypa*vaux(293,m) + rywp*vaux(293,m+1) + c4*3.0d0 
      vaux(376,m) = rzpa*vaux(293,m) + rzwp*vaux(293,m+1)
      vaux(377,m) = rzpa*vaux(294,m) + rzwp*vaux(294,m+1) + c4 
      vaux(378,m) = rypa*vaux(296,m) + rywp*vaux(296,m+1)
      vaux(379,m) = rzpa*vaux(296,m) + rzwp*vaux(296,m+1) + c5 
      vaux(380,m) = rypa*vaux(297,m) + rywp*vaux(297,m+1) + c6*4.0d0 
      vaux(381,m) = rzpa*vaux(297,m) + rzwp*vaux(297,m+1)
      vaux(382,m) = rzpa*vaux(298,m) + rzwp*vaux(298,m+1) + c6 
      vaux(383,m) = rypa*vaux(300,m) + rywp*vaux(300,m+1) + c7 
      vaux(384,m) = rypa*vaux(301,m) + rywp*vaux(301,m+1)
      vaux(385,m) = rzpa*vaux(301,m) + rzwp*vaux(301,m+1) + c7*4.0d0 
      vaux(386,m) = rypa*vaux(302,m) + rywp*vaux(302,m+1) + c8*5.0d0 
      vaux(387,m) = rzpa*vaux(302,m) + rzwp*vaux(302,m+1)
      vaux(388,m) = rzpa*vaux(303,m) + rzwp*vaux(303,m+1) + c8 
      vaux(389,m) = rzpa*vaux(304,m) + rzwp*vaux(304,m+1) + c9 
      vaux(390,m) = rypa*vaux(306,m) + rywp*vaux(306,m+1) + c10 
      vaux(391,m) = rypa*vaux(307,m) + rywp*vaux(307,m+1)
      vaux(392,m) = rzpa*vaux(307,m) + rzwp*vaux(307,m+1) + c10*5.0d0 
      vaux(393,m) = rxpa*vaux(315,m) + rxwp*vaux(315,m+1) + c20*4.0d0 
      vaux(394,m) = rzpa*vaux(308,m) + rzwp*vaux(308,m+1)
      vaux(395,m) = rzpa*vaux(309,m) + rzwp*vaux(309,m+1) + c11 
      vaux(396,m) = rzpa*vaux(310,m) + rzwp*vaux(310,m+1) + c12 
      vaux(397,m) = rypa*vaux(312,m) + rywp*vaux(312,m+1) + c13 
      vaux(398,m) = rypa*vaux(313,m) + rywp*vaux(313,m+1) + c14 
      vaux(399,m) = rypa*vaux(314,m) + rywp*vaux(314,m+1)
      vaux(400,m) = rxpa*vaux(322,m) + rxwp*vaux(322,m+1) + c23*4.0d0 
      vaux(401,m) = rxpa*vaux(323,m) + rxwp*vaux(323,m+1) + c24*3.0d0 
      vaux(402,m) = rzpa*vaux(315,m) + rzwp*vaux(315,m+1)
      vaux(403,m) = rzpa*vaux(316,m) + rzwp*vaux(316,m+1) + c15 
      vaux(404,m) = rzpa*vaux(317,m) + rzwp*vaux(317,m+1) + c16 
      vaux(405,m) = rzpa*vaux(318,m) + rzwp*vaux(318,m+1) + c17 
      vaux(406,m) = rypa*vaux(320,m) + rywp*vaux(320,m+1) + c18 
      vaux(407,m) = rypa*vaux(321,m) + rywp*vaux(321,m+1) + c19 
      vaux(408,m) = rypa*vaux(322,m) + rywp*vaux(322,m+1)
      vaux(409,m) = rxpa*vaux(331,m) + rxwp*vaux(331,m+1) + c25*3.0d0 
      vaux(410,m) = rxpa*vaux(332,m) + rxwp*vaux(332,m+1) + c26 
      vaux(411,m) = rzpa*vaux(323,m) + rzwp*vaux(323,m+1)
      vaux(412,m) = rzpa*vaux(324,m) + rzwp*vaux(324,m+1) + c20 
      vaux(413,m) = rzpa*vaux(325,m) + rzwp*vaux(325,m+1) + c21 
      vaux(414,m) = rxpa*vaux(336,m) + rxwp*vaux(336,m+1) + c27 
      vaux(415,m) = rxpa*vaux(337,m) + rxwp*vaux(337,m+1) + c28 
      vaux(416,m) = rypa*vaux(329,m) + rywp*vaux(329,m+1) + c22 
      vaux(417,m) = rypa*vaux(330,m) + rywp*vaux(330,m+1) + c23 
      vaux(418,m) = rypa*vaux(331,m) + rywp*vaux(331,m+1)
      vaux(419,m) = rxpa*vaux(341,m) + rxwp*vaux(341,m+1) + c29 
      vaux(420,m) = rxpa*vaux(342,m) + rxwp*vaux(342,m+1) + c30 
      vaux(421,m) = rzpa*vaux(332,m) + rzwp*vaux(332,m+1)
      vaux(422,m) = rzpa*vaux(333,m) + rzwp*vaux(333,m+1) + c24 
      vaux(423,m) = rxpa*vaux(345,m) + rxwp*vaux(345,m+1) + c33 
      vaux(424,m) = rxpa*vaux(346,m) + rxwp*vaux(346,m+1) + c34 
      vaux(425,m) = rxpa*vaux(347,m) + rxwp*vaux(347,m+1) + c35 
      vaux(426,m) = rxpa*vaux(348,m) + rxwp*vaux(348,m+1) + c36 
      vaux(427,m) = rxpa*vaux(349,m) + rxwp*vaux(349,m+1) + c37 
      vaux(428,m) = rypa*vaux(340,m) + rywp*vaux(340,m+1) + c25 
      vaux(429,m) = rypa*vaux(341,m) + rywp*vaux(341,m+1)
      vaux(430,m) = rxpa*vaux(352,m) + rxwp*vaux(352,m+1) + c40 
      vaux(431,m) = rxpa*vaux(353,m) + rxwp*vaux(353,m+1)
      vaux(432,m) = rxpa*vaux(354,m) + rxwp*vaux(354,m+1)
      vaux(433,m) = rxpa*vaux(355,m) + rxwp*vaux(355,m+1)
      vaux(434,m) = rxpa*vaux(356,m) + rxwp*vaux(356,m+1)
      vaux(435,m) = rxpa*vaux(357,m) + rxwp*vaux(357,m+1)
      vaux(436,m) = rxpa*vaux(358,m) + rxwp*vaux(358,m+1)
      vaux(437,m) = rxpa*vaux(359,m) + rxwp*vaux(359,m+1)
      vaux(438,m) = rxpa*vaux(360,m) + rxwp*vaux(360,m+1)
      vaux(439,m) = rxpa*vaux(361,m) + rxwp*vaux(361,m+1)
      vaux(440,m) = rxpa*vaux(362,m) + rxwp*vaux(362,m+1)
      vaux(441,m) = rxpa*vaux(363,m) + rxwp*vaux(363,m+1)
      vaux(442,m) = rxpa*vaux(364,m) + rxwp*vaux(364,m+1)
      vaux(443,m) = rypa*vaux(353,m) + rywp*vaux(353,m+1) + c30*11.0d0
      vaux(444,m) = rzpa*vaux(353,m) + rzwp*vaux(353,m+1)
      vaux(445,m) = rzpa*vaux(354,m) + rzwp*vaux(354,m+1) + c30 
      vaux(446,m) = rzpa*vaux(355,m) + rzwp*vaux(355,m+1) + c31 
      vaux(447,m) = rzpa*vaux(356,m) + rzwp*vaux(356,m+1) + c32 
      vaux(448,m) = rzpa*vaux(357,m) + rzwp*vaux(357,m+1) + c33*4.0d0 
      vaux(449,m) = rzpa*vaux(358,m) + rzwp*vaux(358,m+1) + c34*5.0d0 
      vaux(450,m) = rypa*vaux(360,m) + rywp*vaux(360,m+1) + c37*4.0d0 
      vaux(451,m) = rypa*vaux(361,m) + rywp*vaux(361,m+1) + c38 
      vaux(452,m) = rypa*vaux(362,m) + rywp*vaux(362,m+1) + c39 
      vaux(453,m) = rypa*vaux(363,m) + rywp*vaux(363,m+1) + c40 
      vaux(454,m) = rypa*vaux(364,m) + rywp*vaux(364,m+1)
      vaux(455,m) = rzpa*vaux(364,m) + rzwp*vaux(364,m+1) + c40*11.0d0
    enddo
    if( lab.eq.13 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,455,vaux,fint,li1,li2)
      else
        gama( 0) = 0.0d0
        gama( 1) = fnket*fdbra 
        gama( 2) = gama( 1) + gama(1) 
        gama( 3) = gama( 2) + gama(1) 
        gama( 4) = gama( 3) + gama(1) 
        gama( 5) = gama( 4) + gama(1) 
        gama( 6) = gama( 5) + gama(1) 
        gama( 7) = gama( 6) + gama(1) 
        gama( 8) = gama( 7) + gama(1) 
        gama( 9) = gama( 8) + gama(1) 
        gama(10) = gama( 9) + gama(1) 
        gama(11) = gama(10) + gama(1) 
        gama(12) = gama(11) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,455,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [13 | * ] --!
    f( 1) = fdbra
    f( 2) = f( 1) + f(1) 
    f( 3) = f( 2) + f(1) 
    f( 4) = f( 3) + f(1) 
    f( 5) = f( 4) + f(1) 
    f( 6) = f( 5) + f(1) 
    f( 7) = f( 6) + f(1) 
    f( 8) = f( 7) + f(1) 
    f( 9) = f( 8) + f(1) 
    f(10) = f( 9) + f(1) 
    f(11) = f(10) + f(1) 
    f(12) = f(11) + f(1) 
    mmax=mmax-1
    do m=1,mmax
      vaux(456,m) = rxpa*vaux(365,m) + rxwp*vaux(365,m+1) + f( 12)*(vaux(287,m) - fnbra*vaux(287,m+1))
      vaux(457,m) = rypa*vaux(365,m) + rywp*vaux(365,m+1)
      vaux(458,m) = rzpa*vaux(365,m) + rzwp*vaux(365,m+1)
      vaux(459,m) = rypa*vaux(366,m) + rywp*vaux(366,m+1) + f(  1)*(vaux(287,m) - fnbra*vaux(287,m+1))
      vaux(460,m) = rzpa*vaux(366,m) + rzwp*vaux(366,m+1)
      vaux(461,m) = rzpa*vaux(367,m) + rzwp*vaux(367,m+1) + f(  1)*(vaux(287,m) - fnbra*vaux(287,m+1))
      vaux(462,m) = rypa*vaux(368,m) + rywp*vaux(368,m+1) + f(  2)*(vaux(288,m) - fnbra*vaux(288,m+1))
      vaux(463,m) = rzpa*vaux(368,m) + rzwp*vaux(368,m+1)
      vaux(464,m) = rypa*vaux(370,m) + rywp*vaux(370,m+1)
      vaux(465,m) = rzpa*vaux(370,m) + rzwp*vaux(370,m+1) + f(  2)*(vaux(289,m) - fnbra*vaux(289,m+1))
      vaux(466,m) = rypa*vaux(371,m) + rywp*vaux(371,m+1) + f(  3)*(vaux(290,m) - fnbra*vaux(290,m+1))
      vaux(467,m) = rzpa*vaux(371,m) + rzwp*vaux(371,m+1)
      vaux(468,m) = rzpa*vaux(372,m) + rzwp*vaux(372,m+1) + f(  1)*(vaux(290,m) - fnbra*vaux(290,m+1))
      vaux(469,m) = rypa*vaux(374,m) + rywp*vaux(374,m+1)
      vaux(470,m) = rzpa*vaux(374,m) + rzwp*vaux(374,m+1) + f(  3)*(vaux(292,m) - fnbra*vaux(292,m+1))
      vaux(471,m) = rypa*vaux(375,m) + rywp*vaux(375,m+1) + f(  4)*(vaux(293,m) - fnbra*vaux(293,m+1))
      vaux(472,m) = rzpa*vaux(375,m) + rzwp*vaux(375,m+1)
      vaux(473,m) = rzpa*vaux(376,m) + rzwp*vaux(376,m+1) + f(  1)*(vaux(293,m) - fnbra*vaux(293,m+1))
      vaux(474,m) = rypa*vaux(378,m) + rywp*vaux(378,m+1) + f(  1)*(vaux(296,m) - fnbra*vaux(296,m+1))
      vaux(475,m) = rypa*vaux(379,m) + rywp*vaux(379,m+1)
      vaux(476,m) = rzpa*vaux(379,m) + rzwp*vaux(379,m+1) + f(  4)*(vaux(296,m) - fnbra*vaux(296,m+1))
      vaux(477,m) = rypa*vaux(380,m) + rywp*vaux(380,m+1) + f(  5)*(vaux(297,m) - fnbra*vaux(297,m+1))
      vaux(478,m) = rzpa*vaux(380,m) + rzwp*vaux(380,m+1)
      vaux(479,m) = rzpa*vaux(381,m) + rzwp*vaux(381,m+1) + f(  1)*(vaux(297,m) - fnbra*vaux(297,m+1))
      vaux(480,m) = rzpa*vaux(382,m) + rzwp*vaux(382,m+1) + f(  2)*(vaux(298,m) - fnbra*vaux(298,m+1))
      vaux(481,m) = rypa*vaux(384,m) + rywp*vaux(384,m+1) + f(  1)*(vaux(301,m) - fnbra*vaux(301,m+1))
      vaux(482,m) = rypa*vaux(385,m) + rywp*vaux(385,m+1)
      vaux(483,m) = rzpa*vaux(385,m) + rzwp*vaux(385,m+1) + f(  5)*(vaux(301,m) - fnbra*vaux(301,m+1))
      vaux(484,m) = rxpa*vaux(393,m) + rxwp*vaux(393,m+1) + f(  5)*(vaux(315,m) - fnbra*vaux(315,m+1))
      vaux(485,m) = rzpa*vaux(386,m) + rzwp*vaux(386,m+1)
      vaux(486,m) = rzpa*vaux(387,m) + rzwp*vaux(387,m+1) + f(  1)*(vaux(302,m) - fnbra*vaux(302,m+1))
      vaux(487,m) = rzpa*vaux(388,m) + rzwp*vaux(388,m+1) + f(  2)*(vaux(303,m) - fnbra*vaux(303,m+1))
      vaux(488,m) = rypa*vaux(390,m) + rywp*vaux(390,m+1) + f(  2)*(vaux(306,m) - fnbra*vaux(306,m+1))
      vaux(489,m) = rypa*vaux(391,m) + rywp*vaux(391,m+1) + f(  1)*(vaux(307,m) - fnbra*vaux(307,m+1))
      vaux(490,m) = rypa*vaux(392,m) + rywp*vaux(392,m+1)
      vaux(491,m) = rxpa*vaux(400,m) + rxwp*vaux(400,m+1) + f(  5)*(vaux(322,m) - fnbra*vaux(322,m+1))
      vaux(492,m) = rxpa*vaux(401,m) + rxwp*vaux(401,m+1) + f(  4)*(vaux(323,m) - fnbra*vaux(323,m+1))
      vaux(493,m) = rzpa*vaux(393,m) + rzwp*vaux(393,m+1)
      vaux(494,m) = rzpa*vaux(394,m) + rzwp*vaux(394,m+1) + f(  1)*(vaux(308,m) - fnbra*vaux(308,m+1))
      vaux(495,m) = rzpa*vaux(395,m) + rzwp*vaux(395,m+1) + f(  2)*(vaux(309,m) - fnbra*vaux(309,m+1))
      vaux(496,m) = rzpa*vaux(396,m) + rzwp*vaux(396,m+1) + f(  3)*(vaux(310,m) - fnbra*vaux(310,m+1))
      vaux(497,m) = rypa*vaux(398,m) + rywp*vaux(398,m+1) + f(  2)*(vaux(313,m) - fnbra*vaux(313,m+1))
      vaux(498,m) = rypa*vaux(399,m) + rywp*vaux(399,m+1) + f(  1)*(vaux(314,m) - fnbra*vaux(314,m+1))
      vaux(499,m) = rypa*vaux(400,m) + rywp*vaux(400,m+1)
      vaux(500,m) = rxpa*vaux(409,m) + rxwp*vaux(409,m+1) + f(  4)*(vaux(331,m) - fnbra*vaux(331,m+1))
      vaux(501,m) = rxpa*vaux(410,m) + rxwp*vaux(410,m+1) + f(  3)*(vaux(332,m) - fnbra*vaux(332,m+1))
      vaux(502,m) = rzpa*vaux(401,m) + rzwp*vaux(401,m+1)
      vaux(503,m) = rzpa*vaux(402,m) + rzwp*vaux(402,m+1) + f(  1)*(vaux(315,m) - fnbra*vaux(315,m+1))
      vaux(504,m) = rzpa*vaux(403,m) + rzwp*vaux(403,m+1) + f(  2)*(vaux(316,m) - fnbra*vaux(316,m+1))
      vaux(505,m) = rzpa*vaux(404,m) + rzwp*vaux(404,m+1) + f(  3)*(vaux(317,m) - fnbra*vaux(317,m+1))
      vaux(506,m) = rypa*vaux(406,m) + rywp*vaux(406,m+1) + f(  3)*(vaux(320,m) - fnbra*vaux(320,m+1))
      vaux(507,m) = rypa*vaux(407,m) + rywp*vaux(407,m+1) + f(  2)*(vaux(321,m) - fnbra*vaux(321,m+1))
      vaux(508,m) = rypa*vaux(408,m) + rywp*vaux(408,m+1) + f(  1)*(vaux(322,m) - fnbra*vaux(322,m+1))
      vaux(509,m) = rypa*vaux(409,m) + rywp*vaux(409,m+1)
      vaux(510,m) = rxpa*vaux(419,m) + rxwp*vaux(419,m+1) + f(  3)*(vaux(341,m) - fnbra*vaux(341,m+1))
      vaux(511,m) = rxpa*vaux(420,m) + rxwp*vaux(420,m+1) + f(  2)*(vaux(342,m) - fnbra*vaux(342,m+1))
      vaux(512,m) = rzpa*vaux(410,m) + rzwp*vaux(410,m+1)
      vaux(513,m) = rzpa*vaux(411,m) + rzwp*vaux(411,m+1) + f(  1)*(vaux(323,m) - fnbra*vaux(323,m+1))
      vaux(514,m) = rzpa*vaux(412,m) + rzwp*vaux(412,m+1) + f(  2)*(vaux(324,m) - fnbra*vaux(324,m+1))
      vaux(515,m) = rxpa*vaux(424,m) + rxwp*vaux(424,m+1) + f(  2)*(vaux(346,m) - fnbra*vaux(346,m+1))
      vaux(516,m) = rxpa*vaux(425,m) + rxwp*vaux(425,m+1) + f(  2)*(vaux(347,m) - fnbra*vaux(347,m+1))
      vaux(517,m) = rxpa*vaux(426,m) + rxwp*vaux(426,m+1) + f(  2)*(vaux(348,m) - fnbra*vaux(348,m+1))
      vaux(518,m) = rypa*vaux(417,m) + rywp*vaux(417,m+1) + f(  2)*(vaux(330,m) - fnbra*vaux(330,m+1))
      vaux(519,m) = rypa*vaux(418,m) + rywp*vaux(418,m+1) + f(  1)*(vaux(331,m) - fnbra*vaux(331,m+1))
      vaux(520,m) = rypa*vaux(419,m) + rywp*vaux(419,m+1)
      vaux(521,m) = rxpa*vaux(430,m) + rxwp*vaux(430,m+1) + f(  2)*(vaux(352,m) - fnbra*vaux(352,m+1))
      vaux(522,m) = rxpa*vaux(431,m) + rxwp*vaux(431,m+1) + f(  1)*(vaux(353,m) - fnbra*vaux(353,m+1))
      vaux(523,m) = rzpa*vaux(420,m) + rzwp*vaux(420,m+1)
      vaux(524,m) = rzpa*vaux(421,m) + rzwp*vaux(421,m+1) + f(  1)*(vaux(332,m) - fnbra*vaux(332,m+1))
      vaux(525,m) = rxpa*vaux(434,m) + rxwp*vaux(434,m+1) + f(  1)*(vaux(356,m) - fnbra*vaux(356,m+1))
      vaux(526,m) = rxpa*vaux(435,m) + rxwp*vaux(435,m+1) + f(  1)*(vaux(357,m) - fnbra*vaux(357,m+1))
      vaux(527,m) = rxpa*vaux(436,m) + rxwp*vaux(436,m+1) + f(  1)*(vaux(358,m) - fnbra*vaux(358,m+1))
      vaux(528,m) = rxpa*vaux(437,m) + rxwp*vaux(437,m+1) + f(  1)*(vaux(359,m) - fnbra*vaux(359,m+1))
      vaux(529,m) = rxpa*vaux(438,m) + rxwp*vaux(438,m+1) + f(  1)*(vaux(360,m) - fnbra*vaux(360,m+1))
      vaux(530,m) = rxpa*vaux(439,m) + rxwp*vaux(439,m+1) + f(  1)*(vaux(361,m) - fnbra*vaux(361,m+1))
      vaux(531,m) = rypa*vaux(429,m) + rywp*vaux(429,m+1) + f(  1)*(vaux(341,m) - fnbra*vaux(341,m+1))
      vaux(532,m) = rypa*vaux(430,m) + rywp*vaux(430,m+1)
      vaux(533,m) = rxpa*vaux(442,m) + rxwp*vaux(442,m+1) + f(  1)*(vaux(364,m) - fnbra*vaux(364,m+1))
      vaux(534,m) = rxpa*vaux(443,m) + rxwp*vaux(443,m+1)
      vaux(535,m) = rxpa*vaux(444,m) + rxwp*vaux(444,m+1)
      vaux(536,m) = rxpa*vaux(445,m) + rxwp*vaux(445,m+1)
      vaux(537,m) = rxpa*vaux(446,m) + rxwp*vaux(446,m+1)
      vaux(538,m) = rxpa*vaux(447,m) + rxwp*vaux(447,m+1)
      vaux(539,m) = rxpa*vaux(448,m) + rxwp*vaux(448,m+1)
      vaux(540,m) = rxpa*vaux(449,m) + rxwp*vaux(449,m+1)
      vaux(541,m) = rxpa*vaux(450,m) + rxwp*vaux(450,m+1)
      vaux(542,m) = rxpa*vaux(451,m) + rxwp*vaux(451,m+1)
      vaux(543,m) = rxpa*vaux(452,m) + rxwp*vaux(452,m+1)
      vaux(544,m) = rxpa*vaux(453,m) + rxwp*vaux(453,m+1)
      vaux(545,m) = rxpa*vaux(454,m) + rxwp*vaux(454,m+1)
      vaux(546,m) = rxpa*vaux(455,m) + rxwp*vaux(455,m+1)
      vaux(547,m) = rypa*vaux(443,m) + rywp*vaux(443,m+1) + f( 12)*(vaux(353,m) - fnbra*vaux(353,m+1))
      vaux(548,m) = rzpa*vaux(443,m) + rzwp*vaux(443,m+1)
      vaux(549,m) = rzpa*vaux(444,m) + rzwp*vaux(444,m+1) + f(  1)*(vaux(353,m) - fnbra*vaux(353,m+1))
      vaux(550,m) = rzpa*vaux(445,m) + rzwp*vaux(445,m+1) + f(  2)*(vaux(354,m) - fnbra*vaux(354,m+1))
      vaux(551,m) = rzpa*vaux(446,m) + rzwp*vaux(446,m+1) + f(  3)*(vaux(355,m) - fnbra*vaux(355,m+1))
      vaux(552,m) = rzpa*vaux(447,m) + rzwp*vaux(447,m+1) + f(  4)*(vaux(356,m) - fnbra*vaux(356,m+1))
      vaux(553,m) = rzpa*vaux(448,m) + rzwp*vaux(448,m+1) + f(  5)*(vaux(357,m) - fnbra*vaux(357,m+1))
      vaux(554,m) = rypa*vaux(450,m) + rywp*vaux(450,m+1) + f(  5)*(vaux(360,m) - fnbra*vaux(360,m+1))
      vaux(555,m) = rypa*vaux(451,m) + rywp*vaux(451,m+1) + f(  4)*(vaux(361,m) - fnbra*vaux(361,m+1))
      vaux(556,m) = rypa*vaux(452,m) + rywp*vaux(452,m+1) + f(  3)*(vaux(362,m) - fnbra*vaux(362,m+1))
      vaux(557,m) = rypa*vaux(453,m) + rywp*vaux(453,m+1) + f(  2)*(vaux(363,m) - fnbra*vaux(363,m+1))
      vaux(558,m) = rypa*vaux(454,m) + rywp*vaux(454,m+1) + f(  1)*(vaux(364,m) - fnbra*vaux(364,m+1))
      vaux(559,m) = rypa*vaux(455,m) + rywp*vaux(455,m+1)
      vaux(560,m) = rzpa*vaux(455,m) + rzwp*vaux(455,m+1) + f( 12)*(vaux(364,m) - fnbra*vaux(364,m+1))
    enddo
    if( lab.eq.14 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,560,vaux,fint,li1,li2)
      else
        gama( 0) = 0.0d0
        gama( 1) = fnket*fdbra 
        gama( 2) = gama( 1) + gama(1) 
        gama( 3) = gama( 2) + gama(1) 
        gama( 4) = gama( 3) + gama(1) 
        gama( 5) = gama( 4) + gama(1) 
        gama( 6) = gama( 5) + gama(1) 
        gama( 7) = gama( 6) + gama(1) 
        gama( 8) = gama( 7) + gama(1) 
        gama( 9) = gama( 8) + gama(1) 
        gama(10) = gama( 9) + gama(1) 
        gama(11) = gama(10) + gama(1) 
        gama(12) = gama(11) + gama(1) 
        gama(13) = gama(12) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,560,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    !-- [14 | * ] --!
    f(13) = f(12) + f(1) 
    mmax=mmax-1
    do m=1,mmax
       vaux(561,m) = rxpa*vaux(456,m) + rxwp*vaux(456,m+1) + f( 13)*(vaux(365,m) - fnbra*vaux(365,m+1))
       vaux(562,m) = rypa*vaux(456,m) + rywp*vaux(456,m+1)
       vaux(563,m) = rzpa*vaux(456,m) + rzwp*vaux(456,m+1)
       vaux(564,m) = rypa*vaux(457,m) + rywp*vaux(457,m+1) + f(  1)*(vaux(365,m) - fnbra*vaux(365,m+1))
       vaux(565,m) = rzpa*vaux(457,m) + rzwp*vaux(457,m+1)
       vaux(566,m) = rzpa*vaux(458,m) + rzwp*vaux(458,m+1) + f(  1)*(vaux(365,m) - fnbra*vaux(365,m+1))
       vaux(567,m) = rypa*vaux(459,m) + rywp*vaux(459,m+1) + f(  2)*(vaux(366,m) - fnbra*vaux(366,m+1))
       vaux(568,m) = rzpa*vaux(459,m) + rzwp*vaux(459,m+1)
       vaux(569,m) = rypa*vaux(461,m) + rywp*vaux(461,m+1)
       vaux(570,m) = rzpa*vaux(461,m) + rzwp*vaux(461,m+1) + f(  2)*(vaux(367,m) - fnbra*vaux(367,m+1))
       vaux(571,m) = rypa*vaux(462,m) + rywp*vaux(462,m+1) + f(  3)*(vaux(368,m) - fnbra*vaux(368,m+1))
       vaux(572,m) = rzpa*vaux(462,m) + rzwp*vaux(462,m+1)
       vaux(573,m) = rzpa*vaux(463,m) + rzwp*vaux(463,m+1) + f(  1)*(vaux(368,m) - fnbra*vaux(368,m+1))
       vaux(574,m) = rypa*vaux(465,m) + rywp*vaux(465,m+1)
       vaux(575,m) = rzpa*vaux(465,m) + rzwp*vaux(465,m+1) + f(  3)*(vaux(370,m) - fnbra*vaux(370,m+1))
       vaux(576,m) = rypa*vaux(466,m) + rywp*vaux(466,m+1) + f(  4)*(vaux(371,m) - fnbra*vaux(371,m+1))
       vaux(577,m) = rzpa*vaux(466,m) + rzwp*vaux(466,m+1)
       vaux(578,m) = rzpa*vaux(467,m) + rzwp*vaux(467,m+1) + f(  1)*(vaux(371,m) - fnbra*vaux(371,m+1))
       vaux(579,m) = rypa*vaux(469,m) + rywp*vaux(469,m+1) + f(  1)*(vaux(374,m) - fnbra*vaux(374,m+1))
       vaux(580,m) = rypa*vaux(470,m) + rywp*vaux(470,m+1)
       vaux(581,m) = rzpa*vaux(470,m) + rzwp*vaux(470,m+1) + f(  4)*(vaux(374,m) - fnbra*vaux(374,m+1))
       vaux(582,m) = rypa*vaux(471,m) + rywp*vaux(471,m+1) + f(  5)*(vaux(375,m) - fnbra*vaux(375,m+1))
       vaux(583,m) = rzpa*vaux(471,m) + rzwp*vaux(471,m+1)
       vaux(584,m) = rzpa*vaux(472,m) + rzwp*vaux(472,m+1) + f(  1)*(vaux(375,m) - fnbra*vaux(375,m+1))
       vaux(585,m) = rzpa*vaux(473,m) + rzwp*vaux(473,m+1) + f(  2)*(vaux(376,m) - fnbra*vaux(376,m+1))
       vaux(586,m) = rypa*vaux(475,m) + rywp*vaux(475,m+1) + f(  1)*(vaux(379,m) - fnbra*vaux(379,m+1))
       vaux(587,m) = rypa*vaux(476,m) + rywp*vaux(476,m+1)
       vaux(588,m) = rzpa*vaux(476,m) + rzwp*vaux(476,m+1) + f(  5)*(vaux(379,m) - fnbra*vaux(379,m+1))
       vaux(589,m) = rypa*vaux(477,m) + rywp*vaux(477,m+1) + f(  6)*(vaux(380,m) - fnbra*vaux(380,m+1))
       vaux(590,m) = rzpa*vaux(477,m) + rzwp*vaux(477,m+1)
       vaux(591,m) = rzpa*vaux(478,m) + rzwp*vaux(478,m+1) + f(  1)*(vaux(380,m) - fnbra*vaux(380,m+1))
       vaux(592,m) = rzpa*vaux(479,m) + rzwp*vaux(479,m+1) + f(  2)*(vaux(381,m) - fnbra*vaux(381,m+1))
       vaux(593,m) = rypa*vaux(481,m) + rywp*vaux(481,m+1) + f(  2)*(vaux(384,m) - fnbra*vaux(384,m+1))
       vaux(594,m) = rypa*vaux(482,m) + rywp*vaux(482,m+1) + f(  1)*(vaux(385,m) - fnbra*vaux(385,m+1))
       vaux(595,m) = rypa*vaux(483,m) + rywp*vaux(483,m+1)
       vaux(596,m) = rzpa*vaux(483,m) + rzwp*vaux(483,m+1) + f(  6)*(vaux(385,m) - fnbra*vaux(385,m+1))
       vaux(597,m) = rxpa*vaux(492,m) + rxwp*vaux(492,m+1) + f(  5)*(vaux(401,m) - fnbra*vaux(401,m+1))
       vaux(598,m) = rzpa*vaux(484,m) + rzwp*vaux(484,m+1)
       vaux(599,m) = rzpa*vaux(485,m) + rzwp*vaux(485,m+1) + f(  1)*(vaux(386,m) - fnbra*vaux(386,m+1))
       vaux(600,m) = rzpa*vaux(486,m) + rzwp*vaux(486,m+1) + f(  2)*(vaux(387,m) - fnbra*vaux(387,m+1))
       vaux(601,m) = rzpa*vaux(487,m) + rzwp*vaux(487,m+1) + f(  3)*(vaux(388,m) - fnbra*vaux(388,m+1))
       vaux(602,m) = rypa*vaux(489,m) + rywp*vaux(489,m+1) + f(  2)*(vaux(391,m) - fnbra*vaux(391,m+1))
       vaux(603,m) = rypa*vaux(490,m) + rywp*vaux(490,m+1) + f(  1)*(vaux(392,m) - fnbra*vaux(392,m+1))
       vaux(604,m) = rypa*vaux(491,m) + rywp*vaux(491,m+1)
       vaux(605,m) = rxpa*vaux(500,m) + rxwp*vaux(500,m+1) + f(  5)*(vaux(409,m) - fnbra*vaux(409,m+1))
       vaux(606,m) = rxpa*vaux(501,m) + rxwp*vaux(501,m+1) + f(  4)*(vaux(410,m) - fnbra*vaux(410,m+1))
       vaux(607,m) = rzpa*vaux(492,m) + rzwp*vaux(492,m+1)
       vaux(608,m) = rzpa*vaux(493,m) + rzwp*vaux(493,m+1) + f(  1)*(vaux(393,m) - fnbra*vaux(393,m+1))
       vaux(609,m) = rzpa*vaux(494,m) + rzwp*vaux(494,m+1) + f(  2)*(vaux(394,m) - fnbra*vaux(394,m+1))
       vaux(610,m) = rzpa*vaux(495,m) + rzwp*vaux(495,m+1) + f(  3)*(vaux(395,m) - fnbra*vaux(395,m+1))
       vaux(611,m) = rypa*vaux(497,m) + rywp*vaux(497,m+1) + f(  3)*(vaux(398,m) - fnbra*vaux(398,m+1))
       vaux(612,m) = rypa*vaux(498,m) + rywp*vaux(498,m+1) + f(  2)*(vaux(399,m) - fnbra*vaux(399,m+1))
       vaux(613,m) = rypa*vaux(499,m) + rywp*vaux(499,m+1) + f(  1)*(vaux(400,m) - fnbra*vaux(400,m+1))
       vaux(614,m) = rypa*vaux(500,m) + rywp*vaux(500,m+1)
       vaux(615,m) = rxpa*vaux(510,m) + rxwp*vaux(510,m+1) + f(  4)*(vaux(419,m) - fnbra*vaux(419,m+1))
       vaux(616,m) = rxpa*vaux(511,m) + rxwp*vaux(511,m+1) + f(  3)*(vaux(420,m) - fnbra*vaux(420,m+1))
       vaux(617,m) = rzpa*vaux(501,m) + rzwp*vaux(501,m+1)
       vaux(618,m) = rzpa*vaux(502,m) + rzwp*vaux(502,m+1) + f(  1)*(vaux(401,m) - fnbra*vaux(401,m+1))
       vaux(619,m) = rzpa*vaux(503,m) + rzwp*vaux(503,m+1) + f(  2)*(vaux(402,m) - fnbra*vaux(402,m+1))
       vaux(620,m) = rzpa*vaux(504,m) + rzwp*vaux(504,m+1) + f(  3)*(vaux(403,m) - fnbra*vaux(403,m+1))
       vaux(621,m) = rxpa*vaux(516,m) + rxwp*vaux(516,m+1) + f(  3)*(vaux(425,m) - fnbra*vaux(425,m+1))
       vaux(622,m) = rypa*vaux(507,m) + rywp*vaux(507,m+1) + f(  3)*(vaux(407,m) - fnbra*vaux(407,m+1))
       vaux(623,m) = rypa*vaux(508,m) + rywp*vaux(508,m+1) + f(  2)*(vaux(408,m) - fnbra*vaux(408,m+1))
       vaux(624,m) = rypa*vaux(509,m) + rywp*vaux(509,m+1) + f(  1)*(vaux(409,m) - fnbra*vaux(409,m+1))
       vaux(625,m) = rypa*vaux(510,m) + rywp*vaux(510,m+1)
       vaux(626,m) = rxpa*vaux(521,m) + rxwp*vaux(521,m+1) + f(  3)*(vaux(430,m) - fnbra*vaux(430,m+1))
       vaux(627,m) = rxpa*vaux(522,m) + rxwp*vaux(522,m+1) + f(  2)*(vaux(431,m) - fnbra*vaux(431,m+1))
       vaux(628,m) = rzpa*vaux(511,m) + rzwp*vaux(511,m+1)
       vaux(629,m) = rzpa*vaux(512,m) + rzwp*vaux(512,m+1) + f(  1)*(vaux(410,m) - fnbra*vaux(410,m+1))
       vaux(630,m) = rzpa*vaux(513,m) + rzwp*vaux(513,m+1) + f(  2)*(vaux(411,m) - fnbra*vaux(411,m+1))
       vaux(631,m) = rxpa*vaux(526,m) + rxwp*vaux(526,m+1) + f(  2)*(vaux(435,m) - fnbra*vaux(435,m+1))
       vaux(632,m) = rxpa*vaux(527,m) + rxwp*vaux(527,m+1) + f(  2)*(vaux(436,m) - fnbra*vaux(436,m+1))
       vaux(633,m) = rxpa*vaux(528,m) + rxwp*vaux(528,m+1) + f(  2)*(vaux(437,m) - fnbra*vaux(437,m+1))
       vaux(634,m) = rxpa*vaux(529,m) + rxwp*vaux(529,m+1) + f(  2)*(vaux(438,m) - fnbra*vaux(438,m+1))
       vaux(635,m) = rypa*vaux(519,m) + rywp*vaux(519,m+1) + f(  2)*(vaux(418,m) - fnbra*vaux(418,m+1))
       vaux(636,m) = rypa*vaux(520,m) + rywp*vaux(520,m+1) + f(  1)*(vaux(419,m) - fnbra*vaux(419,m+1))
       vaux(637,m) = rypa*vaux(521,m) + rywp*vaux(521,m+1)
       vaux(638,m) = rxpa*vaux(533,m) + rxwp*vaux(533,m+1) + f(  2)*(vaux(442,m) - fnbra*vaux(442,m+1))
       vaux(639,m) = rxpa*vaux(534,m) + rxwp*vaux(534,m+1) + f(  1)*(vaux(443,m) - fnbra*vaux(443,m+1))
       vaux(640,m) = rzpa*vaux(522,m) + rzwp*vaux(522,m+1)
       vaux(641,m) = rzpa*vaux(523,m) + rzwp*vaux(523,m+1) + f(  1)*(vaux(420,m) - fnbra*vaux(420,m+1))
       vaux(642,m) = rxpa*vaux(537,m) + rxwp*vaux(537,m+1) + f(  1)*(vaux(446,m) - fnbra*vaux(446,m+1))
       vaux(643,m) = rxpa*vaux(538,m) + rxwp*vaux(538,m+1) + f(  1)*(vaux(447,m) - fnbra*vaux(447,m+1))
       vaux(644,m) = rxpa*vaux(539,m) + rxwp*vaux(539,m+1) + f(  1)*(vaux(448,m) - fnbra*vaux(448,m+1))
       vaux(645,m) = rxpa*vaux(540,m) + rxwp*vaux(540,m+1) + f(  1)*(vaux(449,m) - fnbra*vaux(449,m+1))
       vaux(646,m) = rxpa*vaux(541,m) + rxwp*vaux(541,m+1) + f(  1)*(vaux(450,m) - fnbra*vaux(450,m+1))
       vaux(647,m) = rxpa*vaux(542,m) + rxwp*vaux(542,m+1) + f(  1)*(vaux(451,m) - fnbra*vaux(451,m+1))
       vaux(648,m) = rxpa*vaux(543,m) + rxwp*vaux(543,m+1) + f(  1)*(vaux(452,m) - fnbra*vaux(452,m+1))
       vaux(649,m) = rypa*vaux(532,m) + rywp*vaux(532,m+1) + f(  1)*(vaux(430,m) - fnbra*vaux(430,m+1))
       vaux(650,m) = rypa*vaux(533,m) + rywp*vaux(533,m+1)
       vaux(651,m) = rxpa*vaux(546,m) + rxwp*vaux(546,m+1) + f(  1)*(vaux(455,m) - fnbra*vaux(455,m+1))
       vaux(652,m) = rxpa*vaux(547,m) + rxwp*vaux(547,m+1)
       vaux(653,m) = rxpa*vaux(548,m) + rxwp*vaux(548,m+1)
       vaux(654,m) = rxpa*vaux(549,m) + rxwp*vaux(549,m+1)
       vaux(655,m) = rxpa*vaux(550,m) + rxwp*vaux(550,m+1)
       vaux(656,m) = rxpa*vaux(551,m) + rxwp*vaux(551,m+1)
       vaux(657,m) = rxpa*vaux(552,m) + rxwp*vaux(552,m+1)
       vaux(658,m) = rxpa*vaux(553,m) + rxwp*vaux(553,m+1)
       vaux(659,m) = rxpa*vaux(554,m) + rxwp*vaux(554,m+1)
       vaux(660,m) = rxpa*vaux(555,m) + rxwp*vaux(555,m+1)
       vaux(661,m) = rxpa*vaux(556,m) + rxwp*vaux(556,m+1)
       vaux(662,m) = rxpa*vaux(557,m) + rxwp*vaux(557,m+1)
       vaux(663,m) = rxpa*vaux(558,m) + rxwp*vaux(558,m+1)
       vaux(664,m) = rxpa*vaux(559,m) + rxwp*vaux(559,m+1)
       vaux(665,m) = rxpa*vaux(560,m) + rxwp*vaux(560,m+1)
       vaux(666,m) = rypa*vaux(547,m) + rywp*vaux(547,m+1) + f( 13)*(vaux(443,m) - fnbra*vaux(443,m+1))
       vaux(667,m) = rzpa*vaux(547,m) + rzwp*vaux(547,m+1)
       vaux(668,m) = rzpa*vaux(548,m) + rzwp*vaux(548,m+1) + f(  1)*(vaux(443,m) - fnbra*vaux(443,m+1))
       vaux(669,m) = rzpa*vaux(549,m) + rzwp*vaux(549,m+1) + f(  2)*(vaux(444,m) - fnbra*vaux(444,m+1))
       vaux(670,m) = rzpa*vaux(550,m) + rzwp*vaux(550,m+1) + f(  3)*(vaux(445,m) - fnbra*vaux(445,m+1))
       vaux(671,m) = rzpa*vaux(551,m) + rzwp*vaux(551,m+1) + f(  4)*(vaux(446,m) - fnbra*vaux(446,m+1))
       vaux(672,m) = rzpa*vaux(552,m) + rzwp*vaux(552,m+1) + f(  5)*(vaux(447,m) - fnbra*vaux(447,m+1))
       vaux(673,m) = rzpa*vaux(553,m) + rzwp*vaux(553,m+1) + f(  6)*(vaux(448,m) - fnbra*vaux(448,m+1))
       vaux(674,m) = rypa*vaux(555,m) + rywp*vaux(555,m+1) + f(  5)*(vaux(451,m) - fnbra*vaux(451,m+1))
       vaux(675,m) = rypa*vaux(556,m) + rywp*vaux(556,m+1) + f(  4)*(vaux(452,m) - fnbra*vaux(452,m+1))
       vaux(676,m) = rypa*vaux(557,m) + rywp*vaux(557,m+1) + f(  3)*(vaux(453,m) - fnbra*vaux(453,m+1))
       vaux(677,m) = rypa*vaux(558,m) + rywp*vaux(558,m+1) + f(  2)*(vaux(454,m) - fnbra*vaux(454,m+1))
       vaux(678,m) = rypa*vaux(559,m) + rywp*vaux(559,m+1) + f(  1)*(vaux(455,m) - fnbra*vaux(455,m+1))
       vaux(679,m) = rypa*vaux(560,m) + rywp*vaux(560,m+1)
       vaux(680,m) = rzpa*vaux(560,m) + rzwp*vaux(560,m+1) + f( 13)*(vaux(455,m) - fnbra*vaux(455,m+1))
    enddo
    if( lab.eq.15 )then
      if( mmax.eq.1 )then
        call nuclear_getbra(nmin,680,vaux,fint,li1,li2)
      else
        gama( 0) = 0.0d0
        gama( 1) = fnket*fdbra 
        gama( 2) = gama( 1) + gama(1) 
        gama( 3) = gama( 2) + gama(1) 
        gama( 4) = gama( 3) + gama(1) 
        gama( 5) = gama( 4) + gama(1) 
        gama( 6) = gama( 5) + gama(1) 
        gama( 7) = gama( 6) + gama(1) 
        gama( 8) = gama( 7) + gama(1) 
        gama( 9) = gama( 8) + gama(1) 
        gama(10) = gama( 9) + gama(1) 
        gama(11) = gama(10) + gama(1) 
        gama(12) = gama(11) + gama(1) 
        gama(13) = gama(12) + gama(1) 
        gama(14) = gama(13) + gama(1) 
        call nuclear_ket(li1,li2,fint,lf2,mmax,nmin,680,lc,lcd,   &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)
      endif
      return
    endif

    write(6,*) ' Improper choice of the angular momentum l'
    write(6,*) ' expected:  <1-15>'
    write(6,*) ' required: ',lab
    write(6,*) ' Abnormal termination in module: interest_osr, class: nuclear, routine: nuclear_bra'
    stop
  end subroutine

! ---------------------------------------------------------------------------------------
! KET ROUTINE
!
! todo: do buducna doprogramuj zmenu dolnej hranice v cykle po n
! todo: do buducna odstran redundancie v clenoch
! ---------------------------------------------------------------------------------------
  subroutine nuclear_ket(li1,li2,fint,lf2,mmax,nmin,nmax,lc,lcd,  &
                         fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq)

    !-- output --!
    integer, intent(out) :: li1,li2 
    real(8), intent(out) :: fint(*)

    !-- input --!
    integer             :: mmax
    integer, intent(in) :: nmin,nmax
    integer, intent(in) :: lc,lcd,lf2
    real(8), intent(in) :: fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq 

    !-- local --!
    integer :: m,n,mmin
    integer :: mox,moy,moz 
    real(8) :: xmo,ymo,zmo 


    !-- [ 1 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax                                                    
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        waux(1,n,m) =      vaux(n,m)
        waux(2,n,m) = rxqc*vaux(n,m) + rxwq*vaux(n,m+1) + xmo*vaux(mox,m+1) 
        waux(3,n,m) = ryqc*vaux(n,m) + rywq*vaux(n,m+1) + ymo*vaux(moy,m+1)
        waux(4,n,m) = rzqc*vaux(n,m) + rzwq*vaux(n,m+1) + zmo*vaux(moz,m+1)
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 4    - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,4,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 2 | * ] --!
    f(1)=fdket
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(1,n,m) - fnket*waux(1,n,m+1))

        waux( 5,n,m) = rxqc*waux(2,n,m) + rxwq*waux(2,n,m+1) + xmo*waux(2,mox,m+1) + c1 
        waux( 6,n,m) = rxqc*waux(3,n,m) + rxwq*waux(3,n,m+1) + xmo*waux(3,mox,m+1)
        waux( 7,n,m) = rxqc*waux(4,n,m) + rxwq*waux(4,n,m+1) + xmo*waux(4,mox,m+1) 
        waux( 8,n,m) = ryqc*waux(3,n,m) + rywq*waux(3,n,m+1) + ymo*waux(3,moy,m+1) + c1
        waux( 9,n,m) = rzqc*waux(3,n,m) + rzwq*waux(3,n,m+1) + zmo*waux(3,moz,m+1)
        waux(10,n,m) = rzqc*waux(4,n,m) + rzwq*waux(4,n,m+1) + zmo*waux(4,moz,m+1) + c1
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 10   - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,10,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 3 | * ] --!
    f(2)=f(1)+f(1)
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(2)*(waux(2,n,m) - fnket*waux(2,n,m+1))
        c2 = f(2)*(waux(3,n,m) - fnket*waux(3,n,m+1))
        c3 = f(2)*(waux(4,n,m) - fnket*waux(4,n,m+1))
        
        waux(11,n,m) = rxqc*waux( 5,n,m) + rxwq*waux( 5,n,m+1) + xmo*waux( 5,mox,m+1) + c1
        waux(12,n,m) = ryqc*waux( 5,n,m) + rywq*waux( 5,n,m+1) + ymo*waux( 5,moy,m+1)
        waux(13,n,m) = rzqc*waux( 5,n,m) + rzwq*waux( 5,n,m+1) + zmo*waux( 5,moz,m+1)
        waux(14,n,m) = rxqc*waux( 8,n,m) + rxwq*waux( 8,n,m+1) + xmo*waux( 8,mox,m+1)
        waux(15,n,m) = rxqc*waux( 9,n,m) + rxwq*waux( 9,n,m+1) + xmo*waux( 9,mox,m+1)
        waux(16,n,m) = rxqc*waux(10,n,m) + rxwq*waux(10,n,m+1) + xmo*waux(10,mox,m+1)
        waux(17,n,m) = ryqc*waux( 8,n,m) + rywq*waux( 8,n,m+1) + ymo*waux( 8,moy,m+1) + c2
        waux(18,n,m) = rzqc*waux( 8,n,m) + rzwq*waux( 8,n,m+1) + zmo*waux( 8,moz,m+1)
        waux(19,n,m) = ryqc*waux(10,n,m) + rywq*waux(10,n,m+1) + ymo*waux(10,moy,m+1)
        waux(20,n,m) = rzqc*waux(10,n,m) + rzwq*waux(10,n,m+1) + zmo*waux(10,moz,m+1) + c3
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 20   - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,20,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 4 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux( 5,n,m) - fnket*waux( 5,n,m+1))
        c2 = f(1)*(waux( 8,n,m) - fnket*waux( 8,n,m+1))
        c3 = f(1)*(waux(10,n,m) - fnket*waux(10,n,m+1))
        
        waux(21,n,m) = rxqc*waux(11,n,m) + rxwq*waux(11,n,m+1) + xmo*waux(11,mox,m+1) + c1*3.0d0
        waux(22,n,m) = ryqc*waux(11,n,m) + rywq*waux(11,n,m+1) + ymo*waux(11,moy,m+1) 
        waux(23,n,m) = rzqc*waux(11,n,m) + rzwq*waux(11,n,m+1) + zmo*waux(11,moz,m+1) 
        waux(24,n,m) = ryqc*waux(12,n,m) + rywq*waux(12,n,m+1) + ymo*waux(12,moy,m+1) + c1
        waux(25,n,m) = rzqc*waux(12,n,m) + rzwq*waux(12,n,m+1) + zmo*waux(12,moz,m+1) 
        waux(26,n,m) = rzqc*waux(13,n,m) + rzwq*waux(13,n,m+1) + zmo*waux(13,moz,m+1) + c1
        waux(27,n,m) = rxqc*waux(17,n,m) + rxwq*waux(17,n,m+1) + xmo*waux(17,mox,m+1) 
        waux(28,n,m) = rxqc*waux(18,n,m) + rxwq*waux(18,n,m+1) + xmo*waux(18,mox,m+1) 
        waux(29,n,m) = rxqc*waux(19,n,m) + rxwq*waux(19,n,m+1) + xmo*waux(19,mox,m+1) 
        waux(30,n,m) = rxqc*waux(20,n,m) + rxwq*waux(20,n,m+1) + xmo*waux(20,mox,m+1) 
        waux(31,n,m) = ryqc*waux(17,n,m) + rywq*waux(17,n,m+1) + ymo*waux(17,moy,m+1) + c2*3.0d0
        waux(32,n,m) = rzqc*waux(17,n,m) + rzwq*waux(17,n,m+1) + zmo*waux(17,moz,m+1) 
        waux(33,n,m) = rzqc*waux(18,n,m) + rzwq*waux(18,n,m+1) + zmo*waux(18,moz,m+1) + c2
        waux(34,n,m) = ryqc*waux(20,n,m) + rywq*waux(20,n,m+1) + ymo*waux(20,moy,m+1) 
        waux(35,n,m) = rzqc*waux(20,n,m) + rzwq*waux(20,n,m+1) + zmo*waux(20,moz,m+1) + c3*3.0d0
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 35   - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,35,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 5 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(11,n,m) - fnket*waux(11,n,m+1))
        c2 = f(1)*(waux(17,n,m) - fnket*waux(17,n,m+1))
        c3 = f(1)*(waux(20,n,m) - fnket*waux(20,n,m+1))
        
        waux(36,n,m) = rxqc*waux(21,n,m) + rxwq*waux(21,n,m+1) + xmo*waux(21,mox,m+1) + c1*4.0d0
        waux(37,n,m) = ryqc*waux(21,n,m) + rywq*waux(21,n,m+1) + ymo*waux(21,moy,m+1) 
        waux(38,n,m) = rzqc*waux(21,n,m) + rzwq*waux(21,n,m+1) + zmo*waux(21,moz,m+1) 
        waux(39,n,m) = ryqc*waux(22,n,m) + rywq*waux(22,n,m+1) + ymo*waux(22,moy,m+1) + c1         
        waux(40,n,m) = rzqc*waux(22,n,m) + rzwq*waux(22,n,m+1) + zmo*waux(22,moz,m+1) 
        waux(41,n,m) = rzqc*waux(23,n,m) + rzwq*waux(23,n,m+1) + zmo*waux(23,moz,m+1) + c1        
        waux(42,n,m) = rxqc*waux(27,n,m) + rxwq*waux(27,n,m+1) + xmo*waux(27,mox,m+1) + c2       
        waux(43,n,m) = rzqc*waux(24,n,m) + rzwq*waux(24,n,m+1) + zmo*waux(24,moz,m+1) 
        waux(44,n,m) = ryqc*waux(26,n,m) + rywq*waux(26,n,m+1) + ymo*waux(26,moy,m+1) 
        waux(45,n,m) = rxqc*waux(30,n,m) + rxwq*waux(30,n,m+1) + xmo*waux(30,mox,m+1) + c3        
        waux(46,n,m) = rxqc*waux(31,n,m) + rxwq*waux(31,n,m+1) + xmo*waux(31,mox,m+1) 
        waux(47,n,m) = rxqc*waux(32,n,m) + rxwq*waux(32,n,m+1) + xmo*waux(32,mox,m+1) 
        waux(48,n,m) = rxqc*waux(33,n,m) + rxwq*waux(33,n,m+1) + xmo*waux(33,mox,m+1) 
        waux(49,n,m) = rxqc*waux(34,n,m) + rxwq*waux(34,n,m+1) + xmo*waux(34,mox,m+1) 
        waux(50,n,m) = rxqc*waux(35,n,m) + rxwq*waux(35,n,m+1) + xmo*waux(35,mox,m+1) 
        waux(51,n,m) = ryqc*waux(31,n,m) + rywq*waux(31,n,m+1) + ymo*waux(31,moy,m+1) + c2*4.0d0
        waux(52,n,m) = rzqc*waux(31,n,m) + rzwq*waux(31,n,m+1) + zmo*waux(31,moz,m+1) 
        waux(53,n,m) = rzqc*waux(32,n,m) + rzwq*waux(32,n,m+1) + zmo*waux(32,moz,m+1) + c2         
        waux(54,n,m) = ryqc*waux(34,n,m) + rywq*waux(34,n,m+1) + ymo*waux(34,moy,m+1) + c3          
        waux(55,n,m) = ryqc*waux(35,n,m) + rywq*waux(35,n,m+1) + ymo*waux(35,moy,m+1) 
        waux(56,n,m) = rzqc*waux(35,n,m) + rzwq*waux(35,n,m+1) + zmo*waux(35,moz,m+1) + c3*4.0d0
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 56   - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,56,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 6 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(21,n,m) - fnket*waux(21,n,m+1))
        c4 = f(2)*(waux(22,n,m) - fnket*waux(22,n,m+1))
        c5 = f(2)*(waux(23,n,m) - fnket*waux(23,n,m+1))
        c6 = f(1)*(waux(24,n,m) - fnket*waux(24,n,m+1))
        c2 = f(1)*(waux(31,n,m) - fnket*waux(31,n,m+1))
        c7 = f(2)*(waux(32,n,m) - fnket*waux(32,n,m+1))
        c3 = f(1)*(waux(35,n,m) - fnket*waux(35,n,m+1))
        
        waux(57,n,m) = rxqc*waux(36,n,m) + rxwq*waux(36,n,m+1) + xmo*waux(36,mox,m+1) + c1*5.0d0
        waux(58,n,m) = ryqc*waux(36,n,m) + rywq*waux(36,n,m+1) + ymo*waux(36,moy,m+1) 
        waux(59,n,m) = rzqc*waux(36,n,m) + rzwq*waux(36,n,m+1) + zmo*waux(36,moz,m+1) 
        waux(60,n,m) = ryqc*waux(37,n,m) + rywq*waux(37,n,m+1) + ymo*waux(37,moy,m+1) + c1 
        waux(61,n,m) = rzqc*waux(37,n,m) + rzwq*waux(37,n,m+1) + zmo*waux(37,moz,m+1) 
        waux(62,n,m) = rzqc*waux(38,n,m) + rzwq*waux(38,n,m+1) + zmo*waux(38,moz,m+1) + c1 
        waux(63,n,m) = ryqc*waux(39,n,m) + rywq*waux(39,n,m+1) + ymo*waux(39,moy,m+1) + c4 
        waux(64,n,m) = rzqc*waux(39,n,m) + rzwq*waux(39,n,m+1) + zmo*waux(39,moz,m+1) 
        waux(65,n,m) = ryqc*waux(41,n,m) + rywq*waux(41,n,m+1) + ymo*waux(41,moy,m+1) 
        waux(66,n,m) = rzqc*waux(41,n,m) + rzwq*waux(41,n,m+1) + zmo*waux(41,moz,m+1) + c5 
        waux(67,n,m) = rxqc*waux(46,n,m) + rxwq*waux(46,n,m+1) + xmo*waux(46,mox,m+1) + c2 
        waux(68,n,m) = rzqc*waux(42,n,m) + rzwq*waux(42,n,m+1) + zmo*waux(42,moz,m+1) 
        waux(69,n,m) = rzqc*waux(43,n,m) + rzwq*waux(43,n,m+1) + zmo*waux(43,moz,m+1) + c6 
        waux(70,n,m) = ryqc*waux(45,n,m) + rywq*waux(45,n,m+1) + ymo*waux(45,moy,m+1) 
        waux(71,n,m) = rxqc*waux(50,n,m) + rxwq*waux(50,n,m+1) + xmo*waux(50,mox,m+1) + c3 
        waux(72,n,m) = rxqc*waux(51,n,m) + rxwq*waux(51,n,m+1) + xmo*waux(51,mox,m+1) 
        waux(73,n,m) = rxqc*waux(52,n,m) + rxwq*waux(52,n,m+1) + xmo*waux(52,mox,m+1) 
        waux(74,n,m) = rxqc*waux(53,n,m) + rxwq*waux(53,n,m+1) + xmo*waux(53,mox,m+1) 
        waux(75,n,m) = rxqc*waux(54,n,m) + rxwq*waux(54,n,m+1) + xmo*waux(54,mox,m+1) 
        waux(76,n,m) = rxqc*waux(55,n,m) + rxwq*waux(55,n,m+1) + xmo*waux(55,mox,m+1) 
        waux(77,n,m) = rxqc*waux(56,n,m) + rxwq*waux(56,n,m+1) + xmo*waux(56,mox,m+1) 
        waux(78,n,m) = ryqc*waux(51,n,m) + rywq*waux(51,n,m+1) + ymo*waux(51,moy,m+1) + c2*5.0d0
        waux(79,n,m) = rzqc*waux(51,n,m) + rzwq*waux(51,n,m+1) + zmo*waux(51,moz,m+1) 
        waux(80,n,m) = rzqc*waux(52,n,m) + rzwq*waux(52,n,m+1) + zmo*waux(52,moz,m+1) + c2 
        waux(81,n,m) = rzqc*waux(53,n,m) + rzwq*waux(53,n,m+1) + zmo*waux(53,moz,m+1) + c7 
        waux(82,n,m) = ryqc*waux(55,n,m) + rywq*waux(55,n,m+1) + ymo*waux(55,moy,m+1) + c3 
        waux(83,n,m) = ryqc*waux(56,n,m) + rywq*waux(56,n,m+1) + ymo*waux(56,moy,m+1) 
        waux(84,n,m) = rzqc*waux(56,n,m) + rzwq*waux(56,n,m+1) + zmo*waux(56,moz,m+1) + c3*5.0d0
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 84   - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,84,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 7 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(36,n,m) - fnket*waux(36,n,m+1))
        c4 = f(2)*(waux(37,n,m) - fnket*waux(37,n,m+1))
        c5 = f(2)*(waux(38,n,m) - fnket*waux(38,n,m+1))
        c6 = f(1)*(waux(39,n,m) - fnket*waux(39,n,m+1))
        c7 = f(1)*(waux(42,n,m) - fnket*waux(42,n,m+1))
        c8 = f(1)*(waux(45,n,m) - fnket*waux(45,n,m+1))
        c9 = f(2)*(waux(46,n,m) - fnket*waux(46,n,m+1))
        c10= f(2)*(waux(50,n,m) - fnket*waux(50,n,m+1))
        c2 = f(1)*(waux(51,n,m) - fnket*waux(51,n,m+1))
        c11= f(2)*(waux(52,n,m) - fnket*waux(52,n,m+1))
        c12= f(2)*(waux(55,n,m) - fnket*waux(55,n,m+1))
        c3 = f(1)*(waux(56,n,m) - fnket*waux(56,n,m+1))
        
        waux( 85,n,m) = rxqc*waux(57,n,m) + rxwq*waux(57,n,m+1) + xmo*waux(57,mox,m+1) + c1*6.0d0 
        waux( 86,n,m) = ryqc*waux(57,n,m) + rywq*waux(57,n,m+1) + ymo*waux(57,moy,m+1) 
        waux( 87,n,m) = rzqc*waux(57,n,m) + rzwq*waux(57,n,m+1) + zmo*waux(57,moz,m+1) 
        waux( 88,n,m) = ryqc*waux(58,n,m) + rywq*waux(58,n,m+1) + ymo*waux(58,moy,m+1) + c1 
        waux( 89,n,m) = rzqc*waux(58,n,m) + rzwq*waux(58,n,m+1) + zmo*waux(58,moz,m+1) 
        waux( 90,n,m) = rzqc*waux(59,n,m) + rzwq*waux(59,n,m+1) + zmo*waux(59,moz,m+1) + c1 
        waux( 91,n,m) = ryqc*waux(60,n,m) + rywq*waux(60,n,m+1) + ymo*waux(60,moy,m+1) + c4 
        waux( 92,n,m) = rzqc*waux(60,n,m) + rzwq*waux(60,n,m+1) + zmo*waux(60,moz,m+1) 
        waux( 93,n,m) = ryqc*waux(62,n,m) + rywq*waux(62,n,m+1) + ymo*waux(62,moy,m+1) 
        waux( 94,n,m) = rzqc*waux(62,n,m) + rzwq*waux(62,n,m+1) + zmo*waux(62,moz,m+1) + c5 
        waux( 95,n,m) = rxqc*waux(67,n,m) + rxwq*waux(67,n,m+1) + xmo*waux(67,mox,m+1) + c9 
        waux( 96,n,m) = rzqc*waux(63,n,m) + rzwq*waux(63,n,m+1) + zmo*waux(63,moz,m+1) 
        waux( 97,n,m) = rzqc*waux(64,n,m) + rzwq*waux(64,n,m+1) + zmo*waux(64,moz,m+1) + c6 
        waux( 98,n,m) = ryqc*waux(66,n,m) + rywq*waux(66,n,m+1) + ymo*waux(66,moy,m+1) 
        waux( 99,n,m) = rxqc*waux(71,n,m) + rxwq*waux(71,n,m+1) + xmo*waux(71,mox,m+1) + c10
        waux(100,n,m) = rxqc*waux(72,n,m) + rxwq*waux(72,n,m+1) + xmo*waux(72,mox,m+1) + c2 
        waux(101,n,m) = rzqc*waux(67,n,m) + rzwq*waux(67,n,m+1) + zmo*waux(67,moz,m+1) 
        waux(102,n,m) = rzqc*waux(68,n,m) + rzwq*waux(68,n,m+1) + zmo*waux(68,moz,m+1) + c7 
        waux(103,n,m) = ryqc*waux(70,n,m) + rywq*waux(70,n,m+1) + ymo*waux(70,moy,m+1) + c8 
        waux(104,n,m) = ryqc*waux(71,n,m) + rywq*waux(71,n,m+1) + ymo*waux(71,moy,m+1) 
        waux(105,n,m) = rxqc*waux(77,n,m) + rxwq*waux(77,n,m+1) + xmo*waux(77,mox,m+1) + c3 
        waux(106,n,m) = rxqc*waux(78,n,m) + rxwq*waux(78,n,m+1) + xmo*waux(78,mox,m+1) 
        waux(107,n,m) = rxqc*waux(79,n,m) + rxwq*waux(79,n,m+1) + xmo*waux(79,mox,m+1) 
        waux(108,n,m) = rxqc*waux(80,n,m) + rxwq*waux(80,n,m+1) + xmo*waux(80,mox,m+1) 
        waux(109,n,m) = rxqc*waux(81,n,m) + rxwq*waux(81,n,m+1) + xmo*waux(81,mox,m+1) 
        waux(110,n,m) = rxqc*waux(82,n,m) + rxwq*waux(82,n,m+1) + xmo*waux(82,mox,m+1) 
        waux(111,n,m) = rxqc*waux(83,n,m) + rxwq*waux(83,n,m+1) + xmo*waux(83,mox,m+1) 
        waux(112,n,m) = rxqc*waux(84,n,m) + rxwq*waux(84,n,m+1) + xmo*waux(84,mox,m+1) 
        waux(113,n,m) = ryqc*waux(78,n,m) + rywq*waux(78,n,m+1) + ymo*waux(78,moy,m+1) + c2*6.0d0 
        waux(114,n,m) = rzqc*waux(78,n,m) + rzwq*waux(78,n,m+1) + zmo*waux(78,moz,m+1) 
        waux(115,n,m) = rzqc*waux(79,n,m) + rzwq*waux(79,n,m+1) + zmo*waux(79,moz,m+1) + c2 
        waux(116,n,m) = rzqc*waux(80,n,m) + rzwq*waux(80,n,m+1) + zmo*waux(80,moz,m+1) + c11 
        waux(117,n,m) = ryqc*waux(82,n,m) + rywq*waux(82,n,m+1) + ymo*waux(82,moy,m+1) + c12 
        waux(118,n,m) = ryqc*waux(83,n,m) + rywq*waux(83,n,m+1) + ymo*waux(83,moy,m+1) + c3 
        waux(119,n,m) = ryqc*waux(84,n,m) + rywq*waux(84,n,m+1) + ymo*waux(84,moy,m+1) 
        waux(120,n,m) = rzqc*waux(84,n,m) + rzwq*waux(84,n,m+1) + zmo*waux(84,moz,m+1) + c3*6.0d0 
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 120  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,120,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 8 | * ]--!
    f(3)=f(2)+f(1)
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux( 57,n,m) - fnket*waux( 57,n,m+1))
        c4 = f(2)*(waux( 58,n,m) - fnket*waux( 58,n,m+1))
        c5 = f(2)*(waux( 59,n,m) - fnket*waux( 59,n,m+1))
        c6 = f(1)*(waux( 60,n,m) - fnket*waux( 60,n,m+1))
        c7 = f(3)*(waux( 62,n,m) - fnket*waux( 62,n,m+1))
        c8 = f(1)*(waux( 63,n,m) - fnket*waux( 63,n,m+1))
        c9 = f(1)*(waux( 66,n,m) - fnket*waux( 66,n,m+1))
        c10= f(1)*(waux( 67,n,m) - fnket*waux( 67,n,m+1))
        c11= f(1)*(waux( 71,n,m) - fnket*waux( 71,n,m+1))
        c12= f(2)*(waux( 72,n,m) - fnket*waux( 72,n,m+1))
        c13= f(2)*(waux( 77,n,m) - fnket*waux( 77,n,m+1))
        c2 = f(1)*(waux( 78,n,m) - fnket*waux( 78,n,m+1))
        c14= f(2)*(waux( 79,n,m) - fnket*waux( 79,n,m+1))
        c15= f(3)*(waux( 80,n,m) - fnket*waux( 80,n,m+1))
        c16= f(1)*(waux( 81,n,m) - fnket*waux( 81,n,m+1))
        c17= f(2)*(waux( 83,n,m) - fnket*waux( 83,n,m+1))
        c3 = f(1)*(waux( 84,n,m) - fnket*waux( 84,n,m+1))
        
        waux(121,n,m) = rxqc*waux( 85,n,m) + rxwq*waux( 85,n,m+1) + xmo*waux( 85,mox,m+1) + c1*7.0d0 
        waux(122,n,m) = ryqc*waux( 85,n,m) + rywq*waux( 85,n,m+1) + ymo*waux( 85,moy,m+1) 
        waux(123,n,m) = rzqc*waux( 85,n,m) + rzwq*waux( 85,n,m+1) + zmo*waux( 85,moz,m+1) 
        waux(124,n,m) = ryqc*waux( 86,n,m) + rywq*waux( 86,n,m+1) + ymo*waux( 86,moy,m+1) + c1 
        waux(125,n,m) = rzqc*waux( 86,n,m) + rzwq*waux( 86,n,m+1) + zmo*waux( 86,moz,m+1) 
        waux(126,n,m) = rzqc*waux( 87,n,m) + rzwq*waux( 87,n,m+1) + zmo*waux( 87,moz,m+1) + c1 
        waux(127,n,m) = ryqc*waux( 88,n,m) + rywq*waux( 88,n,m+1) + ymo*waux( 88,moy,m+1) + c4 
        waux(128,n,m) = rzqc*waux( 88,n,m) + rzwq*waux( 88,n,m+1) + zmo*waux( 88,moz,m+1) 
        waux(129,n,m) = ryqc*waux( 90,n,m) + rywq*waux( 90,n,m+1) + ymo*waux( 90,moy,m+1) 
        waux(130,n,m) = rzqc*waux( 90,n,m) + rzwq*waux( 90,n,m+1) + zmo*waux( 90,moz,m+1) + c5 
        waux(131,n,m) = ryqc*waux( 91,n,m) + rywq*waux( 91,n,m+1) + ymo*waux( 91,moy,m+1) + c6*3.0d0
        waux(132,n,m) = rzqc*waux( 91,n,m) + rzwq*waux( 91,n,m+1) + zmo*waux( 91,moz,m+1) 
        waux(133,n,m) = rzqc*waux( 92,n,m) + rzwq*waux( 92,n,m+1) + zmo*waux( 92,moz,m+1) + c6 
        waux(134,n,m) = ryqc*waux( 94,n,m) + rywq*waux( 94,n,m+1) + ymo*waux( 94,moy,m+1) 
        waux(135,n,m) = rzqc*waux( 94,n,m) + rzwq*waux( 94,n,m+1) + zmo*waux( 94,moz,m+1) + c7 
        waux(136,n,m) = rxqc*waux(100,n,m) + rxwq*waux(100,n,m+1) + xmo*waux(100,mox,m+1) + c12 
        waux(137,n,m) = rzqc*waux( 95,n,m) + rzwq*waux( 95,n,m+1) + zmo*waux( 95,moz,m+1) 
        waux(138,n,m) = rzqc*waux( 96,n,m) + rzwq*waux( 96,n,m+1) + zmo*waux( 96,moz,m+1) + c8 
        waux(139,n,m) = ryqc*waux( 98,n,m) + rywq*waux( 98,n,m+1) + ymo*waux( 98,moy,m+1) + c9 
        waux(140,n,m) = ryqc*waux( 99,n,m) + rywq*waux( 99,n,m+1) + ymo*waux( 99,moy,m+1) 
        waux(141,n,m) = rxqc*waux(105,n,m) + rxwq*waux(105,n,m+1) + xmo*waux(105,mox,m+1) + c13 
        waux(142,n,m) = rxqc*waux(106,n,m) + rxwq*waux(106,n,m+1) + xmo*waux(106,mox,m+1) + c2 
        waux(143,n,m) = rzqc*waux(100,n,m) + rzwq*waux(100,n,m+1) + zmo*waux(100,moz,m+1) 
        waux(144,n,m) = rzqc*waux(101,n,m) + rzwq*waux(101,n,m+1) + zmo*waux(101,moz,m+1) + c10 
        waux(145,n,m) = rxqc*waux(109,n,m) + rxwq*waux(109,n,m+1) + xmo*waux(109,mox,m+1) + c16 
        waux(146,n,m) = ryqc*waux(104,n,m) + rywq*waux(104,n,m+1) + ymo*waux(104,moy,m+1) + c11 
        waux(147,n,m) = ryqc*waux(105,n,m) + rywq*waux(105,n,m+1) + ymo*waux(105,moy,m+1) 
        waux(148,n,m) = rxqc*waux(112,n,m) + rxwq*waux(112,n,m+1) + xmo*waux(112,mox,m+1) + c3 
        waux(149,n,m) = rxqc*waux(113,n,m) + rxwq*waux(113,n,m+1) + xmo*waux(113,mox,m+1) 
        waux(150,n,m) = rxqc*waux(114,n,m) + rxwq*waux(114,n,m+1) + xmo*waux(114,mox,m+1) 
        waux(151,n,m) = rxqc*waux(115,n,m) + rxwq*waux(115,n,m+1) + xmo*waux(115,mox,m+1) 
        waux(152,n,m) = rxqc*waux(116,n,m) + rxwq*waux(116,n,m+1) + xmo*waux(116,mox,m+1) 
        waux(153,n,m) = rxqc*waux(117,n,m) + rxwq*waux(117,n,m+1) + xmo*waux(117,mox,m+1) 
        waux(154,n,m) = rxqc*waux(118,n,m) + rxwq*waux(118,n,m+1) + xmo*waux(118,mox,m+1) 
        waux(155,n,m) = rxqc*waux(119,n,m) + rxwq*waux(119,n,m+1) + xmo*waux(119,mox,m+1) 
        waux(156,n,m) = rxqc*waux(120,n,m) + rxwq*waux(120,n,m+1) + xmo*waux(120,mox,m+1) 
        waux(157,n,m) = ryqc*waux(113,n,m) + rywq*waux(113,n,m+1) + ymo*waux(113,moy,m+1) + c2*7.0d0 
        waux(158,n,m) = rzqc*waux(113,n,m) + rzwq*waux(113,n,m+1) + zmo*waux(113,moz,m+1) 
        waux(159,n,m) = rzqc*waux(114,n,m) + rzwq*waux(114,n,m+1) + zmo*waux(114,moz,m+1) + c2 
        waux(160,n,m) = rzqc*waux(115,n,m) + rzwq*waux(115,n,m+1) + zmo*waux(115,moz,m+1) + c14 
        waux(161,n,m) = rzqc*waux(116,n,m) + rzwq*waux(116,n,m+1) + zmo*waux(116,moz,m+1) + c15 
        waux(162,n,m) = ryqc*waux(118,n,m) + rywq*waux(118,n,m+1) + ymo*waux(118,moy,m+1) + c17 
        waux(163,n,m) = ryqc*waux(119,n,m) + rywq*waux(119,n,m+1) + ymo*waux(119,moy,m+1) + c3 
        waux(164,n,m) = ryqc*waux(120,n,m) + rywq*waux(120,n,m+1) + ymo*waux(120,moy,m+1) 
        waux(165,n,m) = rzqc*waux(120,n,m) + rzwq*waux(120,n,m+1) + zmo*waux(120,moz,m+1) + c3*7.0d0 
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 165  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,165,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 9 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux( 85,n,m) - fnket*waux( 85,n,m+1))
        c2 = f(2)*(waux( 86,n,m) - fnket*waux( 86,n,m+1))
        c3 = f(2)*(waux( 87,n,m) - fnket*waux( 87,n,m+1))
        c4 = f(1)*(waux( 88,n,m) - fnket*waux( 88,n,m+1))
        c5 = f(3)*(waux( 90,n,m) - fnket*waux( 90,n,m+1))
        c6 = f(1)*(waux( 91,n,m) - fnket*waux( 91,n,m+1))
        c7 = f(1)*(waux( 94,n,m) - fnket*waux( 94,n,m+1))
        c8 = f(1)*(waux( 95,n,m) - fnket*waux( 95,n,m+1))
        c9 = f(2)*(waux( 96,n,m) - fnket*waux( 96,n,m+1))
        c10= f(1)*(waux( 99,n,m) - fnket*waux( 99,n,m+1))
        c11= f(1)*(waux(100,n,m) - fnket*waux(100,n,m+1))
        c12= f(1)*(waux(105,n,m) - fnket*waux(105,n,m+1))
        c13= f(2)*(waux(106,n,m) - fnket*waux(106,n,m+1))
        c14= f(2)*(waux(112,n,m) - fnket*waux(112,n,m+1))
        c15= f(1)*(waux(113,n,m) - fnket*waux(113,n,m+1))
        c16= f(2)*(waux(114,n,m) - fnket*waux(114,n,m+1))
        c17= f(3)*(waux(115,n,m) - fnket*waux(115,n,m+1))
        c18= f(1)*(waux(116,n,m) - fnket*waux(116,n,m+1))
        c19= f(1)*(waux(117,n,m) - fnket*waux(117,n,m+1))
        c20= f(3)*(waux(118,n,m) - fnket*waux(118,n,m+1))
        c21= f(2)*(waux(119,n,m) - fnket*waux(119,n,m+1))
        c22= f(1)*(waux(120,n,m) - fnket*waux(120,n,m+1))
        
        waux(166,n,m) = rxqc*waux(121,n,m) + rxwq*waux(121,n,m+1) + xmo*waux(121,mox,m+1) + c1*8.0d0 
        waux(167,n,m) = ryqc*waux(121,n,m) + rywq*waux(121,n,m+1) + ymo*waux(121,moy,m+1) 
        waux(168,n,m) = rzqc*waux(121,n,m) + rzwq*waux(121,n,m+1) + zmo*waux(121,moz,m+1) 
        waux(169,n,m) = ryqc*waux(122,n,m) + rywq*waux(122,n,m+1) + ymo*waux(122,moy,m+1) + c1 
        waux(170,n,m) = rzqc*waux(122,n,m) + rzwq*waux(122,n,m+1) + zmo*waux(122,moz,m+1) 
        waux(171,n,m) = rzqc*waux(123,n,m) + rzwq*waux(123,n,m+1) + zmo*waux(123,moz,m+1) + c1 
        waux(172,n,m) = ryqc*waux(124,n,m) + rywq*waux(124,n,m+1) + ymo*waux(124,moy,m+1) + c2 
        waux(173,n,m) = rzqc*waux(124,n,m) + rzwq*waux(124,n,m+1) + zmo*waux(124,moz,m+1) 
        waux(174,n,m) = ryqc*waux(126,n,m) + rywq*waux(126,n,m+1) + ymo*waux(126,moy,m+1) 
        waux(175,n,m) = rzqc*waux(126,n,m) + rzwq*waux(126,n,m+1) + zmo*waux(126,moz,m+1) + c3 
        waux(176,n,m) = ryqc*waux(127,n,m) + rywq*waux(127,n,m+1) + ymo*waux(127,moy,m+1) + c4*3.0d0
        waux(177,n,m) = rzqc*waux(127,n,m) + rzwq*waux(127,n,m+1) + zmo*waux(127,moz,m+1) 
        waux(178,n,m) = rzqc*waux(128,n,m) + rzwq*waux(128,n,m+1) + zmo*waux(128,moz,m+1) + c4 
        waux(179,n,m) = ryqc*waux(130,n,m) + rywq*waux(130,n,m+1) + ymo*waux(130,moy,m+1) 
        waux(180,n,m) = rzqc*waux(130,n,m) + rzwq*waux(130,n,m+1) + zmo*waux(130,moz,m+1) + c5 
        waux(181,n,m) = rxqc*waux(136,n,m) + rxwq*waux(136,n,m+1) + xmo*waux(136,mox,m+1) + c11*3.0d0 
        waux(182,n,m) = rzqc*waux(131,n,m) + rzwq*waux(131,n,m+1) + zmo*waux(131,moz,m+1) 
        waux(183,n,m) = rzqc*waux(132,n,m) + rzwq*waux(132,n,m+1) + zmo*waux(132,moz,m+1) + c6 
        waux(184,n,m) = ryqc*waux(134,n,m) + rywq*waux(134,n,m+1) + ymo*waux(134,moy,m+1) + c7 
        waux(185,n,m) = ryqc*waux(135,n,m) + rywq*waux(135,n,m+1) + ymo*waux(135,moy,m+1) 
        waux(186,n,m) = rxqc*waux(141,n,m) + rxwq*waux(141,n,m+1) + xmo*waux(141,mox,m+1) + c12*3.0d0 
        waux(187,n,m) = rxqc*waux(142,n,m) + rxwq*waux(142,n,m+1) + xmo*waux(142,mox,m+1) + c13 
        waux(188,n,m) = rzqc*waux(136,n,m) + rzwq*waux(136,n,m+1) + zmo*waux(136,moz,m+1) 
        waux(189,n,m) = rzqc*waux(137,n,m) + rzwq*waux(137,n,m+1) + zmo*waux(137,moz,m+1) + c8 
        waux(190,n,m) = rzqc*waux(138,n,m) + rzwq*waux(138,n,m+1) + zmo*waux(138,moz,m+1) + c9 
        waux(191,n,m) = ryqc*waux(140,n,m) + rywq*waux(140,n,m+1) + ymo*waux(140,moy,m+1) + c10 
        waux(192,n,m) = ryqc*waux(141,n,m) + rywq*waux(141,n,m+1) + ymo*waux(141,moy,m+1) 
        waux(193,n,m) = rxqc*waux(148,n,m) + rxwq*waux(148,n,m+1) + xmo*waux(148,mox,m+1) + c14 
        waux(194,n,m) = rxqc*waux(149,n,m) + rxwq*waux(149,n,m+1) + xmo*waux(149,mox,m+1) + c15 
        waux(195,n,m) = rzqc*waux(142,n,m) + rzwq*waux(142,n,m+1) + zmo*waux(142,moz,m+1) 
        waux(196,n,m) = rzqc*waux(143,n,m) + rzwq*waux(143,n,m+1) + zmo*waux(143,moz,m+1) + c11 
        waux(197,n,m) = rxqc*waux(152,n,m) + rxwq*waux(152,n,m+1) + xmo*waux(152,mox,m+1) + c18 
        waux(198,n,m) = rxqc*waux(153,n,m) + rxwq*waux(153,n,m+1) + xmo*waux(153,mox,m+1) + c19 
        waux(199,n,m) = ryqc*waux(147,n,m) + rywq*waux(147,n,m+1) + ymo*waux(147,moy,m+1) + c12 
        waux(200,n,m) = ryqc*waux(148,n,m) + rywq*waux(148,n,m+1) + ymo*waux(148,moy,m+1) 
        waux(201,n,m) = rxqc*waux(156,n,m) + rxwq*waux(156,n,m+1) + xmo*waux(156,mox,m+1) + c22 
        waux(202,n,m) = rxqc*waux(157,n,m) + rxwq*waux(157,n,m+1) + xmo*waux(157,mox,m+1) 
        waux(203,n,m) = rxqc*waux(158,n,m) + rxwq*waux(158,n,m+1) + xmo*waux(158,mox,m+1) 
        waux(204,n,m) = rxqc*waux(159,n,m) + rxwq*waux(159,n,m+1) + xmo*waux(159,mox,m+1) 
        waux(205,n,m) = rxqc*waux(160,n,m) + rxwq*waux(160,n,m+1) + xmo*waux(160,mox,m+1) 
        waux(206,n,m) = rxqc*waux(161,n,m) + rxwq*waux(161,n,m+1) + xmo*waux(161,mox,m+1) 
        waux(207,n,m) = rxqc*waux(162,n,m) + rxwq*waux(162,n,m+1) + xmo*waux(162,mox,m+1) 
        waux(208,n,m) = rxqc*waux(163,n,m) + rxwq*waux(163,n,m+1) + xmo*waux(163,mox,m+1) 
        waux(209,n,m) = rxqc*waux(164,n,m) + rxwq*waux(164,n,m+1) + xmo*waux(164,mox,m+1) 
        waux(210,n,m) = rxqc*waux(165,n,m) + rxwq*waux(165,n,m+1) + xmo*waux(165,mox,m+1) 
        waux(211,n,m) = ryqc*waux(157,n,m) + rywq*waux(157,n,m+1) + ymo*waux(157,moy,m+1) + c15*8.0d0 
        waux(212,n,m) = rzqc*waux(157,n,m) + rzwq*waux(157,n,m+1) + zmo*waux(157,moz,m+1) 
        waux(213,n,m) = rzqc*waux(158,n,m) + rzwq*waux(158,n,m+1) + zmo*waux(158,moz,m+1) + c15 
        waux(214,n,m) = rzqc*waux(159,n,m) + rzwq*waux(159,n,m+1) + zmo*waux(159,moz,m+1) + c16 
        waux(215,n,m) = rzqc*waux(160,n,m) + rzwq*waux(160,n,m+1) + zmo*waux(160,moz,m+1) + c17 
        waux(216,n,m) = ryqc*waux(162,n,m) + rywq*waux(162,n,m+1) + ymo*waux(162,moy,m+1) + c20 
        waux(217,n,m) = ryqc*waux(163,n,m) + rywq*waux(163,n,m+1) + ymo*waux(163,moy,m+1) + c21 
        waux(218,n,m) = ryqc*waux(164,n,m) + rywq*waux(164,n,m+1) + ymo*waux(164,moy,m+1) + c22 
        waux(219,n,m) = ryqc*waux(165,n,m) + rywq*waux(165,n,m+1) + ymo*waux(165,moy,m+1) 
        waux(220,n,m) = rzqc*waux(165,n,m) + rzwq*waux(165,n,m+1) + zmo*waux(165,moz,m+1) + c22*8.0d0 
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 220  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,220,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [ 10 | * ]--!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(121,n,m) - fnket*waux(121,n,m+1))
        c2 = f(2)*(waux(122,n,m) - fnket*waux(122,n,m+1))
        c3 = f(2)*(waux(123,n,m) - fnket*waux(123,n,m+1))
        c4 = f(1)*(waux(124,n,m) - fnket*waux(124,n,m+1))
        c5 = f(3)*(waux(126,n,m) - fnket*waux(126,n,m+1))
        c6 = f(1)*(waux(127,n,m) - fnket*waux(127,n,m+1))
        c7 = f(1)*(waux(130,n,m) - fnket*waux(130,n,m+1))
        c8 = f(1)*(waux(131,n,m) - fnket*waux(131,n,m+1))
        c9 = f(2)*(waux(132,n,m) - fnket*waux(132,n,m+1))
        c10= f(1)*(waux(135,n,m) - fnket*waux(135,n,m+1))
        c11= f(1)*(waux(136,n,m) - fnket*waux(136,n,m+1))
        c12= f(2)*(waux(137,n,m) - fnket*waux(137,n,m+1))
        c13= f(2)*(waux(140,n,m) - fnket*waux(140,n,m+1))
        c14= f(1)*(waux(141,n,m) - fnket*waux(141,n,m+1))
        c15= f(1)*(waux(142,n,m) - fnket*waux(142,n,m+1))
        c16= f(1)*(waux(148,n,m) - fnket*waux(148,n,m+1))
        c17= f(2)*(waux(149,n,m) - fnket*waux(149,n,m+1))
        c18= f(2)*(waux(156,n,m) - fnket*waux(156,n,m+1))
        c19= f(1)*(waux(157,n,m) - fnket*waux(157,n,m+1))
        c20= f(2)*(waux(158,n,m) - fnket*waux(158,n,m+1))
        c21= f(3)*(waux(159,n,m) - fnket*waux(159,n,m+1))
        c22= f(1)*(waux(160,n,m) - fnket*waux(160,n,m+1))
        c23= f(1)*(waux(161,n,m) - fnket*waux(161,n,m+1))
        c24= f(1)*(waux(162,n,m) - fnket*waux(162,n,m+1))
        c25= f(3)*(waux(163,n,m) - fnket*waux(163,n,m+1))
        c26= f(2)*(waux(164,n,m) - fnket*waux(164,n,m+1))
        c27= f(1)*(waux(165,n,m) - fnket*waux(165,n,m+1))
        
        waux(221,n,m) = rxqc*waux(166,n,m) + rxwq*waux(166,n,m+1) + xmo*waux(166,mox,m+1) + c1*9.0d0 
        waux(222,n,m) = ryqc*waux(166,n,m) + rywq*waux(166,n,m+1) + ymo*waux(166,moy,m+1)
        waux(223,n,m) = rzqc*waux(166,n,m) + rzwq*waux(166,n,m+1) + zmo*waux(166,moz,m+1)
        waux(224,n,m) = ryqc*waux(167,n,m) + rywq*waux(167,n,m+1) + ymo*waux(167,moy,m+1) + c1 
        waux(225,n,m) = rzqc*waux(167,n,m) + rzwq*waux(167,n,m+1) + zmo*waux(167,moz,m+1)
        waux(226,n,m) = rzqc*waux(168,n,m) + rzwq*waux(168,n,m+1) + zmo*waux(168,moz,m+1) + c1 
        waux(227,n,m) = ryqc*waux(169,n,m) + rywq*waux(169,n,m+1) + ymo*waux(169,moy,m+1) + c2 
        waux(228,n,m) = rzqc*waux(169,n,m) + rzwq*waux(169,n,m+1) + zmo*waux(169,moz,m+1)
        waux(229,n,m) = ryqc*waux(171,n,m) + rywq*waux(171,n,m+1) + ymo*waux(171,moy,m+1)
        waux(230,n,m) = rzqc*waux(171,n,m) + rzwq*waux(171,n,m+1) + zmo*waux(171,moz,m+1) + c3 
        waux(231,n,m) = ryqc*waux(172,n,m) + rywq*waux(172,n,m+1) + ymo*waux(172,moy,m+1) + c4*3.0d0 
        waux(232,n,m) = rzqc*waux(172,n,m) + rzwq*waux(172,n,m+1) + zmo*waux(172,moz,m+1)
        waux(233,n,m) = rzqc*waux(173,n,m) + rzwq*waux(173,n,m+1) + zmo*waux(173,moz,m+1) + c4 
        waux(234,n,m) = ryqc*waux(175,n,m) + rywq*waux(175,n,m+1) + ymo*waux(175,moy,m+1)
        waux(235,n,m) = rzqc*waux(175,n,m) + rzwq*waux(175,n,m+1) + zmo*waux(175,moz,m+1) + c5 
        waux(236,n,m) = ryqc*waux(176,n,m) + rywq*waux(176,n,m+1) + ymo*waux(176,moy,m+1) + c6*4.0d0 
        waux(237,n,m) = rzqc*waux(176,n,m) + rzwq*waux(176,n,m+1) + zmo*waux(176,moz,m+1)
        waux(238,n,m) = rzqc*waux(177,n,m) + rzwq*waux(177,n,m+1) + zmo*waux(177,moz,m+1) + c6 
        waux(239,n,m) = ryqc*waux(179,n,m) + rywq*waux(179,n,m+1) + ymo*waux(179,moy,m+1) + c7 
        waux(240,n,m) = ryqc*waux(180,n,m) + rywq*waux(180,n,m+1) + ymo*waux(180,moy,m+1)
        waux(241,n,m) = rzqc*waux(180,n,m) + rzwq*waux(180,n,m+1) + zmo*waux(180,moz,m+1) + c7*4.0d0 
        waux(242,n,m) = rxqc*waux(187,n,m) + rxwq*waux(187,n,m+1) + xmo*waux(187,mox,m+1) + c15*3.0d0 
        waux(243,n,m) = rzqc*waux(181,n,m) + rzwq*waux(181,n,m+1) + zmo*waux(181,moz,m+1)
        waux(244,n,m) = rzqc*waux(182,n,m) + rzwq*waux(182,n,m+1) + zmo*waux(182,moz,m+1) + c8 
        waux(245,n,m) = rzqc*waux(183,n,m) + rzwq*waux(183,n,m+1) + zmo*waux(183,moz,m+1) + c9 
        waux(246,n,m) = ryqc*waux(185,n,m) + rywq*waux(185,n,m+1) + ymo*waux(185,moy,m+1) + c10 
        waux(247,n,m) = ryqc*waux(186,n,m) + rywq*waux(186,n,m+1) + ymo*waux(186,moy,m+1)
        waux(248,n,m) = rxqc*waux(193,n,m) + rxwq*waux(193,n,m+1) + xmo*waux(193,mox,m+1) + c16*3.0d0 
        waux(249,n,m) = rxqc*waux(194,n,m) + rxwq*waux(194,n,m+1) + xmo*waux(194,mox,m+1) + c17 
        waux(250,n,m) = rzqc*waux(187,n,m) + rzwq*waux(187,n,m+1) + zmo*waux(187,moz,m+1)
        waux(251,n,m) = rzqc*waux(188,n,m) + rzwq*waux(188,n,m+1) + zmo*waux(188,moz,m+1) + c11 
        waux(252,n,m) = rzqc*waux(189,n,m) + rzwq*waux(189,n,m+1) + zmo*waux(189,moz,m+1) + c12 
        waux(253,n,m) = ryqc*waux(191,n,m) + rywq*waux(191,n,m+1) + ymo*waux(191,moy,m+1) + c13 
        waux(254,n,m) = ryqc*waux(192,n,m) + rywq*waux(192,n,m+1) + ymo*waux(192,moy,m+1) + c14 
        waux(255,n,m) = ryqc*waux(193,n,m) + rywq*waux(193,n,m+1) + ymo*waux(193,moy,m+1)
        waux(256,n,m) = rxqc*waux(201,n,m) + rxwq*waux(201,n,m+1) + xmo*waux(201,mox,m+1) + c18 
        waux(257,n,m) = rxqc*waux(202,n,m) + rxwq*waux(202,n,m+1) + xmo*waux(202,mox,m+1) + c19 
        waux(258,n,m) = rzqc*waux(194,n,m) + rzwq*waux(194,n,m+1) + zmo*waux(194,moz,m+1)
        waux(259,n,m) = rzqc*waux(195,n,m) + rzwq*waux(195,n,m+1) + zmo*waux(195,moz,m+1) + c15 
        waux(260,n,m) = rxqc*waux(205,n,m) + rxwq*waux(205,n,m+1) + xmo*waux(205,mox,m+1) + c22 
        waux(261,n,m) = rxqc*waux(206,n,m) + rxwq*waux(206,n,m+1) + xmo*waux(206,mox,m+1) + c23 
        waux(262,n,m) = rxqc*waux(207,n,m) + rxwq*waux(207,n,m+1) + xmo*waux(207,mox,m+1) + c24 
        waux(263,n,m) = ryqc*waux(200,n,m) + rywq*waux(200,n,m+1) + ymo*waux(200,moy,m+1) + c16 
        waux(264,n,m) = ryqc*waux(201,n,m) + rywq*waux(201,n,m+1) + ymo*waux(201,moy,m+1)
        waux(265,n,m) = rxqc*waux(210,n,m) + rxwq*waux(210,n,m+1) + xmo*waux(210,mox,m+1) + c27 
        waux(266,n,m) = rxqc*waux(211,n,m) + rxwq*waux(211,n,m+1) + xmo*waux(211,mox,m+1)
        waux(267,n,m) = rxqc*waux(212,n,m) + rxwq*waux(212,n,m+1) + xmo*waux(212,mox,m+1)
        waux(268,n,m) = rxqc*waux(213,n,m) + rxwq*waux(213,n,m+1) + xmo*waux(213,mox,m+1)
        waux(269,n,m) = rxqc*waux(214,n,m) + rxwq*waux(214,n,m+1) + xmo*waux(214,mox,m+1)
        waux(270,n,m) = rxqc*waux(215,n,m) + rxwq*waux(215,n,m+1) + xmo*waux(215,mox,m+1)
        waux(271,n,m) = rxqc*waux(216,n,m) + rxwq*waux(216,n,m+1) + xmo*waux(216,mox,m+1)
        waux(272,n,m) = rxqc*waux(217,n,m) + rxwq*waux(217,n,m+1) + xmo*waux(217,mox,m+1)
        waux(273,n,m) = rxqc*waux(218,n,m) + rxwq*waux(218,n,m+1) + xmo*waux(218,mox,m+1)
        waux(274,n,m) = rxqc*waux(219,n,m) + rxwq*waux(219,n,m+1) + xmo*waux(219,mox,m+1)
        waux(275,n,m) = rxqc*waux(220,n,m) + rxwq*waux(220,n,m+1) + xmo*waux(220,mox,m+1)
        waux(276,n,m) = ryqc*waux(211,n,m) + rywq*waux(211,n,m+1) + ymo*waux(211,moy,m+1) + c19*9.0d0 
        waux(277,n,m) = rzqc*waux(211,n,m) + rzwq*waux(211,n,m+1) + zmo*waux(211,moz,m+1)
        waux(278,n,m) = rzqc*waux(212,n,m) + rzwq*waux(212,n,m+1) + zmo*waux(212,moz,m+1) + c19 
        waux(279,n,m) = rzqc*waux(213,n,m) + rzwq*waux(213,n,m+1) + zmo*waux(213,moz,m+1) + c20 
        waux(280,n,m) = rzqc*waux(214,n,m) + rzwq*waux(214,n,m+1) + zmo*waux(214,moz,m+1) + c21 
        waux(281,n,m) = rzqc*waux(215,n,m) + rzwq*waux(215,n,m+1) + zmo*waux(215,moz,m+1) + c22*4.0d0 
        waux(282,n,m) = ryqc*waux(217,n,m) + rywq*waux(217,n,m+1) + ymo*waux(217,moy,m+1) + c25 
        waux(283,n,m) = ryqc*waux(218,n,m) + rywq*waux(218,n,m+1) + ymo*waux(218,moy,m+1) + c26 
        waux(284,n,m) = ryqc*waux(219,n,m) + rywq*waux(219,n,m+1) + ymo*waux(219,moy,m+1) + c27 
        waux(285,n,m) = ryqc*waux(220,n,m) + rywq*waux(220,n,m+1) + ymo*waux(220,moy,m+1)
        waux(286,n,m) = rzqc*waux(220,n,m) + rzwq*waux(220,n,m+1) + zmo*waux(220,moz,m+1) + c27*9.0d0 
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 286  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,286,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [11 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(166,n,m) - fnket*waux(166,n,m+1))
        c2 = f(2)*(waux(167,n,m) - fnket*waux(167,n,m+1))
        c3 = f(2)*(waux(168,n,m) - fnket*waux(168,n,m+1))
        c4 = f(1)*(waux(169,n,m) - fnket*waux(169,n,m+1))
        c5 = f(3)*(waux(171,n,m) - fnket*waux(171,n,m+1))
        c6 = f(1)*(waux(172,n,m) - fnket*waux(172,n,m+1))
        c7 = f(1)*(waux(175,n,m) - fnket*waux(175,n,m+1))
        c8 = f(1)*(waux(176,n,m) - fnket*waux(176,n,m+1))
        c9 = f(2)*(waux(177,n,m) - fnket*waux(177,n,m+1))
        c10= f(1)*(waux(180,n,m) - fnket*waux(180,n,m+1))
        c11= f(1)*(waux(181,n,m) - fnket*waux(181,n,m+1))
        c12= f(2)*(waux(182,n,m) - fnket*waux(182,n,m+1))
        c13= f(2)*(waux(185,n,m) - fnket*waux(185,n,m+1))
        c14= f(1)*(waux(186,n,m) - fnket*waux(186,n,m+1))
        c15= f(1)*(waux(187,n,m) - fnket*waux(187,n,m+1))
        c16= f(2)*(waux(188,n,m) - fnket*waux(188,n,m+1))
        c17= f(2)*(waux(192,n,m) - fnket*waux(192,n,m+1))
        c18= f(1)*(waux(193,n,m) - fnket*waux(193,n,m+1))
        c19= f(1)*(waux(194,n,m) - fnket*waux(194,n,m+1))
        c20= f(1)*(waux(201,n,m) - fnket*waux(201,n,m+1))
        c21= f(2)*(waux(202,n,m) - fnket*waux(202,n,m+1))
        c22= f(2)*(waux(206,n,m) - fnket*waux(206,n,m+1))
        c23= f(2)*(waux(210,n,m) - fnket*waux(210,n,m+1))
        c24= f(1)*(waux(211,n,m) - fnket*waux(211,n,m+1))
        c25= f(2)*(waux(212,n,m) - fnket*waux(212,n,m+1))
        c26= f(3)*(waux(213,n,m) - fnket*waux(213,n,m+1))
        c27= f(1)*(waux(214,n,m) - fnket*waux(214,n,m+1))
        c28= f(1)*(waux(215,n,m) - fnket*waux(215,n,m+1))
        c29= f(1)*(waux(216,n,m) - fnket*waux(216,n,m+1))
        c30= f(1)*(waux(217,n,m) - fnket*waux(217,n,m+1))
        c31= f(3)*(waux(218,n,m) - fnket*waux(218,n,m+1))
        c32= f(2)*(waux(219,n,m) - fnket*waux(219,n,m+1))
        c33= f(1)*(waux(220,n,m) - fnket*waux(220,n,m+1))
        
        waux(287,n,m) = rxqc*waux(221,n,m) + rxwq*waux(221,n,m+1) + xmo*waux(221,mox,m+1) + c1*10.0d0
        waux(288,n,m) = ryqc*waux(221,n,m) + rywq*waux(221,n,m+1) + ymo*waux(221,moy,m+1)
        waux(289,n,m) = rzqc*waux(221,n,m) + rzwq*waux(221,n,m+1) + zmo*waux(221,moz,m+1)
        waux(290,n,m) = ryqc*waux(222,n,m) + rywq*waux(222,n,m+1) + ymo*waux(222,moy,m+1) + c1 
        waux(291,n,m) = rzqc*waux(222,n,m) + rzwq*waux(222,n,m+1) + zmo*waux(222,moz,m+1)
        waux(292,n,m) = rzqc*waux(223,n,m) + rzwq*waux(223,n,m+1) + zmo*waux(223,moz,m+1) + c1
        waux(293,n,m) = ryqc*waux(224,n,m) + rywq*waux(224,n,m+1) + ymo*waux(224,moy,m+1) + c2
        waux(294,n,m) = rzqc*waux(224,n,m) + rzwq*waux(224,n,m+1) + zmo*waux(224,moz,m+1)
        waux(295,n,m) = ryqc*waux(226,n,m) + rywq*waux(226,n,m+1) + ymo*waux(226,moy,m+1)
        waux(296,n,m) = rzqc*waux(226,n,m) + rzwq*waux(226,n,m+1) + zmo*waux(226,moz,m+1) + c3 
        waux(297,n,m) = ryqc*waux(227,n,m) + rywq*waux(227,n,m+1) + ymo*waux(227,moy,m+1) + c4*3.0d0
        waux(298,n,m) = rzqc*waux(227,n,m) + rzwq*waux(227,n,m+1) + zmo*waux(227,moz,m+1)
        waux(299,n,m) = rzqc*waux(228,n,m) + rzwq*waux(228,n,m+1) + zmo*waux(228,moz,m+1) + c4 
        waux(300,n,m) = ryqc*waux(230,n,m) + rywq*waux(230,n,m+1) + ymo*waux(230,moy,m+1)
        waux(301,n,m) = rzqc*waux(230,n,m) + rzwq*waux(230,n,m+1) + zmo*waux(230,moz,m+1) + c5 
        waux(302,n,m) = ryqc*waux(231,n,m) + rywq*waux(231,n,m+1) + ymo*waux(231,moy,m+1) + c6*4.0d0 
        waux(303,n,m) = rzqc*waux(231,n,m) + rzwq*waux(231,n,m+1) + zmo*waux(231,moz,m+1)
        waux(304,n,m) = rzqc*waux(232,n,m) + rzwq*waux(232,n,m+1) + zmo*waux(232,moz,m+1) + c6 
        waux(305,n,m) = ryqc*waux(234,n,m) + rywq*waux(234,n,m+1) + ymo*waux(234,moy,m+1) + c7 
        waux(306,n,m) = ryqc*waux(235,n,m) + rywq*waux(235,n,m+1) + ymo*waux(235,moy,m+1)
        waux(307,n,m) = rzqc*waux(235,n,m) + rzwq*waux(235,n,m+1) + zmo*waux(235,moz,m+1) + c7*4.0d0 
        waux(308,n,m) = rxqc*waux(242,n,m) + rxwq*waux(242,n,m+1) + xmo*waux(242,mox,m+1) + c15*4.0d0 
        waux(309,n,m) = rzqc*waux(236,n,m) + rzwq*waux(236,n,m+1) + zmo*waux(236,moz,m+1)
        waux(310,n,m) = rzqc*waux(237,n,m) + rzwq*waux(237,n,m+1) + zmo*waux(237,moz,m+1) + c8 
        waux(311,n,m) = rzqc*waux(238,n,m) + rzwq*waux(238,n,m+1) + zmo*waux(238,moz,m+1) + c9 
        waux(312,n,m) = ryqc*waux(240,n,m) + rywq*waux(240,n,m+1) + ymo*waux(240,moy,m+1) + c10 
        waux(313,n,m) = ryqc*waux(241,n,m) + rywq*waux(241,n,m+1) + ymo*waux(241,moy,m+1)
        waux(314,n,m) = rxqc*waux(248,n,m) + rxwq*waux(248,n,m+1) + xmo*waux(248,mox,m+1) + c18*4.0d0 
        waux(315,n,m) = rxqc*waux(249,n,m) + rxwq*waux(249,n,m+1) + xmo*waux(249,mox,m+1) + c19*3.0d0 
        waux(316,n,m) = rzqc*waux(242,n,m) + rzwq*waux(242,n,m+1) + zmo*waux(242,moz,m+1)
        waux(317,n,m) = rzqc*waux(243,n,m) + rzwq*waux(243,n,m+1) + zmo*waux(243,moz,m+1) + c11 
        waux(318,n,m) = rzqc*waux(244,n,m) + rzwq*waux(244,n,m+1) + zmo*waux(244,moz,m+1) + c12 
        waux(319,n,m) = ryqc*waux(246,n,m) + rywq*waux(246,n,m+1) + ymo*waux(246,moy,m+1) + c13 
        waux(320,n,m) = ryqc*waux(247,n,m) + rywq*waux(247,n,m+1) + ymo*waux(247,moy,m+1) + c14 
        waux(321,n,m) = ryqc*waux(248,n,m) + rywq*waux(248,n,m+1) + ymo*waux(248,moy,m+1)
        waux(322,n,m) = rxqc*waux(256,n,m) + rxwq*waux(256,n,m+1) + xmo*waux(256,mox,m+1) + c20*3.0d0 
        waux(323,n,m) = rxqc*waux(257,n,m) + rxwq*waux(257,n,m+1) + xmo*waux(257,mox,m+1) + c21 
        waux(324,n,m) = rzqc*waux(249,n,m) + rzwq*waux(249,n,m+1) + zmo*waux(249,moz,m+1)
        waux(325,n,m) = rzqc*waux(250,n,m) + rzwq*waux(250,n,m+1) + zmo*waux(250,moz,m+1) + c15 
        waux(326,n,m) = rzqc*waux(251,n,m) + rzwq*waux(251,n,m+1) + zmo*waux(251,moz,m+1) + c16 
        waux(327,n,m) = rxqc*waux(261,n,m) + rxwq*waux(261,n,m+1) + xmo*waux(261,mox,m+1) + c22 
        waux(328,n,m) = ryqc*waux(254,n,m) + rywq*waux(254,n,m+1) + ymo*waux(254,moy,m+1) + c17 
        waux(329,n,m) = ryqc*waux(255,n,m) + rywq*waux(255,n,m+1) + ymo*waux(255,moy,m+1) + c18 
        waux(330,n,m) = ryqc*waux(256,n,m) + rywq*waux(256,n,m+1) + ymo*waux(256,moy,m+1)
        waux(331,n,m) = rxqc*waux(265,n,m) + rxwq*waux(265,n,m+1) + xmo*waux(265,mox,m+1) + c23 
        waux(332,n,m) = rxqc*waux(266,n,m) + rxwq*waux(266,n,m+1) + xmo*waux(266,mox,m+1) + c24 
        waux(333,n,m) = rzqc*waux(257,n,m) + rzwq*waux(257,n,m+1) + zmo*waux(257,moz,m+1)
        waux(334,n,m) = rzqc*waux(258,n,m) + rzwq*waux(258,n,m+1) + zmo*waux(258,moz,m+1) + c19 
        waux(335,n,m) = rxqc*waux(269,n,m) + rxwq*waux(269,n,m+1) + xmo*waux(269,mox,m+1) + c27 
        waux(336,n,m) = rxqc*waux(270,n,m) + rxwq*waux(270,n,m+1) + xmo*waux(270,mox,m+1) + c28 
        waux(337,n,m) = rxqc*waux(271,n,m) + rxwq*waux(271,n,m+1) + xmo*waux(271,mox,m+1) + c29 
        waux(338,n,m) = rxqc*waux(272,n,m) + rxwq*waux(272,n,m+1) + xmo*waux(272,mox,m+1) + c30 
        waux(339,n,m) = ryqc*waux(264,n,m) + rywq*waux(264,n,m+1) + ymo*waux(264,moy,m+1) + c20 
        waux(340,n,m) = ryqc*waux(265,n,m) + rywq*waux(265,n,m+1) + ymo*waux(265,moy,m+1)
        waux(341,n,m) = rxqc*waux(275,n,m) + rxwq*waux(275,n,m+1) + xmo*waux(275,mox,m+1) + c33 
        waux(342,n,m) = rxqc*waux(276,n,m) + rxwq*waux(276,n,m+1) + xmo*waux(276,mox,m+1)
        waux(343,n,m) = rxqc*waux(277,n,m) + rxwq*waux(277,n,m+1) + xmo*waux(277,mox,m+1)
        waux(344,n,m) = rxqc*waux(278,n,m) + rxwq*waux(278,n,m+1) + xmo*waux(278,mox,m+1)
        waux(345,n,m) = rxqc*waux(279,n,m) + rxwq*waux(279,n,m+1) + xmo*waux(279,mox,m+1)
        waux(346,n,m) = rxqc*waux(280,n,m) + rxwq*waux(280,n,m+1) + xmo*waux(280,mox,m+1)
        waux(347,n,m) = rxqc*waux(281,n,m) + rxwq*waux(281,n,m+1) + xmo*waux(281,mox,m+1)
        waux(348,n,m) = rxqc*waux(282,n,m) + rxwq*waux(282,n,m+1) + xmo*waux(282,mox,m+1)
        waux(349,n,m) = rxqc*waux(283,n,m) + rxwq*waux(283,n,m+1) + xmo*waux(283,mox,m+1)
        waux(350,n,m) = rxqc*waux(284,n,m) + rxwq*waux(284,n,m+1) + xmo*waux(284,mox,m+1)
        waux(351,n,m) = rxqc*waux(285,n,m) + rxwq*waux(285,n,m+1) + xmo*waux(285,mox,m+1)
        waux(352,n,m) = rxqc*waux(286,n,m) + rxwq*waux(286,n,m+1) + xmo*waux(286,mox,m+1)
        waux(353,n,m) = ryqc*waux(276,n,m) + rywq*waux(276,n,m+1) + ymo*waux(276,moy,m+1) + c24*10.0d0
        waux(354,n,m) = rzqc*waux(276,n,m) + rzwq*waux(276,n,m+1) + zmo*waux(276,moz,m+1)
        waux(355,n,m) = rzqc*waux(277,n,m) + rzwq*waux(277,n,m+1) + zmo*waux(277,moz,m+1) + c24 
        waux(356,n,m) = rzqc*waux(278,n,m) + rzwq*waux(278,n,m+1) + zmo*waux(278,moz,m+1) + c25 
        waux(357,n,m) = rzqc*waux(279,n,m) + rzwq*waux(279,n,m+1) + zmo*waux(279,moz,m+1) + c26 
        waux(358,n,m) = rzqc*waux(280,n,m) + rzwq*waux(280,n,m+1) + zmo*waux(280,moz,m+1) + c27*4.0d0 
        waux(359,n,m) = ryqc*waux(282,n,m) + rywq*waux(282,n,m+1) + ymo*waux(282,moy,m+1) + c30*4.0d0 
        waux(360,n,m) = ryqc*waux(283,n,m) + rywq*waux(283,n,m+1) + ymo*waux(283,moy,m+1) + c31 
        waux(361,n,m) = ryqc*waux(284,n,m) + rywq*waux(284,n,m+1) + ymo*waux(284,moy,m+1) + c32 
        waux(362,n,m) = ryqc*waux(285,n,m) + rywq*waux(285,n,m+1) + ymo*waux(285,moy,m+1) + c33 
        waux(363,n,m) = ryqc*waux(286,n,m) + rywq*waux(286,n,m+1) + ymo*waux(286,moy,m+1)
        waux(364,n,m) = rzqc*waux(286,n,m) + rzwq*waux(286,n,m+1) + zmo*waux(286,moz,m+1) + c33*10.0d0
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 364  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,364,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [12 | * ] --!
    mmax=mmax-1
    do m=1,mmax
      do n=1,nmax
        mox =      mo(1,n)
        moy =      mo(3,n)
        moz =      mo(5,n)
        xmo = gama(mo(2,n))
        ymo = gama(mo(4,n))
        zmo = gama(mo(6,n))

        c1 = f(1)*(waux(221,n,m) - fnket*waux(221,n,m+1))
        c2 = f(2)*(waux(222,n,m) - fnket*waux(222,n,m+1))
        c3 = f(2)*(waux(223,n,m) - fnket*waux(223,n,m+1))
        c4 = f(1)*(waux(224,n,m) - fnket*waux(224,n,m+1))
        c5 = f(3)*(waux(226,n,m) - fnket*waux(226,n,m+1))
        c6 = f(1)*(waux(227,n,m) - fnket*waux(227,n,m+1))
        c7 = f(1)*(waux(230,n,m) - fnket*waux(230,n,m+1))
        c8 = f(1)*(waux(231,n,m) - fnket*waux(231,n,m+1))
        c9 = f(2)*(waux(232,n,m) - fnket*waux(232,n,m+1))
        c10= f(1)*(waux(235,n,m) - fnket*waux(235,n,m+1))
        c11= f(1)*(waux(236,n,m) - fnket*waux(236,n,m+1))
        c12= f(2)*(waux(237,n,m) - fnket*waux(237,n,m+1))
        c13= f(2)*(waux(240,n,m) - fnket*waux(240,n,m+1))
        c14= f(1)*(waux(241,n,m) - fnket*waux(241,n,m+1))
        c15= f(1)*(waux(242,n,m) - fnket*waux(242,n,m+1))
        c16= f(2)*(waux(243,n,m) - fnket*waux(243,n,m+1))
        c17= f(3)*(waux(244,n,m) - fnket*waux(244,n,m+1))
        c18= f(2)*(waux(247,n,m) - fnket*waux(247,n,m+1))
        c19= f(1)*(waux(248,n,m) - fnket*waux(248,n,m+1))
        c20= f(1)*(waux(249,n,m) - fnket*waux(249,n,m+1))
        c21= f(2)*(waux(250,n,m) - fnket*waux(250,n,m+1))
        c22= f(2)*(waux(255,n,m) - fnket*waux(255,n,m+1))
        c23= f(1)*(waux(256,n,m) - fnket*waux(256,n,m+1))
        c24= f(1)*(waux(257,n,m) - fnket*waux(257,n,m+1))
        c25= f(1)*(waux(265,n,m) - fnket*waux(265,n,m+1))
        c26= f(2)*(waux(266,n,m) - fnket*waux(266,n,m+1))
        c27= f(2)*(waux(270,n,m) - fnket*waux(270,n,m+1))
        c28= f(2)*(waux(271,n,m) - fnket*waux(271,n,m+1))
        c29= f(2)*(waux(275,n,m) - fnket*waux(275,n,m+1))
        c30= f(1)*(waux(276,n,m) - fnket*waux(276,n,m+1))
        c31= f(2)*(waux(277,n,m) - fnket*waux(277,n,m+1))
        c32= f(3)*(waux(278,n,m) - fnket*waux(278,n,m+1))
        c33= f(1)*(waux(279,n,m) - fnket*waux(279,n,m+1))
        c34= f(1)*(waux(280,n,m) - fnket*waux(280,n,m+1))
        c35= f(1)*(waux(281,n,m) - fnket*waux(281,n,m+1))
        c36= f(1)*(waux(282,n,m) - fnket*waux(282,n,m+1))
        c37= f(1)*(waux(283,n,m) - fnket*waux(283,n,m+1))
        c38= f(3)*(waux(284,n,m) - fnket*waux(284,n,m+1))
        c39= f(2)*(waux(285,n,m) - fnket*waux(285,n,m+1))
        c40= f(1)*(waux(286,n,m) - fnket*waux(286,n,m+1))
        
        waux(365,n,m) = rxqc*waux(287,n,m) + rxwq*waux(287,n,m+1) + xmo*waux(287,mox,m+1) + c1*11.0d0
        waux(366,n,m) = ryqc*waux(287,n,m) + rywq*waux(287,n,m+1) + ymo*waux(287,moy,m+1)
        waux(367,n,m) = rzqc*waux(287,n,m) + rzwq*waux(287,n,m+1) + zmo*waux(287,moz,m+1)
        waux(368,n,m) = ryqc*waux(288,n,m) + rywq*waux(288,n,m+1) + ymo*waux(288,moy,m+1) + c1 
        waux(369,n,m) = rzqc*waux(288,n,m) + rzwq*waux(288,n,m+1) + zmo*waux(288,moz,m+1)
        waux(370,n,m) = rzqc*waux(289,n,m) + rzwq*waux(289,n,m+1) + zmo*waux(289,moz,m+1) + c1
        waux(371,n,m) = ryqc*waux(290,n,m) + rywq*waux(290,n,m+1) + ymo*waux(290,moy,m+1) + c2 
        waux(372,n,m) = rzqc*waux(290,n,m) + rzwq*waux(290,n,m+1) + zmo*waux(290,moz,m+1)
        waux(373,n,m) = ryqc*waux(292,n,m) + rywq*waux(292,n,m+1) + ymo*waux(292,moy,m+1)
        waux(374,n,m) = rzqc*waux(292,n,m) + rzwq*waux(292,n,m+1) + zmo*waux(292,moz,m+1) + c3 
        waux(375,n,m) = ryqc*waux(293,n,m) + rywq*waux(293,n,m+1) + ymo*waux(293,moy,m+1) + c4*3.0d0 
        waux(376,n,m) = rzqc*waux(293,n,m) + rzwq*waux(293,n,m+1) + zmo*waux(293,moz,m+1)
        waux(377,n,m) = rzqc*waux(294,n,m) + rzwq*waux(294,n,m+1) + zmo*waux(294,moz,m+1) + c4 
        waux(378,n,m) = ryqc*waux(296,n,m) + rywq*waux(296,n,m+1) + ymo*waux(296,moy,m+1)
        waux(379,n,m) = rzqc*waux(296,n,m) + rzwq*waux(296,n,m+1) + zmo*waux(296,moz,m+1) + c5 
        waux(380,n,m) = ryqc*waux(297,n,m) + rywq*waux(297,n,m+1) + ymo*waux(297,moy,m+1) + c6*4.0d0 
        waux(381,n,m) = rzqc*waux(297,n,m) + rzwq*waux(297,n,m+1) + zmo*waux(297,moz,m+1)
        waux(382,n,m) = rzqc*waux(298,n,m) + rzwq*waux(298,n,m+1) + zmo*waux(298,moz,m+1) + c6 
        waux(383,n,m) = ryqc*waux(300,n,m) + rywq*waux(300,n,m+1) + ymo*waux(300,moy,m+1) + c7 
        waux(384,n,m) = ryqc*waux(301,n,m) + rywq*waux(301,n,m+1) + ymo*waux(301,moy,m+1)
        waux(385,n,m) = rzqc*waux(301,n,m) + rzwq*waux(301,n,m+1) + zmo*waux(301,moz,m+1) + c7*4.0d0 
        waux(386,n,m) = ryqc*waux(302,n,m) + rywq*waux(302,n,m+1) + ymo*waux(302,moy,m+1) + c8*5.0d0 
        waux(387,n,m) = rzqc*waux(302,n,m) + rzwq*waux(302,n,m+1) + zmo*waux(302,moz,m+1)
        waux(388,n,m) = rzqc*waux(303,n,m) + rzwq*waux(303,n,m+1) + zmo*waux(303,moz,m+1) + c8 
        waux(389,n,m) = rzqc*waux(304,n,m) + rzwq*waux(304,n,m+1) + zmo*waux(304,moz,m+1) + c9 
        waux(390,n,m) = ryqc*waux(306,n,m) + rywq*waux(306,n,m+1) + ymo*waux(306,moy,m+1) + c10 
        waux(391,n,m) = ryqc*waux(307,n,m) + rywq*waux(307,n,m+1) + ymo*waux(307,moy,m+1)
        waux(392,n,m) = rzqc*waux(307,n,m) + rzwq*waux(307,n,m+1) + zmo*waux(307,moz,m+1) + c10*5.0d0 
        waux(393,n,m) = rxqc*waux(315,n,m) + rxwq*waux(315,n,m+1) + xmo*waux(315,mox,m+1) + c20*4.0d0 
        waux(394,n,m) = rzqc*waux(308,n,m) + rzwq*waux(308,n,m+1) + zmo*waux(308,moz,m+1)
        waux(395,n,m) = rzqc*waux(309,n,m) + rzwq*waux(309,n,m+1) + zmo*waux(309,moz,m+1) + c11 
        waux(396,n,m) = rzqc*waux(310,n,m) + rzwq*waux(310,n,m+1) + zmo*waux(310,moz,m+1) + c12 
        waux(397,n,m) = ryqc*waux(312,n,m) + rywq*waux(312,n,m+1) + ymo*waux(312,moy,m+1) + c13 
        waux(398,n,m) = ryqc*waux(313,n,m) + rywq*waux(313,n,m+1) + ymo*waux(313,moy,m+1) + c14 
        waux(399,n,m) = ryqc*waux(314,n,m) + rywq*waux(314,n,m+1) + ymo*waux(314,moy,m+1)
        waux(400,n,m) = rxqc*waux(322,n,m) + rxwq*waux(322,n,m+1) + xmo*waux(322,mox,m+1) + c23*4.0d0 
        waux(401,n,m) = rxqc*waux(323,n,m) + rxwq*waux(323,n,m+1) + xmo*waux(323,mox,m+1) + c24*3.0d0 
        waux(402,n,m) = rzqc*waux(315,n,m) + rzwq*waux(315,n,m+1) + zmo*waux(315,moz,m+1)
        waux(403,n,m) = rzqc*waux(316,n,m) + rzwq*waux(316,n,m+1) + zmo*waux(316,moz,m+1) + c15 
        waux(404,n,m) = rzqc*waux(317,n,m) + rzwq*waux(317,n,m+1) + zmo*waux(317,moz,m+1) + c16 
        waux(405,n,m) = rzqc*waux(318,n,m) + rzwq*waux(318,n,m+1) + zmo*waux(318,moz,m+1) + c17 
        waux(406,n,m) = ryqc*waux(320,n,m) + rywq*waux(320,n,m+1) + ymo*waux(320,moy,m+1) + c18 
        waux(407,n,m) = ryqc*waux(321,n,m) + rywq*waux(321,n,m+1) + ymo*waux(321,moy,m+1) + c19 
        waux(408,n,m) = ryqc*waux(322,n,m) + rywq*waux(322,n,m+1) + ymo*waux(322,moy,m+1)
        waux(409,n,m) = rxqc*waux(331,n,m) + rxwq*waux(331,n,m+1) + xmo*waux(331,mox,m+1) + c25*3.0d0 
        waux(410,n,m) = rxqc*waux(332,n,m) + rxwq*waux(332,n,m+1) + xmo*waux(332,mox,m+1) + c26 
        waux(411,n,m) = rzqc*waux(323,n,m) + rzwq*waux(323,n,m+1) + zmo*waux(323,moz,m+1)
        waux(412,n,m) = rzqc*waux(324,n,m) + rzwq*waux(324,n,m+1) + zmo*waux(324,moz,m+1) + c20 
        waux(413,n,m) = rzqc*waux(325,n,m) + rzwq*waux(325,n,m+1) + zmo*waux(325,moz,m+1) + c21 
        waux(414,n,m) = rxqc*waux(336,n,m) + rxwq*waux(336,n,m+1) + xmo*waux(336,mox,m+1) + c27 
        waux(415,n,m) = rxqc*waux(337,n,m) + rxwq*waux(337,n,m+1) + xmo*waux(337,mox,m+1) + c28 
        waux(416,n,m) = ryqc*waux(329,n,m) + rywq*waux(329,n,m+1) + ymo*waux(329,moy,m+1) + c22 
        waux(417,n,m) = ryqc*waux(330,n,m) + rywq*waux(330,n,m+1) + ymo*waux(330,moy,m+1) + c23 
        waux(418,n,m) = ryqc*waux(331,n,m) + rywq*waux(331,n,m+1) + ymo*waux(331,moy,m+1)
        waux(419,n,m) = rxqc*waux(341,n,m) + rxwq*waux(341,n,m+1) + xmo*waux(341,mox,m+1) + c29 
        waux(420,n,m) = rxqc*waux(342,n,m) + rxwq*waux(342,n,m+1) + xmo*waux(342,mox,m+1) + c30 
        waux(421,n,m) = rzqc*waux(332,n,m) + rzwq*waux(332,n,m+1) + zmo*waux(332,moz,m+1)
        waux(422,n,m) = rzqc*waux(333,n,m) + rzwq*waux(333,n,m+1) + zmo*waux(333,moz,m+1) + c24 
        waux(423,n,m) = rxqc*waux(345,n,m) + rxwq*waux(345,n,m+1) + xmo*waux(345,mox,m+1) + c33 
        waux(424,n,m) = rxqc*waux(346,n,m) + rxwq*waux(346,n,m+1) + xmo*waux(346,mox,m+1) + c34 
        waux(425,n,m) = rxqc*waux(347,n,m) + rxwq*waux(347,n,m+1) + xmo*waux(347,mox,m+1) + c35 
        waux(426,n,m) = rxqc*waux(348,n,m) + rxwq*waux(348,n,m+1) + xmo*waux(348,mox,m+1) + c36 
        waux(427,n,m) = rxqc*waux(349,n,m) + rxwq*waux(349,n,m+1) + xmo*waux(349,mox,m+1) + c37 
        waux(428,n,m) = ryqc*waux(340,n,m) + rywq*waux(340,n,m+1) + ymo*waux(340,moy,m+1) + c25 
        waux(429,n,m) = ryqc*waux(341,n,m) + rywq*waux(341,n,m+1) + ymo*waux(341,moy,m+1)
        waux(430,n,m) = rxqc*waux(352,n,m) + rxwq*waux(352,n,m+1) + xmo*waux(352,mox,m+1) + c40 
        waux(431,n,m) = rxqc*waux(353,n,m) + rxwq*waux(353,n,m+1) + xmo*waux(353,mox,m+1)
        waux(432,n,m) = rxqc*waux(354,n,m) + rxwq*waux(354,n,m+1) + xmo*waux(354,mox,m+1)
        waux(433,n,m) = rxqc*waux(355,n,m) + rxwq*waux(355,n,m+1) + xmo*waux(355,mox,m+1)
        waux(434,n,m) = rxqc*waux(356,n,m) + rxwq*waux(356,n,m+1) + xmo*waux(356,mox,m+1)
        waux(435,n,m) = rxqc*waux(357,n,m) + rxwq*waux(357,n,m+1) + xmo*waux(357,mox,m+1)
        waux(436,n,m) = rxqc*waux(358,n,m) + rxwq*waux(358,n,m+1) + xmo*waux(358,mox,m+1)
        waux(437,n,m) = rxqc*waux(359,n,m) + rxwq*waux(359,n,m+1) + xmo*waux(359,mox,m+1)
        waux(438,n,m) = rxqc*waux(360,n,m) + rxwq*waux(360,n,m+1) + xmo*waux(360,mox,m+1)
        waux(439,n,m) = rxqc*waux(361,n,m) + rxwq*waux(361,n,m+1) + xmo*waux(361,mox,m+1)
        waux(440,n,m) = rxqc*waux(362,n,m) + rxwq*waux(362,n,m+1) + xmo*waux(362,mox,m+1)
        waux(441,n,m) = rxqc*waux(363,n,m) + rxwq*waux(363,n,m+1) + xmo*waux(363,mox,m+1)
        waux(442,n,m) = rxqc*waux(364,n,m) + rxwq*waux(364,n,m+1) + xmo*waux(364,mox,m+1)
        waux(443,n,m) = ryqc*waux(353,n,m) + rywq*waux(353,n,m+1) + ymo*waux(353,moy,m+1) + c30*11.0d0
        waux(444,n,m) = rzqc*waux(353,n,m) + rzwq*waux(353,n,m+1) + zmo*waux(353,moz,m+1)
        waux(445,n,m) = rzqc*waux(354,n,m) + rzwq*waux(354,n,m+1) + zmo*waux(354,moz,m+1) + c30 
        waux(446,n,m) = rzqc*waux(355,n,m) + rzwq*waux(355,n,m+1) + zmo*waux(355,moz,m+1) + c31 
        waux(447,n,m) = rzqc*waux(356,n,m) + rzwq*waux(356,n,m+1) + zmo*waux(356,moz,m+1) + c32 
        waux(448,n,m) = rzqc*waux(357,n,m) + rzwq*waux(357,n,m+1) + zmo*waux(357,moz,m+1) + c33*4.0d0 
        waux(449,n,m) = rzqc*waux(358,n,m) + rzwq*waux(358,n,m+1) + zmo*waux(358,moz,m+1) + c34*5.0d0 
        waux(450,n,m) = ryqc*waux(360,n,m) + rywq*waux(360,n,m+1) + ymo*waux(360,moy,m+1) + c37*4.0d0 
        waux(451,n,m) = ryqc*waux(361,n,m) + rywq*waux(361,n,m+1) + ymo*waux(361,moy,m+1) + c38 
        waux(452,n,m) = ryqc*waux(362,n,m) + rywq*waux(362,n,m+1) + ymo*waux(362,moy,m+1) + c39 
        waux(453,n,m) = ryqc*waux(363,n,m) + rywq*waux(363,n,m+1) + ymo*waux(363,moy,m+1) + c40 
        waux(454,n,m) = ryqc*waux(364,n,m) + rywq*waux(364,n,m+1) + ymo*waux(364,moy,m+1)
        waux(455,n,m) = rzqc*waux(364,n,m) + rzwq*waux(364,n,m+1) + zmo*waux(364,moz,m+1) + c40*11.0d0
      enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 455  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,455,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [13 | * ] --!
    f( 1) = fdket
    f( 2) = f( 1) + f(1) 
    f( 3) = f( 2) + f(1) 
    f( 4) = f( 3) + f(1) 
    f( 5) = f( 4) + f(1) 
    f( 6) = f( 5) + f(1) 
    f( 7) = f( 6) + f(1) 
    f( 8) = f( 7) + f(1) 
    f( 9) = f( 8) + f(1) 
    f(10) = f( 9) + f(1) 
    f(11) = f(10) + f(1) 
    f(12) = f(11) + f(1) 
    mmax=mmax-1
    do m=1,mmax
    do n=1,nmax
    mox =      mo(1,n)
    moy =      mo(3,n)
    moz =      mo(5,n)
    xmo = gama(mo(2,n))
    ymo = gama(mo(4,n))
    zmo = gama(mo(6,n))
    waux(456,n,m)=rxqc*waux(365,n,m)+rxwq*waux(365,n,m+1)+xmo*waux(365,mox,m+1)+f(12)*(waux(287,n,m)-fnket*waux(287,n,m+1))
    waux(457,n,m)=ryqc*waux(365,n,m)+rywq*waux(365,n,m+1)+ymo*waux(365,moy,m+1)
    waux(458,n,m)=rzqc*waux(365,n,m)+rzwq*waux(365,n,m+1)+zmo*waux(365,moz,m+1)
    waux(459,n,m)=ryqc*waux(366,n,m)+rywq*waux(366,n,m+1)+ymo*waux(366,moy,m+1)+f( 1)*(waux(287,n,m)-fnket*waux(287,n,m+1))
    waux(460,n,m)=rzqc*waux(366,n,m)+rzwq*waux(366,n,m+1)+zmo*waux(366,moz,m+1)
    waux(461,n,m)=rzqc*waux(367,n,m)+rzwq*waux(367,n,m+1)+zmo*waux(367,moz,m+1)+f( 1)*(waux(287,n,m)-fnket*waux(287,n,m+1))
    waux(462,n,m)=ryqc*waux(368,n,m)+rywq*waux(368,n,m+1)+ymo*waux(368,moy,m+1)+f( 2)*(waux(288,n,m)-fnket*waux(288,n,m+1))
    waux(463,n,m)=rzqc*waux(368,n,m)+rzwq*waux(368,n,m+1)+zmo*waux(368,moz,m+1)
    waux(464,n,m)=ryqc*waux(370,n,m)+rywq*waux(370,n,m+1)+ymo*waux(370,moy,m+1)
    waux(465,n,m)=rzqc*waux(370,n,m)+rzwq*waux(370,n,m+1)+zmo*waux(370,moz,m+1)+f( 2)*(waux(289,n,m)-fnket*waux(289,n,m+1))
    waux(466,n,m)=ryqc*waux(371,n,m)+rywq*waux(371,n,m+1)+ymo*waux(371,moy,m+1)+f( 3)*(waux(290,n,m)-fnket*waux(290,n,m+1))
    waux(467,n,m)=rzqc*waux(371,n,m)+rzwq*waux(371,n,m+1)+zmo*waux(371,moz,m+1)
    waux(468,n,m)=rzqc*waux(372,n,m)+rzwq*waux(372,n,m+1)+zmo*waux(372,moz,m+1)+f( 1)*(waux(290,n,m)-fnket*waux(290,n,m+1))
    waux(469,n,m)=ryqc*waux(374,n,m)+rywq*waux(374,n,m+1)+ymo*waux(374,moy,m+1)
    waux(470,n,m)=rzqc*waux(374,n,m)+rzwq*waux(374,n,m+1)+zmo*waux(374,moz,m+1)+f( 3)*(waux(292,n,m)-fnket*waux(292,n,m+1))
    waux(471,n,m)=ryqc*waux(375,n,m)+rywq*waux(375,n,m+1)+ymo*waux(375,moy,m+1)+f( 4)*(waux(293,n,m)-fnket*waux(293,n,m+1))
    waux(472,n,m)=rzqc*waux(375,n,m)+rzwq*waux(375,n,m+1)+zmo*waux(375,moz,m+1)
    waux(473,n,m)=rzqc*waux(376,n,m)+rzwq*waux(376,n,m+1)+zmo*waux(376,moz,m+1)+f( 1)*(waux(293,n,m)-fnket*waux(293,n,m+1))
    waux(474,n,m)=ryqc*waux(378,n,m)+rywq*waux(378,n,m+1)+ymo*waux(378,moy,m+1)+f( 1)*(waux(296,n,m)-fnket*waux(296,n,m+1))
    waux(475,n,m)=ryqc*waux(379,n,m)+rywq*waux(379,n,m+1)+ymo*waux(379,moy,m+1)
    waux(476,n,m)=rzqc*waux(379,n,m)+rzwq*waux(379,n,m+1)+zmo*waux(379,moz,m+1)+f( 4)*(waux(296,n,m)-fnket*waux(296,n,m+1))
    waux(477,n,m)=ryqc*waux(380,n,m)+rywq*waux(380,n,m+1)+ymo*waux(380,moy,m+1)+f( 5)*(waux(297,n,m)-fnket*waux(297,n,m+1))
    waux(478,n,m)=rzqc*waux(380,n,m)+rzwq*waux(380,n,m+1)+zmo*waux(380,moz,m+1)
    waux(479,n,m)=rzqc*waux(381,n,m)+rzwq*waux(381,n,m+1)+zmo*waux(381,moz,m+1)+f( 1)*(waux(297,n,m)-fnket*waux(297,n,m+1))
    waux(480,n,m)=rzqc*waux(382,n,m)+rzwq*waux(382,n,m+1)+zmo*waux(382,moz,m+1)+f( 2)*(waux(298,n,m)-fnket*waux(298,n,m+1))
    waux(481,n,m)=ryqc*waux(384,n,m)+rywq*waux(384,n,m+1)+ymo*waux(384,moy,m+1)+f( 1)*(waux(301,n,m)-fnket*waux(301,n,m+1))
    waux(482,n,m)=ryqc*waux(385,n,m)+rywq*waux(385,n,m+1)+ymo*waux(385,moy,m+1)
    waux(483,n,m)=rzqc*waux(385,n,m)+rzwq*waux(385,n,m+1)+zmo*waux(385,moz,m+1)+f( 5)*(waux(301,n,m)-fnket*waux(301,n,m+1))
    waux(484,n,m)=rxqc*waux(393,n,m)+rxwq*waux(393,n,m+1)+xmo*waux(393,mox,m+1)+f( 5)*(waux(315,n,m)-fnket*waux(315,n,m+1))
    waux(485,n,m)=rzqc*waux(386,n,m)+rzwq*waux(386,n,m+1)+zmo*waux(386,moz,m+1)
    waux(486,n,m)=rzqc*waux(387,n,m)+rzwq*waux(387,n,m+1)+zmo*waux(387,moz,m+1)+f( 1)*(waux(302,n,m)-fnket*waux(302,n,m+1))
    waux(487,n,m)=rzqc*waux(388,n,m)+rzwq*waux(388,n,m+1)+zmo*waux(388,moz,m+1)+f( 2)*(waux(303,n,m)-fnket*waux(303,n,m+1))
    waux(488,n,m)=ryqc*waux(390,n,m)+rywq*waux(390,n,m+1)+ymo*waux(390,moy,m+1)+f( 2)*(waux(306,n,m)-fnket*waux(306,n,m+1))
    waux(489,n,m)=ryqc*waux(391,n,m)+rywq*waux(391,n,m+1)+ymo*waux(391,moy,m+1)+f( 1)*(waux(307,n,m)-fnket*waux(307,n,m+1))
    waux(490,n,m)=ryqc*waux(392,n,m)+rywq*waux(392,n,m+1)+ymo*waux(392,moy,m+1)
    waux(491,n,m)=rxqc*waux(400,n,m)+rxwq*waux(400,n,m+1)+xmo*waux(400,mox,m+1)+f( 5)*(waux(322,n,m)-fnket*waux(322,n,m+1))
    waux(492,n,m)=rxqc*waux(401,n,m)+rxwq*waux(401,n,m+1)+xmo*waux(401,mox,m+1)+f( 4)*(waux(323,n,m)-fnket*waux(323,n,m+1))
    waux(493,n,m)=rzqc*waux(393,n,m)+rzwq*waux(393,n,m+1)+zmo*waux(393,moz,m+1)
    waux(494,n,m)=rzqc*waux(394,n,m)+rzwq*waux(394,n,m+1)+zmo*waux(394,moz,m+1)+f( 1)*(waux(308,n,m)-fnket*waux(308,n,m+1))
    waux(495,n,m)=rzqc*waux(395,n,m)+rzwq*waux(395,n,m+1)+zmo*waux(395,moz,m+1)+f( 2)*(waux(309,n,m)-fnket*waux(309,n,m+1))
    waux(496,n,m)=rzqc*waux(396,n,m)+rzwq*waux(396,n,m+1)+zmo*waux(396,moz,m+1)+f( 3)*(waux(310,n,m)-fnket*waux(310,n,m+1))
    waux(497,n,m)=ryqc*waux(398,n,m)+rywq*waux(398,n,m+1)+ymo*waux(398,moy,m+1)+f( 2)*(waux(313,n,m)-fnket*waux(313,n,m+1))
    waux(498,n,m)=ryqc*waux(399,n,m)+rywq*waux(399,n,m+1)+ymo*waux(399,moy,m+1)+f( 1)*(waux(314,n,m)-fnket*waux(314,n,m+1))
    waux(499,n,m)=ryqc*waux(400,n,m)+rywq*waux(400,n,m+1)+ymo*waux(400,moy,m+1)
    waux(500,n,m)=rxqc*waux(409,n,m)+rxwq*waux(409,n,m+1)+xmo*waux(409,mox,m+1)+f( 4)*(waux(331,n,m)-fnket*waux(331,n,m+1))
    waux(501,n,m)=rxqc*waux(410,n,m)+rxwq*waux(410,n,m+1)+xmo*waux(410,mox,m+1)+f( 3)*(waux(332,n,m)-fnket*waux(332,n,m+1))
    waux(502,n,m)=rzqc*waux(401,n,m)+rzwq*waux(401,n,m+1)+zmo*waux(401,moz,m+1)
    waux(503,n,m)=rzqc*waux(402,n,m)+rzwq*waux(402,n,m+1)+zmo*waux(402,moz,m+1)+f( 1)*(waux(315,n,m)-fnket*waux(315,n,m+1))
    waux(504,n,m)=rzqc*waux(403,n,m)+rzwq*waux(403,n,m+1)+zmo*waux(403,moz,m+1)+f( 2)*(waux(316,n,m)-fnket*waux(316,n,m+1))
    waux(505,n,m)=rzqc*waux(404,n,m)+rzwq*waux(404,n,m+1)+zmo*waux(404,moz,m+1)+f( 3)*(waux(317,n,m)-fnket*waux(317,n,m+1))
    waux(506,n,m)=ryqc*waux(406,n,m)+rywq*waux(406,n,m+1)+ymo*waux(406,moy,m+1)+f( 3)*(waux(320,n,m)-fnket*waux(320,n,m+1))
    waux(507,n,m)=ryqc*waux(407,n,m)+rywq*waux(407,n,m+1)+ymo*waux(407,moy,m+1)+f( 2)*(waux(321,n,m)-fnket*waux(321,n,m+1))
    waux(508,n,m)=ryqc*waux(408,n,m)+rywq*waux(408,n,m+1)+ymo*waux(408,moy,m+1)+f( 1)*(waux(322,n,m)-fnket*waux(322,n,m+1))
    waux(509,n,m)=ryqc*waux(409,n,m)+rywq*waux(409,n,m+1)+ymo*waux(409,moy,m+1)
    waux(510,n,m)=rxqc*waux(419,n,m)+rxwq*waux(419,n,m+1)+xmo*waux(419,mox,m+1)+f( 3)*(waux(341,n,m)-fnket*waux(341,n,m+1))
    waux(511,n,m)=rxqc*waux(420,n,m)+rxwq*waux(420,n,m+1)+xmo*waux(420,mox,m+1)+f( 2)*(waux(342,n,m)-fnket*waux(342,n,m+1))
    waux(512,n,m)=rzqc*waux(410,n,m)+rzwq*waux(410,n,m+1)+zmo*waux(410,moz,m+1)
    waux(513,n,m)=rzqc*waux(411,n,m)+rzwq*waux(411,n,m+1)+zmo*waux(411,moz,m+1)+f( 1)*(waux(323,n,m)-fnket*waux(323,n,m+1))
    waux(514,n,m)=rzqc*waux(412,n,m)+rzwq*waux(412,n,m+1)+zmo*waux(412,moz,m+1)+f( 2)*(waux(324,n,m)-fnket*waux(324,n,m+1))
    waux(515,n,m)=rxqc*waux(424,n,m)+rxwq*waux(424,n,m+1)+xmo*waux(424,mox,m+1)+f( 2)*(waux(346,n,m)-fnket*waux(346,n,m+1))
    waux(516,n,m)=rxqc*waux(425,n,m)+rxwq*waux(425,n,m+1)+xmo*waux(425,mox,m+1)+f( 2)*(waux(347,n,m)-fnket*waux(347,n,m+1))
    waux(517,n,m)=rxqc*waux(426,n,m)+rxwq*waux(426,n,m+1)+xmo*waux(426,mox,m+1)+f( 2)*(waux(348,n,m)-fnket*waux(348,n,m+1))
    waux(518,n,m)=ryqc*waux(417,n,m)+rywq*waux(417,n,m+1)+ymo*waux(417,moy,m+1)+f( 2)*(waux(330,n,m)-fnket*waux(330,n,m+1))
    waux(519,n,m)=ryqc*waux(418,n,m)+rywq*waux(418,n,m+1)+ymo*waux(418,moy,m+1)+f( 1)*(waux(331,n,m)-fnket*waux(331,n,m+1))
    waux(520,n,m)=ryqc*waux(419,n,m)+rywq*waux(419,n,m+1)+ymo*waux(419,moy,m+1)
    waux(521,n,m)=rxqc*waux(430,n,m)+rxwq*waux(430,n,m+1)+xmo*waux(430,mox,m+1)+f( 2)*(waux(352,n,m)-fnket*waux(352,n,m+1))
    waux(522,n,m)=rxqc*waux(431,n,m)+rxwq*waux(431,n,m+1)+xmo*waux(431,mox,m+1)+f( 1)*(waux(353,n,m)-fnket*waux(353,n,m+1))
    waux(523,n,m)=rzqc*waux(420,n,m)+rzwq*waux(420,n,m+1)+zmo*waux(420,moz,m+1)
    waux(524,n,m)=rzqc*waux(421,n,m)+rzwq*waux(421,n,m+1)+zmo*waux(421,moz,m+1)+f( 1)*(waux(332,n,m)-fnket*waux(332,n,m+1))
    waux(525,n,m)=rxqc*waux(434,n,m)+rxwq*waux(434,n,m+1)+xmo*waux(434,mox,m+1)+f( 1)*(waux(356,n,m)-fnket*waux(356,n,m+1))
    waux(526,n,m)=rxqc*waux(435,n,m)+rxwq*waux(435,n,m+1)+xmo*waux(435,mox,m+1)+f( 1)*(waux(357,n,m)-fnket*waux(357,n,m+1))
    waux(527,n,m)=rxqc*waux(436,n,m)+rxwq*waux(436,n,m+1)+xmo*waux(436,mox,m+1)+f( 1)*(waux(358,n,m)-fnket*waux(358,n,m+1))
    waux(528,n,m)=rxqc*waux(437,n,m)+rxwq*waux(437,n,m+1)+xmo*waux(437,mox,m+1)+f( 1)*(waux(359,n,m)-fnket*waux(359,n,m+1))
    waux(529,n,m)=rxqc*waux(438,n,m)+rxwq*waux(438,n,m+1)+xmo*waux(438,mox,m+1)+f( 1)*(waux(360,n,m)-fnket*waux(360,n,m+1))
    waux(530,n,m)=rxqc*waux(439,n,m)+rxwq*waux(439,n,m+1)+xmo*waux(439,mox,m+1)+f( 1)*(waux(361,n,m)-fnket*waux(361,n,m+1))
    waux(531,n,m)=ryqc*waux(429,n,m)+rywq*waux(429,n,m+1)+ymo*waux(429,moy,m+1)+f( 1)*(waux(341,n,m)-fnket*waux(341,n,m+1))
    waux(532,n,m)=ryqc*waux(430,n,m)+rywq*waux(430,n,m+1)+ymo*waux(430,moy,m+1)
    waux(533,n,m)=rxqc*waux(442,n,m)+rxwq*waux(442,n,m+1)+xmo*waux(442,mox,m+1)+f( 1)*(waux(364,n,m)-fnket*waux(364,n,m+1))
    waux(534,n,m)=rxqc*waux(443,n,m)+rxwq*waux(443,n,m+1)+xmo*waux(443,mox,m+1)
    waux(535,n,m)=rxqc*waux(444,n,m)+rxwq*waux(444,n,m+1)+xmo*waux(444,mox,m+1)
    waux(536,n,m)=rxqc*waux(445,n,m)+rxwq*waux(445,n,m+1)+xmo*waux(445,mox,m+1)
    waux(537,n,m)=rxqc*waux(446,n,m)+rxwq*waux(446,n,m+1)+xmo*waux(446,mox,m+1)
    waux(538,n,m)=rxqc*waux(447,n,m)+rxwq*waux(447,n,m+1)+xmo*waux(447,mox,m+1)
    waux(539,n,m)=rxqc*waux(448,n,m)+rxwq*waux(448,n,m+1)+xmo*waux(448,mox,m+1)
    waux(540,n,m)=rxqc*waux(449,n,m)+rxwq*waux(449,n,m+1)+xmo*waux(449,mox,m+1)
    waux(541,n,m)=rxqc*waux(450,n,m)+rxwq*waux(450,n,m+1)+xmo*waux(450,mox,m+1)
    waux(542,n,m)=rxqc*waux(451,n,m)+rxwq*waux(451,n,m+1)+xmo*waux(451,mox,m+1)
    waux(543,n,m)=rxqc*waux(452,n,m)+rxwq*waux(452,n,m+1)+xmo*waux(452,mox,m+1)
    waux(544,n,m)=rxqc*waux(453,n,m)+rxwq*waux(453,n,m+1)+xmo*waux(453,mox,m+1)
    waux(545,n,m)=rxqc*waux(454,n,m)+rxwq*waux(454,n,m+1)+xmo*waux(454,mox,m+1)
    waux(546,n,m)=rxqc*waux(455,n,m)+rxwq*waux(455,n,m+1)+xmo*waux(455,mox,m+1)
    waux(547,n,m)=ryqc*waux(443,n,m)+rywq*waux(443,n,m+1)+ymo*waux(443,moy,m+1)+f(12)*(waux(353,n,m)-fnket*waux(353,n,m+1))
    waux(548,n,m)=rzqc*waux(443,n,m)+rzwq*waux(443,n,m+1)+zmo*waux(443,moz,m+1)
    waux(549,n,m)=rzqc*waux(444,n,m)+rzwq*waux(444,n,m+1)+zmo*waux(444,moz,m+1)+f( 1)*(waux(353,n,m)-fnket*waux(353,n,m+1))
    waux(550,n,m)=rzqc*waux(445,n,m)+rzwq*waux(445,n,m+1)+zmo*waux(445,moz,m+1)+f( 2)*(waux(354,n,m)-fnket*waux(354,n,m+1))
    waux(551,n,m)=rzqc*waux(446,n,m)+rzwq*waux(446,n,m+1)+zmo*waux(446,moz,m+1)+f( 3)*(waux(355,n,m)-fnket*waux(355,n,m+1))
    waux(552,n,m)=rzqc*waux(447,n,m)+rzwq*waux(447,n,m+1)+zmo*waux(447,moz,m+1)+f( 4)*(waux(356,n,m)-fnket*waux(356,n,m+1))
    waux(553,n,m)=rzqc*waux(448,n,m)+rzwq*waux(448,n,m+1)+zmo*waux(448,moz,m+1)+f( 5)*(waux(357,n,m)-fnket*waux(357,n,m+1))
    waux(554,n,m)=ryqc*waux(450,n,m)+rywq*waux(450,n,m+1)+ymo*waux(450,moy,m+1)+f( 5)*(waux(360,n,m)-fnket*waux(360,n,m+1))
    waux(555,n,m)=ryqc*waux(451,n,m)+rywq*waux(451,n,m+1)+ymo*waux(451,moy,m+1)+f( 4)*(waux(361,n,m)-fnket*waux(361,n,m+1))
    waux(556,n,m)=ryqc*waux(452,n,m)+rywq*waux(452,n,m+1)+ymo*waux(452,moy,m+1)+f( 3)*(waux(362,n,m)-fnket*waux(362,n,m+1))
    waux(557,n,m)=ryqc*waux(453,n,m)+rywq*waux(453,n,m+1)+ymo*waux(453,moy,m+1)+f( 2)*(waux(363,n,m)-fnket*waux(363,n,m+1))
    waux(558,n,m)=ryqc*waux(454,n,m)+rywq*waux(454,n,m+1)+ymo*waux(454,moy,m+1)+f( 1)*(waux(364,n,m)-fnket*waux(364,n,m+1))
    waux(559,n,m)=ryqc*waux(455,n,m)+rywq*waux(455,n,m+1)+ymo*waux(455,moy,m+1)
    waux(560,n,m)=rzqc*waux(455,n,m)+rzwq*waux(455,n,m+1)+zmo*waux(455,moz,m+1)+f(12)*(waux(364,n,m)-fnket*waux(364,n,m+1))
    enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 560  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,560,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    !-- [14 | * ] --!
    f(13) = f(12) + f(1) 
    mmax=mmax-1
    do m=1,mmax
    do n=1,nmax
    mox =      mo(1,n)
    moy =      mo(3,n)
    moz =      mo(5,n)
    xmo = gama(mo(2,n))
    ymo = gama(mo(4,n))
    zmo = gama(mo(6,n))
    waux(561,n,m)=rxqc*waux(456,n,m)+rxwq*waux(456,n,m+1)+xmo*waux(456,mox,m+1)+f(13)*(waux(365,n,m)-fnket*waux(365,n,m+1))
    waux(562,n,m)=ryqc*waux(456,n,m)+rywq*waux(456,n,m+1)+ymo*waux(456,moy,m+1)
    waux(563,n,m)=rzqc*waux(456,n,m)+rzwq*waux(456,n,m+1)+zmo*waux(456,moz,m+1)
    waux(564,n,m)=ryqc*waux(457,n,m)+rywq*waux(457,n,m+1)+ymo*waux(457,moy,m+1)+f( 1)*(waux(365,n,m)-fnket*waux(365,n,m+1))
    waux(565,n,m)=rzqc*waux(457,n,m)+rzwq*waux(457,n,m+1)+zmo*waux(457,moz,m+1)
    waux(566,n,m)=rzqc*waux(458,n,m)+rzwq*waux(458,n,m+1)+zmo*waux(458,moz,m+1)+f( 1)*(waux(365,n,m)-fnket*waux(365,n,m+1))
    waux(567,n,m)=ryqc*waux(459,n,m)+rywq*waux(459,n,m+1)+ymo*waux(459,moy,m+1)+f( 2)*(waux(366,n,m)-fnket*waux(366,n,m+1))
    waux(568,n,m)=rzqc*waux(459,n,m)+rzwq*waux(459,n,m+1)+zmo*waux(459,moz,m+1)
    waux(569,n,m)=ryqc*waux(461,n,m)+rywq*waux(461,n,m+1)+ymo*waux(461,moy,m+1)
    waux(570,n,m)=rzqc*waux(461,n,m)+rzwq*waux(461,n,m+1)+zmo*waux(461,moz,m+1)+f( 2)*(waux(367,n,m)-fnket*waux(367,n,m+1))
    waux(571,n,m)=ryqc*waux(462,n,m)+rywq*waux(462,n,m+1)+ymo*waux(462,moy,m+1)+f( 3)*(waux(368,n,m)-fnket*waux(368,n,m+1))
    waux(572,n,m)=rzqc*waux(462,n,m)+rzwq*waux(462,n,m+1)+zmo*waux(462,moz,m+1)
    waux(573,n,m)=rzqc*waux(463,n,m)+rzwq*waux(463,n,m+1)+zmo*waux(463,moz,m+1)+f( 1)*(waux(368,n,m)-fnket*waux(368,n,m+1))
    waux(574,n,m)=ryqc*waux(465,n,m)+rywq*waux(465,n,m+1)+ymo*waux(465,moy,m+1)
    waux(575,n,m)=rzqc*waux(465,n,m)+rzwq*waux(465,n,m+1)+zmo*waux(465,moz,m+1)+f( 3)*(waux(370,n,m)-fnket*waux(370,n,m+1))
    waux(576,n,m)=ryqc*waux(466,n,m)+rywq*waux(466,n,m+1)+ymo*waux(466,moy,m+1)+f( 4)*(waux(371,n,m)-fnket*waux(371,n,m+1))
    waux(577,n,m)=rzqc*waux(466,n,m)+rzwq*waux(466,n,m+1)+zmo*waux(466,moz,m+1)
    waux(578,n,m)=rzqc*waux(467,n,m)+rzwq*waux(467,n,m+1)+zmo*waux(467,moz,m+1)+f( 1)*(waux(371,n,m)-fnket*waux(371,n,m+1))
    waux(579,n,m)=ryqc*waux(469,n,m)+rywq*waux(469,n,m+1)+ymo*waux(469,moy,m+1)+f( 1)*(waux(374,n,m)-fnket*waux(374,n,m+1))
    waux(580,n,m)=ryqc*waux(470,n,m)+rywq*waux(470,n,m+1)+ymo*waux(470,moy,m+1)
    waux(581,n,m)=rzqc*waux(470,n,m)+rzwq*waux(470,n,m+1)+zmo*waux(470,moz,m+1)+f( 4)*(waux(374,n,m)-fnket*waux(374,n,m+1))
    waux(582,n,m)=ryqc*waux(471,n,m)+rywq*waux(471,n,m+1)+ymo*waux(471,moy,m+1)+f( 5)*(waux(375,n,m)-fnket*waux(375,n,m+1))
    waux(583,n,m)=rzqc*waux(471,n,m)+rzwq*waux(471,n,m+1)+zmo*waux(471,moz,m+1)
    waux(584,n,m)=rzqc*waux(472,n,m)+rzwq*waux(472,n,m+1)+zmo*waux(472,moz,m+1)+f( 1)*(waux(375,n,m)-fnket*waux(375,n,m+1))
    waux(585,n,m)=rzqc*waux(473,n,m)+rzwq*waux(473,n,m+1)+zmo*waux(473,moz,m+1)+f( 2)*(waux(376,n,m)-fnket*waux(376,n,m+1))
    waux(586,n,m)=ryqc*waux(475,n,m)+rywq*waux(475,n,m+1)+ymo*waux(475,moy,m+1)+f( 1)*(waux(379,n,m)-fnket*waux(379,n,m+1))
    waux(587,n,m)=ryqc*waux(476,n,m)+rywq*waux(476,n,m+1)+ymo*waux(476,moy,m+1)
    waux(588,n,m)=rzqc*waux(476,n,m)+rzwq*waux(476,n,m+1)+zmo*waux(476,moz,m+1)+f( 5)*(waux(379,n,m)-fnket*waux(379,n,m+1))
    waux(589,n,m)=ryqc*waux(477,n,m)+rywq*waux(477,n,m+1)+ymo*waux(477,moy,m+1)+f( 6)*(waux(380,n,m)-fnket*waux(380,n,m+1))
    waux(590,n,m)=rzqc*waux(477,n,m)+rzwq*waux(477,n,m+1)+zmo*waux(477,moz,m+1)
    waux(591,n,m)=rzqc*waux(478,n,m)+rzwq*waux(478,n,m+1)+zmo*waux(478,moz,m+1)+f( 1)*(waux(380,n,m)-fnket*waux(380,n,m+1))
    waux(592,n,m)=rzqc*waux(479,n,m)+rzwq*waux(479,n,m+1)+zmo*waux(479,moz,m+1)+f( 2)*(waux(381,n,m)-fnket*waux(381,n,m+1))
    waux(593,n,m)=ryqc*waux(481,n,m)+rywq*waux(481,n,m+1)+ymo*waux(481,moy,m+1)+f( 2)*(waux(384,n,m)-fnket*waux(384,n,m+1))
    waux(594,n,m)=ryqc*waux(482,n,m)+rywq*waux(482,n,m+1)+ymo*waux(482,moy,m+1)+f( 1)*(waux(385,n,m)-fnket*waux(385,n,m+1))
    waux(595,n,m)=ryqc*waux(483,n,m)+rywq*waux(483,n,m+1)+ymo*waux(483,moy,m+1)
    waux(596,n,m)=rzqc*waux(483,n,m)+rzwq*waux(483,n,m+1)+zmo*waux(483,moz,m+1)+f( 6)*(waux(385,n,m)-fnket*waux(385,n,m+1))
    waux(597,n,m)=rxqc*waux(492,n,m)+rxwq*waux(492,n,m+1)+xmo*waux(492,mox,m+1)+f( 5)*(waux(401,n,m)-fnket*waux(401,n,m+1))
    waux(598,n,m)=rzqc*waux(484,n,m)+rzwq*waux(484,n,m+1)+zmo*waux(484,moz,m+1)
    waux(599,n,m)=rzqc*waux(485,n,m)+rzwq*waux(485,n,m+1)+zmo*waux(485,moz,m+1)+f( 1)*(waux(386,n,m)-fnket*waux(386,n,m+1))
    waux(600,n,m)=rzqc*waux(486,n,m)+rzwq*waux(486,n,m+1)+zmo*waux(486,moz,m+1)+f( 2)*(waux(387,n,m)-fnket*waux(387,n,m+1))
    waux(601,n,m)=rzqc*waux(487,n,m)+rzwq*waux(487,n,m+1)+zmo*waux(487,moz,m+1)+f( 3)*(waux(388,n,m)-fnket*waux(388,n,m+1))
    waux(602,n,m)=ryqc*waux(489,n,m)+rywq*waux(489,n,m+1)+ymo*waux(489,moy,m+1)+f( 2)*(waux(391,n,m)-fnket*waux(391,n,m+1))
    waux(603,n,m)=ryqc*waux(490,n,m)+rywq*waux(490,n,m+1)+ymo*waux(490,moy,m+1)+f( 1)*(waux(392,n,m)-fnket*waux(392,n,m+1))
    waux(604,n,m)=ryqc*waux(491,n,m)+rywq*waux(491,n,m+1)+ymo*waux(491,moy,m+1)
    waux(605,n,m)=rxqc*waux(500,n,m)+rxwq*waux(500,n,m+1)+xmo*waux(500,mox,m+1)+f( 5)*(waux(409,n,m)-fnket*waux(409,n,m+1))
    waux(606,n,m)=rxqc*waux(501,n,m)+rxwq*waux(501,n,m+1)+xmo*waux(501,mox,m+1)+f( 4)*(waux(410,n,m)-fnket*waux(410,n,m+1))
    waux(607,n,m)=rzqc*waux(492,n,m)+rzwq*waux(492,n,m+1)+zmo*waux(492,moz,m+1)
    waux(608,n,m)=rzqc*waux(493,n,m)+rzwq*waux(493,n,m+1)+zmo*waux(493,moz,m+1)+f( 1)*(waux(393,n,m)-fnket*waux(393,n,m+1))
    waux(609,n,m)=rzqc*waux(494,n,m)+rzwq*waux(494,n,m+1)+zmo*waux(494,moz,m+1)+f( 2)*(waux(394,n,m)-fnket*waux(394,n,m+1))
    waux(610,n,m)=rzqc*waux(495,n,m)+rzwq*waux(495,n,m+1)+zmo*waux(495,moz,m+1)+f( 3)*(waux(395,n,m)-fnket*waux(395,n,m+1))
    waux(611,n,m)=ryqc*waux(497,n,m)+rywq*waux(497,n,m+1)+ymo*waux(497,moy,m+1)+f( 3)*(waux(398,n,m)-fnket*waux(398,n,m+1))
    waux(612,n,m)=ryqc*waux(498,n,m)+rywq*waux(498,n,m+1)+ymo*waux(498,moy,m+1)+f( 2)*(waux(399,n,m)-fnket*waux(399,n,m+1))
    waux(613,n,m)=ryqc*waux(499,n,m)+rywq*waux(499,n,m+1)+ymo*waux(499,moy,m+1)+f( 1)*(waux(400,n,m)-fnket*waux(400,n,m+1))
    waux(614,n,m)=ryqc*waux(500,n,m)+rywq*waux(500,n,m+1)+ymo*waux(500,moy,m+1)
    waux(615,n,m)=rxqc*waux(510,n,m)+rxwq*waux(510,n,m+1)+xmo*waux(510,mox,m+1)+f( 4)*(waux(419,n,m)-fnket*waux(419,n,m+1))
    waux(616,n,m)=rxqc*waux(511,n,m)+rxwq*waux(511,n,m+1)+xmo*waux(511,mox,m+1)+f( 3)*(waux(420,n,m)-fnket*waux(420,n,m+1))
    waux(617,n,m)=rzqc*waux(501,n,m)+rzwq*waux(501,n,m+1)+zmo*waux(501,moz,m+1)
    waux(618,n,m)=rzqc*waux(502,n,m)+rzwq*waux(502,n,m+1)+zmo*waux(502,moz,m+1)+f( 1)*(waux(401,n,m)-fnket*waux(401,n,m+1))
    waux(619,n,m)=rzqc*waux(503,n,m)+rzwq*waux(503,n,m+1)+zmo*waux(503,moz,m+1)+f( 2)*(waux(402,n,m)-fnket*waux(402,n,m+1))
    waux(620,n,m)=rzqc*waux(504,n,m)+rzwq*waux(504,n,m+1)+zmo*waux(504,moz,m+1)+f( 3)*(waux(403,n,m)-fnket*waux(403,n,m+1))
    waux(621,n,m)=rxqc*waux(516,n,m)+rxwq*waux(516,n,m+1)+xmo*waux(516,mox,m+1)+f( 3)*(waux(425,n,m)-fnket*waux(425,n,m+1))
    waux(622,n,m)=ryqc*waux(507,n,m)+rywq*waux(507,n,m+1)+ymo*waux(507,moy,m+1)+f( 3)*(waux(407,n,m)-fnket*waux(407,n,m+1))
    waux(623,n,m)=ryqc*waux(508,n,m)+rywq*waux(508,n,m+1)+ymo*waux(508,moy,m+1)+f( 2)*(waux(408,n,m)-fnket*waux(408,n,m+1))
    waux(624,n,m)=ryqc*waux(509,n,m)+rywq*waux(509,n,m+1)+ymo*waux(509,moy,m+1)+f( 1)*(waux(409,n,m)-fnket*waux(409,n,m+1))
    waux(625,n,m)=ryqc*waux(510,n,m)+rywq*waux(510,n,m+1)+ymo*waux(510,moy,m+1)
    waux(626,n,m)=rxqc*waux(521,n,m)+rxwq*waux(521,n,m+1)+xmo*waux(521,mox,m+1)+f( 3)*(waux(430,n,m)-fnket*waux(430,n,m+1))
    waux(627,n,m)=rxqc*waux(522,n,m)+rxwq*waux(522,n,m+1)+xmo*waux(522,mox,m+1)+f( 2)*(waux(431,n,m)-fnket*waux(431,n,m+1))
    waux(628,n,m)=rzqc*waux(511,n,m)+rzwq*waux(511,n,m+1)+zmo*waux(511,moz,m+1)
    waux(629,n,m)=rzqc*waux(512,n,m)+rzwq*waux(512,n,m+1)+zmo*waux(512,moz,m+1)+f( 1)*(waux(410,n,m)-fnket*waux(410,n,m+1))
    waux(630,n,m)=rzqc*waux(513,n,m)+rzwq*waux(513,n,m+1)+zmo*waux(513,moz,m+1)+f( 2)*(waux(411,n,m)-fnket*waux(411,n,m+1))
    waux(631,n,m)=rxqc*waux(526,n,m)+rxwq*waux(526,n,m+1)+xmo*waux(526,mox,m+1)+f( 2)*(waux(435,n,m)-fnket*waux(435,n,m+1))
    waux(632,n,m)=rxqc*waux(527,n,m)+rxwq*waux(527,n,m+1)+xmo*waux(527,mox,m+1)+f( 2)*(waux(436,n,m)-fnket*waux(436,n,m+1))
    waux(633,n,m)=rxqc*waux(528,n,m)+rxwq*waux(528,n,m+1)+xmo*waux(528,mox,m+1)+f( 2)*(waux(437,n,m)-fnket*waux(437,n,m+1))
    waux(634,n,m)=rxqc*waux(529,n,m)+rxwq*waux(529,n,m+1)+xmo*waux(529,mox,m+1)+f( 2)*(waux(438,n,m)-fnket*waux(438,n,m+1))
    waux(635,n,m)=ryqc*waux(519,n,m)+rywq*waux(519,n,m+1)+ymo*waux(519,moy,m+1)+f( 2)*(waux(418,n,m)-fnket*waux(418,n,m+1))
    waux(636,n,m)=ryqc*waux(520,n,m)+rywq*waux(520,n,m+1)+ymo*waux(520,moy,m+1)+f( 1)*(waux(419,n,m)-fnket*waux(419,n,m+1))
    waux(637,n,m)=ryqc*waux(521,n,m)+rywq*waux(521,n,m+1)+ymo*waux(521,moy,m+1)
    waux(638,n,m)=rxqc*waux(533,n,m)+rxwq*waux(533,n,m+1)+xmo*waux(533,mox,m+1)+f( 2)*(waux(442,n,m)-fnket*waux(442,n,m+1))
    waux(639,n,m)=rxqc*waux(534,n,m)+rxwq*waux(534,n,m+1)+xmo*waux(534,mox,m+1)+f( 1)*(waux(443,n,m)-fnket*waux(443,n,m+1))
    waux(640,n,m)=rzqc*waux(522,n,m)+rzwq*waux(522,n,m+1)+zmo*waux(522,moz,m+1)
    waux(641,n,m)=rzqc*waux(523,n,m)+rzwq*waux(523,n,m+1)+zmo*waux(523,moz,m+1)+f( 1)*(waux(420,n,m)-fnket*waux(420,n,m+1))
    waux(642,n,m)=rxqc*waux(537,n,m)+rxwq*waux(537,n,m+1)+xmo*waux(537,mox,m+1)+f( 1)*(waux(446,n,m)-fnket*waux(446,n,m+1))
    waux(643,n,m)=rxqc*waux(538,n,m)+rxwq*waux(538,n,m+1)+xmo*waux(538,mox,m+1)+f( 1)*(waux(447,n,m)-fnket*waux(447,n,m+1))
    waux(644,n,m)=rxqc*waux(539,n,m)+rxwq*waux(539,n,m+1)+xmo*waux(539,mox,m+1)+f( 1)*(waux(448,n,m)-fnket*waux(448,n,m+1))
    waux(645,n,m)=rxqc*waux(540,n,m)+rxwq*waux(540,n,m+1)+xmo*waux(540,mox,m+1)+f( 1)*(waux(449,n,m)-fnket*waux(449,n,m+1))
    waux(646,n,m)=rxqc*waux(541,n,m)+rxwq*waux(541,n,m+1)+xmo*waux(541,mox,m+1)+f( 1)*(waux(450,n,m)-fnket*waux(450,n,m+1))
    waux(647,n,m)=rxqc*waux(542,n,m)+rxwq*waux(542,n,m+1)+xmo*waux(542,mox,m+1)+f( 1)*(waux(451,n,m)-fnket*waux(451,n,m+1))
    waux(648,n,m)=rxqc*waux(543,n,m)+rxwq*waux(543,n,m+1)+xmo*waux(543,mox,m+1)+f( 1)*(waux(452,n,m)-fnket*waux(452,n,m+1))
    waux(649,n,m)=ryqc*waux(532,n,m)+rywq*waux(532,n,m+1)+ymo*waux(532,moy,m+1)+f( 1)*(waux(430,n,m)-fnket*waux(430,n,m+1))
    waux(650,n,m)=ryqc*waux(533,n,m)+rywq*waux(533,n,m+1)+ymo*waux(533,moy,m+1)
    waux(651,n,m)=rxqc*waux(546,n,m)+rxwq*waux(546,n,m+1)+xmo*waux(546,mox,m+1)+f( 1)*(waux(455,n,m)-fnket*waux(455,n,m+1))
    waux(652,n,m)=rxqc*waux(547,n,m)+rxwq*waux(547,n,m+1)+xmo*waux(547,mox,m+1)
    waux(653,n,m)=rxqc*waux(548,n,m)+rxwq*waux(548,n,m+1)+xmo*waux(548,mox,m+1)
    waux(654,n,m)=rxqc*waux(549,n,m)+rxwq*waux(549,n,m+1)+xmo*waux(549,mox,m+1)
    waux(655,n,m)=rxqc*waux(550,n,m)+rxwq*waux(550,n,m+1)+xmo*waux(550,mox,m+1)
    waux(656,n,m)=rxqc*waux(551,n,m)+rxwq*waux(551,n,m+1)+xmo*waux(551,mox,m+1)
    waux(657,n,m)=rxqc*waux(552,n,m)+rxwq*waux(552,n,m+1)+xmo*waux(552,mox,m+1)
    waux(658,n,m)=rxqc*waux(553,n,m)+rxwq*waux(553,n,m+1)+xmo*waux(553,mox,m+1)
    waux(659,n,m)=rxqc*waux(554,n,m)+rxwq*waux(554,n,m+1)+xmo*waux(554,mox,m+1)
    waux(660,n,m)=rxqc*waux(555,n,m)+rxwq*waux(555,n,m+1)+xmo*waux(555,mox,m+1)
    waux(661,n,m)=rxqc*waux(556,n,m)+rxwq*waux(556,n,m+1)+xmo*waux(556,mox,m+1)
    waux(662,n,m)=rxqc*waux(557,n,m)+rxwq*waux(557,n,m+1)+xmo*waux(557,mox,m+1)
    waux(663,n,m)=rxqc*waux(558,n,m)+rxwq*waux(558,n,m+1)+xmo*waux(558,mox,m+1)
    waux(664,n,m)=rxqc*waux(559,n,m)+rxwq*waux(559,n,m+1)+xmo*waux(559,mox,m+1)
    waux(665,n,m)=rxqc*waux(560,n,m)+rxwq*waux(560,n,m+1)+xmo*waux(560,mox,m+1)
    waux(666,n,m)=ryqc*waux(547,n,m)+rywq*waux(547,n,m+1)+ymo*waux(547,moy,m+1)+f(13)*(waux(443,n,m)-fnket*waux(443,n,m+1))
    waux(667,n,m)=rzqc*waux(547,n,m)+rzwq*waux(547,n,m+1)+zmo*waux(547,moz,m+1)
    waux(668,n,m)=rzqc*waux(548,n,m)+rzwq*waux(548,n,m+1)+zmo*waux(548,moz,m+1)+f( 1)*(waux(443,n,m)-fnket*waux(443,n,m+1))
    waux(669,n,m)=rzqc*waux(549,n,m)+rzwq*waux(549,n,m+1)+zmo*waux(549,moz,m+1)+f( 2)*(waux(444,n,m)-fnket*waux(444,n,m+1))
    waux(670,n,m)=rzqc*waux(550,n,m)+rzwq*waux(550,n,m+1)+zmo*waux(550,moz,m+1)+f( 3)*(waux(445,n,m)-fnket*waux(445,n,m+1))
    waux(671,n,m)=rzqc*waux(551,n,m)+rzwq*waux(551,n,m+1)+zmo*waux(551,moz,m+1)+f( 4)*(waux(446,n,m)-fnket*waux(446,n,m+1))
    waux(672,n,m)=rzqc*waux(552,n,m)+rzwq*waux(552,n,m+1)+zmo*waux(552,moz,m+1)+f( 5)*(waux(447,n,m)-fnket*waux(447,n,m+1))
    waux(673,n,m)=rzqc*waux(553,n,m)+rzwq*waux(553,n,m+1)+zmo*waux(553,moz,m+1)+f( 6)*(waux(448,n,m)-fnket*waux(448,n,m+1))
    waux(674,n,m)=ryqc*waux(555,n,m)+rywq*waux(555,n,m+1)+ymo*waux(555,moy,m+1)+f( 5)*(waux(451,n,m)-fnket*waux(451,n,m+1))
    waux(675,n,m)=ryqc*waux(556,n,m)+rywq*waux(556,n,m+1)+ymo*waux(556,moy,m+1)+f( 4)*(waux(452,n,m)-fnket*waux(452,n,m+1))
    waux(676,n,m)=ryqc*waux(557,n,m)+rywq*waux(557,n,m+1)+ymo*waux(557,moy,m+1)+f( 3)*(waux(453,n,m)-fnket*waux(453,n,m+1))
    waux(677,n,m)=ryqc*waux(558,n,m)+rywq*waux(558,n,m+1)+ymo*waux(558,moy,m+1)+f( 2)*(waux(454,n,m)-fnket*waux(454,n,m+1))
    waux(678,n,m)=ryqc*waux(559,n,m)+rywq*waux(559,n,m+1)+ymo*waux(559,moy,m+1)+f( 1)*(waux(455,n,m)-fnket*waux(455,n,m+1))
    waux(679,n,m)=ryqc*waux(560,n,m)+rywq*waux(560,n,m+1)+ymo*waux(560,moy,m+1)
    waux(680,n,m)=rzqc*waux(560,n,m)+rzwq*waux(560,n,m+1)+zmo*waux(560,moz,m+1)+f(13)*(waux(455,n,m)-fnket*waux(455,n,m+1))
    enddo
    enddo
    if( mmax.eq.1 )then
      mmin = max(1,(lc*(lc-1)*(lc+1))/6+1)
      li2  = 680  - mmin + 1
      li1  = nmax - nmin + 1
      call nuclear_getket(mmin,680,nmin,nmax,waux,fint,max(li2,lf2))
      return
    endif

    write(6,*) ' Improper choice of the angular momentum l'
    write(6,*) ' expected:  <1-15>'
    write(6,*) ' required: ',lcd
    write(6,*) ' Abnormal termination in module: interest_osr, class: nuclear, routine: nuclear_ket'
    stop
  end subroutine

! ---------------------------------------------------------------------------------------
! downward recursion: x < 12.0
! ---------------------------------------------------------------------------------------
  subroutine nuclear_boysxina(nval,xval,xdif,fnx,fnxtab)  
    integer :: nval, k, nmax
    real(8) :: xval, xdif, fnx(*),fnxtab(*)
    real(8) :: expx, xval2, sum, fact, denom 


    expx  = dexp(-xval)
    xval2 = xval + xval
    sum   = fnxtab(nval)
    fact  = 1.0d0

    do k=1,6
      fact = -fact*xdif/dfloat(k)
      sum  =  sum + fnxtab(nval+k)*fact
    enddo
    fnx(nval) = sum

    nmax  = nval - 1
    denom = dfloat(nmax+nmax+1)
    do k=nmax,1,-1
      denom  = denom - 2.0d0
      fnx(k) = (xval2*fnx(k+1) + expx)/denom
    enddo
  end subroutine

! ---------------------------------------------------------------------------------------
! upward recursion: near asymptotic region
! ---------------------------------------------------------------------------------------
  subroutine nuclear_boysxinb(nval,xval,fnx)  
    integer :: nval, n
    real(8) :: xval, fnx(*)
    real(8) :: rxval2, expx, gval, fact, denom
    real(8), parameter :: gfac0 =  0.4999489092d0,&
                          gfac1 = -0.2473631686d0,&
                          gfac2 =  0.3211809090d0,&
                          gfac3 = -0.3811559346d0 

    rxval2 = 1.0d0/(xval + xval)
    expx   = dexp(-xval)
    gval   = gfac0 + (gfac1 + (gfac2 + gfac3/xval)/xval)/xval
    fnx(1) = 0.5d0*dsqrt(pi/xval) - expx*gval/xval

    fact  = expx*rxval2
    denom = -1.0d0
    do n=2,nval,1
      denom  = denom + 2.0d0
      fnx(n) = denom*rxval2*fnx(n-1) - fact
    enddo
  end subroutine

! ---------------------------------------------------------------------------------------
! asymptotic region (54)
! ---------------------------------------------------------------------------------------
  subroutine nuclear_boysxinc(nval,xval,fnx) 
    integer :: nval, n
    real(8) :: rxval2, denom
    real(8) :: xval, fnx(*)

    rxval2 = 1.0d0/(xval + xval)
    fnx(1) = 0.5d0*dsqrt(pi/xval)

    denom = -1.0d0
    do n=2,nval,1
      denom  = denom + 2.0d0
      fnx(n) = denom*rxval2*fnx(n-1)
    enddo
  end subroutine

! ---------------------------------------------------------------------------------------
! accumulate bra results 
! ---------------------------------------------------------------------------------------
  subroutine nuclear_getbra(iini,iend,a,b,li1,li2)
    integer :: i,iini,iend,li1,li2
    real(8) :: a(*),b(*)
 
    li1=0 
    do i=iini,iend
      li1=li1+1
      b(li1)=a(i)
    enddo
    li2=1
  end subroutine

! ---------------------------------------------------------------------------------------
! accumulate ket results 
! ---------------------------------------------------------------------------------------
  subroutine nuclear_getket(jini,jend,iini,iend,a,b,ldb)
    integer :: ldb
    integer :: i,ii,iini,iend
    integer :: j,jj,jini,jend
    real(8) :: a(ncc2,*),b(ldb,*)
 
    ii=0
    do i=iini,iend
      ii=ii+1
      jj=0
      do j=jini,jend
        jj=jj+1 
        b(jj,ii)=a(j,i)
      enddo
    enddo
    if(jj.gt.ldb)then
      write(6,*) ' Error: (li2>lf2)' 
      write(6,*) ' Abnormal termination in module: interest_osr, class: nuclear, routine: getket'
      stop
    endif
  end subroutine

! ---------------------------------------------------------------------------------------
! module initialization 
! ---------------------------------------------------------------------------------------
  subroutine interest_osr_initialize
    integer :: istart, istop, i, j, ii, m1, m2
    integer :: ier, l, lx, ly, lval, mval, nval, madr
    real(8) :: denom, dmax, rdmax, xval, xval2, term, sum, expx
    integer, allocatable :: lin(:), min(:), nin(:)

    allocate( lin(ncc2), min(ncc2), nin(ncc2), stat=ier )
    if( ier.ne.0 )then
      write(6,*) ' Error in allocation of: lin(:), min(:), nin(:)'
      stop       ' Stop in interest_osr_initialize'
    endif

    madr = 0
    do l=0,(lmx2-1),1
      do lx=0,l
        lval = l - lx
        do ly=0,lx
          mval = lx - ly
          nval = ly
          madr      = madr + 1
          lin(madr) = lval
          min(madr) = mval
          nin(madr) = nval
        enddo
      enddo
    enddo

    m1=1
    m2=0
    mo(1,1)=1
    mo(2,1)=lin(1)
    do i=2,lmx2
      ii=(i*(i-1))/2
      do j=1,ii
        m2=m2+1 
        m1=m1+1
        mo(1,m1)=m2
        mo(2,m1)=lin(m1)
      enddo
      do j=1,i
        m1=m1+1
        mo(1,m1)=m2
        mo(2,m1)=lin(m1)
      enddo
    enddo

    m1=1
    m2=0
    mo(3,1)=1
    mo(4,1)=min(1)
    do i=2,lmx2
      m1=m1+1
      mo(3,m1)=(m2+1)
      mo(4,m1)=min(m1)
      do ii=1,(i-1)
        do j=1,ii
          m2=m2+1 
          m1=m1+1
          mo(3,m1)=m2
          mo(4,m1)=min(m1)
        enddo
        m1=m1+1
        mo(3,m1)=m2
        mo(4,m1)=min(m1)
      enddo
    enddo

    m1=1
    m2=0
    mo(5,1)=1
    mo(6,1)=nin(1)
    do i=2,lmx2
      m1=m1+1
      mo(5,m1)=(m2+1)
      mo(6,m1)=nin(m1)
      do ii=1,(i-1)
        m1=m1+1
        mo(5,m1)=(m2+1)
        mo(6,m1)=nin(m1)
        do j=1,ii
          m2=m2+1 
          m1=m1+1
          mo(5,m1)=m2
          mo(6,m1)=nin(m1)
        enddo
      enddo
    enddo

    deallocate( lin )
    deallocate( min )
    deallocate( nin )

!   m2 = 0
!   do i=1,lmx1
!     m1 = m2 + 1      
!     m2 = m2 + i*(i+1)/2 
!     write(6,*)
!     write(6,'(a,i2)') 'Ang. momentum: ',i
!     write(6,*)
!     write(6,'(21(7x,i2))') (ii,ii=m1,m2) 
!     write(6,'(21(7x,i2))') (mo(1,ii),ii=m1,m2) 
!     write(6,'(21(7x,i2))') (mo(3,ii),ii=m1,m2) 
!     write(6,'(21(7x,i2))') (mo(5,ii),ii=m1,m2) 
!     write(6,*)
!     write(6,'(21(7x,i2))') (mo(2,ii),ii=m1,m2) 
!     write(6,'(21(7x,i2))') (mo(4,ii),ii=m1,m2) 
!     write(6,'(21(7x,i2))') (mo(6,ii),ii=m1,m2) 
!   enddo

    istop = lmxx
    denom = 1.0d0 
    do i=1,lmxx  
      fnxtab(i) = 1.0d0/denom 
      denom     = denom + 2.0d0
    enddo 

    dmax  = dfloat(lmxx+lmxx-1)
    rdmax = 1/dmax
    do i = 1,120 
      xval   = 0.1d0*dfloat(i)   
      xval2  = xval + xval
      denom  = dmax
      term   = rdmax
      sum    = term
      do ii  = 2,200
        denom = denom + 2.0d0
        term  = term*xval2/denom
        sum   = sum + term
        if( term.le.1.0d-15 )exit
      enddo 
      istart        = istop + 1
      istop         = istop + lmxx 
      expx          = dexp(-xval)
      fnxtab(istop) = expx*sum
      denom         = dmax 

      do j=istop,istart+1,-1   
        denom       = denom - 2.0d0 
        fnxtab(j-1) = (fnxtab(j)*xval2 + expx)/denom 
      enddo
    enddo
  end subroutine
! ---------------------------------------------------------------------------------------
end module
