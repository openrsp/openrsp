! -------------------------------------------------------------------------
!
! Library:      InteRest
!
! File:         module_interest.f90 
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
  subroutine interest_initialize()

    use module_interest_osr
    implicit none

    if( is_interest_initialized )then
      write(6,*)
      write(6,'(2x,a)')'Integral module InteRest was already initialized'
      write(6,*)
    else
      write(6,*)
      write(6,'(2x,a)')'Integral module InteRest is initialized'
      call print_interest_git_revision_info
      write(6,*)
      call interest_osr_initialize()
      is_interest_initialized = .true.
    endif

  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_overlap(rkb,                                            &
                              fint_inp,na,nb,nrkb,                            &
                              la_inp,alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp,&
                              lb_inp,beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp,&
                              ncentr,atomic_data                              )
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nrkb
    real(8), intent(out) :: fint_inp(*)
    !-- input --!
    logical, intent(in) :: rkb
    integer, intent(in) :: ncentr 
    integer, intent(in) :: la_inp,lb_inp
    real(8), intent(in) :: atomic_data(ncentr,5)
    real(8), intent(in) :: alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp
    real(8), intent(in) :: beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp
    !-- local --!
    integer :: n1
    integer :: nint
    integer :: index
    integer :: la,lb
    real(8) :: fint(441*25)      !fixme
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: alpha,ax,ay,az,anorm
    real(8) :: beta, bx,by,bz,bnorm
    real(8) :: pexp,rexp,rxpa,rypa,rzpa,cntr


    !> AB-pair indexation
    if( la_inp >= lb_inp )then
      la = la_inp    
      lb = lb_inp   
      ax = ax_inp   
      ay = ay_inp   
      az = az_inp    
      bx = bx_inp   
      by = by_inp   
      bz = bz_inp   
      index = 12
      alpha = alpha_inp
      anorm = anorm_inp
      beta  = beta_inp  
      bnorm = bnorm_inp
    else
      la = lb_inp    
      lb = la_inp   
      ax = bx_inp   
      ay = by_inp   
      az = bz_inp    
      bx = ax_inp   
      by = ay_inp   
      bz = az_inp   
      index = 21
      alpha = beta_inp
      anorm = bnorm_inp
      beta  = alpha_inp  
      bnorm = anorm_inp
    endif


    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    rrab = rxab*rxab + ryab*ryab + rzab*rzab

    pexp = alpha + beta
    if( rrab.lt.1.d-12 )then
      rxpa = 0.0d0 
      rypa = 0.0d0
      rzpa = 0.0d0
      rxab = 0.0d0  
      ryab = 0.0d0  
      rzab = 0.0d0  
      cntr = anorm*bnorm 
    else
      rxpa = ( alpha*ax + beta*bx )/pexp - ax
      rypa = ( alpha*ay + beta*by )/pexp - ay
      rzpa = ( alpha*az + beta*bz )/pexp - az
      rexp = alpha*beta/pexp
      cntr = anorm*bnorm*dexp(-rexp*rrab)
    endif

    na = ncc(la)
    nb = ncc(lb)
    if( .not. rkb )then
      nrkb = 1
      nint = na*nb
      call interest_osr_class_overlap(n1,fint,cntr,(la  ),(la+lb-1),pexp,rxpa,rypa,rzpa)
      call interest_hrr_bra(fint,nint,1,la,lb,rxab,ryab,rzab)
      select case( index )
        case( 12 )
                  fint_inp(1:nint)=fint(1:nint) 
        case( 21 )
                  na = ncc(lb)
                  nb = ncc(la)
                  call transpose_ab( fint, fint_inp, ncc(la), ncc(lb) ) 
      end select 
    else
      nrkb = 2
      nint = 5*na*nb
      call interest_osr_class_overlap(n1,fint,cntr,(la-1),(la+lb+1),pexp,rxpa,rypa,rzpa)
      call interest_hrr_bra_rkb(fint,max(n1,nint),1,la,lb,rxab,ryab,rzab,alpha,beta)
      select case( index )
        case( 12 )
                  fint_inp(1:nint)=fint(1:nint) 
        case( 21 )
                  na = ncc(lb)
                  nb = ncc(la)
                  call transpose_ab_rkb( fint, fint_inp, ncc(la), ncc(lb) ) 
      end select 
    endif

  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_nuclear_attraction_point_nucleus(rkb,                                            &
                                                       fint_inp,na,nb,nrkb,                            &
                                                       la_inp,alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp,&
                                                       lb_inp,beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp,&
                                                       ncentr,atomic_data                              )
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nrkb
    real(8), intent(out) :: fint_inp(*)
    !-- input --!
    logical, intent(in) :: rkb
    integer, intent(in) :: ncentr 
    integer, intent(in) :: la_inp,lb_inp
    real(8), intent(in) :: atomic_data(ncentr,5)
    real(8), intent(in) :: alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp
    real(8), intent(in) :: beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp
    !-- local --!
    integer :: nint
    integer :: index
    integer :: la,lb
    integer :: n1,n2,iat
    real(8) :: fint(441*25)      !fixme
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: alpha,ax,ay,az,anorm
    real(8) :: beta, bx,by,bz,bnorm
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: cntr1,cntr,tval,fnbra,fdbra
    real(8) :: cx,cy,cz,px,py,pz,rxcp,rycp,rzcp
    real(8), dimension(ncc2) :: fintl


    !> AB-pair indexation
    if( la_inp >= lb_inp )then
      la = la_inp    
      lb = lb_inp   
      ax = ax_inp   
      ay = ay_inp   
      az = az_inp    
      bx = bx_inp   
      by = by_inp   
      bz = bz_inp   
      index = 12
      alpha = alpha_inp
      anorm = anorm_inp
      beta  = beta_inp  
      bnorm = bnorm_inp
    else
      la = lb_inp    
      lb = la_inp   
      ax = bx_inp   
      ay = by_inp   
      az = bz_inp    
      bx = ax_inp   
      by = ay_inp   
      bz = az_inp   
      index = 21
      alpha = beta_inp
      anorm = bnorm_inp
      beta  = alpha_inp  
      bnorm = anorm_inp
    endif


    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    rrab = rxab*rxab + ryab*ryab + rzab*rzab

    pexp = alpha + beta
    if( rrab.lt.1.d-12 )then
      px    = ax 
      py    = ay
      pz    = az
      rxpa  = 0.0d0 
      rypa  = 0.0d0
      rzpa  = 0.0d0
      rxab  = 0.0d0  
      ryab  = 0.0d0  
      rzab  = 0.0d0  
      cntr1 = 2.0d0*(pi/pexp)*anorm*bnorm 
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rexp  = alpha*beta/pexp
      cntr1 = 2.0d0*(pi/pexp)*anorm*bnorm*dexp(-rexp*rrab)
    endif

    na    = ncc(la)
    nb    = ncc(lb)
    fnbra = 1.0d0 
    fdbra = 1.0d0/(2.0d0*pexp)  

    do iat=1,ncentr
      if( atomic_data(iat,1) <= 0.0d0 )cycle
      cx= atomic_data(iat,2)
      cy= atomic_data(iat,3)
      cz= atomic_data(iat,4)
      rxcp = cx - px
      rycp = cy - py
      rzcp = cz - pz
      cntr = atomic_data(iat,1)*cntr1 
      tval = pexp*(rxcp*rxcp + rycp*rycp + rzcp*rzcp)
      if( .not.rkb )then
        call interest_osr_class_nuclear(n1,n2,fintl,cntr,tval,(la  ),(la+lb-1),       &
                                        1,fnbra,fdbra,rxpa,rypa,rzpa,rxcp,rycp,rzcp,  &
                                        1,1,1,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      else
        call interest_osr_class_nuclear(n1,n2,fintl,cntr,tval,(la-1),(la+lb+1),       &
                                        1,fnbra,fdbra,rxpa,rypa,rzpa,rxcp,rycp,rzcp,  &
                                        1,1,1,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      if( iat==1 ) fint(1:n1) = 0.0d0
      fint(1:n1) = fint(1:n1) - fintl(1:n1)                 
    enddo  

    if( .not.rkb )then
      nrkb = 1
      nint = na*nb
      call interest_hrr_bra(fint,nint,1,la,lb,rxab,ryab,rzab)
      select case( index )
        case( 12 )
                  fint_inp(1:nint)=fint(1:nint) 
        case( 21 )
                  na = ncc(lb)
                  nb = ncc(la)
                  call transpose_ab( fint, fint_inp, ncc(la), ncc(lb) ) 
      end select 
    else
      nrkb = 5
      nint = 5*na*nb
      call interest_hrr_bra_rkb(fint,max(n1,nint),1,la,lb,rxab,ryab,rzab,alpha,beta)
      select case( index )
        case( 12 )
                  fint_inp(1:nint)=fint(1:nint) 
        case( 21 )
                  na = ncc(lb)
                  nb = ncc(la)
                  call transpose_ab_rkb( fint, fint_inp, ncc(la), ncc(lb) ) 
      end select 
    endif

  end subroutine

! --------------------------------------------------------------------------
!>
  subroutine interest_dipole(rkb,                    &
                             fintx,finty,fintz,nint, &
                             la,alpha,ax,ay,az,anorm,&
                             lb,beta, bx,by,bz,bnorm,&
                             gx,gy,gz                )
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: nint
    real(8), intent(out) :: fintx(*)
    real(8), intent(out) :: finty(*)
    real(8), intent(out) :: fintz(*)
    !-- input --!
    logical, intent(in) :: rkb 
    integer, intent(in) :: la,lb
    real(8), intent(in) :: gx,gy,gz 
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    !-- local --!
    integer :: n1
    real(8) :: rxpg,rypg,rzpg
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: pexp,rexp,rxpa,rypa,rzpa,cntr


    if( .not.rkb )then
      nint = ncc(la)*ncc(lb)
    else
      nint = ncc(la)*ncc(lb)*5
    endif
    fintx(1:nint) = 0.0d0
    finty(1:nint) = 0.0d0
    fintz(1:nint) = 0.0d0

    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    pexp = alpha + beta
    rrab = rxab*rxab + ryab*ryab + rzab*rzab


    if( rrab.lt.1.d-12 )then
      rxpa = 0.0d0 
      rypa = 0.0d0
      rzpa = 0.0d0
      rxab = 0.0d0  
      ryab = 0.0d0  
      rzab = 0.0d0  
      cntr = anorm*bnorm 
    else
      rxpa = ( alpha*ax + beta*bx )/pexp - ax
      rypa = ( alpha*ay + beta*by )/pexp - ay
      rzpa = ( alpha*az + beta*bz )/pexp - az
      rexp = alpha*beta/pexp
      cntr = anorm*bnorm*dexp(-rexp*rrab)
    endif

    rxpg = rxpa + ax - gx 
    rypg = rypa + ay - gy 
    rzpg = rzpa + az - gz 

    if( .not.rkb )then
      call interest_osr_class_dipole(n1,fintx,finty,fintz,cntr,(la  ),(la+lb-1),pexp,rxpa,rypa,rzpa,rxpg,rypg,rzpg)
      call interest_hrr_bra(fintx,nint,1,la,lb,rxab,ryab,rzab)
      call interest_hrr_bra(finty,nint,1,la,lb,rxab,ryab,rzab)
      call interest_hrr_bra(fintz,nint,1,la,lb,rxab,ryab,rzab)
    else
      stop 'fix rkb case => stop'
      call interest_osr_class_dipole(n1,fintx,finty,fintz,cntr,(la-1),(la+lb+1),pexp,rxpa,rypa,rzpa,rxpg,rypg,rzpg)
      call interest_hrr_bra_rkb(fintx,nint,1,la,lb,rxab,ryab,rzab,alpha,beta)
      call interest_hrr_bra_rkb(finty,nint,1,la,lb,rxab,ryab,rzab,alpha,beta)
      call interest_hrr_bra_rkb(fintz,nint,1,la,lb,rxab,ryab,rzab,alpha,beta)
    endif

  end subroutine

! -------------------------------------------------------------------------
!>
  subroutine transpose_ab( A, B, na, nb ) 
    !> input
    integer, intent(in) :: na, nb
    real(8), intent(in) :: A(na,nb)
    !> output
    real(8), intent(out) :: B(nb,na)
    !> local
    integer :: i,j
    !> 
    do j=1,nb
      do i=1,na
        B(j,i) = A(i,j) 
      enddo
    enddo
  end subroutine

! -------------------------------------------------------------------------
!>

  subroutine transpose_ab_rkb( A, B, na, nb ) 
    !> input
    integer, intent(in) :: na, nb
    real(8), intent(in) :: A(5,na,nb)
    !> output
    real(8), intent(out) :: B(5,nb,na)
    !> local
    integer :: i,j
    !> 
    do j=1,nb
      do i=1,na
        B(1,j,i) =  A(1,i,j) 
        B(2,j,i) =  A(2,i,j) 
        B(3,j,i) = -A(3,i,j) 
        B(4,j,i) = -A(4,i,j) 
        B(5,j,i) = -A(5,i,j) 
      enddo
    enddo
  end subroutine

! -------------------------------------------------------------------------
