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
  subroutine interest_overlap(rkb,                    &
                              fint,na,nb,nrkb,        &
                              la,alpha,ax,ay,az,anorm,&
                              lb,beta, bx,by,bz,bnorm,&
                              ncentr,atomic_data      )
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nrkb
    real(8), intent(out) :: fint(*)
    !-- input --!
    integer, intent(in) :: la,lb
    integer, intent(in) :: ncentr 
    logical, intent(in) :: rkb 
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    real(8), intent(in) :: atomic_data(ncentr,5) 
    !-- local --!
    integer :: n1
    integer :: nint
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: pexp,rexp,rxpa,rypa,rzpa,cntr


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
    else
      nrkb = 2
      nint = 5*na*nb
      call interest_osr_class_overlap(n1,fint,cntr,(la-1),(la+lb+1),pexp,rxpa,rypa,rzpa)
      call interest_hrr_bra_rkb(fint,max(n1,nint),1,la,lb,rxab,ryab,rzab,alpha,beta)
    endif

  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_nuclear_attraction_point_nucleus(rkb,                    &
                                                       fint,na,nb,nrkb,        &
                                                       la,alpha,ax,ay,az,anorm,&
                                                       lb,beta, bx,by,bz,bnorm,&
                                                       ncentr,atomic_data      )
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nrkb
    real(8), intent(out) :: fint(*)
    !-- input --!
    logical, intent(in) :: rkb 
    integer, intent(in) :: la,lb
    integer, intent(in) :: ncentr 
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    real(8), intent(in) :: atomic_data(ncentr,5) 
    !-- local --!
    integer :: nint
    integer :: n1,n2,iat
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: cx,cy,cz,px,py,pz,rxcp,rycp,rzcp
    real(8) :: cntr1,cntr,tval,fnbra,fdbra
    real(8), dimension(ncc2) :: fintl


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
    else
      nrkb = 5
      nint = 5*na*nb
      call interest_hrr_bra_rkb(fint,max(n1,nint),1,la,lb,rxab,ryab,rzab,alpha,beta)
    endif
  end subroutine

! --------------------------------------------------------------------------
!>on output: (5,c,d,5,a,b)
  subroutine interest_eri(rkb,                                            &
                          factor,fint,nint,                               &
                          la,alpha,ax,ay,az,anorm,lb,beta ,bx,by,bz,bnorm,&
                          lc,gamma,cx,cy,cz,cnorm,ld,delta,dx,dy,dz,dnorm )

    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: nint
    real(8), intent(out) :: fint(*)
    !-- input --!
    logical, intent(in) :: rkb 
    real(8), intent(in) :: factor 
    integer, intent(in) :: la,lb,lc,ld
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    real(8), intent(in) :: gamma,cx,cy,cz,cnorm
    real(8), intent(in) :: delta,dx,dy,dz,dnorm
    !-- local --!
    integer :: n1,n2,nab,ncd
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: rxcd,rycd,rzcd,rrcd
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: qexp,sexp,rxqc,ryqc,rzqc
    real(8) :: px,py,pz,qx,qy,qz
    real(8) :: psq,ppq,rho,rxpq,rypq,rzpq
    real(8) :: rxwp,rywp,rzwp,rxwq,rywq,rzwq
    real(8) :: cntr,cntab,cntcd,tval,fnbra,fdbra,fnket,fdket
    real(8), parameter :: pi52 = 2.0d0*dsqrt(pi*pi*pi*pi*pi)



    !-- AB-pair --!
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
      cntab = anorm*bnorm 
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rexp  = alpha*beta/pexp
      cntab = anorm*bnorm*dexp(-rexp*rrab)
    endif


    !-- CD-pair --!
    rxcd = cx - dx
    rycd = cy - dy
    rzcd = cz - dz
    rrcd = rxcd*rxcd + rycd*rycd + rzcd*rzcd

    qexp = gamma + delta
    if( rrcd.lt.1.d-12 )then
      qx    = cx 
      qy    = cy
      qz    = cz
      rxqc  = 0.0d0 
      ryqc  = 0.0d0
      rzqc  = 0.0d0
      rxcd  = 0.0d0  
      rycd  = 0.0d0  
      rzcd  = 0.0d0  
      cntcd = cnorm*dnorm 
    else
      qx    = ( gamma*cx + delta*dx )/qexp
      qy    = ( gamma*cy + delta*dy )/qexp
      qz    = ( gamma*cz + delta*dz )/qexp
      rxqc  = qx - cx
      ryqc  = qy - cy
      rzqc  = qz - cz
      sexp  = gamma*delta/qexp
      cntcd = cnorm*dnorm*dexp(-sexp*rrcd)
    endif


    !-- ABCD --!
    psq = pexp+qexp
    ppq = pexp*qexp
    rho = ppq/psq

    rxpq = px - qx
    rypq = py - qy
    rzpq = pz - qz
    tval = rho*((rxpq*rxpq)+(rypq*rypq)+(rzpq*rzpq))
    cntr = pi52*(1.0d0/dsqrt(psq))*(1.0d0/ppq)*cntab*cntcd*factor

    fnbra = qexp/psq
    fnket = pexp/psq 
    fdbra = 1.0d0/(2.0d0*pexp)
    fdket = 1.0d0/(2.0d0*qexp)
    rxwp  = ( pexp*px + qexp*qx )/psq - px
    rywp  = ( pexp*py + qexp*qy )/psq - py
    rzwp  = ( pexp*pz + qexp*qz )/psq - pz
    rxwq  = ( pexp*px + qexp*qx )/psq - qx
    rywq  = ( pexp*py + qexp*qy )/psq - qy
    rzwq  = ( pexp*pz + qexp*qz )/psq - qz


    if( .not.rkb )then
      nab  = ncc(la)*ncc(lb)
      ncd  = ncc(lc)*ncc(ld)
      call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                         &
                                     (la  ),(la+lb-1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                     (lc  ),(lc+ld-1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
      call interest_hrr_ket(fint,ncd,nab,la,lb,rxab,ryab,rzab,n2)
      call interest_hrr_bra(fint,ncd,nab,lc,ld,rxcd,rycd,rzcd) 
    else
      nab  = ncc(la)*ncc(lb)*5
      ncd  = ncc(lc)*ncc(ld)*5
      call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                         &
                                     (la-1),(la+lb+1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                     (lc-1),(lc+ld+1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
      call interest_hrr_ket_rkb(fint,max(n2,ncd),max(n1,nab),la,lb,rxab,ryab,rzab,alpha,beta,n2)
      call interest_hrr_bra_rkb(fint,max(n2,ncd),       nab ,lc,ld,rxcd,rycd,rzcd,gamma,delta) 
    endif
    nint = nab*ncd 

  end subroutine

! --------------------------------------------------------------------------
!>on output: (5,c,d,5,a,b)
  subroutine interest_eri_test(rkb,factor,fint_inp,nint,                       &
                               la_inp,alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp,&
                               lb_inp,beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp,&
                               lc_inp,gamma_inp,cx_inp,cy_inp,cz_inp,cnorm_inp,&
                               ld_inp,delta_inp,dx_inp,dy_inp,dz_inp,dnorm_inp )

    use module_interest_osr
    use module_interest_hrr
    implicit none

    !> output
    integer, intent(out) :: nint
    real(8), intent(out) :: fint_inp(*)
    !> input
    logical, intent(in) :: rkb 
    real(8), intent(in) :: factor 
    integer, intent(in) :: la_inp,lb_inp,lc_inp,ld_inp
    real(8), intent(in) :: alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp
    real(8), intent(in) :: beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp
    real(8), intent(in) :: gamma_inp,cx_inp,cy_inp,cz_inp,cnorm_inp
    real(8), intent(in) :: delta_inp,dx_inp,dy_inp,dz_inp,dnorm_inp
    !> local
    integer :: index
    integer :: la,lb,lc,ld
    real(8) :: fint(441*441*25)      !fixme: this value should reconsidered!!!!!!
    real(8) :: alpha,ax,ay,az,anorm
    real(8) :: beta, bx,by,bz,bnorm
    real(8) :: gamma,cx,cy,cz,cnorm
    real(8) :: delta,dx,dy,dz,dnorm
    integer :: n1,n2,nab,ncd
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: rxcd,rycd,rzcd,rrcd
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: qexp,sexp,rxqc,ryqc,rzqc
    real(8) :: px,py,pz,qx,qy,qz
    real(8) :: psq,ppq,rho,rxpq,rypq,rzpq
    real(8) :: rxwp,rywp,rzwp,rxwq,rywq,rzwq
    real(8) :: cntr,cntab,cntcd,tval,fnbra,fdbra,fnket,fdket
    real(8), parameter :: pi52 = 2.0d0*dsqrt(pi*pi*pi*pi*pi)


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
      index = 1200
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
      index = 2100
      alpha = beta_inp
      anorm = bnorm_inp
      beta  = alpha_inp  
      bnorm = anorm_inp
    endif

    !> CD-pair indexation
    if( lc_inp >= ld_inp )then
      lc = lc_inp    
      ld = ld_inp   
      cx = cx_inp   
      cy = cy_inp   
      cz = cz_inp    
      dx = dx_inp   
      dy = dy_inp   
      dz = dz_inp   
      gamma = gamma_inp
      cnorm = cnorm_inp
      delta = delta_inp  
      dnorm = dnorm_inp
      index = index +34 
    else
      lc = ld_inp    
      ld = lc_inp   
      cx = dx_inp   
      cy = dy_inp   
      cz = dz_inp    
      dx = cx_inp   
      dy = cy_inp   
      dz = cz_inp   
      gamma = delta_inp
      cnorm = dnorm_inp
      delta = gamma_inp  
      dnorm = cnorm_inp
      index = index +43 
    endif
      
    !> AB-pair
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
      cntab = anorm*bnorm 
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rexp  = alpha*beta/pexp
      cntab = anorm*bnorm*dexp(-rexp*rrab)
    endif

    !> CD-pair
    rxcd = cx - dx
    rycd = cy - dy
    rzcd = cz - dz
    rrcd = rxcd*rxcd + rycd*rycd + rzcd*rzcd

    qexp = gamma + delta
    if( rrcd.lt.1.d-12 )then
      qx    = cx 
      qy    = cy
      qz    = cz
      rxqc  = 0.0d0 
      ryqc  = 0.0d0
      rzqc  = 0.0d0
      rxcd  = 0.0d0  
      rycd  = 0.0d0  
      rzcd  = 0.0d0  
      cntcd = cnorm*dnorm 
    else
      qx    = ( gamma*cx + delta*dx )/qexp
      qy    = ( gamma*cy + delta*dy )/qexp
      qz    = ( gamma*cz + delta*dz )/qexp
      rxqc  = qx - cx
      ryqc  = qy - cy
      rzqc  = qz - cz
      sexp  = gamma*delta/qexp
      cntcd = cnorm*dnorm*dexp(-sexp*rrcd)
    endif


    !> ABCD
    psq = pexp+qexp
    ppq = pexp*qexp
    rho = ppq/psq

    rxpq = px - qx
    rypq = py - qy
    rzpq = pz - qz
    tval = rho*((rxpq*rxpq)+(rypq*rypq)+(rzpq*rzpq))
    cntr = pi52*(1.0d0/dsqrt(psq))*(1.0d0/ppq)*cntab*cntcd*factor

    fnbra = qexp/psq
    fnket = pexp/psq 
    fdbra = 1.0d0/(2.0d0*pexp)
    fdket = 1.0d0/(2.0d0*qexp)
    rxwp  = ( pexp*px + qexp*qx )/psq - px
    rywp  = ( pexp*py + qexp*qy )/psq - py
    rzwp  = ( pexp*pz + qexp*qz )/psq - pz
    rxwq  = ( pexp*px + qexp*qx )/psq - qx
    rywq  = ( pexp*py + qexp*qy )/psq - qy
    rzwq  = ( pexp*pz + qexp*qz )/psq - qz

    if( .not.rkb )then
      nab  = ncc(la)*ncc(lb)
      ncd  = ncc(lc)*ncc(ld)
      nint = nab*ncd 
      call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                         &
                                     (la  ),(la+lb-1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                     (lc  ),(lc+ld-1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
      call interest_hrr_ket(fint,ncd,nab,la,lb,rxab,ryab,rzab,n2)
      call interest_hrr_bra(fint,ncd,nab,lc,ld,rxcd,rycd,rzcd) 
      select case( index )
        case( 1234 ); fint_inp(1:nint)=fint(1:nint) 
        case( 2134 ); call transpose_ab( fint, fint_inp, ncc(la), ncc(lb), ncd ) 
        case( 1243 ); call transpose_cd( fint, fint_inp, ncc(lc), ncc(ld), nab ) 
        case( 2143 ); fint_inp(1:nint)=fint(1:nint)
                      call transpose_ab( fint_inp, fint,     ncc(la), ncc(lb), ncd )
                      call transpose_cd( fint,     fint_inp, ncc(lc), ncc(ld), nab )
      end select 
    else
      nab  = ncc(la)*ncc(lb)*5
      ncd  = ncc(lc)*ncc(ld)*5
      nint = nab*ncd 
      call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                         &
                                     (la-1),(la+lb+1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                     (lc-1),(lc+ld+1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
      call interest_hrr_ket_rkb(fint,max(n2,ncd),max(n1,nab),la,lb,rxab,ryab,rzab,alpha,beta,n2)
      call interest_hrr_bra_rkb(fint,max(n2,ncd),       nab ,lc,ld,rxcd,rycd,rzcd,gamma,delta) 
      select case( index )
        case( 1234 ); fint_inp(1:nint)=fint(1:nint) 
        case( 2134 ); call transpose_ab_rkb( fint, fint_inp, ncc(la), ncc(lb), ncd ) 
        case( 1243 ); call transpose_cd_rkb( fint, fint_inp, ncc(lc), ncc(ld), nab ) 
        case( 2143 ); fint_inp(1:nint)=fint(1:nint)
                      call transpose_ab_rkb( fint_inp, fint,     ncc(la), ncc(lb), ncd )
                      call transpose_cd_rkb( fint,     fint_inp, ncc(lc), ncc(ld), nab )
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
  subroutine transpose_ab( A, B, na, nb, ncd ) 
    implicit none

    !> input
    integer, intent(in) :: na, nb, ncd 
    real(8), intent(in) :: A(ncd,na,nb)
    !> output
    real(8), intent(out) :: B(ncd,nb,na)
    !> local
    integer :: i,j,k 


    do j=1,nb
      do i=1,na
        do k=1,ncd
          B(k,j,i) = A(k,i,j) 
        enddo
      enddo
    enddo
  
  end subroutine

! -------------------------------------------------------------------------
  subroutine transpose_cd( A, B, nc, nd, nab ) 
    implicit none

    !> input
    integer, intent(in) :: nc, nd, nab 
    real(8), intent(in) :: A(nc,nd,nab)
    !> output
    real(8), intent(out) :: B(nd,nc,nab)
    !> local
    integer :: i,j,k 


    do k=1,nab
      do j=1,nd
        do i=1,nc
          B(j,i,k) = A(i,j,k) 
        enddo
      enddo
    enddo
  
  end subroutine

! -------------------------------------------------------------------------
  subroutine transpose_ab_rkb( A, B, na, nb, ncd ) 
    implicit none

    !> input
    integer, intent(in) :: na, nb, ncd 
    real(8), intent(in) :: A(ncd,5,na,nb)
    !> output
    real(8), intent(out) :: B(ncd,5,nb,na)
    !> local
    integer :: i,j,k 


    do j=1,nb
      do i=1,na
        do k=1,ncd
          B(k,1,j,i) =  A(k,1,i,j) 
          B(k,2,j,i) =  A(k,2,i,j) 
          B(k,3,j,i) = -A(k,3,i,j) 
          B(k,4,j,i) = -A(k,4,i,j) 
          B(k,5,j,i) = -A(k,5,i,j) 
        enddo
      enddo
    enddo
  
  end subroutine

! -------------------------------------------------------------------------
  subroutine transpose_cd_rkb( A, B, nc, nd, nab ) 
    implicit none

    !> input
    integer, intent(in) :: nc, nd, nab 
    real(8), intent(in) :: A(5,nc,nd,nab)
    !> output
    real(8), intent(out) :: B(5,nd,nc,nab)
    !> local
    integer :: i,j,k 


    do k=1,nab
      do j=1,nd
        do i=1,nc
          B(1,j,i,k) =  A(1,i,j,k) 
          B(2,j,i,k) =  A(2,i,j,k) 
          B(3,j,i,k) = -A(3,i,j,k) 
          B(4,j,i,k) = -A(4,i,j,k) 
          B(5,j,i,k) = -A(5,i,j,k) 
        enddo
      enddo
    enddo
  
  end subroutine
! -------------------------------------------------------------------------
