! -------------------------------------------------------------------------
!
! Library:      InteRest
!
! File:         module_interest_eri.f90 
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
!>
  subroutine interest_eri_basic(factor,fint_inp,nint,                          &
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
    real(8), intent(in) :: factor 
    integer, intent(in) :: la_inp,lb_inp,lc_inp,ld_inp
    real(8), intent(in) :: alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp
    real(8), intent(in) :: beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp
    real(8), intent(in) :: gamma_inp,cx_inp,cy_inp,cz_inp,cnorm_inp
    real(8), intent(in) :: delta_inp,dx_inp,dy_inp,dz_inp,dnorm_inp
    !> local
    integer :: index
    integer :: la,lb,lc,ld
    real(8) :: fint(441*441)  !todo
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
    real(8), parameter :: pi52 = 0.34986836655249725d2 


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

    nab  = ncc(la)*ncc(lb)
    ncd  = ncc(lc)*ncc(ld)
    nint = nab*ncd 
    call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                      &
                                    la,(la+lb-1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                    lc,(lc+ld-1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
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


  contains
! --------------------------------------------------------------------------
!>
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
!>
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
end subroutine
