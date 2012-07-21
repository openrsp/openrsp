! ------------------------------------------------------------------------------------
!
! Program:      Dirac
!
! Module:       interest_eri_diff
!
! Description:  - interfacing basic data between Dirac & InteRest program
!               - providing four-center integrals and their analytical derivatives
!
! Contains:
!
! Licensing:    dirac_copyright_start
!                     Copyright (c) 2010 by the authors of DIRAC.
!                     All Rights Reserved.
!
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!               dirac_copyright_end
!
! Author:       Michal  Repisky (michal.repisky@uit.no)
!
! Logs:
!
! Comments:
!
! ------------------------------------------------------------------------------------
MODULE MODULE_INTEREST_ERI_DIFF

  implicit none

  public initialize_interest_eri_diff
  public interest_eri_diff

  private

  type type_atom
       real(8) :: charge
       real(8) :: coordinate_x
       real(8) :: coordinate_y
       real(8) :: coordinate_z
       real(8) :: gnu_exponent
  end  type

  type type_gto
       integer :: index
       integer :: offset
       integer :: lvalue
       integer :: origin
       integer :: sdegen
       integer :: cdegen
       real(8) :: exponent(1)
       real(8) :: coefficient(1)
  end  type

  integer, save :: shell_end(2) = 0
  integer, save :: shell_start(2) = 0

! type(type_constant)             :: constant
  type(type_atom),    allocatable :: atom(:)
  type(type_gto),     allocatable :: gto(:)

  logical :: is_initialized = .false.

CONTAINS
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE initialize_interest_eri_diff()

#ifdef PRG_DIRAC
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbdhf.h"
#include "nucdata.h"

    !> local variables
    integer :: i
    integer :: j
    integer :: ij
    integer :: ier
    type( type_gto ) :: tmpgto

    if (is_initialized) return

    !> initialize InteRest integral package
    call interest_initialize()

    !> test if the large component basis is uncontracted
    if( nplrg.ne.nlarge )then
      write(6,*) ' Error: contracted basis set found => stop!'
      write(6,*) ' Note:  uncontracted basis sets are only supported'
      stop
    endif

    !> interface molecular data
    allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
    do i=1,size(atom)
      atom(i) = type_atom( charge(i), &
                           cord(1,i), &
                           cord(2,i), &
                           cord(3,i), &
                           gnuexp(i)  )
    enddo

    !> interface basis set data
    !> fixme: contracted basis sets
    allocate( gto(nlrgsh+nsmlsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
    i =0
    ij=0
    do while( i < size(gto) )
      do j=1,nrco(i+1)
        i=i+1
        gto(i) = type_gto( ij,                      &
                           0,                       &
                           nhkt(i),                 &
                           ncent(i),                &
                           (2*nhkt(i)-1),           &
                           (nhkt(i)*(nhkt(i)+1))/2, &
                           priexp(i),               &
                           priccf(i,j)              )
        !> note: cartesian indexation
        ij=ij+(nhkt(i)*(nhkt(i)+1))/2
      enddo
    enddo

    !> reorder basis functions to the angular momentum
    !> ascending order to increase the integral performance
    !> set the AO-vector offsets and shell-degeneracy
!   do i=1,size(gto)-1
!     do j=(i+1),size(gto)
!       if( gto(j)%lvalue >= gto(i)%lvalue )cycle
!       tmpgto = gto(i)
!       gto(i) = gto(j)
!       gto(j) = tmpgto
!     enddo
!   enddo

    do i=1,size(gto)-1
      gto(i+1)%offset = gto(i)%offset + gto(i)%cdegen
    enddo

    !> temporary print ...
    do i=1,size(gto)
      write(6,'(1x,5(a,i3,3x),3(a,d15.5,3x))') ' i = ',i,                             &
                                               ' index  = ',     gto(i)%index,        &
                                               ' offset = ',     gto(i)%offset,       &
                                               ' lvalue = ',     gto(i)%lvalue,       &
                                               ' origin = ',     gto(i)%origin,       &
                                               ' exponent = ',   gto(i)%exponent(1),  &
                                               ' coefficient = ',gto(i)%coefficient(1)
    enddo

    shell_start(1) = 1
    shell_end(1)   = nlrgsh
    shell_start(2) = shell_end(1) + 1
    shell_end(2)   = shell_end(1) + nsmlsh

    write(6,*)
    write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
    write(6,'(2x,a,i5  )') 'Total number of spherical basis functions:',sum(gto(:)%sdegen)
    write(6,'(2x,a,i5,a)') 'Total number of cartesian basis functions:',sum(gto(:)%cdegen),' (used for calculation)'
    write(6,*)

    is_initialized = .true.
#endif /* ifdef PRG_DIRAC */

  END SUBROUTINE

! ------------------------------------------------------------------------------------
!> note: on input G, P are assumed to be in cartesian GTOs
  SUBROUTINE interest_eri_diff(ndim, Gmat, Pmat, iblocks)

    integer, intent(in)  :: ndim
    real(8), intent(in)  :: Pmat(ndim, ndim, 4)
    real(8), intent(out) :: Gmat(ndim, ndim, 4)
    integer, intent(in)  :: iblocks(4)

#ifdef PRG_DIRAC
    !> local
    integer :: i, j, k, l
    integer :: li, ni, oi
    integer :: lj, nj, oj
    integer :: lk, nk, ok
    integer :: ll, nl, ol
    integer :: nint, lmax, ibas
    integer :: nij, nkl, ij, kl
    integer :: iz
    integer :: i_start, i_end
    integer :: j_start, j_end
    integer :: k_start, k_end
    integer :: l_start, l_end
    real(8) :: ei, ci, xi, yi, zi
    real(8) :: ej, cj, xj, yj, zj
    real(8) :: ek, ck, xk, yk, zk
    real(8) :: el, cl, xl, yl, zl

    real(8) :: gout(441*441)  !todo: better definition
    real(8) :: fij, fkl, fijkl


    !> InteRest interface
    interface
      subroutine interest_eri_basic(factor,fint,nint,       &
                                    la,alpha,ax,ay,az,anorm,&
                                    lb,beta ,bx,by,bz,bnorm,&
                                    lc,gamma,cx,cy,cz,cnorm,&
                                    ld,delta,dx,dy,dz,dnorm )

        integer, intent(out) :: nint
        real(8), intent(out) :: fint(*)
        real(8), intent(in ) :: factor
        integer, intent(in ) :: la,lb,lc,ld
        real(8), intent(in ) :: alpha,ax,ay,az,anorm
        real(8), intent(in ) :: beta, bx,by,bz,bnorm
        real(8), intent(in ) :: gamma,cx,cy,cz,cnorm
        real(8), intent(in ) :: delta,dx,dy,dz,dnorm
      end subroutine
    end interface

    i_start = shell_start(iblocks(1))
    i_end   = shell_end(iblocks(1))

    j_start = shell_start(iblocks(2))
    j_end   = shell_end(iblocks(2))

    k_start = shell_start(iblocks(3))
    k_end   = shell_end(iblocks(3))

    l_start = shell_start(iblocks(4))
    l_end   = shell_end(iblocks(4))

       iloop: do i = i_start, i_end
         li =       gto(i)%lvalue   !s=1, p=2, d=3 (l+1)
         oi =       gto(i)%offset   !example: two p functions, offsets are [0, 3]
         ni =       gto(i)%cdegen   !cartesian degenracy s=1, p=3, d=6, f=10
         ei =       gto(i)%exponent(1)
         xi = atom( gto(i)%origin )%coordinate_x
         yi = atom( gto(i)%origin )%coordinate_y
         zi = atom( gto(i)%origin )%coordinate_z
         ci =       gto(i)%coefficient(1)

         jloop: do j = j_start, j_end
           lj =       gto(j)%lvalue
           oj =       gto(j)%offset
           nj =       gto(j)%cdegen
           ej =       gto(j)%exponent(1)
           xj = atom( gto(j)%origin )%coordinate_x
           yj = atom( gto(j)%origin )%coordinate_y
           zj = atom( gto(j)%origin )%coordinate_z
           cj =       gto(j)%coefficient(1)

           kloop: do k = k_start, k_end
             lk =       gto(k)%lvalue
             ok =       gto(k)%offset
             nk =       gto(k)%cdegen
             ek =       gto(k)%exponent(1)
             xk = atom( gto(k)%origin )%coordinate_x
             yk = atom( gto(k)%origin )%coordinate_y
             zk = atom( gto(k)%origin )%coordinate_z
             ck =       gto(k)%coefficient(1)

             lloop: do l = l_start, l_end
               ll   =       gto(l)%lvalue
               nl   =       gto(l)%cdegen
               ol   =       gto(l)%offset
               el   =       gto(l)%exponent(1)
               xl   = atom( gto(l)%origin )%coordinate_x
               yl   = atom( gto(l)%origin )%coordinate_y
               zl   = atom( gto(l)%origin )%coordinate_z
               cl   =       gto(l)%coefficient(1)
               fijkl= 1.0d0

               ! fijkl is scaling constant
               ! gout  is output (integrals)
               ! [ij|kl] [electron1|electron2]
               ! gout = (ncc(k), ncc(l), ncc(i), ncc(j))
               ! example: [sp|df]: [6, 10, 1, 3] this is the layout in mem
               ! limitation: up to h functions (incl)

               !> call InteRest library routine for a given batch
               call interest_eri_basic(fijkl,gout,nint,  &
                                       li,ei,xi,yi,zi,ci,&
                                       lj,ej,xj,yj,zj,cj,&
                                       lk,ek,xk,yk,zk,ck,&
                                       ll,el,xl,yl,zl,cl )

!              !> process integral batch with the density matrix
               call process_dG( ni, nj, nk, nl, oi, oj, ok, ol, gout, &
                                Pmat, Gmat, ndim, ndim, .true., 1.0d0             )

             enddo lloop
           enddo kloop
         enddo jloop
       enddo iloop
#endif /* ifdef PRG_DIRAC */

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!
#ifdef PRG_DIRAC
  SUBROUTINE process_dG( ic, jc, kc, lc, io, jo, ko, lo, gout, dP, dG, n1, n2, doK, scaleK )

    !> input
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    logical, intent(in) :: doK
    real(8), intent(in) :: scaleK
    integer, intent(in) :: ic, jc
    integer, intent(in) :: kc, lc
    integer, intent(in) :: io, jo, ko, lo
    real(8), intent(in) :: gout(kc,lc,ic,jc)
    real(8), intent(in) :: dP(n1,n1,*)

    !> output
    real(8), intent(out) :: dG(n2,n2,*)

    !> local
    integer :: n, ns
    integer :: i, j, k, l
    real(8) :: dPkl, dPkj(4)
    integer :: ibas, jbas, kbas, lbas


    !> coulomb part
      do l=1,lc
        lbas=lo+l
        do k=1,kc
          kbas=ko+k
          dPkl=dP(kbas,lbas,1)*2.0d0
          do j=1,jc
            jbas=jo+j
            do i=1,ic
              ibas=io+i
              dG(ibas,jbas,1) = dG(ibas,jbas,1) + gout(k,l,i,j)*dPkl
            enddo
          enddo
        enddo
      enddo

    !> exchange part
    if( doK )then
        do l=1,lc
          lbas=lo+l
          do k=1,kc
            kbas=ko+k
            do j=1,jc
              jbas=jo+j
              !P is symmetric
              dPkj(1) = dP(jbas,kbas,1)*scaleK
              dPkj(2) = dP(jbas,kbas,2)*scaleK
              dPkj(3) = dP(jbas,kbas,3)*scaleK
              dPkj(4) = dP(jbas,kbas,4)*scaleK
              do i=1,ic
                ibas=io+i
                dG(ibas,lbas,1) = dG(ibas,lbas,1) - gout(k,l,i,j)*dPkj(1)
                dG(ibas,lbas,2) = dG(ibas,lbas,2) - gout(k,l,i,j)*dPkj(2)
                dG(ibas,lbas,3) = dG(ibas,lbas,3) - gout(k,l,i,j)*dPkj(3)
                dG(ibas,lbas,4) = dG(ibas,lbas,4) - gout(k,l,i,j)*dPkj(4)
              enddo
            enddo
          enddo
        enddo
    endif

  END SUBROUTINE
#endif /* ifdef PRG_DIRAC */

end module
