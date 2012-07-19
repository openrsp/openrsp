! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       interest_interface 
!
! Description:  - interfacing basic data between Dirac & InteRest program 
!               - providing two- and four-center integrals in RKB/RMB balanced basis
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
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     
!                            
! ------------------------------------------------------------------------------------

MODULE INTEREST_INTERFACE

  implicit none

  public initialize_interest 
  public process_overlap_integrals 
  public process_two_electron_integrals
  public process_kinetic_energy_integrals 
  public process_nuclear_attraction_integrals

  private

  !> matrix dimensions
     !> note: by definition nrow, ncol refer to the number of cartesian GTOs
     integer :: nshl = 1
     integer :: nrow = 0
     integer :: ncol = 0
     integer :: nr_shells = 0 
     real(8) :: prescreening_treshold = 1.d-15
  !> end

  !> definition of "constant" data type 
     type type_constant
          real(8) :: alpha2      !=(1/c)^2 
     end  type
     type(type_constant) :: constant 
  !> end

  !> definition of "symmetry" data type 
     !todo: later implementation 
     type type_symmetry
          integer :: dummy 
     end  type
     type(type_symmetry) :: symmetry 
  !> end

  !> definition of "atom" data type 
     type type_atom
          real(8) :: charge 
          real(8) :: coordinate_x
          real(8) :: coordinate_y
          real(8) :: coordinate_z
          real(8) :: gnu_exponent 
     end  type
     type(type_atom), allocatable :: atom(:)
  !> end

  !> definition of "gto" data type
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
     type(type_gto), allocatable, target :: gto(:)
  !> end

CONTAINS
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE initialize_interest()

  use codata

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
    !> end
  

    !> initialize InteRest integral package 
       call interest_initialize()
    !> end

    !> test if the large component basis is uncontracted
       if( nplrg.ne.nlarge )then
         write(6,*) ' Error: contracted basis set found => stop!'
         write(6,*) ' Note:  uncontracted basis sets are only supported'
         stop
       endif
    !> end
  
    !> interface constants
       constant = type_constant( alpha2 )
    !> end

    !> interface molecular data
       allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), &
                              cord(1,i), &
                              cord(2,i), &
                              cord(3,i), &
                              gnuexp(i)  )
       enddo
    !> end
  
    !> interface basis set data
    !> fixme: contracted basis sets
       allocate( gto(nlrgsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
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
    !> end

    !> reorder basis functions to the angular momentum 
    !> ascending order to increase the integral performance
    !> set the AO-vector offsets and shell-degeneracy
       do i=1,size(gto)-1
         do j=(i+1),size(gto)
           if( gto(j)%lvalue >= gto(i)%lvalue )cycle
           tmpgto = gto(i)
           gto(i) = gto(j)
           gto(j) = tmpgto
         enddo
       enddo
       do i=1,size(gto)-1
         gto(i+1)%offset = gto(i)%offset + gto(i)%cdegen
       enddo
    !> end

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
      !write(6,*) 
      !write(6,'(2x,a,i5)') 'Basis set informations:'
      !write(6,'(2x,a,i5)') '-----------------------'
      !write(6,*) 
      !write(6,*) 
      !write(6,'(2x,a)')    'label    atom   charge   prim   cont   L-basis   RKB     S-basis'
      !write(6,'(2x,70a1)') ('-',i=1,70)
      !write(6,'(2x,70a1)') ('-',i=1,70)

       nrow      = sum( gto(:)%cdegen)
       ncol      = sum( gto(:)%cdegen)
       nr_shells = size(gto)

       write(6,*) 
       write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
       write(6,'(2x,a,i5  )') 'Total number of spherical basis functions:',sum(gto(:)%sdegen)
       write(6,'(2x,a,i5,a)') 'Total number of cartesian basis functions:',sum(gto(:)%cdegen),' (used for calculation)'
       write(6,*) 
    !> end

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_overlap_integrals( S ) 

    !> input
       real(8) :: S(2*nrow,2*ncol) 
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Smat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_overlap 
    !> end


    allocate( Smat( nrow,ncol,2 ) ); Smat = 0.0d0
    call two_center_integrals( interest_overlap, Smat, nrow )
    call symmetrize_matrix(                                 A = Smat(1,1,1), ida = nrow )  !Re(A)
    call symmetrize_matrix( factor = constant%alpha2/4.0d0, A = Smat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        S(i     ,j     ) = Smat(i,j,1) 
        S(i+nrow,j+ncol) = Smat(i,j,2) 
      enddo
    enddo
    deallocate( Smat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>
  SUBROUTINE process_kinetic_energy_integrals( T ) 

    !> input
       real(8) :: T(2*nrow,2*ncol) 
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Tmat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_overlap 
    !> end


    allocate( Tmat( nrow,ncol,2 ) ); Tmat = 0.0d0
    call two_center_integrals( interest_overlap, Tmat, nrow )
    call symmetrize_matrix( factor = 0.5d0, A = Tmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    do j=1,ncol
      do i=1,nrow
        T(i+nrow,j     ) = T(i+nrow,j     ) + Tmat(i,j,2)
        T(i     ,j+ncol) = T(i     ,j+ncol) + Tmat(i,j,2)
        T(i+nrow,j+ncol) = T(i+nrow,j+ncol) - Tmat(i,j,2)
      enddo
    enddo
    deallocate( Tmat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE process_nuclear_attraction_integrals( V ) 
 
    !> input
    !> note: we assume nz=4 !!!
       real(8) :: V(2*nrow,2*ncol,4)
    !> end

    !> local
       integer              :: i, j 
       real(8), allocatable :: Vmat(:,:,:)
    !> end

    !> external InteRest integral routines
       external :: interest_nuclear_attraction_point_nucleus 
    !> end


    allocate( Vmat( nrow,ncol,5 ) ); Vmat = 0.0d0
    call two_center_integrals( interest_nuclear_attraction_point_nucleus, Vmat, nrow )
    call symmetrize_matrix    (                                 A = Vmat(1,1,1), ida = nrow )  !Re(A)
    call symmetrize_matrix    ( factor = constant%alpha2/4.0d0, A = Vmat(1,1,2), ida = nrow )  !Re(A) = (xx + yy + zz)
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,3), ida = nrow )  !Im(A) = (xy - yx)      
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,4), ida = nrow )  !Re(B) = (zx - xz)      
    call antisymmetrize_matrix( factor = constant%alpha2/4.0d0, A = Vmat(1,1,5), ida = nrow )  !Im(B) = (yz - zy)     
    do j=1,ncol
      do i=1,nrow
        V(i     ,j     ,1) = V(i     ,j     ,1) + Vmat(i,j,1)
        V(i+nrow,j+ncol,1) = V(i+nrow,j+ncol,1) + Vmat(i,j,2)
        V(i+nrow,j+ncol,2) = V(i+nrow,j+ncol,2) + Vmat(i,j,3)
        V(i+nrow,j+ncol,3) = V(i+nrow,j+ncol,3) + Vmat(i,j,4)
        V(i+nrow,j+ncol,4) = V(i+nrow,j+ncol,4) + Vmat(i,j,5)
      enddo
    enddo
    deallocate( Vmat )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> 
  SUBROUTINE two_center_integrals( interest_routine, final_integrals, ida ) 

    !> input
       interface 
         subroutine interest_routine(rkb,                    &
                                     fint,na,nb,nrkb,        &
                                     la,alpha,ax,ay,az,anorm,&
                                     lb,beta, bx,by,bz,bnorm,&
                                     ncentr,atomic_data      )
       
           integer, intent(out) :: na 
           integer, intent(out) :: nb 
           logical, intent(in)  :: rkb
           integer, intent(out) :: nrkb
           integer, intent(in)  :: la,lb
           integer, intent(in)  :: ncentr
           real(8), intent(out) :: fint(*)
           real(8), intent(in)  :: alpha,ax,ay,az,anorm
           real(8), intent(in)  :: beta, bx,by,bz,bnorm
           real(8), intent(in)  :: atomic_data(ncentr,5)
         end subroutine
       end interface
    !> end

    !> input
       integer, intent(in)  :: ida           !todo: better definition
       logical              :: rkb = .true.  !todo: better definition
    !> end

    !> output
       real(8) :: final_integrals(*)
    !> end

    !> local
       integer :: nrkb
       integer :: i, ni
       integer :: j, nj
       real(8) :: gout(5*441)  !todo: better definition 
    !> end


    do j=1,nr_shells
      do i=j,nr_shells
        call interest_routine(rkb,                                &
                              gout, ni, nj, nrkb,                 &
                              !> bra-function
                              gto(i)%lvalue,                      &
                              gto(i)%exponent(1),                 &
                              atom( gto(i)%origin )%coordinate_x, &
                              atom( gto(i)%origin )%coordinate_y, &
                              atom( gto(i)%origin )%coordinate_z, &
                              gto(i)%coefficient(1),              & 
                              !> ket-function
                              gto(j)%lvalue,                      &
                              gto(j)%exponent(1),                 &
                              atom( gto(j)%origin )%coordinate_x, &
                              atom( gto(j)%origin )%coordinate_y, &
                              atom( gto(j)%origin )%coordinate_z, &
                              gto(j)%coefficient(1),              &
                              !> atomic data 
                              size(atom),                         &
                              (/atom(1:size(atom))%charge,        &
                                atom(1:size(atom))%coordinate_x,  &
                                atom(1:size(atom))%coordinate_y,  &
                                atom(1:size(atom))%coordinate_z,  &
                                atom(1:size(atom))%gnu_exponent/) )

        call store_batch(gout, ni, nj,  nrkb,       &
                         final_integrals, ida, ida, &
                         gto(i)%offset,             &
                         gto(j)%offset              )

      enddo
    enddo    

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> Copy batch A -> matrix B 
  SUBROUTINE store_batch( A, ida, jda, nrkb, B, idb, jdb, ioffset, joffset )

    !> input 
       integer :: nrkb
       integer :: ioffset 
       integer :: joffset 
       integer :: ida, jda 
       integer :: idb, jdb
       real(8) :: A(5,ida,jda)
    !> end

    !> output
       real(8) :: B(idb,jdb,nrkb)
    !> end

    !> local
       integer :: i
       integer :: j
       integer :: k
    !> end


    do k=1,nrkb
      do j=1,jda
        do i=1,ida
          B(ioffset+i,joffset+j,k) = A(k,i,j) 
        enddo
      enddo
    enddo
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> note: on input G, P are assumed to be in cartesian GTOs 
  SUBROUTINE process_two_electron_integrals( Gmat, Pmat, intflg ) 

    !> input
    !> note:  we assume nz=4 !!!
    !> fixme: the case when nz!=4
       real(8), intent(in) :: Pmat(2*nrow,2*ncol,4)
       integer, intent(in) :: intflg
       interface 
         subroutine interest_eri(rkb,factor,fint,nint,   &
                                 la,alpha,ax,ay,az,anorm,&
                                 lb,beta ,bx,by,bz,bnorm,&
                                 lc,gamma,cx,cy,cz,cnorm,&
                                 ld,delta,dx,dy,dz,dnorm )
       
           integer, intent(out) :: nint
           real(8), intent(out) :: fint(*)
           logical, intent(in ) :: rkb 
           real(8), intent(in ) :: factor 
           integer, intent(in ) :: la,lb,lc,ld
           real(8), intent(in ) :: alpha,ax,ay,az,anorm
           real(8), intent(in ) :: beta, bx,by,bz,bnorm
           real(8), intent(in ) :: gamma,cx,cy,cz,cnorm
           real(8), intent(in ) :: delta,dx,dy,dz,dnorm
         end subroutine
       end interface
    !> end

    !> output
       real(8), intent(out) :: Gmat(2*nrow,2*ncol,4)
    !> end

    !> local
       integer :: i, j, k, l
       integer :: li, ni, oi   
       integer :: lj, nj, oj   
       integer :: lk, nk, ok   
       integer :: ll, nl, ol   
       integer :: nint, lmax, ibas
       integer :: nij, nkl, ij, kl
       integer :: nr_integrals_total
       integer :: nr_integrals_pased
       real(8) :: ei, ci, xi, yi, zi 
       real(8) :: ej, cj, xj, yj, zj 
       real(8) :: ek, ck, xk, yk, zk 
       real(8) :: el, cl, xl, yl, zl 

       real(8) :: screened 
       real(8) :: max_ijij
       real(8) :: max_dens
       real(8) :: max_dpij
       real(8) :: max_dpki
       real(8) :: max_dpkj
       real(8) :: gout(441*441)  !todo: better definition 
       real(8) :: fij, fkl, fijkl

       real(8), allocatable :: dPi (:)
       real(8), allocatable :: dPj (:)  
       real(8), allocatable :: dGi (:)
       real(8), allocatable :: dGj (:)  
       real(8), allocatable :: dPij(:)
       real(8), allocatable :: dGij(:)

       real(8), allocatable :: P(:,:)
       real(8), allocatable :: G(:,:)  

       integer              :: record_length
       real(8), allocatable :: dG_guess(:,:,:)
       logical              :: lbit
    !> end

    !  lbit(intflg,i2typ)=T means to calculate this integral class i2typ = 1(LL),2(SL),3(SS),4(GT);
    
   
    nr_integrals_total = 0   
    nr_integrals_pased = 0

    !> allocate working/temporary fields
       allocate( P       ( nrow,ncol             ) )
       allocate( G       ( nrow,ncol             ) )
       allocate( dPij    ( nshl*21*21            ) ) !todo: better definition 
       allocate( dGij    ( nshl*21*21            ) ) !todo: better definition 
       allocate( dPi     ( nshl*21*nrow          ) ) !todo: better definition 
       allocate( dPj     ( nshl*21*nrow          ) ) !todo: better definition 
       allocate( dGi     ( nshl*21*nrow          ) ) !todo: better definition 
       allocate( dGj     ( nshl*21*nrow          ) ) !todo: better definition 
       allocate( dG_guess( 2,nr_shells,nr_shells ) )
    !> end

    !> forward resorting of the density matrix elements
       G(:,:)=0.0d0
       do j=1,ncol
         do i=1,nrow
           P(i,j) = 0.5d0*Pmat(i,j,1)  !Re(A,LL)
         enddo
       enddo
       call forward_ao_resorting( P )
    !> end

    !> precalculate maximum shell pair contributions to the Fock matrix
       call max_shell_pair_contributions( P, nrow, dG_guess )
    !> end


iloop: do i=1,nr_shells
         ij = (i-1)*nr_shells
         li =       gto(i)%lvalue
         oi =       gto(i)%offset
         ni =       gto(i)%cdegen
         ei =       gto(i)%exponent(1)
         ci =       gto(i)%coefficient(1)
         xi = atom( gto(i)%origin )%coordinate_x
         yi = atom( gto(i)%origin )%coordinate_y
         zi = atom( gto(i)%origin )%coordinate_z
     
         !> K-preparation
            ibas = oi+ni
            call get_dP_exchange( P, nrow, i, ibas,  dPi, dGi )
         !> end
       

jloop:   do j=1,i
           ij = ij+1             
           lj =       gto(j)%lvalue
           oj =       gto(j)%offset
           nj =       gto(j)%cdegen
           ej =       gto(j)%exponent(1)
           cj =       gto(j)%coefficient(1)
           xj = atom( gto(j)%origin )%coordinate_x
           yj = atom( gto(j)%origin )%coordinate_y
           zj = atom( gto(j)%origin )%coordinate_z
           fij= 1.0d0;if( j==i )fij=0.5d0
       
           !> prescreening
              max_ijij = dG_guess(1,j,i)
              max_dpij = dG_guess(2,j,i)
           !> end

           !> K-preparation
              call get_dP_exchange( P, nrow, j, ibas,  dPj, dGj )
           !> end

           !> J-preparation
              call get_dP_coulomb( P, nrow, i, j, dPij, dGij )
           !> end


kloop:     do k=1,i
             kl = (k-1)*nr_shells
             lk =       gto(k)%lvalue
             ok =       gto(k)%offset
             nk =       gto(k)%cdegen
             ek =       gto(k)%exponent(1)
             ck =       gto(k)%coefficient(1)
             xk = atom( gto(k)%origin )%coordinate_x
             yk = atom( gto(k)%origin )%coordinate_y
             zk = atom( gto(k)%origin )%coordinate_z
       
             !> prescreening
                max_dpkj = dG_guess(2,k,j)
                max_dpki = dG_guess(2,k,i)
             !> end


             lmax=k
             if( k==i )lmax=j
lloop:       do l=1,lmax
               kl   = kl+1
               ll   =       gto(l)%lvalue
               nl   =       gto(l)%cdegen
               ol   =       gto(l)%offset
               el   =       gto(l)%exponent(1)
               cl   =       gto(l)%coefficient(1)
               xl   = atom( gto(l)%origin )%coordinate_x
               yl   = atom( gto(l)%origin )%coordinate_y
               zl   = atom( gto(l)%origin )%coordinate_z
               fijkl= fij;if(  k==l  )fijkl=0.5d0*fijkl
                          if( ij==kl )fijkl=0.5d0*fijkl
       
               !> prescreening
                  nr_integrals_total = nr_integrals_total + ni*nj*nk*nl
                  max_dens = max(max_dpij,max_dpki,max_dpkj,dG_guess(2,l,j),dG_guess(2,l,i),dG_guess(2,l,k))
                  if( max_dens*max_ijij*dG_guess(1,l,k) < prescreening_treshold )cycle
                  nr_integrals_pased = nr_integrals_pased + ni*nj*nk*nl
               !> end          

               !> call InteRest library routine for a given batch
                  call interest_eri(.false.,           &
                                    fijkl,gout,nint,   &
                                    li,ei,xi,yi,zi,ci, &
                                    lj,ej,xj,yj,zj,cj, &
                                    lk,ek,xk,yk,zk,ck, &
                                    ll,el,xl,yl,zl,cl  )
               !> end
       
               !> process integral batch with the density matrix
                  call process_dG( ni, nj, nk, nl, oi, oj, ok, ol, gout,            &
                                   P, G, nrow, dPij, dGij, dPi, dPj, dGi, dGj, ibas )
               !> end 

             enddo lloop
           enddo kloop
       
           !> add J-contribution
              call add_dG_coulomb( G, nrow, dGij, ni, nj, oi, oj )
           !> end

           !> add K-contribution
              call add_dG_exchange( G, nrow, j, ibas, dGj )
           !> end

         enddo jloop

         !> add K-contribution
            call add_dG_exchange( G, nrow, i, ibas, dGi )
         !> end
       
       enddo iloop
    !> end

    !> make full G matrix 
                    call symmetrize_matrix( A = G, ida = nrow )
       if( nshl==2 )call symmetrize_matrix( A = G, ida = nrow )
    !> end


    !> prescreening
       screened = 100.d0*(1.0d0-real(nr_integrals_pased)/real(nr_integrals_total))  
    !> end

    !> backward resorting of the Fock matrix
       do j=1,ncol
         do i=1,nrow
           Gmat(i,j,1) = G(i,j)  !Re(A,LL)
         enddo
       enddo
    !> end

    deallocate( P        )
    deallocate( G        )
    deallocate( dPi      )
    deallocate( dPj      )
    deallocate( dGi      )
    deallocate( dGj      )
    deallocate( dPij     )
    deallocate( dGij     )
    deallocate( dG_guess )

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!
  SUBROUTINE get_dP_coulomb( dP, idP, ishell, jshell, dP_block, dG_block )

    !> input
       integer, intent(in) :: idP 
       integer, intent(in) :: ishell 
       integer, intent(in) :: jshell 
       real(8), intent(in) :: dP(idP,idP,nshl)
    !> end

    !> output
       real(8), intent(out) :: dP_block(*)
       real(8), intent(out) :: dG_block(*)
    !> end

    !> local
       integer :: n
       integer :: ij, jbas
       integer :: i, io, ic
       integer :: j, jo, jc
    !> end


    ij = 0
    io = gto(ishell)%offset
    ic = gto(ishell)%cdegen
    jo = gto(jshell)%offset
    jc = gto(jshell)%cdegen

    do n=1,nshl
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          dG_block(ij+i) = 0.0d0
          dP_block(ij+i) = dP(io+i,jbas,n) !note that we use here: P(ab) = P(ba)
        enddo
        ij=ij+ic
      enddo
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!
  SUBROUTINE get_dP_exchange( dP, idP, ishell, jbas, dP_block, dG_block )

    !> input
       integer, intent(in) :: idP 
       integer, intent(in) :: jbas
       integer, intent(in) :: ishell 
       real(8), intent(in) :: dP(idP,idP,nshl)
    !> end

    !> output
       real(8), intent(out) :: dP_block(*)
       real(8), intent(out) :: dG_block(*)
    !> end

    !> local
       integer :: n
       integer :: i, j, ij
       integer :: ibas, io, ic
    !> end


    ij = 0
    io = gto(ishell)%offset
    ic = gto(ishell)%cdegen

    do n=1,nshl
      do i=1,ic
        ibas = io+i
        do j=1,jbas
          dG_block(ij+j) = 0.0d0
          dP_block(ij+j) = dP(j,ibas,n) !vertical stripe (:,ibas) recorder from the UT of P=CC^{dag}
        enddo
        ij=ij+jbas
      enddo
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!
  SUBROUTINE process_dG( ic, jc, kc, lc, io, jo, ko, lo, gout, dP, dG, &
                         idP, dPij, dGij, dPi, dPj, dGi, dGj, nbas     )

    !> input
       integer, intent(in) :: idP 
       integer, intent(in) :: nbas
       integer, intent(in) :: ic, jc
       integer, intent(in) :: kc, lc
       integer, intent(in) :: io, jo, ko, lo
       real(8), intent(in) :: gout(kc,lc,ic,jc)
       real(8), intent(in) :: dPi (nbas,ic, nshl)
       real(8), intent(in) :: dPj (nbas,jc, nshl)
       real(8), intent(in) :: dPij(ic,  jc, nshl)
       real(8), intent(in) :: dP  (idP, idP,nshl)
    !> end

    !> output
       real(8), intent(out) :: dGij(ic,jc)
       real(8), intent(out) :: dGi (nbas,ic,nshl)
       real(8), intent(out) :: dGj (nbas,jc,nshl)
       real(8), intent(out) :: dG  (idP,idP,nshl)
    !> end

    !> local
       integer :: n
       integer :: i, j, k, l
       real(8) :: dPli, dPlj
       real(8) :: dGli, dGlj
       real(8) :: dPkl, dGkl
       real(8) :: coulomb_factor
       integer :: ibas, jbas, kbas, lbas
    !> end


                coulomb_factor = 4.0d0/real(nshl)
    if( ko==lo )coulomb_factor = 2.0d0*coulomb_factor

    !> coulomb part
       do n=1,nshl
         do l=1,lc
           lbas=lo+l
           do k=1,kc
             kbas=ko+k
             dGkl=0.0d0
             dPkl=dP(kbas,lbas,n)
             do j=1,jc
               do i=1,ic
                 dGkl      = dGkl      + gout(k,l,i,j)*dPij(i,j,n)
                 dGij(i,j) = dGij(i,j) + gout(k,l,i,j)*dPkl
               enddo
             enddo
             dG(kbas,lbas,n) = dG(kbas,lbas,n) + coulomb_factor*dGkl
           enddo
         enddo
       enddo
    !> end


    !> exchange part 
       do n=1,nshl
         do j=1,jc
           jbas=jo+j
           do i=1,ic
             ibas=io+i
             do l=1,lc
               lbas=lo+l 
               dGli=0.0d0
               dGlj=0.0d0
               dPlj=dPj(lbas,j,n)
               dPli=dPi(lbas,i,n)
               do k=1,kc
                 kbas=ko+k 
                             dGli          = dGli          - gout(k,l,i,j)*dPj(kbas,j,n)
                             dGi(kbas,i,n) = dGi(kbas,i,n) - gout(k,l,i,j)*dPlj
                 if( io==ko )dGi(ibas,k,n) = dGi(ibas,k,n) - gout(k,l,i,j)*dPlj
                             dGlj          = dGlj          - gout(k,l,i,j)*dPi(kbas,i,n)
                             dGj(kbas,j,n) = dGj(kbas,j,n) - gout(k,l,i,j)*dPli
                 if( jo==ko )dGj(jbas,k,n) = dGj(jbas,k,n) - gout(k,l,i,j)*dPli
               enddo
                           dGi(lbas,i,n) = dGi(lbas,i,n) + dGli
               if( io==lo )dGi(ibas,l,n) = dGi(ibas,l,n) + dGli 
                           dGj(lbas,j,n) = dGj(lbas,j,n) + dGlj 
               if( jo==lo )dGj(jbas,l,n) = dGj(jbas,l,n) + dGlj
             enddo
           enddo
         enddo
       enddo
    !> end

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> add the coulomb contribution to the G-matrix lower triangle  
  SUBROUTINE add_dG_coulomb( dG, idG, dG_block, ic, jc, io, jo )

    !> input
       integer, intent(in) :: idG 
       integer, intent(in) :: ic, jc
       integer, intent(in) :: io, jo
       real(8), intent(in) :: dG_block(ic,jc)
    !> end

    !> output
       real(8), intent(out) :: dG(idG,idG,nshl)
    !> end

    !> local
       integer :: n
       integer :: i, j, ibas, jbas
       real(8) :: coulomb_factor
    !> end

                coulomb_factor = 4.0d0/real(nshl)
    if( io==jo )coulomb_factor = 2.0d0*coulomb_factor


    do n=1,nshl
      do j=1,jc
        jbas=jo+j
        do i=1,ic
          ibas=io+i
          dG(ibas,jbas,n) = dG(ibas,jbas,n) + coulomb_factor*dG_block(i,j)
        enddo
      enddo
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> add the exchange contribution to the G-matrix lower triangle  
  SUBROUTINE add_dG_exchange( dG, idG, jshell, nbas, dG_block )

    !> input
       integer, intent(in) :: idG 
       integer, intent(in) :: nbas
       integer, intent(in) :: jshell 
       real(8), intent(in) :: dG_block(*)
    !> end

    !> output
       real(8), intent(out) :: dG(idG,idG,nshl)
    !> end

    !> local
       integer :: i, j, ij, jbas
       integer :: n, jo, jc, min
    !> end


    ij  = 0
    jo  = gto(jshell)%offset
    jc  = gto(jshell)%cdegen
    min = jo+jc 

    do n=1,nshl
      do j=1,jc
        jbas=jo+j
        do i=1,min
          dG(jbas,i,n) = dG(jbas,i,n) + dG_block(ij+i)
        enddo
        do i=(min+1),nbas
          dG(i,jbas,n) = dG(i,jbas,n) + dG_block(ij+i)
        enddo
        ij=ij+nbas
      enddo
    enddo

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!> calculate the maximum value of diagonal (i.e. [ij|ij]) integrals and differential density matrix elements over 
!> all pairs of basis functions within a given shell pair in order to estimate the shell pair contribution 
!> to the Fock matrix
  SUBROUTINE max_shell_pair_contributions( dP, idP, dG )

    !> input
       integer :: idP
       real(8) :: dP(idP,idP,nshl)
    !> end

    !> output
       real(8) :: dG(2,size(gto),size(gto))
    !> end

    !> local 
       interface
         subroutine interest_eri(rkb,                                            &
                                 factor,fint,nint,                               &
                                 la,alpha,ax,ay,az,anorm,lb,beta ,bx,by,bz,bnorm,&
                                 lc,gamma,cx,cy,cz,cnorm,ld,delta,dx,dy,dz,dnorm )

           integer, intent(out) :: nint
           real(8), intent(out) :: fint(*)
           logical, intent(in ) :: rkb
           real(8), intent(in ) :: factor
           integer, intent(in ) :: la,lb,lc,ld
           real(8), intent(in ) :: alpha,ax,ay,az,anorm
           real(8), intent(in ) :: beta, bx,by,bz,bnorm
           real(8), intent(in ) :: gamma,cx,cy,cz,cnorm
           real(8), intent(in ) :: delta,dx,dy,dz,dnorm
         end subroutine
       end interface

       integer :: n, nint
       integer :: i, ii, io, ic
       integer :: j, jj, jo, jc

       real(8) :: maximum
       real(8) :: current
       real(8) :: gout(441*441)   !todo: rozumnejsia definicia
    !> end


    !> initiate dG
       dG = 1.0d-20
    !> end

    !> find max{dsqrt(abs[ij|ij])}
iloop: do i=1,size(gto)
jloop:   do j=1,i
           !> call InteRest library routine for a given batch
              call interest_eri(.false.,                            &
                                1.0d0,gout,nint,                    &
                                !> bra-function(1)
                                gto(i)%lvalue,                      &
                                gto(i)%exponent(1),                 &
                                atom( gto(i)%origin )%coordinate_x, &
                                atom( gto(i)%origin )%coordinate_y, &
                                atom( gto(i)%origin )%coordinate_z, &
                                gto(i)%coefficient(1),              &
                                !> ket-function(1)
                                gto(j)%lvalue,                      &
                                gto(j)%exponent(1),                 &
                                atom( gto(j)%origin )%coordinate_x, &
                                atom( gto(j)%origin )%coordinate_y, &
                                atom( gto(j)%origin )%coordinate_z, &
                                gto(j)%coefficient(1),              &
                                !> bra-function(2)
                                gto(i)%lvalue,                      &
                                gto(i)%exponent(1),                 &
                                atom( gto(i)%origin )%coordinate_x, &
                                atom( gto(i)%origin )%coordinate_y, &
                                atom( gto(i)%origin )%coordinate_z, &
                                gto(i)%coefficient(1),              &
                                !> ket-function(2)
                                gto(j)%lvalue,                      &
                                gto(j)%exponent(1),                 &
                                atom( gto(j)%origin )%coordinate_x, &
                                atom( gto(j)%origin )%coordinate_y, &
                                atom( gto(j)%origin )%coordinate_z, &
                                gto(j)%coefficient(1)               )
           !> end

           maximum = dG(1,j,i) 
           do n=1,nint
             current = dsqrt(abs(gout(n)))
             if( current > maximum ) maximum = current
           enddo
           dG(1,j,i) = maximum
           dG(1,i,j) = maximum
         enddo jloop
       enddo iloop
    !> end

    !> find max{abs(dP(i,j))}
       do n=1,nshl
         do j=1,size(gto)
           jo = gto(j)%offset
           jc = gto(j)%cdegen
           do i=j,size(gto)
             io = gto(i)%offset
             ic = gto(i)%cdegen
             maximum = dG(2,i,j) 
             do jj=1,jc
               do ii=1,ic
                 current = abs( dP(io+ii,jo+jj,n) )
                 if( current > maximum ) maximum = current
               enddo
             enddo
             dG(2,i,j) = maximum
             dG(2,j,i) = maximum
           enddo
         enddo
       enddo
    !> end

  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input: lower triangle, on output: square symmetric matrix 
  SUBROUTINE symmetrize_matrix( factor, A, ida )

    !> input
    integer,           intent(in) :: ida
    real(8), optional, intent(in) :: factor

    !> output
    real(8), intent(inout) :: A(ida,ida)

    !> local
    integer :: i
    integer :: j


    if( present(factor) )then
      do j=1,ida
        do i=j,ida
          A(i,j) = A(i,j)*factor
          A(j,i) = A(i,j) 
        enddo
      enddo
    else
      do j=1,ida
        do i=j,ida
          A(j,i) = A(i,j) 
        enddo
      enddo
    endif

    call backward_ao_resorting( A )
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input: lower triangle, on output: square antisymmetric matrix 
  SUBROUTINE antisymmetrize_matrix( factor, A, ida )

    !> input
    integer,           intent(in) :: ida
    real(8), optional, intent(in) :: factor

    !> output
    real(8), intent(inout) :: A(ida,ida)

    !> local
    integer :: i
    integer :: j


    if( present(factor) )then
      do j=1,ida
        do i=j,ida
          A(i,j) =  A(i,j)*factor
          A(j,i) = -A(i,j) 
        enddo
      enddo
    else
      do j=1,ida
        do i=j,ida
          A(j,i) = -A(i,j) 
        enddo
      enddo
    endif

    call backward_ao_resorting( A )
     
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input:  the matrix in "input"  AO order 
!>on output: the matrix in integral AO order 
!>fixme:     clsed-shell and LL-only case!!!
  SUBROUTINE forward_ao_resorting( A )

    !> input/ouput
    real(8), intent(inout) :: A(nrow,ncol) 

    !> local
    real(8), allocatable :: B(:,:)
    integer :: i, ii, ic, io_new, io_old
    integer :: j, jj, jc, jo_new, jo_old


    allocate( B(nrow,ncol) )

    B = A
    do j=1,nr_shells
      jc     = gto(j)%cdegen
      jo_new = gto(j)%offset
      jo_old = gto(j)%index 
      do i=1,nr_shells
        ic     = gto(i)%cdegen
        io_new = gto(i)%offset
        io_old = gto(i)%index 
        do jj=1,jc
          do ii=1,ic
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)  !Re(A,LL) 
          enddo
        enddo
      enddo
    enddo
        
  END SUBROUTINE
! ------------------------------------------------------------------------------------
!>on input:  the matrix in "integral" AO order 
!>on output: the matrix in "input"    AO order 
!>fixme:     clsed-shell and LL-only case!!!
  SUBROUTINE backward_ao_resorting( A )

    !> input/ouput
    real(8), intent(inout) :: A(nrow,ncol) 

    !> local
    real(8), allocatable :: B(:,:)
    integer :: i, ii, ic, io_new, io_old
    integer :: j, jj, jc, jo_new, jo_old


    allocate( B(nrow,ncol) )

    B = A
    do j=1,nr_shells
      jc     = gto(j)%cdegen
      jo_new = gto(j)%index 
      jo_old = gto(j)%offset
      do i=1,nr_shells
        ic     = gto(i)%cdegen
        io_new = gto(i)%index 
        io_old = gto(i)%offset
        do jj=1,jc
          do ii=1,ic
            A(io_new+ii,jo_new+jj) = B(io_old+ii,jo_old+jj)  !Re(A,LL) 
          enddo
        enddo
      enddo
    enddo

    deallocate( B )
        
  END SUBROUTINE
! ------------------------------------------------------------------------------------
END MODULE
