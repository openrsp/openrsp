module interface_basis

   use basis_set, only: cgto
   use matrix_lowlevel, only : matrix, mat_init
   use eri_contractions, only: ctr_arg, set_eri_contractions_xfac
   use eri_basis_loops, only: unopt_geodiff_loop
   implicit none

#ifndef VAR_LSDALTON
   !dependent on dalton/dirac commonblocks
   public interface_basis_init
#endif
   public interface_basis_init_general
   public interface_basis_main_general
   public interface_basis_xfac_general
   public interface_basis_init_argument_general
   public get_nr_ao
   public interface_basis_finalize

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   type(cgto), pointer, public :: basis_large(:)
   real(8), allocatable        :: exp_and_ctr_large(:)
   type(cgto), pointer, public :: basis_small(:)
   real(8), allocatable        :: exp_and_ctr_small(:)
   type(ctr_arg) :: ctr_arg_item(1)
!  non-allocatables
   integer :: nr_ao

contains

  !A proper interface where the required info is provided as primitive
  !Arguments to the subroutine (NOT through common block specific for the 
  !host program). TK
  subroutine interface_basis_init_general(nbast,nr_blocks_large,nEXPSIZE,nCCSIZE,nCC,&
       & ANGMOM,CCINDEX,ICHARGE,IBASIS,nPrim,nCont,start_exponents,start_CC,&
       & CENTER,exponents,ContC,MAXMOM,AMOM)
     implicit none
     integer,intent(in) :: nbast,nr_blocks_large,nEXPSIZE,nCCSIZE,nCC,MAXMOM
     integer,intent(in) :: ANGMOM(nr_blocks_large),CCINDEX(nr_blocks_large)
     integer,intent(in) :: ICHARGE(nr_blocks_large),IBASIS(nr_blocks_large),AMOM(nCC)
     integer,intent(in) :: nPrim(nCC),nCont(nCC),start_exponents(nCC),start_CC(nCC)
     real(8),intent(in) :: CENTER(3,nr_blocks_large)
     real(8),target :: exponents(nEXPSIZE),ContC(nCCSIZE)
     integer :: i,maxm,nelms,istartCC,istart,nrow,ncont1,j
     real(8) :: rescal(MAXMOM+1),factor,TMPContC(nCCSIZE)
     ! coefficient rescaling factors
     maxm = 0
     rescal(maxm+1) = 2*sqrt(acos(-1d0))
     do I=1,MAXMOM
        maxm = maxm + 1
        rescal(maxm+1) = rescal(maxm) / sqrt(2*maxm+1.d0)
     end do
     do i=1,nCC
        nelms = nPrim(CCINDEX(I))*nCont(CCINDEX(I))
        istartCC = start_CC(CCINDEX(I))
        factor = rescal(amom(i)+1)
        do j=istartCC+1,istartCC+nelms
           TMPContC(j) = factor*ContC(j)
        enddo
     enddo
     do i=1,nCCSIZE
        ContC(i) = TMPContC(i)
     enddo

     nr_ao = nbast
     nullify(basis_large)
     nullify(basis_small)
     allocate(basis_large(nr_blocks_large))
     allocate(exp_and_ctr_large(1)) !not used.
     do I = 1,nr_blocks_large
        basis_large(I)%CENT(1) = CENTER(1,I)
        basis_large(I)%CENT(2) = CENTER(2,I)
        basis_large(I)%CENT(3) = CENTER(3,I)
        basis_large(I)%MOM = ANGMOM(I)
        istart = start_exponents(CCINDEX(I))
        nrow = nPrim(CCINDEX(I))
        basis_large(I)%exp => exponents(istart+1:istart+nrow)
        istartCC = start_CC(CCINDEX(I))
        ncont1 = nCont(CCINDEX(I))
        nelms = nrow*nCont1
        basis_large(I)%nbas   = nCont1 * (2*basis_large(I)%mom + 1)
        call point(ContC(istartCC+1:istartCC+nelms), basis_large(I)%ctr,nrow,nCont1)
        basis_large(I)%CHARGE = ICHARGE(I)
        basis_large(I)%iBAS = IBASIS(I)
        basis_large(I)%iCENT = -1 ! I do not think this should be used
     enddo
     is_initialized = .true.
  contains
    ! point to re-ranked array
    subroutine point(ctr, ctr_pt,ne,nc)
     implicit none
      integer,  intent(in) :: ne,nc
      real(8), target,  intent(in)  :: ctr(nc,ne)
      real(8), pointer, intent(out) :: ctr_pt(:,:)
      ctr_pt => ctr
    end subroutine
   end subroutine

   subroutine interface_basis_init_argument_general(Dfull,average,ncor,geo,nbast)
     implicit none
     integer,intent(in) :: ncor,geo,nbast
!     type(matrix),intent(in) :: Dmat
!     type(matrix),intent(inout) :: Gmat
     real(8),intent(in) :: Dfull(nbast,nbast)
     real(8),target,intent(inout)  :: average(ncor)
     ctr_arg_item(1)%geo = geo
     ctr_arg_item(1)%comp = -huge(1)
     ctr_arg_item(1)%ncor = ncor
     nullify(ctr_arg_item(1)%dens)
     allocate(ctr_arg_item(1)%dens)
     call mat_init(ctr_arg_item(1)%dens, nbast, nbast, &
                   .false., .false., .false., .false., .false.)
     CALL DCOPY(NBAST*NBAST,Dfull,1,ctr_arg_item(1)%dens%elms,1)
     !
     nullify(ctr_arg_item(1)%fock_or_dens)
     allocate(ctr_arg_item(1)%fock_or_dens)
     call mat_init(ctr_arg_item(1)%fock_or_dens, nbast, nbast, &
                   .false., .false., .false., .false., .false.)
     CALL DCOPY(NBAST*NBAST,Dfull,1,ctr_arg_item(1)%fock_or_dens%elms,1)
!     ctr_arg_item(1)%dens = Dmat
!     ctr_arg_item(1)%fock_or_dens = Gmat
     ctr_arg_item(1)%average => average
   end subroutine 

   subroutine interface_basis_main_general()
     implicit none
     call unopt_geodiff_loop(basis_large, basis_small, ctr_arg_item)
   end subroutine 

   subroutine interface_basis_xfac_general(ExchangeFactor)
     implicit none
     real(8) :: ExchangeFactor
     call set_eri_contractions_xfac(ExchangeFactor)
   end subroutine 


#ifndef VAR_LSDALTON
   !dependent on dalton/dirac commonblocks
   subroutine interface_basis_init()

      integer :: nr_blocks_large  = 0
      integer :: nr_exp_ctr_large = 0
      integer :: nr_blocks_small  = 0
      integer :: nr_exp_ctr_small = 0

      integer :: i

! need MXSHEL
#include "maxorb.h"
! uses NLRGSH, NSMLSH
#include "shells.h"

#ifdef PRG_DALTON
! uses nbast
#include "inforb.h"
#endif

#ifdef PRG_DIRAC
! uses ntbas(0)
#include "dcbbas.h"
#endif

#ifdef PRG_DALTON
      nr_ao = nbast
#endif

#ifdef PRG_DIRAC
      nr_ao = ntbas(0)
#endif

      nullify(basis_large)
      nullify(basis_small)

      call shells_find_sizes(nr_blocks_large,  &
                             nr_exp_ctr_large, &
                             1,                &
                             nlrgsh)

      allocate(basis_large(nr_blocks_large))
      allocate(exp_and_ctr_large(nr_exp_ctr_large))

      call shells_to_type_cgto(nr_blocks_large,   &
                               nr_exp_ctr_large,  &
                               exp_and_ctr_large, &
                               basis_large,       &
                               1,                 &
                               nlrgsh)

#ifdef PRG_DIRAC
      call shells_find_sizes(nr_blocks_small,  &
                             nr_exp_ctr_small, &
                             nlrgsh + 1,       &
                             nlrgsh + nsmlsh)

      allocate(basis_small(nr_blocks_small))
      allocate(exp_and_ctr_small(nr_exp_ctr_small))

      call shells_to_type_cgto(nr_blocks_small,   &
                               nr_exp_ctr_small,  &
                               exp_and_ctr_small, &
                               basis_small,       &
                               nlrgsh + 1,        &
                               nlrgsh + nsmlsh)
#endif

      is_initialized = .true.

   end subroutine
#endif

   subroutine interface_basis_finalize()

      deallocate(basis_large)
      nullify(basis_large)
      deallocate(exp_and_ctr_large)

#ifdef PRG_DIRAC
      deallocate(basis_small)
      nullify(basis_small)
      deallocate(exp_and_ctr_small)
#endif

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_basis'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   integer function get_nr_ao()
      call check_if_interface_is_initialized()
      get_nr_ao = nr_ao
   end function

#ifndef VAR_LSDALTON
   !dependent on dalton/dirac commonblocks
  !> Count the number of contracted Gaussian-type orbital shells
  !> \param ncgto, and number of exponents and contraction coefficents
  !> \param nectr, so that memory for data structures can be allocated
  subroutine SHELLS_find_sizes(ncgto,   &
                               nectr,   &
                               i_start, &
                               i_end)

    integer, intent(out) :: ncgto
    integer, intent(out) :: nectr
    integer, intent(in)  :: i_start
    integer, intent(in)  :: i_end

    ! need MXSHEL
#include "maxorb.h"
    ! need NBCH NUCO NRCO
#include "shells.h"
    logical haveit(MXSHEL)
    integer i, j
    ! count the number of cgto blocks and the number of cgto
    ! exponents (from +1 below) and contraction coefficients
    ncgto = 0
    nectr = 0
    haveit(:) = .false.
    do i = i_start, i_end
       ! if not the first contacted in this block, skip
       if (NUMCF(i) /= 1) cycle
       ! index of AO shell block
       j = NBCH(i)
       if (j <= 0 .or. j > MXSHEL) &
          call quit('SHELLS_sizes error: unexpected NBCH, <=0 or >MXSHEL')
       ! count as cgto
       ncgto = ncgto + 1
       ! if first occurance of this block, count exponents and
       ! contraction coefficients
       if (haveit(j)) cycle
       haveit(j) = .true.
       nectr = nectr + NUCO(i)*(NRCO(i)+1)
    end do
  end subroutine
#endif

#ifndef VAR_LSDALTON
  !dependent on dalton/dirac commonblocks
  subroutine SHELLS_to_type_cgto(ncgto,   &
                                 nectr,   &
                                 ectr,    &
                                 bas,     &
                                 i_start, &
                                 i_end)
    use basis_set, only: cgto
    integer,         intent(in)  :: ncgto
    integer,         intent(in)  :: nectr
    type(cgto),      intent(out) :: bas(ncgto)
    real(8), target, intent(out) :: ectr(nectr)
    integer,         intent(in)  :: i_start
    integer,         intent(in)  :: i_end
    ! need MXSHEL
#include "maxorb.h"
    ! need NBCH NUCO NRCO NUMCF CENT NHKT NSTRT
#include "shells.h"
    ! need MXCONT
#include "aovec.h"
    ! need PRIEXP PRICCF
#include "primit.h"
    ! need MAXQNM
#include "maxmom.h"
    integer alreadyis(MXSHEL)
    real(8) rescal(MXQNM+1)
    integer i, j, k, l, m, ne, nc, maxm
    ! coefficient rescaling factors
    maxm = 0
    rescal(maxm+1) = 2*sqrt(acos(-1d0))
    k = 0 !index in bas
    l = 0 !offset in ectr
    alreadyis(:) = 0
    do i = i_start, i_end
       if (NBCH(i) <= 0 .or. NBCH(i) > MXSHEL) &
          call quit('SHELLS_sizes error: unexpected NBCH<=0 or >MXSHEL')
       ! if not the first contracted in this block, skip
       if (NUMCF(i) /= 1) cycle
       ! index of AO shell block
       j = NBCH(i)
       ! fill cgto entry
       k = k + 1
       ne = NUCO(i) !number of exponents
       nc = NRCO(i) !number of contracted
       bas(k)%cent   = CENT(i,:,1)
       bas(k)%mom    = NHKT(i) - 1
       bas(k)%charge = huge(1) !ajt FIXME wrong
       bas(k)%icent  = NCENT(i)
       bas(k)%nbas   = nc * (2*bas(k)%mom + 1)
       bas(k)%ibas   = NSTRT(i)
       if (alreadyis(j) /= 0) then
          bas(k)%exp => bas(alreadyis(j))%exp(:)
          bas(k)%ctr => bas(alreadyis(j))%ctr(:,:)
       else
          alreadyis(j) = k
          ! copy exponents from PRIEXP to ectr
          j = JSTRT(i)
          ectr(l+1:l+ne) = PRIEXP(j+1:j+ne)
          bas(k)%exp => ectr(l+1:l+ne)
          l = l + ne
          ! copy and rescale contraction coefs from PRICCF to ectr
          do while (maxm < bas(k)%mom)
             maxm = maxm + 1
             rescal(maxm+1) = rescal(maxm) / sqrt(2*maxm+1d0)
          end do
          ectr(l+1:l+ne*nc) = rescal(bas(k)%mom+1) * reshape( &
                              transpose(PRICCF(j+1:j+ne,:nc)),(/ne*nc/))
          call point(ectr(l+1:l+ne*nc), bas(k)%ctr)
          l = l + ne*nc
       end if
    end do
  contains
    ! point to re-ranked array
    subroutine point(ctr, ctr_pt)
      real(8), target,  intent(in)  :: ctr(nc,ne)
      real(8), pointer, intent(out) :: ctr_pt(:,:)
      ctr_pt => ctr
    end subroutine
  end subroutine
#endif

end module
