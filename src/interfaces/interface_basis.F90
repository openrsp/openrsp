module interface_basis

   use basis_set, only: cgto

   implicit none

   public interface_basis_init
   public interface_basis_finalize
   public get_nr_ao

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   type(cgto), pointer, public :: basis_large(:)
   real(8), allocatable        :: exp_and_ctr_large(:)
#ifdef PRG_DIRAC
   type(cgto), pointer, public :: basis_small(:)
   real(8), allocatable        :: exp_and_ctr_small(:)
#endif

!  non-allocatables
   integer :: nr_ao

contains

   subroutine interface_basis_init()

      integer :: nr_blocks_large  = 0
      integer :: nr_exp_ctr_large = 0

#ifdef PRG_DIRAC
      integer :: nr_blocks_small  = 0
      integer :: nr_exp_ctr_small = 0
#endif

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

      call shells_find_sizes(nr_blocks_large,  &
                             nr_exp_ctr_large, &
                             1,                &
                             nlrgsh)

      nullify(basis_large)
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

      nullify(basis_small)
      allocate(basis_small(nr_blocks_small))
      allocate(exp_and_ctr_small(nr_exp_ctr_small))

      call shells_to_type_cgto(nr_blocks_small,   &
                               nr_exp_ctr_small,  &
                               exp_and_ctr_small, &
                               basis_small,       &
                               nlrgsh + 1,        &
                               nlrgsh + nsmlsh)

      print *, 'debug: large component basis set'
      print *, '--------------------------------'
      print *, 'nr blocks: ', nr_blocks_large
      do i = 1, nr_blocks_large
         print *, i, basis_large(i)%mom, basis_large(i)%nbas
         print *, 'ibas: ', basis_large(i)%ibas
         print *, 'exp:  ', basis_large(i)%exp
         print *, 'ctr:  ', basis_large(i)%ctr
      end do

      print *, 'debug: small component basis set'
      print *, '--------------------------------'
      print *, 'nr blocks: ', nr_blocks_small
      do i = 1, nr_blocks_small
         print *, i, basis_small(i)%mom, basis_small(i)%nbas
         print *, 'ibas: ', basis_small(i)%ibas
         print *, 'exp:  ', basis_small(i)%exp
         print *, 'ctr:  ', basis_small(i)%ctr
      end do
#endif

      is_initialized = .true.

   end subroutine

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

end module
