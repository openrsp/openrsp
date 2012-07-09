module interface_basis

   use basis_set, only: cgto

   implicit none

   public interface_basis_init
   public interface_basis_finalize
   public get_nr_ao

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   type(cgto), pointer, public :: interface_basis_pointer(:)
   real(8), allocatable :: exp_and_ctr(:)

!  non-allocatables
   integer :: nr_ao

contains

   subroutine interface_basis_init()

      integer              :: num_cgto_blocks
      integer              :: num_exp_and_ctr

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

      call shells_find_sizes(num_cgto_blocks, num_exp_and_ctr)

      nullify(interface_basis_pointer)
      allocate(interface_basis_pointer(num_cgto_blocks))
      allocate(exp_and_ctr(num_exp_and_ctr))

      call shells_to_type_cgto(num_cgto_blocks, &
                               num_exp_and_ctr, &
                               exp_and_ctr,     &
                               interface_basis_pointer)

!     do not deallocate!!! otherwise exponents will point to nirvana
!     deallocate(exp_and_ctr)

      is_initialized = .true.

   end subroutine

   subroutine interface_basis_finalize()

      deallocate(interface_basis_pointer)
      nullify(interface_basis_pointer)

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
  subroutine SHELLS_find_sizes(ncgto, nectr)
    integer, intent(out) :: ncgto, nectr
    ! need MXSHEL
#include "maxorb.h"
    ! need NLRGSH NBCH NUCO NRCO
#include "shells.h"
    logical haveit(MXSHEL)
    integer i, j
    ! count the number of cgto blocks and the number of cgto
    ! exponents (from +1 below) and contraction coefficients
    ncgto = 0
    nectr = 0
    haveit(:) = .false.
    do i = 1, NLRGSH
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

  subroutine SHELLS_to_type_cgto(ncgto, nectr, ectr, bas)
    !radovan: fixme ugly dependency on basis_set
    use basis_set, only: cgto
    integer,          intent(in)  :: ncgto, nectr
    type(cgto),       intent(out) :: bas(ncgto)
    real(8), target, intent(out) :: ectr(nectr)
    ! need MXSHEL
#include "maxorb.h"
    ! need NLRGSH NBCH NUCO NRCO NUMCF CENT NHKT NSTRT
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
    do i = 1, NLRGSH
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
