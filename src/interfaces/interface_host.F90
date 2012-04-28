module interface_host

!  interface_host_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_host_init

   implicit none

   public interface_host_init
   public interface_host_finalize

   public get_nr_ao
   public get_nr_atoms
   public get_print_unit
   public get_input_unit
   public get_nuc_charge
   public get_nuc_xyz
   public get_is_ks_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   integer :: nr_ao
   integer :: nr_atoms
   integer :: print_unit
   integer :: input_unit
   logical :: is_ks_calculation

!  allocatables, deallocate them in interface_host_finalize
   real(8), allocatable :: nuc_charge(:)
   real(8), allocatable :: nuc_xyz(:, :)

contains

   subroutine interface_host_init()

#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"

      nr_ao      = nbast
      nr_atoms   = natoms
      print_unit = lupri
      input_unit = lucmd
      is_ks_calculation = dodft

      allocate(nuc_charge(nr_atoms))
      nuc_charge = charge(:nr_atoms)

      allocate(nuc_xyz(3, nr_atoms))
      nuc_xyz = cord(:, :nr_atoms)

      is_initialized = .true.

   end subroutine

   subroutine interface_host_finalize()
!     here deallocate everything that is allocated
      if (allocated(nuc_charge)) deallocate(nuc_charge)
      if (allocated(nuc_xyz))    deallocate(nuc_xyz)

      is_initialized = .false.
   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_host'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   integer function get_nr_ao()
      call check_if_interface_is_initialized()
      get_nr_ao = nr_ao
   end function

   integer function get_nr_atoms()
      call check_if_interface_is_initialized()
      get_nr_atoms = nr_atoms
   end function

   integer function get_print_unit()
      call check_if_interface_is_initialized()
      get_print_unit = print_unit
   end function

   integer function get_input_unit()
      call check_if_interface_is_initialized()
      get_input_unit = input_unit
   end function

   real(8) function get_nuc_charge(i)
      integer, intent(in) :: i
      call check_if_interface_is_initialized()
      get_nuc_charge = nuc_charge(i)
   end function

   real(8) function get_nuc_xyz(i, j)
      integer, intent(in) :: i
      integer, intent(in) :: j
      call check_if_interface_is_initialized()
      get_nuc_xyz = nuc_xyz(i, j)
   end function

   logical function get_is_ks_calculation()
      call check_if_interface_is_initialized()
      get_is_ks_calculation = is_ks_calculation
   end function

end module
