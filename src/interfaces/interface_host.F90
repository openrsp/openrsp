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

   private

   logical :: is_initialized = .false.

   integer :: nr_ao
   integer :: nr_atoms
   integer :: print_unit
   integer :: input_unit

contains

   subroutine interface_host_init()

#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"

      nr_ao      = nbast
      nr_atoms   = natoms
      print_unit = lupri
      input_unit = lucmd

      is_initialized = .true.

   end subroutine

   subroutine interface_host_finalize()
!     here deallocate everything that is allocated
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

end module
