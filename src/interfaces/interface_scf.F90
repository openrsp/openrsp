module interface_scf

!  interface_scf_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_scf_init

   implicit none

   public interface_scf_init
   public interface_scf_finalize

   public get_is_restricted_scf_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_restricted_scf_calculation

contains

   subroutine interface_scf_init()

#include "inforb.h"

      is_restricted_scf_calculation = (nasht == 0)

      is_initialized = .true.

   end subroutine

   subroutine interface_scf_finalize()

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_scf'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   logical function get_is_restricted_scf_calculation()
      call check_if_interface_is_initialized()
      get_is_restricted_scf_calculation = is_restricted_scf_calculation
   end function

end module
