module interface_io

!  interface_io_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_io_init

   implicit none

   public interface_io_init
   public interface_io_finalize

   public get_print_unit
   public get_input_unit

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   integer :: print_unit
   integer :: input_unit

contains

#ifdef VAR_LSDALTON
  !A proper interface that recieves info as primitve arguments to subroutine
  !and is therefore independent on the Host program. TK
   subroutine interface_io_init(lupri,lucmd)
     implicit none
     integer,intent(in) :: lupri,lucmd
      print_unit = lupri
      input_unit = lucmd
      is_initialized = .true.
   end subroutine
#else
  !Crappy interrface that uses host specific commen blocks. 
   subroutine interface_io_init()
#include "priunit.h"
      print_unit = lupri
      input_unit = lucmd
      is_initialized = .true.
   end subroutine
#endif

   subroutine interface_io_finalize()
      is_initialized = .false.
   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_io'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   integer function get_print_unit()
      call check_if_interface_is_initialized()
      get_print_unit = print_unit
   end function

   integer function get_input_unit()
      call check_if_interface_is_initialized()
      get_input_unit = input_unit
   end function

end module
