module interface_pcm

!  interface_pcm_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_pcm_init

   implicit none

   public interface_pcm_init
   public interface_pcm_finalize

   public get_is_pcm_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_pcm_calculation

contains

   subroutine interface_pcm_init(l)

      logical, intent(in) :: l

      is_pcm_calculation = l

      is_initialized = .true.

   end subroutine

   subroutine interface_pcm_finalize()

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_pcm'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   logical function get_is_pcm_calculation()
      call check_if_interface_is_initialized()
      get_is_pcm_calculation = is_pcm_calculation
   end function

end module
