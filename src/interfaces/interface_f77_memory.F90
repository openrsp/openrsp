module interface_f77_memory

   implicit none

   public interface_f77_memory_init
   public interface_f77_memory_finalize

   public f77_memory_select
   public f77_memory_deselect

   private

   real(8), pointer, public :: f77_work(:)

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   integer :: f77_work_len  !length of the whole array
   integer :: f77_work_next !position of the non-used work array
   integer :: f77_work_left !amount of the left work array

contains

   subroutine interface_f77_memory_init(work_len, work)

      integer, intent(in) :: work_len
      real(8), target     :: work(:)

      f77_work_len  = work_len
      f77_work_next = 1
      f77_work_left = work_len

      f77_work => work

      is_initialized = .true.

   end subroutine

   subroutine interface_f77_memory_finalize()

      f77_work_len  = 0
      f77_work_next = 0
      f77_work_left = 0

      nullify(f77_work)

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_f77_memory'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

!  integer function get_nr_ao()
!     call check_if_interface_is_initialized()
!     get_nr_ao = nr_ao
!  end function

   subroutine f77_memory_select(work_len, work)

      integer, intent(in) :: work_len
      real(8), pointer    :: work(:)

      call check_if_interface_is_initialized()

      if (work_len > f77_work_left) then
         print *, 'error: work_len > f77_work_left in f77_memory_select'
         stop 1
      else
         f77_work_next = f77_work_next + work_len
         f77_work_left = f77_work_left - work_len
         work => f77_work
      end if

   end subroutine

   subroutine f77_memory_deselect(work_len, work)

      integer, intent(in) :: work_len
      real(8), pointer    :: work(:)

      call check_if_interface_is_initialized()

      f77_work_next = f77_work_next - work_len
      f77_work_left = f77_work_left + work_len
      nullify(work)

   end subroutine

end module
