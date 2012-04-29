module interface_f77_memory

   implicit none

   public interface_f77_memory_init
   public interface_f77_memory_finalize

   public f77_memory_select
   public f77_memory_deselect

   public get_f77_memory_total
   public get_f77_memory_next
   public get_f77_memory_left

   public set_f77_memory_next

   private

   real(8), pointer, public :: f77_memory(:)

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   integer :: f77_memory_total !length of the whole array
   integer :: f77_memory_next  !position of the non-used work array
   integer :: f77_memory_left  !amount of the left work array

contains

   subroutine interface_f77_memory_init(work_len, work)

      integer, intent(in) :: work_len
      real(8), target     :: work(:)

      f77_memory_total = work_len
      f77_memory_next  = 1
      f77_memory_left  = work_len

      f77_memory => work

      is_initialized = .true.

   end subroutine

   subroutine interface_f77_memory_finalize()

      f77_memory_total = 0
      f77_memory_next  = 0
      f77_memory_left  = 0

      nullify(f77_memory)

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_f77_memory'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   subroutine f77_memory_select(work_len, work)

      integer, intent(in) :: work_len
      real(8), pointer    :: work(:)

      call check_if_interface_is_initialized()

      if (work_len > f77_memory_left) then
         print *, 'error: work_len > f77_memory_left in f77_memory_select'
         stop 1
      else
         f77_memory_next = f77_memory_next + work_len
         f77_memory_left = f77_memory_left - work_len
         work => f77_memory
      end if

   end subroutine

   subroutine f77_memory_deselect(work_len, work)

      integer, intent(in) :: work_len
      real(8), pointer    :: work(:)

      call check_if_interface_is_initialized()

      f77_memory_next = f77_memory_next - work_len
      f77_memory_left = f77_memory_left + work_len
      nullify(work)

   end subroutine

   integer function get_f77_memory_total()
      call check_if_interface_is_initialized()
      get_f77_memory_total = f77_memory_total
   end function

   integer function get_f77_memory_next()
      call check_if_interface_is_initialized()
      get_f77_memory_next = f77_memory_next
   end function

   integer function get_f77_memory_left()
      call check_if_interface_is_initialized()
      get_f77_memory_left = f77_memory_left
   end function

   subroutine set_f77_memory_next(i)

      integer, intent(in) :: i

      call check_if_interface_is_initialized()

      f77_memory_next = i
      f77_memory_left = f77_memory_total - i + 1

   end subroutine


end module
