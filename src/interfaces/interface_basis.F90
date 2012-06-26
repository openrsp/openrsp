module interface_basis

   implicit none

   public interface_basis_init
   public interface_basis_finalize

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

contains

   subroutine interface_basis_init()

      is_initialized = .true.

   end subroutine

   subroutine interface_basis_finalize()

      is_initialized = .false.

   end subroutine

end module
