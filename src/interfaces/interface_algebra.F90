module interface_algebra

   implicit none

   public get_nz

   private

   integer :: nz = 1

contains

   integer function get_nz()
      get_nz = nz
   end function

end module
