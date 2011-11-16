!> @file
!> Contains module rsp_functions

!> The response functions module is OpenRSP's front-end and calculates
!> response functions. Remember that module rsp_backend must be correctly
!> set up before calling rsp_function
module rsp_functions

  use matrix_genop
  use rsp_contribs
  
  implicit none

  public rsp_function
  public rsp_func_tensor
  private

contains

  !> a single response function component
  function rsp_function(fields)
    type(rsp_field), intent(in) :: fields(:)
    complex(8)                  :: rsp_function
    rsp_function = 0
  end function

  !> a response function tensor
  subroutine rsp_func_tensor(n, fields, dims, rsp)
    integer,         intent(in)  :: n
    type(rsp_field), intent(in)  :: fields(n)
    integer,         intent(in)  :: dims(n)
    complex(8),      intent(out) :: rsp(*) !any-rank, size=product(dims)
    rsp(:product(dims)) = 0
  end subroutine

end module
