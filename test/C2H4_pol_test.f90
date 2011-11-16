program test
  use assertion     !subroutines assert*, assertion_threshold
  use rsp_backend   !in openrsp@repo this is the file-loading testing backend
  use rsp_functions !library's front end
  implicit none
  ! value of the XX-component of the negative polarizability
  complex(8) val
  ! this will instruct the backend to load files from subdirectory C2H4_exp 
  call rsp_backend_set_case('C2H4_exp')
  ! rsp_function function returns a single component
  val = rsp_function((/rsp_field('EL  ',1,(0d0,0d0)), &
                       rsp_field('EL  ',1,(0d0,0d0))/))
  ! assert(a,b,msg) compares a and b, and either prints '.' indicating success,
  ! or msg indicating failure. Threshold can be modified
  call assert(val, (-1.054432d0,0d0), '-porarizability_xx wrong')
end program
