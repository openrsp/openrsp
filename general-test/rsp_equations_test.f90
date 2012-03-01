program test
  use assertion     !subroutines assert*, assertion_threshold
  use rsp_backend   !which in openrsp@repo is the file-loading testing backend
  use rsp_contribs  !integrals/Fock
  use rsp_equations !for rsp_dens
  implicit none
  ! the response density and the residual contracted from it
  type(matrix) Df, R
  ! this will instruct the backend to load files from subdirectory C2H4_exp 
  call rsp_backend_set_case('H2O2_exp')
  call rsp_dens((/rsp_field('EL  ',1,1,(0d0,0d0))/), Df)
  ! from Df, Hf and G(Df), contract right hand side manually,
  ...
  ! assert matrix R is small
  call assert_zero(R, 'Response equation for EL_x was not solved correctly')
end program
