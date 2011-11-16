! non-copyrighted template

!> @file
!> Contains module rsp_backend

!> The backend module provides openrsp with access to the host program's
!> self-consistent density, Fock and overlap matrices, to the number of
!> atoms in the molecule, whether there is E_xc, and relays calls for
!> 1el-, 2el-, and XC-integrals, nuclear potential contributions,
!> and calls to the response solver
module rsp_backend

  use matrix_genop

  implicit none

  ! public oneint?
  public twoint
  ! public xcint?
  ! public nucpot?
  ! public oneint_arg?
  public twoint_arg
  ! public xcint_arg?
  private


  !> Type used for vectorizing and 'storing' calls to two-electron program.
  !> To prevent tampering, this should be moved into another, common file
  type twoint_arg
     !> order/index of geometry
     integer               :: geo
     !> order/index of London
     integer               :: lon
     !> references to matrices involed. If Fa present and not Db Gb(Da) is
     !> added to Fa. Should be properly allocated
     type(matrix), pointer :: Da, Db, Fa
     !> if Db present and not Fa, average trG(Da)Db goes here
     real(8)               :: ave
  end type

contains

  !> Run two-electron program to calculate a certain batch of
  !> Fock matrices and/or averages. Successive args must come in 'order'
  subroutine twoint(args)
    type(twoint_arg), intent(inout) :: args(:)
    
  end subroutine

end module
