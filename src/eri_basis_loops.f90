! Copyright 2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module eri_basis_loops

!> 
module eri_basis_loops

  use basis_set,        only: cgto
  use cgto_diff_eri,    only: geodiff_eri
  use eri_contractions, only: ctr_arg
  implicit none
  public unopt_geodiff_loop
  private

contains

  subroutine unopt_geodiff_loop(basis, ctrs)
    type(cgto), intent(in)  :: basis(:)
    type(ctr_arg)           :: ctrs(:)
    integer i, n, a, b, c, d, dorder
    ! initialize args' averages to zero and determine highest order
    dorder = 0
    do i = 1, size(ctrs)
       if (associated(ctrs(i)%average)) &
          ctrs(i)%average = 0
       dorder = max(dorder, ctrs(i)%geo)
    end do
    ! loop, severly unoptimized
    do d = 1, size(basis)
       do c = 1, size(basis)
          do b = 1, size(basis)
             do a = 1, size(basis)
                call geodiff_eri(dorder, basis(a), basis(b), &
                                         basis(c), basis(d), ctrs)
             end do
          end do
       end do
    end do
  end subroutine

end module
