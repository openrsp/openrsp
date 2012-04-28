! ajt This is a backend, so should carry the same license as the
!     host project. (on the other hand, the 'template', file-loading
!     testing backend to be created in OpenRsp upstream, should be
!     un-licensed, so to make deployment of the library as easy as possible)

!> @file Contains module rsp_backend

module rsp_backend

  ! this file should eventually overtake dalton_ifc's function
  use dalton_ifc, only: DIPNUC_ifc, &
                        GRADNN_ifc, &
                        HESSNN_ifc

  use interface_host

  implicit none

  public rsp_backend_setup
  public rsp_backend_finalize
  public nuclear_potential
  public is_ks_calculation

  !--------- private -----------
  real(8), save, private, pointer :: nuc_charges(:)  !nuclear charges (natom)
  real(8), save, private, pointer :: nuc_coords(:,:) !nuclear coordinates (3,natom)

  ! needs to be down here, or doxygen will think
  private !type contents are private (bug in doxygen, maybe)
  
contains


  !> (derivatives of) nuclear repulsion and nuclei--field interaction
  subroutine nuclear_potential(derv, ncor, field, ncomp, nucpot)
    !> order of geometry derivative requested
    integer,      intent(in)  :: derv
    !> number of geometry coordinates
    integer,      intent(in)  :: ncor
    !> one external potential label, or 'NONE'
    character(4), intent(in)  :: field
    !> number of components of the field (must be all, ie. 1/3/6)
    integer,      intent(in)  :: ncomp
    !> output tensor
    real(8),      intent(out) :: nucpot(ncor**derv * ncomp)
    !-------------------------------------------------------
    integer i, na
    na = get_nr_atoms()
    if (ncor /= 3*na) &
       call quit('rsp_backend nuclear_potential error: expected ncor = 3*natom')
    if (.not.associated(nuc_charges) .and. &
        .not.associated(nuc_coords)) call rsp_backend_setup
    if ((field == 'NONE' .and. ncomp /= 1) .or. &
        (field == 'EL  ' .and. ncomp /= 3)) &
       call quit('rsp_backend nuclear_potential error: expected ncomp = 1 or 3')
    ! switch
    if (derv == 0 .and. field == 'NONE') then
       call quit('rsp_ error: unperturbed nuc.rep. not implemented')
    else if (derv == 0 .and. field == 'EL  ') then
       call DIPNUC_ifc(nucpot)
       nucpot = -nucpot !dip = -dE/dF
    else if (derv == 1 .and. field == 'NONE') then
       call GRADNN_ifc(na, nucpot)
    else if (derv == 2 .and. field == 'NONE') then
       call HESSNN_ifc(na, nucpot)
    else if (derv == 3 .and. field == 'NONE') then
       call cubicff_nuc(na, nuc_charges, nuc_coords, nucpot)
    else if (derv == 4 .and. field == 'NONE') then
       call quarticff_nuc(na, nuc_charges, nuc_coords, nucpot)
    else
       print *,'rsp_backend nuclear_potential error: not implented, derv =', &
               derv, '  field = ', field
       call quit('rsp_backend nuclear_potential error: not implented')
    end if
  end subroutine


  !> Nuclear repulsion contribution to cubic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  subroutine cubicff_nuc(na, chg, cor, cub)
     integer, intent(in)  :: na
     real(8), intent(in)  :: chg(na)
     real(8), intent(in)  :: cor(3,na)
     real(8), intent(out) :: cub(3,na,3,na,3,na)
     real(8) r(3), rrr(3,3,3), trace(3)
     integer i, j, k, l
     cub = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct triple tensor prod with traces removed
           r = cor(:,j) - cor(:,i)
           rrr = reshape((/((r(:)*r(k)*r(l),k=1,3),l=1,3)/),(/3,3,3/))
           trace = rrr(:,1,1) + rrr(:,2,2) + rrr(:,3,3)
           do k = 1, 3
              rrr(:,k,k) = rrr(:,k,k) - trace/5
              rrr(k,:,k) = rrr(k,:,k) - trace/5
              rrr(k,k,:) = rrr(k,k,:) - trace/5
           end do
           ! apply scale factor 15 Qi Qj / r^7
           rrr = rrr * (15 * chg(i) * chg(j) / sum(r**2)**(7/2d0))
           cub(:,i,:,i,:,i) = cub(:,i,:,i,:,i) + rrr !iii
           cub(:,i,:,i,:,j) = -rrr !iij
           cub(:,i,:,j,:,i) = -rrr !iji
           cub(:,j,:,i,:,i) = -rrr !jii
           cub(:,i,:,j,:,j) =  rrr !ijj
           cub(:,j,:,i,:,j) =  rrr !jij
           cub(:,j,:,j,:,i) =  rrr !jji
           cub(:,j,:,j,:,j) = cub(:,j,:,j,:,j) - rrr !jjj
        end do
     end do
  end subroutine


  !> Nuclear repulsion contribution to quartic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  !> Note: Dimensions 7 and 8 of 'qua' have been joined, to comply with
  !> fortran standard (max 7 dims)
  subroutine quarticff_nuc(na, chg, cor, qua)
     integer, intent(in)  :: na
     real(8), intent(in)  :: chg(na)
     real(8), intent(in)  :: cor(3,na)
     real(8), intent(out) :: qua(3,na,3,na,3,na,3*na)
     real(8) r(3), rrrr(3,3,3,3), trace(3,3), trace2
     integer i, j, k, l, m
     qua = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct quadruple tensor product with traces removed
           r = cor(:,j) - cor(:,i)
           rrrr = reshape((/(((r(:)*r(k)*r(l)*r(m), &
                          k=1,3),l=1,3),m=1,3)/),(/3,3,3,3/))
           trace  = rrrr(:,:,1,1) + rrrr(:,:,2,2) + rrrr(:,:,3,3)
           trace2 = trace(1,1) + trace(2,2) + trace(3,3)
           do k = 1, 3
              trace(k,k) = trace(k,k) - trace2/10
           end do
           do k = 1, 3
              rrrr(:,:,k,k) = rrrr(:,:,k,k) - trace/7
              rrrr(:,k,:,k) = rrrr(:,k,:,k) - trace/7
              rrrr(k,:,:,k) = rrrr(k,:,:,k) - trace/7
              rrrr(:,k,k,:) = rrrr(:,k,k,:) - trace/7
              rrrr(k,:,k,:) = rrrr(k,:,k,:) - trace/7
              rrrr(k,k,:,:) = rrrr(k,k,:,:) - trace/7
           end do
           ! apply scale factor 105 Qi Qj / r^9
           rrrr = rrrr * (105 * chg(i) * chg(j) / sum(r**2)**(9/2d0))
           qua(:,i,:,i,:,i,3*i-2:3*i) = qua(:,i,:,i,:,i,3*i-2:3*i) + rrrr !iiii
           qua(:,i,:,i,:,i,3*j-2:3*j) = -rrrr !iiij
           qua(:,i,:,i,:,j,3*i-2:3*i) = -rrrr !iiji
           qua(:,i,:,j,:,i,3*i-2:3*i) = -rrrr !ijii
           qua(:,j,:,i,:,i,3*i-2:3*i) = -rrrr !jiii
           qua(:,i,:,i,:,j,3*j-2:3*j) =  rrrr !iijj
           qua(:,i,:,j,:,i,3*j-2:3*j) =  rrrr !ijij
           qua(:,j,:,i,:,i,3*j-2:3*j) =  rrrr !jiij
           qua(:,i,:,j,:,j,3*i-2:3*i) =  rrrr !ijji
           qua(:,j,:,i,:,j,3*i-2:3*i) =  rrrr !jiji
           qua(:,j,:,j,:,i,3*i-2:3*i) =  rrrr !jjii
           qua(:,i,:,j,:,j,3*j-2:3*j) = -rrrr !ijjj
           qua(:,j,:,i,:,j,3*j-2:3*j) = -rrrr !jijj
           qua(:,j,:,j,:,i,3*j-2:3*j) = -rrrr !jjij
           qua(:,j,:,j,:,j,3*i-2:3*i) = -rrrr !jjji
           qua(:,j,:,j,:,j,3*j-2:3*j) = qua(:,j,:,j,:,j,3*j-2:3*j) + rrrr !jjjj
        end do
     end do
  end subroutine


  !> fetch whatever needed from COMMON blocks, etc.
  subroutine rsp_backend_setup()
    implicit integer (i,m-n)
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
    na = get_nr_atoms()
    allocate(nuc_charges(na))
    nuc_charges = CHARGE(:na)
    allocate(nuc_coords(3,na))
    nuc_coords = CORD(:,:na)
  end subroutine


  subroutine rsp_backend_finalize()
    deallocate(nuc_charges)
    deallocate(nuc_coords)
  end subroutine

  function is_ks_calculation()
     logical :: is_ks_calculation
#include "maxorb.h"
#include "priunit.h"
#include "mxcent.h"
#include "infinp.h"
     is_ks_calculation = dodft
  end function
  
end module
