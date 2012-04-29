module interface_molecule

!  interface_molecule_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_molecule_init

   implicit none

   public interface_molecule_init
   public interface_molecule_finalize

   public get_nr_ao
   public get_nr_atoms
   public get_nuc_name
   public get_nuc_charge
   public get_nuc_xyz
   public get_is_ks_calculation
   public get_dipole_origin

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   integer :: nr_ao
   integer :: nr_atoms
   logical :: is_ks_calculation !fixme move to interface_xc
   real(8) :: dipole_origin(3)

!  allocatables, deallocate them in interface_molecule_finalize
   character(4), allocatable :: nuc_name(:)
   real(8),      allocatable :: nuc_charge(:)
   real(8),      allocatable :: nuc_xyz(:, :)

contains

   subroutine interface_molecule_init()

#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "orgcom.h"

      nr_ao      = nbast
      nr_atoms   = natoms
      is_ks_calculation = dodft

      allocate(nuc_name(nr_atoms))
      nuc_name = namn(:nr_atoms)

      allocate(nuc_charge(nr_atoms))
      nuc_charge = charge(:nr_atoms)

      allocate(nuc_xyz(3, nr_atoms))
      nuc_xyz = cord(:, :nr_atoms)

      dipole_origin = diporg

      is_initialized = .true.

   end subroutine

   subroutine interface_molecule_finalize()

!     here deallocate everything that is allocated
      if (allocated(nuc_name))   deallocate(nuc_name)
      if (allocated(nuc_charge)) deallocate(nuc_charge)
      if (allocated(nuc_xyz))    deallocate(nuc_xyz)

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_molecule'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   integer function get_nr_ao()
      call check_if_interface_is_initialized()
      get_nr_ao = nr_ao
   end function

   integer function get_nr_atoms()
      call check_if_interface_is_initialized()
      get_nr_atoms = nr_atoms
   end function

   character(4) function get_nuc_name(i)
      integer, intent(in) :: i
      call check_if_interface_is_initialized()
      get_nuc_name = nuc_name(i)
   end function

   real(8) function get_nuc_charge(i)
      integer, intent(in) :: i
      call check_if_interface_is_initialized()
      get_nuc_charge = nuc_charge(i)
   end function

   real(8) function get_nuc_xyz(i, j)
      integer, intent(in) :: i
      integer, intent(in) :: j
      call check_if_interface_is_initialized()
      get_nuc_xyz = nuc_xyz(i, j)
   end function

   logical function get_is_ks_calculation()
      call check_if_interface_is_initialized()
      get_is_ks_calculation = is_ks_calculation
   end function

   function get_dipole_origin()
      real(8) :: get_dipole_origin(3)
      call check_if_interface_is_initialized()
      get_dipole_origin = dipole_origin
   end function

end module
