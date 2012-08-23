module interface_molecule

!  interface_molecule_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_molecule_init

   implicit none

   public interface_molecule_init
   public interface_molecule_finalize

   public get_nr_atoms
   public get_nuc_name
   public get_nuc_charge
   public get_nuc_isotope
   public get_nuc_xyz
   public get_dipole_origin

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   integer :: nr_atoms
   real(8) :: dipole_origin(3)

!  allocatables, deallocate them in interface_molecule_finalize
   character(4), allocatable :: nuc_name(:)
   real(8),      allocatable :: nuc_charge(:)
   integer,      allocatable :: nuc_isotope(:)
   real(8),      allocatable :: nuc_xyz(:, :)

contains

#ifdef VAR_LSDALTON
  !A proper interface that is independent of the host program. TK
   subroutine interface_molecule_init(natoms,namn,charge,cord,diporg,isotop)
     implicit none
     integer      :: natoms
     character(4) :: namn(natoms)
     real(8)      :: charge(natoms),cord(3,natoms),diporg(3)
     integer      :: isotop(natoms)
     nr_atoms = natoms
     allocate(nuc_name(nr_atoms))
     nuc_name = namn(:nr_atoms)
     allocate(nuc_charge(nr_atoms))
     nuc_charge = charge(:nr_atoms)
     allocate(nuc_isotope(nr_atoms))
     nuc_isotope = isotop(:nr_atoms)
     allocate(nuc_xyz(3, nr_atoms))
     nuc_xyz = cord(:, :nr_atoms)
     dipole_origin = diporg
     is_initialized = .true.
   end subroutine
#else
   subroutine interface_molecule_init()

#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "orgcom.h"

      nr_atoms = natoms

      allocate(nuc_name(nr_atoms))
      nuc_name = namn(:nr_atoms)

      allocate(nuc_charge(nr_atoms))
      nuc_charge = charge(:nr_atoms)

      allocate(nuc_isotope(nr_atoms))
#ifdef PRG_DIRAC
      ! isotop array is not available in DIRAC
      ! this should be program independent
      ! for the moment use workaround
      nuc_isotope = 1
#else
      nuc_isotope = isotop(:nr_atoms)
#endif

      allocate(nuc_xyz(3, nr_atoms))
      nuc_xyz = cord(:, :nr_atoms)

      dipole_origin = diporg

      is_initialized = .true.

   end subroutine
#endif

   subroutine interface_molecule_finalize()

!     here deallocate everything that is allocated
      if (allocated(nuc_name))    deallocate(nuc_name)
      if (allocated(nuc_charge))  deallocate(nuc_charge)
      if (allocated(nuc_isotope)) deallocate(nuc_isotope)
      if (allocated(nuc_xyz))     deallocate(nuc_xyz)

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_molecule'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

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

   integer function get_nuc_isotope(i)
      integer, intent(in) :: i
      call check_if_interface_is_initialized()
      get_nuc_isotope = nuc_isotope(i)
   end function

   real(8) function get_nuc_xyz(i, j)
      integer, intent(in) :: i
      integer, intent(in) :: j
      call check_if_interface_is_initialized()
      get_nuc_xyz = nuc_xyz(i, j)
   end function

   function get_dipole_origin()
      real(8) :: get_dipole_origin(3)
      call check_if_interface_is_initialized()
      get_dipole_origin = dipole_origin
   end function

end module
