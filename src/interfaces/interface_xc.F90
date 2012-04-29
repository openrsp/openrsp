module interface_xc

!  interface_xc_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_xc_init

   use matrix_defop

   implicit none

   public interface_xc_init
   public interface_xc_finalize

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT
   public di_get_geomDeriv_FxD_DFT
   public di_get_geomDeriv_GxD_DFT
   public di_get_geomDeriv_molgrad_DFT

   public get_is_ks_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_ks_calculation

contains

   subroutine interface_xc_init()

#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"

      is_ks_calculation = dodft

      is_initialized = .true.

   end subroutine

   subroutine interface_xc_finalize()

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_xc'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   logical function get_is_ks_calculation()
      call check_if_interface_is_initialized()
      get_is_ks_calculation = is_ks_calculation
   end function

   !> \brief calculates the magnetic derivative of the F^xc matrix defined in Equation XX of
   !> \author Bin Gao
   !> \date 2009-12-10
   !> \param D
   !> \return Fx contains the magnetic derivative
   subroutine di_get_MagDeriv_FxD_DFT( Fx, D )
     type(matrix), intent(in) :: D
     type(matrix), intent(inout) :: Fx(3)
     call QUIT( 'di_get_MagDeriv_FxD_DFT is not implemented!' )
   end subroutine


   !> \brief calculates the magnetic derivative of the G^xc matrix defined in Equation XX of
   !> \author Bin Gao
   !> \date 2009-12-10
   !> \param D
   !> \param A
   !> \return Gx contains the magnetic derivative
   subroutine di_get_MagDeriv_GxD_DFT( D, A, Gx )
     type(matrix), intent(in) :: D
     type(matrix), intent(in) :: A
     type(matrix), intent(inout) :: Gx(3)
     call QUIT( 'di_get_MagDeriv_GxD_DFT is not implemented!' )
   end subroutine


   !> \brief calculates the exchange correlation part of the geometric derivative of the F^ks matrix
   !> \author Bin Gao
   !> \date 2009-12-10
   !> \param natoms is the number of atoms
   !> \param D
   !> \param B
   !> \return gradient contains the exchange correlation part of the geometric derivative of the F^ks matrix
   subroutine di_get_geomDeriv_FxD_DFT( gradient, natoms, D, B )
     type(matrix), intent(in) :: D
     type(matrix), intent(in) :: B
     real(8), intent(out) :: gradient( 3*natoms )
     integer, intent(in) :: natoms
     call QUIT( 'di_get_geomDeriv_FxD_DFT is not implemented!' )
   end subroutine


   !> \brief calculates the exchange correlation part of the geometric derivative of the F^ks matrix
   !> \author Bin Gao
   !> \date 2009-12-10
   !> \param natoms is the number of atoms
   !> \param D
   !> \return gradient contains the exchange correlation part of the geometric derivative of the F^ks matrix
   subroutine di_get_geomDeriv_GxD_DFT( gradient, natoms, D, A, B )
     type(matrix), intent(in) :: D
     type(matrix), intent(in) :: A
     type(matrix), intent(in) :: B
     real(8), intent(out) :: gradient( 3*natoms )
     integer, intent(in) :: natoms
     call QUIT( 'di_get_geomDeriv_GxD_DFT is not implemented!' )
   end subroutine


   !> \brief calculates the exchange correlation part of the geometric derivative of the F^ks matrix
   !> \author Bin Gao
   !> \date 2009-12-10
   !> \param natoms is the number of atoms
   !> \param D
   !> \return gradient contains the exchange correlation part of the geometric derivative of the F^ks matrix
   subroutine di_get_geomDeriv_molgrad_DFT( gradient, natoms, D )
     type(matrix), intent(in) :: D
     real(8), intent(out) :: gradient( 3*natoms )
     integer, intent(in) :: natoms
     call QUIT( 'di_get_geomDeriv_molgrad_DFT is not implemented!' )
   end subroutine

end module
