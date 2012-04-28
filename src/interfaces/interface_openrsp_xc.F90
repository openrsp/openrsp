module interface_openrsp_xc

   use matrix_defop

   implicit none

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT
   public di_get_geomDeriv_FxD_DFT
   public di_get_geomDeriv_GxD_DFT
   public di_get_geomDeriv_molgrad_DFT

   private

contains

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
