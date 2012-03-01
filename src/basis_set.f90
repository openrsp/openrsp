! Copyright 2012-     ?
!           2009-2011 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3

!> @file
!> Contains module basis_set

!> This module contains type(s) describing sets of basis functions
!> as well as a routine calculating coefficients for the transformation
!> from Cartesian to spherical (solid-harmonic) basis
module basis_set

  implicit none
  public cgto
  public cart_to_spher_coef

  !> Derived type describing a shell-tuple of contracted
  !> Gaussian-type orbitals (basis functions)
  type cgto
    !sequence
    !> coordinates of center
    real(8) :: cent(3)
    !> anugular momentum
    integer :: mom
    !> exponents in exp(-exp(i)|r-cent|^2)
    real(8), pointer :: exp(:)
    !> contraction coefficients ctr(nrad,nexp)
    real(8), pointer :: ctr(:,:)
    !-----  additional info, for convenience  -----
    !> charge of nucleus
    integer :: charge
    !> index in molecule(:) of cent
    integer :: icent
    !> number of AOs in shell == nrad * (2*mom+1)
    integer :: nbas
    !> offset of shell in AO basis enumeration
    integer :: ibas
  end type

  private

contains

  !> Cartesian to spherical (real solid-harmonic) contraction
  !> coefficients. Spherical indices ordered im(Yl)...Y0...re(Yl)
  function cart_to_spher_coef(amom) result(coef)
    !> momentum of Cartesians to contract
    !> and angular momentum of the resulting sphericals
    integer, intent(in) :: amom !angular momentum
    !> resulting contraction coefficients (mostly zeros)
    real(8)             :: coef(2*amom+1, (amom+1)*(amom+2)/2)
    !------------------------------------
    real(8), parameter :: pi = acos(-1d0)
    ! Legendre coefficient, normalization constant
    real(8) :: legendre, invnorm
    ! binomial coefficients for (x+iy)^i (xx+yy)^j
    integer :: binomial(amom+1), binomixy(amom+1)
    integer :: i, j, k, m, mxy, my
    ! zero all, before filling in non-zeros
    coef(:,:) = 0
    ! initialize to normalization factor for m=0
    invnorm = sqrt((2*amom+1d0)/pi)/2
    ! initialize binomial coefficients
    binomial(1)  = 1
    binomial(2:) = 0
    ! loop over spherical indices (+ and - together)
    do m = 0, amom
       ! calculate next normalization factor (except m==0)
       if (m/=0) invnorm = invnorm &
                         * sqrt((amom+1d0-m)*(amom+m)) / (2*m)
       ! increment binomial coefficients
       binomial(2:m+1) = binomial(2:m+1) + binomial(1:m)
       ! ++--copy binomial coefficients into binomixy (for recursion)
       binomixy(:m+1) = (/((-1)**(i/2),i=0,m)/) &
                      * binomial(:m+1)
       binomixy(m+2:) = 0
       ! first Legendre coefficient is always 1
       legendre = 1
       ! loop over contributing rows (with mx+my=mxy)
       do mxy = m, amom, 2
          ! offset of row
          i = (amom+1)*(amom+2)/2 - (mxy+1)*(mxy+2)/2
          ! increment Legendre coefficient
          if (mxy/=m) legendre = -(legendre &
                               * (amom-mxy+1) &
                               * (amom-mxy+2)) &
                               / ((mxy-m) * (mxy+m))
          ! increment +-binomial coefficients
          if (mxy/=m) binomixy(3:mxy+1) = binomixy(3:mxy+1) &
                                        + binomixy(1:mxy-1)
          ! coefficients for real part Ym
          do my=0, mxy, 2
             coef(amom+1+m,i+my+1) = invnorm * legendre &
                                   * binomixy(my+1)
          end do
          ! coefficients for imaginary part Y-m
          if (m==0) cycle !they are zero m==0
          do my=1, mxy, 2
             coef(amom+1-m,i+my+1) = invnorm * legendre &
                                   * binomixy(my+1)
          end do
       end do
       ! for m other than m=0 there are both real and imag parts,
       ! and the normalization constants are equal and sqrt(2)
       ! times the normalization constant for the complex
       if (m==0) invnorm = invnorm * sqrt(2d0)
    end do
    ! if p shell, reorder yzx as xyz
    if (amom==1) then
       coef(1,1) = coef(3,1);  coef(3,1)=0
       coef(2,2) = coef(1,2);  coef(1,2)=0
       coef(3,3) = coef(2,3);  coef(2,3)=0
    end if
  end function

end module
