  module rsp_mag_prop

  contains

  function D_efishg(setup_i, Effff)

    implicit none

    integer :: setup_i
    complex(8), dimension(3,3,3,3) :: Effff
    real(8), dimension(3,3) :: A
    complex(8) :: D_efishg

    ! Define A
    A(1,1) = 0.5
    A(1,2) = 2.0
    A(1,3) = 1.0
    A(2,1) = 0.5
    A(2,2) = 1.0
    A(2,3) = -2.0
    A(3,1) = 0.25
    A(3,2) = 3.0
    A(3,3) = -1.0

    D_efishg = ( A(setup_i,1)  /  225.0 ) * R_efishg(setup_i, A, Effff)**(2.0)

  end function


  function R_efishg(setup_i, A, Effff)

    implicit none

    integer :: setup_i, lambda, mu
    complex(8), dimension(3,3,3,3) :: Effff
    real(8), dimension(3,3) :: A
    complex(8) :: R_efishg

    R_efishg = 0.0

    do lambda = 1, 3
       do mu = 1, 3

          R_efishg = R_efishg + A(setup_i, 2) * Effff(lambda, lambda, mu, mu) + &
                                A(setup_i, 3) * Effff(lambda, mu, mu, lambda) 

       end do
    end do

  end function


  function M_efishg(setup_i, Effff, Effmfww, Effmfw2w)

    implicit none

    integer :: setup_i, nu, tau
    complex(8), dimension(3,3,3,3) :: Effff, Effmfww, Effmfw2w
    real(8), dimension(3,3) :: A
    real(8), dimension(3,6) :: B
    complex(8) :: M_efishg

    ! Define A
    A(1,1) = 0.5
    A(1,2) = 2.0
    A(1,3) = 1.0
    A(2,1) = 0.5
    A(2,2) = 1.0
    A(2,3) = -2.0
    A(3,1) = 0.25
    A(3,2) = 3.0
    A(3,3) = -1.0

    ! Define B

    B(1,1) = -1.0
    B(1,2) = 1.0
    B(1,3) = -4.0
    B(1,4) = 1.0
    B(1,5) = -2.0
    B(1,6) = -1.0
    B(2,1) = 1.0
    B(2,2) = 4.0
    B(2,3) = -1.0
    B(2,4) = -1.0
    B(2,5) =  1.0
    B(2,6) = -2.0
    B(3,1) = 0.5
    B(3,2) = 3.0
    B(3,3) = 3.0
    B(3,4) = -2.0
    B(3,5) = 3.0
    B(3,6) = -1.0




    M_efishg = 0.0

    do nu = 1, 3
       do tau = 1, 3

          M_efishg = M_efishg + ( B(setup_i,1)  /  225.0 ) * &
                     R_efishg(setup_i, A, Effff) * &
                   ( B(setup_i,2) * Effmfww(nu, nu, tau, tau) + &
                     B(setup_i,3) * Effmfww(nu, tau, nu, tau) + &
                     B(setup_i,4) * Effmfww(nu, tau, tau, nu) + &
                     B(setup_i,5) * Effmfw2w(tau, nu, nu, tau) + &
                     B(setup_i,6) * Effmfw2w(tau, tau, nu, nu) )

       end do
    end do

  end function


  function L_efishg(setup_i, Effff, Effqfww, Effqfw2w)

    implicit none

    integer :: setup_i, nu, tau, pi, rho, sigma, pisigma_offset, nusigma_offset
    complex(8), dimension(3,3,3,3) :: Effff
    complex(8), dimension(3,3,6,3) :: Effqfww, Effqfw2w
    real(8), dimension(3,3) :: A
    real(8), dimension(3,5) :: C
    complex(8) :: L_efishg


    ! Define A
    A(1,1) = 0.5
    A(1,2) = 2.0
    A(1,3) = 1.0
    A(2,1) = 0.5
    A(2,2) = 1.0
    A(2,3) = -2.0
    A(3,1) = 0.25
    A(3,2) = 3.0
    A(3,3) = -1.0

    ! Define C
    C(1,1) = -1.0
    C(1,2) = 1.0
    C(1,3) = 0.0
    C(1,4) = -1.0
    C(1,5) = 0.0
    C(2,1) = -1.0
    C(2,2) = 1.0
    C(2,3) = -1.0
    C(2,4) = 2.0
    C(2,5) = 2.0
    C(3,1) = -0.5
    C(3,2) = 2.0
    C(3,3) = -1.0
    C(3,4) = 1.0
    C(3,5) = 2.0




    L_efishg = 0.0

    do nu = 1, 3
       do tau = 1, 3
          do pi = 1, 3
             do rho = 1, 3
                do sigma = 1, 3

                   pisigma_offset = elgr_offset(pi, sigma)
                   nusigma_offset = elgr_offset(nu, sigma)

                   L_efishg = L_efishg + ( C(setup_i,1)  /  225.0 ) * &
                                           R_efishg(setup_i, A, Effff ) * &
                              ( W_efishg_1(setup_i, C, nu, tau, pi, rho, sigma) * &
                                Effqfww(nu, tau, pisigma_offset, rho) + &
                                W_efishg_2(setup_i, C, nu, tau, pi, rho, sigma) * &
                                Effqfw2w(pi, tau, nusigma_offset, rho) )

                end do
             end do
          end do
       end do
    end do

  end function


  function W_efishg_1(setup_i, C, a, b, g, d, e)

    implicit none

    integer :: setup_i, a, b, g, d, e
    real(8), dimension(3,5) :: C
    real(8) :: W_efishg_1


    W_efishg_1 = C(setup_i, 2) * eps(a, b, g) * delta(d, e) + &
                 C(setup_i, 3) * eps(a, b, d) * delta(g, e) + &
                 C(setup_i, 4) * eps(a, g, d) * delta(b, e)

  end function


  function W_efishg_2(setup_i, C, a, b, g, d, e)

    implicit none

    integer :: setup_i, a, b, g, d, e
    real(8), dimension(3,5) :: C
    real(8) :: W_efishg_2

    W_efishg_2 = C(setup_i, 5) * eps(a, b, d) * delta(g, e)

  end function


  function delta(a, b)

    implicit none

    integer :: a, b, delta

    if (a == b) then

       delta = 1.0

    else

       delta = 0.0

    end if

  end function


  ! Levi-Civita
  function eps(a, b, c)

    implicit none

    integer :: a, b, c, eps

    eps = (a - b)*(b - c)*(c - a)/2

  end function


  function elgr_offset(a, b)

    implicit none

    integer :: a, b, elgr_offset

    elgr_offset = 0

    if ( (a == 1) .and. (b == 1) ) then
       elgr_offset = 1

    else if ( (a == 1) .and. (b == 2) ) then
       elgr_offset = 2
    
    else if ( (a == 1) .and. (b == 3) ) then
       elgr_offset = 3

    else if ( (a == 2) .and. (b == 1) ) then
       elgr_offset = 2

    else if ( (a == 2) .and. (b == 2) ) then
       elgr_offset = 4

    else if ( (a == 2) .and. (b == 3) ) then
       elgr_offset = 5

    else if ( (a == 3) .and. (b == 1) ) then
       elgr_offset = 3

    else if ( (a == 3) .and. (b == 2) ) then
       elgr_offset = 5

    else if ( (a == 3) .and. (b == 3) ) then
       elgr_offset = 6

    end if

  end function

end module rsp_mag_prop
