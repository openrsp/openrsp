
module vib_pv_contribs

contains

!-------------------------------------------------------------------------------
!
! Contains functions for transforming derivatives from Cartesian to normal
! coordinates. The various functions do this task for various combinations of
! "ways" in the properties being transformed (i.e. dipole moments being a 
! one-way term and polarizability tensors being two-way terms), and the degree 
! of differentiation of the quantity under consideration.
!
! The functions are named with this in mind - for instance, trans_cartnc_1w1d
! transforms the gradient of a one-way term (such as the molecular dipole
! moment tensor), while trans_cartnc_3w2d transforms the Hessian of a 3-way
! term (such as the molecular 1st hyperpolarizability tensor).
!
!-------------------------------------------------------------------------------

  function trans_cartnc_0w1d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, pp
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8)						:: tmp ! Temporary storage
    complex(8), 	dimension(n_c)				:: A ! The tensor being transformed
    complex(8), 	dimension(n_nm)				:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

! Loop over normal coordinates

    do p = 1, n_nm

! Loop over Cartesian coordinates

      do pp = 1, n_c
	
! Add to a temporary variable element for element the Cartesian element projected on the normal coordinate

	tmp = tmp + A(pp) * T(pp, p)
	
      end do

! Put in result tensor

      B(p) = tmp

    end do      
      
  end function

! The succeeding functions follow the same procedure, but also taking into account increased tensor ranks
! due to increased levels of differentiation and/or increased ranks of the tensors being differentiated.

  function trans_cartnc_0w2d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, q, pp, qp
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8)						:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm

	do qp = 1, n_c

	  tmp = tmp + A(p, qp) * T(qp, q)

	end do

	Bt(p, q) = tmp

      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm

	do pp = 1, n_c

	  tmp = tmp + Bt(pp, q) * T(pp, q)

	end do

	B(p, q) = tmp

      end do
    end do       
      
  end function

  function trans_cartnc_0w3d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, r, pp, qp, rp
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8)							:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, n_c)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, n_nm)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm

	  do rp = 1, n_c

	    tmp = tmp + A(p, q, rp) * T(rp, r)

	  end do

	  Bt(p, q, r) = tmp

	end do
      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do qp = 1, n_c

	    tmp = tmp + Bt(p, qp, r) * T(qp, q)

	  end do

	  A(p, q, r) = tmp

	end do
      end do
    end do


    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm

	  do pp = 1, n_c

	    tmp = tmp + A(pp, q, r) * T(pp, p)

	  end do

	  B(p, q, r) = tmp

	end do
      end do
    end do


  end function

  function trans_cartnc_0w4d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, r, s, pp, qp, rp, sp
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8)							:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, n_c, n_c)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, n_nm, n_nm)		:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  do s = 1, n_nm

	    do sp = 1, n_c

	    tmp = tmp + A(p, q, r, sp) * T(sp, s)

	    end do

	    Bt(p, q, r, s) = tmp

	    end do
	  end do
	end do
      end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  do s = 1, n_nm

	    do rp = 1, n_c

	      tmp = tmp + A(p, q, rp, s) * T(rp, r)

	    end do

	    Bt(p, q, r, s) = tmp

	  end do
	end do
      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  do s = 1, n_nm
	  
	    do qp = 1, n_c

	      tmp = tmp + Bt(p, qp, r, s) * T(qp, q)

	    end do

	    A(p, q, r, s) = tmp

	  end do
	end do
      end do
    end do


    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  do s = 1, n_nm

	    do pp = 1, n_c

	      tmp = tmp + A(pp, q, r, s) * T(pp, p)

	    end do

	    B(p, q, r, s) = tmp

	  end do
	end do
      end do
    end do


  end function



  function trans_cartnc_1w1d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, pp, i
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, 3)			:: A ! The tensor being transformed
    complex(8), 	dimension(n_nm, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm

      do i = 1, 3
	tmp(i) = 0.0
      end do

      do pp = 1, n_c
	do i = 1, 3

	  tmp(i) = tmp(i) + A(pp, i) * T(pp, p)

	end do
      end do

      B(p, :) = tmp(:)

    end do      
      
  end function

  function trans_cartnc_2w1d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, pp, i, j
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, 3, 3)			:: A ! The tensor being transformed
    complex(8), 	dimension(n_nm, 3, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm

      do i = 1, 3
	do j = 1, 3
	  tmp(i, j) = 0.0
	end do
      end do


      do pp = 1, n_c
	do i = 1, 3
	  do j = 1, 3

	    tmp(i, j) = tmp(i, j) + A(pp, i, j) * T(pp, p)

	  end do
	end do
      end do

      B(p, :, :) = tmp(:, :)

    end do      
      
  end function

  function trans_cartnc_3w1d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, pp, i, j, k
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, 3, 3, 3)				:: A ! The tensor being transformed
    complex(8), 	dimension(n_nm, 3, 3, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 		:: T ! The matrix of transformation coefficients

    do p = 1, n_nm

      do i = 1, 3
	do j = 1, 3
	  do k = 1, 3
	    tmp(i, j, k) = 0.0
	  end do
	end do
      end do


      do pp = 1, n_c
	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3

	      tmp(i, j, k) = tmp(i, j, k ) + A(pp, i, j, k) * T(pp, p)

	    end do
	  end do
	end do
      end do

      B(p, :, :, :) = tmp(:, :, :)

    end do      
      
  end function


  function trans_cartnc_1w2d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, q, pp, qp, i
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, 3)		:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do i = 1, 3
	  tmp(i) = 0.0
	end do

	do qp = 1, n_c
	  do i = 1, 3

	    tmp(i) = tmp(i) + A(p, qp, i) * T(qp, q)

	  end do
	end do

	Bt(p, q, :) = tmp(:)

      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  tmp(i) = 0.0
	end do

	do pp = 1, n_c
	  do i = 1, 3

	    tmp(i) = tmp(i) + Bt(pp, q, i) * T(pp, q)

	  end do
	end do

	B(p, q, :) = tmp(:)

      end do
    end do       
      
  end function

  function trans_cartnc_2w2d(n_c, n_nm, A, T) result(B)

    implicit none

    integer						:: p, q, pp, qp, i, j
    integer, 	intent(in) 				:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, 3, 3)		:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, 3, 3)		:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3	
	    tmp(i, j) = 0.0	
	  end do
	end do

	do qp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3

	      tmp(i, j) = tmp(i, j) + A(p, qp, i, j) * T(qp, q)

	    end do
	  end do
	end do

	Bt(p, q, :, :) = tmp(:, :)

      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3	
	    tmp(i, j) = 0.0	
	  end do
	end do

	do pp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3

	      tmp(i, j) = tmp(i, j) + Bt(pp, q, i, j) * T(pp, q)

	    end do
	  end do
	end do

	B(p, q, :, :) = tmp(:, :)

      end do
    end do       
      
  end function


  function trans_cartnc_3w2d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, pp, qp, i, j, k
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, 3, 3, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, 3, 3, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      tmp(i, j, k) = 0.0
	    end do
	  end do
	end do

	do qp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
    
		tmp(i, j, k) = tmp(i, j, k) + A(p, qp, i, j, k) * T(qp, q)

	      end do
	    end do
	  end do
	end do

	Bt(p, q, :, : , :) = tmp(:, :, :)

      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      tmp(i, j, k) = 0.0
	    end do
	  end do
	end do

	do pp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		tmp(i, j, k) = tmp(i, j, k) + Bt(pp, q, i, j, k) * T(pp, q)

	      end do
	    end do
	  end do
	end do

	B(p, q, :, :, :) = tmp(:, :, :)

      end do
    end do       
      
  end function



  function trans_cartnc_1w3d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, r, pp, qp, rp, i
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3)					:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, n_c, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, n_nm, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    tmp(i) = 0.0
	  end do
	  	

	  do rp = 1, n_c
	    do i = 1, 3

	      tmp(i) = tmp(i) + A(p, q, rp, i) * T(rp, r)

	    end do
	  end do

	  Bt(p, q, r, :) = tmp(:)

	end do
      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    tmp(i) = 0.0
	  end do
	  	

	  do qp = 1, n_c
	    do i = 1, 3

	      tmp(i) = tmp(i) + Bt(p, qp, r, i) * T(qp, q)

	    end do
	  end do

	  A(p, q, r, :) = tmp(:)

	end do
      end do
    end do


    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    tmp(i) = 0.0
	  end do
	  	

	  do pp = 1, n_c
	    do i = 1, 3

	      tmp(i) = tmp(i) + A(pp, q, r, i) * T(pp, p)

	    end do
	  end do

	  B(p, q, r, :) = tmp(:)

	end do
      end do
    end do


  end function



  function trans_cartnc_2w3d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, r, pp, qp, rp, i, j
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3)					:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, n_c, 3, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, n_nm, 3, 3)		:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      tmp(i, j) = 0.0
	    end do
	  end do
	  	

	  do rp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3

		tmp(i, j) = tmp(i, j) + A(p, q, rp, i, j) * T(rp, r)

	      end do
	    end do
	  end do

	  Bt(p, q, r, :, :) = tmp(:, :)

	end do
      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      tmp(i, j) = 0.0
	    end do
	  end do
	  	

	  do qp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3

		tmp(i, j) = tmp(i, j) + Bt(p, qp, r, i, j) * T(qp, q)

	      end do
	    end do
	  end do

	  A(p, q, r, :, :) = tmp(:, :)

	end do
      end do
    end do


    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      tmp(i, j) = 0.0
	    end do
	  end do
	  	

	  do pp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3

		tmp(i, j) = tmp(i, j) + A(pp, q, r, i, j) * T(pp, p)

	      end do
	    end do
	  end do

	  B(p, q, r, :, :) = tmp(:, :)

	end do
      end do
    end do


  end function


  function trans_cartnc_3w3d(n_c, n_nm, A, T) result(B)

    implicit none

    integer								:: p, q, r, pp, qp, rp, i, j, k
    integer, 	intent(in) 						:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3, 3)					:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, n_c, 3, 3, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, n_nm, 3, 3, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 		:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		tmp(i, j, k) = 0.0
	      end do
	    end do
	  end do
	  	

	  do rp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  tmp(i, j, k) = tmp(i, j, k) + A(p, q, rp, i, j, k) * T(rp, r)

		end do
	      end do
	    end do
	  end do

	  Bt(p, q, r, :, :, :) = tmp(:, :, :)

	end do
      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		tmp(i, j, k) = 0.0
	      end do
	    end do
	  end do
	  	

	  do qp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  tmp(i, j, k) = tmp(i, j, k) + Bt(p, qp, r, i, j, k) * T(qp, q)

		end do
	      end do
	    end do
	  end do

	  A(p, q, r, :, :, :) = tmp(:, :, :)

	end do
      end do
    end do


    do p = 1, n_nm
      do q = 1, n_nm
	do r = 1, n_nm
	  
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		tmp(i, j, k) = 0.0
	      end do
	    end do
	  end do
	  	

	  do pp = 1, n_c
	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		tmp(i, j, k) = tmp(i, j, k) + A(pp, q, r, i, j, k) * T(pp, p)

		end do
	      end do
	    end do
	  end do

	  B(p, q, r, :, :, :) = tmp(:, :, :)

	end do
      end do
    end do


  end function


  function trans_cartnc_4w1d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, pp, i, j, k, m
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3, 3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, 3, 3, 3, 3)			:: A ! The tensor being transformed
    complex(8), 	dimension(n_nm, 3, 3, 3, 3)			:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in) 		:: T ! The matrix of transformation coefficients

    do p = 1, n_nm

      do i = 1, 3
	do j = 1, 3
	  do k = 1, 3
	    do m = 1, 3
	      tmp(i, j, k, m) = 0.0
	    end do
	  end do
	end do
      end do


      do pp = 1, n_c
	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      do m = 1, 3

		tmp(i, j, k, m) = tmp(i, j, k, m) + A(pp, i, j, k, m) * T(pp, p)

	      end do
	    end do
	  end do
	end do
      end do

      B(p, :, :, :, :) = tmp(:, :, :, :)

    end do      
      
  end function


  function trans_cartnc_4w2d(n_c, n_nm, A, T) result(B)

    implicit none

    integer							:: p, q, pp, qp, i, j, k, m
    integer, 	intent(in) 					:: n_c, n_nm ! Number of Cartesian and normal coordinates, resp.
    complex(8),	dimension(3, 3, 3, 3)				:: tmp ! Temporary storage
    complex(8), 	dimension(n_c, n_c, 3, 3, 3, 3)			:: A, Bt ! The tensor being transformed and temporary storage
    complex(8), 	dimension(n_nm, n_nm, 3, 3, 3, 3)		:: B ! The resulting NC transformed tensor
    real(8), 	dimension(n_c, n_nm), intent(in)	 	:: T ! The matrix of transformation coefficients

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      do m = 1, 3
		tmp(i, j, k, m) = 0.0
	      end do
	    end do
	  end do
	end do

	do qp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do m = 1, 3
    
		  tmp(i, j, k, m) = tmp(i, j, k, m) + A(p, qp, i, j, k, m) * T(qp, q)

		end do
	      end do
	    end do
	  end do
	end do

	Bt(p, q, :, :, :, :) = tmp(:, :, :, :)

      end do
    end do

    do p = 1, n_nm
      do q = 1, n_nm

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      do m = 1, 3
		tmp(i, j, k, m) = 0.0
	      end do
	    end do
	  end do
	end do

	do pp = 1, n_c
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do m = 1, 3

		  tmp(i, j, k, m) = tmp(i, j, k, m) + Bt(pp, q, i, j, k, m) * T(pp, q)
	  
		end do
	      end do
	    end do
	  end do
	end do

	B(p, q, :, :, :, :) = tmp(:, :, :, :)

      end do
    end do       
      
  end function


!-------------------------------------------------------------------------------
!
! Contains a few small combinatoric functions (factorization and permutation
! list generation).
!
!-------------------------------------------------------------------------------



  recursive function fact(n) result(factorial)

  ! Takes a number and returns its factorial

    implicit none
    integer :: n, factorial

    if (n > 1) then
  
      factorial = n*fact(n-1)

    else

      factorial = 1

    end if

  end function  

  recursive function make_perm(n, factn, seed) result(perm_ind_list)

    ! Takes a list of indices and returns all permutations
    ! of the elements of the list.

    implicit none
    integer, intent(in) :: n, factn
    integer :: i, j
    integer, dimension(n), intent(in) :: seed
    integer, dimension(factn,n) :: perm_ind_list

    if (n > 2) then

      do j = 1, size(seed)

	perm_ind_list(fact(n - 1) * (j - 1) + 1: fact(n - 1) * (j), 1) = seed(j) 

	if (j == 1) then

	  perm_ind_list(fact(n - 1) * (j - 1) + 1: fact(n - 1) * (j), 2:size(seed)) = &
	  make_perm(n - 1, fact(n - 1), seed(2:size(seed)))

	elseif (j == size(seed)) then

	  perm_ind_list(fact(n - 1) * (j - 1) + 1: fact(n - 1) * (j), 2:size(seed)) = &
	  make_perm(n - 1, fact(n - 1), seed(1:size(seed) - 1))

	else

	  perm_ind_list(fact(n - 1) * (j - 1) + 1: fact(n - 1) * (j), 2:size(seed)) = &
	  make_perm(n - 1, fact(n - 1), ((/seed(1:j -1), seed(j+1:size(seed))/)))

	end if

      end do

    elseif(n == 2) then

      perm_ind_list(1, :) = (/seed(1), seed(2)/)
      perm_ind_list(2, :) = (/seed(2), seed(1)/)

    else

      write(*,*) 'Error in make_perm: n is too small'

    end if

  end function

  function alpha_pv(n_nm, nm_e, opt_e, dm_1d, dm_2d, dm_3d, fcc, fcq, rst_elec, rst_mech, rst_order)
!--------------------------------------------------------------------------------
! Calculates and returns the pure vibrational contribution to the molecular
! polarizability tensor from various derivatives taken with 
! respect to the normal coordinates of the molecule. Derivatives above 
! first order are taken as optional, should their implementation be 
! unavailable. 
!
! Functions referred to as "Bishop/Luis/Kirtman" terms reproduce terms presented
! in Journal of Chemical Physics, Vol. 108, No. 24 (1998), 10013-17.
! 
!--------------------------------------------------------------------------------

    implicit none

      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i,j
      integer, optional,					intent(in) 	:: rst_elec ! Restriction on order of electrical anharmonicity
									    ! If defined: Restricted to [rst_elec]th order, otherwise max order
      integer, optional,					intent(in) 	:: rst_mech ! Restriction on order of mechanical anharmonicity
									    ! If defined: Restricted to [rst_mech]th order, otherwist max order
      integer, optional,					intent(in) 	:: rst_order ! Restriction on total order of anharmonicity
									    ! If defined: Restricted to [rst_order]th order, otherwist max order
      integer, dimension(2, 2)							:: perm
      real(8), dimension(2)							:: opt_e ! Photon energies, negative sum first, then incident
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 2)						:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 2)						:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3), optional, 		intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3), optional,	intent(in) 	:: dm_3d ! 3rd                   "
      complex(8), dimension(n_nm, n_nm, n_nm), optional, 		intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm), optional,	intent(in) 	:: fcq   ! Quartic force constants
      complex(8), dimension(3, 3)							:: alpha_pv


    ! Precalculate terms for use in Bishop/Luis/Kirtman term functions

    nm_e_inv = invvec(n_nm, nm_e)

    lambda_a = make_lambda_a(2, n_nm, opt_e, nm_e)
    lambda_ab = make_lambda_ab(2, n_nm, opt_e, nm_e)

    do i = 1, 3
      do j = 1, 3

	alpha_pv(i, j) = 0.0
 
      end do
    end do
    
    ! Decide between various configurations of anharmonicity order restrictions, do calculation accordingly

    if (present(rst_order)) then

      if (rst_order == 1) then

      ! Take terms up to and including 1st total order of anharmonicity

	write(*,*) 'Warning: Due to vanishing terms, no precision is gained above 0th order.'

	do i = 1, 3
	  do j = 1, 3

	    ! Make permutation indices. While a trivial matter here,
	    ! the function performing this is called nonetheless, for
	    ! consistency with the corresponding task for beta and gamma.

	    perm = make_perm(2, fact(2), (/i, j/))

	    ! Add Bishop/Luis/Kirtman terms to the array to be returned

	    alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	  end do
	end do

      end if

      if (rst_order == 2) then

      ! Take terms up to and including 2nd total order of anharmonicity

	do i = 1, 3
	  do j = 1, 3

	    perm = make_perm(2, fact(2), (/i, j/))

	    alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
	    pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d) + &
	    pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc) + &
	    pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)

	  end do
	end do


      end if


    elseif (present(rst_elec) .OR. present(rst_mech)) then

      if (present(rst_elec) .AND. present(rst_mech)) then

	if (rst_elec == 0) then

	  if (rst_mech == 0) then

	  ! 0th order for electric, 0th order for mechanical

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	      end do
	    end do


	  end if



	  if (rst_mech == 1) then

	  ! 0th order for electric, 1st order for mechanical

	    write(*,*) 'Warning: Due to vanishing terms, no precision is gained above 0th order.'
	    
	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	      end do
	    end do

	  end if
	  

	  if (rst_mech == 2) then

	  ! 0th order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)

	      end do
	    end do

	  end if
	end if


	if (rst_elec == 1) then

	  if (rst_mech == 0) then

	  ! 1st order for electric, 0th order for mechanical

	    write(*,*) 'Warning: Due to vanishing terms, no precision is gained above 0th order.'

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	      end do
	    end do

	  end if


	  if (rst_mech == 1) then

	  ! 1st order for electric, 1st order for mechanical

	    do i = 1, 3
 
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc)

	      end do
	    end do


	  end if


	  if (rst_mech == 2) then

	  ! 1st order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc) + &
		pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)

	      end do
	    end do

	  end if
	end if

	if (rst_elec == 2) then

	  if (rst_mech == 0) then

	  ! 2nd order for electric, 0th order for mechanical


	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d)

	      end do
	    end do

	  end if


	  if (rst_mech == 1) then

	  ! 2nd order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d) + &
		pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc)

	      end do
	    end do

	   end if




	  if (rst_mech == 2) then

	  ! 2nd order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3

		perm = make_perm(2, fact(2), (/i, j/))

		alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
		pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d) + &
		pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc) + &
		pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)

	      end do
	    end do


	  end if
	end if

      elseif (present(rst_elec)) then

	if (rst_elec == 0) then

	! 0th order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	    end do
	  end do


	end if


	if (rst_elec == 1) then

	! 1st order for electric, 0th order for mechanical

	  write(*,*) 'Warning: Due to vanishing terms, no precision is gained above 0th order.'

	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	    end do
	  end do



	end if


	if (rst_elec == 2) then

	! 2nd order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
	      pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d)

	    end do
	  end do

	end if


      elseif (present(rst_mech)) then

	if (rst_mech == 0) then

	! 0th order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	    end do
	  end do

	end if


	if (rst_mech == 1) then

	! 0th order for electric, 1st order for mechanical

	  write(*,*) 'Warning: Due to vanishing terms, no precision is gained above 0th order.'

	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	    end do
	  end do


	end if


	if (rst_mech == 2) then

	! 0th order for electric, 2nd order for mechanical


	  do i = 1, 3
	    do j = 1, 3

	      perm = make_perm(2, fact(2), (/i, j/))

	      alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d) + &
	      pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)

	    end do
	  end do


	end if

      end if

    else

    ! No order specification, so do 0th order

      do i = 1, 3
	do j = 1, 3

	  perm = make_perm(2, fact(2), (/i, j/))

	  alpha_pv(i, j) = alpha_pv(i, j) + pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)

	end do
      end do

    end if

! write(*,*) 'Calculated alpha PV'
! 
! do i = 1, 3
! 
! write(*,*)  real(alpha_pv(i,:))
! 
! end do

  end function


  function beta_pv(n_nm, nm_e, opt_e, dm_1d, dm_2d, dm_3d, po_1d, po_2d, po_3d, fcc, fcq, rst_elec, rst_mech, rst_order)
!-------------------------------------------------------------------------------
! Calculates and returns the pure vibrational contribution to the molecular
! 1st hyperpolarizability tensor from various derivatives taken with 
! respect to the normal coordinates of the molecule. Derivatives above 
! first order are taken as optional, should their implementation be 
! unavailable. 
!
! Functions labeled as "Bishop/Luis/Kirtman" terms reproduce terms presented in
! Journal of Chemical Physics, Vol. 108, No. 24 (1998), 10013-17.
! 
!-------------------------------------------------------------------------------

    implicit none

      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i, j, k
      integer, optional,					intent(in) 	:: rst_elec ! Restriction on order of electrical anharmonicity
									    ! If defined: Restricted to [rst_elec]th order, otherwise max order
      integer, optional,					intent(in) 	:: rst_mech ! Restriction on order of mechanical anharmonicity
									    ! If defined: Restricted to [rst_mech]th order, otherwist max order
      integer, optional,					intent(in) 	:: rst_order ! Restriction on total order of anharmonicity
									    ! If defined: Restricted to [rst_order]th order, otherwist max order
      integer, dimension(6, 3)							:: perm
      real(8), dimension(3)							:: opt_e ! Photon energies, negative sum first, then incident
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 3)						:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 3)						:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3), optional, 		intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3), optional,	intent(in) 	:: dm_3d ! 3rd                   "
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3), optional,		intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3), optional,	intent(in) 	:: po_3d ! 3rd   
      complex(8), dimension(n_nm, n_nm, n_nm), optional, 		intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm), optional,	intent(in) 	:: fcq   ! Quartic force constants
      complex(8), dimension(3, 3, 3)						:: beta_pv


    ! Precalculate terms for use in Bishop/Luis/Kirtman term functions

    nm_e_inv = invvec(n_nm, nm_e)

    lambda_a = make_lambda_a(3, n_nm, opt_e, nm_e)
    lambda_ab = make_lambda_ab(3, n_nm, opt_e, nm_e)


    do i = 1, 3
      do j = 1, 3
	do k = 1, 3

	  beta_pv(i, j, k) = 0.0
 
	end do
      end do
    end do


    ! Decide between various configurations of anharmonicity order restrictions, do calculation accordingly

    if (present(rst_order)) then

      if (rst_order == 1) then

      ! Take terms up to and including 1st total order of anharmonicity

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3

	      perm = make_perm(3, fact(3), (/i, j, k/))

	      beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
	      lambda_a, dm_1d, po_1d) + &
	      pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
	      pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc)
      
	    end do
	  end do
	end do


      end if

      if (rst_order == 2) then

      ! Take terms up to and including 2nd total order of anharmonicity

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3

	      perm = make_perm(3, fact(3), (/i, j, k/))

	      beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, lambda_a, dm_1d, po_1d) + &
	      pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
	      pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
	      pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
	      dm_2d, dm_3d, po_1d, po_2d, po_3d) + &
	      pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, po_1d, po_2d, fcc) + &
	      pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)

	    end do
	  end do
	end do

      end if


    elseif (present(rst_elec) .OR. present(rst_mech)) then

      if (present(rst_elec) .AND. present(rst_mech)) then

	if (rst_elec == 0) then

	  if (rst_mech == 0) then

	  ! 0th order for electric, 0th order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  ! Make permutation indices
	      
		  perm = make_perm(3, fact(3), (/i, j, k/))

		  ! Add Bishop/Luis/Kirtman terms to the array to be returned
  
		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d)

		end do
	      end do
	    end do

	  end if



	  if (rst_mech == 1) then

	  ! 0th order for electric, 1st order for mechanical

	    do i = 1, 3
 	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm,&
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc)
		  
		end do
	      end do
	    end do

	  end if




	  if (rst_mech == 2) then

	  ! 0th order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		  pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)

		end do
	      end do
	    end do

	  end if
	end if


	if (rst_elec == 1) then

	  if (rst_mech == 0) then

	  ! 1st order for electric, 0th order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d)
		  
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 1) then

	  ! 1st order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		  pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, po_1d, po_2d, fcc)
		  
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 2) then

	  ! 1st order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		  pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, po_1d, po_2d, fcc) + &
		  pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)

		end do
	      end do
	    end do

	  end if
	end if

	if (rst_elec == 2) then

	  if (rst_mech == 0) then

	  ! 2nd order for electric, 0th order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		  pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, dm_3d, po_1d, po_2d, po_3d)
		  
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 1) then

	  ! 2nd order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		  lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		  pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, dm_3d, po_1d, po_2d, po_3d) + &
		  pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, po_1d, po_2d, fcc)

		end do
	      end do
	    end do  

	  end if


	  if (rst_mech == 2) then

	  ! 2nd order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3

		  perm = make_perm(3, fact(3), (/i, j, k/))

		  beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, lambda_a, dm_1d, po_1d) + &
		  pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		  pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		  pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		  dm_2d, dm_3d, po_1d, po_2d, po_3d) + &
		  pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, po_1d, po_2d, fcc) + &
		  pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)

		end do
	      end do
	    end do

	  end if
	end if

      elseif (present(rst_elec)) then

	if (rst_elec == 0) then

	! 0th order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		lambda_a, dm_1d, po_1d)

	      end do
	    end do
	  end do

	end if


	if (rst_elec == 1) then

	! 1st order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		lambda_a, dm_1d, po_1d) + &
		pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d)
		  
	      end do
	    end do
	  end do

	end if


	if (rst_elec == 2) then

	! 2nd order for electric, 0th order for mechanical
	  
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		lambda_a, dm_1d, po_1d) + &
		pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d) + &
		pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		dm_2d, dm_3d, po_1d, po_2d, po_3d)
		  
	      end do
	    end do
	  end do

	end if


      elseif (present(rst_mech)) then

	if (rst_mech == 0) then

	! 0th order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		lambda_a, dm_1d, po_1d)

	      end do
	    end do
	  end do

	end if


	if (rst_mech == 1) then

	! 0th order for electric, 1st order for mechanical

	  do i = 1, 3
 	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm,&
		lambda_a, dm_1d, po_1d) + &
		pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc)
		  
	      end do
	    end do
	  end do


	end if


	if (rst_mech == 2) then

	! 0th order for electric, 2nd order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3

		perm = make_perm(3, fact(3), (/i, j, k/))

		beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
		lambda_a, dm_1d, po_1d) + &
		pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc) + &
		pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)

	      end do
	    end do
	  end do

	end if

      end if

    else

    ! No order specification, so do 0th order
      do i = 1, 3
	do j = 1, 3
	  do k = 1, 3

	    perm = make_perm(3, fact(3), (/i, j, k/))

	    beta_pv(i, j, k) = beta_pv(i, j, k) + pv_dippol_00(n_nm, perm, &
	    lambda_a, dm_1d, po_1d)

	  end do
	end do
      end do

    end if


! write(*,*) 'Calculated beta PV'
! 
! do i = 1, 3
! do j = 1, 3
! 
! write(*,*)  real(beta_pv(i,j,:))
! 
! end do
! write(*,*) ' '
! end do

  end function


  function gamma_pv(n_nm, nm_e, opt_e, dm_1d, dm_2d, dm_3d, po_1d, po_2d, po_3d, hp_1d, hp_2d, hp_3d, fcc, fcq, &
  rst_elec, rst_mech, rst_order)
!-------------------------------------------------------------------------------
! Calculates and returns the pure vibrational contribution to the molecular
! 2nd hyperpolarizability tensor from various derivatives taken with 
! respect to the normal coordinates of the molecule. Derivatives above 
! first order are taken as optional, should their implementation be 
! unavailable. 
!
! Functions labeled as "Bishop/Luis/Kirtman" terms reproduce terms presented in
! Journal of Chemical Physics, Vol. 108, No. 24 (1998), 10013-17.
! 
!-------------------------------------------------------------------------------
    implicit none
      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l
      integer, optional,					intent(in) 	:: rst_elec ! Restriction on order of electrical anharmonicity
									    ! If defined: Restricted to [rst_elec]th order, otherwise max order
      integer, optional,					intent(in) 	:: rst_mech ! Restriction on order of mechanical anharmonicity
									    ! If defined: Restricted to [rst_mech]th order, otherwist max order
      integer, optional,					intent(in) 	:: rst_order ! Restriction on total order of anharmonicity
									    ! If defined: Restricted to [rst_order]th order, otherwist max order
      integer, dimension(24, 4)							:: perm
      
      real(8), dimension(4)							:: opt_e ! Photon energies, negative sum first, then incident
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4)						:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4)						:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4)						:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4, 4)					:: lambda_abij ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3), optional, 		intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3), optional,	intent(in) 	:: dm_3d ! 3rd                   "
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3), optional,		intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3), optional,	intent(in) 	:: po_3d ! 3rd   
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st            "              1st hyperpolarizability
      complex(8), dimension(n_nm, n_nm, 3, 3, 3), optional,	intent(in) 	:: hp_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3, 3), optional,	intent(in) 	:: hp_3d ! 3rd                   "
      complex(8), dimension(n_nm, n_nm, n_nm), optional, 		intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm), optional,	intent(in) 	:: fcq   ! Quartic force constants
      complex(8), dimension(3, 3, 3, 3)						:: gamma_pv


    ! Precalculate terms for use in Bishop/Luis/Kirtman term functions

    nm_e_inv = invvec(n_nm, nm_e)

    lambda_a = make_lambda_a(4, n_nm, opt_e, nm_e)
    lambda_ab = make_lambda_ab(4, n_nm, opt_e, nm_e)
    lambda_aij = make_lambda_aij(4, n_nm, opt_e, nm_e)
    lambda_abij = make_lambda_abij(4, n_nm, opt_e, nm_e)

    do i = 1, 3
      do j = 1, 3
	do k = 1, 3
	  do l = 1, 3
	  
	    gamma_pv(i, j, k, l) = 0.0
    
	  end do
	end do
      end do
    end do



    ! Decide between various configurations of anharmonicity order restrictions, do calculation accordingly

    if (present(rst_order)) then

      if (rst_order == 1) then

      ! Take terms up to and including 1st total order of anharmonicity

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      do l = 1, 3

		perm = make_perm(4, fact(4), (/i, j, k, l/))

		gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc)

	      end do
	    end do
	  end do
	end do

      end if

      if (rst_order == 2) then

      ! Take terms up to and including 2nd total order of anharmonicity

	do i = 1, 3
	  do j = 1, 3
	    do k = 1, 3
	      do l = 1, 3

		perm = make_perm(4, fact(4), (/i, j, k, l/))

		gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
		po_1d, po_2d, po_3d) + &
		pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc) + &
		pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq) + &
		pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, & 
		dm_3d, hp_1d, hp_2d, hp_3d) + &
		pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, hp_1d, hp_2d, fcc) + &
		pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq) + &
		pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d) + &
		pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc) + &
		pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)

	      end do
	    end do
	  end do
	end do

      end if


    elseif (present(rst_elec) .OR. present(rst_mech)) then

      if (present(rst_elec) .AND. present(rst_mech)) then

	if (rst_elec == 0) then

	  if (rst_mech == 0) then

	  ! 0th order for electric, 0th order for mechanical


	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    ! Make permutation indices

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    ! Add Bishop/Luis/Kirtman terms to the array to be returned

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d)
		    
		  end do
		end do
	      end do
	    end do

	  end if



	  if (rst_mech == 1) then

	  ! 0th order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc)
		    
		  end do
		end do
	      end do
	    end do

	  end if




	  if (rst_mech == 2) then

	  ! 0th order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		    pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq) + &
		    pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq) + &
		    pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)

		  end do
		end do
	      end do
	    end do
  
	  end if

	end if

	if (rst_elec == 1) then

	  if (rst_mech == 0) then

	  ! 1st order for electric, 0th order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d)

		  end do
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 1) then

	  ! 1st order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		    pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc) + &
		    pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		    dm_2d, hp_1d, hp_2d, fcc) + &
		    pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc)
		
		  end do
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 2) then

	  ! 1st order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

! Some strange line remained here and was removed. Check here early if anything goes wrong.

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		    pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc) + &
		    pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq) + &
		    pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, hp_1d, hp_2d, fcc) + &
		    pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq) + &
		    pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc) + &
		    pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)

		  end do
		end do
	      end do
	    end do

	  end if
	end if

	if (rst_elec == 2) then

	  if (rst_mech == 0) then

	  ! 2nd order for electric, 0th order for mechanical


	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		    pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
		    po_1d, po_2d, po_3d) + &
		    pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, & 
		    dm_3d, hp_1d, hp_2d, hp_3d) + &
		    pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d)

		  end do
		end do
	      end do
	    end do


	  end if


	  if (rst_mech == 1) then

	  ! 2nd order for electric, 1st order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		    pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		    pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
		    po_1d, po_2d, po_3d) + &
		    pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc) + &
		    pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, & 
		    dm_3d, hp_1d, hp_2d, hp_3d) + &
		    pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, &
		    dm_2d, hp_1d, hp_2d, fcc) + &
		    pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d) + &
		    pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc)

		  end do
		end do
	      end do
	    end do

	  end if


	  if (rst_mech == 2) then

	  ! 2nd order for electric, 2nd order for mechanical

	    do i = 1, 3
	      do j = 1, 3
		do k = 1, 3
		  do l = 1, 3

		    perm = make_perm(4, fact(4), (/i, j, k, l/))

		    gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		    pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		    pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		    pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		    pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
		    po_1d, po_2d, po_3d) + &
		    pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc) + &
		    pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq) + &
		    pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, & 
		    dm_3d, hp_1d, hp_2d, hp_3d) + &
		    pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, hp_1d, hp_2d, fcc) + &
		    pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq) + &
		    pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d) + &
		    pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc) + &
		    pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)

		  end do
		end do
	      end do
	    end do

	  end if
	end if

      elseif (present(rst_elec)) then

	if (rst_elec == 0) then

	! 0th order for electric, 0th order for mechanical


	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d)
		    
		end do
	      end do
	    end do
	  end do


	end if


	if (rst_elec == 1) then

	! 1st order for electric, 0th order for mechanical


	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		  pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d)

		end do
	      end do
	    end do
	  end do


	end if


	if (rst_elec == 2) then

	! 2nd order for electric, 0th order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		  pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d) + &
		  pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
		  po_1d, po_2d, po_3d) + &
		  pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, & 
		  dm_3d, hp_1d, hp_2d, hp_3d) + &
		  pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d)
		  

		end do
	      end do
	    end do
	  end do

	end if


      elseif (present(rst_mech)) then

	if (rst_mech == 0) then

	! 0th order for electric, 0th order for mechanical


	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d)
		    
		end do
	      end do
	    end do
	  end do


	end if


	if (rst_mech == 1) then

	! 0th order for electric, 1st order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		  pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc)
		    
		end do
	      end do
	    end do
	  end do

	end if


	if (rst_mech == 2) then

	! 0th order for electric, 2nd order for mechanical

	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
		do l = 1, 3

		  perm = make_perm(4, fact(4), (/i, j, k, l/))

		  gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
		  pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
		  pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d) + &
		  pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc) + &
		  pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq) + &
		  pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq) + &
		  pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)

		end do
	      end do
	    end do
	  end do

	end if

      end if

    else

    ! No order specification, so do 0th order

      do i = 1, 3
	do j = 1, 3
	  do k = 1, 3
	    do l = 1, 3

	      perm = make_perm(4, fact(4), (/i, j, k, l/))

	      gamma_pv(i, j, k, l) = gamma_pv(i, j, k, l) + &
	      pv_polsq_00(n_nm, perm, lambda_aij, po_1d) + &
	      pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d)
		    
	    end do
	  end do
	end do
      end do

    end if


! write(*,*) 'Calculated gamma PV'
! 
! do i = 1, 3
! do j = 1, 3
! do k = 1, 3
! 
! write(*,*)  real(gamma_pv(i,j,k,:))
! 
! end do 
! write(*,*) ' '
! end do
! write(*,*) ' '
! end do

  end function


!-------------------------------------------------------------------------------
!
! Contains "Bishop/Luis/Kirtman" terms reproducing the terms presented
! in Journal of Chemical Physics, Vol. 108, No. 24 (1998), 10013-17.
! These functions are called from pvalpha.f90, pvbeta.f90 and pvgamma.f90.
! 
!-------------------------------------------------------------------------------


  function pv_dipsq_00(n_nm, perm, lambda_a, dm_1d)
      ! Bishop/Luis/Kirtman's mu^2 (0,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, m
      integer, dimension(2, 2)							:: perm, permopt
      
      real(8), dimension(n_nm, 2),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8)									:: pv_dipsq_00 


      permopt = make_perm(2, fact(2), (/1, 2/))

      pv_dipsq_00 = 0.0

! Looping over normal modes with i, and permutations with m

    do i=1, n_nm
      do m=1, 2

	pv_dipsq_00 = pv_dipsq_00 + dm_1d(i, perm(m, 1)) * dm_1d(i, perm(m, 2)) * lambda_a(i, permopt(m, 1))

      end do
    end do

! Multiplying the coefficient in front of the summation

    pv_dipsq_00 = 0.5 * pv_dipsq_00

  end function

! The rest of these functions follow the same pattern

  function pv_dipsq_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d)
      ! Bishop/Luis/Kirtman's mu^2 (2,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(2, 2)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 2),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 2),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),	 		intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3),			intent(in) 	:: dm_3d ! 3rd                   "
      complex(8)									:: pv_dipsq_20 


      permopt = make_perm(2, fact(2), (/1, 2/))

      pv_dipsq_20 = 0.0

    hbar = 1.0
    
    do i=1, n_nm
      do j=1, n_nm
	do m=1, 2

	  pv_dipsq_20 = pv_dipsq_20 + dm_2d(i, j, perm(m, 1)) * dm_2d(i, j, perm(m, 2)) * &
	  (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) + &
	  dm_3d(i, i, j, perm(m, 1)) * dm_1d(i, perm(m, 2)) * nm_e_inv(i) * lambda_a(j, permopt(m, 1)) + &
	  dm_1d(i, perm(m, 1)) * dm_3d(i, i, j, perm(m, 2)) * nm_e_inv(i) * lambda_a(j, permopt(m, 1))
	
	end do
      end do
    end do

    pv_dipsq_20 = hbar * 0.125 * pv_dipsq_20



  end function

  function pv_dipsq_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, fcc)
      ! Bishop/Luis/Kirtman's mu^2 (1,1) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(2, 2)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 2),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 2),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_dipsq_11 

      permopt = make_perm(2, fact(2), (/1, 2/))

      pv_dipsq_11 = 0.0

    hbar = 1.0
    
    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm

	  do m=1, 2

! write(*,*) 'i j k ', i, j, k
! write(*,*) 'added ', fcc(i, j, k) * dm_2d(i, j, perm(m, 1)) * dm_1d(k, perm(m, 2)) * &
! 	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) 
	    
! write(*,*) 'added also', fcc(j, k, k) * dm_2d(i, j, perm(m, 1)) * dm_1d(i, perm(m, 2)) * (nm_e_inv(j)**(2.0)) * &
! 	    nm_e_inv(k) * lambda_a(i, permopt(m, 1))
	
	    pv_dipsq_11 = pv_dipsq_11 + fcc(i, j, k) * dm_2d(i, j, perm(m, 1)) * dm_1d(k, perm(m, 2)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) + &
	    fcc(j, k, k) * dm_2d(i, j, perm(m, 1)) * dm_1d(i, perm(m, 2)) * (nm_e_inv(j)**(2.0)) * &
	    nm_e_inv(k) * lambda_a(i, permopt(m, 1))
	   
	  end do
	end do
      end do
    end do

    pv_dipsq_11 = -hbar * 0.25 * pv_dipsq_11

  end function

  function pv_dipsq_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, fcc, fcq)
      ! Bishop/Luis/Kirtman's mu^2 (0,2) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m
      integer, dimension(2, 2)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 2),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 2),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm),		intent(in) 	:: fcq   ! Quartic force constants
      complex(8)									:: pv_dipsq_02 

      permopt = make_perm(2, fact(2), (/1, 2/))

      pv_dipsq_02 = 0.0

    hbar = 1.0
    
    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 2
	
	    pv_dipsq_02 = pv_dipsq_02 + nm_e_inv(i) * fcq(i, i, j, k) * dm_1d(j, perm(m, 1)) * &
	    dm_1d(k, perm(m, 2)) * lambda_a(j, permopt(m, 1)) * lambda_a(k, permopt(m, 1))
	      
	      do l=1, n_nm

		pv_dipsq_02 = pv_dipsq_02 - fcc(i, i, j) * fcc(j, k, l) * dm_1d(k, perm(m, 1)) * &
		dm_1d(l, perm(m, 2)) * (nm_e_inv(j)**(2.0)) * lambda_a(k, permopt(m, 1)) * &
		lambda_a(l, permopt(m, 1))  + 2.0 * fcc(i, j, k) * fcc(i, j, l) * dm_1d(k, perm(m, 1)) * &
		dm_1d(l, perm(m, 2)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) * &
		lambda_a(l, permopt(m, 1))


	      end do
	   
	  end do
	end do
      end do
    end do

    pv_dipsq_02 = -hbar * 0.125 * pv_dipsq_02



  end function

  function pv_dipcu_10(n_nm, perm, lambda_a, dm_1d, dm_2d)
      ! Bishop/Luis/Kirtman's mu^3 (1,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),	 		intent(in) 	:: dm_2d ! 2nd                   "
      complex(8)									:: pv_dipcu_10 

      permopt = make_perm(3, fact(3), (/1, 2, 3/))

      
      pv_dipcu_10 = 0.0

    do i=1, n_nm
      do j=1, n_nm
	do m=1, 6
	
	    pv_dipcu_10 = pv_dipcu_10 + dm_1d(i, perm(m, 1)) * dm_2d(i, j, perm(m,2)) * &
	    dm_1d(j, perm(m, 3)) * lambda_a(i, permopt(m, 1)) * lambda_a(j, permopt(m, 3))
	   
	end do
      end do
    end do

    pv_dipcu_10 = 0.5 * pv_dipcu_10
  
  end function

  function pv_dipcu_01(n_nm, perm, lambda_a, dm_1d, fcc)
      ! Bishop/Luis/Kirtman's mu^3 (0,1) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_dipcu_01

      permopt = make_perm(3, fact(3), (/1, 2, 3/))
      
    
      pv_dipcu_01 = 0.0

    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 6
	
	    pv_dipcu_01 = pv_dipcu_01 + fcc(i, j, k) * dm_1d(i, perm(m, 1)) * dm_1d(j, perm(m, 2)) * &
	    dm_1d(k, perm(m, 3)) * lambda_a(i, permopt(m, 1)) * lambda_a(j, permopt(m, 2)) * &
	    lambda_a(k, permopt(m, 3))
	   
	  end do
	end do
      end do
    end do

    pv_dipcu_01 = -1.0/(6.0) * pv_dipcu_01


  end function

  function pv_dipsqpol_10(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, po_1d, po_2d)
      ! Bishop/Luis/Kirtman's mu^2 alpha (1,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8)									:: pv_dipsqpol_10 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_dipsqpol_10 = 0.0

    do i=1, n_nm
      do j=1, n_nm
	do m= 1, 24
	
	    pv_dipsqpol_10 = pv_dipsqpol_10 + dm_1d(i, perm(m, 1)) * po_2d(i, j, perm(m, 2), perm(m, 3)) * &
	    dm_1d(j, perm(m, 4)) * lambda_a(i, permopt(m, 1)) * lambda_a(j, permopt(m, 4)) + 2.0 * &
	    dm_1d(i, perm(m, 1)) * dm_2d(i, j, perm(m, 2)) * po_1d(j, perm(m, 3), perm(m, 4)) * &
	    lambda_a(i, permopt(m, 1)) * lambda_aij(j, permopt(m, 3), permopt(m, 4))
	   
	end do
      end do
    end do

    pv_dipsqpol_10 = 0.25 * pv_dipsqpol_10
        
  end function

  function pv_dipsqpol_01(n_nm, perm, lambda_a, lambda_aij, dm_1d, po_1d, fcc)
      ! Bishop/Luis/Kirtman's mu^2 alpha (0,1) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(24, 4)							:: perm, permopt
      
      !real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_dipsqpol_01 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_dipsqpol_01 = 0.0


    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_dipsqpol_01 = pv_dipsqpol_01 + fcc(i, j, k) * dm_1d(i, perm(m, 1)) * &
	    dm_1d(j, perm(m, 2)) * po_1d(k, perm(m, 3), perm(m, 4)) * lambda_a(i, permopt(m, 1)) * &
	    lambda_a(j, permopt(m, 2)) * lambda_aij(k, permopt(m,3), permopt(m, 4))
	   
	  end do
	end do
      end do
    end do

    pv_dipsqpol_01 = -0.25 * pv_dipsqpol_01
        
  end function

  function pv_dipqu_20(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, dm_3d)
      ! Bishop/Luis/Kirtman's mu^4 (2,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3),			intent(in) 	:: dm_3d ! 3rd                   "
      complex(8)									:: pv_dipqu_20 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_dipqu_20 = 0.0

        
    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_dipqu_20 = pv_dipqu_20 + 3.0 * dm_1d(i, perm(m, 1)) * dm_2d(i, j, perm(m, 2)) * &
	    dm_2d(j, k, perm(m, 3)) * dm_1d(k, perm(m, 4)) * lambda_a(i, permopt(m, 1)) * &
	    lambda_aij(j, permopt(m, 3), permopt(m, 4)) * lambda_a(k, permopt(m, 4)) + &
	    dm_3d(i, j, k, perm(m, 1)) * dm_1d(i, perm(m, 2)) * dm_1d(j, perm(m, 3)) * &
	    dm_1d(k, perm(m, 4)) * lambda_a(i, permopt(m, 2)) * lambda_a(j, permopt(m, 3)) * &
	    lambda_a(k, permopt(m, 4))
	   
	  end do
	end do
      end do
    end do

    pv_dipqu_20 = (1.0/6.0) * pv_dipqu_20


  end function

  function pv_dipqu_11(n_nm, perm, lambda_a, lambda_aij, dm_1d, dm_2d, fcc)
      ! Bishop/Luis/Kirtman's mu^4 (1,1) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_dipqu_11 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_dipqu_11 = 0.0

        
    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do l=1, n_nm
	    do m=1, 4
	      
	      pv_dipqu_11 = pv_dipqu_11 + fcc(i, j, k) * dm_1d(i, perm(m, 1)) * dm_1d(j, perm(m, 2)) * &
	      dm_2d(k, l, perm(m, 3)) * dm_1d(l, perm(m, 4)) * lambda_a(i, permopt(m, 1)) * &
	      lambda_a(j, permopt(m, 2)) * lambda_aij(k, permopt(m, 3), permopt(m, 4)) * lambda_a(l, permopt(m, 4))
	  
	    end do
	  end do
	end do
      end do
    end do

    pv_dipqu_11 = -0.5 * pv_dipqu_11

  end function

  function pv_dipqu_02(n_nm, perm, lambda_a, lambda_aij, dm_1d, fcc, fcq)
      ! Bishop/Luis/Kirtman's mu^4 (0,2) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m, n
      integer, dimension(24, 4)							:: perm, permopt
      
      !real(8), dimension(n_nm),				intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm),		intent(in) 	:: fcq   ! Quartic force constants
      complex(8)									:: pv_dipqu_02 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_dipqu_02 = 0.0

        
    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do l=1, n_nm
	    do m=1, 24
	      
	      pv_dipqu_02 = pv_dipqu_02 + fcq(i, j, k, l) * dm_1d(i, perm(m, 1)) * &
	      dm_1d(j, perm(m, 2)) * dm_1d(k, perm(m, 3)) * dm_1d(l, perm(m, 4)) * &
	      lambda_a(i, permopt(m, 1)) * lambda_a(j, permopt(m, 2)) * lambda_a(k, permopt(m, 3)) * &
	      lambda_a(l, permopt(m, 4))
	      
	      do n=1, n_nm
		
		pv_dipqu_02 = pv_dipqu_02 - 3.0 * fcc(i, j, k) * fcc(k, l, n) * dm_1d(i, perm(m, 1)) * &
		dm_1d(j, perm(m, 2))  * dm_1d(l, perm(m, 3)) * dm_1d(n, perm(m, 4)) * &
		lambda_a(i, permopt(m, 1)) * lambda_a(j, permopt(m, 2)) * & 
		lambda_aij(k, permopt(m, 3), permopt(m, 4)) * lambda_a(l, permopt(m, 3)) * &
		lambda_a(n, permopt(m, 4))

	      end do
	  
	    end do
	  end do
	end do
      end do
    end do

    pv_dipqu_02 = (-1.0/24.0) * pv_dipqu_02

  end function


  function pv_dippol_00(n_nm, perm, lambda_a, dm_1d, po_1d)
      ! Bishop/Luis/Kirtman's mu alpha (0,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8)									:: pv_dippol_00 

      permopt = make_perm(3, fact(3), (/1, 2, 3/))

      pv_dippol_00 = 0.0


    do i=1, n_nm
      do m=1, 6
	
	pv_dippol_00 = pv_dippol_00 + dm_1d(i, perm(m, 1)) * po_1d(i, perm(m, 2), perm(m, 3)) * &
	lambda_a(i, permopt(m, 1))

      end do
    end do

    pv_dippol_00 = 0.5 * pv_dippol_00
        
  end function

  function pv_dippol_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d, po_1d, po_2d, po_3d)
      ! Bishop/Luis/Kirtman's mu alpha (2,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 3),				intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3),			intent(in) 	:: dm_3d ! 3rd                   "
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3),		intent(in) 	:: po_3d ! 3rd                   "
      complex(8)									:: pv_dippol_20 

      permopt = make_perm(3, fact(3), (/1, 2, 3/))

      pv_dippol_20 = 0.0


    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do m=1, 6
	
	  pv_dippol_20 = pv_dippol_20 + dm_2d(i, j, perm(m, 1)) * po_2d(i, j, perm(m, 2), perm(m, 3)) * &
	  (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) + &
	  dm_3d(i, i, j, perm(m, 1)) * po_1d(i, perm(m, 2), perm(m, 3)) * nm_e_inv(i) * &
	  lambda_a(j, permopt(m, 1)) + dm_1d(i, perm(m, 1)) * & 
	  po_3d(i, i, j, perm(m, 2), perm(m, 3)) * nm_e_inv(i) * lambda_a(j, permopt(m, 1))
	
	end do
      end do
    end do

    pv_dippol_20 = hbar * 0.125 * pv_dippol_20

        
  end function

  function pv_dippol_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, po_1d, po_2d, fcc)
      ! Bishop/Luis/Kirtman's mu alpha (1,1) term
      ! This implementation may be erroneous, please check again if anything goes wrong

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 3),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_dippol_11 

      permopt = make_perm(3, fact(3), (/1, 2, 3/))

      pv_dippol_11 = 0.0

    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 6
	
	    pv_dippol_11 = pv_dippol_11 + fcc(i, j, k) * dm_2d(i, j, perm(m, 1)) * po_1d(k, perm(m, 2), perm(m, 3)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) + &    
	    fcc(i, j, k) * po_2d(i, j, perm(m, 2), perm(m, 3)) * dm_1d(k, perm(m, 1)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) + &
	    fcc(j, k, k) * dm_2d(i, j, perm(m, 1)) * po_1d(i, perm(m, 2), perm(m, 3)) * (nm_e_inv(j)**(2.0)) * &
	    nm_e_inv(k) * lambda_a(i, permopt(m, 1)) + &	    
	    fcc(j, k, k) * po_2d(i, j, perm(m, 2), perm(m, 3)) * dm_1d(i, perm(m, 1)) * (nm_e_inv(j)**(2.0)) * &
	    nm_e_inv(k) * lambda_a(i, permopt(m, 1))
	   
	  end do
	end do
      end do
    end do

    pv_dippol_11 = -hbar * 0.125 * pv_dippol_11
        
  end function
  
  function pv_dippol_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, po_1d, fcc, fcq)
      ! Bishop/Luis/Kirtman's mu alpha (0,2) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m
      integer, dimension(6, 3)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 3),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 3),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st            "              polarizability
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm),		intent(in) 	:: fcq   ! Quartic force constants
      complex(8)									:: pv_dippol_02 

      permopt = make_perm(3, fact(3), (/1, 2, 3/))

      pv_dippol_02 = 0.0

    hbar = 1.0


    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 6
	
	    pv_dippol_02 = pv_dippol_02 + nm_e_inv(i) * fcq(i, i, j, k) * dm_1d(j, perm(m, 1)) * &
	    po_1d(k, perm(m, 2), perm(m, 3)) * lambda_a(j, permopt(m, 1)) * lambda_a(k, permopt(m, 1))
	      
	      do l=1, n_nm

		pv_dippol_02 = pv_dippol_02 - fcc(i, i, j) * fcc(j, k, l) * dm_1d(k, perm(m, 1)) * &
		po_1d(l, perm(m, 2), perm(m, 3)) * (nm_e_inv(j)**(2.0)) * lambda_a(k, permopt(m, 1)) * &
		lambda_a(l, permopt(m, 1))  + 2.0 * fcc(i, j, k) * fcc(i, j, l) * dm_1d(k, perm(m, 1)) * &
		po_1d(l, perm(m, 2), perm(m, 3)) * lambda_ab(i, j, permopt(m, 1)) * &
		lambda_a(k, permopt(m, 1)) * lambda_a(l, permopt(m, 1))

	      end do
	   
	  end do
	end do
      end do
    end do

    pv_dippol_02 = -hbar * 0.125 * pv_dippol_02
        
  end function



  function pv_polsq_00(n_nm, perm, lambda_aij, po_1d)
      ! Bishop/Luis/Kirtman's alpha^2 (0,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, m
      integer, dimension(24, 4)							:: perm, permopt
      
      !real(8), dimension(n_nm),				intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st normal mode derivative of polarizability
      complex(8)									:: pv_polsq_00 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_polsq_00 = 0.0

        
    do i=1, n_nm
      do m=1, 24
	
	pv_polsq_00 = pv_polsq_00 + po_1d(i, perm(m, 1), perm(m, 2)) * &
	po_1d(i, perm(m, 3), perm(m, 4)) * lambda_aij(i, permopt(m, 3), permopt(m, 4))

      end do
    end do

    pv_polsq_00 = 0.125 * pv_polsq_00

  end function

  function pv_polsq_20(n_nm, perm, lambda_a, lambda_ab, lambda_aij, lambda_abij, nm_e_inv, &
  po_1d, po_2d, po_3d)
      ! Bishop/Luis/Kirtman's alpha^2 (2,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4, 4),			intent(in)	:: lambda_abij ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st normal mode derivative of polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3),		intent(in) 	:: po_3d ! 3rd   
      complex(8)									:: pv_polsq_20 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_polsq_20 = 0.0


    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do m=1, 24
	
	  pv_polsq_20 = pv_polsq_20 + po_2d(i, j, perm(m, 1), perm(m, 2)) * &
	  po_2d(i, j, perm(m, 3), perm(m, 4)) * &
	  (nm_e_inv(i) + nm_e_inv(j)) * lambda_abij(i, j, permopt(m, 3), permopt(m, 4)) + &
	  po_3d(i, i, j, perm(m, 1), perm(m, 2)) * po_1d(i, perm(m, 3), perm(m, 4)) * & 
	  nm_e_inv(i) * lambda_aij(j, permopt(m, 3), permopt(m, 4)) + &
	  po_1d(i, perm(m, 1), perm(m, 2)) * po_3d(i, i, j, perm(m, 3), perm(m, 4)) * &
	  nm_e_inv(i) * lambda_aij(j, permopt(m, 3), permopt(m, 4))
	
	end do
      end do
    end do

    pv_polsq_20 = hbar * 0.03125 * pv_polsq_20

        
  end function

  function pv_polsq_11(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, po_2d, fcc)
      ! Bishop/Luis/Kirtman's alpha^2 (1,1) term
      ! This implementation may be erroneous, please check again if anything goes wrong

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4, 4),			intent(in)	:: lambda_abij ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st normal mode derivative of polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_polsq_11 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_polsq_11 = 0.0

    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_polsq_11 = pv_polsq_11 + fcc(i, j, k) * po_2d(i, j, perm(m, 1), perm(m, 2)) * &
	    po_1d(k, perm(m, 3), perm(m, 4)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_abij(i, j, permopt(m, 3), permopt(m, 4)) *  & 
	    lambda_aij(k, permopt(m, 3), permopt(m, 4)) + &
	    fcc(j, k, k) * po_2d(i, j, perm(m, 1), perm(m, 2)) * po_1d(i, perm(m, 3), perm(m, 4)) * &
	    (nm_e_inv(j)**(2.0)) *  nm_e_inv(k) * lambda_aij(i, permopt(m, 3), permopt(m, 4))
	   
	  end do
	end do
      end do
    end do

    pv_polsq_11 = -hbar * 0.0625 * pv_polsq_11
        
  end function

  function pv_polsq_02(n_nm, perm, lambda_aij, lambda_abij, nm_e_inv, po_1d, fcc, fcq)
      ! Bishop/Luis/Kirtman's alpha^2 (0,2) term
      ! Original article is inconsistent with subscripts for lambda
      ! The term lambda_alpha is taken to be lambda_d

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4, 4),				intent(in)	:: lambda_aij  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4, 4),			intent(in)	:: lambda_abij ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3, 3),				intent(in) 	:: po_1d ! 1st normal mode derivative of polarizability
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm),		intent(in) 	:: fcq   ! Quartic force constants
      complex(8)									:: pv_polsq_02 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_polsq_02 = 0.0

    hbar = 1.0


    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_polsq_02 = pv_polsq_02 + nm_e_inv(i) * fcq(i, i, j, k) * po_1d(j, perm(m, 1), perm(m, 2)) * &
	    po_1d(k, perm(m, 3), perm(m, 4)) * lambda_aij(j, permopt(m, 3), permopt(m, 4)) * &
	    lambda_aij(k, permopt(m, 3), permopt(m, 4))
	      
	      do l=1, n_nm

		pv_polsq_02 = pv_polsq_02 - fcc(i, i, j) * fcc(j, k, l) *  & 
		po_1d(k, perm(m, 1), perm(m, 2)) * po_1d(l, perm(m, 3), perm(m, 4)) * &
		(nm_e_inv(j)**(2.0)) * lambda_aij(k, permopt(m, 3), permopt(m, 4)) * &
		lambda_aij(l, permopt(m, 3), permopt(m, 4))  + 2.0 * fcc(i, j, k) * fcc(i, j, l) * &
		po_1d(k, perm(m, 1), perm(m, 2)) * po_1d(l, perm(m, 3), perm(m, 4)) * &
		lambda_abij(i, j, permopt(m, 3), permopt(m, 4)) * &
		lambda_aij(k, permopt(m, 3), permopt(m, 4)) * lambda_aij(l, permopt(m, 3), permopt(m, 4))

	      end do
	   
	  end do
	end do
      end do
    end do

    pv_polsq_02 = -hbar * 0.03125 * pv_polsq_02
        
  end function


  function pv_diphyp_00(n_nm, perm, lambda_a, dm_1d, hp_1d)
      ! Bishop/Luis/Kirtman's mu beta (0,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st            "              1st hyperpolarizability
      complex(8)									:: pv_diphyp_00 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_diphyp_00 = 0.0

        
    do i=1, n_nm
      do m=1, 24
	
	pv_diphyp_00 = pv_diphyp_00 + dm_1d(i, perm(m, 1)) * &
	hp_1d(i, perm(m, 2), perm(m, 3), perm(m, 4)) * lambda_a(i, permopt(m, 1))

      end do
    end do

    pv_diphyp_00 = (1.0/6.0) * pv_diphyp_00

  end function

  function pv_diphyp_20(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, dm_3d, hp_1d, hp_2d, hp_3d)
      ! Bishop/Luis/Kirtman's mu beta (2,0) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3),			intent(in) 	:: dm_3d ! 3rd                   "
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st            "              1st hyperpolarizability
      complex(8), dimension(n_nm, n_nm, 3, 3, 3),			intent(in) 	:: hp_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm, 3, 3, 3),		intent(in) 	:: hp_3d ! 3rd                   "
      complex(8)									:: pv_diphyp_20 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_diphyp_20 = 0.0


    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do m=1, 24
	
	  pv_diphyp_20 = pv_diphyp_20 + dm_2d(i, j, perm(m, 1)) * & 
	  hp_2d(i, j, perm(m, 2), perm(m, 3), perm(m, 4)) * (nm_e_inv(i) + nm_e_inv(j)) * &
	  lambda_ab(i, j, permopt(m, 1)) + dm_3d(i, i, j, perm(m, 1)) * &
	  hp_1d(i, perm(m, 2), perm(m, 3), perm(m, 4)) * nm_e_inv(i) * lambda_a(j, permopt(m, 1)) + &
	  dm_1d(i, perm(m, 1)) * hp_3d(i, i, j, perm(m, 2), perm(m, 3), perm(m, 4)) * &
	  nm_e_inv(i) * lambda_a(j, permopt(m, 1))
	
	end do
      end do
    end do

    pv_diphyp_20 = hbar * (1.0/24.0) * pv_diphyp_20

        
  end function

  function pv_diphyp_11(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, dm_2d, hp_1d, hp_2d, fcc)
      ! Bishop/Luis/Kirtman's mu beta (1,1) term
      ! This implementation may be erroneous, please check again if anything goes wrong
      ! Nothing seemed to go wrong in testing, so assume it's OK for now

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, n_nm, 3),			intent(in) 	:: dm_2d ! 2nd                   "
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st            "              1st hyperpolarizability
      complex(8), dimension(n_nm, n_nm, 3, 3, 3),			intent(in) 	:: hp_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8)									:: pv_diphyp_11 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_diphyp_11 = 0.0

    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_diphyp_11 = pv_diphyp_11 + fcc(i, j, k) * dm_2d(i, j, perm(m, 1)) *  &
	    hp_1d(k, perm(m, 2), perm(m, 3), perm(m, 4)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) + &
	    fcc(i, j, k) * hp_2d(i, j, perm(m, 2), perm(m, 3), perm(m, 4)) * &
	    dm_1d(k, perm(m, 1)) * &
	    (nm_e_inv(i) + nm_e_inv(j)) * lambda_ab(i, j, permopt(m, 1)) * lambda_a(k, permopt(m, 1)) + &
	    fcc(j, k, k) * dm_2d(i, j, perm(m, 1)) * hp_1d(i, perm(m, 2), perm(m, 3), perm(m, 4)) * & 
	    (nm_e_inv(j)**(2.0)) * nm_e_inv(k) * lambda_a(i, permopt(m, 1)) + &
	    fcc(j, k, k) * hp_2d(i, j, perm(m, 2), perm(m, 3), perm(m, 4)) * & 
	    dm_1d(i, perm(m, 1)) * (nm_e_inv(j)**(2.0)) * nm_e_inv(k) * lambda_a(i, permopt(m, 1))

	   
	  end do
	end do
      end do
    end do

    pv_diphyp_11 = -hbar * (1.0/24.0) * pv_diphyp_11
        
  end function

  function pv_diphyp_02(n_nm, perm, lambda_a, lambda_ab, nm_e_inv, dm_1d, hp_1d, fcc, fcq)
      ! Bishop/Luis/Kirtman's mu beta (0,2) term

      implicit none
      integer,							intent(in)	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, m
      integer, dimension(24, 4)							:: perm, permopt
      
      real(8)									:: hbar
      real(8), dimension(n_nm),					intent(in)	:: nm_e_inv ! Transition energies of normal modes and inverse
      real(8), dimension(n_nm, 4),				intent(in)	:: lambda_a  ! Combination terms for Bishop/Luis/Kirtman terms
      real(8), dimension(n_nm, n_nm, 4),			intent(in)	:: lambda_ab ! Combination terms for Bishop/Luis/Kirtman terms

      complex(8), dimension(n_nm, 3), 				intent(in) 	:: dm_1d ! 1st normal mode derivative of dipole moment
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st            "              1st hyperpolarizability
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(n_nm, n_nm, n_nm, n_nm),		intent(in) 	:: fcq   ! Quartic force constants
      complex(8)									:: pv_diphyp_02 

      permopt = make_perm(4, fact(4), (/1, 2, 3, 4/))

      pv_diphyp_02 = 0.0


    hbar = 1.0

    do i=1, n_nm
      do j=1, n_nm
	do k=1, n_nm
	  do m=1, 24
	
	    pv_diphyp_02 = pv_diphyp_02 + nm_e_inv(i) * fcq(i, i, j, k) * dm_1d(j, perm(m, 1)) * &
	    hp_1d(k, perm(m, 2), perm(m, 3), perm(m, 4)) * lambda_a(j, permopt(m, 1)) * &
	    lambda_a(k, permopt(m, 1))
	      
	      do l=1, n_nm

		pv_diphyp_02 = pv_diphyp_02 - fcc(i, i, j) * fcc(j, k, l) * dm_1d(k, perm(m, 1)) * &
		hp_1d(l, perm(m, 2), perm(m, 3), perm(m, 4)) * (nm_e_inv(j)**(2.0)) * &
		lambda_a(k, permopt(m, 1)) * &
		lambda_a(l, permopt(m, 1))  + 2.0 * fcc(i, j, k) * fcc(i, j, l) * dm_1d(k, perm(m, 1)) * &
		hp_1d(l, perm(m, 2), perm(m, 3), perm(m, 4)) * lambda_ab(i, j, permopt(m, 1)) * & 
		lambda_a(k, permopt(m, 1)) * &
		lambda_a(l, permopt(m, 1))

	      end do
	   
	  end do
	end do
      end do
    end do

    pv_diphyp_02 = -hbar * (1.0/24.0) * pv_diphyp_02
        
  end function


!-------------------------------------------------------------------------------
!
! Contains various functions used in the calculation of the pure vibrational
! contribution to the molecular polarizability and 1st/2nd hyperpolarizability.
! The functions in this file is called by the various functions calculating
! the PV quantities.
!
!-------------------------------------------------------------------------------

  function make_lambda_a(n_o, n_nm, opt_e, nm_e)

    ! Creates a combination term of incident energies and vibrational excitation energies

    integer, intent(in):: n_o, n_nm
    integer :: i, k
    real(8), dimension(n_o), intent(in) :: opt_e
    real(8), dimension(n_nm), intent(in) :: nm_e
    real(8), dimension(n_nm, n_o) :: make_lambda_a

    do i = 1, n_nm
      do k = 1, n_o

	make_lambda_a(i, k) = 1.0/((nm_e(i))**(2.0) - (opt_e(k))**(2.0))

      end do
    end do

  end function


  function make_lambda_ab(n_o, n_nm, opt_e, nm_e)

    ! Creates a combination term of incident energies and vibrational excitation energies

    integer, intent(in):: n_o, n_nm
    integer :: i, j, k
    real(8), dimension(n_o), intent(in) :: opt_e
    real(8), dimension(n_nm), intent(in) :: nm_e
    real(8), dimension(n_nm, n_nm, n_o) :: make_lambda_ab

    do i = 1, n_nm
      do j = 1, n_nm
	do k = 1, n_o

	  make_lambda_ab(i, j, k) = 1.0/((nm_e(i) + nm_e(j))**(2.0) - (opt_e(k))**(2.0))

	end do
      end do
    end do

  end function


  function make_lambda_aij(n_o, n_nm, opt_e, nm_e)

    ! Creates a combination term of incident energies and vibrational excitation energies

    integer, intent(in):: n_o, n_nm
    integer :: i, k, l
    real(8), dimension(n_o), intent(in) :: opt_e
    real(8), dimension(n_nm), intent(in) :: nm_e
    real(8), dimension(n_nm, n_o, n_o) :: make_lambda_aij

    do i = 1, n_nm
      do k = 1, n_o
	do l = 1, n_o

	  make_lambda_aij(i, k, l) = 1.0/((nm_e(i))**(2.0) - (opt_e(k) + opt_e(l))**(2.0))

	end do
      end do
    end do

  end function

  function make_lambda_abij(n_o, n_nm, opt_e, nm_e)

    ! Creates a combination term of incident energies and vibrational excitation energies

    integer, intent(in):: n_o, n_nm
    integer :: i, j, k, l
    real(8), dimension(n_o), intent(in) :: opt_e
    real(8), dimension(n_nm), intent(in) :: nm_e
    real(8), dimension(n_nm, n_nm, n_o, n_o) :: make_lambda_abij

    do i = 1, n_nm
      do j = 1, n_nm
	do k = 1, n_o
	  do l = 1, n_o

	    make_lambda_abij(i, j, k, l) = 1.0 / & 
	    ((nm_e(i) + nm_e(j))**(2.0) - (opt_e(k) + opt_e(l))**(2.0))
	
	  end do
	end do
      end do
    end do

  end function


  function invvec(n, vec)

    ! Takes a vector and creates the inverse element by element.
  
    implicit none

    integer :: n, i
    real(8), dimension(n), intent(in) :: vec
    real(8), dimension(n) :: invvec

    do i = 1, n
      
      if (vec(i) == 0.0) then
  
	invvec(i) = 0.0

      else
	
	invvec(i) = 1.0/vec(i)

      end if

    end do

  end function



  function alpha_zpva(n_nm, nm_e, opt_e, po_1d, po_2d, fcc)
!--------------------------------------------------------------------------------
!
! Calculates and returns the ZPVA correction contribution to the molecular
! polarizability tensor from various derivatives taken with 
! respect to the normal coordinates of the molecule. Derivatives above 
! first order are taken as optional, should their implementation be 
! unavailable. The ZPVA correction is calculated only to first order in both
! electrical and mechanical anharmonicity. The next nonvanishing order is the 
! third order, requiring fourth derivatives of electrical properties and quintic
! force constants.
!
! Functions referred to as "Quinet/Kirtman/Champagne" terms reproduce terms 
! presented in Journal of Chemical Physics, Vol. 118, No. 2 (2003), 505-513.
! An adapted form of these terms is used here.
! 
!--------------------------------------------------------------------------------

    implicit none
      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i, j, p, q
      real(8)									:: hbar, fsum
      real(8), dimension(2)							:: opt_e
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      complex(8), dimension(n_nm, 3, 3), 				intent(in) 	:: po_1d ! 1st normal mode derivative of molecular polarizability
      complex(8), dimension(n_nm, n_nm, 3, 3),			intent(in) 	:: po_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm),			intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(3, 3)							:: alpha_zpva


    hbar = 1.0

    ! Precalculate terms for use in Quinet/Kirtman/Champagne term functions

    nm_e_inv = invvec(n_nm, nm_e)

    ! Loop over tensor elements and normal modes

    do i = 1, 3
      do j = 1, 3
	do p = 1, n_nm

	  ! Add electrical anharmonicity term
	  
	  alpha_zpva(i, j) = alpha_zpva(i, j) + po_2d(p, p, i, j) * nm_e_inv(p)

	  ! Prepare a contribution to the mechanical anharmonicity term

	  fsum = 0.0

	  do q = 1, n_nm

	    fsum = fsum - fcc(p, q, q) * nm_e_inv(q)

	  end do

	  ! Add mechanical anharmonicity term

	  alpha_zpva(i, j) = alpha_zpva(i, j) + po_1d(p, i, j) * ((nm_e_inv(p))**2.0)

	end do
	
	alpha_zpva(i, j) = hbar * 0.25 * alpha_zpva(i, j)

      end do
    end do


  end function


  function beta_zpva(n_nm, nm_e, opt_e, hp_1d, hp_2d, fcc)
!--------------------------------------------------------------------------------
!
! Calculates and returns the ZPVA correction contribution to the molecular
! 1st hyperpolarizability tensor.
!
! Functions referred to as "Quinet/Kirtman/Champagne" terms reproduce terms 
! presented in Journal of Chemical Physics, Vol. 118, No. 2 (2003), 505-513.
! An adapted form of these terms is used here.
! 
!--------------------------------------------------------------------------------

    implicit none
      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i, j, k, p, q
      real(8)									:: hbar, fsum
      real(8), dimension(2)							:: opt_e
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      complex(8), dimension(n_nm, 3, 3, 3),			intent(in) 	:: hp_1d ! 1st normal mode derivative of mol. 1st hyperpolarizability
      complex(8), dimension(n_nm, n_nm, 3, 3, 3),  		intent(in) 	:: hp_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm), 	 		intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(3, 3, 3)						:: beta_zpva


    hbar = 1.0

    ! Precalculate terms for use in Quinet/Kirtman/Champagne term functions

    nm_e_inv = invvec(n_nm, nm_e)

    ! Loop over tensor elements and normal modes

    do i = 1, 3
      do j = 1, 3
	do k = 1, 3
	  do p = 1, n_nm

	    ! Add electrical anharmonicity term
	  
	    beta_zpva(i, j, k) = beta_zpva(i, j, k) + hp_2d(p, p, i, j, k) * nm_e_inv(p)

	    ! Prepare a contribution to the mechanical anharmonicity term

	    fsum = 0.0

	    do q = 1, n_nm

	      fsum = fsum - fcc(p, q, q) * nm_e_inv(q)

	    end do

	    ! Add mechanical anharmonicity term

	    beta_zpva(i, j, k) = beta_zpva(i, j, k) + hp_1d(p, i, j, k) * ((nm_e_inv(p))**2.0)

	  end do
	
	  beta_zpva(i, j, k) = hbar * 0.25 * beta_zpva(i, j, k)

	end do
      end do
    end do


  end function


  function gamma_zpva(n_nm, nm_e, opt_e, hp2_1d, hp2_2d, fcc)
!--------------------------------------------------------------------------------
!
! Calculates and returns the ZPVA correction contribution to the molecular
! 1st hyperpolarizability tensor.
!
! Functions referred to as "Quinet/Kirtman/Champagne" terms reproduce terms 
! presented in Journal of Chemical Physics, Vol. 118, No. 2 (2003), 505-513.
! An adapted form of these terms is used here.
! 
!--------------------------------------------------------------------------------

    implicit none
      integer,							intent(in) 	:: n_nm ! Number of normal modes
      integer									:: i, j, k, l, p, q
      real(8)									:: hbar, fsum
      real(8), dimension(2)							:: opt_e
      real(8), dimension(n_nm)							:: nm_e, nm_e_inv ! Transition energies of normal modes and inverse
      complex(8), dimension(n_nm, 3, 3, 3, 3),			intent(in) 	:: hp2_1d ! 1st normal mode derivative of mol. 2nd hyperpolarizability
      complex(8), dimension(n_nm, n_nm, 3, 3, 3, 3),  		intent(in) 	:: hp2_2d ! 2nd                   "
      complex(8), dimension(n_nm, n_nm, n_nm), 	 		intent(in) 	:: fcc   ! Cubic force constants
      complex(8), dimension(3, 3, 3, 3)						:: gamma_zpva


    hbar = 1.0

    ! Precalculate terms for use in Quinet/Kirtman/Champagne term functions

    nm_e_inv = invvec(n_nm, nm_e)

    ! Loop over tensor elements and normal modes

    do i = 1, 3
      do j = 1, 3
	do k = 1, 3
	  do l = 1, 3
	    do p = 1, n_nm

	      ! Add electrical anharmonicity term
	  
	      gamma_zpva(i, j, k, l) = gamma_zpva(i, j, k, l) + &
	      hp2_2d(p, p, i, j, k, l) * nm_e_inv(p)

	      ! Prepare a contribution to the mechanical anharmonicity term

	      fsum = 0.0

	      do q = 1, n_nm

		fsum = fsum - fcc(p, q, q) * nm_e_inv(q)

	      end do

	      ! Add mechanical anharmonicity term

	      gamma_zpva(i, j, k, l) = gamma_zpva(i, j, k, l) + &
	      hp2_1d(p, i, j, k, l) * ((nm_e_inv(p))**2.0)

	    end do
	
	    gamma_zpva(i, j, k, l) = hbar * 0.25 * gamma_zpva(i, j, k, l)

	  end do
	end do
      end do
    end do


  end function


  function hpsq_ravg_vv(hp_1d_tnm)
!--------------------------------------------------------------------------------
!
! Calculates and returns the combination of hyperpolarizability derivatives
! associated with the VV (vertically-vertically plane polarized incident light) 
! hyper-Raman intensity.
!
! This expression is found in J. Chem. Phys. 124, 244312 (2006) by Quinet, O.,
! Champagne, B., and Rodriguez, V.
!
!--------------------------------------------------------------------------------

    implicit none
      integer :: i, j, k
      double precision :: hpsq_ravg_vv
      double precision, dimension(3,3,3) :: hp_1d_tnm ! hp_1d_tnm is hyperpolarizability derivatives for 'this normal mode'

    ! Loop over Cartesian axes and make the appropriate combinations

    hpsq_ravg_vv = 0.0

    do i = 1, 3

      hpsq_ravg_vv = hpsq_ravg_vv + 1.0/7.0 * (( hp_1d_tnm(i,i,i) )**(2.0))

      do j = 1, 3

	if (.NOT.(j == i)) then

	  hpsq_ravg_vv = hpsq_ravg_vv + 1.0/35.0 * (4.0 * (( hp_1d_tnm(i,i,j) )**(2.0)) + &
	  2.0 * hp_1d_tnm(i,i,i) * hp_1d_tnm(i,j,j) + 4.0 * hp_1d_tnm(j,i,i) * hp_1d_tnm(i,i,j) + &
	  4.0 * hp_1d_tnm(i,i,i) * hp_1d_tnm(j,j,i) + (( hp_1d_tnm(j,i,i) )**(2.0)) )
	  
	  do k = 1, 3

	    if (.NOT.((k == i) .OR. (k == j))) then

	    hpsq_ravg_vv = hpsq_ravg_vv + 1.0/105.0 * (4.0 * hp_1d_tnm(i,i,j) * hp_1d_tnm(j,k,k) + &
	    hp_1d_tnm(j,i,i) * hp_1d_tnm(j,k,k) + 4.0 * hp_1d_tnm(i,i,j) * hp_1d_tnm(k,k,j) + &
	    2.0 * (( hp_1d_tnm(i,j,k) )**(2.0)) + 4.0 * hp_1d_tnm(i,j,k) * hp_1d_tnm(j,i,k))

	    end if

	  end do

	end if
 
      end do
    end do

  end function

  function hpsq_ravg_hv(hp_1d_tnm)
!--------------------------------------------------------------------------------
!
! Calculates and returns the combination of hyperpolarizability derivatives
! associated with the HV (horizontally-vertically plane polarized incident light)
! hyper-Raman intensity.
!
! This expression is found in J. Chem. Phys. 124, 244312 (2006) by Quinet, O.,
! Champagne, B., and Rodriguez, V.
!
!--------------------------------------------------------------------------------

    implicit none

      integer :: i, j, k
      double precision :: hpsq_ravg_hv
      double precision, dimension(3,3,3) :: hp_1d_tnm ! hp_1d_tnm is hyperpolarizability derivatives for 'this normal mode'

    ! Loop over Cartesian axes and make the appropriate combinations

    hpsq_ravg_hv = 0.0

    do i = 1, 3

      hpsq_ravg_hv = hpsq_ravg_hv + 1.0/35.0 * (( hp_1d_tnm(i,i,i) )**(2.0))

      do j = 1, 3

	if (.NOT.(j == i)) then

	  hpsq_ravg_hv = hpsq_ravg_hv + 4.0/105.0 * hp_1d_tnm(i,i,i) * hp_1d_tnm(i,j,j) - &
	  2.0/35.0 * hp_1d_tnm(i,i,i) * hp_1d_tnm(j,j,i) + 8.0/105.0 * (( hp_1d_tnm(i,i,j) )**(2.0)) + &
	  3.0/35.0 * (( hp_1d_tnm(i,j,j) )**(2.0)) - 2.0/35.0 * hp_1d_tnm(i,i,j) * hp_1d_tnm(j,i,i)
	  	  
	  do k = 1, 3

	    if (.NOT.((k == i) .OR. (k == j))) then

	    hpsq_ravg_hv = hpsq_ravg_hv + 1.0/35.0 * hp_1d_tnm(i,j,j) * hp_1d_tnm(i,k,k) - &
	    2.0/105.0 * hp_1d_tnm(i,i,k) * hp_1d_tnm(j,j,k) - &
	    2.0/105.0 * hp_1d_tnm(i,i,j) * hp_1d_tnm(j,k,k) + 2.0/35.0 * (( hp_1d_tnm(i,j,k) )**(2.0)) - &
	    2.0/105.0 * hp_1d_tnm(i,j,k) * hp_1d_tnm(j,i,k)

	    end if

	  end do

	end if
 
      end do
    end do

  end function

  function hyp_raman_vv(n_nm, temp, opt_e, nm_e, hp_1d)
!--------------------------------------------------------------------------------
!
! Calculates and returns the hyper-Raman intensity associated with 
! the VV (vertically-vertically plane polarized incident light) scenario.
!
! This expression is adapted from one found in J. Chem. Phys. 124, 244312 (2006)
! by Quinet, O., Champagne, B., and Rodriguez, V.
!
!--------------------------------------------------------------------------------

    implicit none

      integer :: n_nm, i
      double precision :: hbar, boltz, temp, nm_e ! H-bar, Boltzmann's constant and temperature
      double precision :: opt_e, hyp_raman_vv
      double precision, dimension(3,3,3) :: hp_1d

    hbar = 1.0
    boltz = 0.0000031668153
    hyp_raman_vv = 0.0

    if (temp == 0.0)  then
    
      hyp_raman_vv = hyp_raman_vv + hbar * ((2.0 * opt_e - nm_e)**(4.0)) / &
      (2.0 * nm_e ) * hpsq_ravg_vv(hp_1d)

    else

      hyp_raman_vv = hyp_raman_vv + hbar * ((2.0 * opt_e - nm_e)**(4.0)) / &
      ((2.0 * nm_e) * (1.0 - exp( - 1.0 * hbar * nm_e / (boltz * temp) ))) * &
      hpsq_ravg_vv(hp_1d)

    end if



  end function


  function hyp_raman_hv(n_nm, temp, opt_e, nm_e, hp_1d)
!--------------------------------------------------------------------------------
!
! Calculates and returns the hyper-Raman intensity associated with 
! the HV (horizontally-vertically plane polarized incident light) scenario.
!
! This expression is adapted from one found in J. Chem. Phys. 124, 244312 (2006)
! by Quinet, O., Champagne, B., and Rodriguez, V.
!
!--------------------------------------------------------------------------------

    implicit none

      integer :: n_nm, i
      double precision :: hbar, boltz, temp, nm_e ! H-bar, Boltzmann's constant and temperature
      double precision :: opt_e, hyp_raman_hv
      double precision, dimension(3,3,3) :: hp_1d

    hbar = 1.0
    boltz = 0.0000031668153
    hyp_raman_hv = 0.0

    if (temp == 0.0)  then

      hyp_raman_hv = hyp_raman_hv + hbar * ((2.0 * opt_e - nm_e)**(4.0)) / &
      (2.0 * nm_e) * hpsq_ravg_hv(hp_1d)

    else

      hyp_raman_hv = hyp_raman_hv + hbar * ((2.0 * opt_e - nm_e)**(4.0)) / &
      ((2.0 * nm_e) * (1.0 - exp( - 1.0 * hbar * nm_e / (boltz * temp) ))) * &
      hpsq_ravg_hv(hp_1d)

    end if

  end function



end module vib_pv_contribs