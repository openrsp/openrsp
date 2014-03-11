module interface_xc

!  interface_xc_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_xc_init

   use matrix_defop, matrix => openrsp_matrix
   use interface_molecule
   use rsp_field_tuple
   use rsp_indices_and_addressing
   use rsp_sdf_caching
   use iso_c_binding

   implicit none

   public interface_xc_init
   public interface_xc_finalize
   public xcint_wakeup_workers

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT

   public rsp_xcave_interface
   public rsp_xcint_interface

   public get_is_ks_calculation
   public openrsp_set_functional

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_ks_calculation

#ifdef VAR_MPI
#include "mpif.h"
#endif

interface xcint_init
   subroutine xcint_init(basis_type,              &
                         nr_centers,              &
                         center_xyz,              &
                         center_element,          &
                         nr_shells,               &
                         shell_center,            &
                         l_quantum_nr,            &
                         nr_primitives_per_shell, &
                         primitive_exp,           &
                         contraction_coef,        &
                         radial_precision,        &
                         angular_order) bind (C, name = "xcint_init")

      use iso_c_binding
      implicit none

      integer(c_int), value :: basis_type
      integer(c_int), value :: nr_centers
      real(c_double)        :: center_xyz(*)
      integer(c_int)        :: center_element(*)
      integer(c_int), value :: nr_shells
      integer(c_int)        :: shell_center(*)
      integer(c_int)        :: l_quantum_nr(*)
      integer(c_int)        :: nr_primitives_per_shell(*)
      real(c_double)        :: primitive_exp(*)
      real(c_double)        :: contraction_coef(*)
      real(c_double), value :: radial_precision
      integer(c_int), value :: angular_order
   end subroutine
end interface

interface xcint_set_functional
   subroutine xcint_set_functional(line, hfx, mu, beta) bind (C, name = "xcint_set_functional")
      use iso_c_binding
      implicit none
      character(c_char) :: line
      real(c_double)    :: hfx
      real(c_double)    :: mu
      real(c_double)    :: beta
   end subroutine
end interface

interface xcint_integrate
   subroutine xcint_integrate(num_dmat,        &
                              dmat,            &
                              fmat,            &
                              energy,          &
                              get_ave,         &
                              geo_derv_order,  &
                              geo_coor,        &
                              num_fields,      &
                              kn_rule,         &
                              force_sequential &
                             ) bind (C, name = "xcint_integrate")
      use iso_c_binding
      implicit none
      integer(c_int), value :: num_dmat
      real(c_double)        :: dmat(*)
      real(c_double)        :: fmat(*)
      real(c_double)        :: energy
      integer(c_int), value :: get_ave
      integer(c_int), value :: geo_derv_order
      integer(c_int)        :: geo_coor(*)
      integer(c_int), value :: num_fields
      integer(c_int)        :: kn_rule(2)
      integer(c_int), value :: force_sequential
   end subroutine
end interface

contains

   subroutine xcint_wakeup_workers()
      integer :: ierr, num_proc
      integer :: iprint = 0
#ifdef VAR_MPI
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"
      call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)
      if (num_proc > 1) then
         CALL MPIXBCAST(XCINT_MPI_WAKEUP_SIGNAL, 1,'INTEGER', MASTER)
         CALL MPIXBCAST(iprint, 1,'INTEGER', MASTER)
      end if
#endif
   end subroutine

   subroutine interface_xc_init()

#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "shells.h"
#include "primit.h"

      integer(c_int), allocatable :: nr_primitives_per_shell(:)
      real(c_double), allocatable :: primitive_exp(:, :)
      real(c_double), allocatable :: contraction_coef(:, :)
      integer(c_int), allocatable :: l_quantum_nr(:)
      integer(c_int), allocatable :: shell_center(:)
      integer                     :: i, j, icount, iprim, ishell, iround
      integer(c_int)              :: basis_type
      integer(c_int)              :: nr_shells
      integer(c_int)              :: nr_centers
      real(c_double), allocatable :: center_xyz(:, :)
      integer(c_int), allocatable :: center_element(:)

      if (is_initialized) return

      nr_shells = kmax
      nr_centers = nucind
      allocate(center_xyz(3, nr_centers))
      center_xyz = cord

      allocate(center_element(nr_centers))
      do i = 1, nr_centers
         center_element(i) = nint(charge(i))
      end do

      allocate(nr_primitives_per_shell(nr_shells))
      do iround = 1, 2
         if (iround == 2) then
            allocate(primitive_exp(maxval(nr_primitives_per_shell), nr_shells))
            allocate(contraction_coef(maxval(nr_primitives_per_shell), nr_shells))
         end if
         do ishell = 1, nr_shells
            i = jstrt(ishell) + 1
            j = jstrt(ishell) + nuco(ishell)
            icount = 0
            do iprim = i, j
               if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
                  icount = icount + 1
                  nr_primitives_per_shell(ishell) = icount
                  if (iround == 2) then
                     contraction_coef(icount, ishell) = priccf(iprim, numcf(ishell))
                     primitive_exp(icount, ishell)    = priexp(iprim)
                  end if
               end if
            end do
         end do
      end do

      allocate(l_quantum_nr(nr_shells))
      allocate(shell_center(nr_shells))
      do ishell = 1, nr_shells
         l_quantum_nr(ishell) = nhkt(ishell) - 1
         shell_center(ishell) = ncent(ishell)
      end do

      ! we default to spherical basis
      basis_type = 1
      do ishell = 1, nr_shells
         if ((l_quantum_nr(ishell) > 1) .and. .not. sphr(ishell)) then
            ! basis is cartesian
            basis_type = 2
            exit
         end if
      end do

#ifdef VAR_MPI
      call xcint_set_mpi_comm(MPI_COMM_WORLD)
#endif

      call xcint_init(basis_type,                                            &
                      nr_centers,                                            &
                      reshape(center_xyz, (/size(center_xyz)/)),             &
                      center_element,                                        &
                      nr_shells,                                             &
                      shell_center,                                          &
                      l_quantum_nr,                                          &
                      nr_primitives_per_shell,                               &
                      reshape(primitive_exp, (/size(primitive_exp)/)),       &
                      reshape(contraction_coef, (/size(contraction_coef)/)), &
                      1.0d-13,                                               &
                      17);

      deallocate(center_xyz)
      deallocate(center_element)
      deallocate(nr_primitives_per_shell)
      deallocate(primitive_exp)
      deallocate(contraction_coef)
      deallocate(l_quantum_nr)
      deallocate(shell_center)

      is_initialized = .true.

   end subroutine

   subroutine openrsp_set_functional(line, hfx, mu, beta)
      character(*),   intent(in)  :: line
      real(c_double), intent(out) :: hfx
      real(c_double), intent(out) :: mu
      real(c_double), intent(out) :: beta
      call xcint_set_functional(line//C_NULL_CHAR, hfx, mu, beta)
      is_ks_calculation = .true.
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





   !> Exchange-correlation perturbed by fields f, averaged over densities D
   !> New routine under development (by MaR)
   subroutine rsp_xcave_interface(mat_dim,       &
                        pert,          &
                        kn,            &
                        num_blks,      &
                        blk_sizes,     &
                        blk_info,      &
                        property_size, &
                        prop,          &
                        D_sdf)

!     --------------------------------------------------------------------------
      integer,       intent(in)    :: mat_dim
      integer, allocatable, dimension(:) :: emptyint
      type(p_tuple), intent(in)    :: pert
      integer,       intent(in)    :: kn(2)
      integer,       intent(in)    :: num_blks
      integer,       intent(in)    :: blk_sizes(num_blks)
      integer,       intent(in)    :: blk_info(num_blks, 3)
      integer,       intent(in)    :: property_size
      complex(8),    intent(inout) :: prop(property_size)
      type(SDF)                    :: D_sdf
!     --------------------------------------------------------------------------
      complex(8)                   :: res(property_size) !fixme, allocate
      integer                      :: i, j, k, l, m, n, p
      integer                      :: maxcomp1, maxcomp2, maxcomp3, maxcomp4, maxcomp5, maxcomp6
      integer                      :: element
      integer                      :: nr_dmat
      integer                      :: nr_atoms
      integer                      :: nr_geo
      integer                      :: nr_el
      integer                      :: ipart
      logical                      :: combination_found
      type(matrix),   allocatable  :: dmat_tuple(:)
      real(c_double)               :: xc_energy
      integer(c_int)              :: get_ave
      integer(c_int)              :: force_sequential
!     --------------------------------------------------------------------------

      get_ave = 1
      force_sequential = 0

      res = 0.0

      if (.not. get_is_ks_calculation()) return

      nr_atoms = get_nr_atoms()

      nr_geo = count(pert%plab == 'GEO ')
      nr_el  = count(pert%plab == 'EL  ')


      combination_found = .false.

      nr_dmat = 2**(nr_geo + nr_el - 1)

      allocate(dmat_tuple(nr_dmat))
      
      allocate(emptyint(0))

      ! MaR: Begin EL only cases
      
      if (nr_geo == 0 .and. nr_el == 1) then
               
      ! No contribution
      
      end if      
      
      if (nr_geo == 0 .and. nr_el == 2) then
         combination_found = .true.
      
      ! No contribution      
      
      end if      
      
      if (nr_geo == 0 .and. nr_el == 3) then
         combination_found = .true.
      
      ! Assumes (k,n) = (0,2)
      ! Then no contribution
            
      end if
      
      if (nr_geo == 0 .and. nr_el == 4) then
         combination_found = .true.
      
      ! Assumes (k,n) = (0,3)
      ! Then no contribution
            
      end if
      
      if (nr_geo == 0 .and. nr_el == 5) then
         combination_found = .true.

         ! Assumes (k,n) = (2,2)


         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))
         
                        
         do i = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            if (p_tuple_compare(p_tuple_getone(pert, 1) , &
                p_tuple_getone(pert, 2))) then
               maxcomp2 = i
            else
               maxcomp2 = 3
            end if

            do j = 1, maxcomp2

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                  p_tuple_getone(pert, 2)), (/i,j/), dmat_tuple(4))

               if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                   p_tuple_getone(pert, 3))) then
                  maxcomp3 = j
               else
                  maxcomp3 = 3
               end if

               do k = 1, maxcomp3

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(5))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                     p_tuple_getone(pert, 3)), (/i,k/), dmat_tuple(6))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(7))

!                   call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
!                                      merge_p_tuple(p_tuple_getone(pert, 3),  &
!                                      p_tuple_getone(pert, 4))), (/i,j,k/), &
!                                      dmat_tuple(8))

                  if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                     p_tuple_getone(pert, 4))) then
                     maxcomp4 = k
                  else
                     maxcomp4 = 3
                  end if

                  do l = 1, maxcomp4
      
                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(8))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                        p_tuple_getone(pert, 4)), (/i,l/), dmat_tuple(9))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(10))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(11))

                     if (p_tuple_compare(p_tuple_getone(pert, 4) , &
                        p_tuple_getone(pert, 5))) then
                        maxcomp5 = l
                     else
                        maxcomp5 = 3
                     end if


                     do m = 1, maxcomp5

                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/m/), dmat_tuple(12))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 5)), (/i,m/), dmat_tuple(13))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 5)), (/j,m/), dmat_tuple(14))
                     
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 5)), (/k,m/), dmat_tuple(15))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5)), (/l,m/), dmat_tuple(16))                                                                
                     
                     
                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/i, j, k, l, m/))

                        call xcint_wakeup_workers()
                        call xcint_integrate(nr_dmat,                      &
                                             (/dmat_tuple(1)%elms,   &
                                               dmat_tuple(2)%elms,   &
                                               dmat_tuple(3)%elms,   &
                                               dmat_tuple(4)%elms,   &
                                               dmat_tuple(5)%elms,   &
                                               dmat_tuple(6)%elms,   &
                                               dmat_tuple(7)%elms,   &
                                               dmat_tuple(8)%elms,   &
                                               dmat_tuple(9)%elms,   &
                                               dmat_tuple(10)%elms,   &
                                               dmat_tuple(11)%elms,   &
                                               dmat_tuple(12)%elms,   &
                                               dmat_tuple(13)%elms,   &
                                               dmat_tuple(14)%elms,   &
                                               dmat_tuple(15)%elms,   &
                                               dmat_tuple(16)%elms/), &
                                               (/0.0d0/),              &
                                               xc_energy,              &
                                               get_ave,                &
                                               0,                      &
                                               emptyint,                   &
                                               5,                      &
                                               kn,                     &
                                               force_sequential)

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do
            
      end if
      
      
      if (nr_geo == 0 .and. nr_el == 6) then
         combination_found = .true.
      
         ! Assumes (k,n) = (2,3)


         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))
         
                        
         do i = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            if (p_tuple_compare(p_tuple_getone(pert, 1) , &
                p_tuple_getone(pert, 2))) then
               maxcomp2 = i
            else
               maxcomp2 = 3
            end if

            do j = 1, maxcomp2

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                  p_tuple_getone(pert, 2)), (/i,j/), dmat_tuple(4))

               if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                   p_tuple_getone(pert, 3))) then
                  maxcomp3 = j
               else
                  maxcomp3 = 3
               end if

               do k = 1, maxcomp3

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(5))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                     p_tuple_getone(pert, 3)), (/i,k/), dmat_tuple(6))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(7))

!                   call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
!                                      merge_p_tuple(p_tuple_getone(pert, 3),  &
!                                      p_tuple_getone(pert, 4))), (/i,j,k/), &
!                                      dmat_tuple(8))

                  if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                     p_tuple_getone(pert, 4))) then
                     maxcomp4 = k
                  else
                     maxcomp4 = 3
                  end if

                  do l = 1, maxcomp4
      
                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(8))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                        p_tuple_getone(pert, 4)), (/i,l/), dmat_tuple(9))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(10))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(11))
                                        
                                        
                                        
                                        
                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 4))), (/j,k,l/), &
                                        dmat_tuple(12))

                     if (p_tuple_compare(p_tuple_getone(pert, 4) , &
                        p_tuple_getone(pert, 5))) then
                        maxcomp5 = l
                     else
                        maxcomp5 = 3
                     end if


                     do m = 1, maxcomp5

                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/m/), dmat_tuple(13))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 5)), (/i,m/), dmat_tuple(14))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 5)), (/j,m/), dmat_tuple(15))
                     
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 5)), (/k,m/), dmat_tuple(16))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5)), (/l,m/), dmat_tuple(17))                                                                


                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 5))), (/j,k,m/), &
                                           dmat_tuple(18))                         
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5))), (/j,l,m/), &
                                           dmat_tuple(19))                                                                    
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5))), (/k,l,m/), &
                                           dmat_tuple(20))
                                           
                                           
                     if (p_tuple_compare(p_tuple_getone(pert, 5) , &
                        p_tuple_getone(pert, 6))) then
                        maxcomp6 = m
                     else
                        maxcomp6 = 3
                     end if


                     do n = 1, maxcomp5                                                                                                           
                     
                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 6), (/n/), dmat_tuple(21))
                     
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 6)), (/i,n/), dmat_tuple(22))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 6)), (/j,n/), dmat_tuple(23))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 6)), (/k,n/), dmat_tuple(24))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 6)), (/l,n/), dmat_tuple(25))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 5),  &
                                           p_tuple_getone(pert, 6)), (/m,n/), dmat_tuple(26))                                           
                     
                     
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 6))), (/j,k,n/), &
                                           dmat_tuple(27))
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 6))), (/j,l,n/), &
                                           dmat_tuple(28))        
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 6))), (/k,l,n/), &
                                           dmat_tuple(29))        
                                           
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           merge_p_tuple(p_tuple_getone(pert, 5),  &
                                           p_tuple_getone(pert, 6))), (/j,m,n/), &
                                           dmat_tuple(30))        
                                           
                       call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           merge_p_tuple(p_tuple_getone(pert, 5),  &
                                           p_tuple_getone(pert, 6))), (/k,m,n/), &
                                           dmat_tuple(31))        
                                           
                       call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           merge_p_tuple(p_tuple_getone(pert, 5),  &
                                           p_tuple_getone(pert, 6))), (/l,m,n/), &
                                           dmat_tuple(32))        
                     
                     
                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/i, j, k, l, m, n/))

                        call xcint_wakeup_workers()
                        call xcint_integrate(nr_dmat,                      &
                                             (/(dmat_tuple(i)%elms, i = 1, nr_dmat)/), &
                                               (/0.0d0/),              &
                                               xc_energy,              &
                                               get_ave,                &
                                               0,                      &
                                               emptyint,                   &
                                               6,                      &
                                               kn,                     &
                                               force_sequential)

                        res(element) = cmplx(xc_energy, 0.0d0)

                        end do
                     end do
                  end do
               end do
            end do
         end do
                 
      end if

      
      
      ! End EL only cases
      
      if (nr_geo == 1 .and. nr_el == 0) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))
         do i = 1, nr_atoms*3
            element = i
            call xcint_wakeup_workers()
            call xcint_integrate(1,                      &
                                 (/dmat_tuple(1)%elms/), &
                                 (/0.0d0/),              &
                                 xc_energy,              &
                                 get_ave,                &
                                 1,                      &
                                 (/i/),                  &
                                 0,                      &
                                 kn,                     &
                                 force_sequential)
            res(element) = cmplx(xc_energy, 0.0d0)
         end do
      end if

      if (nr_geo == 2 .and. nr_el == 0) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3
            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xcint_wakeup_workers()
               call xcint_integrate(2,                      &
                                    (/dmat_tuple(1)%elms,   &
                                      dmat_tuple(2)%elms/), &
                                    (/0.0d0/),              &
                                    xc_energy,              &
                                    get_ave,                &
                                    2,                      &
                                    (/i, j/),               &
                                    0,                      &
                                    kn,                     &
                                    force_sequential)

               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         end do
      end if

      if (nr_geo == 1 .and. nr_el == 1) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do j = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

            do i = 1, nr_atoms*3

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xcint_wakeup_workers()
               call xcint_integrate(2,                      &
                                    (/dmat_tuple(1)%elms,   &
                                      dmat_tuple(2)%elms/), &
                                    (/0.0d0/),              &
                                    xc_energy,              &
                                    get_ave,                &
                                    1,                      &
                                    (/i/),                  &
                                    1,                      &
                                    kn,                     &
                                    force_sequential)

               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         end do
      end if

      if (nr_geo == 3 .and. nr_el == 0) then
         combination_found = .true.


         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do k = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(4,                      &
                                       (/dmat_tuple(1)%elms,   &
                                         dmat_tuple(2)%elms,   &
                                         dmat_tuple(3)%elms,   &
                                         dmat_tuple(4)%elms/), &
                                       (/0.0d0/),              &
                                       xc_energy,              &
                                       get_ave,                &
                                       3,                      &
                                       (/i, j, k/),            &
                                       0,                      &
                                       kn,                     &
                                       force_sequential)

                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         end do
      end if

      if (nr_geo == 2 .and. nr_el == 1) then
         combination_found = .true.

         ! Set up for (k,n) = (1,1)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do k = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

            do j = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do i = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))


                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(4,                      &
                                       (/dmat_tuple(1)%elms,   &
                                         dmat_tuple(2)%elms,   &
                                         dmat_tuple(3)%elms,   &
                                         dmat_tuple(4)%elms/), &
                                       (/0.0d0/),              &
                                       xc_energy,              &
                                       get_ave,                &
                                       2,                      &
                                       (/i, j/),               &
                                       1,                      &
                                       kn,                     &
                                       force_sequential)
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         end do
      end if

      if (nr_geo == 1 .and. nr_el == 2) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         ipart = 0
         do k = 1, 3

            if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                p_tuple_getone(pert, 3))) then
               maxcomp2 = k
            else
               maxcomp2 = 3
            end if

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/k/), dmat_tuple(2))

            do j = 1, maxcomp2

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/j/), dmat_tuple(3))
            call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                               p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(4))

               do i = 1, nr_atoms*3
                  ipart = ipart + 1
                  print *, 'computing XC ave GFF contribution', ipart, 'out of', 3*maxcomp2*nr_atoms*3

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(4,                      &
                                       (/dmat_tuple(1)%elms,   &
                                         dmat_tuple(2)%elms,   &
                                         dmat_tuple(3)%elms,   &
                                         dmat_tuple(4)%elms/), &
                                       (/0.0d0/),              &
                                       xc_energy,              &
                                       get_ave,                &
                                       1,                      &
                                       (/i/),                  &
                                       2,                      &
                                       kn,                     &
                                       force_sequential)

                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         end do
         print *, 'done computing XC ave GFF contribution'
      end if

      if (nr_geo == 4 .and. nr_el == 0) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do l = 1, nr_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(6))

            do k = 1, l

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                  p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(8))

               do j = 1, k

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))
               
                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(5))
               
                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(7))

                  do i = 1, j

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

                     element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                               blk_info, blk_sizes, (/i, j, k, l/))

                     call xcint_wakeup_workers()
                     call xcint_integrate(8,                      &
                                          (/dmat_tuple(1)%elms,   &
                                            dmat_tuple(2)%elms,   &
                                            dmat_tuple(3)%elms,   &
                                            dmat_tuple(4)%elms,   &
                                            dmat_tuple(5)%elms,   &
                                            dmat_tuple(6)%elms,   &
                                            dmat_tuple(7)%elms,   &
                                            dmat_tuple(8)%elms/), &
                                          (/0.0d0/),              &
                                          xc_energy,              &
                                          get_ave,                &
                                          4,                      &
                                          (/i, j, k, l/),         &
                                          0,                      &
                                          kn,                     &
                                          force_sequential)

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 3 .and. nr_el == 1) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do l = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(6))

            do k = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                  p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(8))

               do j = 1, k

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(5))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(7))

                  do i = 1, j

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

                     element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                               blk_info, blk_sizes, (/i, j, k, l/))

                     call xcint_wakeup_workers()
                     call xcint_integrate(8,                      &
                                          (/dmat_tuple(1)%elms,   &
                                            dmat_tuple(2)%elms,   &
                                            dmat_tuple(3)%elms,   &
                                            dmat_tuple(4)%elms,   &
                                            dmat_tuple(5)%elms,   &
                                            dmat_tuple(6)%elms,   &
                                            dmat_tuple(7)%elms,   &
                                            dmat_tuple(8)%elms/), &
                                          (/0.0d0/),              &
                                          xc_energy,              &
                                          get_ave,                &
                                          3,                      &
                                          (/i, j, k/),            &
                                          1,                      &
                                          kn,                     &
                                          force_sequential)

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 2 .and. nr_el == 2) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do l = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(6))

            if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                p_tuple_getone(pert, 4))) then
               maxcomp2 = l
            else
               maxcomp2 = 3
            end if

            do k = 1, maxcomp2

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                  p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(8))

               do j = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(5))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(7))

                  do i = 1, j

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

                     element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                               blk_info, blk_sizes, (/i, j, k, l/))

                     call xcint_wakeup_workers()
                     call xcint_integrate(8,                      &
                                          (/dmat_tuple(1)%elms,   &
                                            dmat_tuple(2)%elms,   &
                                            dmat_tuple(3)%elms,   &
                                            dmat_tuple(4)%elms,   &
                                            dmat_tuple(5)%elms,   &
                                            dmat_tuple(6)%elms,   &
                                            dmat_tuple(7)%elms,   &
                                            dmat_tuple(8)%elms/), &
                                          (/0.0d0/),              &
                                          xc_energy,              &
                                          get_ave,                &
                                          2,                      &
                                          (/i, j/),               &
                                          2,                      &
                                          kn,                     &
                                          force_sequential)

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 1 .and. nr_el == 3) then
         combination_found = .true.

         ! The array of density matrices should have length 4

         ! Get the unperturbed D and put it in position 1 in dmat_tuple
         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, 3

            ! Get Df w.r.t. the first electric field (perturbation #2 in the GFFF tuple) at index i
            ! and put it in position 2 in dmat_tuple
            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/i/), dmat_tuple(2))

            ! Check if electric field 1 and 2 are equivalent (same frequency)
            ! If so, do "nonredundant loop" for field 2, if not, do all 3 components of both field 1 and 2
            if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                p_tuple_getone(pert, 3))) then
               maxcomp2 = i
            else
               maxcomp2 = 3
            end if

            do j = 1, maxcomp2

               ! Df w.r.t. the second electric field - position 3
               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/j/), dmat_tuple(3))

               ! Dff w.r.t. fields 1 and 2 - position 4
               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 3)), (/i,j/), dmat_tuple(4))

               ! Check if fields 2 and 3 are equivalent and adjust loop extent accordingly
               ! We shouldn't need to check fields 1 and 3 because the perturbations should be sorted already
               ! If not sorted, we would have risked missing some nonredundancy loop extent limitation by
               ! not comparing fields 1 and 3
               if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                   p_tuple_getone(pert, 4))) then
                  maxcomp3 = j
               else
                  maxcomp3 = 3
               end if

               do k = 1, maxcomp3

               ! Df w.r.t. field 3 - position 5
               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/k/), dmat_tuple(5))

               ! Dff w.r.t. fields 1 and 3 - position 6
               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 4)), (/i,k/), dmat_tuple(6))

               ! Dff w.r.t. fields 1 and 2 - position 7
               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                  p_tuple_getone(pert, 4)), (/j,k/), dmat_tuple(7))

               ! Dff w.r.t. fields 1, 2 and 3 - position 8
               call sdf_getdata_s(D_sdf, p_tuple_remove_first(pert), (/i,j,k/), &
                                  dmat_tuple(8))

                  do l = 1, nr_atoms*3

                     ! Get the offset in the response tensor
                     element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                               blk_info, blk_sizes, (/l, i, j, k/))

                     call xcint_wakeup_workers()
                     call xcint_integrate(8,                      &
                                          (/dmat_tuple(1)%elms,   &
                                            dmat_tuple(2)%elms,   &
                                            dmat_tuple(3)%elms,   &
                                            dmat_tuple(4)%elms,   &
                                            dmat_tuple(5)%elms,   &
                                            dmat_tuple(6)%elms,   &
                                            dmat_tuple(7)%elms,   &
                                            dmat_tuple(8)%elms/), &
                                          (/0.0d0/),              &
                                          xc_energy,              &
                                          get_ave,                &
                                          1,                      &
                                          (/l/),                  &
                                          3,                      &
                                          kn,                     &
                                          force_sequential)

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 1 .and. nr_el == 4) then
         combination_found = .true.

         ! Using (k,n) = (0,4)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/i/), dmat_tuple(2))

            if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                p_tuple_getone(pert, 3))) then
               maxcomp2 = i
            else
               maxcomp2 = 3
            end if

            do j = 1, maxcomp2

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                  p_tuple_getone(pert, 3)), (/i,j/), dmat_tuple(4))

               if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                   p_tuple_getone(pert, 4))) then
                  maxcomp3 = j
               else
                  maxcomp3 = 3
               end if

               do k = 1, maxcomp3

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/k/), dmat_tuple(5))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 4)), (/i,k/), dmat_tuple(6))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                     p_tuple_getone(pert, 4)), (/j,k/), dmat_tuple(7))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     merge_p_tuple(p_tuple_getone(pert, 3),  &
                                     p_tuple_getone(pert, 4))), (/i,j,k/), &
                                     dmat_tuple(8))

                  if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                     p_tuple_getone(pert, 4))) then
                     maxcomp4 = k
                  else
                     maxcomp4 = 3
                  end if

                  do l = 1, maxcomp4

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/l/), dmat_tuple(9))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 5)), (/i,l/), dmat_tuple(10))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 5)), (/j,l/), dmat_tuple(11))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 5))), (/i,j,l/), &
                                        dmat_tuple(12))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                        p_tuple_getone(pert, 5)), (/k,l/), dmat_tuple(13))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        merge_p_tuple(p_tuple_getone(pert, 4),  &
                                        p_tuple_getone(pert, 5))), (/i,k,l/), &
                                        dmat_tuple(14))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        merge_p_tuple(p_tuple_getone(pert, 4),  &
                                        p_tuple_getone(pert, 5))), (/j,k,l/), &
                                        dmat_tuple(15))

                     call sdf_getdata_s(D_sdf, p_tuple_remove_first(pert), (/i,j,k,l/), &
                                        dmat_tuple(16))

                     do m = 1, nr_atoms*3

                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/m, i, j, k, l/))

                        print *, 'ERROR: implement g=1 f=4 in openrsp/xcave'
                        stop 1

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 2 .and. nr_el == 3) then
          combination_found = .true.
      
          ! Using (k,n) = (1,3)
      
          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))
      
          do i = 1, nr_atoms*3
      
             call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))
      
             do j = 1, i
      
                call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))
      
                do k = 1, 3
      
                   call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))
      
                   call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                      p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(5))
      
                   if (p_tuple_compare(p_tuple_getone(pert, 3) , &
                       p_tuple_getone(pert, 4))) then
                      maxcomp2 = k
                   else
                      maxcomp2 = 3
                   end if
      
                   do l = 1, maxcomp2
      
                      call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(6))
      
                      call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                         p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(7))
      
                      call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                         p_tuple_getone(pert, 3)), (/k,l/), dmat_tuple(8))
      
                      call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                         merge_p_tuple(p_tuple_getone(pert, 3),  &
                                         p_tuple_getone(pert, 4))), (/j,k,l/), &
                                         dmat_tuple(9))
      
                      if (p_tuple_compare(p_tuple_getone(pert, 4) , &
                         p_tuple_getone(pert, 5))) then
                         maxcomp3 = l
                      else
                         maxcomp3 = 3
                      end if
      
                      do m = 1, maxcomp3
      
                         call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/m/), dmat_tuple(10))
      
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                            p_tuple_getone(pert, 5)), (/j,m/), dmat_tuple(11))
      
                         ! radovan: FIXME note to myself to consider changing dmats 12 and 13
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                            p_tuple_getone(pert, 5)), (/k,m/), dmat_tuple(12))
      
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                            p_tuple_getone(pert, 5)), (/l,m/), dmat_tuple(13))
      
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                            merge_p_tuple(p_tuple_getone(pert, 3),  &
                                            p_tuple_getone(pert, 5))), (/j,k,m/), &
                                            dmat_tuple(14))
      
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                            merge_p_tuple(p_tuple_getone(pert, 4),  &
                                            p_tuple_getone(pert, 5))), (/j,l,m/), &
                                            dmat_tuple(15))
      
                         call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                            merge_p_tuple(p_tuple_getone(pert, 4),  &
                                            p_tuple_getone(pert, 5))), (/k,l,m/), &
                                            dmat_tuple(16))
      
                         element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                   blk_info, blk_sizes, (/i, j, k, l, m/))
      
                         print *, 'ERROR: implement g=2 f=3 in openrsp/xcave'
                         stop 1
      
                         res(element) = cmplx(xc_energy, 0.0d0)
      
                      end do
                   end do
                end do
             end do
          end do
      end if

      if (nr_geo == 3 .and. nr_el == 2) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                  p_tuple_getone(pert, 2)), (/i,j/), dmat_tuple(4))

               do k = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(5))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                     p_tuple_getone(pert, 3)), (/i,k/), dmat_tuple(6))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(7))

                  do l = 1, 3

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(8))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                        p_tuple_getone(pert, 4)), (/i,l/), dmat_tuple(9))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(10))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(11))

                     if (p_tuple_compare(p_tuple_getone(pert, 4) , &
                        p_tuple_getone(pert, 5))) then
                        maxcomp3 = l
                     else
                        maxcomp3 = 3
                     end if

                     do m = 1, maxcomp3

                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/m/), dmat_tuple(12))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 5)), (/i,m/), dmat_tuple(13))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 5)), (/j,m/), dmat_tuple(14))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 5)), (/k,m/), dmat_tuple(15))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5)), (/l,m/), dmat_tuple(16))

                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/i, j, k, l, m/))

                        print *, 'ERROR: implement g=3 f=2 in openrsp/xcave'
                        stop 1

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do
      end if

      if (nr_geo == 4 .and. nr_el == 1) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                  p_tuple_getone(pert, 2)), (/i,j/), dmat_tuple(4))

               do k = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(5))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                     p_tuple_getone(pert, 3)), (/i,k/), dmat_tuple(6))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                     p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(7))

                  do l = 1, k

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(8))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                        p_tuple_getone(pert, 4)), (/i,l/), dmat_tuple(9))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(10))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                        p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(11))

                     do m = 1, 3

                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 5), (/m/), dmat_tuple(12))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 5)), (/i,m/), dmat_tuple(13))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 5)), (/j,m/), dmat_tuple(14))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 5)), (/k,m/), dmat_tuple(15))

                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 4),  &
                                           p_tuple_getone(pert, 5)), (/l,m/), dmat_tuple(16))

                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/i, j, k, l, m/))

                        print *, 'ERROR: implement g=4 f=1 in openrsp/xcave'
                        stop 1

                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do
         end do
      end if

      deallocate(dmat_tuple)

      if (.not. combination_found) then
         print *, 'ERROR: unsupported geo-el combination in rsp_xcave_interface', &
                  nr_geo, &
                  nr_el
         stop 1
      end if

      do i = 1, property_size
         prop(i) = prop(i) + res(i)
      end do

   end subroutine

















   !> Exchange-correlation perturbed by fields 'field', contracted over
   !> densities 'D', added to Fock matrices 'F'
   subroutine rsp_xcint_interface(pert_labels, D, F, Fg, Fgg)

!     ---------------------------------------------------------------------------
      character(4), intent(in)              :: pert_labels(:)
      type(matrix), intent(in)              :: D(:)
      type(matrix), intent(inout), optional :: F
      type(matrix), intent(inout), optional :: Fg(:)
      type(matrix), intent(inout), optional :: Fgg(:, :)
!     ---------------------------------------------------------------------------
      integer              :: icenter
      integer              :: ioff
      integer              :: imat, i, j
      integer              :: mat_dim
      integer              :: nr_atoms
      integer(c_int)       :: nr_dmat
      real(c_double), allocatable :: xc_dmat(:)
      real(c_double), allocatable :: xc_fmat(:)
      real(c_double)              :: xc_energy
      integer(c_int)       :: get_ave
      integer(c_int)       :: force_sequential
!     ---------------------------------------------------------------------------

      get_ave = 0
      force_sequential = 0

      if (.not. get_is_ks_calculation()) then
         return
      end if

      nr_atoms = get_nr_atoms()
      nr_dmat  = size(D)

      mat_dim = D(1)%nrow
      allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
      xc_dmat = 0.0d0
      do imat = 1, nr_dmat
         call daxpy(mat_dim*mat_dim, 1.0d0, D(imat)%elms, 1, xc_dmat((imat-1)*mat_dim*mat_dim + 1), 1)
      end do

      allocate(xc_fmat(mat_dim*mat_dim))

      if (present(F)) then
         call xcint_wakeup_workers()
         call xcint_integrate(nr_dmat,   &
                              xc_dmat,   &
                              xc_fmat,   &
                              xc_energy, &
                              get_ave,   &
                              0,         &
                              (/0/),     &
                              nr_dmat-1, &
                              (/0, 0/),  &
                              force_sequential)
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms, 1)
      end if

      if (present(Fg)) then
         do i = 1, nr_atoms*3
            call xcint_wakeup_workers()
            call xcint_integrate(nr_dmat,   &
                                 xc_dmat,   &
                                 xc_fmat,   &
                                 xc_energy, &
                                 get_ave,   &
                                 1,         &
                                 (/i/),     &
                                 nr_dmat-1, &
                                 (/0, 0/),  &
                                 force_sequential)
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         do i = 1, nr_atoms*3
            do j = 1, i
               call xcint_wakeup_workers()
               call xcint_integrate(nr_dmat,   &
                                    xc_dmat,   &
                                    xc_fmat,   &
                                    xc_energy, &
                                    get_ave,   &
                                    2,         &
                                    (/i, j/),  &
                                    nr_dmat-1, &
                                    (/0, 0/),  &
                                    force_sequential)
               call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fgg(i, j)%elms, 1)
            end do
         end do
      end if

      deallocate(xc_dmat)
      deallocate(xc_fmat)

   end subroutine

end module
