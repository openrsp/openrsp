module interface_xc

!  interface_xc_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_xc_init

   use matrix_defop
   use xcint_integrator
   use interface_molecule
   use rsp_field_tuple
   use rsp_indices_and_addressing
   use rsp_sdf_caching

   implicit none

   public interface_xc_init
   public interface_xc_finalize

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT

   public rsp_xcave_interface
   public rsp_xcint_interface

   public get_is_ks_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_ks_calculation

contains

   subroutine set_is_ks_calculation()

#include "maxorb.h"
#include "infinp.h"

      is_ks_calculation = dodft

   end subroutine

   subroutine interface_xc_init()

      use xcint_interface

#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "shells.h"
#include "primit.h"

      integer, allocatable :: nr_primitives_per_shell(:)
      real(8), allocatable :: primitive_exp(:, :)
      real(8), allocatable :: contraction_coef(:, :)
      integer, allocatable :: l_quantum_nr(:)
      integer, allocatable :: shell_center(:)
      integer              :: i, j, icount, iprim, ishell, iround
      logical              :: is_spherical

      allocate(nr_primitives_per_shell(kmax))
      do iround = 1, 2
         if (iround == 2) then
            allocate(primitive_exp(maxval(nr_primitives_per_shell), kmax))
            allocate(contraction_coef(maxval(nr_primitives_per_shell), kmax))
         end if
         do ishell = 1, kmax
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

      allocate(l_quantum_nr(kmax))
      allocate(shell_center(kmax))
      do ishell = 1, kmax
         l_quantum_nr(ishell) = nhkt(ishell) - 1
         shell_center(ishell) = ncent(ishell)
      end do

      is_spherical = .false.
      do ishell = 1, kmax
         if (sphr(ishell)) then
            is_spherical = .true.
         end if
      end do

      call xcint_interface_init(nr_centers=nucind,                                           &
                                center_charge=charge,                                        &
                                center_xyz=cord,                                             &
                                nr_shells=kmax,                                              &
                                l_quantum_nr=l_quantum_nr,                                   &
                                shell_center=shell_center,                                   &
                                nr_primitives_per_shell=nr_primitives_per_shell,             &
                                max_nr_primitives_per_shell=maxval(nr_primitives_per_shell), &
                                primitive_exp=primitive_exp,                                 &
                                contraction_coef=contraction_coef,                           &
                                is_spherical=is_spherical)

      deallocate(nr_primitives_per_shell)
      deallocate(primitive_exp)
      deallocate(contraction_coef)
      deallocate(l_quantum_nr)
      deallocate(shell_center)

      call set_is_ks_calculation()

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
      integer                      :: i, j, k, l, m, p
      integer                      :: maxcomp1, maxcomp2, maxcomp3, maxcomp4
      integer                      :: element
      integer                      :: nr_atoms
      integer                      :: nr_geo
      integer                      :: nr_el
      logical                      :: combination_found
      real(8)                      :: xc_energy
      type(matrix), allocatable, dimension(:) :: dmat_tuple
!     --------------------------------------------------------------------------

      res = 0.0

#ifndef PRG_DIRAC
      if (.not. get_is_ks_calculation()) return

      nr_atoms = get_nr_atoms()

      nr_geo = count(pert%plab == 'GEO ')
      nr_el  = count(pert%plab == 'EL  ')

      combination_found = .false.

      if (nr_geo == 0 .and. nr_el == 1) then
         combination_found = .true.
         ! nothing to do
      end if

      if (nr_geo == 0 .and. nr_el == 2) then
         combination_found = .true.
         ! nothing to do
      end if

      if (nr_geo == 0 .and. nr_el == 3) then
         combination_found = .true.
         ! nothing to do
      end if

      if (nr_geo == 0 .and. nr_el == 4) then
         combination_found = .true.
         ! nothing to do
      end if

      if (nr_geo == 0 .and. nr_el == 5) then
         combination_found = .true.
         ! nothing to do
      end if

      if (nr_geo == 1 .and. nr_el == 0) then
         combination_found = .true.

         allocate(dmat_tuple(1))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3
            element = i
            call xc_integrate(                   &
                    mat_dim=mat_dim,             &
                    nr_dmat=1,                   &
                    dmat=(/dmat_tuple(1)%elms/), &
                    energy=xc_energy,            &
                    get_ave=.true.,              &
                    fmat=(/0.0d0/),              &
                    geo_coor=(/i/),              &
                    pert_labels=(/pert%plab/),   &
                    kn=kn                        &
                 )
            res(element) = cmplx(xc_energy, 0.0d0)
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 2 .and. nr_el == 0) then
         combination_found = .true.

         allocate(dmat_tuple(2))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3

          ! call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xc_integrate(                   &
                       mat_dim=mat_dim,             &
                       nr_dmat=2,                   &
                       dmat=(/dmat_tuple(1)%elms,   &
                              dmat_tuple(2)%elms/), &
                       energy=xc_energy,            &
                       get_ave=.true.,              &
                       fmat=(/0.0d0/),              &
                       geo_coor=(/i, j/),           &
                    pert_labels=(/pert%plab/),   &
                       kn=kn                        &
                    )

               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 1 .and. nr_el == 1) then
         combination_found = .true.

         allocate(dmat_tuple(2))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do j = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

            do i = 1, nr_atoms*3

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xc_integrate(                   &
                       mat_dim=mat_dim,             &
                       nr_dmat=2,                   &
                       dmat=(/dmat_tuple(1)%elms,   &
                              dmat_tuple(2)%elms/), &
                       energy=xc_energy,            &
                       get_ave=.true.,              &
                       fmat=(/0.0d0/),              &
                       geo_coor=(/i/),              &
                    pert_labels=(/pert%plab/),   &
                       kn=kn                        &
                    )
               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 3 .and. nr_el == 0) then
         combination_found = .true.

         allocate(dmat_tuple(4))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, nr_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do k = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xc_integrate(                &
                          mat_dim=mat_dim,          &
                          nr_dmat=4,                &
                        dmat=(/dmat_tuple(1)%elms,  &
                              dmat_tuple(2)%elms,   &
                              dmat_tuple(3)%elms,   &
                              dmat_tuple(4)%elms/), &
                          energy=xc_energy,         &
                          get_ave=.true.,           &
                          fmat=(/0.0d0/),           &
                          geo_coor=(/i, j, k/),     &
                    pert_labels=(/pert%plab/),   &
                          kn=kn                     &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 2 .and. nr_el == 1) then
         combination_found = .true.

         ! Set up for (k,n) = (1,1)
         ! Almost sure of ordering in dmat_tuple

         allocate(dmat_tuple(4))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do k = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

            do j = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do i = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))


                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))


                  call xc_integrate(                 &
                          mat_dim=mat_dim,           &
                          nr_dmat=4,                 &
                          dmat=(/dmat_tuple(1)%elms, &
                              dmat_tuple(2)%elms,    &
                              dmat_tuple(3)%elms,    &
                              dmat_tuple(4)%elms/),  &
                          energy=xc_energy,          &
                          get_ave=.true.,            &
                          fmat=(/0.0d0/),            &
                          geo_coor=(/i, j/),         &
                    pert_labels=(/pert%plab/),   &
                          kn=kn                      &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)

               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 1 .and. nr_el == 2) then
         combination_found = .true.

         allocate(dmat_tuple(4))

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

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

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xc_integrate(                 &
                          mat_dim=mat_dim,           &
                          nr_dmat=4,                 &
                          dmat=(/dmat_tuple(1)%elms, &
                              dmat_tuple(2)%elms,    &
                              dmat_tuple(3)%elms,    &
                              dmat_tuple(4)%elms/),  &
                          energy=xc_energy,          &
                          get_ave=.true.,            &
                          fmat=(/0.0d0/),            &
                          geo_coor=(/i/),            &
                    pert_labels=(/pert%plab/),   &
                          kn=kn                      &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 4 .and. nr_el == 0) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         allocate(dmat_tuple(8))

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

                     call xc_integrate(                   &
                             mat_dim=mat_dim,             &
                             nr_dmat=8,                   &
                             dmat=(/dmat_tuple(1)%elms,   &
                                    dmat_tuple(2)%elms,   &
                                    dmat_tuple(3)%elms,   &
                                    dmat_tuple(4)%elms,   &
                                    dmat_tuple(5)%elms,   &
                                    dmat_tuple(6)%elms,   &
                                    dmat_tuple(7)%elms,   &
                                    dmat_tuple(8)%elms/), &
                             energy=xc_energy,            &
                             get_ave=.true.,              &
                             fmat=(/0.0d0/),              &
                             geo_coor=(/i, j, k, l/),     &
                             pert_labels=(/pert%plab/),   &
                             kn=kn                        &
                          )

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 3 .and. nr_el == 1) then
         combination_found = .true.

            ! Using (k,n) = (1,2)

            allocate(dmat_tuple(8))

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

                        call xc_integrate(              &
                                mat_dim=mat_dim,        &
                                nr_dmat=8,              &
                             dmat=(/dmat_tuple(1)%elms, &
                                 dmat_tuple(2)%elms,    &
                                 dmat_tuple(3)%elms,    &
                                 dmat_tuple(4)%elms,    &
                                 dmat_tuple(5)%elms,    &
                                 dmat_tuple(6)%elms,    &
                                 dmat_tuple(7)%elms,    &
                                 dmat_tuple(8)%elms/),  &
                                energy=xc_energy,       &
                                get_ave=.true.,         &
                                fmat=(/0.0d0/),         &
                                geo_coor=(/i, j, k/),         &
                       pert_labels=(/pert%plab/),   &
                                kn=kn                   &
                             )

                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do

            deallocate(dmat_tuple)
      end if

      if (nr_geo == 2 .and. nr_el == 2) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         allocate(dmat_tuple(8))

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

                     call xc_integrate(              &
                             mat_dim=mat_dim,        &
                             nr_dmat=8,              &
                          dmat=(/dmat_tuple(1)%elms, &
                              dmat_tuple(2)%elms,    &
                              dmat_tuple(3)%elms,    &
                              dmat_tuple(4)%elms,    &
                              dmat_tuple(5)%elms,    &
                              dmat_tuple(6)%elms,    &
                              dmat_tuple(7)%elms,    &
                              dmat_tuple(8)%elms/),  &
                             energy=xc_energy,       &
                             get_ave=.true.,         &
                             fmat=(/0.0d0/),         &
                             geo_coor=(/i, j/),         &
                             pert_labels=(/pert%plab/),   &
                             kn=kn                   &
                          )

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if

      if (nr_geo == 1 .and. nr_el == 3) then
         combination_found = .true.

         ! The array of density matrices should have length 4
         allocate(dmat_tuple(8))

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

                     call xc_integrate(              &
                             mat_dim=mat_dim,        &
                             nr_dmat=8,              &
                          dmat=(/dmat_tuple(1)%elms, &
                              dmat_tuple(2)%elms,    &
                              dmat_tuple(3)%elms,    &
                              dmat_tuple(4)%elms,    &
                              dmat_tuple(5)%elms,    &
                              dmat_tuple(6)%elms,    &
                              dmat_tuple(7)%elms,    &
                              dmat_tuple(8)%elms/),  &
                             energy=xc_energy,       &
                             get_ave=.true.,         &
                             fmat=(/0.0d0/),         &
                             geo_coor=(/l/),         &
                    pert_labels=(/pert%plab/),   &
                             kn=kn                   &
                          )

                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)
      end if



      if (nr_geo == 1 .and. nr_el == 4) then
         combination_found = .true.

         ! Using (k,n) = (0,4)

         allocate(dmat_tuple(16))

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

                        call xc_integrate(              &
                                mat_dim=mat_dim,        &
                                nr_dmat=size(dmat_tuple),    &
                                dmat=(/(dmat_tuple(p)%elms, p = 1, size(dmat_tuple))/),  &
                                energy=xc_energy,       &
                                get_ave=.true.,         &
                                fmat=(/0.0d0/),         &
                                geo_coor=(/m/),         &
                                pert_labels=(/pert%plab/),   &
                                kn=kn                   &
                             )

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do


         deallocate(dmat_tuple)
      end if


     if (nr_geo == 2 .and. nr_el == 3) then
         ! radovan: note to myself to consider changing dmats 12 and 13
         combination_found = .true.

         ! Using (k,n) = (1,3)

         allocate(dmat_tuple(16))

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

                        call xc_integrate(              &
                                mat_dim=mat_dim,        &
                                nr_dmat=size(dmat_tuple),    &
                                dmat=(/(dmat_tuple(p)%elms, p = 1, size(dmat_tuple))/),  &
                                energy=xc_energy,       &
                                get_ave=.true.,         &
                                fmat=(/0.0d0/),         &
                                geo_coor=(/i, j/),         &
                                pert_labels=(/pert%plab/),   &
                                kn=kn                   &
                             )

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)

      end if



     if (nr_geo == 3 .and. nr_el == 2) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         allocate(dmat_tuple(16))

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

                        call xc_integrate(              &
                                mat_dim=mat_dim,        &
                                nr_dmat=size(dmat_tuple),    &
                                dmat=(/(dmat_tuple(p)%elms, p = 1, size(dmat_tuple))/),  &
                                energy=xc_energy,       &
                                get_ave=.true.,         &
                                fmat=(/0.0d0/),         &
                                geo_coor=(/i, j, k/),         &
                                pert_labels=(/pert%plab/),   &
                                kn=kn                   &
                             )

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)

      end if


     if (nr_geo == 4 .and. nr_el == 1) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         allocate(dmat_tuple(16))

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

                        call xc_integrate(              &
                                mat_dim=mat_dim,        &
                                nr_dmat=size(dmat_tuple),    &
                                dmat=(/(dmat_tuple(p)%elms, p = 1, size(dmat_tuple))/),  &
                                energy=xc_energy,       &
                                get_ave=.true.,         &
                                fmat=(/0.0d0/),         &
                                geo_coor=(/i, j, k, l/),         &
                                pert_labels=(/pert%plab/),   &
                                kn=kn                   &
                             )

                        res(element) = cmplx(xc_energy, 0.0d0)

                     end do
                  end do
               end do
            end do
         end do

         deallocate(dmat_tuple)

      end if




      if (.not. combination_found) then
         print *, 'ERROR: unsupported geo-el combination in rsp_xcave_interface', &
                  nr_geo, &
                  nr_el
         stop 1
      end if
#endif /* ifdef PRG_DIRAC */

      prop = prop + res

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
      integer              :: nr_dmat
      real(8), allocatable :: xc_dmat(:)
      real(8), allocatable :: xc_fmat(:)
      real(8)              :: xc_energy
!     ---------------------------------------------------------------------------

#ifndef PRG_DIRAC
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
         call xc_integrate(                  &
                           mat_dim=mat_dim,  &
                           nr_dmat=nr_dmat,  &
                           dmat=xc_dmat,     &
                           energy=xc_energy, &
                           get_ave=.false.,  &
                           fmat=xc_fmat,     &
                           geo_coor=(/0/),   &
                           pert_labels=pert_labels, &
                           kn=(/0, 0/)       &
                          )
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms, 1)
      end if

      if (present(Fg)) then
         do i = 1, nr_atoms*3
            call xc_integrate(                  &
                              mat_dim=mat_dim,  &
                              nr_dmat=nr_dmat,  &
                              dmat=xc_dmat,     &
                              energy=xc_energy, &
                              get_ave=.false.,  &
                              fmat=xc_fmat,     &
                              geo_coor=(/i/),   &
                           pert_labels=pert_labels, &
                              kn=(/0, 0/)       &
                             )
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         do i = 1, nr_atoms*3
            do j = 1, i
               call xc_integrate(                   &
                                 mat_dim=mat_dim,   &
                                 nr_dmat=nr_dmat,   &
                                 dmat=xc_dmat,      &
                                 energy=xc_energy,  &
                                 get_ave=.false.,   &
                                 fmat=xc_fmat,      &
                                 geo_coor=(/i, j/), &
                           pert_labels=pert_labels, &
                                 kn=(/0, 0/)        &
                                )
               call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fgg(i, j)%elms, 1)
            end do
         end do
      end if

      deallocate(xc_dmat)
      deallocate(xc_fmat)
#endif

   end subroutine

end module

   subroutine external_rsp_xcint_interface(nr_ao, dmat, fmat, xc_energy)

      use xcint_integrator

      integer      :: nr_ao
      real(8)      :: dmat(nr_ao, nr_ao)
      real(8)      :: fmat(*)
      real(8)      :: xc_energy

      integer      :: i, j, k

      real(8), allocatable :: xcmat(:, :)
      allocate(xcmat(nr_ao, nr_ao))
      xcmat = 0.0d0

      xc_energy = 0.0d0

#ifndef PRG_DIRAC
      call xc_integrate(                  &
                        mat_dim=nr_ao,    &
                        nr_dmat=1,        &
                        dmat=0.5d0*dmat,  &
                        energy=xc_energy, &
                        get_ave=.false.,  &
                        fmat=xcmat,       &
                        geo_coor=(/0/),   &
                        pert_labels=(/'NONE'/), &
                        kn=(/0, 0/)       &
                       )

      k = 0
      do i = 1, nr_ao
         do j = 1, i
            k = k + 1
            fmat(k) = fmat(k) + xcmat(j, i)
         end do
      end do
#endif /* ifndef PRG_DIRAC */

!     print *, 'raboof xc energy', xc_energy
      deallocate(xcmat)

   end subroutine
