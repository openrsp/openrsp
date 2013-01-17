module interface_xc

!  interface_xc_init is the ONLY routine
!  that is program specific
!  you are NOT allowed to introduce anything specific
!  to a host program outside of interface_xc_init

   use matrix_defop
#ifndef PRG_DIRAC
   use xcint_main
#endif
   use interface_molecule
   use rsp_field_tuple
   use rsp_indices_and_addressing
   use rsp_sdf_caching

   implicit none

   public interface_xc_init
   public interface_xc_finalize

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT

   ! New xcave routine
   public rsp_xcave_new

   public rsp_xcave
   public rsp_xcint

   public get_is_ks_calculation

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_ks_calculation

contains

   subroutine interface_xc_init()
#ifdef VAR_LSDALTON
     STOP 'interface_xc_init'
#else
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"

      is_ks_calculation = dodft

      is_initialized = .true.
#endif
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
   subroutine rsp_xcave_new(mat_dim, pert, num_blks, blk_sizes, blk_info, property_size, prop, D_sdf)

!     ---------------------------------------------------------------------------
      integer                            :: num_blks, property_size
      type(p_tuple)           :: pert
      type(SDF) :: D_sdf
      complex(8),   dimension(property_size) :: prop, res
      type(matrix), allocatable, dimension(:) :: dmat_tuple
!     ---------------------------------------------------------------------------
      integer                            :: i, j, k, l, maxcomp1, maxcomp2, maxcomp3
      integer                            :: element
      integer                            :: mat_dim
      integer, dimension(num_blks) :: blk_sizes
      integer, dimension(num_blks,3) :: blk_info
      integer                            :: nr_atoms
      real(8)                            :: xc_energy
      integer                            :: nr_pert_geo
      integer                            :: nr_pert_el
!     ---------------------------------------------------------------------------

      res = 0.0

      nr_atoms = get_nr_atoms()
!       mat_dim  = D%nrow

      nr_pert_geo = 0
      nr_pert_el  = 0
      do i = 1, pert%n_perturbations
         select case (pert%plab(i))
            case ('EL  ')
               nr_pert_el = nr_pert_el + 1
            case ('GEO ')
               nr_pert_geo = nr_pert_geo + 1
            case default
               print *, 'Error: rsp_xcave: unknown perturbation type', pert%plab(i)
               stop 1
         end select
      end do

      if (nr_pert_el > 0) then
         res(1:((nr_atoms*3)**nr_pert_geo)*(3*nr_pert_el)) = 0.0d0
      else
         res(1:((nr_atoms*3)**nr_pert_geo)) = 0.0d0
      end if
      if (.not. get_is_ks_calculation()) return


      if (pert%n_perturbations == 1) then

         if (count(pert%plab == 'GEO ') == 1) then

            allocate(dmat_tuple(1))

            call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

            do i = 1, nr_atoms*3
               element = i
               call xc_integrate(             &
                       mat_dim=mat_dim,       &
                       nr_dmat=1,             &
                       dmat=(/dmat_tuple(1)%elms/),       &
                       energy=xc_energy,      &
                       get_ave=.true.,        &
                       fmat=(/0.0d0/),        &
                       geo_coor=(/i/)         &
                    )
               res(element) = cmplx(xc_energy, 0.0d0)
            end do

            deallocate(dmat_tuple)

         else if (count(pert%plab == 'EL  ') == 1) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         else

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         end if
 

      else if (pert%n_perturbations == 2) then


         if (count(pert%plab == 'GEO ') == 2) then

            allocate(dmat_tuple(3))

            call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

            do i = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

               do j = 1, i

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j/))

                  call xc_integrate(                  &
                          mat_dim=mat_dim,            &
                          nr_dmat=3,                  &
                          dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms,  &
                                 dmat_tuple(3)%elms/),&
                          energy=xc_energy,           &
                          get_ave=.true.,             &
                          fmat=(/0.0d0/),             &
                          geo_coor=(/i, j/)           &
                       )

                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do

            deallocate(dmat_tuple)

            
         else if ((count(pert%plab == 'GEO ') == 1) .AND. &
                (count(pert%plab == 'EL  ') == 1)) then

            allocate(dmat_tuple(2))

            call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

            do j = 1, 3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

               do i = 1, nr_atoms*3

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j/))

                  call xc_integrate(        &
                          mat_dim=mat_dim,  &
                          nr_dmat=2,  &
                          dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms/),&
                          energy=xc_energy, &
                          get_ave=.true.,   &
                          fmat=(/0.0d0/),   &
                          geo_coor=(/i/)    &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do

            deallocate(dmat_tuple)

         else if (count(pert%plab == 'EL  ') == 2) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         else

            write(*,*) 'WARNING: rsp_xcave: Unknown perturbation tuple:', pert%plab

         end if


      else if (pert%n_perturbations == 3) then

         if (count(pert%plab == 'GEO ') == 3) then

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

                     call xc_integrate(           &
                             mat_dim=mat_dim,     &
                             nr_dmat=4,           &
                           dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms,  &
                                 dmat_tuple(3)%elms,  &
                                 dmat_tuple(4)%elms/),&
                             energy=xc_energy,    &
                             get_ave=.true.,      &
                             fmat=(/0.0d0/),      &
                             geo_coor=(/i, j, k/) &
                          )
                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do

            deallocate(dmat_tuple)
            
         else if ((count(pert%plab == 'GEO ') == 2) .AND. &
                (count(pert%plab == 'EL  ') == 1)) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

            
         else if ((count(pert%plab == 'GEO ') == 1) .AND. &
                (count(pert%plab == 'EL  ') == 2)) then

            allocate(dmat_tuple(4))

            call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))




            do k = 1, 3

               if (p_tuple_compare(p_tuple_getone(pert, 2) , &
                   p_tuple_getone(pert, 3)) .eqv. .TRUE.) then
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
! write(*,*) 'At i j k', i, j, k, ', got offset:', element

! write(*,*) 'dmat 1', dmat_tuple(1)%elms
! write(*,*) 'dmat 2', dmat_tuple(2)%elms
! write(*,*) 'dmat 3', dmat_tuple(3)%elms
! write(*,*) 'dmat 4', dmat_tuple(4)%elms

                     call xc_integrate(        &
                             mat_dim=mat_dim,  &
                             nr_dmat=4,  &
                             dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms,  &
                                 dmat_tuple(3)%elms,  &
                                 dmat_tuple(4)%elms/),&
                             energy=xc_energy, &
                             get_ave=.true.,   &
                             fmat=(/0.0d0/),   &
                             geo_coor=(/i/)    &
                          )
                     res(element) = cmplx(xc_energy, 0.0d0)
! write(*,*) 'result:', res(element)
                  end do
               end do
            end do

            deallocate(dmat_tuple)

         else if (count(pert%plab == 'EL  ') == 3) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         else

            write(*,*) 'WARNING: rsp_xcave: Unknown perturbation tuple:', pert%plab

         end if



      else if (pert%n_perturbations == 4) then


         if (count(pert%plab == 'GEO ') == 4) then


            allocate(dmat_tuple(11))

            call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

            do i = 1, nr_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

               do j = 1, i

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

                  call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                     p_tuple_getone(pert, 2)), (/i,j/), dmat_tuple(6))

                  do k = 1, j

                     call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(3))

                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                        p_tuple_getone(pert, 3)), (/i,k/), dmat_tuple(7))
                     call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                        p_tuple_getone(pert, 3)), (/j,k/), dmat_tuple(9))


                     do l = 1, k

                        call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(4))


                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 1),  &
                                           p_tuple_getone(pert, 4)), (/i,l/), dmat_tuple(8))
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 2),  &
                                           p_tuple_getone(pert, 4)), (/j,l/), dmat_tuple(10))
                        call sdf_getdata_s(D_sdf, merge_p_tuple(p_tuple_getone(pert, 3),  &
                                           p_tuple_getone(pert, 4)), (/k,l/), dmat_tuple(11))


                        element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                  blk_info, blk_sizes, (/i, j, k, l/))

                        call xc_integrate(               &
                                mat_dim=mat_dim,         &
                                nr_dmat=11,              &
                             dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms,  &
                                 dmat_tuple(3)%elms,  &
                                 dmat_tuple(4)%elms,  &
                                 dmat_tuple(5)%elms,  &
                                 dmat_tuple(6)%elms,  &
                                 dmat_tuple(7)%elms,  &
                                 dmat_tuple(8)%elms,  &
                                 dmat_tuple(9)%elms,  &
                                 dmat_tuple(10)%elms,  &
                                 dmat_tuple(11)%elms/),&
                                energy=xc_energy,        &
                                get_ave=.true.,          &
                                fmat=(/0.0d0/),          &
                                geo_coor=(/i, j, k, l/)  &
                             )
                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do
            
            deallocate(dmat_tuple)

         else if ((count(pert%plab == 'GEO ') == 3) .AND. &
                (count(pert%plab == 'EL  ') == 1)) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         else if ((count(pert%plab == 'GEO ') == 2) .AND. &
                (count(pert%plab == 'EL  ') == 2)) then

            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab
            
         else if ((count(pert%plab == 'GEO ') == 1) .AND. &
                (count(pert%plab == 'EL  ') == 3)) then

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
                   p_tuple_getone(pert, 3)) .eqv. .TRUE.) then
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
                      p_tuple_getone(pert, 4)) .eqv. .TRUE.) then
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

! write(*,*) 'At i j k l', i, j, k, l, ', got offset:', element

                        call xc_integrate(        &
                                mat_dim=mat_dim,  &
                                nr_dmat=8,  &
                             dmat=(/dmat_tuple(1)%elms,  &
                                 dmat_tuple(2)%elms,  &
                                 dmat_tuple(3)%elms,  &
                                 dmat_tuple(4)%elms,  &
                                 dmat_tuple(5)%elms,  &
                                 dmat_tuple(6)%elms,  &
                                 dmat_tuple(7)%elms,  &
                                 dmat_tuple(8)%elms/),&
                                energy=xc_energy, &
                                get_ave=.true.,   &
                                fmat=(/0.0d0/),   &
                                geo_coor=(/l/)    &
                             )

! (The result is added to the property tensor at the very end of the routine)
                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do

            deallocate(dmat_tuple)

         else if (count(pert%plab == 'EL  ') == 4) then


            write(*,*) 'WARNING: rsp_xcave: Unsupported perturbation tuple:', pert%plab

         else

            write(*,*) 'WARNING: rsp_xcave: Unknown perturbation tuple:', pert%plab

         end if

      else

         write(*,*) 'WARNING: rsp_xcave: Unsupported number of perturbations:', &
                    pert%n_perturbations

      end if

      prop = prop + res


   end subroutine























   !> Exchange-correlation perturbed by fields f, averaged over densities D
   subroutine rsp_xcave(pert, res, D, Dg, Dgg, Df, Dff, Dfff)

!     ---------------------------------------------------------------------------
      character(*), intent(in)           :: pert
      complex(8),   intent(out)          :: res(*)
      type(matrix), intent(in)           :: D
      type(matrix), intent(in), optional :: Dg(:)
      type(matrix), intent(in), optional :: Dgg(:, :)
      type(matrix), intent(in), optional :: Df(:)
      type(matrix), intent(in), optional :: Dff(:, :)
      type(matrix), intent(in), optional :: Dfff(:, :, :)
!     ---------------------------------------------------------------------------
      integer                            :: i, j, k, l
      integer                            :: element
      integer                            :: mat_dim
      integer                            :: nr_atoms
      integer                            :: nr_dmat
      real(8)                            :: xc_energy
      real(8), allocatable               :: xc_dmat(:)
      integer                            :: nr_pert_geo
      integer                            :: nr_pert_el
!     ---------------------------------------------------------------------------

      nr_atoms = get_nr_atoms()
      mat_dim  = D%nrow

      nr_pert_geo = 0
      nr_pert_el  = 0
      do i = 1, len(pert)
         select case (pert(i:i))
            case ('f')
               nr_pert_el = nr_pert_el + 1
            case ('g')
               nr_pert_geo = nr_pert_geo + 1
            case default
               print *, 'error: unknown pert component'
               stop 1
         end select
      end do

      if (nr_pert_el > 0) then
         res(1:((nr_atoms*3)**nr_pert_geo)*(3*nr_pert_el)) = 0.0d0
      else
         res(1:((nr_atoms*3)**nr_pert_geo)) = 0.0d0
      end if
      if (.not. get_is_ks_calculation()) return

      select case (pert)
         case ('g')
            do i = 1, nr_atoms*3
               element = i
               call xc_integrate(             &
                       mat_dim=mat_dim,       &
                       nr_dmat=1,             &
                       dmat=(/D%elms/),       &
                       energy=xc_energy,      &
                       get_ave=.true.,        &
                       fmat=(/0.0d0/),        &
                       geo_coor=(/i/)         &
                    )
               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         case ('gg')
            do i = 1, nr_atoms*3
               do j = 1, i
                  element = i &
                          + (j-1)*nr_atoms*3
                  call xc_integrate(                 &
                          mat_dim=mat_dim,           &
                          nr_dmat=3,                 &
                          dmat=(/D%elms,             &
                                 Dg(i)%elms,         &
                                 Dg(j)%elms/),       &
                          energy=xc_energy,          &
                          get_ave=.true.,            &
                          fmat=(/0.0d0/),            &
                          geo_coor=(/i, j/)          &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         case ('gf')
            nr_dmat = 3
            allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
            xc_dmat = 0.0d0
            call dcopy(mat_dim*mat_dim, D%elms, 1, xc_dmat(1), 1)
            do j = 1, 3
               call dcopy(mat_dim*mat_dim, Df(j)%elms, 1, xc_dmat(mat_dim*mat_dim*1 + 1), 1)
               do i = 1, nr_atoms*3
                  element = i &
                          + (j-1)*nr_atoms*3
                  call xc_integrate(        &
                          mat_dim=mat_dim,  &
                          nr_dmat=nr_dmat,  &
                          dmat=xc_dmat,     &
                          energy=xc_energy, &
                          get_ave=.true.,   &
                          fmat=(/0.0d0/),   &
                          geo_coor=(/i/)    &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         case ('gff')
            nr_dmat = 4
            allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
            xc_dmat = 0.0d0
            call dcopy(mat_dim*mat_dim, D%elms, 1, xc_dmat(1), 1)
            do k = 1, 3
               call dcopy(mat_dim*mat_dim, Df(k)%elms, 1, xc_dmat(mat_dim*mat_dim*1 + 1), 1)
               do j = 1, k
                  call dcopy(mat_dim*mat_dim, Df(j)%elms, 1, xc_dmat(mat_dim*mat_dim*2 + 1), 1)
                  call dcopy(mat_dim*mat_dim, Dff(j, k)%elms, 1, xc_dmat(mat_dim*mat_dim*3 + 1), 1)



                  do i = 1, nr_atoms*3
                     element = i                &
                             + (j-1)*nr_atoms*3 &
                             + (k-1)*(nr_atoms*3)*3
write(*,*) 'At i j k', i, j, k, ', got offset:', element

write(*,*) 'dmat 1', D%elms
write(*,*) 'dmat 2', Df(k)%elms
write(*,*) 'dmat 3', Df(j)%elms
write(*,*) 'dmat 4', Dff(j,k)%elms


                     call xc_integrate(        &
                             mat_dim=mat_dim,  &
                             nr_dmat=nr_dmat,  &
                             dmat=xc_dmat,     &
                             energy=xc_energy, &
                             get_ave=.true.,   &
                             fmat=(/0.0d0/),   &
                             geo_coor=(/i/)    &
                          )
                     res(element) = cmplx(xc_energy, 0.0d0)
write(*,*) 'result:', res(element)
                  end do
               end do
            end do
write(*,*) 'size', element
write(*,*) 'end result:', real(res(1:element))

         case ('gfff')
            nr_dmat = 8
            allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
            xc_dmat = 0.0d0
            call dcopy(mat_dim*mat_dim, D%elms, 1, xc_dmat(1), 1)
            do i = 1, 3
               call dcopy(mat_dim*mat_dim, Df(i)%elms, 1, xc_dmat(mat_dim*mat_dim*1 + 1), 1)
               do j = 1, i
                  call dcopy(mat_dim*mat_dim, Df(j)%elms, 1, xc_dmat(mat_dim*mat_dim*2 + 1), 1)
                  call dcopy(mat_dim*mat_dim, Dff(i, j)%elms, 1, xc_dmat(mat_dim*mat_dim*3 + 1), 1)
                  do k = 1, j
                     call dcopy(mat_dim*mat_dim, Df(k)%elms, 1, xc_dmat(mat_dim*mat_dim*4 + 1), 1)
                     call dcopy(mat_dim*mat_dim, Dff(i, k)%elms, 1, xc_dmat(mat_dim*mat_dim*5 + 1), 1)
                     call dcopy(mat_dim*mat_dim, Dff(j, k)%elms, 1, xc_dmat(mat_dim*mat_dim*6 + 1), 1)
                     call dcopy(mat_dim*mat_dim, Dfff(i, j, k)%elms, 1, xc_dmat(mat_dim*mat_dim*7 + 1), 1)
                     do l = 1, nr_atoms*3
                        element = i                    &
                                + (j-1)*nr_atoms*3     &
                                + (k-1)*(nr_atoms*3)*3 &
                                + (l-1)*(nr_atoms*3)*3*3
                        call xc_integrate(        &
                                mat_dim=mat_dim,  &
                                nr_dmat=nr_dmat,  &
                                dmat=xc_dmat,     &
                                energy=xc_energy, &
                                get_ave=.true.,   &
                                fmat=(/0.0d0/),   &
                                geo_coor=(/l/)    &
                             )
                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do
         case ('ggg')
            do i = 1, nr_atoms*3
               do j = 1, i
                  do k = 1, j
                     element = i                &
                             + (j-1)*nr_atoms*3 &
                             + (k-1)*(nr_atoms*3)**2
                     call xc_integrate(           &
                             mat_dim=mat_dim,     &
                             nr_dmat=4,           &
                             dmat=(/D%elms,       &
                                    Dg(i)%elms,   &
                                    Dg(j)%elms,   &
                                    Dg(k)%elms/), &
                             energy=xc_energy,    &
                             get_ave=.true.,      &
                             fmat=(/0.0d0/),      &
                             geo_coor=(/i, j, k/) &
                          )
                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         case ('gggg')
            do i = 1, nr_atoms*3
               do j = 1, i
                  do k = 1, j
                     do l = 1, k
                        element = i                     &
                                + (j-1)*nr_atoms*3      &
                                + (k-1)*(nr_atoms*3)**2 &
                                + (l-1)*(nr_atoms*3)**3
                        call xc_integrate(               &
                                mat_dim=mat_dim,         &
                                nr_dmat=11,              &
                                dmat=(/D%elms,           &
                                       Dg(i)%elms,       &
                                       Dg(j)%elms,       &
                                       Dg(k)%elms,       &
                                       Dg(l)%elms,       &
                                       Dgg(i, j)%elms,   &
                                       Dgg(i, k)%elms,   &
                                       Dgg(i, l)%elms,   &
                                       Dgg(j, k)%elms,   &
                                       Dgg(j, l)%elms,   &
                                       Dgg(k, l)%elms/), &
                                energy=xc_energy,        &
                                get_ave=.true.,          &
                                fmat=(/0.0d0/),          &
                                geo_coor=(/i, j, k, l/)  &
                             )
                        res(element) = cmplx(xc_energy, 0.0d0)
                     end do
                  end do
               end do
            end do
         case default
            print *, 'error: perturbation not implemented in xcave'
            stop 1
      end select

      if (allocated(xc_dmat)) deallocate(xc_dmat)

   end subroutine


   !> Exchange-correlation perturbed by fields 'field', contracted over
   !> densities 'D', added to Fock matrices 'F'
   subroutine rsp_xcint(D, F, Fg, Fgg)

!     ---------------------------------------------------------------------------
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
                           geo_coor=(/0/)    &
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
                              geo_coor=(/i/)    &
                             )
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         do i = 1, nr_atoms*3
            do j = 1, i
               call xc_integrate(                  &
                                 mat_dim=mat_dim,  &
                                 nr_dmat=nr_dmat,  &
                                 dmat=xc_dmat,     &
                                 energy=xc_energy, &
                                 get_ave=.false.,  &
                                 fmat=xc_fmat,     &
                                 geo_coor=(/i, j/) &
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

   subroutine external_rsp_xcint(nr_ao, dmat, fmat, xc_energy)

      use interface_ao_specific
      use xcint_main

      integer      :: nr_ao
      real(8)      :: dmat(nr_ao, nr_ao)
      real(8)      :: fmat(*)
      real(8)      :: xc_energy

      integer      :: i, j, k

      real(8), allocatable :: xcmat(:, :)
      allocate(xcmat(nr_ao, nr_ao))
      xcmat = 0.0d0

      xc_energy = 0.0d0

      call interface_ao_write()
      call xc_integrate(                  &
                        mat_dim=nr_ao,    &
                        nr_dmat=1,        &
                        dmat=0.5d0*dmat,  &
                        energy=xc_energy, &
                        get_ave=.false.,  &
                        fmat=xcmat,       &
                        geo_coor=(/0/)    &
                       )

      k = 0
      do i = 1, nr_ao
         do j = 1, i
            k = k + 1
            fmat(k) = fmat(k) + xcmat(j, i)
         end do
      end do

!     print *, 'raboof xc energy', xc_energy
      deallocate(xcmat)

   end subroutine
