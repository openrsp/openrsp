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

   implicit none

   public interface_xc_init
   public interface_xc_finalize

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT

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
   subroutine rsp_xcave(pert, res, D, Dg, Dgg, Df, Dff)

!     ---------------------------------------------------------------------------
      character(*), intent(in)           :: pert
      complex(8),   intent(out)          :: res(*)
      type(matrix), intent(in)           :: D
      type(matrix), intent(in), optional :: Dg(:)
      type(matrix), intent(in), optional :: Dgg(:, :)
      type(matrix), intent(in), optional :: Df(:)
      type(matrix), intent(in), optional :: Dff(:, :)
!     ---------------------------------------------------------------------------
      integer                            :: i, j, k, l
      integer                            :: element
      integer                            :: mat_dim
      integer                            :: nr_atoms
      integer                            :: nr_dmat
      real(8)                            :: xc_energy
      real(8), allocatable               :: xc_dmat(:)
!     ---------------------------------------------------------------------------

      nr_atoms = get_nr_atoms()
      mat_dim  = D%nrow

      select case (pert)
         case ('g')
            res(1:(nr_atoms*3)**len(pert)) = 0.0d0
            if (.not. get_is_ks_calculation()) return
            do i = 1, nr_atoms*3
               element = i
               call xc_integrate(                &
                       xc_mat_dim=mat_dim,       &
                       xc_nr_dmat=1,             &
                       xc_dmat=(/D%elms_alpha/), &
                       xc_nr_geo_pert=1,         &
                       xc_nr_fld_pert=0,         &
                       xc_energy=xc_energy,      &
                       xc_geo_coor=(/i/)         &
                    )
               res(element) = cmplx(xc_energy, 0.0d0)
            end do
         case ('gg')
            res(1:(nr_atoms*3)**len(pert)) = 0.0d0
            if (.not. get_is_ks_calculation()) return
            do i = 1, nr_atoms*3
               do j = 1, i
                  element = i &
                          + (j-1)*nr_atoms*3
                  call xc_integrate(                    &
                          xc_mat_dim=mat_dim,           &
                          xc_nr_dmat=3,                 &
                          xc_dmat=(/D%elms_alpha,       &
                                    Dg(i)%elms_alpha,   &
                                    Dg(j)%elms_alpha/), &
                          xc_nr_geo_pert=2,             &
                          xc_nr_fld_pert=0,             &
                          xc_energy=xc_energy,          &
                          xc_geo_coor=(/i, j/)          &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         case ('gf')
!           res(1:(nr_atoms*3)**len(pert)) = 0.0d0
            if (.not. get_is_ks_calculation()) return
            nr_dmat = 3
            allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
            xc_dmat = 0.0d0
            call dcopy(mat_dim*mat_dim, D%elms_alpha, 1, xc_dmat(1), 1)
            do j = 1, 3
               call dcopy(mat_dim*mat_dim, Df(j)%elms_alpha, 1, xc_dmat(mat_dim*mat_dim*2 + 1), 1)
               do i = 1, nr_atoms*3
                  element = i &
                          + (j-1)*nr_atoms*3
                  call xc_integrate(           &
                          xc_mat_dim=mat_dim,  &
                          xc_nr_dmat=nr_dmat,  &
                          xc_dmat=xc_dmat,     &
                          xc_nr_geo_pert=1,    &
                          xc_nr_fld_pert=1,    &
                          xc_energy=xc_energy, &
                          xc_geo_coor=(/i, 0/) &
                       )
                  res(element) = cmplx(xc_energy, 0.0d0)
               end do
            end do
         case ('gff')
            ! nothing done so far
         case ('ggg')
            res(1:(nr_atoms*3)**len(pert)) = 0.0d0
            if (.not. get_is_ks_calculation()) return
            do i = 1, nr_atoms*3
               do j = 1, i
                  do k = 1, j
                     element = i                &
                             + (j-1)*nr_atoms*3 &
                             + (k-1)*(nr_atoms*3)**2
                     call xc_integrate(                    &
                             xc_mat_dim=mat_dim,           &
                             xc_nr_dmat=4,                 &
                             xc_dmat=(/D%elms_alpha,       &
                                       Dg(i)%elms_alpha,   &
                                       Dg(j)%elms_alpha,   &
                                       Dg(k)%elms_alpha/), &
                             xc_nr_geo_pert=3,             &
                             xc_nr_fld_pert=0,             &
                             xc_energy=xc_energy,          &
                             xc_geo_coor=(/i, j, k/)       &
                          )
                     res(element) = cmplx(xc_energy, 0.0d0)
                  end do
               end do
            end do
         case ('gggg')
            res(1:(nr_atoms*3)**len(pert)) = 0.0d0
            if (.not. get_is_ks_calculation()) return
            do i = 1, nr_atoms*3
               do j = 1, i
                  do k = 1, j
                     do l = 1, k
                        element = i                     &
                                + (j-1)*nr_atoms*3      &
                                + (k-1)*(nr_atoms*3)**2 &
                                + (l-1)*(nr_atoms*3)**3
                        call xc_integrate(                        &
                                xc_mat_dim=mat_dim,               &
                                xc_nr_dmat=11,                    &
                                xc_dmat=(/D%elms_alpha,           &
                                          Dg(i)%elms_alpha,       &
                                          Dg(j)%elms_alpha,       &
                                          Dg(k)%elms_alpha,       &
                                          Dg(l)%elms_alpha,       &
                                          Dgg(i, j)%elms_alpha,   &
                                          Dgg(i, k)%elms_alpha,   &
                                          Dgg(i, l)%elms_alpha,   &
                                          Dgg(j, k)%elms_alpha,   &
                                          Dgg(j, l)%elms_alpha,   &
                                          Dgg(k, l)%elms_alpha/), &
                                xc_nr_geo_pert=4,                 &
                                xc_nr_fld_pert=0,                 &
                                xc_energy=xc_energy,              &
                                xc_geo_coor=(/i, j, k, l/)        &
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
         call daxpy(mat_dim*mat_dim, 1.0d0, D(imat)%elms_alpha, 1, xc_dmat((imat-1)*mat_dim*mat_dim + 1), 1)
      end do

      allocate(xc_fmat(mat_dim*mat_dim))

      if (present(F)) then
         call xc_integrate(                     &
                           xc_mat_dim=mat_dim,  &
                           xc_nr_dmat=nr_dmat,  &
                           xc_dmat=xc_dmat,     &
                           xc_nr_geo_pert=0,    &
                           xc_nr_fld_pert=0,    &
                           xc_energy=xc_energy, &
                           xc_fmat=xc_fmat      &
                          )
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms_alpha, 1)
      end if

      if (present(Fg)) then
         do i = 1, nr_atoms*3
            call xc_integrate(                     &
                              xc_mat_dim=mat_dim,  &
                              xc_nr_dmat=nr_dmat,  &
                              xc_dmat=xc_dmat,     &
                              xc_nr_geo_pert=1,    &
                              xc_nr_fld_pert=0,    &
                              xc_energy=xc_energy, &
                              xc_fmat=xc_fmat,     &
                              xc_geo_coor=(/i/)    &
                             )
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fg(i)%elms_alpha, 1)
         end do
      end if

      if (present(Fgg)) then
         do i = 1, nr_atoms*3
            do j = 1, i
               call xc_integrate(                     &
                                 xc_mat_dim=mat_dim,  &
                                 xc_nr_dmat=nr_dmat,  &
                                 xc_dmat=xc_dmat,     &
                                 xc_nr_geo_pert=2,    &
                                 xc_nr_fld_pert=0,    &
                                 xc_energy=xc_energy, &
                                 xc_fmat=xc_fmat,     &
                                 xc_geo_coor=(/i, j/) &
                                )
               call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fgg(i, j)%elms_alpha, 1)
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
      call xc_integrate(                     &
                        xc_mat_dim=nr_ao,    &
                        xc_nr_dmat=1,        &
                        xc_dmat=0.5d0*dmat,  &
                        xc_nr_geo_pert=0,    &
                        xc_nr_fld_pert=0,    &
                        xc_energy=xc_energy, &
                        xc_fmat=xcmat        &
                       )

      k = 0
      do i = 1, nr_ao
         do j = 1, i
            k = k + 1
            fmat(k) = fmat(k) + xcmat(j, i)
         end do
      end do

      print *, 'raboof xc energy', xc_energy
      deallocate(xcmat)

   end subroutine
