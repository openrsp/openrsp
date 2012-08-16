module interface_interest

   use openrsp_cfg

   implicit none

   public interest_get_int
   public interest_get_ave

   private

   type type_atom
      real(8) :: charge
      real(8) :: coordinate_x
      real(8) :: coordinate_y
      real(8) :: coordinate_z
      real(8) :: gnu_exponent
   end type

   type type_gto
      integer :: idx
      integer :: offset
      integer :: lvalue
      integer :: origin
      integer :: sdegen
      integer :: cdegen
      real(8) :: ex(1)
      real(8) :: coefficient(1)
      integer :: center
   end type

   integer, save :: shell_start(2) = 0
   integer, save :: shell_end(2) = 0

   type(type_atom), allocatable :: atom(:)
   type(type_gto),  allocatable :: gto(:)

   logical :: interface_is_initialized = .false.

   integer, parameter   :: maxl = 10
   integer, allocatable :: imat_u(:, :)
   integer, allocatable :: imat_d(:, :, :)
   real(8), allocatable :: imat_f(:, :, :)

   integer :: nr_centers

contains

   subroutine initialize_interest_eri_diff()

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"

      integer        :: i
      integer        :: j
      integer        :: ij
      integer        :: ier
      type(type_gto) :: tmpgto

      if (interface_is_initialized) return

      call init_arrays()

      !> initialize InteRest integral package
      call interest_initialize()

      !> test if the large component basis is uncontracted
      if (nplrg /= nlarge) then
         print *, ' Error: contracted basis set found => stop!'
         print *, ' Note:  uncontracted basis sets are only supported'
         stop 1
      endif

      !> interface molecular data
      allocate(atom(nucind), stat=ier)
      if (ier /= 0) stop 'Error in allocation: atom(:)'
      do i = 1, size(atom)
         atom(i) = type_atom(charge(i), &
                             cord(1,i), &
                             cord(2,i), &
                             cord(3,i), &
                             gnuexp(i))
      enddo

      !> interface basis set data
      !> fixme: contracted basis sets
      allocate(gto(nlrgsh+nsmlsh), stat=ier)
      if (ier /= 0) stop 'Error in allocation: gto(:)'
      i  = 0
      ij = 0
      do while (i < size(gto))
        do j = 1, nrco(i+1)
           i = i + 1
           gto(i) = type_gto(ij,                      &
                             0,                       &
                             nhkt(i)-1,               &
                             ncent(i),                &
                             (2*nhkt(i)-1),           &
                             (nhkt(i)*(nhkt(i)+1))/2, &
                             priexp(i),               &
                             priccf(i,j),             &
                             ncent(i))
           !> note: cartesian indexation
           ij = ij + (nhkt(i)*(nhkt(i)+1))/2
        end do
      end do

      do i = 1, size(gto)-1
         gto(i+1)%offset = gto(i)%offset + gto(i)%cdegen
      enddo

      ! get total number of centers
      nr_centers = 0
      do i = 1, size(gto)
         if (gto(i)%center > nr_centers) nr_centers = nr_centers + 1
      enddo

      shell_start(1) = 1
      shell_end(1)   = nlrgsh
      shell_start(2) = shell_end(1) + 1
      shell_end(2)   = shell_end(1) + nsmlsh

      print *, 'total number of basis function shells:    ', size(gto)
      print *, 'total number of spherical basis functions:', sum(gto(:)%sdegen)
      print *, 'total number of cartesian basis functions:', sum(gto(:)%cdegen)

      interface_is_initialized = .true.

   end subroutine

   subroutine interest_get_int(ndim, dmat, gmat, order, only_icoor)

      integer, intent(in)  :: ndim
      real(8), intent(out) :: gmat(ndim, ndim, *)
      real(8), intent(in)  :: dmat(ndim, ndim, *)
      integer, intent(in)  :: order
      integer, intent(in)  :: only_icoor

      call interest_eri_diff_block(ndim,           &
                                   dmat,           &
                                   gmat,           &
                                   order,          &
                                   (/1, 1, 1, 1/), &
                                   (/0.0d0/),      &
                                   .false.,        &
                                   only_icoor)
#ifdef PRG_DIRAC
      if (.not. openrsp_cfg_skip_llss) then
         call interest_eri_diff_block(ndim,           &
                                      dmat,           &
                                      gmat,           &
                                      order,          &
                                      (/1, 1, 2, 2/), &
                                      (/0.0d0/),      &
                                      .false.,        &
                                      only_icoor)
         call interest_eri_diff_block(ndim,           &
                                      dmat,           &
                                      gmat,           &
                                      order,          &
                                      (/2, 2, 1, 1/), &
                                      (/0.0d0/),      &
                                      .false.,        &
                                      only_icoor)
      end if
#endif

   end subroutine

   subroutine interest_get_ave(ndim, dmat1, dmat2, order, ave)

      integer, intent(in)  :: ndim
      real(8), intent(in)  :: dmat1(ndim, ndim, *)
      real(8), intent(in)  :: dmat2(ndim, ndim, *)
      integer, intent(in)  :: order
      real(8), intent(out) :: ave(*)

      ave(1) = 0.0d0

      call interest_eri_diff_block(ndim,           &
                                   dmat1,          &
                                   dmat2,          &
                                   order,          &
                                   (/1, 1, 1, 1/), &
                                   ave,            &
                                   .true.,         &
                                   0)
#ifdef PRG_DIRAC
      if (.not. openrsp_cfg_skip_llss) then
         call interest_eri_diff_block(ndim,           &
                                      dmat1,          &
                                      dmat2,          &
                                      order,          &
                                      (/1, 1, 2, 2/), &
                                      ave,            &
                                      .true.,         &
                                      0)
         call interest_eri_diff_block(ndim,           &
                                      dmat1,          &
                                      dmat2,          &
                                      order,          &
                                      (/2, 2, 1, 1/), &
                                      ave,            &
                                      .true.,         &
                                      0)
      end if
#endif

   end subroutine

   subroutine interest_eri_diff_block(ndim,    &
                                      dmat,    &
                                      gmat,    &
                                      order,   &
                                      iblocks, &
                                      ave,     &
                                      get_ave, &
                                      only_icoor)

      integer, intent(in)    :: ndim
      real(8)                :: dmat(ndim, ndim, *)
      real(8)                :: gmat(ndim, ndim, *)
      integer, intent(in)    :: order
      integer, intent(in)    :: iblocks(4)
      real(8)                :: ave(*)
      logical, intent(in)    :: get_ave
      integer, intent(in)    :: only_icoor

      integer, parameter :: max_nr_integrals = 194481 !fixme hardcoded
      integer, parameter :: max_ave_length   = 100    !fixme hardcoded

      real(8) :: gint(max_nr_integrals)
      real(8) :: g_u(max_nr_integrals)
      real(8) :: g_d(max_nr_integrals)
      real(8) :: ex(4), coef(4), xyz(3, 4)
      integer :: ang(4), deg(4), off(4)
      real(8) :: f

      integer :: ii, ij, ik, il
      integer :: cent(4)
      integer :: icent, ixyz
      integer :: nr_integrals
      integer :: nr_elements
      integer :: icent_start, icent_end
      integer :: ixyz_start,  ixyz_end

      ! make sure it is initialized
      ! it does not cost anything
      ! routine will return if already initialized
      call initialize_interest_eri_diff()

      if (get_ave) then
         icent_start = 1
         icent_end   = nr_centers
         ixyz_start  = 1
         ixyz_end    = 3
      else
         icent_start = 1 + (only_icoor-1)/3
         ixyz_start  = 1 + mod((only_icoor-1), 3)
         icent_end   = icent_start
         ixyz_end    = ixyz_start
      end if

      i_function_loop: do ii = shell_start(iblocks(1)), shell_end(iblocks(1))
         ang(1)    =      gto(ii)%lvalue
         off(1)    =      gto(ii)%offset
         deg(1)    =      gto(ii)%cdegen
         ex(1)     =      gto(ii)%ex(1)
         coef(1)   =      gto(ii)%coefficient(1)
         cent(1)   =      gto(ii)%center
         xyz(1, 1) = atom(gto(ii)%origin)%coordinate_x
         xyz(2, 1) = atom(gto(ii)%origin)%coordinate_y
         xyz(3, 1) = atom(gto(ii)%origin)%coordinate_z

      j_function_loop: do ij = shell_start(iblocks(2)), shell_end(iblocks(2))
         ang(2)    =      gto(ij)%lvalue
         off(2)    =      gto(ij)%offset
         deg(2)    =      gto(ij)%cdegen
         ex(2)     =      gto(ij)%ex(1)
         coef(2)   =      gto(ij)%coefficient(1)
         cent(2)   =      gto(ij)%center
         xyz(1, 2) = atom(gto(ij)%origin)%coordinate_x
         xyz(2, 2) = atom(gto(ij)%origin)%coordinate_y
         xyz(3, 2) = atom(gto(ij)%origin)%coordinate_z

      k_function_loop: do ik = shell_start(iblocks(3)), shell_end(iblocks(3))
         ang(3)    =      gto(ik)%lvalue
         off(3)    =      gto(ik)%offset
         deg(3)    =      gto(ik)%cdegen
         ex(3)     =      gto(ik)%ex(1)
         coef(3)   =      gto(ik)%coefficient(1)
         cent(3)   =      gto(ik)%center
         xyz(1, 3) = atom(gto(ik)%origin)%coordinate_x
         xyz(2, 3) = atom(gto(ik)%origin)%coordinate_y
         xyz(3, 3) = atom(gto(ik)%origin)%coordinate_z

      l_function_loop: do il = shell_start(iblocks(4)), shell_end(iblocks(4))
         ang(4)    =      gto(il)%lvalue
         off(4)    =      gto(il)%offset
         deg(4)    =      gto(il)%cdegen
         ex(4)     =      gto(il)%ex(1)
         coef(4)   =      gto(il)%coefficient(1)
         cent(4)   =      gto(il)%center
         xyz(1, 4) = atom(gto(il)%origin)%coordinate_x
         xyz(2, 4) = atom(gto(il)%origin)%coordinate_y
         xyz(3, 4) = atom(gto(il)%origin)%coordinate_z

         select case (order)

         case (0)
!-------------------------------------------------------------------------------
!           order 0
!-------------------------------------------------------------------------------

            call get_integrals(gint, ang, ex, coef, xyz)
            call contract_integrals(deg,     &
                                    deg,     &
                                    off,     &
                                    1, 2, 3, 4, &
                                    gint,    &
                                    ndim,    &
                                    dmat,    &
                                    gmat,    &
                                    get_ave, &
                                    ave,     &
                                    1.0d0)

         case (1)
!-------------------------------------------------------------------------------
!           order 1
!-------------------------------------------------------------------------------

            icent_loop: do icent = icent_start, icent_end

!              ijkl is in memory (k, l, i, j)
!              we will always differentiate on the slowest index

               !fixme: if ave and unp D, avoid multiple calls and scale by 4.0

               call first_order(1, 2, 3, 4, &
                                ixyz_start, &
                                ixyz_end,   &
                                deg,        &
                                off,        &
                                ang,        &
                                coef,       &
                                ex,         &
                                xyz,        &
                                cent,       &
                                g_u,        &
                                g_d,        &
                                gint,       &
                                dmat,       &
                                gmat,       &
                                ave,        &
                                get_ave,    &
                                icent,      &
                                ndim)

               call first_order(2, 1, 4, 3, &
                                ixyz_start, &
                                ixyz_end,   &
                                deg,        &
                                off,        &
                                ang,        &
                                coef,       &
                                ex,         &
                                xyz,        &
                                cent,       &
                                g_u,        &
                                g_d,        &
                                gint,       &
                                dmat,       &
                                gmat,       &
                                ave,        &
                                get_ave,    &
                                icent,      &
                                ndim)

               call first_order(3, 4, 1, 2, &
                                ixyz_start, &
                                ixyz_end,   &
                                deg,        &
                                off,        &
                                ang,        &
                                coef,       &
                                ex,         &
                                xyz,        &
                                cent,       &
                                g_u,        &
                                g_d,        &
                                gint,       &
                                dmat,       &
                                gmat,       &
                                ave,        &
                                get_ave,    &
                                icent,      &
                                ndim)

               call first_order(4, 3, 2, 1, &
                                ixyz_start, &
                                ixyz_end,   &
                                deg,        &
                                off,        &
                                ang,        &
                                coef,       &
                                ex,         &
                                xyz,        &
                                cent,       &
                                g_u,        &
                                g_d,        &
                                gint,       &
                                dmat,       &
                                gmat,       &
                                ave,        &
                                get_ave,    &
                                icent,      &
                                ndim)

            end do icent_loop

         case default
!-------------------------------------------------------------------------------
!           order too high
!-------------------------------------------------------------------------------
            print *, 'error: order too high in interest interface'
            stop 1

         end select

      end do l_function_loop
      end do k_function_loop
      end do j_function_loop
      end do i_function_loop

   end subroutine

   subroutine get_integrals(gint, &
                            l,    &
                            e,    &
                            c,    &
                            xyz)

      real(8), intent(inout) :: gint(*)
      integer, intent(in)    :: l(4)
      real(8), intent(in)    :: e(4)
      real(8), intent(in)    :: c(4)
      real(8), intent(in)    :: xyz(3, 4)

      real(8) :: f
      integer :: nr_integrals

      !> InteRest interface
      interface
         subroutine interest_eri_basic(f, gint, nr_integrals,  &
                                       li, ei, xi, yi, zi, ci, &
                                       lj, ej, xj, yj, zj, cj, &
                                       lk, ek, xk, yk, zk, ck, &
                                       ll, el, xl, yl, zl, cl)

           real(8), intent(in)  :: f
           real(8), intent(out) :: gint(*)
           integer, intent(out) :: nr_integrals
           integer, intent(in)  :: li, lj, lk, ll
           real(8), intent(in)  :: ei, ej, ek, el
           real(8), intent(in)  :: xi, xj, xk, xl
           real(8), intent(in)  :: yi, yj, yk, yl
           real(8), intent(in)  :: zi, zj, zk, zl
           real(8), intent(in)  :: ci, cj, ck, cl
         end subroutine
      end interface

      f = 1.0d0

      !> call InteRest library routine for a given batch
      call interest_eri_basic(f,                                                   &
                              gint,                                                &
                              nr_integrals,                                        &
                              l(1)+1, e(1), xyz(1, 1), xyz(2, 1), xyz(3, 1), c(1), &
                              l(2)+1, e(2), xyz(1, 2), xyz(2, 2), xyz(3, 2), c(2), &
                              l(3)+1, e(3), xyz(1, 3), xyz(2, 3), xyz(3, 3), c(3), &
                              l(4)+1, e(4), xyz(1, 4), xyz(2, 4), xyz(3, 4), c(4))

   end subroutine

   subroutine contract_integrals(n,              &
                                 m,              &
                                 o,              &
                                 pi, pj, pk, pl, &
                                 gint,           &
                                 ndim,           &
                                 dmat,           &
                                 gmat,           &
                                 get_ave,        &
                                 average,        &
                                 scale_exchange)

      integer, intent(in)    :: n(4)
      integer, intent(in)    :: m(4)
      integer, intent(in)    :: o(4)
      integer, intent(in)    :: pi, pj, pk, pl
      real(8), intent(in)    :: gint(m(3), m(4), m(1), m(2))
      integer, intent(in)    :: ndim
      real(8), intent(in)    :: dmat(ndim, ndim, *)
      real(8), intent(out)   :: gmat(ndim, ndim, *)
      logical, intent(in)    :: get_ave
      real(8), intent(inout) :: average(*)
      real(8), intent(in)    :: scale_exchange

      integer :: i, j, k, l
      integer :: idx(4)
      integer :: bas(4)
      real(8) :: pkl, pkj(4), g

      ! coulomb
      do l = 1, n(4)
         bas(4) = o(4) + l
         do k = 1, n(3)
            bas(3) = o(3) + k
            pkl = 2.0d0*dmat(bas(3), bas(4), 1)
            do j = 1, n(2)
               bas(2) = o(2) + j
               do i = 1, n(1)
                  bas(1) = o(1) + i
                  idx(1) = i
                  idx(2) = j
                  idx(3) = k
                  idx(4) = l
                  g = gint(idx(pk), idx(pl), idx(pi), idx(pj))
                  if (get_ave) then
                     average(1) = average(1) + gmat(bas(1), bas(2), 1)*g*pkl
                  else
                     gmat(bas(1), bas(2), 1) = gmat(bas(1), bas(2), 1) + g*pkl
                  end if
               end do
            end do
         end do
      end do

      ! exchange
      if (dabs(scale_exchange) < 1.0d-20) return
      do l = 1, n(4)
         bas(4) = o(4) + l
         do k = 1, n(3)
            bas(3) = o(3) + k
            do j = 1, n(2)
               bas(2) = o(2) + j
               pkj(1) = dmat(bas(2), bas(3), 1)
#ifdef PRG_DIRAC
               pkj(2) = dmat(bas(2), bas(3), 2)
               pkj(3) = dmat(bas(2), bas(3), 3)
               pkj(4) = dmat(bas(2), bas(3), 4)
#endif
               pkj = scale_exchange*pkj
               do i = 1, n(1)
                  bas(1) = o(1) + i
                  idx(1) = i
                  idx(2) = j
                  idx(3) = k
                  idx(4) = l
                  g = gint(idx(pk), idx(pl), idx(pi), idx(pj))
                  if (get_ave) then
                     average(1) = average(1) - gmat(bas(1), bas(4), 1)*g*pkj(1)
#ifdef PRG_DIRAC
                     average(1) = average(1) - gmat(bas(1), bas(4), 2)*g*pkj(2)
                     average(1) = average(1) - gmat(bas(1), bas(4), 3)*g*pkj(3)
                     average(1) = average(1) - gmat(bas(1), bas(4), 4)*g*pkj(4)
#endif
                  else
                     gmat(bas(1), bas(4), 1) = gmat(bas(1), bas(4), 1) - g*pkj(1)
#ifdef PRG_DIRAC
                     gmat(bas(1), bas(4), 2) = gmat(bas(1), bas(4), 2) - g*pkj(2)
                     gmat(bas(1), bas(4), 3) = gmat(bas(1), bas(4), 3) - g*pkj(3)
                     gmat(bas(1), bas(4), 4) = gmat(bas(1), bas(4), 4) - g*pkj(4)
#endif
                  end if
               end do
            end do
         end do
      end do

  end subroutine

   subroutine init_arrays()

!     integer :: il, ia, ib, ii, ij, ik, ip
      integer :: i, j, k, l, m
      integer :: ndim

      ndim = (maxl+1)*(maxl+2)/2

      ! construct up matrix

      if (allocated(imat_u)) deallocate(imat_u)
      allocate(imat_u(ndim, ndim))
      imat_u = 1

      ! first row is just like this:
      ! 1 2 3 4 5 6 ...
      do i = 1, ndim
         imat_u(i, 1) = i
      end do

      ! construct rows 2 to ndim
      ! they look like this:
      ! 1 2 2 3 3 3 4 4 4 4 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! 1 2 2 3 3 3 4 4 4 4 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! 1 2 2 3 3 3 4 4 4 4 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! 1 1 1 1 1 1 1 1 1 1 ...
      ! ...
      m = 2
      do i = 1, maxl
         l = 0
         do j = 1, maxl+1
            do k = 1, j
               l = l + 1
               imat_u(l, m) = j
            end do
         end do
         m = m + i + 1
      end do

      ! construct rows 2 to ndim
      ! by adding row N to row N-1
      do i = 2, ndim
         do j = 1, ndim
            imat_u(j, i) = imat_u(j, i-1) + imat_u(j, i)
         end do
      end do

      ! construct down matrices

      if (allocated(imat_d)) deallocate(imat_d)
      if (allocated(imat_f)) deallocate(imat_f)
      allocate(imat_d(ndim, 3, maxl))
      allocate(imat_f(ndim, 3, maxl))
      imat_d = 0
      imat_f = 0.0d0

      do i = 1, maxl
         do k = 1, i*(i+1)/2
            imat_d(k, 1, i) = k
         end do
         k = 0
         l = 0
         do m = 1, i
            do j = 1, m
               k = k + 1
               l = l + 1
               imat_d(l+1, 2, i) = k
               imat_d(l+2, 3, i) = k
               imat_f(k,   1, i) = real(i - m + 1)
               imat_f(l+1, 2, i) = real(m - j + 1)
               imat_f(l+2, 3, i) = real(j)
            end do
            l = l + 1
         end do
      end do

   end subroutine

   subroutine first_order(pi, pj, pk, pl, &
                          ixyz_start,     &
                          ixyz_end,       &
                          deg,            &
                          off,            &
                          ang,            &
                          coef,           &
                          ex,             &
                          xyz,            &
                          cent,           &
                          g_u,            &
                          g_d,            &
                          gint,           &
                          dmat,           &
                          gmat,           &
                          ave,            &
                          get_ave,        &
                          icent,          &
                          ndim)

      integer :: pi, pj, pk, pl
      integer :: ixyz_start
      integer :: ixyz_end
      integer :: deg(4)
      integer :: off(4)
      integer :: ang(4)
      real(8) :: coef(4)
      real(8) :: ex(4)
      real(8) :: xyz(3, 4)
      integer :: cent(4)
      real(8) :: g_u(*)
      real(8) :: g_d(*)
      real(8) :: gint(*)
      real(8) :: dmat(*)
      real(8) :: gmat(*)
      real(8) :: ave(*)
      logical :: get_ave
      integer :: icent
      integer :: ndim

      integer :: ixyz
      integer :: nr_elements
      integer :: icoor
      integer :: ideg
      integer :: i
      real(8) :: f

      if (cent(pj) /= icent) return

      call get_integrals(g_u,                                                        &
                         (/     ang(pi),      ang(pj)+1,    ang(pk),      ang(pl)/), &
                         (/      ex(pi),       ex(pj),       ex(pk),       ex(pl)/), &
                         (/    coef(pi),     coef(pj),     coef(pk),     coef(pl)/), &
                         (/xyz(1:3, pi), xyz(1:3, pj), xyz(1:3, pk), xyz(1:3, pl)/))
      if (ang(pj) > 0) then
         call get_integrals(g_d,                                                        &
                            (/     ang(pi),      ang(pj)-1,    ang(pk),      ang(pl)/), &
                            (/      ex(pi),       ex(pj),       ex(pk),       ex(pl)/), &
                            (/    coef(pi),     coef(pj),     coef(pk),     coef(pl)/), &
                            (/xyz(1:3, pi), xyz(1:3, pj), xyz(1:3, pk), xyz(1:3, pl)/))
      end if

      nr_elements = deg(pk)*deg(pl)*deg(pi)

      do ixyz = ixyz_start, ixyz_end

         do ideg = 1, deg(pj)
            do i = 1, nr_elements
               gint(i + nr_elements*(ideg - 1)) = 2.0d0*ex(pj)*g_u(i + nr_elements*(imat_u(ideg, ixyz) - 1))
            end do
         end do

         if (ang(pj) > 0) then
            do ideg = 1, deg(pj)
               f = imat_f(ideg, ixyz, ang(pj))
               if (f > 0.0d0) then
                  do i = 1, nr_elements
                     gint(i + nr_elements*(ideg - 1)) = gint(i + nr_elements*(ideg - 1)) &
                                                      - f*g_d(i + nr_elements*(imat_d(ideg, ixyz, ang(pj)) - 1))
                  end do
               end if
            end do
         end if

         if (get_ave) then
            icoor = (icent-1)*3 + ixyz
         end if

         call contract_integrals(deg, &
                                 (/deg(pi), deg(pj), deg(pk), deg(pl)/), &
                                 off, &
                                 pi, pj, pk, pl, &
                                 gint,                                   &
                                 ndim,                                   &
                                 dmat,                                   &
                                 gmat,                                   &
                                 get_ave,                                &
                                 ave(icoor),                             &
                                 1.0d0)

      end do

   end subroutine

end module
