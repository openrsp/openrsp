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

   integer, allocatable :: cdeg(:)
   integer, allocatable :: ijk_to_ic(:, :, :)
   integer, allocatable :: ic_to_ijk(:, :, :)

   integer, parameter   :: maxl = 3
   integer, allocatable :: imat_u(:, :)
   integer, allocatable :: imat_d(:, :, :)
   integer, allocatable :: imat_f(:, :, :)

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

      real(8) :: gint(max_nr_integrals, max_ave_length)
      real(8) :: g_u(max_nr_integrals)
      real(8) :: g_d(max_nr_integrals)
      real(8) :: ex(4), coef(4), xyz(3, 4), cent(4)
      integer :: ang(4), ang_temp(4), deg(4), off(4)
      integer :: ciw(4)
      integer :: cr(4), cir(4)
      integer :: ci, cj, ck, cl
      integer :: ir, iw

      integer :: ii, ij, ik, il
      integer :: ifun, ic
      integer :: icent, ixyz
      integer :: ijk(3), ijk_u(3), ijk_d(3)
      integer :: nr_integrals
      integer :: icoor
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
                                    off,     &
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
            gint(1:deg(1)*deg(2)*deg(3)*deg(4), 1:3*nr_centers) = 0.0d0

            icent_loop: do icent = icent_start, icent_end
               ifun_loop: do ifun = 1, 4
                  if (cent(ifun) == icent) then

                     ang_temp = ang
                     ang_temp(ifun) = ang_temp(ifun) + 1
                     call get_integrals(g_u, ang_temp, ex, coef, xyz)
                     if (ang(ifun) > 0) then
                        ang_temp = ang
                        ang_temp(ifun) = ang_temp(ifun) - 1
                        call get_integrals(g_d, ang_temp, ex, coef, xyz)
                     end if

                     ! the quintuple loop below can be optimized
                     ! but i will do that later after having
                     ! higher order derivatives
                     ! now i prefer compact and readable code
                     iw = 0
                     do cj = 1, deg(2)
                        ciw(2) = cj
                        do ci = 1, deg(1)
                           ciw(1) = ci
                           do cl = 1, deg(4)
                              ciw(4) = cl
                              do ck = 1, deg(3)
                                 ciw(3) = ck
                                 iw = iw + 1

                                 ijk = get_ijk(ang(ifun), ciw(ifun))

                                 do ixyz = ixyz_start, ixyz_end

                                    if (get_ave) then
                                       icoor = (icent-1)*3 + ixyz
                                    else
                                       icoor = 1
                                    end if

                                    ijk_u       = ijk
                                    ijk_u(ixyz) = ijk_u(ixyz) + 1
                                    ijk_d       = ijk
                                    ijk_d(ixyz) = ijk_d(ixyz) - 1

                                    cr        = deg
                                    cr(ifun)  = cdeg(ang(ifun)+1)
                                    cir       = ciw
                                    cir(ifun) = get_ic(ijk_u)

                                    ir = (cir(2)-1)*cr(3)*cr(4)*cr(1) &
                                       + (cir(1)-1)*cr(3)*cr(4)       &
                                       + (cir(4)-1)*cr(3)             &
                                       +  cir(3)

                                    gint(iw, icoor) = gint(iw, icoor) + g_u(ir)*2.0d0*ex(ifun)

                                    if (ijk_d(ixyz) < 0) cycle

                                    cr        = deg
                                    cr(ifun)  = cdeg(ang(ifun)-1)
                                    cir       = ciw
                                    cir(ifun) = get_ic(ijk_d)

                                    ir = (cir(2)-1)*cr(3)*cr(4)*cr(1) &
                                       + (cir(1)-1)*cr(3)*cr(4)       &
                                       + (cir(4)-1)*cr(3)             &
                                       +  cir(3)

                                    gint(iw, icoor) = gint(iw, icoor) - g_d(ir)*ijk(ixyz)

                                 end do
                              end do
                           end do
                        end do
                     end do

                  end if
               end do ifun_loop
            end do icent_loop

            if (get_ave) then
               do icoor = 1, 3*nr_centers
                  call contract_integrals(deg,            &
                                          off,            &
                                          gint(1, icoor), &
                                          ndim,           &
                                          dmat,           &
                                          gmat,           &
                                          get_ave,        &
                                          ave(icoor),     &
                                          1.0d0)
               end do
            else
               call contract_integrals(deg,     &
                                       off,     &
                                       gint,    &
                                       ndim,    &
                                       dmat,    &
                                       gmat,    &
                                       get_ave, &
                                       ave,     &
                                       1.0d0)
            end if

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

   subroutine contract_integrals(n,        &
                                 o,        &
                                 gint,     &
                                 ndim,     &
                                 dmat,     &
                                 gmat,     &
                                 get_ave,  &
                                 average,  &
                                 scale_exchange)

      integer, intent(in)    :: n(4)
      integer, intent(in)    :: o(4)
      real(8), intent(in)    :: gint(n(3), n(4), n(1), n(2), *)
      integer, intent(in)    :: ndim
      real(8), intent(in)    :: dmat(ndim, ndim, *)
      real(8), intent(out)   :: gmat(ndim, ndim, *)
      logical, intent(in)    :: get_ave
      real(8), intent(inout) :: average(*)
      real(8), intent(in)    :: scale_exchange

      integer :: i, j, k, l
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
                  g = gint(k, l, i, j, 1)
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
                  g = gint(k, l, i, j, 1)
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

      integer :: il, ia, ib, ii, ij, ik, ip
      integer :: i, j, k, l, m
      integer :: ndim

      if (allocated(ijk_to_ic)) deallocate(ijk_to_ic)
      if (allocated(ic_to_ijk)) deallocate(ic_to_ijk)
      if (allocated(cdeg)) deallocate(cdeg)

      allocate(ijk_to_ic(0:maxl, 0:maxl, 0:maxl))
      allocate(ic_to_ijk(0:maxl, (maxl + 1)*(maxl + 2)/2, 3))
      allocate(cdeg(0:maxl))

      ijk_to_ic = 0
      ic_to_ijk = 0
      cdeg = 0

      do il = 0, maxl
         cdeg(il) = ((il + 1)*(il + 2))/2
         ip = 0
         do ia = 1, il + 1
            do ib = 1, ia
               ip = ip + 1
               ii = il + 1 - ia
               ij = ia - ib
               ik = ib - 1
               ijk_to_ic(ii, ij, ik) = ip
               ic_to_ijk(il, ip, 1)  = ii
               ic_to_ijk(il, ip, 2)  = ij
               ic_to_ijk(il, ip, 3)  = ik
            end do
         end do
      end do

      ndim = (maxl+1)*(maxl+2)/2

      ! construct up matrix

      allocate(imat_u(ndim, ndim))
      imat_u = 1

      do i = 1, ndim
         imat_u(i, 1) = i
      end do

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

      do i = 2, ndim
         do j = 1, ndim
            imat_u(j, i) = imat_u(j, i-1) + imat_u(j, i)
         end do
      end do

      ! construct down matrices

      allocate(imat_d(ndim, 3, maxl))
      allocate(imat_f(ndim, 3, maxl))
      imat_d = 0
      imat_f = 0

      do l = 1, maxl
         do k = 1, l*(l+1)/2
            imat_d(k, 1, l) = k
         end do
      end do

      do i = 1, maxl
         k = 0
         l = 0
         do m = 1, i
            do j = 1, m
               k = k + 1
               l = l + 1
               imat_d(l+1, 2, i) = k
               imat_d(l+2, 3, i) = k
            end do
            l = l + 1
         end do
      end do

      do i = 1, maxl
         l = 0
         do m = 1, i
            do j = 1, m
               l = l + 1
               imat_f(l, 1, i) = i - m + 1
            end do
         end do
      end do

      do i = 1, maxl
         k = 0
         l = 0
         do m = 1, i
            do j = 1, m
               k = k + 1
               l = l + 1
               imat_f(l+1, 2, i) = m - j + 1
               imat_f(l+2, 3, i) = j
            end do
            l = l + 1
         end do
      end do

!     do i = 1, maxl
!        print *, 'raboof down', i
!        write(*, '(a, 10i3)') 'raboof', imat_d(1:ndim, 1, i)
!        write(*, '(a, 10i3)') 'raboof', imat_d(1:ndim, 2, i)
!        write(*, '(a, 10i3)') 'raboof', imat_d(1:ndim, 3, i)
!     end do

!     do i = 1, maxl
!        print *, 'raboof fac', i
!        write(*, '(a, 10i3)') 'raboof', imat_f(1:ndim, 1, i)
!        write(*, '(a, 10i3)') 'raboof', imat_f(1:ndim, 2, i)
!        write(*, '(a, 10i3)') 'raboof', imat_f(1:ndim, 3, i)
!     end do

   end subroutine

   function get_ijk(il, ic)
      integer :: il, ic
      integer :: get_ijk(3)
      get_ijk(1:3) = ic_to_ijk(il, ic, 1:3)
   end function

   function get_ic(ijk)
      integer :: get_ic, ijk(3)
      get_ic = ijk_to_ic(ijk(1), ijk(2), ijk(3))
   end function

end module
