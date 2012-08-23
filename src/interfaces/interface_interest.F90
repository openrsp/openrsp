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
   integer, allocatable :: ic_to_ijk(:, :, :)
   integer, allocatable :: ijk_to_ic(:, :, :)
   integer, allocatable :: ijk(:, :, :)
   integer, allocatable :: ijk_temp(:, :, :)
   real(8), allocatable :: prefactor(:, :)
   integer, allocatable :: cdeg(:)

   integer :: nr_centers
      integer, parameter :: max_nr_integrals = 194481 !fixme hardcoded

contains

#ifdef VAR_LSDALTON
   subroutine initialize_interest_eri_diff(nplrg,nucind,nlarge,natom,nlrgsh,&
        & nsmlsh,charge,cord,gnuexp,priexp,nrco,nhkt,priccf,ncent)
     implicit none
     integer,intent(in) :: nplrg,nucind,nlarge,natom,nlrgsh,nsmlsh,nhkt(natom)
     real(8),intent(in) :: charge(natom),cord(3,natom),gnuexp(natom)
     real(8),intent(in) :: priexp(natom),priccf(:,:)
     integer,intent(in) :: nrco(:),ncent(natom)
!
     integer        :: i,j,ij,ier
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
         atom(i) = type_atom(charge(i),cord(1,i),cord(2,i),cord(3,i),gnuexp(i))
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
           gto(i) = type_gto(ij,0,nhkt(i)-1,ncent(i),(2*nhkt(i)-1), &
                & (nhkt(i)*(nhkt(i)+1))/2, priexp(i), priccf(i,j), ncent(i))
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
#else
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
#endif

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

      integer, parameter :: max_ave_length   = 100    !fixme hardcoded

      real(8) :: gint(max_nr_integrals, 100) !fixme hardcoded
      real(8) :: gint_temp(max_nr_integrals)
      real(8) :: ex(4), coef(4), xyz(3, 4)
      integer :: ang(4), deg(4), off(4)
      real(8) :: f

      integer :: ii, ij, ik, il
      integer :: cent(4)
      integer :: icent, ixyz
      integer :: nr_integrals
      integer :: icent_start, icent_end
      integer :: ixyz_start,  ixyz_end
      integer :: icoor
      integer :: ifun

      ! make sure it is initialized
      ! it does not cost anything
      ! routine will return if already initialized
#ifdef VAR_LSDALTON
      STOP 'LSDALTON require more input arguments call to initialize_interest_eri_diff'
!      call initialize_interest_eri_diff()
#else
      call initialize_interest_eri_diff()
#endif

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

            icent_loop: do icent = icent_start, icent_end

               !fixme: if ave and unp D, avoid multiple calls and scale by 4.0
               do ifun = 1, 4
                  if (cent(ifun) /= icent) cycle

                  if (get_ave) then
                     do icoor = 1, 3*nr_centers
                        gint(1:deg(1)*deg(2)*deg(3)*deg(4), icoor) = 0.0d0
                     end do
                  else
                     gint(1:deg(1)*deg(2)*deg(3)*deg(4), 1) = 0.0d0
                  end if

                  call first_order(ifun, &
                                   ixyz_start, &
                                   ixyz_end,   &
                                   deg,        &
                                   off,        &
                                   ang,        &
                                   coef,       &
                                   ex,         &
                                   xyz,        &
                                   cent,       &
                                   gint_temp,        &
                                   gint,       &
                                   dmat,       &
                                   gmat,       &
                                   ave,        &
                                   get_ave,    &
                                   icent,      &
                                   ndim)

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
                     call contract_integrals(deg,        &
                                             off,        &
                                             gint(1, 1), &
                                             ndim,       &
                                             dmat,       &
                                             gmat,       &
                                             get_ave,    &
                                             ave,        &
                                             1.0d0)
                  end if
               end do

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
                                 o,              &
                                 gint,           &
                                 ndim,           &
                                 dmat,           &
                                 gmat,           &
                                 get_ave,        &
                                 average,        &
                                 scale_exchange)

      integer, intent(in)    :: n(4)
      integer, intent(in)    :: o(4)
      real(8), intent(in)    :: gint(n(3), n(4), n(1), n(2))
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
                  g = gint(k, l, i, j)
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
                  g = gint(k, l, i, j)
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
      integer :: i, j, k, l, m, n
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

      if (allocated(ic_to_ijk)) deallocate(ic_to_ijk)
      if (allocated(ijk_to_ic)) deallocate(ijk_to_ic)
      if (allocated(ijk))       deallocate(ijk)
      if (allocated(ijk_temp))     deallocate(ijk_temp)
      if (allocated(prefactor))    deallocate(prefactor)
      if (allocated(cdeg))    deallocate(cdeg)

      allocate(ic_to_ijk(ndim, 3, 0:maxl))
      allocate(ijk_to_ic(0:maxl, 0:maxl, 0:maxl))
      allocate(ijk(ndim, 3, 4))
      allocate(ijk_temp(ndim, 3, 4))
      allocate(prefactor(ndim, 4))
      allocate(cdeg(0:maxl))

      ic_to_ijk = 0
      ijk_to_ic = 0

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

      do i = 0, maxl
         k = 0
         l = 0
         do m = 1, i
            do j = 1, m
               k = k + 1
               l = l + 1
               ic_to_ijk(k,   1, i) = i - m + 1
               ic_to_ijk(l+1, 2, i) = m - j + 1
               ic_to_ijk(l+2, 3, i) = j
               ijk_to_ic(i - m, m - j, j - 1) = k
            end do
            l = l + 1
         end do
         cdeg(i) = (i+1)*(i+2)/2
      end do

!     do l = 0, maxl
!        print *, 'raboof l', l
!        n = (l+1)*(l+2)/2
!        do k = 1, n
!           print *, ic_to_ijk(k, 1:3, l), ijk_to_ic(ic_to_ijk(k, 1, l), ic_to_ijk(k, 2, l), ic_to_ijk(k, 3, l))
!        end do
!     end do

   end subroutine

   subroutine first_order(ifd, &
                          ixyz_start,     &
                          ixyz_end,       &
                          deg,            &
                          off,            &
                          ang,            &
                          coef,           &
                          ex,             &
                          xyz,            &
                          cent,           &
                          gint_temp,            &
                          gint,           &
                          dmat,           &
                          gmat,           &
                          ave,            &
                          get_ave,        &
                          icent,          &
                          ndim)

      integer :: ifd
      integer :: ixyz_start
      integer :: ixyz_end
      integer :: deg(4)
      integer :: off(4)
      integer :: ang(4)
      real(8) :: coef(4)
      real(8) :: ex(4)
      real(8) :: xyz(3, 4)
      integer :: cent(4)
      real(8) :: gint_temp(*)
      real(8) :: gint(max_nr_integrals, *)
      real(8) :: dmat(*)
      real(8) :: gmat(*)
      real(8) :: ave(*)
      logical :: get_ave
      integer :: icent
      integer :: ndim

      integer :: ixyz, jxyz
      integer :: icoor
      integer :: ideg
      integer :: ifun, jfun, kfun, lfun
      integer :: i
      integer :: ang_temp(4)

      do ifun = 1, 4
         do ixyz = 1, 3
            ijk(1:deg(ifun), ixyz, ifun) = ic_to_ijk(1:deg(ifun), ixyz, ang(ifun))
         end do
      end do

!     up contribution

      ang_temp = ang
      ang_temp(ifd) = ang_temp(ifd) + 1
      call get_integrals(gint_temp, ang_temp, ex, coef, xyz)

      do ifun = 1, 4
         prefactor(1:deg(ifun), ifun) = 1.0d0
      end do
      do ideg = 1, deg(ifd)
         prefactor(ideg, ifd) = ex(ifd)*2.0d0
      end do

      do ixyz = ixyz_start, ixyz_end

         if (get_ave) then
            icoor = (icent-1)*3 + ixyz
         else
            icoor = 1
         end if

         do ifun = 1, 4
            do jxyz = 1, 3
               ijk_temp(1:deg(ifun), jxyz, ifun) = ijk(1:deg(ifun), jxyz, ifun)
            end do
         end do
         do ideg = 1, deg(ifd)
            ijk_temp(ideg, ixyz, ifd) = ijk_temp(ideg, ixyz, ifd) + 1
         end do

         call add_integrals(gint(1, icoor), &
                            gint_temp,      &
                            deg,            &
                            ang_temp,       &
                            ijk_temp,       &
                            prefactor)

      end do

!     down contribution

      ang_temp = ang
      if (ang_temp(ifd) > 0) then
         ang_temp(ifd) = ang_temp(ifd) - 1
         call get_integrals(gint_temp, ang_temp, ex, coef, xyz)

         do ixyz = ixyz_start, ixyz_end

            if (get_ave) then
               icoor = (icent-1)*3 + ixyz
            else
               icoor = 1
            end if

            do ifun = 1, 4
               prefactor(1:deg(ifun), ifun) = 1.0d0
            end do

            do ifun = 1, 4
               do jxyz = 1, 3
                  ijk_temp(1:deg(ifun), jxyz, ifun) = ijk(1:deg(ifun), jxyz, ifun)
               end do
            end do
            do ideg = 1, deg(ifd)
               if (ijk_temp(ideg, ixyz, ifd) > 0) then
                  prefactor(ideg, ifd) = -real(ijk_temp(ideg, ixyz, ifd))
               else
                  prefactor(ideg, ifd) = 0.0d0
               end if
               ijk_temp(ideg, ixyz, ifd) = ijk_temp(ideg, ixyz, ifd) - 1
            end do

            call add_integrals(gint(1, icoor), &
                               gint_temp,      &
                               deg,            &
                               ang_temp,       &
                               ijk_temp,       &
                               prefactor)

         end do
      end if

   end subroutine

   subroutine add_integrals(g_out,   &
                            g_in,    &
                            deg_out, &
                            ang_in,  &
                            ijk_in,  &
                            prefactor_in)

!     ijkl is in memory (k, l, i, j)

      real(8), intent(inout) :: g_out(*)
      real(8), intent(in)    :: g_in(*)
      integer, intent(in)    :: deg_out(4)
      integer, intent(in)    :: ang_in(4)
      integer, intent(in)    :: ijk_in(:, :, :)
      real(8), intent(in)    :: prefactor_in(:, :)

      integer :: pr,   pw
      integer :: ifun, jfun, kfun, lfun
      integer :: ir,   jr,   kr,   lr
      integer :: iw,   jw,   kw,   lw
      real(8) :: fi,   fj,   fk,   fl
      real(8) :: ldri, ldrj, ldrk, ldrl
      real(8) :: ldwi, ldwj, ldwk, ldwl
      real(8) :: ofri, ofrj, ofrk, ofrl
      real(8) :: ofwi, ofwj, ofwk, ofwl

      ifun = 1
      jfun = 2
      kfun = 3
      lfun = 4

      ldrk = 0
      ldrl = cdeg(ang_in(kfun))
      ldri = cdeg(ang_in(lfun))*ldrl
      ldrj = cdeg(ang_in(ifun))*ldri

      ldwk = 0
      ldwl = deg_out(kfun)
      ldwi = deg_out(lfun)*ldwl
      ldwj = deg_out(ifun)*ldwi

      do jw = 1, deg_out(jfun)
         fj = prefactor_in(jw, jfun)
         if (fj == 0.0d0) cycle
         jr = ijk_to_ic(ijk_in(jw, 1, jfun), ijk_in(jw, 2, jfun), ijk_in(jw, 3, jfun))
         ofrj = (jr - 1)*ldrj
         ofwj = (jw - 1)*ldwj
         do iw = 1, deg_out(ifun)
            fi = prefactor_in(iw, ifun)*fj
            if (fi == 0.0d0) cycle
            ir = ijk_to_ic(ijk_in(iw, 1, ifun), ijk_in(iw, 2, ifun), ijk_in(iw, 3, ifun))
            ofri = (ir - 1)*ldri + ofrj
            ofwi = (iw - 1)*ldwi + ofwj
            do lw = 1, deg_out(lfun)
               fl = prefactor_in(lw, lfun)*fi
               if (fl == 0.0d0) cycle
               lr = ijk_to_ic(ijk_in(lw, 1, lfun), ijk_in(lw, 2, lfun), ijk_in(lw, 3, lfun))
               ofrl = (lr - 1)*ldrl + ofri
               ofwl = (lw - 1)*ldwl + ofwi
               do kw = 1, deg_out(kfun)
                  fk = prefactor_in(kw, kfun)*fl
                  if (fk == 0.0d0) cycle
                  kr = ijk_to_ic(ijk_in(kw, 1, kfun), ijk_in(kw, 2, kfun), ijk_in(kw, 3, kfun))
                  pr = ofrl + kr
                  pw = ofwl + kw
                  g_out(pw) = g_out(pw) + fk*g_in(pr)
               end do
            end do
         end do
      end do

   end subroutine

end module
