module interface_interest

   implicit none

   public interest_eri_diff

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
   end type

   integer, save :: shell_start(2) = 0
   integer, save :: shell_end(2) = 0

   type(type_atom), allocatable :: atom(:)
   type(type_gto),  allocatable :: gto(:)

   logical :: interface_is_initialized = .false.

contains

   subroutine initialize_interest_eri_diff()

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#ifdef PRG_DIRAC
#include "dcbdhf.h"
#include "nucdata.h"
#endif

      integer        :: i
      integer        :: j
      integer        :: ij
      integer        :: ier
      type(type_gto) :: tmpgto

      if (interface_is_initialized) return

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
                             nhkt(i),                 &
                             ncent(i),                &
                             (2*nhkt(i)-1),           &
                             (nhkt(i)*(nhkt(i)+1))/2, &
                             priexp(i),               &
                             priccf(i,j))
           !> note: cartesian indexation
           ij = ij + (nhkt(i)*(nhkt(i)+1))/2
        end do
      end do

      do i = 1, size(gto)-1
         gto(i+1)%offset = gto(i)%offset + gto(i)%cdegen
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

   subroutine interest_eri_diff(ndim, Gmat, Pmat)

      integer, intent(in)  :: ndim
      real(8), intent(in)  :: Pmat(ndim, ndim, 4)
      real(8), intent(out) :: Gmat(ndim, ndim, 4)

      call interest_eri_diff_block(ndim, Gmat, Pmat, (/1, 1, 1, 1/))
#ifdef PRG_DIRAC
      call interest_eri_diff_block(ndim, Gmat, Pmat, (/1, 1, 2, 2/))
      call interest_eri_diff_block(ndim, Gmat, Pmat, (/2, 2, 1, 1/))
#endif

   end subroutine

   subroutine interest_eri_diff_block(ndim, Gmat, Pmat, iblocks)

      integer, intent(in)  :: ndim
      real(8), intent(in)  :: Pmat(ndim, ndim, 4)
      real(8), intent(out) :: Gmat(ndim, ndim, 4)
      integer, intent(in)  :: iblocks(4)

      integer, parameter :: max_nr_integrals = 194481 !fixme hardcoded

      real(8) :: gout(max_nr_integrals)
      real(8) :: f
      real(8) :: e(4)
      real(8) :: c(4)
      real(8) :: x(4)
      real(8) :: y(4)
      real(8) :: z(4)
      integer :: l(4)
      integer :: n(4)
      integer :: o(4)

      integer :: ii, ij, ik, il
      integer :: nr_integrals

      !> InteRest interface
      interface
         subroutine interest_eri_basic(factor, fint, nr_integrals, &
                                       li, ei, xi, yi, zi, ci,     &
                                       lj, ej, xj, yj, zj, cj,     &
                                       lk, ek, xk, yk, zk, ck,     &
                                       ll, el, xl, yl, zl, cl)

           real(8), intent(in)  :: factor
           real(8), intent(out) :: fint(*)
           integer, intent(out) :: nr_integrals
           integer, intent(in)  :: li, lj, lk, ll
           real(8), intent(in)  :: ei, ej, ek, el
           real(8), intent(in)  :: xi, xj, xk, xl
           real(8), intent(in)  :: yi, yj, yk, yl
           real(8), intent(in)  :: zi, zj, zk, zl
           real(8), intent(in)  :: ci, cj, ck, cl
         end subroutine
      end interface

      ! make sure it is initialized
      ! it does not cost anything
      ! routine will return if already initialized
      call initialize_interest_eri_diff()

      f = 1.0d0

      iloop: do ii = shell_start(iblocks(1)), shell_end(iblocks(1))
         l(1) =      gto(ii)%lvalue
         o(1) =      gto(ii)%offset
         n(1) =      gto(ii)%cdegen
         e(1) =      gto(ii)%ex(1)
         x(1) = atom(gto(ii)%origin)%coordinate_x
         y(1) = atom(gto(ii)%origin)%coordinate_y
         z(1) = atom(gto(ii)%origin)%coordinate_z
         c(1) =      gto(ii)%coefficient(1)

         jloop: do ij = shell_start(iblocks(2)), shell_end(iblocks(2))
            l(2) =      gto(ij)%lvalue
            o(2) =      gto(ij)%offset
            n(2) =      gto(ij)%cdegen
            e(2) =      gto(ij)%ex(1)
            x(2) = atom(gto(ij)%origin)%coordinate_x
            y(2) = atom(gto(ij)%origin)%coordinate_y
            z(2) = atom(gto(ij)%origin)%coordinate_z
            c(2) =      gto(ij)%coefficient(1)

            kloop: do ik = shell_start(iblocks(3)), shell_end(iblocks(3))
               l(3) =      gto(ik)%lvalue
               o(3) =      gto(ik)%offset
               n(3) =      gto(ik)%cdegen
               e(3) =      gto(ik)%ex(1)
               x(3) = atom(gto(ik)%origin)%coordinate_x
               y(3) = atom(gto(ik)%origin)%coordinate_y
               z(3) = atom(gto(ik)%origin)%coordinate_z
               c(3) =      gto(ik)%coefficient(1)

               lloop: do il = shell_start(iblocks(4)), shell_end(iblocks(4))
                  l(4) =      gto(il)%lvalue
                  o(4) =      gto(il)%offset
                  n(4) =      gto(il)%cdegen
                  e(4) =      gto(il)%ex(1)
                  x(4) = atom(gto(il)%origin)%coordinate_x
                  y(4) = atom(gto(il)%origin)%coordinate_y
                  z(4) = atom(gto(il)%origin)%coordinate_z
                  c(4) =      gto(il)%coefficient(1)

                  ! f is scaling constant
                  ! gout  is output (integrals)
                  ! [ij|kl] [electron1|electron2]
                  ! gout = (ncc(k), ncc(l), ncc(i), ncc(j))
                  ! example: [sp|df]: [6, 10, 1, 3] this is the layout in mem
                  ! limitation: up to h functions (incl)

                  !> call InteRest library routine for a given batch
                  call interest_eri_basic(f,                                  &
                                          gout,                               &
                                          nr_integrals,                       &
                                          l(1), e(1), x(1), y(1), z(1), c(1), &
                                          l(2), e(2), x(2), y(2), z(2), c(2), &
                                          l(3), e(3), x(3), y(3), z(3), c(3), &
                                          l(4), e(4), x(4), y(4), z(4), c(4) )

                  call process_dG(n,    &
                                  o,    &
                                  gout, &
                                  ndim, &
                                  Pmat, &
                                  Gmat, &
                                  1.0d0)

               end do lloop
            end do kloop
         end do jloop
      end do iloop

   end subroutine

   subroutine process_dG(n,    &
                         o,    &
                         gout, &
                         ndim, &
                         Pmat, &
                         Gmat, &
                         scale_exchange)

      integer, intent(in)  :: n(4)
      integer, intent(in)  :: o(4)
      real(8), intent(in)  :: gout(n(3), n(4), n(1), n(2))
      integer, intent(in)  :: ndim
      real(8), intent(in)  :: Pmat(ndim, ndim, *)
      real(8), intent(out) :: Gmat(ndim, ndim, *)
      real(8), intent(in)  :: scale_exchange

      integer :: i, j, k, l
      integer :: bas(4)
      real(8) :: pkl, pkj(4)

      ! coulomb
      do l = 1, n(4)
         bas(4) = o(4) + l
         do k = 1, n(3)
            bas(3) = o(3) + k
            pkl = 2.0d0*Pmat(bas(3), bas(4), 1)
            do j = 1, n(2)
               bas(2) = o(2) + j
               do i = 1, n(1)
                  bas(1) = o(1) + i
                  Gmat(bas(1), bas(2), 1) = Gmat(bas(1), bas(2), 1) &
                                          + gout(k, l, i, j)*pkl
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
               pkj(1) = scale_exchange*Pmat(bas(2), bas(3), 1)
               pkj(2) = scale_exchange*Pmat(bas(2), bas(3), 2)
               pkj(3) = scale_exchange*Pmat(bas(2), bas(3), 3)
               pkj(4) = scale_exchange*Pmat(bas(2), bas(3), 4)
               do i = 1, n(1)
                  bas(1) = o(1) + i
                  Gmat(bas(1), bas(4), 1) = Gmat(bas(1), bas(4), 1) &
                                          - gout(k, l, i, j)*pkj(1)
#ifdef PRG_DIRAC
                  Gmat(bas(1), bas(4), 2) = Gmat(bas(1), bas(4), 2) &
                                          - gout(k, l, i, j)*pkj(2)
                  Gmat(bas(1), bas(4), 3) = Gmat(bas(1), bas(4), 3) &
                                          - gout(k, l, i, j)*pkj(3)
                  Gmat(bas(1), bas(4), 4) = Gmat(bas(1), bas(4), 4) &
                                          - gout(k, l, i, j)*pkj(4)
#endif
               end do
            end do
         end do
      end do

  end subroutine

end module
