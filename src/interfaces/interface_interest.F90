module interface_interest

   !fixme: * if ave and unp D, avoid loop over ifun and scale by 4.0
   !       * parallel code relies on the fact that common blocks are broadcast
   !         before this module is called, better would be to sync basis set info here

   use openrsp_cfg

   implicit none

   public interest_get_int
   public interest_mpi_wake_up
   public interest_mpi_launch_slave_process

   private

#ifdef VAR_MPI
#include "mpif.h"
#endif

   ! max angular momentum (s = 0, p = 1, ...)
   integer, parameter :: maxl = 7
   ! max nr of integrals ((maxl+1)*(maxl+2)/2)**4
   integer, parameter :: max_nr_integrals = 1679616
   ! 3*3*3
   integer, parameter :: max_ave_length   = 27

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

   integer, allocatable :: ic_to_ijk(:, :, :)
   integer, allocatable :: ijk_to_ic(:, :, :)
   integer, allocatable :: ijk_list(:, :, :)
   real(8), allocatable :: prefactor_list(:, :)
   integer, allocatable :: cdeg(:)

   integer :: nr_centers

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
      call init_permanent_arrays()
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
      integer        :: ierr
      integer        :: proc_rank
      type(type_gto) :: tmpgto

      if (interface_is_initialized) return

#ifdef VAR_MPI
      call mpi_comm_rank(MPI_COMM_WORLD, proc_rank, ierr)
#else
      proc_rank = 0
#endif

      call init_permanent_arrays()

      !> initialize InteRest integral package
      call interest_initialize()

      !> test if the large component basis is uncontracted
      if (nplrg /= nlarge) then
         print *, ' Error: contracted basis set found => stop!'
         print *, ' Note:  uncontracted basis sets are only supported'
         stop 1
      endif

      !> interface molecular data
      allocate(atom(nucind), stat=ierr)
      if (ierr /= 0) stop 'Error in allocation: atom(:)'
      do i = 1, size(atom)
         atom(i) = type_atom(charge(i), &
                             cord(1,i), &
                             cord(2,i), &
                             cord(3,i), &
                             gnuexp(i))
      enddo

      !> interface basis set data
      !> fixme: contracted basis sets
      allocate(gto(nlrgsh+nsmlsh), stat=ierr)
      if (ierr /= 0) stop 'Error in allocation: gto(:)'
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

      if (proc_rank == 0) then
         print *, 'total number of basis function shells:    ', size(gto)
         print *, 'total number of spherical basis functions:', sum(gto(:)%sdegen)
         print *, 'total number of cartesian basis functions:', sum(gto(:)%cdegen)
      end if

      interface_is_initialized = .true.

   end subroutine
#endif

   subroutine interest_mpi_wake_up()

#ifdef VAR_MPI
      integer :: nr_proc
      integer :: ierr

      call mpi_comm_size(MPI_COMM_WORLD, nr_proc, ierr)
      if (nr_proc > 1) then
         call dirac_parctl(8)
      end if
#endif

   end subroutine

   subroutine interest_mpi_launch_slave_process()

      integer :: ndim
      real(8) :: array_dmat(1)
      real(8) :: array_gmat(1)
      integer :: order
      integer :: only_icoor
      integer :: ave_length
      real(8) :: array_ave(1)

      call interest_get_int(ndim, array_dmat, array_gmat, order, only_icoor, ave_length, array_ave)

   end subroutine

   subroutine interest_get_int(ndim, array_dmat, array_gmat, order, only_icoor, ave_length, array_ave)

      integer                      :: ndim
      real(8), target              :: array_dmat(*)
      real(8), target              :: array_gmat(*)
      integer                      :: order
      integer                      :: only_icoor
      integer                      :: ave_length
      real(8), target              :: array_ave(*)

      logical                      :: get_ave
      integer                      :: nr_proc
      integer                      :: proc_rank
      integer                      :: ierr

      real(8), pointer             :: dmat(:)
      real(8), pointer             :: gmat(:)
      real(8), pointer             :: ave(:)

      real(8), allocatable, target :: dmat_container(:)
      real(8), allocatable, target :: gmat_container(:)
      real(8), allocatable, target :: ave_container(:)

#ifdef VAR_MPI
      call mpi_comm_size(MPI_COMM_WORLD, nr_proc, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, proc_rank, ierr)
#else
      nr_proc   = 1
      proc_rank = 0
#endif

#ifdef VAR_MPI
      if (nr_proc > 1) then
         call mpi_bcast(ndim,       1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(order,      1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(only_icoor, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(ave_length, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
      end if
#endif

      nullify(dmat)
      nullify(gmat)
      nullify(ave)

      if (proc_rank == 0) then
         dmat => array_dmat(1:ndim*ndim*4)
         gmat => array_gmat(1:ndim*ndim*4)
         ave  => array_ave(1:ave_length)
      else
         allocate(dmat_container(ndim*ndim*4))
         allocate(gmat_container(ndim*ndim*4))
         gmat_container = 0.0d0
         allocate(ave_container(ave_length))
         dmat => dmat_container
         gmat => gmat_container
         ave  => ave_container
      end if

      if (ave_length > 0) then
         ave(1:ave_length) = 0.0d0
         get_ave = .true.
      else
         get_ave = .false.
      end if

#ifdef VAR_MPI
      if (nr_proc > 1) then
         call mpi_bcast(dmat, ndim*ndim*4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         if (get_ave) then
            call mpi_bcast(gmat, ndim*ndim*4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         end if
         call mpi_bcast(openrsp_cfg_skip_llss, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(openrsp_cfg_skip_ssss, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
#endif

      call interest_eri_diff_block(ndim,           &
                                   dmat,           &
                                   gmat,           &
                                   order,          &
                                   (/1, 1, 1, 1/), &
                                   ave,            &
                                   get_ave,        &
                                   only_icoor)

#ifdef PRG_DIRAC
      if (.not. openrsp_cfg_skip_llss) then
         call interest_eri_diff_block(ndim,           &
                                      dmat,           &
                                      gmat,           &
                                      order,          &
                                      (/1, 1, 2, 2/), &
                                      ave,            &
                                      get_ave,        &
                                      only_icoor)
         call interest_eri_diff_block(ndim,           &
                                      dmat,           &
                                      gmat,           &
                                      order,          &
                                      (/2, 2, 1, 1/), &
                                      ave,            &
                                      get_ave,        &
                                      only_icoor)
      end if

      if (.not. openrsp_cfg_skip_ssss) then
         call interest_eri_diff_block(ndim,           &
                                      dmat,           &
                                      gmat,           &
                                      order,          &
                                      (/2, 2, 2, 2/), &
                                      ave,            &
                                      get_ave,        &
                                      only_icoor)
      end if
#endif

#ifdef VAR_MPI
      if (nr_proc > 1) then
         if (get_ave) then
            if (proc_rank == 0) then
               call mpi_reduce(MPI_IN_PLACE, ave, ave_length, MPI_REAL8, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            else
               call mpi_reduce(ave, MPI_IN_PLACE, ave_length, MPI_REAL8, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            end if
         else
            if (proc_rank == 0) then
               call mpi_reduce(MPI_IN_PLACE, gmat, ndim*ndim*4, MPI_REAL8, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            else
               call mpi_reduce(gmat, MPI_IN_PLACE, ndim*ndim*4, MPI_REAL8, mpi_sum, 0, MPI_COMM_WORLD, ierr)
            end if
         end if
      end if
#endif

      nullify(dmat)
      nullify(gmat)
      nullify(ave)

      if (proc_rank /= 0) then
         deallocate(dmat_container)
         deallocate(gmat_container)
         deallocate(ave_container)
      end if

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

      real(8) :: gint(max_nr_integrals, max_ave_length)
      real(8) :: gint_temp(max_nr_integrals)
      real(8) :: ex(4), coef(4), xyz(3, 4)
      integer :: ang(4), deg(4), off(4)
      real(8) :: f

      integer :: ii, ij, ik, il
      integer :: cent(4)
      integer :: icent, jcent, kcent
      integer :: nr_integrals
      integer :: icent_start, icent_end
      integer :: ixyz_start,  ixyz_end
      integer :: i, j, k, l
      integer :: icoor, jcoor, kcoor, lcoor
      integer :: ixyz,  jxyz,  kxyz,  lxyz
      integer :: ifun,  jfun,  kfun,  lfun
      integer :: ang_temp(4)
      integer :: ideg
      logical :: recalc_integrals
      integer :: nr_proc
      integer :: proc_rank
      integer :: ierr

      ! make sure it is initialized
      ! it does not cost anything
      ! routine will return if already initialized
#ifdef VAR_LSDALTON
      STOP 'LSDALTON require more input arguments call to initialize_interest_eri_diff'
!      call initialize_interest_eri_diff()
#else
      call initialize_interest_eri_diff()
#endif

#ifdef VAR_MPI
      call mpi_comm_size(MPI_COMM_WORLD, nr_proc, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, proc_rank, ierr)
#else
      nr_proc   = 1
      proc_rank = 0
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

         ! this is the parallelization
         ! node 1 will do all quadruples with ii=1
         ! node 2 will do all quadruples with ii=2
         ! node 3 will do all quadruples with ii=3
         ! ...
         if (mod(ii, nr_proc) /= proc_rank) cycle

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

            nr_integrals = deg(1)*deg(2)*deg(3)*deg(4)
            call init_ijk(deg, ang)

            do icent = icent_start, icent_end

               ! zero out integrals
               if (get_ave) then
                  do ixyz = 1, 3
                     do i = 1, nr_integrals
                        gint(i, ixyz) = 0.0d0
                     end do
                  end do
               else
                  do i = 1, nr_integrals
                     gint(i, 1) = 0.0d0
                  end do
               end if

               do ifun = 1, 4
                  if (cent(ifun) /= icent) cycle

                  ! u contribution
                  ang_temp = ang
                  ang_temp(ifun) = ang_temp(ifun) + 1
                  recalc_integrals = .true.
                  do ixyz = ixyz_start, ixyz_end
                     if (get_ave) then
                        k = ixyz
                     else
                        k = 1
                     end if
                     call init_prefactors(deg)
                     call fac_u(deg, ex, ifun, ixyz)
                     call ijk_u(deg, ifun, ixyz)
                     call get_integral_contribution(gint,      &
                                                    gint_temp, &
                                                    deg,       &
                                                    ang_temp,  &
                                                    ex,        &
                                                    coef,      &
                                                    xyz,       &
                                                    k,         &
                                                    recalc_integrals)
                     call ijk_d(deg, ifun, ixyz)
                  end do

                  ! d contribution
                  ang_temp = ang
                  ang_temp(ifun) = ang_temp(ifun) - 1
                  recalc_integrals = .true.
                  do ixyz = ixyz_start, ixyz_end
                     if (get_ave) then
                        k = ixyz
                     else
                        k = 1
                     end if
                     call init_prefactors(deg)
                     call fac_d(deg, ex, ifun, ixyz)
                     call ijk_d(deg, ifun, ixyz)
                     call get_integral_contribution(gint,      &
                                                    gint_temp, &
                                                    deg,       &
                                                    ang_temp,  &
                                                    ex,        &
                                                    coef,      &
                                                    xyz,       &
                                                    k,         &
                                                    recalc_integrals)
                     call ijk_u(deg, ifun, ixyz)
                  end do

               end do

               if (get_ave) then
                  do ixyz = 1, 3
                     k = (icent - 1)*3 + ixyz
                     call contract_integrals(deg,           &
                                             off,           &
                                             gint(1, ixyz), &
                                             ndim,          &
                                             dmat,          &
                                             gmat,          &
                                             get_ave,       &
                                             ave(k),        &
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

            end do

         case (2)
!-------------------------------------------------------------------------------
!           order 2
!-------------------------------------------------------------------------------

            nr_integrals = deg(1)*deg(2)*deg(3)*deg(4)
            call init_ijk(deg, ang)

            do icent = 1, nr_centers
               do jcent = icent, nr_centers

                  ! zero out integrals
                  k = 0
                  do ixyz = 1, 3
                     icoor = (icent-1)*3 + ixyz
                     do jxyz = 1, 3
                        jcoor = (jcent-1)*3 + jxyz
                        if (jcoor >= icoor) then
                           k = k + 1
                           do i = 1, nr_integrals
                              gint(i, k) = 0.0d0
                           end do
                        end if
                     end do
                  end do

                  do ifun = 1, 4
                     if (cent(ifun) /= icent) cycle
                     do jfun = 1, 4
                        if (cent(jfun) /= jcent) cycle

                        ! uu contribution
                        recalc_integrals = .true.
                        ang_temp = ang
                        ang_temp(ifun) = ang_temp(ifun) + 1
                        ang_temp(jfun) = ang_temp(jfun) + 1
                        k = 0
                        do ixyz = 1, 3
                           icoor = (icent-1)*3 + ixyz
                           do jxyz = 1, 3
                              jcoor = (jcent-1)*3 + jxyz
                              if (jcoor >= icoor) then
                                 k = k + 1
                                 call init_prefactors(deg)
                                 call fac_u(deg, ex, ifun, ixyz)
                                 call ijk_u(deg, ifun, ixyz)
                                 call fac_u(deg, ex, jfun, jxyz)
                                 call ijk_u(deg, jfun, jxyz)
                                 call get_integral_contribution(gint,      &
                                                                gint_temp, &
                                                                deg,       &
                                                                ang_temp,  &
                                                                ex,        &
                                                                coef,      &
                                                                xyz,       &
                                                                k,         &
                                                                recalc_integrals)
                                 call ijk_d(deg, ifun, ixyz)
                                 call ijk_d(deg, jfun, jxyz)
                              end if
                           end do
                        end do

                        ! ud contribution
                        recalc_integrals = .true.
                        ang_temp = ang
                        ang_temp(ifun) = ang_temp(ifun) + 1
                        ang_temp(jfun) = ang_temp(jfun) - 1
                        k = 0
                        do ixyz = 1, 3
                           icoor = (icent-1)*3 + ixyz
                           do jxyz = 1, 3
                              jcoor = (jcent-1)*3 + jxyz
                              if (jcoor >= icoor) then
                                 k = k + 1
                                 call init_prefactors(deg)
                                 call fac_u(deg, ex, ifun, ixyz)
                                 call ijk_u(deg, ifun, ixyz)
                                 call fac_d(deg, ex, jfun, jxyz)
                                 call ijk_d(deg, jfun, jxyz)
                                 call get_integral_contribution(gint,      &
                                                                gint_temp, &
                                                                deg,       &
                                                                ang_temp,  &
                                                                ex,        &
                                                                coef,      &
                                                                xyz,       &
                                                                k,         &
                                                                recalc_integrals)
                                 call ijk_d(deg, ifun, ixyz)
                                 call ijk_u(deg, jfun, jxyz)
                              end if
                           end do
                        end do

                        ! du contribution
                        if (ifun /= jfun) then
                           ! if ifun == jfun then du raw integrals
                           ! are the same as ud, in this case do not recalculate
                           recalc_integrals = .true.
                        end if
                        ang_temp = ang
                        ang_temp(ifun) = ang_temp(ifun) - 1
                        ang_temp(jfun) = ang_temp(jfun) + 1
                        k = 0
                        do ixyz = 1, 3
                           icoor = (icent-1)*3 + ixyz
                           do jxyz = 1, 3
                              jcoor = (jcent-1)*3 + jxyz
                              if (jcoor >= icoor) then
                                 k = k + 1
                                 call init_prefactors(deg)
                                 call fac_d(deg, ex, ifun, ixyz)
                                 call ijk_d(deg, ifun, ixyz)
                                 call fac_u(deg, ex, jfun, jxyz)
                                 call ijk_u(deg, jfun, jxyz)
                                 call get_integral_contribution(gint,      &
                                                                gint_temp, &
                                                                deg,       &
                                                                ang_temp,  &
                                                                ex,        &
                                                                coef,      &
                                                                xyz,       &
                                                                k,         &
                                                                recalc_integrals)
                                 call ijk_u(deg, ifun, ixyz)
                                 call ijk_d(deg, jfun, jxyz)
                              end if
                           end do
                        end do

                        ! dd contribution
                        recalc_integrals = .true.
                        ang_temp = ang
                        ang_temp(ifun) = ang_temp(ifun) - 1
                        ang_temp(jfun) = ang_temp(jfun) - 1
                        k = 0
                        do ixyz = 1, 3
                           icoor = (icent-1)*3 + ixyz
                           do jxyz = 1, 3
                              jcoor = (jcent-1)*3 + jxyz
                              if (jcoor >= icoor) then
                                 k = k + 1
                                 call init_prefactors(deg)
                                 call fac_d(deg, ex, ifun, ixyz)
                                 call ijk_d(deg, ifun, ixyz)
                                 call fac_d(deg, ex, jfun, jxyz)
                                 call ijk_d(deg, jfun, jxyz)
                                 call get_integral_contribution(gint,      &
                                                                gint_temp, &
                                                                deg,       &
                                                                ang_temp,  &
                                                                ex,        &
                                                                coef,      &
                                                                xyz,       &
                                                                k,         &
                                                                recalc_integrals)
                                 call ijk_u(deg, ifun, ixyz)
                                 call ijk_u(deg, jfun, jxyz)
                              end if
                           end do
                        end do

                     end do
                  end do

                  ! contract integrals
                  k = 0
                  do ixyz = 1, 3
                     icoor = (icent-1)*3 + ixyz
                     do jxyz = 1, 3
                        jcoor = (jcent-1)*3 + jxyz
                        if (jcoor >= icoor) then
                           k = k + 1
                           l = (icoor - 1)*nr_centers*3 &
                             +  jcoor
                           call contract_integrals(deg,        &
                                                   off,        &
                                                   gint(1, k), &
                                                   ndim,       &
                                                   dmat,       &
                                                   gmat,       &
                                                   get_ave,    &
                                                   ave(l),     &
                                                   1.0d0)
                        end if
                     end do
                  end do

               end do
            end do

         case (3)
!-------------------------------------------------------------------------------
!           order 3
!-------------------------------------------------------------------------------

            nr_integrals = deg(1)*deg(2)*deg(3)*deg(4)
            call init_ijk(deg, ang)

            do icent = 1, nr_centers
               do jcent = icent, nr_centers
                  do kcent = jcent, nr_centers

                     ! zero out integrals
                     k = 0
                     do ixyz = 1, 3
                        icoor = (icent-1)*3 + ixyz
                        do jxyz = 1, 3
                           jcoor = (jcent-1)*3 + jxyz
                           do kxyz = 1, 3
                              kcoor = (kcent-1)*3 + kxyz
                              if (jcoor >= icoor) then
                                 if (kcoor >= jcoor) then
                                    k = k + 1
                                    do i = 1, nr_integrals
                                       gint(i, k) = 0.0d0
                                    end do
                                 end if
                              end if
                           end do
                        end do
                     end do

                     do ifun = 1, 4
                        if (cent(ifun) /= icent) cycle
                        do jfun = 1, 4
                           if (cent(jfun) /= jcent) cycle
                           do kfun = 1, 4
                              if (cent(kfun) /= kcent) cycle

                              ! uuu contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) + 1
                              ang_temp(jfun) = ang_temp(jfun) + 1
                              ang_temp(kfun) = ang_temp(kfun) + 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_u(deg, ex, ifun, ixyz)
                                             call ijk_u(deg, ifun, ixyz)
                                             call fac_u(deg, ex, jfun, jxyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call fac_u(deg, ex, kfun, kxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_d(deg, ifun, ixyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! duu contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) - 1
                              ang_temp(jfun) = ang_temp(jfun) + 1
                              ang_temp(kfun) = ang_temp(kfun) + 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_d(deg, ex, ifun, ixyz)
                                             call ijk_d(deg, ifun, ixyz)
                                             call fac_u(deg, ex, jfun, jxyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call fac_u(deg, ex, kfun, kxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_u(deg, ifun, ixyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! udu contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) + 1
                              ang_temp(jfun) = ang_temp(jfun) - 1
                              ang_temp(kfun) = ang_temp(kfun) + 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_u(deg, ex, ifun, ixyz)
                                             call ijk_u(deg, ifun, ixyz)
                                             call fac_d(deg, ex, jfun, jxyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call fac_u(deg, ex, kfun, kxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_d(deg, ifun, ixyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! uud contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) + 1
                              ang_temp(jfun) = ang_temp(jfun) + 1
                              ang_temp(kfun) = ang_temp(kfun) - 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_u(deg, ex, ifun, ixyz)
                                             call ijk_u(deg, ifun, ixyz)
                                             call fac_u(deg, ex, jfun, jxyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call fac_d(deg, ex, kfun, kxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_d(deg, ifun, ixyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! udd contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) + 1
                              ang_temp(jfun) = ang_temp(jfun) - 1
                              ang_temp(kfun) = ang_temp(kfun) - 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_u(deg, ex, ifun, ixyz)
                                             call ijk_u(deg, ifun, ixyz)
                                             call fac_d(deg, ex, jfun, jxyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call fac_d(deg, ex, kfun, kxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_d(deg, ifun, ixyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! dud contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) - 1
                              ang_temp(jfun) = ang_temp(jfun) + 1
                              ang_temp(kfun) = ang_temp(kfun) - 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_d(deg, ex, ifun, ixyz)
                                             call ijk_d(deg, ifun, ixyz)
                                             call fac_u(deg, ex, jfun, jxyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call fac_d(deg, ex, kfun, kxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_u(deg, ifun, ixyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! ddu contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) - 1
                              ang_temp(jfun) = ang_temp(jfun) - 1
                              ang_temp(kfun) = ang_temp(kfun) + 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_d(deg, ex, ifun, ixyz)
                                             call ijk_d(deg, ifun, ixyz)
                                             call fac_d(deg, ex, jfun, jxyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call fac_u(deg, ex, kfun, kxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_u(deg, ifun, ixyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                              ! ddd contribution
                              recalc_integrals = .true.
                              ang_temp = ang
                              ang_temp(ifun) = ang_temp(ifun) - 1
                              ang_temp(jfun) = ang_temp(jfun) - 1
                              ang_temp(kfun) = ang_temp(kfun) - 1
                              k = 0
                              do ixyz = 1, 3
                                 icoor = (icent-1)*3 + ixyz
                                 do jxyz = 1, 3
                                    jcoor = (jcent-1)*3 + jxyz
                                    do kxyz = 1, 3
                                       kcoor = (kcent-1)*3 + kxyz
                                       if (jcoor >= icoor) then
                                          if (kcoor >= jcoor) then
                                             k = k + 1
                                             call init_prefactors(deg)
                                             call fac_d(deg, ex, ifun, ixyz)
                                             call ijk_d(deg, ifun, ixyz)
                                             call fac_d(deg, ex, jfun, jxyz)
                                             call ijk_d(deg, jfun, jxyz)
                                             call fac_d(deg, ex, kfun, kxyz)
                                             call ijk_d(deg, kfun, kxyz)
                                             call get_integral_contribution(gint,      &
                                                                            gint_temp, &
                                                                            deg,       &
                                                                            ang_temp,  &
                                                                            ex,        &
                                                                            coef,      &
                                                                            xyz,       &
                                                                            k,         &
                                                                            recalc_integrals)
                                             call ijk_u(deg, ifun, ixyz)
                                             call ijk_u(deg, jfun, jxyz)
                                             call ijk_u(deg, kfun, kxyz)
                                          end if
                                       end if
                                    end do
                                 end do
                              end do

                           end do
                        end do
                     end do

                     ! contract integrals
                     k = 0
                     do ixyz = 1, 3
                        icoor = (icent-1)*3 + ixyz
                        do jxyz = 1, 3
                           jcoor = (jcent-1)*3 + jxyz
                           do kxyz = 1, 3
                              kcoor = (kcent-1)*3 + kxyz
                              if (jcoor >= icoor) then
                                 if (kcoor >= jcoor) then
                                    k = k + 1
                                    l = (icoor - 1)*nr_centers*3*nr_centers*3 &
                                      + (jcoor - 1)*nr_centers*3              &
                                      +  kcoor
                                    call contract_integrals(deg,        &
                                                            off,        &
                                                            gint(1, k), &
                                                            ndim,       &
                                                            dmat,       &
                                                            gmat,       &
                                                            get_ave,    &
                                                            ave(l),     &
                                                            1.0d0)
                                 end if
                              end if
                           end do
                        end do
                     end do

                  end do
               end do
            end do

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

      use module_interest_eri

      real(8), intent(inout) :: gint(*)
      integer, intent(in)    :: l(4)
      real(8), intent(in)    :: e(4)
      real(8), intent(in)    :: c(4)
      real(8), intent(in)    :: xyz(3, 4)

      real(8) :: f

      f = 1.0d0

      !> call InteRest library routine for a given batch
      call interest_eri(f,                                                   &
                        gint,                                                &
                        l(1)+1, e(1), xyz(1, 1), xyz(2, 1), xyz(3, 1), c(1), &
                        l(2)+1, e(2), xyz(1, 2), xyz(2, 2), xyz(3, 2), c(2), &
                        l(3)+1, e(3), xyz(1, 3), xyz(2, 3), xyz(3, 3), c(3), &
                        l(4)+1, e(4), xyz(1, 4), xyz(2, 4), xyz(3, 4), c(4))

   end subroutine

   subroutine contract_integrals(n,       &
                                 o,       &
                                 gint,    &
                                 ndim,    &
                                 dmat,    &
                                 gmat,    &
                                 get_ave, &
                                 average, &
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

   subroutine init_permanent_arrays()

      integer :: i, j, k, l, m, n
      integer :: ndim

      ndim = (maxl+1)*(maxl+2)/2

      if (allocated(ic_to_ijk)) deallocate(ic_to_ijk)
      if (allocated(ijk_to_ic)) deallocate(ijk_to_ic)
      if (allocated(ijk_list))       deallocate(ijk_list)
      if (allocated(prefactor_list))    deallocate(prefactor_list)
      if (allocated(cdeg))    deallocate(cdeg)

      allocate(ic_to_ijk(ndim, 3, 0:maxl))
      allocate(ijk_to_ic(0:maxl, 0:maxl, 0:maxl))
      allocate(ijk_list(ndim, 3, 4))
      allocate(prefactor_list(ndim, 4))
      allocate(cdeg(0:maxl))

      ic_to_ijk = 0
      ijk_to_ic = 0

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

   subroutine init_ijk(deg, ang)

      integer, intent(in) :: deg(4)
      integer, intent(in) :: ang(4)

      integer             :: ifun
      integer             :: ixyz
      integer             :: iang
      integer             :: ideg
      integer             :: ideg_max

      do ifun = 1, 4
         ideg_max = deg(ifun)
         iang = ang(ifun)
         do ixyz = 1, 3
            do ideg = 1, ideg_max
               ijk_list(ideg, ixyz, ifun) = ic_to_ijk(ideg, ixyz, iang)
            end do
         end do
      end do

   end subroutine

   subroutine init_prefactors(deg)

      integer, intent(in) :: deg(4)

      integer             :: ifun
      integer             :: ideg
      integer             :: ideg_max

      do ifun = 1, 4
         ideg_max = deg(ifun)
         do ideg = 1, ideg_max
            prefactor_list(ideg, ifun) = 1.0d0
         end do
      end do

   end subroutine

   subroutine fac_u(deg,  &
                    ex,   &
                    ifun, &
                    ixyz)

      integer, intent(in)    :: deg(4)
      real(8), intent(in)    :: ex(4)
      integer, intent(in)    :: ifun
      integer, intent(in)    :: ixyz

      integer                :: ideg
      integer                :: ideg_max
      real(8)                :: f

      f = 2.0d0*ex(ifun)

      ideg_max = deg(ifun)
      do ideg = 1, ideg_max
         prefactor_list(ideg, ifun) = f*prefactor_list(ideg, ifun)
      end do

   end subroutine

   subroutine ijk_u(deg, ifun, ixyz)

      integer, intent(in) :: deg(4)
      integer, intent(in) :: ifun
      integer, intent(in) :: ixyz

      integer             :: ideg
      integer             :: ideg_max

      ideg_max = deg(ifun)
      do ideg = 1, ideg_max
         ijk_list(ideg, ixyz, ifun) = ijk_list(ideg, ixyz, ifun) + 1
      end do

   end subroutine

   subroutine ijk_d(deg, ifun, ixyz)

      integer, intent(in) :: deg(4)
      integer, intent(in) :: ifun
      integer, intent(in) :: ixyz

      integer             :: ideg
      integer             :: ideg_max

      ideg_max = deg(ifun)
      do ideg = 1, ideg_max
         ijk_list(ideg, ixyz, ifun) = ijk_list(ideg, ixyz, ifun) - 1
      end do

   end subroutine

   subroutine fac_d(deg,  &
                    ex,   &
                    ifun, &
                    ixyz)

      integer, intent(in)    :: deg(4)
      real(8), intent(in)    :: ex(4)
      integer, intent(in)    :: ifun
      integer, intent(in)    :: ixyz

      integer                :: ideg
      integer                :: ideg_max

      ideg_max = deg(ifun)
      do ideg = 1, ideg_max
         if (ijk_list(ideg, ixyz, ifun) > 0) then
            prefactor_list(ideg, ifun) = -real(ijk_list(ideg, ixyz, ifun))*prefactor_list(ideg, ifun)
         else
            prefactor_list(ideg, ifun) = 0.0d0
         end if
      end do

   end subroutine

   logical function all_ang_nonnegative(ang)

      integer, intent(in) :: ang(4)

      integer             :: ifun

      all_ang_nonnegative = .false.

      do ifun = 1, 4
         if (ang(ifun) < 0) return
      end do

      all_ang_nonnegative = .true.

   end function

   subroutine get_integral_contribution(gint,      &
                                        gint_temp, &
                                        deg,       &
                                        ang_temp,  &
                                        ex,        &
                                        coef,      &
                                        xyz,       &
                                        icoor,     &
                                        recalc_integrals)

      real(8), intent(inout) :: gint(max_nr_integrals, max_ave_length)
      real(8), intent(inout) :: gint_temp(max_nr_integrals)
      integer, intent(in)    :: deg(4)
      integer, intent(in)    :: ang_temp(4)
      real(8), intent(in)    :: ex(4)
      real(8), intent(in)    :: coef(4)
      real(8), intent(in)    :: xyz(3, 4)
      integer, intent(in)    :: icoor
      logical, intent(inout) :: recalc_integrals

      if (all_ang_nonnegative(ang_temp)) then
         if (recalc_integrals) then
            call get_integrals(gint_temp, ang_temp, ex, coef, xyz)
            recalc_integrals = .false.
         end if
         call add_integrals(gint(1, icoor), &
                            gint_temp,      &
                            deg,            &
                            ang_temp,       &
                            ijk_list,       &
                            prefactor_list)
      end if

   end subroutine

end module
