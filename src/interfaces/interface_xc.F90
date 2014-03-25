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
   use xcint_fortran_interface

   implicit none

   public interface_xc_init
   public interface_xc_finalize
   public xcint_wakeup_workers

   public di_get_MagDeriv_FxD_DFT
   public di_get_MagDeriv_GxD_DFT

   public rsp_xcave_interface
   public rsp_xcave_interface_new
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

   integer :: ierr

contains

   subroutine xcint_wakeup_workers()
      integer :: num_proc
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

   integer(c_int) function fortran_stdout_function(string) bind(c)
      character(kind=c_char, len=1), intent(in) :: string(*)
      integer(c_int) :: i, string_len
#include "priunit.h"
      i = 1
      do while (.true.)
         if (string(i) == c_null_char) then
            string_len = i - 1
            exit
         end if
         i = i + 1
      end do
      write(lupri, *) string(1:string_len)
      fortran_stdout_function = 0
   end function

   subroutine interface_xc_init()

#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "shells.h"
#include "primit.h"

      integer(c_int), allocatable :: num_primitives_per_shell(:)
      real(c_double), allocatable :: primitive_exp(:, :)
      real(c_double), allocatable :: contraction_coef(:, :)
      integer(c_int), allocatable :: l_quantum_num(:)
      integer(c_int), allocatable :: shell_center(:)
      integer                     :: i, j, icount, iprim, ishell, iround
      integer(c_int)              :: basis_type
      integer(c_int)              :: num_shells
      integer(c_int)              :: num_centers
      real(c_double), allocatable :: center_xyz(:, :)
      integer(c_int), allocatable :: center_element(:)

      type(c_funptr) :: stdout_function

      if (is_initialized) return

      num_shells = kmax
      num_centers = nucind
      allocate(center_xyz(3, num_centers))
      center_xyz = cord

      allocate(center_element(num_centers))
      do i = 1, num_centers
         center_element(i) = nint(charge(i))
      end do

      allocate(num_primitives_per_shell(num_shells))
      do iround = 1, 2
         if (iround == 2) then
            allocate(primitive_exp(maxval(num_primitives_per_shell), num_shells))
            allocate(contraction_coef(maxval(num_primitives_per_shell), num_shells))
         end if
         do ishell = 1, num_shells
            i = jstrt(ishell) + 1
            j = jstrt(ishell) + nuco(ishell)
            icount = 0
            do iprim = i, j
               if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
                  icount = icount + 1
                  num_primitives_per_shell(ishell) = icount
                  if (iround == 2) then
                     contraction_coef(icount, ishell) = priccf(iprim, numcf(ishell))
                     primitive_exp(icount, ishell)    = priexp(iprim)
                  end if
               end if
            end do
         end do
      end do

      allocate(l_quantum_num(num_shells))
      allocate(shell_center(num_shells))
      do ishell = 1, num_shells
         l_quantum_num(ishell) = nhkt(ishell) - 1
         shell_center(ishell) = ncent(ishell)
      end do

      ! we default to spherical basis
      basis_type = XCINT_BASIS_SPHERICAL
      do ishell = 1, num_shells
         if ((l_quantum_num(ishell) > 1) .and. .not. sphr(ishell)) then
            ! basis is cartesian
            basis_type = XCINT_BASIS_CARTESIAN
            exit
         end if
      end do

#ifdef VAR_MPI
      ierr = xcint_set_mpi_comm(MPI_COMM_WORLD)
#endif

      stdout_function = c_funloc(fortran_stdout_function)
      call xcint_set_stdout_function(stdout_function)
      call xcint_set_stderr_function(stdout_function)

      call xcint_set_basis(basis_type,                                      &
                           num_centers,                                     &
                           reshape(center_xyz, (/size(center_xyz)/)),       &
                           center_element,                                  &
                           num_shells,                                      &
                           shell_center,                                    &
                           l_quantum_num,                                   &
                           num_primitives_per_shell,                        &
                           reshape(primitive_exp, (/size(primitive_exp)/)), &
                           reshape(contraction_coef, (/size(contraction_coef)/)))

      deallocate(center_xyz)
      deallocate(center_element)
      deallocate(num_primitives_per_shell)
      deallocate(primitive_exp)
      deallocate(contraction_coef)
      deallocate(l_quantum_num)
      deallocate(shell_center)

      is_initialized = .true.

   end subroutine

   subroutine openrsp_set_functional(line, hfx, mu, beta)
      character(*),   intent(in)  :: line
      real(c_double), intent(out) :: hfx
      real(c_double), intent(out) :: mu
      real(c_double), intent(out) :: beta
      ierr = xcint_set_functional(line//C_NULL_CHAR, hfx, mu, beta)
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





  recursive subroutine rsp_xcave_setup_dmat_perts(pert, sofar, kn, rec_prog, &
                       enc_len, dmat_tuple_len, pert_ids, enc_perts, dmat_perts)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer :: i, j, sofar, dmat_tuple_len, rec_prog, enc_len
    integer, dimension(2) :: kn
    logical :: dmat_already
    integer, dimension(dmat_tuple_len) :: pert_ids
    type(p_tuple), dimension(dmat_tuple_len) :: dmat_perts
    type(p_tuple), dimension(enc_len) :: enc_perts

!         write(*,*) 'called:', pert%pid

    dmat_already = .FALSE.

    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then (at final recursion level) get perturbed F, D, S
    if (pert%n_perturbations > 1) then

       call make_p_tuple_subset(pert, psub)

       do i = size(psub), 1, -1

          dmat_already = .FALSE.

          do j = 1, sofar


             if (psub(i)%n_perturbations == dmat_perts(j)%n_perturbations) then

                if (pid_compare(psub(i)%n_perturbations, psub(i)%pid, dmat_perts(j)%pid)) then

                   dmat_already = .TRUE.

                end if

             end if


          end do

          if (.NOT.(dmat_already)) then

             call rsp_xcave_setup_dmat_perts(psub(i), sofar, kn, rec_prog, enc_len, &
                  dmat_tuple_len, pert_ids, enc_perts, dmat_perts)

          end if

       end do

    end if

    dmat_already = .FALSE.

    do j = 1, rec_prog

       if (pert%n_perturbations == enc_perts(j)%n_perturbations) then

          if (pid_compare(pert%n_perturbations, pert%pid, enc_perts(j)%pid)) then

             dmat_already = .TRUE.

          end if

       end if


    end do

    if (.NOT.(dmat_already)) then

       rec_prog = rec_prog + 1
       call p_tuple_p1_cloneto_p2(pert, enc_perts(rec_prog))

       if (.NOT.(kn_skip(pert%n_perturbations, pert%pid, kn))) then

             sofar = sofar + 1
             call p_tuple_p1_cloneto_p2(pert, dmat_perts(sofar))
             pert_ids(sofar) = rec_prog
!              write(*,*) 'pert id', dmat_perts(sofar)%pid

       end if

    end if

  end subroutine


  subroutine rsp_xcave_interface_new(pert, kn, D,  prop_size, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(:), allocatable :: dmat_perts, enc_perts
    integer, dimension(2) :: kn
    type(SDF) :: D
! MaR: Uncomment this, add to argument list and handle below when callback functionality is ready
!     external :: rsp_get_xcave
    integer :: prop_size, element, ind_dmat_perts, rec_prog, enc_length
    complex(8), dimension(prop_size) :: prop, res
    type(p_tuple) :: emptypert
    type(matrix), dimension(:), allocatable :: dmat_tuple
    integer :: i, j, k, dmat_length, num_blks
    integer, dimension(:), allocatable :: blk_sizes, one_ind, pert_ids
    integer, dimension(:,:), allocatable :: blk_info, indices
    real(c_double)               :: xc_energy
    integer(c_int)              :: get_ave
    integer(c_int)              :: force_sequential

    if (.not. get_is_ks_calculation()) return

    res = 0.0

    dmat_length = 2**(pert%n_perturbations - 1)
    enc_length = 2**pert%n_perturbations

!         dmat_length = pert%n_perturbations * (pert%n_perturbations + 1) / 2

    emptypert = get_emptypert()

    allocate(dmat_perts(dmat_length))
    allocate(enc_perts(enc_length))
    allocate(pert_ids(dmat_length))
    allocate(dmat_tuple(dmat_length))

    call p_tuple_p1_cloneto_p2(emptypert, dmat_perts(1))
    call sdf_getdata_s(D, emptypert, (/1/), dmat_tuple(1))
    pert_ids(1) = 0

    ind_dmat_perts = 1
    rec_prog = 0

    call rsp_xcave_setup_dmat_perts(pert, ind_dmat_perts, kn, rec_prog, &
         enc_length, dmat_length, pert_ids, enc_perts, dmat_perts)



    num_blks = get_num_blks(pert)
    allocate(blk_info(num_blks, 3))
    allocate(blk_sizes(num_blks))
    blk_info = get_blk_info(num_blks, pert)
    blk_sizes = get_triangular_sizes(num_blks, blk_info(1:num_blks, 2), &
                                     blk_info(1:num_blks, 3))

    allocate(indices(prop_size, pert%n_perturbations))
    allocate(one_ind(pert%n_perturbations))


    call make_triangulated_indices(num_blks, blk_info, prop_size, indices)

    write(*,*) 'These p tuples make up the dmat_tuple array', pert_ids

    write(*,*) 'The perturbation ids of these tuples are as follows (blank line means unperturbed):'

    do j = 1, dmat_length

       write(*,*) dmat_perts(j)%pid

    end do

    do i = 1, prop_size

       ! First perturbation is always unperturbed D, handled above
       do j = 2, dmat_length

          do k = 1, dmat_perts(j)%n_perturbations

             one_ind(k) = indices(i, dmat_perts(j)%pid(k))

          end do

          call sdf_getdata_s(D, dmat_perts(j), one_ind(1:dmat_perts(j)%n_perturbations), dmat_tuple(j))

       end do

       element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                        blk_info, blk_sizes, indices(i,:))

       call xcint_wakeup_workers()

       ! Assumes GEO perturbations come first to set up geo indices
       ! If that handling is moved to host program then no assumptions are necessary
       ! This routine is not sufficiently general and should be upgraded
       ! Choice of kn is passed but should not be relevant anymore
       ! Unsure if dmat tuple is correct w.r.t. what XCInt expects

! Commented out until updated

!        call xcint_integrate(XCINT_MODE_RKS, dmat_length,                 &
!                             (/(dmat_tuple(j)%elms, j = 1, dmat_length)/), &
!                             (/0.0d0/),                   &
!                             xc_energy,                   &
!                             get_ave,                     &
!                             count(pert%plab == 'GEO '),  &
!                             (/ (indices(i, j), j = 1, count(pert%plab == 'GEO ')) /), &
!                             count(pert%plab == 'EL  '),  &
!                             kn,                          &
!                             force_sequential)
!
!        res(element) = cmplx(xc_energy, 0.0d0)

    end do

    prop = prop + res

    deallocate(one_ind)
    deallocate(indices)
    deallocate(blk_info)
    deallocate(blk_sizes)
    deallocate(dmat_perts)
    deallocate(enc_perts)
    deallocate(pert_ids)
    deallocate(dmat_tuple)

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
      integer                      :: num_dmat
      integer                      :: num_atoms
      integer                      :: num_geo
      integer                      :: num_el
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

      num_atoms = get_nr_atoms()

      num_geo = count(pert%plab == 'GEO ')
      num_el  = count(pert%plab == 'EL  ')


      combination_found = .false.

      num_dmat = 2**(num_geo + num_el - 1)

      allocate(dmat_tuple(num_dmat))

      allocate(emptyint(0))

      ! MaR: Begin EL only cases

      if (num_geo == 0 .and. num_el == 1) then

      ! No contribution

      end if

      if (num_geo == 0 .and. num_el == 2) then
         combination_found = .true.

      ! No contribution

      end if

      if (num_geo == 0 .and. num_el == 3) then
         combination_found = .true.

      ! Assumes (k,n) = (0,2)
      ! Then no contribution

      end if

      if (num_geo == 0 .and. num_el == 4) then
         combination_found = .true.

      ! Assumes (k,n) = (0,3)
      ! Then no contribution

      end if

      if (num_geo == 0 .and. num_el == 5) then
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
                        call xcint_integrate(XCINT_MODE_RKS, num_dmat,                      &
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


      if (num_geo == 0 .and. num_el == 6) then
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
                        call xcint_integrate(XCINT_MODE_RKS, num_dmat,                      &
                                             (/(dmat_tuple(i)%elms, i = 1, num_dmat)/), &
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

      if (num_geo == 1 .and. num_el == 0) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))
         do i = 1, num_atoms*3
            element = i
            call xcint_wakeup_workers()
            call xcint_integrate(XCINT_MODE_RKS, 1,                      &
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

      if (num_geo == 2 .and. num_el == 0) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, num_atoms*3
            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xcint_wakeup_workers()
               call xcint_integrate(XCINT_MODE_RKS, 2,                      &
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

      if (num_geo == 1 .and. num_el == 1) then
         combination_found = .true.

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do j = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(2))

            do i = 1, num_atoms*3

               element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                         blk_info, blk_sizes, (/i, j/))

               call xcint_wakeup_workers()
               call xcint_integrate(XCINT_MODE_RKS, 2,                      &
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

      if (num_geo == 3 .and. num_el == 0) then
         combination_found = .true.


         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, num_atoms*3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))

            do j = 1, i

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do k = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(XCINT_MODE_RKS, 4,                      &
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

      if (num_geo == 2 .and. num_el == 1) then
         combination_found = .true.

         ! Set up for (k,n) = (1,1)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do k = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 3), (/k/), dmat_tuple(4))

            do j = 1, num_atoms*3

               call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 2), (/j/), dmat_tuple(3))

               do i = 1, j

                  call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 1), (/i/), dmat_tuple(2))


                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(XCINT_MODE_RKS, 4,                      &
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

      if (num_geo == 1 .and. num_el == 2) then
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

               do i = 1, num_atoms*3
                  ipart = ipart + 1
                  print *, 'computing XC ave GFF contribution', ipart, 'out of', 3*maxcomp2*num_atoms*3

                  element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                            blk_info, blk_sizes, (/i, j, k/))

                  call xcint_wakeup_workers()
                  call xcint_integrate(XCINT_MODE_RKS, 4,                      &
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

      if (num_geo == 4 .and. num_el == 0) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do l = 1, num_atoms*3

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
                     call xcint_integrate(XCINT_MODE_RKS, 8,                      &
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

      if (num_geo == 3 .and. num_el == 1) then
         combination_found = .true.

         ! Using (k,n) = (1,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do l = 1, 3

            call sdf_getdata_s(D_sdf, p_tuple_getone(pert, 4), (/l/), dmat_tuple(6))

            do k = 1, num_atoms*3

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
                     call xcint_integrate(XCINT_MODE_RKS, 8,                      &
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

      if (num_geo == 2 .and. num_el == 2) then
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

               do j = 1, num_atoms*3

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
                     call xcint_integrate(XCINT_MODE_RKS, 8,                      &
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

      if (num_geo == 1 .and. num_el == 3) then
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

                  do l = 1, num_atoms*3

                     ! Get the offset in the response tensor
                     element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                               blk_info, blk_sizes, (/l, i, j, k/))

                     call xcint_wakeup_workers()
                     call xcint_integrate(XCINT_MODE_RKS, 8,                      &
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

      if (num_geo == 1 .and. num_el == 4) then
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

                     do m = 1, num_atoms*3

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

      if (num_geo == 2 .and. num_el == 3) then
          combination_found = .true.

          ! Using (k,n) = (1,3)

          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

          do i = 1, num_atoms*3

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

      if (num_geo == 3 .and. num_el == 2) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, num_atoms*3

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

      if (num_geo == 4 .and. num_el == 1) then
         combination_found = .true.

         ! Using (k,n) = (2,2)

         call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), dmat_tuple(1))

         do i = 1, num_atoms*3

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
                  num_geo, &
                  num_el
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
      integer              :: num_atoms
      integer(c_int)       :: num_dmat
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

      num_atoms = get_nr_atoms()
      num_dmat  = size(D)

      mat_dim = D(1)%nrow
      allocate(xc_dmat(mat_dim*mat_dim*num_dmat))
      xc_dmat = 0.0d0
      do imat = 1, num_dmat
         call daxpy(mat_dim*mat_dim, 1.0d0, D(imat)%elms, 1, xc_dmat((imat-1)*mat_dim*mat_dim + 1), 1)
      end do

      allocate(xc_fmat(mat_dim*mat_dim))

      if (present(F)) then
         call xcint_wakeup_workers()
         call xcint_integrate(XCINT_MODE_RKS, num_dmat,   &
                              xc_dmat,   &
                              xc_fmat,   &
                              xc_energy, &
                              get_ave,   &
                              0,         &
                              (/0/),     &
                              num_dmat-1, &
                              (/0, 0/),  &
                              force_sequential)
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms, 1)
      end if

      if (present(Fg)) then
         do i = 1, num_atoms*3
            call xcint_wakeup_workers()
            call xcint_integrate(XCINT_MODE_RKS, num_dmat,   &
                                 xc_dmat,   &
                                 xc_fmat,   &
                                 xc_energy, &
                                 get_ave,   &
                                 1,         &
                                 (/i/),     &
                                 num_dmat-1, &
                                 (/0, 0/),  &
                                 force_sequential)
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         do i = 1, num_atoms*3
            do j = 1, i
               call xcint_wakeup_workers()
               call xcint_integrate(XCINT_MODE_RKS, num_dmat,   &
                                    xc_dmat,   &
                                    xc_fmat,   &
                                    xc_energy, &
                                    get_ave,   &
                                    2,         &
                                    (/i, j/),  &
                                    num_dmat-1, &
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
