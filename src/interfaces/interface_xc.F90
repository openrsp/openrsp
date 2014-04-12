module interface_xc

   use, intrinsic :: iso_c_binding

   use matrix_defop, matrix => openrsp_matrix
   use interface_molecule, only: get_nr_atoms
   use rsp_field_tuple, only: p_tuple,             &
                              pid_compare,         &
                              get_emptypert,       &
                              make_p_tuple_subset, &
                              p_tuple_p1_cloneto_p2
   use rsp_indices_and_addressing, only: kn_skip,                &
                                         get_num_blks,           &
                                         get_blk_info,           &
                                         get_triangular_sizes,   &
                                         get_triang_blks_offset, &
                                         make_triangulated_indices
   use rsp_sdf_caching, only: SDF, &
                              sdf_getdata_s
   use xcint_fortran_interface

   implicit none

   public interface_xc_init
   public interface_xc_finalize

   public xcint_wakeup_workers
   public rsp_xcave_interface
   public rsp_xcint_interface
   public is_ks_calculation
   public openrsp_set_functional

   private

   ! if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

   ! true if this is a KS calculation
   logical :: this_is_ks_calculation = .false.

contains

   subroutine xcint_wakeup_workers()
#ifdef VAR_MPI

   ! local variables
      integer :: num_proc
      integer :: iprint
      integer :: ierr

#include "mpif.h"
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"

      call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)
      if (num_proc > 1) then
         CALL MPIXBCAST(XCINT_MPI_WAKEUP_SIGNAL, 1,'INTEGER', MASTER)
         iprint = 0
         CALL MPIXBCAST(iprint, 1,'INTEGER', MASTER)
      end if
#endif
   end subroutine

   integer(c_int) function fortran_stdout_function(string) bind(c)

   ! input
      ! string to print (c_null_char-terminated)
      character(kind=c_char, len=1), intent(in) :: string(*)

   ! local variables
      integer(c_int) :: i
      integer(c_int) :: string_len

   ! returns
      ! 0 upon success

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
#ifdef VAR_MPI
#include "mpif.h"
#endif

      integer(c_int), allocatable :: shell_num_primitives(:)
      real(c_double), allocatable :: primitive_exp(:)
      real(c_double), allocatable :: contraction_coef(:)
      integer(c_int), allocatable :: l_quantum_num(:)
      integer(c_int), allocatable :: shell_center(:)
      integer                     :: i, j, icount, iprim, ishell, iround, n, ixyz, icenter
      integer(c_int)              :: basis_type
      integer(c_int)              :: num_shells
      integer(c_int)              :: num_centers
      real(c_double), allocatable :: center_xyz(:)
      integer(c_int), allocatable :: center_element(:)

      type(c_funptr) :: stdout_function

      if (is_initialized) return

      num_shells = kmax
      num_centers = nucind
      allocate(center_xyz(3*num_centers))
      n = 1
      do icenter = 1, num_centers
         do ixyz = 1, 3
            center_xyz(n) = cord(ixyz, icenter)
            n = n + 1
         end do
      end do

      allocate(center_element(num_centers))
      do i = 1, num_centers
         center_element(i) = nint(charge(i))
      end do

      allocate(shell_num_primitives(num_shells))
      n = 0
      do iround = 1, 2
         if (iround == 2) then
            allocate(primitive_exp(n))
            allocate(contraction_coef(n))
            n = 0
         end if
         do ishell = 1, num_shells
            i = jstrt(ishell) + 1
            j = jstrt(ishell) + nuco(ishell)
            icount = 0
            do iprim = i, j
               if (dabs(priccf(iprim, numcf(ishell))) > tiny(0.0d0)) then
                  icount = icount + 1
                  shell_num_primitives(ishell) = icount
                  n = n + 1
                  if (iround == 2) then
                     contraction_coef(n) = priccf(iprim, numcf(ishell))
                     primitive_exp(n)    = priexp(iprim)
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
      call xcint_set_mpi_comm(MPI_COMM_WORLD)
#endif

      stdout_function = c_funloc(fortran_stdout_function)
      call xcint_set_stdout_function(stdout_function)
      call xcint_set_stderr_function(stdout_function)

      call xcint_set_basis(basis_type,           &
                           num_centers,          &
                           center_xyz,           &
                           center_element,       &
                           num_shells,           &
                           shell_center,         &
                           l_quantum_num,        &
                           shell_num_primitives, &
                           primitive_exp,        &
                           contraction_coef)

      deallocate(center_xyz)
      deallocate(center_element)
      deallocate(shell_num_primitives)
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
      integer :: ierr
      ierr = xcint_set_functional(line//C_NULL_CHAR, hfx, mu, beta)
      this_is_ks_calculation = .true.
   end subroutine

   subroutine interface_xc_finalize()
      is_initialized = .false.
      this_is_ks_calculation = .false.
   end subroutine

   subroutine check_if_initialized()
      if (.not. is_initialized) then
         print *, 'ERROR: you try to access interface_xc'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   logical function is_ks_calculation()
      call check_if_initialized()
      is_ks_calculation = this_is_ks_calculation
   end function


  recursive subroutine rsp_xcave_setup_dmat_perts(pert, sofar, kn, rec_prog, &
                       enc_len, dmat_tuple_len, pert_ids, enc_perts, dmat_perts)

    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer :: i, j, sofar, dmat_tuple_len, rec_prog, enc_len
    integer, dimension(2) :: kn
    logical :: dmat_already
    integer, dimension(dmat_tuple_len) :: pert_ids
    type(p_tuple), dimension(dmat_tuple_len) :: dmat_perts
    type(p_tuple), dimension(enc_len) :: enc_perts

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

       end if

    end if

  end subroutine


  subroutine rsp_xcave_interface(pert, kn, D,  prop_size, prop)

    type(p_tuple) :: pert
    type(p_tuple), dimension(:), allocatable :: dmat_perts, enc_perts
    integer, dimension(2) :: kn
    type(SDF) :: D
    integer :: prop_size, element, ind_dmat_perts, rec_prog, enc_length
    complex(8), dimension(prop_size) :: prop, res
    type(p_tuple) :: emptypert
    type(matrix), dimension(:), allocatable :: dmat_tuple
    integer :: i, j, k, dmat_length, num_blks
    integer, dimension(:), allocatable :: blk_sizes, one_ind, pert_ids
    integer, dimension(:,:), allocatable :: blk_info, indices

    real(c_double)              :: xc_energy(1)
    real(c_double)              :: num_electrons
    real(c_double)              :: dummy_real(1)
    integer(c_int)              :: num_pert
    integer(c_int), allocatable :: xc_pert(:)
    integer(c_int), allocatable :: comp(:)
    integer(c_int), allocatable :: dmat_to_pert(:)
    integer(c_int), allocatable :: dmat_to_comp(:)

    if (.not. is_ks_calculation()) return

    res = 0.0

    dmat_length = 2**(pert%n_perturbations - 1)
    enc_length = 2**pert%n_perturbations

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

    num_pert = count(pert%plab == 'GEO ') + count(pert%plab == 'EL  ')
    allocate(xc_pert(num_pert))
    allocate(comp(2*num_pert))
    do j = 1, num_pert
       if (pert%plab(j) == 'GEO ') then
          xc_pert(j) = XCINT_PERT_GEO
       end if
       if (pert%plab(j) == 'EL  ') then
          xc_pert(j) = XCINT_PERT_EL
       end if
    end do

    allocate(dmat_to_pert(dmat_length))
    allocate(dmat_to_comp(dmat_length))
    do j = 1, dmat_length
       dmat_to_pert(j) = pert_ids(j)
       dmat_to_comp(j) = 1 ! FIXME not used
    end do

    do i = 1, prop_size

       comp = 1
       k = 1
       do j = 1, count(pert%plab == 'GEO ')
          comp(k)   = indices(i, j)
          comp(k+1) = indices(i, j)
          k = k + 2
       end do

       ! First perturbation is always unperturbed D, handled above
       do j = 2, dmat_length
          do k = 1, dmat_perts(j)%n_perturbations
             one_ind(k) = indices(i, dmat_perts(j)%pid(k))
          end do
          call sdf_getdata_s(D, dmat_perts(j), one_ind(1:dmat_perts(j)%n_perturbations), dmat_tuple(j))
       end do

       element = get_triang_blks_offset(num_blks, pert%n_perturbations, &
                                        blk_info, blk_sizes, indices(i,:))

       ! Assumes GEO perturbations come first to set up geo indices
       ! If that handling is moved to host program then no assumptions are necessary

       call xcint_wakeup_workers()
       call xcint_integrate(XCINT_MODE_RKS,                               &
                            num_pert,                                     &
                            xc_pert,                                      &
                            comp,                                         &
                            dmat_length,                                  &
                            dmat_to_pert,                                 &
                            dmat_to_comp,                                 &
                            (/(dmat_tuple(j)%elms, j = 1, dmat_length)/), &
                            1,                                            &
                            xc_energy,                                    &
                            0,                                            &
                            dummy_real,                                   &
                            num_electrons)
       res(element) = cmplx(xc_energy(1), 0.0d0)

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
    deallocate(xc_pert)
    deallocate(comp)
    deallocate(dmat_to_pert)
    deallocate(dmat_to_comp)

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
      integer(c_int)              :: num_dmat
      integer(c_int), allocatable :: pert(:)
      integer(c_int), allocatable :: comp(:)
      integer(c_int), allocatable :: dmat_to_pert(:)
      integer(c_int), allocatable :: dmat_to_comp(:)
      real(c_double), allocatable :: xc_dmat(:)
      real(c_double), allocatable :: xc_mat(:)
      real(c_double)              :: xc_energy(1)
      real(c_double)              :: num_electrons
!     ---------------------------------------------------------------------------

      if (.not. is_ks_calculation()) then
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

      allocate(xc_mat(mat_dim*mat_dim))
      allocate(dmat_to_pert(num_dmat))
      allocate(dmat_to_comp(num_dmat))
      dmat_to_pert = 1 ! FIXME not used
      dmat_to_comp = 1

      if (present(F)) then
         allocate(pert(num_dmat-1))
         allocate(comp(2*(num_dmat-1)))
         pert = XCINT_PERT_EL
         comp = 0
         call xcint_wakeup_workers()
         call xcint_integrate(XCINT_MODE_RKS, &
                              num_dmat-1,     &
                              pert,           &
                              comp,           &
                              num_dmat,       &
                              dmat_to_pert,   &
                              dmat_to_comp,   &
                              xc_dmat,        &
                              0,              &
                              xc_energy,      &
                              1,              &
                              xc_mat,         &
                              num_electrons)
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, F%elms, 1)
      end if

      if (present(Fg)) then
         allocate(pert(1+num_dmat-1))
         allocate(comp(2*(1+num_dmat-1)))
         pert = XCINT_PERT_EL
         pert(1) = XCINT_PERT_GEO
         comp = 0
         do i = 1, num_atoms*3
            comp(1) = i
            comp(2) = i
            call xcint_wakeup_workers()
            call xcint_integrate(XCINT_MODE_RKS, &
                                 1+num_dmat-1,   &
                                 pert,           &
                                 comp,           &
                                 num_dmat,       &
                              dmat_to_pert,   &
                              dmat_to_comp,   &
                                 xc_dmat,        &
                                 0,              &
                                 xc_energy,      &
                                 1,              &
                                 xc_mat,         &
                                 num_electrons)
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         allocate(pert(2+num_dmat-1))
         allocate(comp(2*(2+num_dmat-1)))
         pert = XCINT_PERT_EL
         pert(1) = XCINT_PERT_GEO
         pert(2) = XCINT_PERT_GEO
         comp = 0
         do i = 1, num_atoms*3
            comp(1) = i
            comp(2) = i
            do j = 1, i
               comp(3) = j
               comp(4) = j
               call xcint_wakeup_workers()
               call xcint_integrate(XCINT_MODE_RKS, &
                                    2+num_dmat-1,   &
                                    pert,           &
                                    comp,           &
                                    num_dmat,       &
                              dmat_to_pert,   &
                              dmat_to_comp,   &
                                    xc_dmat,        &
                                    0,              &
                                    xc_energy,      &
                                    1,              &
                                    xc_mat,         &
                                    num_electrons)
               call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, Fgg(i, j)%elms, 1)
            end do
         end do
      end if

      deallocate(xc_dmat)
      deallocate(xc_mat)
      deallocate(dmat_to_pert)
      deallocate(dmat_to_comp)
      if (allocated(pert)) deallocate (pert)
      if (allocated(comp)) deallocate (comp)

   end subroutine

end module
