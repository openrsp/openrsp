module interface_xc

   use, intrinsic :: iso_c_binding

   use matrix_defop, matrix => openrsp_matrix
   use interface_molecule, only: get_num_atoms
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
   use openrsp_cfg

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

   subroutine check_if_initialized()
      if (.not. is_initialized) then
         print *, 'ERROR: you try to access interface_xc'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine


   subroutine xcint_wakeup_workers()
#ifdef VAR_MPI

   ! local
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


   subroutine interface_xc_finalize()
      is_initialized = .false.
      this_is_ks_calculation = .false.
   end subroutine


   logical function is_ks_calculation()
      call check_if_initialized()
      is_ks_calculation = this_is_ks_calculation
   end function


   integer(c_int) function fortran_stdout_function(string) bind(c)

   ! input
      ! string to print (c_null_char-terminated)
      character(kind=c_char, len=1), intent(in) :: string(*)

   ! local
      integer(c_int) :: i
      integer(c_int) :: string_len
      character(80) :: line

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
      do i = 1, 80
         line(i:i) = ' '
      end do
      if (string_len < 81) then
         ! at the end we remove newline
         do i = 1, string_len-1
            line(i:i) = string(i)
         end do
      else
         do i = 1, 80
            line(i:i) = string(i)
         end do
      end if
      write(lupri, '(a80)') line
      fortran_stdout_function = 0

   end function


   subroutine interface_xc_init()

   ! local
      real(c_double), allocatable :: primitive_exp(:)
      real(c_double), allocatable :: contraction_coef(:)
      real(c_double), allocatable :: center_xyz(:)
      integer(c_int), allocatable :: shell_num_primitives(:)
      integer(c_int), allocatable :: l_quantum_num(:)
      integer(c_int), allocatable :: shell_center(:)
      integer(c_int), allocatable :: center_element(:)
      integer(c_int)              :: basis_type
      integer(c_int)              :: num_shells
      integer(c_int)              :: num_centers
      integer(c_int)              :: ierr
      integer                     :: i
      integer                     :: j
      integer                     :: icount
      integer                     :: iprim
      integer                     :: ishell
      integer                     :: iround
      integer                     :: n
      integer                     :: ixyz
      integer                     :: icenter

#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "shells.h"
#include "primit.h"
#ifdef VAR_MPI
#include "mpif.h"
#endif

      if (is_initialized) return

      num_shells  = kmax
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
         shell_center(ishell)  = ncent(ishell)
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

      if (openrsp_cfg_use_xcint_grid) then
         ! this generates XCint's internal grid
         ierr = xcint_generate_grid(openrsp_cfg_radint,   &
                                    openrsp_cfg_angmin,   &
                                    openrsp_cfg_angint,   &
                                    num_centers,          &
                                    center_xyz,           &
                                    center_element,       &
                                    num_shells,           &
                                    shell_center,         &
                                    l_quantum_num,        &
                                    shell_num_primitives, &
                                    primitive_exp)
         if (ierr /= 0) then
            print *, 'OpenRSP ERROR: problem in xcint_generate_grid'
            stop 1
         end if
      end if

      deallocate(primitive_exp)
      deallocate(contraction_coef)
      deallocate(center_xyz)
      deallocate(shell_num_primitives)
      deallocate(l_quantum_num)
      deallocate(shell_center)
      deallocate(center_element)

      is_initialized = .true.

   end subroutine


   subroutine openrsp_set_functional(line, hfx, mu, beta)

   ! input
      character(*), intent(in) :: line

   ! output
      real(c_double), intent(out) :: hfx
      real(c_double), intent(out) :: mu
      real(c_double), intent(out) :: beta

   ! local
      integer        :: ierr
      type(c_funptr) :: stdout_function

      stdout_function = c_funloc(fortran_stdout_function)
      call xcint_set_stdout_function(stdout_function)
      call xcint_set_stderr_function(stdout_function)

      call xcint_print_splash()

      ierr = xcint_set_functional(line//C_NULL_CHAR, hfx, mu, beta)
      if (ierr /= 0) then
         print *, 'OpenRSP ERROR: problem in xcint_set_functional'
         stop 1
      end if

      this_is_ks_calculation = .true.

   end subroutine


   recursive subroutine rsp_xcave_setup_dmat_perts(pert,           &
                                                   sofar,          &
                                                   kn,             &
                                                   rec_prog,       &
                                                   enc_len,        &
                                                   dmat_tuple_len, &
                                                   pert_ids,       &
                                                   enc_perts,      &
                                                   dmat_perts)
   ! input
      type(p_tuple), intent(in) :: pert
      integer,       intent(in) :: kn(2)
      integer,       intent(in) :: enc_len
      integer,       intent(in) :: dmat_tuple_len
      type(p_tuple), intent(in) :: enc_perts(enc_len)
      type(p_tuple), intent(in) :: dmat_perts(dmat_tuple_len)

   ! output
      integer, intent(inout) :: sofar
      integer, intent(inout) :: rec_prog
      integer, intent(inout) :: pert_ids(dmat_tuple_len)

   ! local
      integer       :: i
      integer       :: j
      logical       :: dmat_already
      type(p_tuple) :: psub(pert%n_perturbations)

      dmat_already = .false.

      ! unless at final recursion level, recurse further
      ! make all size (n - 1) subsets of the perturbations and recurse
      ! then (at final recursion level) get perturbed F, D, S
      if (pert%n_perturbations > 1) then

         call make_p_tuple_subset(pert, psub)

         do i = size(psub), 1, -1

            dmat_already = .false.
            do j = 1, sofar
               if (psub(i)%n_perturbations == dmat_perts(j)%n_perturbations) then
                  if (pid_compare(psub(i)%n_perturbations, psub(i)%pid, dmat_perts(j)%pid)) then
                     dmat_already = .true.
                  end if
               end if
            end do

            if (.not. dmat_already) then
               call rsp_xcave_setup_dmat_perts(psub(i),        &
                                               sofar,          &
                                               kn,             &
                                               rec_prog,       &
                                               enc_len,        &
                                               dmat_tuple_len, &
                                               pert_ids,       &
                                               enc_perts,      &
                                               dmat_perts)
            end if
         end do
      end if

      dmat_already = .false.
      do j = 1, rec_prog
         if (pert%n_perturbations == enc_perts(j)%n_perturbations) then
            if (pid_compare(pert%n_perturbations, pert%pid, enc_perts(j)%pid)) then
               dmat_already = .true.
            end if
         end if
      end do

      if (.not. dmat_already) then
         rec_prog = rec_prog + 1
         call p_tuple_p1_cloneto_p2(pert, enc_perts(rec_prog))

         if (.not. kn_skip(pert%n_perturbations, pert%pid, kn)) then
            sofar = sofar + 1
            call p_tuple_p1_cloneto_p2(pert, dmat_perts(sofar))
            pert_ids(sofar) = rec_prog
         end if
      end if

   end subroutine


   subroutine rsp_xcave_interface(pert, kn, D, prop_size, prop)
   ! assumes GEO perturbations come first to set up geo indices
   ! if that handling is moved to host program then no assumptions are necessary

   ! input
      type(p_tuple), intent(in) :: pert
      integer,       intent(in) :: kn(2)
      type(SDF),     intent(in) :: D
      integer,       intent(in) :: prop_size

   ! output
      complex(8), intent(inout) :: prop(prop_size)

   ! local
      real(c_double)              :: xc_energy(1)
      real(c_double)              :: num_electrons
      real(c_double)              :: dummy_real(1)
      integer(c_int)              :: num_pert
      integer(c_int), allocatable :: xc_pert(:)
      integer(c_int), allocatable :: comp(:)
      integer(c_int), allocatable :: dmat_to_pert(:)
      integer(c_int), allocatable :: dmat_to_comp(:)
      type(p_tuple),  allocatable :: dmat_perts(:)
      type(p_tuple),  allocatable :: enc_perts(:)
      type(matrix),   allocatable :: dmat_tuple(:)
      integer,        allocatable :: blk_sizes(:)
      integer,        allocatable :: one_ind(:)
      integer,        allocatable :: pert_ids(:)
      integer,        allocatable :: blk_info(:, :)
      integer,        allocatable :: indices(:, :)
      type(p_tuple)               :: emptypert
      integer                     :: element
      integer                     :: i
      integer                     :: j
      integer                     :: k
      integer                     :: num_blocks
      integer                     :: dmat_length
      integer                     :: ind_dmat_perts
      integer                     :: rec_prog
      integer                     :: enc_length

      call check_if_initialized()
      if (.not. is_ks_calculation()) return

      dmat_length = 2**(pert%n_perturbations - 1)
      enc_length  = 2**pert%n_perturbations

      allocate(dmat_perts(dmat_length))
      allocate(enc_perts(enc_length))
      allocate(pert_ids(dmat_length))
      allocate(dmat_tuple(dmat_length))

      emptypert = get_emptypert()
      call p_tuple_p1_cloneto_p2(emptypert, dmat_perts(1))
      call sdf_getdata_s(D, emptypert, (/1/), dmat_tuple(1))
      pert_ids(1) = 0
      ind_dmat_perts = 1
      rec_prog = 0

      call rsp_xcave_setup_dmat_perts(pert,           &
                                      ind_dmat_perts, &
                                      kn,             &
                                      rec_prog,       &
                                      enc_length,     &
                                      dmat_length,    &
                                      pert_ids,       &
                                      enc_perts,      &
                                      dmat_perts)

      num_blocks = get_num_blks(pert)
      allocate(blk_info(num_blocks, 3))
      allocate(blk_sizes(num_blocks))
      blk_info = get_blk_info(num_blocks, pert)
      blk_sizes = get_triangular_sizes(num_blocks,                &
                                       blk_info(1:num_blocks, 2), &
                                       blk_info(1:num_blocks, 3))

      allocate(indices(prop_size, pert%n_perturbations))
      allocate(one_ind(pert%n_perturbations))

      call make_triangulated_indices(num_blocks, blk_info, prop_size, indices)

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

         k = 1
         comp = 1
         do j = 1, count(pert%plab == 'GEO ')
            comp(k)   = indices(i, j)
            comp(k+1) = indices(i, j)
            k = k + 2
         end do

         ! first perturbation is always unperturbed D, handled above
         do j = 2, dmat_length
            do k = 1, dmat_perts(j)%n_perturbations
               one_ind(k) = indices(i, dmat_perts(j)%pid(k))
            end do
            call sdf_getdata_s(D,                                        &
                               dmat_perts(j),                            &
                               one_ind(1:dmat_perts(j)%n_perturbations), &
                               dmat_tuple(j))
         end do

         element = get_triang_blks_offset(num_blocks,           &
                                          pert%n_perturbations, &
                                          blk_info,             &
                                          blk_sizes,            &
                                          indices(i,:))

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
         prop(element) = prop(element) + cmplx(xc_energy(1), 0.0d0)

      end do

      deallocate(xc_pert)
      deallocate(comp)
      deallocate(dmat_to_pert)
      deallocate(dmat_to_comp)
      deallocate(dmat_perts)
      deallocate(enc_perts)
      deallocate(dmat_tuple)
      deallocate(blk_sizes)
      deallocate(one_ind)
      deallocate(pert_ids)
      deallocate(blk_info)
      deallocate(indices)

   end subroutine


   subroutine rsp_xcint_interface(pert_labels, D, F, Fg, Fgg)

   ! input
      character(4), intent(in) :: pert_labels(:)
      type(matrix), intent(in) :: D(:)

   ! output
      type(matrix), intent(inout), optional :: F
      type(matrix), intent(inout), optional :: Fg(:)
      type(matrix), intent(inout), optional :: Fgg(:, :)

   ! local
      integer                     :: num_pert
      integer                     :: i, j, n
      integer                     :: mat_dim
      integer                     :: num_atoms
      integer(c_int)              :: num_dmat
      real(c_double)              :: xc_energy(1)
      real(c_double)              :: num_electrons
      real(c_double), allocatable :: xc_mat(:)
      integer(c_int), allocatable :: pert(:)
      integer(c_int), allocatable :: comp(:)
      integer(c_int), allocatable :: dmat_to_pert(:)
      integer(c_int), allocatable :: dmat_to_comp(:)

      call check_if_initialized()
      if (.not. is_ks_calculation()) return

      num_atoms = get_num_atoms()
      num_dmat  = size(D)

      num_pert = 0
      if (present(F))   num_pert =     num_dmat - 1
      if (present(Fg))  num_pert = 1 + num_dmat - 1
      if (present(Fgg)) num_pert = 2 + num_dmat - 1

      allocate(pert(num_pert))
      allocate(comp(2*num_pert))
      comp = 0
      pert = XCINT_PERT_EL

      allocate(dmat_to_pert(num_dmat))
      allocate(dmat_to_comp(num_dmat))
      dmat_to_pert = 1 ! FIXME not used
      dmat_to_comp = 1 ! FIXME not used

      mat_dim = D(1)%nrow
      allocate(xc_mat(mat_dim*mat_dim))

      if (present(F)) then
         call xcint_wakeup_workers()
         call xcint_integrate(XCINT_MODE_RKS,                   &
                              num_pert,                         &
                              pert,                             &
                              comp,                             &
                              num_dmat,                         &
                              dmat_to_pert,                     &
                              dmat_to_comp,                     &
                              (/(D(n)%elms, n = 1, num_dmat)/), &
                              0,                                &
                              xc_energy,                        &
                              1,                                &
                              xc_mat,                           &
                              num_electrons)
         call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, F%elms, 1)
      end if

      if (present(Fg)) then
         pert(1) = XCINT_PERT_GEO
         do i = 1, num_atoms*3
            comp(1) = i
            comp(2) = i
            call xcint_wakeup_workers()
            call xcint_integrate(XCINT_MODE_RKS,                   &
                                 num_pert,                         &
                                 pert,                             &
                                 comp,                             &
                                 num_dmat,                         &
                                 dmat_to_pert,                     &
                                 dmat_to_comp,                     &
                                 (/(D(n)%elms, n = 1, num_dmat)/), &
                                 0,                                &
                                 xc_energy,                        &
                                 1,                                &
                                 xc_mat,                           &
                                 num_electrons)
            call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, Fg(i)%elms, 1)
         end do
      end if

      if (present(Fgg)) then
         pert(1) = XCINT_PERT_GEO
         pert(2) = XCINT_PERT_GEO
         do i = 1, num_atoms*3
            comp(1) = i
            comp(2) = i
            do j = 1, i
               comp(3) = j
               comp(4) = j
               call xcint_wakeup_workers()
               call xcint_integrate(XCINT_MODE_RKS,                   &
                                    num_pert,                         &
                                    pert,                             &
                                    comp,                             &
                                    num_dmat,                         &
                                    dmat_to_pert,                     &
                                    dmat_to_comp,                     &
                                    (/(D(n)%elms, n = 1, num_dmat)/), &
                                    0,                                &
                                    xc_energy,                        &
                                    1,                                &
                                    xc_mat,                           &
                                    num_electrons)
               call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, Fgg(i, j)%elms, 1)
            end do
         end do
      end if

      deallocate(pert)
      deallocate(comp)
      deallocate(dmat_to_pert)
      deallocate(dmat_to_comp)
      deallocate(xc_mat)

   end subroutine

end module
