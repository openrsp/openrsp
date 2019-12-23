! Copyright 2015 Magnus Ringholm
!
!> @file Contains module rsp_general

!> General response routines. This module organizes, computes and prints
!> response function tensors.
module rsp_general

  use rsp_field_tuple, only: p_tuple,                &
                             p_tuple_standardorder,  &
                             p_tuple_remove_first,   &
                             p_tuple_getone,         &
                             p_tuple_extend,         &
                             get_emptypert,          &
                             empty_p_tuple,          &
                             p_tuples_standardorder, &
                             merge_p_tuple,          &
                             p1_cloneto_p2
  use rsp_indices_and_addressing, only: mem_manager,                  &
                                        mem_set_status,               &
                                        mem_incr,                     &
                                        mem_decr,                     &
                                        mem_enough,                   &
                                        mem_exceed,                   &  
                                        get_blk_info,                 &
                                        get_triangular_sizes,         &
                                        get_triangulated_size,        &
                                        get_num_blks,                 &
                                        kn_skip,                      &
                                        get_triang_blks_offset,       &
                                        get_triang_blks_tuple_offset, &
                                        nc_only,                      &
                                        nc_onlysmall,                 &
                                        get_ncarray,                  &
                                        make_triangulated_indices,    &
                                        make_triangulated_tuples_indices
  use rsp_perturbed_matrices, only: derivative_superstructure_getsize, &
                                    derivative_superstructure,         &
                                    rsp_get_matrix_zeta,          &
                                    rsp_get_matrix_lambda,        &
                                    rsp_get_matrix_z,             &
                                    rsp_get_matrix_w,             &
                                    rsp_get_matrix_y
  use rsp_perturbed_sdf, only: rsp_fds, &
                               rsp_xc_wrapper
  use rsp_property_caching, only: contrib_cache_outer,                 &
                                  contrib_cache,                       &
                                  contrib_cache_initialize,            &
                                  contrib_cache_outer_allocate,        &
                                  contrib_cache_outer_add_element,     &
                                  contrib_cache_add_element,           &
                                  contrib_cache_already,               &
                                  contrib_cache_getdata,               &
                                  contrib_cache_getdata_outer,               &
                                  contrib_cache_allocate, &
                                  contrib_cache_retrieve, &
                                  contrib_cache_outer_retrieve, &
                                  contrib_cache_outer_store, &
                                  contrib_cache_store, &
                                  contrib_cache_locate, &
                                  mat_scal_store, &
                                  mat_scal_retrieve, &
                                  rs_check, &
                                  prog_incr, &
                                  prog_init
                                  
  
  use rsp_sdf_caching
  
  use qcmatrix_f

  implicit none

  public openrsp_get_property
!   public openrsp_get_residue
  public print_rsp_tensor
  public print_rsp_tensor_stdout
  public print_rsp_tensor_stdout_tr
  public rsp_xc_wrapper

  private

  real(8) :: time_start
  real(8) :: time_end

  contains
  
  ! Main routine called by host program to calculate properties
  !
  ! Arguments:
  !
  ! n_props: Number of properties (one property may have several frequency configurations)
  ! np: For each property: Number of perturbations (i.e. the order of each property)
  ! pert_dims: For each perturbation across all properties: Dimensionality of the perturbation
  ! pert_first_comp: For each perturbation across all properties: Which is the first component (not in use)
  ! pert_labels: For each perturbation across all properties: Perturbation label
  ! n_freq_cfgs: For each property: Number of frequency configurations for that property
  ! pert_freqs: For each perturbation across all properties and freq. cfgs.: Perturbation frequency
  ! kn_rules: For each property: Choice of (k,n) truncation rule (represented as k for each property)
  ! F_unpert, S_unpert, D_unpert: The unperturbed Fock, overlap and density matrices, respectively
  ! get_rsp_sol: Callback routine for response equation solution
  ! get_ovl_mat: Callback routine for perturbed overlap matrices
  ! get_ovl_exp: Callback routine for Pulay S*W type terms
  ! get_1el_mat: Callback routine for 1-electron perturbed Fock matrix contributions
  ! get_1el_exp: Callback routine for 1-electron contributions to response properties
  ! get_2el_mat: Callback routine for 2-electron perturbed Fock matrix contributions
  ! get_2el_exp: Callback routine for 2-electron contributions to response properties
  ! get_xc_mat: Callback routine for exchange-correlation contributions to perturbed Fock matrices
  ! get_xc_exp: Callback routine for exchange-correlation contributions to response properties
  ! rsp_tensor: Array to hold the response properties upon calculation
  ! file_id: Custom filename for printing response properties (currently not in use)
  ! mem_calibrate_arg: Flag: OpenRSP is to be run in memory calibration mode (must then give next args.)
  ! max_mat_mem: Max # of matrices to be created by OpenRSP
  ! mem_result: Optional (for mem. calibration mode): Calibration result
  

  subroutine openrsp_get_property(n_props, np, pert_dims, pert_first_comp, &
                                 pert_labels, n_freq_cfgs, pert_freqs, &
                                 kn_rules, F_unpert, S_unpert, D_unpert, &
                                 get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, &
                                 get_1el_mat, get_1el_exp, get_2el_mat, get_2el_exp, &
                                 get_xc_mat, get_xc_exp, out_print, r_flag_in, write_threshold, &
                                 rsp_tensor_size, rsp_tensor, residue_order, file_id, mem_calibrate, &
                                 max_mat, mem_result, residue_spec_pert, size_rsi_1, &
                                 residue_spec_index, exenerg, Xf_unpert)
    implicit none

    
    integer(kind=QINT), intent(in) :: n_props
    
    integer(kind=QINT), dimension(n_props), intent(in) :: np, n_freq_cfgs
    integer(kind=QINT), dimension(sum(np)), intent(in) :: pert_dims, pert_first_comp
    character(4), dimension(sum(np)), intent(in) :: pert_labels
    
    character(256) :: filename
    integer :: i, j, k, l, m, n
    integer :: dum_ind
    integer :: len_fds, len_x
    integer, dimension(sum(n_freq_cfgs)) :: prop_sizes, num_blks
    integer(kind=QINT), intent(in), dimension(n_props) :: kn_rules
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    logical :: lfreq_match
    
    integer, allocatable, dimension(:) :: blk_sizes
    integer, allocatable, dimension(:,:) :: blk_info
    complex(8), dimension(dot_product(np, n_freq_cfgs)), intent(in) :: pert_freqs
    
    integer(kind=QINT) num_perts
    real :: timing_start, timing_end
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), dimension(1) :: empty_pert

    external :: get_rsp_sol, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, get_nucpot
    external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
    external :: out_print
    character(len=2047) :: out_str
    integer(kind=QINT), intent(in) :: rsp_tensor_size
    complex(8), dimension(rsp_tensor_size) :: rsp_tensor
    type(QcMat) :: S_unpert, D_unpert, F_unpert ! NOTE: Make optional to exclude in mem. calibration mode

    type(QcMat), dimension(1) :: S_unpert_arr, D_unpert_arr, F_unpert_arr
    type(QcMat), dimension(1) :: Xf_unpert_arr

    type(contrib_cache_outer), allocatable, dimension(:) :: S, D, F, Xf
    integer :: kn(2)
    character(30) :: fmt_str
    real, parameter :: xtiny=1.0d-8
    
    integer(kind=QINT), intent(in) :: residue_order
    
    character, optional, dimension(20) :: file_id
    
    logical, optional :: mem_calibrate
    integer, optional :: max_mat, mem_result
    logical :: xf_was_allocated
    type(mem_manager) :: mem_mgr
    
    integer(kind=QINT), intent(in), optional :: residue_spec_pert(residue_order)
    integer(kind=QINT), intent(in), optional :: size_rsi_1
    integer(kind=QINT), dimension(*), optional :: residue_spec_index
    complex(8), dimension(*), optional :: exenerg
    type(QcMat), optional, dimension(*) :: Xf_unpert 
    
    integer, allocatable, dimension(:,:) :: indices
    real(kind=QREAL) :: write_threshold
    integer :: p
    integer(kind=QINT) :: r_flag_in
    
    ! Restarting data
    
    ! To be connected to API: Flag to determine the restarting setup
    ! Meaning:
    ! 0: Do not load/use any existing restarting data and do not save any new restarting data
    ! (UNUSED) 1: Load and use all existing restarting data but do not save any new restarting data
    ! (UNUSED) 2: Do not load/use any existing restarting data but save all new restarting data 
    ! (overwriting any existing restarting data)
    ! 3: Use any existing restarting data and extend existing restarting data with all new restarting data
    ! Host program must tell setup to use - no default choice in OpenRSP core
    integer :: r_flag
    
    logical :: r_exist, sdf_retrieved
    integer, dimension(3) :: rs_info, rs_calibrate_save
    integer, dimension(3) :: prog_info

    call QcMatInit(S_unpert_arr(1))
    call QcMatAEqB(S_unpert_arr(1), S_unpert)
    call QcMatInit(D_unpert_arr(1))
    call QcMatAEqB(D_unpert_arr(1), D_unpert)    
    call QcMatInit(F_unpert_arr(1))
    call QcMatAEqB(F_unpert_arr(1), F_unpert)    


    xf_was_allocated = .FALSE.
    
    r_flag = r_flag_in
    
    if (present(mem_calibrate)) then
    
       mem_mgr%calibrate = mem_calibrate

       if (mem_mgr%calibrate) then
       
          if (.NOT.(present(mem_result))) then

             write(out_str, *) 'ERROR: Result holder "mem_result" must be given for memory calibration run'
             call out_print(out_str, -1)

             return
             
          end if
          
       end if
          
    else
    
       mem_mgr%calibrate = .FALSE.
       
    end if
    
    
    if (present(max_mat)) then
       
       mem_mgr%max_mat = max_mat
       mem_mgr%remain = max_mat
       mem_mgr%limited = .TRUE.
       
    else
       
       mem_mgr%limited = .FALSE.
       
    end if
    
    

    ! Start progress counter
    
    prog_info = (/0,0,0/)
    call prog_init(rs_info, r_flag)
    
    ! Circumvent restarting mechanism for calibration run
    ! by pretending that last run did not progress beyond the beginning
    if (mem_mgr%calibrate) then
    
       rs_calibrate_save = rs_info
       rs_info = (/0,0,0/)
    
    end if
    
    call prog_incr(prog_info, r_flag, 1)       

    ! Present calculation and initialize perturbation tuple datatypes and
    ! associated size/indexing information

    write(out_str, *) ' '
    call out_print(out_str, 1)
    
    if (residue_order == 0) then
    
       write(out_str, *) 'OpenRSP lib called for (non-residue) response property calculation'
       call out_print(out_str, 1)
    
    elseif (residue_order == 1) then

       write(out_str, *) 'OpenRSP lib called for response property single residue calculation'
       call out_print(out_str, 1)
    
    elseif (residue_order == 2) then
    
       write(out_str, *) 'ERROR: OpenRSP lib called for double residue calculation'
       call out_print(out_str, 0)
       write(out_str, *) 'The only currently supported residues are single residues'
       call out_print(out_str, 0)
       write(out_str, *) 'Cannot proceed with calculation: Exiting OpenRSP library'
       call out_print(out_str, -1)
    
       return
    
    else
    
       write(out_str, *) 'ERROR: OpenRSP lib called for unsupported order of residue'
       call out_print(out_str, 0)
       write(out_str, *) 'The only currently supported residues are single residues'
       call out_print(out_str, 0)
       write(out_str, *) 'Cannot proceed with calculation: Exiting OpenRSP library'
       call out_print(out_str, -1)
    
       return
    
    end if
    
    write(out_str, *) ' '
    call out_print(out_str, 1)
    
    
    if (n_props == 1) then
       write(out_str, *) 'Calculating one property'
       call out_print(out_str, 1)
    else  
       write(out_str, *) 'Calculating', n_props, 'properties'
       call out_print(out_str, 1)
    end if

    write(out_str, *) ' '
    call out_print(out_str, 1)
   
    k = 1
    
    do i = 1, n_props

       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Property', i, 'is order', np(i)
       call out_print(out_str, 1)
       write(out_str, *) 'The choice of k, n is', kn_rules(i), 'and', np(i) - 1 - kn_rules(i)
       call out_print(out_str, 1)
       write(out_str, *) 'The number of components for each perturbation is: ', &
                         pert_dims(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
       call out_print(out_str, 1)
       write(out_str, *) 'The perturbation labels are: ', pert_labels(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
       call out_print(out_str, 1)
       write(out_str, *) 'Number of frequency configurations:', n_freq_cfgs(i)
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       
       if (residue_order > 0) then

          write(out_str, *) 'Note: One perturbation in this tuple plays the role of residue placeholder'
          call out_print(out_str, 1)
          
       end if
       
    
       do j = 1, n_freq_cfgs(i)
       
          kn_rule(k, 1) = kn_rules(i)
          kn_rule(k, 2) = np(i) - 1 - kn_rules(i)
          
          if ((kn_rule(k, 1) - kn_rule(k, 2) > 1)) then

             write(out_str, *) 'ERROR: Invalid choice of (k,n)'
             call out_print(out_str, 0)
             write(out_str, *) 'Valid choices for k are integers between and including 0 and ', (np(i) - 1)/2
             call out_print(out_str, 0)
             write(out_str, *) 'Valid choices of n are such that k + n =', np(i) - 1
             call out_print(out_str, 0)
             write(out_str, *) 'Cannot proceed with calculation: Exiting OpenRSP lib'
             call out_print(out_str, -1)
                          
             return
 
          end if
         
          p_tuples(k)%npert = np(i)
          allocate(p_tuples(k)%pdim(np(i)))
          allocate(p_tuples(k)%plab(np(i)))
          allocate(p_tuples(k)%pid(np(i)))
          allocate(p_tuples(k)%freq(np(i)))
          
          p_tuples%do_residues = residue_order
          
          if (residue_order > 0) then
          
             ! DaF: Initialization of residue-relevant parts of the perturbation tuple:
             allocate(p_tuples(k)%part_of_residue(p_tuples(k)%npert,residue_order))
             allocate(p_tuples(k)%exenerg(residue_order))
             allocate(p_tuples(k)%states(residue_order))
             p_tuples(k)%part_of_residue = .false.
             ! DaF: So far we only accept residue calculation where excitation energies
             ! match single perturbation frequencies, and not frequency sums
             ! Loop over the number of perturbations which contribute to residualization
             do m = 1, residue_order
                p_tuples(k)%exenerg(m) = exenerg(m)
                ! MaR: NOTE that the RHS of the below line used the (uninitialized)
                ! 'residualization' variable. For now set to zero as this line will not
                ! be relevant until possible use for double residues.
                p_tuples(k)%states(m) = 0

                do l = 1, max(residue_spec_pert(1),residue_spec_pert(residue_order))
                
                   ! MaR: Treating residue_spec_index as collapsed, revisit if ordering issues
                
                   p_tuples(k)%part_of_residue(residue_spec_index( &
                   (m - 1) * max(residue_spec_pert(1),residue_spec_pert(residue_order)) + l ) , m) = .true.
                   
                   
                end do
             end do
          
          end if
          
          p_tuples(k)%pdim = pert_dims(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
          p_tuples(k)%plab = pert_labels(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
          p_tuples(k)%pid = (/(m, m = 1, np(i))/)
          p_tuples(k)%freq = pert_freqs(dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + &
          1 + (j - 1)*np(i):dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + (j)*np(i))
          
          if (residue_order > 0) then
          
             ! DaF: Does the residualized perturbation have a proper frequency?
             lfreq_match = .false.
             do l = 1, p_tuples(k)%npert
                if (dabs(dabs(dble(p_tuples(k)%exenerg(1))) - dabs(dble(p_tuples(k)%freq(l)))).gt.xtiny) then
                   lfreq_match = .true.
                end if
             end do
             
             if (.not.lfreq_match) then
             
                write(out_str, *) 'ERROR: No perturbation frequencies matched excitation energy'
                call out_print(out_str, 0)
                write(out_str, *) 'The residue calculation is therefore indeterminate'
                call out_print(out_str, 0)
                write(out_str, *) 'Cannot proceed with calculation: Exiting OpenRSP library'
                call out_print(out_str, -1)

                return             
                         
             end if
          
          end if
          
          ! MaR: NOTE: Here sorting tuples in standard order: Must sort back to get correct order 
          ! for returned tensor: Original pids will give key to re-sort
          
          p_tuples(k) = p_tuple_standardorder(p_tuples(k))
          p_tuples(k)%pid = (/(m, m = 1, np(i))/)
       
          write(out_str, *) 'Frequency configuration', j
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'Frequencies (real part):', (/(real(p_tuples(k)%freq(m)), m = 1, np(i))/)
          call out_print(out_str, 1)
          write(out_str, *) 'Frequencies (imag. part):', (/(aimag(p_tuples(k)%freq(m)), m = 1, np(i))/)
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
          
          
          ! MaR: Rounding down nearly zero frequencies to zero
          ! This is done to avoid size and indexing issues for cases where a perturbation which should
          ! be zero was wrongly set to nearly zero due to float precision issues
          
          do m = 1, np(i)
          
             if (abs(real(p_tuples(k)%freq(m))) < xtiny) then
             
               write(out_str, *) 'NOTE: The real part of frequency', m, 'of property ', i, 'is zero or nearly zero'
               call out_print(out_str, 1)
               write(out_str, *) 'The value is ', real(p_tuples(k)%freq(m))
               call out_print(out_str, 1)
               write(out_str, *) 'Manually set to zero to avoid tensor size and indexing issues'
               call out_print(out_str, 1)
               write(out_str, *) ' '
               call out_print(out_str, 1)
             
               p_tuples(k)%freq(m) = cmplx(0.0d0, aimag(p_tuples(k)%freq(m)), kind=8)
             
             end if
             
             if (abs(aimag(p_tuples(k)%freq(m))) < xtiny) then
             
               write(out_str, *) 'NOTE: The imaginary part of frequency', m, 'of property ', i, 'is zero or nearly zero'
               call out_print(out_str, 1)
               write(out_str, *) 'The value is ', aimag(p_tuples(k)%freq(m))
               call out_print(out_str, 1)
               write(out_str, *) 'Manually set to zero to avoid tensor size and indexing issues'
               call out_print(out_str, 1)
               write(out_str, *) ' '
               call out_print(out_str, 1)
             
               p_tuples(k)%freq(m) = cmplx(real(p_tuples(k)%freq(m)), 0.0d0, kind=8)
             
             end if
          
          
          end do
          
          
       
                    
          num_blks(k) = get_num_blks(p_tuples(k))
       
          write(out_str, *) 'Number of blocks:', num_blks(k)
          call out_print(out_str, 2)
       
          allocate(blk_info(num_blks(k), 3))
          allocate(blk_sizes(num_blks(k)))
          blk_info = get_blk_info(num_blks(k), p_tuples(k))
          
          write(out_str, *) 'Block info:', blk_info
          call out_print(out_str, 2)
          
          blk_sizes = get_triangular_sizes(num_blks(k), blk_info(1:num_blks(k), 2), &
                                           blk_info(1:num_blks(k), 3))

          write(out_str, *) 'Block sizes:', blk_sizes
          call out_print(out_str, 2)
          
          prop_sizes(k) = get_triangulated_size(num_blks(k), blk_info)
          
          write(out_str, *) 'Property size for this frequency configuration:', prop_sizes(k)
          call out_print(out_str, 1)
          
          deallocate(blk_info)
          deallocate(blk_sizes)
       
          k = k + 1
       
      
       end do

              
    end do

        
    
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (mem_mgr%calibrate) then
       
       call mem_incr(mem_mgr, 3 + residue_order, p=prog_info)
    
    else
    
       ! Check if this stage passed previously and if so, then retrieve and skip execution
    
       sdf_retrieved = .FALSE.
       if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then

          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'S, D, F initialization stage was completed'
          call out_print(out_str, 1)
          write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)


!           call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE.)
!           call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE.)
!           call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE.)
          
          if (residue_order > 0) then
       
!              allocate(Xf)
!              call contrib_cache_outer_retrieve(Xf, 'OPENRSP_Xf_CACHE', .FALSE.)
          
          end if
          
          sdf_retrieved = .TRUE.
      
       else
       
          ! Set up S, D, F data structures
       
         call contrib_cache_outer_allocate(S)
         call contrib_cache_outer_allocate(D)
         call contrib_cache_outer_allocate(F)

         call empty_p_tuple(empty_pert(1))

         call contrib_cache_outer_add_element(size(S), S, .TRUE., 1, empty_pert, &
              data_size = 1, data_mat=S_unpert)
         call contrib_cache_outer_add_element(size(D), D, .TRUE., 1, empty_pert, &
              data_size = 1, data_mat=D_unpert)
         call contrib_cache_outer_add_element(size(F), F, .TRUE., 1, empty_pert, &
              data_size = 1, data_mat=F_unpert)
               
!          call contrib_cache_outer_store(S, 'OPENRSP_S_CACHE', r_flag)
!          call contrib_cache_outer_store(D, 'OPENRSP_D_CACHE', r_flag)
!          call contrib_cache_outer_store(F, 'OPENRSP_F_CACHE', r_flag)
         
         if (residue_order > 0) then

    call QcMatInit(Xf_unpert_arr(1))
    call QcMatAEqB(Xf_unpert_arr(1), Xf_unpert(1))    
       
             call contrib_cache_outer_allocate(Xf)
             xf_was_allocated = .TRUE.
             call contrib_cache_outer_add_element(size(Xf), Xf, .TRUE., 1, empty_pert, &
                  data_size = 1, data_mat=Xf_unpert(1))
!              call contrib_cache_outer_store(Xf,'OPENRSP_Xf_CACHE', r_flag)
          
          end if
       
       end if
       
    end if
    

    rsp_tensor(1:sum(prop_sizes)) = 0.0
    
    if (.NOT.(mem_exceed(mem_mgr))) then
    
       ! Calculate properties

       write(out_str, *) 'Starting clock: About to start property calculations'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       call cpu_time(timing_start)

       if (residue_order > 0) then
       
          len_fds = size(F)
          len_x = size(Xf)
       
          call get_prop(n_props, n_freq_cfgs, p_tuples, kn_rule, len_fds, F, D, S, get_rsp_sol, &
                        get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                        get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, out_print, &
                        prop_sizes, rsp_tensor, prog_info, rs_info, r_flag, sdf_retrieved, &
                        mem_mgr, Xf=Xf)
                        
       else
       
          len_fds = size(F)
       
          call get_prop(n_props, n_freq_cfgs, p_tuples, kn_rule, len_fds, F, D, S, get_rsp_sol, &
                        get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                        get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, out_print, &
                        prop_sizes, rsp_tensor, prog_info, rs_info, r_flag, sdf_retrieved, &
                        mem_mgr)
                        
       end if

       call cpu_time(timing_end)

       write(out_str, *) 'Clock stopped: Property calculations finished'
       call out_print(out_str, 1)
       write(out_str, *) 'Time spent:',  timing_end - timing_start, ' seconds'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
   
    end if
    
    
    if (mem_mgr%calibrate) then
    
       call mem_decr(mem_mgr, 3)
       mem_result = mem_mgr%status
       
       ! Put back restart information the way it was
       if (.NOT.(all(rs_calibrate_save == (/0,0,0/)))) then
          
          rs_calibrate_save(1) = rs_calibrate_save(1) - 1
          call prog_incr(rs_calibrate_save, r_flag, 1)
          
       end if
    
    else
    
    
       ! MaR: FIXME: Adapt to include residues
    
       write(out_str, *) 'Writing response tensors to file'
       call out_print(out_str, 2)
    
       ! Print tensors to standardized output file
       ! NOTE: This routine is placed here due to access to index/addressing routines
       ! Should be moved to API level once index/addressing routines are abstracted
       
       open(unit=260, file='rsp_tensor', status='replace', action='write') 
       
       write(260,*) 'VERSION'
       write(260,*) '1'
       write(260,*) 'NUM_PROPERTIES'
       write(260,*) n_props

   
       
       k = 1
       p = 0
       
       do i = 1, n_props
       
          write(260,*) 'NEW_PROPERTY'
          write(260,*) 'ORDER'
          write(260,*) p_tuples(k)%npert
          write(260,*) 'NUM_FREQ_CFGS'
          write(260,*) n_freq_cfgs(i)
          
          write(260,*) 'OPERATORS'
          do j = 1, p_tuples(k)%npert
          
             write(260,*) p_tuples(k)%plab(j)          
          
          end do
          
          
          write(260,*) 'NUM_COMPONENTS'
          do j = 1, p_tuples(k)%npert
          
             write(260,*) p_tuples(k)%pdim(j)          
          
          end do
         
          write(260,*) 'FREQUENCIES'
          
         
          do j = 1, n_freq_cfgs(i)
          
             write(260,*) 'CONFIGURATION'
             
             do n = 1, p_tuples(k)%npert
             
                write(260,*) real(p_tuples(k)%freq(n))
             
             end do
             
         
          
             k = k + 1
         
          end do
          
          k = k - n_freq_cfgs(i)
         
          write(260,*) 'VALUES'
         
         ! FIXME: SOMETHING MAY BE OFF ABOUT THE INDICES: NOT ALL VALUES OF THE LAST PROPERTY ARE WRITTEN
         ! MaR: Update: Now likely fixed and above comment was forgotten, keep in case further problems
         
          do j = 1, n_freq_cfgs(i)
          
             write(260,*) 'CONFIGURATION'
             
             ! Get indices, write index-value pairs
             allocate(blk_info(num_blks(k), 3))
             allocate(blk_sizes(num_blks(k)))
             blk_info = get_blk_info(num_blks(k), p_tuples(k))

             blk_sizes = get_triangular_sizes(num_blks(k), blk_info(1:num_blks(k), 2), &
                                              blk_info(1:num_blks(k), 3))

             
             allocate(indices(product(blk_sizes), sum(blk_info(:,2))))
             
             call make_triangulated_indices(num_blks(k), blk_info, &
                  product(blk_sizes), indices)
                  
                  
             do n = 1, size(indices, 1)
             
                ! NOTE: TENSOR ELEMENTS WITH ABSOLUTE VALUE BELOW write_threshold WILL NOT BE OUTPUT
                if (abs(real(rsp_tensor(p + n))) >= write_threshold) then
                
                   write(260,*) indices(n,:)
                   write(260,*) real(rsp_tensor(p + n))
                
                end if
             
             
             end do
             
             deallocate(indices)
             deallocate(blk_info)
             deallocate(blk_sizes)
             
             p = p + prop_sizes(k)
             k = k + 1
          
          end do
          
       end do
       
       close(260)
       
    end if
    
    
    deallocate(F)
    deallocate(D)
    deallocate(S)
    if (xf_was_allocated) then
       deallocate(Xf)
    end if
    
    write(out_str, *) 'OpenRSP library: Normal termination, returning...'
    call out_print(out_str, 1)

  end subroutine
    
    
    
    
    
    
   
  ! Main property calculation routine - Get perturbed F, D, S and then calculate the properties
  subroutine get_prop(n_props, n_freq_cfgs, p_tuples, kn_rule, len_fds, F, D, S, get_rsp_sol, &
                      get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                      get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, out_print, &
                      prop_sizes, props, prog_info, rs_info, r_flag, sdf_retrieved, &
                      mem_mgr, Xf)

    implicit none

    type(mem_manager) :: mem_mgr
    logical :: traverse_end, sdf_retrieved, contrib_retrieved, props_retrieved
    integer :: n_props, i, j, k
    integer :: r_flag
    integer :: len_fds
    integer, dimension(3) :: prog_info, rs_info
    integer, dimension(n_props) :: n_freq_cfgs
    integer, dimension(sum(n_freq_cfgs)) :: prop_sizes
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    integer :: len_cache
    type(p_tuple) :: emptypert
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), dimension(2) :: emptyp_tuples
    external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
    external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
    complex(8), dimension(*) :: props
    
    type(contrib_cache_outer), allocatable, dimension(:) :: F, D, S
    type(contrib_cache_outer), optional, allocatable, dimension(:) :: Xf
    
    type(contrib_cache), allocatable, dimension(:) :: contribution_cache
    
    external :: out_print
    character(len=1048576) :: out_str
    
    call empty_p_tuple(emptypert)
    emptyp_tuples = (/emptypert, emptypert/)

    call prog_incr(prog_info, r_flag, 1)
  
    ! Check if this stage passed previously and if so, then retrieve and skip execution
  
    contrib_retrieved = .FALSE.
    props_retrieved = .FALSE.
  
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Perturbed overlap/density/Fock matrix stage was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       if(.NOT.(sdf_retrieved)) then
       
!           call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE.)
!           call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE.)
!           call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE.)
          
          if (present(Xf)) then
       
!           call contrib_cache_outer_retrieve(Xf, 'OPENRSP_Xf_CACHE', .FALSE.)
          
          end if
                 
       end if
       
       sdf_retrieved = .TRUE.
  
    else
  
       ! Get all necessary F, D, S derivatives
     
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Calculating perturbed overlap/density/Fock matrices'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)

       call cpu_time(time_start)
       
       
       if (present(Xf)) then
       
          call rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, size(D), F, D, S, &
                       get_rsp_sol, get_ovl_mat, get_1el_mat, &
                       get_2el_mat, get_xc_mat, out_print, .TRUE., &
                       prog_info, rs_info, r_flag, sdf_retrieved, mem_mgr, Xf=Xf)
                       
                       
       else
       
          call rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, size(D), F, D, S, &
                       get_rsp_sol, get_ovl_mat, get_1el_mat, &
                       get_2el_mat, get_xc_mat, out_print, .TRUE., &
                       prog_info, rs_info, r_flag, sdf_retrieved, mem_mgr)
       
       end if
                    
       call cpu_time(time_end)
       
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Finished calculation of perturbed overlap/density/Fock matrices'
       call out_print(out_str, 1)
       write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       
       if (mem_exceed(mem_mgr)) then
       
          return
          
       end if
       
    end if
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'HF energy-type contribution identification was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       if (.NOT.(contrib_retrieved)) then
         
!           call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
    
    else
    
       ! For each property and freq. cfg.: Recurse to identify HF energy-type contributions, store in cache
    
       call contrib_cache_allocate(contribution_cache)
    
       len_cache = size(contribution_cache)
    
       k = 1

       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)
       
             
             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Identifying HF-energy type contributions'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
             
             len_cache = size(contribution_cache)

             call cpu_time(time_start)
             call rsp_energy_recurse(p_tuples(k), p_tuples(k)%npert, kn_rule(k,:), 1, (/emptypert/), &
                  0, size(D), D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, out_print, .TRUE., &
                  len_cache, contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
             call cpu_time(time_end)

write(*,*) 'stage 1'

             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Finished identifying HF energy-type contributions'
             call out_print(out_str, 1)
             write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
             
             k = k + 1
       
          end do
       
       end do
    
       len_cache = size(contribution_cache)
    
!        call contrib_cache_store(len_cache, contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
    
    end if
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'HF energy-type contribution calculation was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       if (.NOT.(contrib_retrieved)) then
       
!           allocate(contribution_cache)
    
!           call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
           
    else
    
       ! Calculate all identified contributions and store in cache
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Calculating HF-energy type contributions'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       call cpu_time(time_start)
    
       len_cache = size(contribution_cache)
    
       ! Traverse linked list and calculate
       do k = 1, len_cache
       
          if (contribution_cache(k)%p_inner%npert == 0) then
          
             cycle
             
          end if   
  
          ! MaR: Don't know where this line comes from, reinstate if needed
          ! if(contribution_cache(k)%p_inner%plab(1).eq.'NUTN') stop 'empty perturbation in rsp_general!'

          write(out_str, *) 'Calculating contribution for inner perturbation tuple with labels:'
          call out_print(out_str, 1)
          write(out_str, *) contribution_cache(k)%p_inner%plab
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
       
          ! Check if this stage passed previously and if so, then skip execution
          if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
          
             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Calculation was completed in previous invocation: Passing to next stage'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
                
             ! Note: No cache retrieval here: In order to get to this position, the
             ! cache would already have been retrieved
          
          else

write(*,*) 'stage a'

             call rsp_energy_calculate(size(D), D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, &
                  out_print, contribution_cache(k), mem_mgr)
             
             if (mem_exceed(mem_mgr)) then
       
                return
          
             end if
          
!              call contrib_cache_store(contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
             
          end if
          
          call prog_incr(prog_info, r_flag, 2)
          
       end do

       call cpu_time(time_end)

       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Finished calculating HF energy-type contributions'
       call out_print(out_str, 1)
       write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
!        call contrib_cache_store(contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
    
    end if
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'HF energy-type contribution assembly was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       if (.NOT.(props_retrieved)) then
       
!           call mat_scal_retrieve(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
          props_retrieved = .TRUE.
          
       end if
           
    else
    
       if (.NOT.(mem_mgr%calibrate)) then
    
          ! For each property and freq. cfg.: Recurse to identify HF energy-type contributions and 
          ! add to the property under consideration
   
          k = 1
    
          do i = 1, n_props
    
             do j = 1, n_freq_cfgs(i)

                write(out_str, *) ' '
                call out_print(out_str, 1)
                write(out_str, *) 'Assembling HF-energy type contributions'
                call out_print(out_str, 1)
                write(out_str, *) ' '
                call out_print(out_str, 1)

                call cpu_time(time_start)
                call rsp_energy_recurse(p_tuples(k), p_tuples(k)%npert, kn_rule(k,:), 1, (/emptypert/), &
                     0, size(D), D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, out_print, .FALSE., &
                     len_cache, contribution_cache, prop_sizes(k), &
                     props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
                call cpu_time(time_end)
                
                write(out_str, *) ' '
                call out_print(out_str, 1)
                write(out_str, *) 'Finished assembling HF energy-type contributions'
                call out_print(out_str, 1)
                write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
                call out_print(out_str, 1)
                write(out_str, *) ' '
                call out_print(out_str, 1)

                write(out_str, *) 'Property sample', props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1: &
                min(sum(prop_sizes(1:k)) - prop_sizes(k) + 100, sum(prop_sizes(1:k))))
                call out_print(out_str, 2)
                
                k = k + 1
       
             end do
        
          end do
       
!           call mat_scal_store(sum(prop_sizes), 'OPENRSP_PROP_CACHE', r_flag, scal=props)
       
       end if

    end if
    
    deallocate(contribution_cache)
    
    contrib_retrieved = .FALSE.
    
    
    ! NEW: XC contribution
    
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
   
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'XC contributions were completed in previous invocation:'
       call out_print(out_str, 1)
       write(out_str, *) 'Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
           
       if (.NOT.(props_retrieved)) then
       
!           call mat_scal_retrieve(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
          props_retrieved = .TRUE.
          
       end if
           
    else
    
       if (.NOT.(mem_mgr%calibrate)) then
    
          ! For each property: Set up XC contribution calculation, calculate and 
          ! add to the property under consideration
   
          k = 1
    
          do i = 1, n_props
          
             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Calculating XC contributions'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)

             call cpu_time(time_start)
             call rsp_xc_wrapper(n_freq_cfgs(i), p_tuples(k:k+n_freq_cfgs(i)-1), kn_rule(k,:), &
                  size(D), D, get_xc_exp, out_print, sum(prop_sizes(k:k+n_freq_cfgs(i) - 1)), mem_mgr, &
                  prop=props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1 : &
                        sum(prop_sizes(1:k+n_freq_cfgs(i) - 1))))
             call cpu_time(time_end)

             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Finished calculating XC contributions'
             call out_print(out_str, 1)
             write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)

             
             write(out_str, *) 'Property sample', props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1: &
             min(sum(prop_sizes(1:k)) - prop_sizes(k) + 100, sum(prop_sizes(1:k))))
             call out_print(out_str, 2)
             
             k = k + n_freq_cfgs(i)
       
          end do
       
!           call mat_scal_store(sum(prop_sizes), 'OPENRSP_PROP_CACHE', r_flag, scal=props)
       
       end if

    end if
    
    
    
    ! END NEW: XC contribution
    
    
    
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Two-factor type contribution identification was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       if (.NOT.(contrib_retrieved)) then
       
!           allocate(contribution_cache)
    
!           call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
           
    else
    
       ! For each property and freq. cfg.: Recurse to identify two-factor contributions and store in cache
    
       k  = 1
    
       call contrib_cache_allocate(contribution_cache)
       len_cache = size(contribution_cache)
    
       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)
          
             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Identifying two-factor contributions'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
          
             call cpu_time(time_start)
          
             call rsp_twofact_recurse(p_tuples(k), &
                  kn_rule(k,:), (/emptypert, emptypert/), out_print, &
                  .TRUE., len_cache, contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
             call cpu_time(time_end)

             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Finished identifying two-factor contributions'
             call out_print(out_str, 1)
             write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
             
             k = k + 1
       
          end do
       
       end do
    
       len_cache = size(contribution_cache)
!        call contrib_cache_store(len_cache, contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
    
    end if
    
    ! Check if this stage passed previously and if so, then retrieve and skip execution
    
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Two-factor type contribution calculation was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
       if (.NOT.(contrib_retrieved)) then
       
!           allocate(contribution_cache)
    
!           call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
       
       
    else
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Calculating two-factor contributions'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)

       call cpu_time(time_start)
    
       len_cache = size(contribution_cache)
    
       ! Traverse linked list and calculate
       do k = 1, len_cache
       
          if (contribution_cache(k)%p_inner%npert == 0) then
          
             cycle
             
          end if   
       
          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'Calculating contribution for factor 1 tuple perturbation labels:'
          call out_print(out_str, 1)
          write(out_str, *) contribution_cache(k)%p_inner%plab
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
       
          ! Check if this stage passed previously and if so, then skip execution
          if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
             
             write(out_str, *) ' '
             call out_print(out_str, 1)
             write(out_str, *) 'Calculation was completed in previous invocation: Passing to next stage'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
          
             ! Note: No cache retrieval here: In order to get to this position, the
             ! cache would already have been retrieved
          
          else

             call rsp_twofact_calculate(size(S), S, size(D), D, size(F), F, get_ovl_exp, &
                                        out_print, contribution_cache(k), mem_mgr)
             
             if (mem_exceed(mem_mgr)) then
       
                return
          
             end if             
             
!              call contrib_cache_store(contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
          
          end if
          
          call prog_incr(prog_info, r_flag, 2)
          
       end do

       call cpu_time(time_end)

       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Finished calculating two-factor contributions'
       call out_print(out_str, 1)
       write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
!        call contrib_cache_store(contribution_cache, r_flag, 'OPENRSP_CONTRIB_CACHE')
       
    end if

    ! Check if this stage passed previously and if so, then retrieve and skip execution
    call prog_incr(prog_info, r_flag, 1)
    
    if (rs_check(prog_info, rs_info, r_flag, lvl=1)) then
    
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'Two-factor type contribution assembly was completed'
       call out_print(out_str, 1)
       write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       if (.NOT.(props_retrieved)) then
       
!           call mat_scal_retrieve(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
          props_retrieved = .TRUE.
       
       end if
       
    else
    
       if (.NOT.(mem_mgr%calibrate)) then
    
          ! For each property and freq. cfg.: Recurse to identify two-factor contributions and 
          ! add to the property under consideration
       
          k = 1
    
          do i = 1, n_props
    
             do j = 1, n_freq_cfgs(i)

             
                write(out_str, *) ' '
                call out_print(out_str, 1)
                write(out_str, *) 'Assembling two-factor contributions'
                call out_print(out_str, 1)
                write(out_str, *) ' '
                call out_print(out_str, 1)
          
                len_cache = size(contribution_cache)
          
                call cpu_time(time_start)
          
                call rsp_twofact_recurse(p_tuples(k), &
                     kn_rule(k,:), (/emptypert, emptypert/), out_print, &
                     .FALSE., len_cache, contribution_cache, prop_sizes(k), &
                     props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
                call cpu_time(time_end)
                
                write(out_str, *) ' '
                call out_print(out_str, 1)
                write(out_str, *) 'Finished assembling two-factor contributions'
                call out_print(out_str, 1)
                write(out_str, *) 'Time spent:', time_end - time_start, 'seconds'
                call out_print(out_str, 1)
                write(out_str, *) ' '
                call out_print(out_str, 1)

                k = k + 1
          
             end do
       
          end do
       
!           call mat_scal_store(sum(prop_sizes), 'OPENRSP_PROP_CACHE', r_flag, scal=props)
       
       end if
  
    end if
  
    deallocate(contribution_cache)
    
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       k = 1
    
       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)
          
             write(out_str, *) ' '
             call out_print(out_str, 2)
             write(out_str, *) 'Property', i, ', freq. config', j
             call out_print(out_str, 2)
             write(out_str, *) props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1: &
                               sum(prop_sizes(1:k)))
             call out_print(out_str, 2)
             write(out_str, *) ' '
             call out_print(out_str, 2)
          
             k = k + 1
          
          end do
       
       end do
       
    end if
    
    call prog_incr(prog_info, r_flag, 1)
   
    
  end subroutine
   
   ! Recurse to identify (dryrun == .TRUE.) or assemble (dryrun == .FALSE.) energy-type contributions
   
   recursive subroutine rsp_energy_recurse(pert, total_num_perturbations, kn, num_p_tuples, &
                                  p_tuples, density_order, len_d, D, get_nucpot, get_1el_exp, &
                                  get_t_exp, get_2el_exp, out_print, dryrun, len_cache, cache, p_size, prop)

    implicit none

    logical :: e_knskip, dryrun, traverse_end
    type(p_tuple) :: pert
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations
    integer :: p_size
    integer :: len_d, len_cache
    logical :: residue_skip
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache), allocatable, dimension(:) :: cache
    complex(8), dimension(p_size), optional :: prop
    external :: get_nucpot, get_1el_exp, get_t_exp, get_2el_exp
    
    external :: out_print
    character(len=2047) :: out_str


    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'


    if (p_tuples(1)%npert == 0) then

       call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples, (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
       density_order, len_d, D, get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, out_print, &
       dryrun, len_cache, cache, p_size=p_size, prop=prop)

    else

       call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations,  &
       kn, num_p_tuples, (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
       p_tuples(2:size(p_tuples))/), density_order, len_D, D, &
       get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, out_print, &
       dryrun, len_cache, cache, p_size=p_size, prop=prop)

    end if
    
       ! 2. Differentiate all of the contraction densities in turn

       ! Find the number of terms

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%npert == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples, t_new, density_order + 1, len_d, D, &
          get_nucpot, get_1el_exp, get_t_exp, get_2el_exp , out_print, &
          dryrun, len_cache, cache, p_size=p_size, prop=prop)

       end do

       ! Since we are only calculating Hartree-Fock type energy terms here,
       ! we don't need to go beyond to perturbed contraction density matrices
       ! (but that is in general needed for XC contribs)
       if (num_p_tuples < 3) then

          ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
          ! a(nother) pert D contraction)

          call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
          density_order + 1, len_d, D, get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, out_print, &
          dryrun, len_cache, cache, p_size=p_size, prop=prop)

       end if

    ! At the final recursion level: Calculate the contrib (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allows it)

    else

       p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)
    
       e_knskip = .FALSE.
       residue_skip = (pert%do_residues.gt.0).and.find_residue_info(p_tuples(1))

       do i = 1, num_p_tuples
 
          if (i > 1) then

             if(kn_skip(p_tuples(i)%npert, p_tuples(i)%pid, kn)) then

                e_knskip = .TRUE.

             end if

          
          elseif (i == 1) then
          
          

          end if

       end do


       if ((e_knskip .EQV. .FALSE.) .AND. (residue_skip .EQV. .FALSE.)) then
       
          if (contrib_cache_already(len_cache, cache, num_p_tuples, p_tuples)) then
          
             if (.NOT.(dryrun)) then
             
                write(out_str, *) 'Cache retrieval: getting contribution'
                call out_print(out_str, 2)

                do i = 1, num_p_tuples

                   if (i == 1) then
    
                      write(out_str, *) 'E', p_tuples(i)%pid
                      call out_print(out_str, 2)
    
                   else 
                   
                      write(out_str, *) 'D', p_tuples(i)%pid
                      call out_print(out_str, 2)
                       
                   end if
    
                end do

                write(out_str, *) ' '
                call out_print(out_str, 2)
             
                ! NOTE (MaR): EVERYTHING MUST BE STANDARD ORDER IN 
                ! THIS CALL (LIKE property_cache_getdata ASSUMES)
                call contrib_cache_getdata(len_cache, cache, num_p_tuples, &
                   p_tuples_standardorder(num_p_tuples, p_tuples), p_size, 0, scal=prop)

             end if

          else

             if (dryrun) then
             
                write(out_str, *) 'Adding newly identified cache element'
                call out_print(out_str, 2)
             
                do i = 1, num_p_tuples

                   if (i == 1) then
                   
                      write(out_str, *) 'E', p_tuples(i)%pid
                      call out_print(out_str, 2)
    
                   else 
                    
                      write(out_str, *) 'D', p_tuples(i)%pid
                      call out_print(out_str, 2)
                       
                   end if
    
                end do
                
                write(out_str, *) ' '
                call out_print(out_str, 2)
             
             
                call contrib_cache_add_element(len_cache, cache, num_p_tuples, &
                     p_tuples_standardorder(num_p_tuples, p_tuples))

             end if

          end if

       end if

    end if

  end subroutine


  subroutine rsp_energy_calculate(len_d, D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, &
             out_print, cache, mem_mgr)

    implicit none

    type(mem_manager) :: mem_mgr
    integer :: mctr, mcurr, miter, msize, octr, mem_track_1, mem_track_2, mem_track
    integer :: dctr, contrib_offset, any_heavy, curr_remain, this_outer_size
    integer, dimension(3) :: curr_pickup, next_pickup
    integer, dimension(2) :: wunit_size, wunit_maxsize
    logical :: intra_pair, mem_done, intra_done
    integer :: cache_offset, i, j, k, m, n, p, offset
    integer :: c1_ctr, c2_ctr, lhs_ctr_1, lhs_ctr_2, rhs_ctr_2
    integer :: total_outer_size_1, total_outer_size_2
    integer :: num_0, num_1, num_pert, tot_num_pert
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
    character(30) :: mat_str, fmt_str
    
    integer :: len_d
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache) :: cache
    
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, LHS_dmat_2, RHS_dmat_2
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext, outer_contract_sizes_2_pair
    integer, allocatable, dimension(:,:) :: outer_contract_sizes_2, blk_sizes
    complex(8), allocatable, dimension(:) :: contrib_0, contrib_1, contrib_2, data_tmp, contrib_2_tmp
    external :: get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp
    
    external :: out_print
    character(len=1048576) :: out_str

    ! To avoid maybe-uninitialized warning
    any_heavy = 0
    
    ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)
    
    write(out_str, *) 'Number of outer contributions for this inner tuple:', cache%num_outer
    call out_print(out_str, 1)
   
    allocate(outer_contract_sizes_1(cache%num_outer))
    allocate(outer_contract_sizes_2(cache%num_outer,2))
   
    ! Traversal: Find number of density matrices for contraction for nuc-nuc, 1-el, 2-el cases
    
    total_outer_size_1 = 0
    total_outer_size_2 = 0
    num_0 = 0
    num_1 = 0
        
    k = 1
    
    do m = 1, size(cache%contribs_outer)
    
       if (cache%contribs_outer(m)%dummy_entry) then
       
          cycle
       
       end if
       
       write(out_str, *) 'Outer contribution', k
       call out_print(out_str, 1)
       
       do i = 1, cache%contribs_outer(m)%num_dmat
       
          write(out_str, *) 'D', cache%contribs_outer(m)%p_tuples(i)%pid
          call out_print(out_str, 1)
                 
       end do
       
       write(out_str, *) ' '
       call out_print(out_str, 1)
       
  
       ! If this is the case, then there are no outer perturbations (no chain rule applications)
       if (cache%contribs_outer(m)%num_dmat == 0) then

          num_0 = 1
          num_1 = num_1 + 1
          outer_contract_sizes_1(k) = 1
          outer_contract_sizes_2(k, :) = (/1,1/)
          
          total_outer_size_1 = total_outer_size_1 + 1
          total_outer_size_2 = total_outer_size_2 + 1
       
       ! One chain rule application: No nuc-nuc terms, only 1-el and 2-el terms
       else if (cache%contribs_outer(m)%num_dmat == 1) then
       
          num_1 = num_1 + 1
       
          outer_contract_sizes_1(k) = cache%contribs_outer(m)%blks_tuple_triang_size(1)
          outer_contract_sizes_2(k, :) = &
          (/cache%contribs_outer(m)%blks_tuple_triang_size(1),1/)
          
          total_outer_size_1 = total_outer_size_1 + &
          cache%contribs_outer(m)%blks_tuple_triang_size(1)
          total_outer_size_2 = total_outer_size_2 + &
          cache%contribs_outer(m)%blks_tuple_triang_size(1)
       
       ! Two chain rule applications: Only 2-el terms
       else if (cache%contribs_outer(m)%num_dmat == 2) then
       
          outer_contract_sizes_1(k) = 0
          outer_contract_sizes_2(k, :) = &
          (/cache%contribs_outer(m)%blks_tuple_triang_size(1), &
          cache%contribs_outer(m)%blks_tuple_triang_size(2)/)
          
          total_outer_size_2 = total_outer_size_2 + &
          cache%contribs_outer(m)%blks_tuple_triang_size(1) * &
          cache%contribs_outer(m)%blks_tuple_triang_size(2)
       
       end if
       
       k = k + 1
    
    end do
    
    ! Make collapsed contraction sizes array for 1-electron routine call
 
    allocate(outer_contract_sizes_1_coll(num_1))
    
    k = 1 
     do i = 1, cache%num_outer
        if (outer_contract_sizes_1(i) > 0) then
           outer_contract_sizes_1_coll(k) = outer_contract_sizes_1(i)
          k = k + 1
        end if
    end do
        
    
    allocate(contrib_0(cache%blks_triang_size))
    contrib_0 = 0.0
    
    
    
    ! Calculate contributions
    
    ! Calculate nuclear-nuclear repulsion contribution
    if (num_0 > 0) then
    
       write(out_str, *) 'Calculating density matrix-independent contribution'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       
       call get_nucpot(num_pert, pert_ext, size(contrib_0), contrib_0)
       
       write(out_str, *) 'Density matrix-independent contribution (sample)', &
       contrib_0(1:min(10,size(contrib_0)))
       call out_print(out_str, 2)
       
    end if
             
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    contrib_1 = 0.0
    
    
    if (mem_enough(mem_mgr, sum(outer_contract_sizes_1(:)))) then
    
       ! Possible to run in non-savings mode
       mem_track = sum(outer_contract_sizes_1(:))

    else if (mem_enough(mem_mgr, 1)) then

       ! Possible to run in savings mode
       call mem_set_status(mem_mgr, 1)
       
       ! Set maximum number of components to do at the same time
       mem_track = mem_mgr%remain
              
    else
    
       ! Not possible to run; flag it and return
       call mem_set_status(mem_mgr, 2)
       return
    
    end if
    
    ! Begin memory savings loop 1
    
    mcurr = 1
    
    do while (mcurr <= sum(outer_contract_sizes_1(:)))
    
       msize = min(mem_track, sum(outer_contract_sizes_1(:)) - mcurr + 1)
       
!        write(*,*) 'Memory loop for 1-el: Starting element is', mcurr, 'and block size is', msize
    
       ! Allocate and set up perturbed density matrices for contractions
       
       call mem_incr(mem_mgr, msize)
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(LHS_dmat_1(msize))
       
          do i = 1, msize
             call QcMatInit(LHS_dmat_1(i))
          end do
       
       end if
       

       k = 1
       lhs_ctr_1 = 1
       mctr = 0
       
       ! Traverse and fetch matrices
       do m = 1, size(cache%contribs_outer)
    
          if (cache%contribs_outer(m)%dummy_entry) then
       
             cycle
       
          end if
              
          ! No chain rule applications
          if (cache%contribs_outer(m)%num_dmat == 0) then
             
             if ((lhs_ctr_1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                if (.NOT.(mem_mgr%calibrate)) then
          
                   call contrib_cache_getdata_outer(size(D), D, 1, (/get_emptypert()/), .FALSE., &
                        1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_1(mctr + 1))
                     
                end if
                
                mctr = mctr + 1
             
             end if
    
          ! One chain rule application
          else if (cache%contribs_outer(m)%num_dmat == 1) then
          
          
             do n = 1, outer_contract_sizes_1(k) 
             
                if ((lhs_ctr_1 + n - 1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
               
                   if (.NOT.(mem_mgr%calibrate)) then
               
                      call contrib_cache_getdata_outer(size(D), D, 1, &
                      (/cache%contribs_outer(m)%p_tuples(1)/), .FALSE., &
                      1, ind_len=size(cache%contribs_outer(m)%indices, 2), &
                      ind_unsorted=cache%contribs_outer(m)%indices(n, :), &
                      mat_sing=LHS_dmat_1(mctr + 1))
                           
                   end if
                           
                   mctr = mctr + 1
                   
                end if
                
             end do
          
          end if
    
          lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
          k = k + 1
   
       end do
       
       
       
       ! Calculate one-electron contributions
       if (num_1 > 0) then
       
          if (.NOT.(mem_mgr%calibrate)) then
          
             write(out_str, *) 'Calculating first-order density matrix-dependent contribution'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
          
             call get_1el_exp(num_pert, pert_ext, msize, &
                              LHS_dmat_1, cache%blks_triang_size*msize, &
                              contrib_1( (mcurr - 1) * cache%blks_triang_size + 1: &
                                         (mcurr + msize - 1) * cache%blks_triang_size))
         
             t_matrix_bra = get_emptypert()
             t_matrix_ket = get_emptypert()
   
          
      ! NOTE: T matrix contributions not reinstated yet, to be done later
      !
      !        call rsp_ovlave_t_matrix_2014(get_ovl_exp, cache%p_inner, cache%p_inner%npert, &
      !                                 t_matrix_bra, t_matrix_ket, msize, &
      !                                 LHS_dmat_1, cache%blks_triang_size*msize, &
      !                                 contrib_1( (mcurr - 1) * cache%blks_triang_size + 1: &
      !                                 (mcurr + msize - 1) * cache%blks_triang_size)
       
       
             write(out_str, *) 'First-order density matrix-dependent contribution (sample)', &
             contrib_1(1:min(10,size(contrib_1)))
             call out_print(out_str, 2)
              
          end if
       
       end if
       

       mcurr = mcurr + msize
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          do i = 1, msize
       
             call QcMatDst(LHS_dmat_1(i))
       
          end do
       
          deallocate(LHS_dmat_1)
          
       end if
       
       call mem_decr(mem_mgr, msize)
       
    end do
       
    ! End memory savings loop 1
       
    
    allocate(contrib_2(cache%blks_triang_size*total_outer_size_2))
    contrib_2 = 0.0
    
    if (mem_enough(mem_mgr, sum(outer_contract_sizes_2(:, 1)) + &
                            sum(outer_contract_sizes_2(:, 2)))) then
    
       ! Possible to run in non-savings mode
       mem_track = sum(outer_contract_sizes_2(:, 1)) + sum(outer_contract_sizes_2(:, 2))
       
    elseif (mem_enough(mem_mgr, 2)) then

       ! Possible to run in savings mode
       call mem_set_status(mem_mgr, 1)
       
       ! Set maximum number of components to do at the same time
       mem_track = mem_mgr%remain
              
    else
    
       ! Not possible to run; flag it and return
       call mem_set_status(mem_mgr, 2)
       return
    
    end if
    
    ! Begin memory savings loop 2
    
    
    mcurr = 1
    curr_pickup = (/1,1,1/)
    next_pickup = (/1,1,1/)
    contrib_offset = 1
    
    mem_done = .FALSE.
    intra_pair = .FALSE.
    intra_done = .TRUE.
    
    do while (.NOT.(mem_done))
    
       ! Find best set of contraction matrices given memory setup:
       
       ! For now, use a relatively straightforward scheme: There are some opportunities to refine, but
       ! they will likely only result in moderate improvements and will be significantly more elaborate
       ! Summary: Fill workunit with all contractions for outer pairs until no more pairs can be fit
       ! If an outer pair does not fit into the workunit, then divide it into workunits and do
       ! all of them in turn, then move to next outer pair
       ! Currently no mixing of "full pairs" and divided workunits - there is some room for
       ! improvement here but it is more difficult to make
       
       curr_remain = mem_track
       
       ! If at start of one outer pair (i.e. not currently in a pair separated into workunits)
       if (.NOT.(intra_pair)) then
       
          wunit_size = (/0,0/)
       
          any_heavy = 0
   
          do i = curr_pickup(1), cache%num_outer
          
             ! Can the whole pair be done?
             if (curr_remain >= outer_contract_sizes_2(i, 1) + outer_contract_sizes_2(i, 2)) then
                          
                ! If yes, add to workload
                curr_remain = curr_remain - (outer_contract_sizes_2(i, 1) + &
                                             outer_contract_sizes_2(i, 2))
                               
                next_pickup = (/i+1, 1, 1/)
                
                wunit_size(1) = wunit_size(1) + outer_contract_sizes_2(i, 1)
                wunit_size(2) = wunit_size(2) + outer_contract_sizes_2(i, 2)
                
                ! If this was the last pair, mark as "done after this unit ends"
                if (i == cache%num_outer) then
                
                   mem_done = .TRUE.
                   
                end if
                                             
 
             ! If not 
             else
             
                ! If this was not the first pair in the work unit, don't add any 
                ! more to the unit and proceed
                if (next_pickup(1) - curr_pickup(1) > 0) then
                
                   exit
                
                ! If this was the first pair, make a workunit separation setup for this pair and 
                ! do the first such piece as this workunit
                else 
                
                   intra_pair = .TRUE.
                   intra_done = .FALSE.
                
                   ! Scheme: Divide half the available memory on the LHS contraction matrices
                   ! and the other half (+ 1 if odd # of matrices available) on the RHS,
                   ! unless LHS or RHS length is smaller than half; if so, do all of such
                   ! and give all remaining matrices to the other
                   
                   ! If enough for all RHS matrices
                   if (outer_contract_sizes_2(i,2) < (mem_track/2 + mod(mem_track,2))  ) then
                   
                   
                   
                      ! Flag to switch increment handling to "right-heavy" mode
                      any_heavy = 1
                   
                      wunit_maxsize = (/mem_track - outer_contract_sizes_2(i,2), &
                                     outer_contract_sizes_2(i,2)/)
                   
                   ! Otherwise, if enough for all LHS matrices
                   elseif (outer_contract_sizes_2(i,1) < mem_track/2) then
                   
                      ! Flag to switch increment handling to "left-heavy" mode
                      any_heavy = 2
                   
                      wunit_maxsize = (/outer_contract_sizes_2(i,1), &
                                     mem_track - outer_contract_sizes_2(i,1)/)
                                     
                   ! Otherwise, the default "half-and-half" scheme applies
                   else
                   
                      wunit_maxsize = (/mem_track/2, mem_track/2 + mod(mem_track,2)/)
                   
                   end if
                   
                   
                   wunit_size = wunit_maxsize
                   
                   ! Set up next pickup point
                   
                   ! If left- or right-heavy
                   if (any_heavy > 0) then
             
                      next_pickup(1) = curr_pickup(1)
                      next_pickup(4 - any_heavy) = 1
                      next_pickup(1 + any_heavy) = curr_pickup(1 + any_heavy) + wunit_maxsize(any_heavy)
          
                   ! Otherwise, do default handling
                   else
                         
                      next_pickup = (/curr_pickup(1), curr_pickup(2), curr_pickup(3) + wunit_maxsize(2)/)
             
                   end if
                   
                   exit
                
                end if
                
             end if
             
          end do
          
       ! If not at start of one outer pair (i.e. if currently in a pair separated into workunits)
       else
       
          wunit_size = wunit_maxsize
       
          ! If left- or right-heavy
          if (any_heavy > 0) then
          
             ! If this was the last unit of the pair:
             if (curr_pickup(any_heavy + 1) + wunit_maxsize(any_heavy) > &
                 outer_contract_sizes_2(curr_pickup(1), any_heavy)) then
             
                ! If this was the last unit of the last pair, mark as "done after this unit ends"
                if (curr_pickup(1) == cache%num_outer) then
                   mem_done = .TRUE.
                end if
                
                intra_done = .TRUE.
                
                ! Set the next pickup to be the start of the next outer pair
                next_pickup = (/curr_pickup(1) + 1, 1, 1/)
                
                wunit_size(any_heavy) = outer_contract_sizes_2(curr_pickup(1), any_heavy) - &
                                        curr_pickup(any_heavy + 1) + 1
                
             ! Otherwise, do only next unit of this pair
             else
             
                next_pickup(1) = curr_pickup(1)
                next_pickup(4 - any_heavy) = 1
                next_pickup(1 + any_heavy) = curr_pickup(1 + any_heavy) + wunit_maxsize(any_heavy)
                
             end if
          
          ! Otherwise, do default handling
          else
          
             ! If this is the last RHS unit of the pair:
             if (curr_pickup(3) + wunit_maxsize(2) > outer_contract_sizes_2(i,2)) then
             
                ! If also the last LHS unit of the pair, then it is the last overall unit of the pair
                if (curr_pickup(2) + wunit_maxsize(1) > outer_contract_sizes_2(i,1)) then
                
                   intra_done = .TRUE.
                
                   next_pickup = (/curr_pickup(1) + 1, 1, 1/)
             
                   ! If this was the last unit of the last pair, mark as "done after this unit ends"
                   if (curr_pickup(1) == cache%num_outer) then
                   
                      mem_done = .TRUE.
                      
                   end if
                   
                ! Otherwise, change to next LHS unit for next iteration
                else 
                
                   next_pickup = (/curr_pickup(1), curr_pickup(2) + wunit_maxsize(1), 1/)
                
                   ! If next LHS unit is the last LHS unit, update the sizes for next workunit
                   ! to make maximum use of memory
                   if (curr_pickup(2) + 2 * wunit_maxsize(1) > outer_contract_sizes_2(i,1)) then
                   
                      wunit_maxsize(1) = outer_contract_sizes_2(i,1) - &
                                      (curr_pickup(2) + wunit_maxsize(1)) + 1
                      wunit_maxsize(2) = min(mem_track - wunit_size(1), outer_contract_sizes_2(i,2))
                      
                   end if
          
                end if
                
                wunit_size(2) = outer_contract_sizes_2(i,2) - curr_pickup(3) + 1
             
             ! Otherwise, do next RHS unit
             else
             
                next_pickup = (/curr_pickup(1), curr_pickup(2), curr_pickup(3) + wunit_maxsize(2)/)
             
             end if
             
          end if   
          
       end if

       
       ! Allocate and set up perturbed density matrices for contractions
    
       call mem_incr(mem_mgr, wunit_size(1))
             
       if (.NOT. mem_mgr%calibrate) then
     
          allocate(LHS_dmat_2(wunit_size(1)))
       
          ! LHS perturbed D for 2-el terms
          do i = 1, size(LHS_dmat_2)
       
             call QcMatInit(LHS_dmat_2(i))
       
          end do
          
       end if
       
       call mem_incr(mem_mgr, wunit_size(2))
       
       if (.NOT. mem_mgr%calibrate) then
       
          allocate(RHS_dmat_2(wunit_size(2)))
       
          ! RHS perturbed D for 2-el terms
          do i = 1, size(RHS_dmat_2)
       
             call QcMatInit(RHS_dmat_2(i))
       
          end do
          
       end if
           
       
! Not sure if my ll to array changes work as intended here
! Next lines are prev version of code kept for reference for potential debugging
!        k = 1
!        
!        ! Cycle cache until at correct outer element
!        do i = 1, curr_pickup(1) - 1
!           outer_next => outer_next%next
!           k = k + 1
!        end do
       
       
! DEC 19: CONTINUE HERE
       
       
       lhs_ctr_2 = 1
       rhs_ctr_2 = 1

       if (.NOT.(intra_pair)) then

          ! FIXME: Loop assumes that there is a dummy entry in position 1 of array
          ! Traverse outer elements and fetch matrices
          
          do i = curr_pickup(1), next_pickup(1) - 1

             if (.NOT.(mem_mgr%calibrate)) then
          
                ! No chain rule applications
                if (cache%contribs_outer(i + 1)%num_dmat == 0) then
           
                   call contrib_cache_getdata_outer(size(D), D, 1, (/get_emptypert()/), .FALSE., &
                        1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_2(lhs_ctr_2))
                   call contrib_cache_getdata_outer(size(D), D, 1, (/get_emptypert()/), .FALSE., &
                        1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))
                        
                ! One chain rule application
                else if (cache%contribs_outer(i + 1)%num_dmat == 1) then
          
                   do m = 1, outer_contract_sizes_2(i, 1) 
                   
                      call contrib_cache_getdata_outer(size(D), D, 1, &
                           (/cache%contribs_outer(i + 1)%p_tuples(1)/), .FALSE., &
                           1, ind_len=size(cache%contribs_outer(i + 1)%indices, 2), &
                           ind_unsorted=cache%contribs_outer(i + 1)%indices(m, :), &
                           mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
                   end do
             
                   call contrib_cache_getdata_outer(size(D), D, 1, (/get_emptypert()/), .FALSE., &
                        1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))
          
                ! Two chain rule applications
                else if (cache%contribs_outer(i + 1)%num_dmat == 2) then
            
                   do m = 1, outer_contract_sizes_2(i, 1) 
              
                      call contrib_cache_getdata_outer(size(D), D, 1, &
                           (/cache%contribs_outer(i + 1)%p_tuples(1)/), .FALSE., &
                           1, ind_len=cache%contribs_outer(i + 1)%p_tuples(1)%npert, &
                           ind_unsorted=cache%contribs_outer(i + 1)%indices(1 + &
                           (m - 1) * cache%contribs_outer(i + 1)%blks_tuple_triang_size(2), &
                           1:cache%contribs_outer(i + 1)%p_tuples(1)%npert), &
                           mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
                   end do
             
                   do n = 1, outer_contract_sizes_2(i, 2) 
             
                      call contrib_cache_getdata_outer(size(D), D, 1, &
                      (/cache%contribs_outer(i + 1)%p_tuples(2)/), .FALSE., &
                           1, ind_len=cache%contribs_outer(i + 1)%p_tuples(2)%npert, &
                           ind_unsorted=cache%contribs_outer(i + 1)%indices(n, &
                           cache%contribs_outer(i + 1)%p_tuples(1)%npert + 1: &
                           cache%contribs_outer(i + 1)%p_tuples(1)%npert + &
                           cache%contribs_outer(i + 1)%p_tuples(2)%npert), &
                           mat_sing=RHS_dmat_2(rhs_ctr_2 + n  - 1))
                   end do          
          
                end if
      
                lhs_ctr_2 = lhs_ctr_2 + outer_contract_sizes_2(i, 1)
                rhs_ctr_2 = rhs_ctr_2 + outer_contract_sizes_2(i, 2)
          
             end if
        
          end do
       
          this_outer_size = dot_product( &
          outer_contract_sizes_2(curr_pickup(1):next_pickup(1) - 1, 1), &
          outer_contract_sizes_2(curr_pickup(1):next_pickup(1) - 1, 2) )
       
          if (.NOT.(mem_mgr%calibrate)) then
          
             write(out_str, *) 'Calculating second-order density matrix-dependent contribution'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
          
       
             ! Calculate two-electron contributions
             call get_2el_exp(num_pert, pert_ext, next_pickup(1) - curr_pickup(1), &
                  outer_contract_sizes_2(curr_pickup(1):next_pickup(1) - 1, 1), LHS_dmat_2, & 
                  outer_contract_sizes_2(curr_pickup(1):next_pickup(1) - 1, 2), RHS_dmat_2, &
                  cache%blks_triang_size*this_outer_size, &               
                  contrib_2(contrib_offset:contrib_offset + cache%blks_triang_size*this_outer_size - 1))
      
             write(out_str, *) 'Second-order density matrix-dependent contribution(sample)', &
             contrib_2(1:min(10,size(contrib_2)))
             call out_print(out_str, 2)
            
          end if
               
          contrib_offset = contrib_offset + this_outer_size

       ! FIXME: I (MaR) don't know if this else condition has ever been tested or encountered
       ! I don't see clearly how this would have worked with the ll scheme and I therefore
       ! don't know how to adapt it for the array cache scheme
       ! I will therefore comment it out and return to it later if relevant
       else
!           
!           if (.NOT.(mem_mgr%calibrate)) then
!           
!              ! No chain rule applications
!              if (outer_next%num_dmat == 0) then
!               
!                 call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
!                      1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_2(1))
!                 call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
!                      1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(1))
!       
!              ! One chain rule application
!              else if (outer_next%num_dmat == 1) then
!              
!                 dctr = 1
!                 
!                 do m = curr_pickup(2), curr_pickup(2) + wunit_size(1) - 1
!                 
!                    call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
!                         1, ind_len=size(outer_next%indices, 2), ind_unsorted=outer_next%indices(m, :), &
!                         mat_sing=LHS_dmat_2(dctr))
!                         
!                    dctr = dctr + 1
!                         
!                 end do
!                 
!                 call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
!                      1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(1))
!              
!              ! Two chain rule applications
!              else if (outer_next%num_dmat == 2) then
!          
!                 dctr = 1
!             
!                 do m = 1, curr_pickup(2), curr_pickup(2) + wunit_size(1) - 1
!                 
!                    call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
!                         1, ind_len=outer_next%p_tuples(1)%npert, &
!                         ind_unsorted=outer_next%indices(1 + &
!                         (m - 1) * outer_next%blks_tuple_triang_size(2), &
!                         1:outer_next%p_tuples(1)%npert), &
!                         mat_sing=LHS_dmat_2(dctr))
!                         
!                    dctr = dctr + 1
!                    
!                 end do
!                 
!                 dctr = 1
!                 
!                 do n = 1, curr_pickup(3), curr_pickup(3) + wunit_size(2) - 1
!                 
!                    call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(2)/), .FALSE., &
!                         1, ind_len=outer_next%p_tuples(2)%npert, &
!                         ind_unsorted=outer_next%indices(n, outer_next%p_tuples(1)%npert + 1: &
!                         outer_next%p_tuples(1)%npert + outer_next%p_tuples(2)%npert), &
!                         mat_sing=RHS_dmat_2(dctr))
!                         
!                    dctr = dctr + 1
!                    
!                 end do
!              
!              end if
!              
!           end if
! 
!           ! Make temporary contrib_2 array to hold return values
!           allocate(contrib_2_tmp(cache%blks_triang_size * wunit_size(1) * wunit_size(2)))
!           
!           if (.NOT.(mem_mgr%calibrate)) then
!           
!              write(out_str, *) 'Calculating second-order density matrix-dependent contribution'
!              call out_print(out_str, 1)
!              write(out_str, *) ' '
!              call out_print(out_str, 1)
!           
!              ! Calculate two-electron contributions
!              call get_2el_exp(num_pert, pert_ext, 1, &
!                   (/wunit_size(1)/), LHS_dmat_2, & 
!                   (/wunit_size(2)/), RHS_dmat_2, &
!                   cache%blks_triang_size * wunit_size(1) * wunit_size(2), &               
!                   contrib_2_tmp)
!                   
!              write(out_str, *) 'Second-order density matrix-dependent contribution (sample)', &
!              contrib_2_tmp(1:min(10,size(contrib_2_tmp)))
!              call out_print(out_str, 2)
!                   
!           
!              ! Put temporary array data into contrib_2 in the appropriate places
!           
!              dctr = 0
!           
!              do m = curr_pickup(2), curr_pickup(2) + wunit_size(1) - 1
!              
!                 do n = curr_pickup(3), curr_pickup(3) + wunit_size(2) - 1
!              
!                    contrib_2(contrib_offset +  &
!                              (m - 1) * outer_contract_sizes_2(k, 2) * cache%blks_triang_size + &
!                              (n - 1) * cache%blks_triang_size : &
!                              contrib_offset + &
!                              (m - 1) * outer_contract_sizes_2(k, 2) * cache%blks_triang_size + &
!                              (n) * cache%blks_triang_size - 1) = &
!                    contrib_2_tmp(cache%blks_triang_size * dctr + 1: &
!                                  cache%blks_triang_size * (dctr + 1))
!              
!                    dctr = dctr + 1
!              
!                 end do
!              
!              end do
!           
!           end if
!           
!           ! Increment contribution offset if this was the last workunit of this outer pair
!           if (next_pickup(1) - curr_pickup(1) > 0) then
!           
!              contrib_offset = contrib_offset + cache%blks_triang_size * &
!                               outer_contract_sizes_2(k, 1) * outer_contract_sizes_2(k, 2) 
!           
!           end if
!           
!           deallocate(contrib_2_tmp)
          
          
          
       end if
       
       
       
       ! Deallocate LHS and RHS contraction matrices
    
       call mem_decr(mem_mgr, wunit_size(1))
             
       if (.NOT. mem_mgr%calibrate) then
     
          do i = 1, size(LHS_dmat_2)
             call QcMatDst(LHS_dmat_2(i))
          end do
          
          deallocate(LHS_dmat_2)
          
       end if
       
       call mem_decr(mem_mgr, wunit_size(2))
       
       if (.NOT. mem_mgr%calibrate) then
       
          do i = 1, size(RHS_dmat_2)
             call QcMatDst(RHS_dmat_2(i))
          end do
          
          deallocate(RHS_dmat_2)
          
       end if
       
       ! Update pickup point for next iteration
       curr_pickup = next_pickup
    
       if (intra_done) then
          intra_pair = .FALSE.
          
       end if
    
    end do
    
    
    
    
    
    ! End memory savings loop 2
    
    
    
    
    
    
    ! Traversal: Add nuc-nuc, 1-el and two-el contributions together
    
    if (.NOT.(mem_mgr%calibrate)) then
     
       k = 1
     
       c1_ctr = 1
       c2_ctr = 1
       
            ! Traverse and fetch matrices
       do m = 1, size(cache%contribs_outer)
    
          if (cache%contribs_outer(m)%dummy_entry) then
       
             cycle
       
          end if
       
          ! Nuc-nuc, one-el and two-el contribution
          if (cache%contribs_outer(m)%num_dmat == 0) then
          
             allocate(cache%contribs_outer(m)%data_scal(cache%blks_triang_size))
             allocate(data_tmp(cache%blks_triang_size))
       
             ! Factor 0.5 for two-el because no chain rule applications
             
             data_tmp = contrib_0
             data_tmp = data_tmp + &
             contrib_1(c1_ctr:c1_ctr + cache%blks_triang_size - 1)
             data_tmp = data_tmp + &
             0.5*contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size - 1)
     
             c1_ctr = c1_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
             c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
     
          ! One-el and two-el contribution
          else if (cache%contribs_outer(m)%num_dmat == 1) then
          
             allocate(cache%contribs_outer(m)%data_scal(cache%blks_triang_size * &
                      outer_contract_sizes_2(k, 1) * &
                      outer_contract_sizes_2(k, 2)))
             allocate(data_tmp(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                      outer_contract_sizes_2(k, 2)))
                      
    
             data_tmp = contrib_1(c1_ctr:c1_ctr + cache%blks_triang_size * &
                               outer_contract_sizes_2(k, 1) - 1)
             
             data_tmp = data_tmp + &
             contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size * &
                               outer_contract_sizes_2(k, 1) - 1)
    
             c1_ctr = c1_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
             c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
                      
          ! Only two-electron contribution
          else if (cache%contribs_outer(m)%num_dmat == 2) then
          
             allocate(cache%contribs_outer(m)%data_scal(cache%blks_triang_size * &
                      outer_contract_sizes_2(k, 1) * &
                      outer_contract_sizes_2(k, 2)))
             allocate(data_tmp(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                      outer_contract_sizes_2(k, 2)))
                      
             data_tmp = contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size * &
                               outer_contract_sizes_2(k, 1) * outer_contract_sizes_2(k, 2) - 1)
          
             c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                      outer_contract_sizes_2(k, 2)
          end if
          
          if (cache%contribs_outer(m)%num_dmat == 0) then
          
          
             cache%contribs_outer(m)%data_scal = data_tmp
             
          else
       
             ! Set up collective block information for indexing
          
             if (cache%p_inner%npert > 0) then
             
                tot_num_pert = cache%p_inner%npert + &
                sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1,  &
                cache%contribs_outer(m)%num_dmat)/))
                      
                allocate(blks_tuple_info(cache%contribs_outer(m)%num_dmat + 1,tot_num_pert, 3))
                allocate(blk_sizes(cache%contribs_outer(m)%num_dmat + 1, tot_num_pert))
             
                blks_tuple_info = 0
                blk_sizes = 0
                   
                do j = 1, cache%contribs_outer(m)%num_dmat + 1
                   
                   if (j == 1) then
                   
                      do n = 1, cache%nblks
                      
                         blks_tuple_info(j, n, :) = cache%blk_info(n, :)
                         
                      end do
                      
                      blk_sizes(j, 1:cache%nblks) = cache%blk_sizes
                   
                   
                   else
                   
                      do n = 1, cache%contribs_outer(m)%nblks_tuple(j - 1)
                   
                         do p = 1, 3
                      
                            blks_tuple_info(j, n, :) = &
                            cache%contribs_outer(m)%blks_tuple_info(j - 1, n, :)
                      
                         end do
                   
                      end do
                   
                      blk_sizes(j, 1:cache%contribs_outer(m)%nblks_tuple(j-1)) = &
                      cache%contribs_outer(m)%blk_sizes(j-1, &
                      1:cache%contribs_outer(m)%nblks_tuple(j-1))
                      
                   end if
                   
                end do
          
                ! Go through elements of data_tmp and store in appropriate cache position
          
                do i = 1, size(cache%contribs_outer(m)%indices, 1)
             
                   do j = 1, size(cache%indices, 1)
      
                      offset = get_triang_blks_tuple_offset(cache%contribs_outer(m)%num_dmat + 1, &
                      cache%p_inner%npert + sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, &
                      n = 1, cache%contribs_outer(m)%num_dmat)/)), &
                      (/cache%nblks, (/(cache%contribs_outer(m)%nblks_tuple(n), n = 1, &
                      cache%contribs_outer(m)%num_dmat) /) /), &
                      (/cache%p_inner%npert, (/(cache%contribs_outer(m)%p_tuples(n)%npert, &
                      n = 1, cache%contribs_outer(m)%num_dmat)/)/), &
                      blks_tuple_info, &
                      blk_sizes, &
                      (/cache%blks_triang_size, &
                      (/(cache%contribs_outer(m)%blks_tuple_triang_size(n), n = 1, &
                      cache%contribs_outer(m)%num_dmat)/)/), &
                      (/cache%indices(j, :), cache%contribs_outer(m)%indices(i, :)/))
                   
!                       write(*,*) 'Index tuple:', (/cache%indices(j, :), cache%contribs_outer(m)%indices(i, :)/)
!                       write(*,*) 'Saving element', j + size(cache%indices, 1) * (i - 1), &
!                       'of data in cache element', offset
!                       write(*,*) 'data is',  data_tmp(j + size(cache%indices, 1) * (i - 1))
                   
                      cache%contribs_outer(m)%data_scal(offset) = data_tmp(j + &
                      size(cache%indices, 1) * (i - 1))
                      
                   end do
             
                end do
                
                deallocate(blk_sizes)
                deallocate(blks_tuple_info)
             
             else
             
                write(out_str, *) 'ERROR: Inner perturbation tuple is empty'
                call out_print(out_str, -1)
                       
             end if
          
          
          end if
          
          deallocate(data_tmp)
      
          k = k + 1
       
       end do
     
    end if   
    
    deallocate(contrib_0)
    deallocate(contrib_1)
    deallocate(contrib_2)
    
    
    deallocate(outer_contract_sizes_1_coll)
    deallocate(outer_contract_sizes_1)
    deallocate(outer_contract_sizes_2)
    
  end subroutine

   ! Recurse to identify (dryrun == .TRUE.) or assemble (dryrun == .FALSE.) two-factor contributions
   recursive subroutine rsp_twofact_recurse(pert, kn, p12, out_print, dryrun, len_cache, cache, p_size, prop)

    implicit none

    logical :: dryrun, lag_eligible
    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    integer :: len_cache
    integer, dimension(2) :: cache_loc
    type(contrib_cache), allocatable, dimension(:) :: cache
    integer ::  i, j
    integer :: hard_offset
    integer, dimension(2) :: kn
    integer :: nblks, block_size, p_size
    integer, allocatable, dimension(:,:,:) :: blk_info
    complex(8), dimension(p_size) :: prop
    
    external :: out_print
    character(len=2047) :: out_str
    
    ! Recurse
    if (pert%npert > 0) then

       call rsp_twofact_recurse(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), out_print, dryrun, &
       len_cache, cache, p_size, prop)
       
       call rsp_twofact_recurse(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), out_print, dryrun, &
       len_cache, cache, p_size, prop)

    ! If at end of recursion, process
    else
    
       p12(1) = p_tuple_standardorder(p12(1))
       p12(2) = p_tuple_standardorder(p12(2))
       
       ! Get size of p12(2) contribution for cache addressing offset
    
       nblks = get_num_blks(p12(1))
       allocate(blk_info(1, nblks, 3))
       blk_info(1,:,:) = get_blk_info(nblks, p12(1))
       block_size = get_triangulated_size(nblks, blk_info)
       deallocate(blk_info)
       
       if (p12(2)%npert > 0) then
       
          nblks = get_num_blks(p12(2))
          allocate(blk_info(1, nblks, 3))
          blk_info(1,:,:) = get_blk_info(nblks, p12(2))
          block_size = block_size * get_triangulated_size(nblks, blk_info)
          deallocate(blk_info)
       
       end if
       
       hard_offset = 0

       ! If not skipping
       if (kn_skip(p12(2)%npert, p12(2)%pid, kn) .EQV. .FALSE.) then

          if (contrib_cache_already(len_cache, cache, 2, p12, n_rule=kn(2))) then

             if (.NOT.(dryrun)) then
             
                write(out_str, *) 'Retrieving Pulay n contribution:'
                call out_print(out_str, 2)
                write(out_str, *) 'S', p12(1)%pid
                call out_print(out_str, 2)
                write(out_str, *) 'W', p12(2)%pid
                call out_print(out_str, 2)
                write(out_str, *) ' '
                call out_print(out_str, 2)
                
                call contrib_cache_getdata(len_cache, cache, 2, p12, p_size, 0, scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
             else
             
                ! If previously identified as Lagrangian type contribution, change flag to
                ! also include Pulay n contribution
                cache_loc = contrib_cache_locate(size(cache), cache, 2, p12, n_rule=kn(2))
                
                if (cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type == 3) then
                
                   cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type = 4
                
                end if
             
             end if
       
          else
          
             if (dryrun) then
             
                write(out_str, *) 'Adding newly identified Pulay n contribution to cache:'
                call out_print(out_str, 2)
                write(out_str, *) 'S', p12(1)%pid
                call out_print(out_str, 2)
                write(out_str, *) 'W', p12(2)%pid
                call out_print(out_str, 2)
                write(out_str, *) ' '
                call out_print(out_str, 2)

                len_cache = size(cache)
                call contrib_cache_add_element(len_cache, cache, 2, p12, n_rule=kn(2))
                
                cache_loc = contrib_cache_locate(size(cache), cache, 2, p12, n_rule=kn(2))
                ! Flag contribution as Pulay n type
                cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type = 1
                
             
             else
             
                write(out_str, *) 'ERROR: Expected to find Pulay contribution but it was not present'
                call out_print(out_str, 0)
                write(out_str, *) 'S', p12(1)%pid
                call out_print(out_str, 0)
                write(out_str, *) 'W', p12(2)%pid
                call out_print(out_str, 0)
                write(out_str, *) ' '
                call out_print(out_str, -1)

             end if

          end if

       end if
       
       ! Contribution must include perturbation "a" to be eligible to be Lagrangian type
       lag_eligible = .FALSE.
       
       do j = 1, p12(1)%npert
       
          if (p12(1)%pid(j) == 1) then
          
             lag_eligible = .TRUE.
          
          end if
       
       end do
       
       ! If not skipping
       if ((kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) .AND. &
           (p12(1)%npert > 0) .AND. lag_eligible) then
           
          len_cache = size(cache) 

          if (contrib_cache_already(len_cache, cache, 2, p12, n_rule=kn(2))) then

             if (.NOT.(dryrun)) then
             
                write(out_str, *) 'Retrieving Lagrange contributions:'
                call out_print(out_str, 2)
                write(out_str, *) 'A', p12(1)%pid
                call out_print(out_str, 2)
                write(out_str, *) 'B', p12(2)%pid
                call out_print(out_str, 2)
                write(out_str, *) ' '
                call out_print(out_str, 2)
                
                len_cache = size(cache)
             
                ! Pulay Lagrange contribution
                call contrib_cache_getdata(len_cache, cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'Prop after Pulay Lagrange retrieval', prop(1:min(12, size(prop)))
                
!                 Idempotency Lagrange contribution
                call contrib_cache_getdata(len_cache, cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'Prop after idempotency Lagrange retrieval', prop(1:min(12, size(prop)))

!                 SCFE Lagrange contribution
                call contrib_cache_getdata(len_cache, cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                
!                 write(*,*) 'Prop after SCFE Lagrange retrieval', prop(1:min(12, size(prop)))
            
             else
             
                ! If previously identified as Pulay n type contribution, change flag to
                ! also include Lagrangian type contribution
                cache_loc = contrib_cache_locate(size(cache), cache, 2, p12, n_rule=kn(2))
                
                if (cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type == 1) then
                
                   cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type = 4
                
                end if
                
             end if
       
          else
          
             if (dryrun) then
             
                write(out_str, *) 'Adding newly identified Lagrange contributions to cache:'
                call out_print(out_str, 2)
                write(out_str, *) 'A', p12(1)%pid
                call out_print(out_str, 2)
                write(out_str, *) 'B', p12(2)%pid
                call out_print(out_str, 2)
                write(out_str, *) ' '
                call out_print(out_str, 2)
                
                len_cache = size(cache)
                call contrib_cache_add_element(len_cache, cache, 2, p12, n_rule=kn(2))
                
                ! Flag contribution as Lagrangian type
                cache_loc = contrib_cache_locate(size(cache), cache, 2, p12, n_rule=kn(2))
                cache(cache_loc(1))%contribs_outer(cache_loc(2))%contrib_type = 3
                
             else
              
                write(out_str, *) 'ERROR: Expected to find Lagrange contributions but they were not present'
                call out_print(out_str, 0)
                write(out_str, *) 'A', p12(1)%pid
                call out_print(out_str, 0)
                write(out_str, *) 'B', p12(2)%pid
                call out_print(out_str, 0)
                write(out_str, *) ' '
                call out_print(out_str, -1)
             
             end if
          
          end if

       end if 

    end if

  end subroutine
  
  ! Calculate two-factor contributions (Pulay n and Pulay, idempotency and SCFE Lagrangian)
  subroutine rsp_twofact_calculate(size_s, S, size_d, D, size_f, F, &
             get_ovl_exp, out_print, cache, mem_mgr)

    implicit none
    
    type(mem_manager) :: mem_mgr
    logical :: mem_done
    integer :: mctr, mcurr, msize, mem_track

    logical :: any_lagrange, select_terms
    integer :: cache_offset, i, j, k, m, n, p, c_ctr, c_snap, lagrange_max_n
    integer :: i_supsize, o_triang_size, offset, tot_num_pert, max_outer_npert
    integer :: o_ctr, size_lagrange, size_pulay_n
    integer :: sstr_incr
    integer, dimension(0) :: nof
    integer, allocatable, dimension(:) :: o_supsize, o_supsize_prime, o_size
    integer, allocatable, dimension(:) :: which_index_is_pid
    integer, allocatable, dimension(:,:) :: blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
    complex(8), allocatable, dimension(:) :: contrib_pulay
    character(30) :: mat_str, fmt_str
    type(p_tuple) :: p_inner
    type(p_tuple), allocatable, dimension(:,:) :: d_struct_inner, d_struct_o, d_struct_o_prime
    
    integer :: size_s
    type(contrib_cache_outer), dimension(size_s) :: S
    integer :: size_d
    type(contrib_cache_outer), dimension(size_s) :: D
    integer :: size_F
    type(contrib_cache_outer), dimension(size_s) :: F
    
    type(contrib_cache) :: cache

    type(QcMat), allocatable, dimension(:) :: Lambda, Zeta, W
    type(QcMat) :: Y, Z, D_unp
           
    integer, allocatable, dimension(:) :: pert_ext
    
    external :: get_ovl_exp
    
    external :: out_print
    character(len=2047) :: out_str

    ! To avoid maybe-uninitialized warning
    c_snap = 0
    offset = 0
    
    max_outer_npert = 0
    
    write(out_str, *) 'Number of outer contributions for this inner tuple:', cache%num_outer
    call out_print(out_str, 1)
    
    
    ! Getting unperturbed D for template use

    call mem_incr(mem_mgr, 1)
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QCMatInit(D_unp)

       call contrib_cache_getdata_outer(size_d, D, 1, (/get_emptypert()/), .FALSE., &
            contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)
            
    end if
    
    ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    ! Unsure about npert argument
    ! Update: Seems OK
    call p_tuple_to_external_tuple(cache%p_inner, cache%p_inner%npert, pert_ext)
    
    i_supsize = 0

    ! Size arrays for derivative superstructures
    allocate(o_supsize(cache%num_outer))
    allocate(o_supsize_prime(cache%num_outer))
    allocate(o_size(cache%num_outer))

    ! Prepare matrices for inner tuple
    
    i_supsize = derivative_superstructure_getsize(p_tuple_remove_first(cache%p_inner), &
                (/cache%p_inner%npert, cache%p_inner%npert/), .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(d_struct_inner(i_supsize, 3))

    sstr_incr = 0
    
    call derivative_superstructure(p_tuple_remove_first(cache%p_inner), &
          (/cache%p_inner%npert, cache%p_inner%npert/), .FALSE., &
          (/get_emptypert(), get_emptypert(), get_emptypert()/), &
          i_supsize, sstr_incr, d_struct_inner)  
          
    sstr_incr = 0

    
    ! Traversal: Determine number of matrices/terms for contraction
    
    any_lagrange = .FALSE.
    lagrange_max_n = 0
    
    size_pulay_n = 0
    size_lagrange = 0
    
    k = 1
    
    do m = 1, size(cache%contribs_outer)
    
       if (cache%contribs_outer(m)%dummy_entry) then
       
          cycle
          
       end if    
    
       if (cache%contribs_outer(m)%p_tuples(1)%npert > max_outer_npert) then
       
          max_outer_npert = cache%contribs_outer(m)%p_tuples(1)%npert
       
       end if
  
       if (cache%contribs_outer(m)%p_tuples(1)%npert == 0) then
       
          o_triang_size = 1
       
       else
       
          o_triang_size = cache%contribs_outer(m)%blks_tuple_triang_size(1)
       
       end if
  
  
       write(out_str, *) 'Outer contribution, type', cache%contribs_outer(m)%contrib_type
       call out_print(out_str, 1)

       do i = 1, cache%contribs_outer(m)%num_dmat
          
          write(out_str, *) 'B', cache%contribs_outer(m)%p_tuples(i)%plab
          call out_print(out_str, 1)
          
       end do
       
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       ! Here and elsewhere: k in kn rule does not matter as long as it 
       ! is > n; it will always > n so used like this
       
       ! Get derivative superstructure sizes
       
       o_supsize(k) = derivative_superstructure_getsize( &
                   cache%contribs_outer(m)%p_tuples(1), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/))
                   
       ! FIXME: Change .FALSE. to .TRUE. below since this is a primed term?
       o_supsize_prime(k) = derivative_superstructure_getsize( &
                   cache%contribs_outer(m)%p_tuples(1), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
       
       o_size(k) = cache%contribs_outer(m)%contrib_type * o_triang_size
       
       ! If contribution includes Pulay n terms, increase size accordingly
       if ((cache%contribs_outer(m)%contrib_type == 1) .OR. &
           (cache%contribs_outer(m)%contrib_type == 4)) then
       
          size_pulay_n = size_pulay_n + o_triang_size
              
       end if
       
       ! If contribution includes Lagrangian type terms, increase size accordingly
       if (cache%contribs_outer(m)%contrib_type >= 3) then
       
          size_lagrange = size_lagrange + 3 * o_triang_size
          
          any_lagrange = .TRUE.
          
          ! Possibly update max n rule for bound in matrix calls
          if (cache%contribs_outer(m)%n_rule > lagrange_max_n) then
          
             lagrange_max_n = cache%contribs_outer(m)%n_rule            
          
          end if
       
       end if
       
       ! Allocate array for answers
       ! The data is ordered consecutively in this order (and entries are omitted if not present):
       ! 1: Pulay n contributions
       ! 2: Pulay Lagrange contributions
       ! 3: SCFE Lagrange contributions
       ! 4: Idempotency Lagrange contributions
      
       allocate(cache%contribs_outer(m)%data_scal(o_size(k) * cache%blks_triang_size))
       
       cache%contribs_outer(m)%data_scal = 0.0
         

       k = k + 1
    
    end do
    
    allocate(which_index_is_pid(cache%p_inner%npert + max_outer_npert))
    
    
    ! Memory management only necessary from here

    if (mem_enough(mem_mgr, size_pulay_n + size_lagrange/3 + 1)) then
    
       ! Possible to run in non-savings mode
       mem_track = size_pulay_n + size_lagrange/3
       
    elseif (mem_enough(mem_mgr, 5)) then
    
       ! Possible to run in savings mode
       call mem_set_status(mem_mgr, 1)
       
       mem_track = mem_mgr%remain - 4
    
    else
    
       ! Not possible to run; flag it and return
       call mem_set_status(mem_mgr, 2)
       return       
    
    end if
    
    ! Pulay size is the number of Pulay n contributions + one third of the
    ! Lagrange type contributions (of which Pulay Lagrange is one of three)
    allocate(contrib_pulay(cache%blks_triang_size * (size_pulay_n + size_lagrange/3)))
    contrib_pulay = 0.0

    
    ! Begin memory management loop 1
    
    mcurr = 1
    
    mem_done = .FALSE.
    
    do while (.NOT.(mem_done))
    
       msize = min(mem_track, size_pulay_n + size_lagrange/3 - mcurr + 1)
       
       ! Flag if this is the last memory iteration
       if (msize + mcurr > size_pulay_n + size_lagrange/3) then
       
          mem_done = .TRUE.
       
       end if
    
       call mem_incr(mem_mgr, msize)
    
       if (.NOT.(mem_mgr%calibrate)) then
     
          allocate(W(msize))
       
          do i = 1, size(W)
       
             call QcMatInit(W(i), D_unp)
             call QcMatZero(W(i))
          
          end do
       
       end if
       
       ! Traversal: Make W matrices and store
       
       o_ctr = 1
       mctr = 0
       
       k = 1
       
       ! Traverse to get W matrices
       do m = 1, size(cache%contribs_outer)
    
          if (cache%contribs_outer(m)%dummy_entry) then
       
             cycle
          
          end if
     
             allocate(d_struct_o(o_supsize(k), 3))
             allocate(d_struct_o_prime(o_supsize_prime(k), 3))
             
             which_index_is_pid = 0
             
             do j = 1, cache%contribs_outer(m)%p_tuples(1)%npert
             
                which_index_is_pid(cache%contribs_outer(m)%p_tuples(1)%pid(j)) = j
             
             end do
             
             sstr_incr = 0
   
             ! Get derivative superstructures according to whether contribution
             ! is n or Lagrange
             if(cache%contribs_outer(m)%p_tuples(1)%npert ==0) then
             
                if (cache%contribs_outer(m)%contrib_type == 1 .OR. &
                    cache%contribs_outer(m)%contrib_type == 4) then
                       
                   call derivative_superstructure(get_emptypert(), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize(k), sstr_incr, d_struct_o)
                   
                   sstr_incr = 0
                
                end if
               
                if (cache%contribs_outer(m)%contrib_type == 3 .OR. &
                    cache%contribs_outer(m)%contrib_type == 4) then
                
                   call derivative_superstructure(get_emptypert(), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .TRUE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize_prime(k), sstr_incr, d_struct_o)
                      
                   sstr_incr = 0
      
                end if
                
                o_triang_size = 1
             
             else
         
                if (cache%contribs_outer(m)%contrib_type == 1 .OR. &
                    cache%contribs_outer(m)%contrib_type == 4) then             
         
                   call derivative_superstructure(cache%contribs_outer(m)%p_tuples(1), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize(k), sstr_incr, d_struct_o)
                   
                   sstr_incr = 0
                
                end if
      
               if (cache%contribs_outer(m)%contrib_type == 3 .OR. &
                   cache%contribs_outer(m)%contrib_type == 4) then
                
                   call derivative_superstructure(cache%contribs_outer(m)%p_tuples(1), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .TRUE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize_prime(k), sstr_incr, d_struct_o_prime)
                   
                   sstr_incr = 0
                   
                end if
     
                o_triang_size = size(cache%contribs_outer(m)%indices, 1)
                
             end if
     
             ! Calculate matrices
             do j = 1, o_triang_size
             
                if ((o_ctr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
                
                   select case (cache%contribs_outer(m)%contrib_type)
                
                   ! Only Pulay n
                   case (1)
                
                      call mem_incr(mem_mgr, 4)
                
                      if (.NOT.(mem_mgr%calibrate)) then
                
                         call rsp_get_matrix_w(o_supsize(k), d_struct_o, cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert, &
                              which_index_is_pid(1:cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert), &
                              size(cache%contribs_outer(m)%indices(j,:)), &
                              cache%contribs_outer(m)%indices(j,:), &
                              size(F), F, size(D), D, size(S), S, W(mctr + 1))
                      
                      end if
                           
                      call mem_decr(mem_mgr, 4)
                
                   ! Only Pulay Lagrange
                   case (3)
                
                      call mem_incr(mem_mgr, 4)
                
                      if (.NOT.(mem_mgr%calibrate)) then                
                
                         call rsp_get_matrix_w(o_supsize_prime(k), d_struct_o_prime, &
                              cache%p_inner%npert + cache%contribs_outer(m)%p_tuples(1)%npert, &
                              which_index_is_pid(1:cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert), & 
                              size(cache%contribs_outer(m)%indices(j,:)), &
                              cache%contribs_outer(m)%indices(j,:), &
                              size(F), F, size(D), D, size(S), S, W(mctr + 1))
                           
                      end if
                           
                      call mem_decr(mem_mgr, 4)                           
                
                   ! Both Pulay n and Lagrange (doing Pulay n in first batch)
                   case (4)

                      call mem_incr(mem_mgr, 4)
                
                      if (.NOT.(mem_mgr%calibrate)) then                     
                   
                         call rsp_get_matrix_w(o_supsize(k), d_struct_o, cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert, &
                              which_index_is_pid(1:cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert), & 
                              size(cache%contribs_outer(m)%indices(j,:)), &
                              cache%contribs_outer(m)%indices(j,:), &
                              size(F), F, size(D), D, size(S), S, W(mctr + 1))
                           
                      end if
                           
                      call mem_decr(mem_mgr, 4)                       
                
                   end select
                
                   mctr = mctr + 1
                
                end if
                     
                o_ctr = o_ctr + 1
          
             end do
             
             ! If contribution type is 4, loop to get Pulay Lagrange matrices
             if (cache%contribs_outer(m)%contrib_type == 4) then
             
                ! Calculate matrices
                do j = 1, o_triang_size
             
                   if ((o_ctr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
                   
                      call mem_incr(mem_mgr, 4)
                
                      if (.NOT.(mem_mgr%calibrate)) then   
                         
                         call rsp_get_matrix_w(o_supsize_prime(k), d_struct_o_prime, &
                              cache%p_inner%npert + cache%contribs_outer(m)%p_tuples(1)%npert, &
                              which_index_is_pid(1:cache%p_inner%npert + &
                              cache%contribs_outer(m)%p_tuples(1)%npert), & 
                              size(cache%contribs_outer(m)%indices(j,:)), &
                              cache%contribs_outer(m)%indices(j,:), &
                              size(F), F, size(D), D, size(S), S, W(mctr + 1))
                              
                      end if
                           
                      call mem_decr(mem_mgr, 4)                                  
                  
                      mctr = mctr + 1
                
                   end if
                     
                   o_ctr = o_ctr + 1
           
                end do
             
             end if
             
          
             deallocate(d_struct_o)
             deallocate(d_struct_o_prime)
                
          k = k + 1
       
       end do
       
       ! Outside traversal: Calculate all Pulay contributions
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          call get_ovl_exp(0, nof, 0, nof, size(pert_ext), pert_ext, size(W), W, &
                           cache%blks_triang_size * msize, &
                           contrib_pulay(cache%blks_triang_size * (mcurr - 1) + 1 : &
                           cache%blks_triang_size * (mcurr - 1 + msize)))
                           
       end if

       if (.NOT.(mem_mgr%calibrate)) then
     
          do i = 1, msize
       
             call QcMatDst(W(i))
          
          end do
          
          deallocate(W)
       
       end if
       
       call mem_decr(mem_mgr, msize)
       
       mcurr = mcurr + msize
    
    end do
    ! End memory management loop 1
   
    contrib_pulay = -2.0 * contrib_pulay
    
    ! Traversal: Store Pulay contributions
    
    k = 1
    o_ctr = 0
    
    ! Traversal: Store Pulay contributions
    do m = 1, size(cache%contribs_outer)
    
       if (cache%contribs_outer(m)%dummy_entry) then
       
          cycle
          
       end if

       c_ctr = 0
    
       if (cache%contribs_outer(m)%p_tuples(1)%npert ==0) then
          
          o_triang_size = 1
          
       else
      
          o_triang_size = size(cache%contribs_outer(m)%indices, 1)
             
       end if
       
       ! Set up block information for indexing
       
       tot_num_pert = cache%p_inner%npert + &
       sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1, &
              cache%contribs_outer(m)%num_dmat)/))
                   
       allocate(blks_tuple_info(cache%contribs_outer(m)%num_dmat + 1,tot_num_pert, 3))
       allocate(blk_sizes(cache%contribs_outer(m)%num_dmat + 1, tot_num_pert))
               
       do j = 1, cache%contribs_outer(m)%num_dmat + 1
          
          if (j == 1) then
             
             do n = 1, cache%nblks
                
                blks_tuple_info(j, n, :) = cache%blk_info(n, :)
             
             end do
             
             blk_sizes(j, 1:cache%nblks) = cache%blk_sizes
          
          else
             
             if (size(cache%contribs_outer(m)%p_tuples) > 0) then
                
                if (cache%contribs_outer(m)%p_tuples(1)%npert > 0) then
                   
                   do n = 1, cache%contribs_outer(m)%nblks_tuple(j - 1)
                      
                      do p = 1, 3
                         
                         blks_tuple_info(j, n, :) = &
                         cache%contribs_outer(m)%blks_tuple_info(j - 1, n, :)
                      
                      end do
                   
                   end do
                   
                   blk_sizes(j, 1:cache%contribs_outer(m)%nblks_tuple(j-1)) = &
                   cache%contribs_outer(m)%blk_sizes(j-1, &
                   1:cache%contribs_outer(m)%nblks_tuple(j-1))
                
                end if
             
             end if
          
          end if
       
       end do
       
        
       ! Store any Pulay n and Lagrange terms
       
       do i = 1, o_triang_size
          
          do j = 1, size(cache%indices, 1)
               
             if (size(cache%contribs_outer(m)%p_tuples) > 0) then
             
                if (cache%contribs_outer(m)%p_tuples(1)%npert > 0) then

                   offset = get_triang_blks_tuple_offset(2, tot_num_pert, &
                   (/cache%nblks, cache%contribs_outer(m)%nblks_tuple(1)/), &
                   (/cache%p_inner%npert, cache%contribs_outer(m)%p_tuples(1)%npert/), &
                   blks_tuple_info, &
                   blk_sizes, &
                   (/cache%blks_triang_size, &
                   cache%contribs_outer(m)%blks_tuple_triang_size(1)/), &
                   (/cache%indices(j, :), cache%contribs_outer(m)%indices(i, :)/) )
                      
                else
                   
                   offset = j
                   
                end if
                   
             end if
                
             if (cache%contribs_outer(m)%contrib_type == 1) then
                
                cache%contribs_outer(m)%data_scal(offset) = contrib_pulay(j + &
                size(cache%indices, 1) * (i - 1) + o_ctr)
                c_ctr = c_ctr + 1
                   
             else if (cache%contribs_outer(m)%contrib_type == 3) then

                cache%contribs_outer(m)%data_scal(offset) = contrib_pulay(j + &
                size(cache%indices, 1) * (i - 1) + o_ctr)
                c_ctr = c_ctr + 1
                   
             else if (cache%contribs_outer(m)%contrib_type == 4) then
                
                cache%contribs_outer(m)%data_scal(offset) = contrib_pulay(j + &
                size(cache%indices, 1) * (i - 1) + o_ctr)
                   
                ! Skip one block to store Pulay Lagrange when both Pulay n and Lagrange
                cache%contribs_outer(m)%data_scal(offset + &
                size(cache%indices, 1) * o_triang_size) = &
                contrib_pulay(j + size(cache%indices, 1) * (i - 1) + &
                size(cache%indices, 1) * o_triang_size  + o_ctr)
                c_ctr = c_ctr + 2
                
                
             end if
                      
          end do
             
       end do

       o_ctr = o_ctr + size(cache%indices, 1) * o_triang_size
              
       deallocate(blk_sizes)
       deallocate(blks_tuple_info)
       
       k = k + 1
       
    end do
    
    deallocate(contrib_pulay)
    
    if (any_lagrange) then
      
       ! Memory management for idempotency and SCFE Lagrange contributions
       if (mem_enough(mem_mgr, 2*cache%blks_triang_size + 6)) then
      
          ! Possible to run in non-savings mode
          mem_track = cache%blks_triang_size
          
       elseif (mem_enough(mem_mgr, 4 + 4)) then
       
          ! Possible to run in savings mode
          call mem_set_status(mem_mgr, 1)
          
          mem_track = (mem_mgr%remain - 6)/2
      
       else
       
          ! Not possible to run; flag it and return
          call mem_set_status(mem_mgr, 2)
          return       
       
       end if
       
       ! Begin memory management loop 2
       
       mcurr = 1
       mem_done = .FALSE.
       
       
       do while (.NOT.(mem_done))
       
          msize = min(mem_track, cache%blks_triang_size - mcurr + 1)
          
          ! Flag if this is the last memory iteration
          if (msize + mcurr > cache%blks_triang_size) then
       
             mem_done = .TRUE.
       
          end if
       
          
!           if (msize < cache%blks_triang_size) then
!              
!              write(*,*) 'Memory loop for SCFE/idem.: Starting element:', mcurr , ', block size:', msize
!         
!           end if
       

          ! FIXME: Change to max order of all properties
          ! Update: May be OK anyway
          ! NOTE: Have another look at this part if results are wrong
      
          which_index_is_pid = 0
             
          do i = 1, cache%p_inner%npert
             
             which_index_is_pid(cache%p_inner%pid(i)) = i
             
          end do

          ! Get zeta and lambda matrices
          
          call mem_incr(mem_mgr, msize)
          call mem_incr(mem_mgr, msize)
          
          if (.NOT.(mem_mgr%calibrate)) then
          
             allocate(Lambda(msize))
             allocate(Zeta(msize))
          
             do i = 1, msize
      
                call QcMatInit(Lambda(i), D_unp)
                call QcMatInit(Zeta(i), D_unp)
                call QcMatZero(Lambda(i))
                call QcMatZero(Zeta(i))

                select_terms = find_residue_info(p_tuple_getone(cache%p_inner,1))
          
                call rsp_get_matrix_zeta(p_tuple_getone(cache%p_inner, 1), (/lagrange_max_n, &
                     lagrange_max_n/), i_supsize, d_struct_inner, maxval(cache%p_inner%pid), &
                     which_index_is_pid(1:maxval(cache%p_inner%pid)), &
                     size(cache%indices(mcurr + i - 1,:)), &
                     cache%indices(mcurr + i - 1,:), size(F), F, size(D), D, &
                     size(S), S, Zeta(i), &
                     select_terms_arg = select_terms)
                     
                call rsp_get_matrix_lambda(p_tuple_getone(cache%p_inner, 1), i_supsize, &
                     d_struct_inner, maxval(cache%p_inner%pid), &
                     which_index_is_pid(1:maxval(cache%p_inner%pid)), &
                     size(cache%indices(mcurr + i - 1,:)), cache%indices(mcurr + i - 1,:), &
                     size(D), D, size(S), S, Lambda(i), select_terms_arg = select_terms)
      
             end do
             
          end if
       
          
         
          
          call mem_incr(mem_mgr, 1)
          call mem_incr(mem_mgr, 1)
          
          if (.NOT.(mem_mgr%calibrate)) then
                   
             call QcMatInit(Z, D_unp)
             call QcMatZero(Z)
             call QcMatInit(Y, D_unp)
             call QcMatZero(Y)
          
          end if
          
          k = 1
          o_ctr = 0

          ! Traversal: Calculate/store idempotency/SCFE contributions
          do m = 1, size(cache%contribs_outer)
    
             if (cache%contribs_outer(m)%dummy_entry) then
       
                cycle
          
             end if
          
             c_ctr = 0
          
             if (cache%contribs_outer(m)%p_tuples(1)%npert ==0) then
                
                o_triang_size = 1
               
             else
            
               o_triang_size = size(cache%contribs_outer(m)%indices, 1)
                   
             end if
             
             ! Set up block information for indexing
             
             tot_num_pert = cache%p_inner%npert + &
             sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1, &
             cache%contribs_outer(m)%num_dmat)/))
                         
             allocate(blks_tuple_info(cache%contribs_outer(m)%num_dmat + 1,tot_num_pert, 3))
             allocate(blk_sizes(cache%contribs_outer(m)%num_dmat + 1, tot_num_pert))
                     
             do j = 1, cache%contribs_outer(m)%num_dmat + 1
                
                if (j == 1) then
                   
                   do n = 1, cache%nblks
                      
                      blks_tuple_info(j, n, :) = cache%blk_info(n, :)
                   
                   end do
                   
                   blk_sizes(j, 1:cache%nblks) = cache%blk_sizes
                
                else
               
                   if (size(cache%contribs_outer(m)%p_tuples) > 0) then
                      
                      if (cache%contribs_outer(m)%p_tuples(1)%npert > 0) then
                         
                         do n = 1, cache%contribs_outer(m)%nblks_tuple(j - 1)
                           
                            do p = 1, 3
                               
                              blks_tuple_info(j, n, :) = &
                              cache%contribs_outer(m)%blks_tuple_info(j - 1, n, :)
                            
                            end do
                         
                         end do
                         
                         blk_sizes(j, 1:cache%contribs_outer(m)%nblks_tuple(j-1)) = &
                         cache%contribs_outer(m)%blk_sizes(j-1, &
                         1:cache%contribs_outer(m)%nblks_tuple(j-1))
                      
                      end if
                  
                   end if
                
                end if
             
             end do
             
             ! Set up counters to store idempotency, SCFE terms in correct positions
             
             if ((cache%contribs_outer(m)%contrib_type == 1) .OR. &
             (cache%contribs_outer(m)%contrib_type == 3)) then
             
                c_snap = o_triang_size * size(cache%indices,1)
      
             else if (cache%contribs_outer(m)%contrib_type == 4) then
      
                c_snap = 2 * o_triang_size * size(cache%indices,1)
                      
             end if
             
             o_ctr = o_ctr + c_snap
            
             ! Calculate and store idempotency and SCFE terms
             if ((cache%contribs_outer(m)%contrib_type == 3) .OR. &
             (cache%contribs_outer(m)%contrib_type == 4)) then
           
                allocate(d_struct_o(o_supsize_prime(k), 3))
        
                
                which_index_is_pid = 0
                
                do j = 1, cache%contribs_outer(m)%p_tuples(1)%npert
                
                   which_index_is_pid(cache%contribs_outer(m)%p_tuples(1)%pid(j)) = j
                 
                end do
                
                sstr_incr = 0
                
                ! Get derivative superstructures
                if(cache%contribs_outer(m)%p_tuples(1)%npert ==0) then
                
                   call derivative_superstructure(get_emptypert(), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .TRUE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize_prime(k), sstr_incr, d_struct_o)
                   
                 else
            
                   call derivative_superstructure(cache%contribs_outer(m)%p_tuples(1), &
                   (/cache%contribs_outer(m)%n_rule, &
                   cache%contribs_outer(m)%n_rule/), .TRUE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                   o_supsize_prime(k), sstr_incr, d_struct_o)
                  
                end if
            
                ! Calculate Y, Z matrices and store
                do i = 1, o_triang_size
        
                   call mem_incr(mem_mgr, 4)
        
                   if (.NOT.(mem_mgr%calibrate)) then        
        
                      call QcMatZero(Y)
                      call rsp_get_matrix_y(o_supsize_prime(k), d_struct_o, cache%p_inner%npert + &
                           cache%contribs_outer(m)%p_tuples(1)%npert, &
                           which_index_is_pid(1:cache%p_inner%npert + &
                           cache%contribs_outer(m)%p_tuples(1)%npert), &
                           size(cache%contribs_outer(m)%indices(i,:)), &
                           cache%contribs_outer(m)%indices(i,:), &
                           size(F), F, size(D), D, size(S), S, Y, select_terms_arg = .FALSE.)
                           
                   end if
                   
                   call mem_decr(mem_mgr, 4)
                   
                   call mem_incr(mem_mgr, 4)
                           
                   if (.NOT.(mem_mgr%calibrate)) then
                   
                      ! NOTE: Rule choice very likely to give correct exclusion but
                      ! have another look if something goes wrong
                      call QcMatZero(Z)
                      call rsp_get_matrix_z(o_supsize_prime(k), d_struct_o, &
                           (/cache%contribs_outer(m)%n_rule, cache%contribs_outer(m)%n_rule/), &
                           cache%p_inner%npert + &
                           cache%contribs_outer(m)%p_tuples(1)%npert, &
                           which_index_is_pid(1:cache%p_inner%npert + &
                           cache%contribs_outer(m)%p_tuples(1)%npert), &
                           size(cache%contribs_outer(m)%indices(i,:)), &
                           cache%contribs_outer(m)%indices(i,:), &
                           size(F), F, size(D), D, size(S), S, Z, select_terms_arg = .FALSE.)
                           
                   end if
                  
                   call mem_decr(mem_mgr,4)
                                         
                   do j = mcurr, mcurr + msize - 1
                 
                      if (size(cache%contribs_outer(m)%p_tuples) > 0) then
                         if (cache%contribs_outer(m)%p_tuples(1)%npert > 0) then
        
                            offset = get_triang_blks_tuple_offset(2, tot_num_pert, &
                            (/cache%nblks, cache%contribs_outer(m)%nblks_tuple(1)/), &
                            (/cache%p_inner%npert, cache%contribs_outer(m)%p_tuples(1)%npert/), &
                            blks_tuple_info, &
                            blk_sizes, &
                            (/cache%blks_triang_size, &
                            cache%contribs_outer(m)%blks_tuple_triang_size(1)/), &
                            (/cache%indices(j, :), cache%contribs_outer(m)%indices(i, :)/))
                         
                         else
                      
                            offset = j
                      
                         end if
                      
                      end if

                      if (.NOT.(mem_mgr%calibrate)) then
                      
!                          write(*,*) 'Saving element', j, 'of data in', offset
                      
                         call QcMatTraceAB(Zeta(j), Z, &
                         cache%contribs_outer(m)%data_scal(c_snap + offset))
                         call QcMatTraceAB(Lambda(j), Y, &
                         cache%contribs_outer(m)%data_scal(c_snap + &
                         cache%blks_triang_size*o_triang_size + offset))
                   
                         cache%contribs_outer(m)%data_scal(c_snap + offset) = &
                         -2.0 * cache%contribs_outer(m)%data_scal(c_snap + offset)
                   
                         cache%contribs_outer(m)%data_scal(c_snap + &
                         cache%blks_triang_size*o_triang_size + offset) = &
                         -2.0 * cache%contribs_outer(m)%data_scal(c_snap + &
                         cache%blks_triang_size*o_triang_size + offset)
                         
                      end if
                
                   end do
                     
                end do
                
                deallocate(d_struct_o)
             
             end if
          
             deallocate(blk_sizes)
             deallocate(blks_tuple_info)
          
             k = k + 1
       
          end do
          
          if (.NOT.(mem_mgr%calibrate)) then
          
             do i = 1, msize
      
                call QcMatDst(Lambda(i))
                call QcMatDst(Zeta(i))
                
             end do
                
             deallocate(Lambda)
             deallocate(Zeta)
          
          end if
          
          call mem_decr(mem_mgr, msize)
          call mem_decr(mem_mgr, msize)
          
          if (.NOT.(mem_mgr%calibrate)) then
                   
             call QcMatDst(Z)
             call QcMatDst(Y)

          end if
          
          call mem_decr(mem_mgr, 1)
          call mem_decr(mem_mgr, 1)          
          
          mcurr = mcurr + msize
          
       end do
       ! End memory management loop 2    
    
    
    end if 
    
    
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QcMatDst(D_unp)
    
    end if
    
    call mem_decr(mem_mgr, 1)
    
    deallocate(which_index_is_pid)
    
    deallocate(d_struct_inner)
    
    deallocate(o_supsize)
    deallocate(o_supsize_prime)
    deallocate(o_size)
    
    
    
    
  end subroutine
  
  
  
  
  

  ! Print tensors recursively
  recursive subroutine print_rsp_tensor(npert, lvl, pdim, prop, offset)

    implicit none

    integer :: npert, i, j, offset, lvl, new_offset
    integer, dimension(npert) :: pdim
    complex(8), dimension(product(pdim)) :: prop

    if (lvl > 1) then

    do i = 1, pdim(npert - lvl + 1)

       new_offset = offset + (i - 1)*product(pdim(npert - lvl + 1:npert))/ &
                                             pdim(npert - lvl + 1)

       call print_rsp_tensor(npert, lvl - 1, pdim, prop, new_offset)

    end do

    open(unit=260, file='rsp_tensor', status='old', action='write', &
         position='append') 
    write(260,*) ' '
    close(260)

    else

    open(unit=260, file='rsp_tensor', status='old', action='write', &
         position='append') 
    write(260,*) real(prop(offset:offset+pdim(npert) - 1))
    close(260)

    end if

  end subroutine

  ! Print tensors recursively (non-redundant storage)
  recursive subroutine print_rsp_tensor_tr(lvl, npert, pdim, ind, &
                       nblks, blk_sizes, blk_info, propsize, prop, print_id, &
                       print_id_human)

    implicit none

    integer :: lvl, npert, nblks, i, propsize, print_id, print_id_human
    integer, dimension(npert) :: pdim
    integer, dimension(npert - 1) :: ind, new_ind
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    complex(8), dimension(pdim(npert)) :: line_for_print
    complex(8), dimension(propsize) :: prop
    character(20) :: format_line

    if (lvl == npert) then

       do i = 1, pdim(npert)

          line_for_print(i) = prop(get_triang_blks_offset(nblks, npert, &
                              blk_info, blk_sizes, (/ ind(:), i  /)))

       end do

       write(print_id, *) real(line_for_print)

       format_line = '(      f20.8)'
       write(format_line(2:7), '(i6)') pdim(npert)
       write(print_id_human, format_line) real(line_for_print)

    else

       new_ind = ind

       do i = 1, pdim(lvl)

          new_ind = ind
          new_ind(lvl) = i

          call print_rsp_tensor_tr(lvl + 1, npert, pdim, new_ind, &
                                   nblks, blk_sizes, blk_info, propsize, prop, print_id, print_id_human)

       end do

       write(print_id, *) ' '
       write(print_id_human, *) ' '

    end if

  end subroutine

  ! Print response tensor (non-redundant storage) to stdout
  recursive subroutine print_rsp_tensor_stdout_tr(lvl, npert, pdim, ind, &
                       nblks, blk_sizes, blk_info, propsize, prop)

    implicit none

    integer :: lvl, npert, nblks, i, propsize
    integer, dimension(npert) :: pdim
    integer, dimension(npert - 1) :: ind, new_ind
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    complex(8), dimension(pdim(npert)) :: line_for_print
    complex(8), dimension(propsize) :: prop

    if (lvl == npert) then

       do i = 1, pdim(npert)

          line_for_print(i) = prop(get_triang_blks_offset(nblks, npert, &
                                   blk_info, blk_sizes, (/ ind(:), i  /)))

       end do

       write(*,*) real(line_for_print)

    else

       new_ind = ind

       do i = 1, pdim(lvl)

          new_ind = ind
          new_ind(lvl) = i

          call print_rsp_tensor_stdout_tr(lvl + 1, npert, pdim, new_ind, &
                                          nblks, blk_sizes, blk_info, propsize, prop)

       end do

       write(*,*) ' '

    end if

  end subroutine

  ! Print response tensor to stdout
  recursive subroutine print_rsp_tensor_stdout(npert, lvl, pdim, prop, offset)

    implicit none

    integer :: npert, i, j, offset, lvl, new_offset
    integer, dimension(npert) :: pdim
    complex(8), dimension(product(pdim)) :: prop

    if (lvl > 1) then

       do i = 1, pdim(npert - lvl + 1)

          new_offset = offset + (i - 1)*product(pdim(npert - lvl + 1:npert))/ &
                                                pdim(npert - lvl + 1)

          call print_rsp_tensor_stdout(npert, lvl - 1, pdim, prop, new_offset)

       end do

       write(*,*) ' '

    else

       write(*,*) real(prop(offset:offset+pdim(npert) - 1))

    end if

  end subroutine


end module
