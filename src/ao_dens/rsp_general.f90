! Copyright 2015 Magnus Ringholm
!
!> @file Contains module rsp_general

!> General response routines. This module organizes, computes and prints
!> response function tensors.
module rsp_general

  use rsp_contribs, only: rsp_ovlave_t_matrix_2014
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
  use rsp_indices_and_addressing, only: get_blk_info,                 &
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
  use rsp_perturbed_sdf, only: rsp_fds
  use rsp_property_caching, only: contrib_cache_outer,                 &
                                  contrib_cache,                       &
                                  contrib_cache_initialize,            &
                                  contrib_cache_next_element,          &
                                  contrib_cache_outer_allocate,     &
                                  contrib_cache_outer_add_element,     &
                                  contrib_cache_outer_next_element,    &
                                  contrib_cache_outer_cycle_first,     &
                                  contrib_cache_cycle_outer,     &
                                  contrib_cache_add_element,           &
                                  contrib_cache_already,               &
                                  contrib_cache_getdata,               &
                                  contrib_cache_getdata_outer,               &
                                  contrib_cache_allocate, &
                                  contrib_cache_retrieve, &
                                  contrib_cache_outer_retrieve, &
                                  contrib_cache_outer_store, &
                                  contrib_cache_store, &
                                  mat_scal_store, &
                                  mat_scal_retrieve, &
                                  rs_check, &
                                  prog_incr, &
                                  prog_init
                                  
  
  use rsp_sdf_caching
  
  use qcmatrix_f

  implicit none

  public openrsp_get_property
  public print_rsp_tensor
  public print_rsp_tensor_stdout
  public print_rsp_tensor_stdout_tr

  private

  real(8) :: time_start
  real(8) :: time_end

  contains
  
  subroutine openrsp_get_property(n_props, np, pert_dims, pert_first_comp, pert_labels, n_freq_cfgs, pert_freqs, &
                                   kn_rules, F_unpert, S_unpert, D_unpert, get_rsp_sol, get_nucpot, &
                                   get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                                   get_2el_mat, get_2el_exp, get_xc_mat, & 
                                   get_xc_exp, id_outp, rsp_tensor, file_id)
    implicit none

    integer(kind=QINT), intent(in) :: n_props
    integer(kind=QINT), dimension(n_props), intent(in) :: np, n_freq_cfgs
    integer(kind=4), intent(in) :: id_outp
    integer(kind=QINT), dimension(sum(np)), intent(in) :: pert_dims, pert_first_comp
    character(4), dimension(sum(np)), intent(in) :: pert_labels
    character(256) :: filename
    integer :: i, j, k, m, n
    integer :: dum_ind
    integer, dimension(sum(n_freq_cfgs)) :: prop_sizes, num_blks
    integer(kind=QINT), intent(in), dimension(n_props) :: kn_rules
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    character, optional, dimension(20) :: file_id
    integer, allocatable, dimension(:) :: blk_sizes
    integer, allocatable, dimension(:,:) :: blk_info
    complex(8), dimension(dot_product(np, n_freq_cfgs)), intent(in) :: pert_freqs
    integer(kind=QINT) num_perts
    real :: timing_start, timing_end
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
    external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
    complex(8), dimension(*) :: rsp_tensor
    type(QcMat) :: S_unpert, D_unpert, F_unpert
    type(contrib_cache_outer), pointer :: S, D, F
    integer :: kn(2)
    logical :: r_exist, sdf_retrieved
    integer, dimension(3) :: rs_info
    integer, dimension(3) :: prog_info
    
    prog_info = (/0,0,0/)
    
    call prog_init(rs_info)
    
    call prog_incr(prog_info, 1)
    
    
    sdf_retrieved = .FALSE.

    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'S, D, F initialization stage was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
    
    
       allocate(S)
       allocate(D)
       allocate(F)
    
       call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE., 260)
       call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE., 260)
       call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE., 260)
    
       sdf_retrieved = .TRUE.
    
    else
    
       ! Set up S, D, F data structures
    
       call contrib_cache_outer_allocate(S)
       call contrib_cache_outer_allocate(D)
       call contrib_cache_outer_allocate(F)
    
       call contrib_cache_outer_add_element(S, .FALSE., 1, (/get_emptypert()/), &
            data_size = 1, data_mat=(/S_unpert/))
       call contrib_cache_outer_add_element(D, .FALSE., 1, (/get_emptypert()/), &
            data_size = 1, data_mat=(/D_unpert/))
       call contrib_cache_outer_add_element(F, .FALSE., 1, (/get_emptypert()/), &
            data_size = 1, data_mat=(/F_unpert/))
            
       call contrib_cache_outer_store(S, 'OPENRSP_S_CACHE')
       call contrib_cache_outer_store(D, 'OPENRSP_D_CACHE')
       call contrib_cache_outer_store(F, 'OPENRSP_F_CACHE')
    
    end if
    
    
    call prog_incr(prog_info, 1)

    ! Present calculation and initialize perturbation tuple datatypes and
    ! associated size/indexing information
    
    write(id_outp,*) ' '
    write(id_outp,*) 'OpenRSP lib called'
    write(id_outp,*) ' '
    
    if (n_props == 1) then
       write(id_outp,*) 'Calculating one property'
    else  
       write(id_outp,*) 'Calculating', n_props, 'properties'
    end if
    write(id_outp,*) ' '
    
   
    k = 1
    
    do i = 1, n_props
    
       write(id_outp,*) 'Property', i, 'is order', np(i)
       write(id_outp,*) 'The choice of k, n is', kn_rules(i), 'and', np(i) - 1 - kn_rules(i)
       write(id_outp,*) ' '
       write(id_outp,*) 'The number of components for each perturbation is:    ', pert_dims(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
       write(id_outp,*) 'The perturbation labels are:                          ', pert_labels(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
       write(id_outp,*) 'Number of frequency configurations:', n_freq_cfgs(i)
       write(id_outp,*) ' '
        
   
       do j = 1, n_freq_cfgs(i)
       
          kn_rule(k, 1) = kn_rules(i)
          kn_rule(k, 2) = np(i) - 1 - kn_rules(i)
          
          if ((kn_rule(k, 1) - kn_rule(k, 2) > 1)) then

             write(id_outp,*) 'ERROR: Invalid choice of (k,n)'
             write(id_outp,*) 'Valid choices for k are integers between and including 0 and ', (np(i) - 1)/2
             write(id_outp,*) 'Valid choices of n are such that k + n =', np(i) - 1
             write(id_outp,*) 'Cannot proceed with calculation: Exiting OpenRSP lib'
             write(id_outp,*) ' '
             return
 
          end if
         
          p_tuples(k)%npert = np(i)
          allocate(p_tuples(k)%pdim(np(i)))
          allocate(p_tuples(k)%plab(np(i)))
          allocate(p_tuples(k)%pid(np(i)))
          allocate(p_tuples(k)%freq(np(i)))
          
          p_tuples(k)%pdim = pert_dims(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
          p_tuples(k)%plab = pert_labels(sum(np(1:i)) - np(i) + 1:sum(np(1:i)))
          p_tuples(k)%pid = (/(m, m = 1, np(i))/)
          p_tuples(k)%freq = pert_freqs(dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + &
          1 + (j - 1)*np(i):dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + (j)*np(i))
       
          write(id_outp,*) 'Frequency configuration', j
          write(id_outp,*) ' '
          write(id_outp,*) 'Frequencies (real part):', (/(real(p_tuples(k)%freq(m)), m = 1, np(i))/)
          write(id_outp,*) 'Frequencies (imag. part):', (/(aimag(p_tuples(k)%freq(m)), m = 1, np(i))/)
          write(id_outp,*) ' '
          
          ! NOTE: ORDERING MAY CHANGE HERE AND COULD NEED REORDERING AFTER CALCULATION
          p_tuples(k) = p_tuple_standardorder(p_tuples(k))
          p_tuples(k)%pid = (/(m, m = 1, np(i))/)
       
          num_blks(k) = get_num_blks(p_tuples(k))
       
          write(*,*) 'Number of blocks:', num_blks(k)
       
          allocate(blk_info(num_blks(k), 3))
          allocate(blk_sizes(num_blks(k)))
          blk_info = get_blk_info(num_blks(k), p_tuples(k))
          write(*,*) 'Block info:', blk_info
          blk_sizes = get_triangular_sizes(num_blks(k), blk_info(1:num_blks(k), 2), &
                                           blk_info(1:num_blks(k), 3))

          write(*,*) 'Block sizes:', blk_sizes
                                           
          prop_sizes(k) = get_triangulated_size(num_blks(k), blk_info)
          
          write(*,*) 'Property size:', prop_sizes(k)
       
          deallocate(blk_info)
          deallocate(blk_sizes)
       
          k = k + 1
       
      
       end do

              
    end do

    rsp_tensor(1:sum(prop_sizes)) = 0.0
    
    ! Calculate properties
    
    write(id_outp,*) 'Starting clock: About to start property calculations'
    write(id_outp,*) ' '

    call cpu_time(timing_start)

    call get_prop(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, get_rsp_sol, &
                  get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                  get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, &
                  id_outp, prop_sizes, rsp_tensor, prog_info, rs_info, sdf_retrieved)

    call cpu_time(timing_end)

    write(id_outp,*) 'Clock stopped: Property calculations finished'
    write(id_outp,*) 'Time spent:',  timing_end - timing_start, ' seconds'
    write(id_outp,*) ' '

    ! Print results
    
    k = 1
    n = 1
    
    do i = 1, n_props
    
       do j = 1, n_freq_cfgs(i)
       
          write(filename, "(A10, I3, A1, I3)") 'rsp_tensor_', i, '_', j
          open(unit=260, file=filename, &
               status='replace', action='write') 
          open(unit=261, file=filename // '_human', &
               status='replace', action='write') 
            
          allocate(blk_info(num_blks(i), 3))
          allocate(blk_sizes(num_blks(i)))
          blk_info = get_blk_info(num_blks(i), p_tuples(k))
          blk_sizes = get_triangular_sizes(num_blks(i), blk_info(1:num_blks(i), 2), &
                                           blk_info(1:num_blks(i), 3))
            
          call print_rsp_tensor_tr(1, p_tuples(k)%npert, p_tuples(k)%pdim, &
          (/ (1, m = 1, (p_tuples(k)%npert - 1) ) /), num_blks(i), blk_sizes, &
          blk_info, prop_sizes(i), rsp_tensor(n:n+prop_sizes(i) - 1), 260, 261)

          deallocate(blk_info)
          deallocate(blk_sizes)
          
          close(260)
          close(261)
                    
          write(*,*) 'Property', i, j, ' was printed to rsp_tensor'
          write(*,*) 'Property (formatted print) was printed to rsp_tensor_human'
          
          k = k + 1
          n = n + prop_sizes(i)
    
       end do
       
       
    end do

    write(*,*) ' '
    write(*,*) 'End of print'

  end subroutine
   
  
  subroutine get_prop(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, get_rsp_sol, &
                  get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
                  get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, &
                  id_outp, prop_sizes, props, prog_info, rs_info, sdf_retrieved)

    implicit none

    
    logical :: traverse_end, sdf_retrieved, contrib_retrieved, props_retrieved
    integer :: n_props, id_outp, i, j, k
    integer, dimension(3) :: prog_info, rs_info
    integer, dimension(n_props) :: n_freq_cfgs
    integer, dimension(sum(n_freq_cfgs)) :: prop_sizes
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    type(p_tuple) :: emptypert
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), dimension(2) :: emptyp_tuples
    external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
    external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
    complex(8), dimension(*) :: props
    type(contrib_cache), pointer :: contribution_cache, cache_next
    type(contrib_cache_outer) :: F, D, S
    
    call empty_p_tuple(emptypert)
    emptyp_tuples = (/emptypert, emptypert/)
  
    call prog_incr(prog_info, 1)
  
    contrib_retrieved = .FALSE.
    props_retrieved = .FALSE.
  
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Perturbed overlap/density/Fock matrix stage was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if(.NOT.(sdf_retrieved)) then
       
          call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE.)
          call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE.)
          call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE.)
          sdf_retrieved = .TRUE.
          
       end if
  
    else
  
       ! Get all necessary F, D, S derivatives
     
       write(id_outp,*) ' '
       write(id_outp,*) 'Calculating perturbed overlap/density/Fock matrices'
       write(id_outp,*) ' '

       call cpu_time(time_start)
        
       call rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, &
                    get_rsp_sol, get_ovl_mat, get_1el_mat, &
                    get_2el_mat, get_xc_mat, .TRUE., id_outp, &
                    prog_info, rs_info, sdf_retrieved)
                     
       call cpu_time(time_end)

       write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
       write(id_outp,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
       write(id_outp,*) ' '
       
    end if
    
    
    call prog_incr(prog_info, 1)
    
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'HF energy-type contribution identification was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(contrib_retrieved)) then
       
          allocate(contribution_cache)
    
          call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
    
    else
    
       ! For each property: Recurse to identify HF energy-type contributions, store in cache
    
       call contrib_cache_allocate(contribution_cache)
    
       k = 1

       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)
       
             write(id_outp,*) ' '
             write(id_outp,*) 'Identifying HF-energy type contributions'
             write(id_outp,*) ' '

             call cpu_time(time_start)
             call rsp_energy_recurse(p_tuples(k), p_tuples(k)%npert, kn_rule(k,:), 1, (/emptypert/), &
                  0, D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, .TRUE., &
                  contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
             call cpu_time(time_end)

             write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
             write(id_outp,*) 'Finished identifying HF energy-type contributions'
             write(id_outp,*) ' '
          
             k = k + 1
       
          end do
       
       end do
    

       call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
    
    end if


    
    call prog_incr(prog_info, 1)
    
    
    
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'HF energy-type contribution calculation was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(contrib_retrieved)) then
       
          allocate(contribution_cache)
    
          call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
           
    else
    
       ! Calculate all identified contributions and store in cache
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Calculating HF-energy type contributions'
       write(id_outp,*) ' '
    
       call cpu_time(time_start)
    
       traverse_end = .FALSE.

       cache_next => contribution_cache

       ! Cycle cache
       do while (cache_next%last .eqv. .FALSE.)
          cache_next => cache_next%next
       end do
       
       cache_next => cache_next%next
       
       if (cache_next%p_inner%npert == 0) then
!           write(*,*) 'cycling dummy'
          cache_next => cache_next%next
       end if
       
       ! Traverse linked list and calculate
       do while (traverse_end .eqv. .FALSE.)
       
          write(*,*) 'Calculating contribution for inner perturbation tuple'
          write(*,*) cache_next%p_inner%plab
          write(*,*) ' '
          
          if (rs_check(prog_info, rs_info, lvl=2)) then
          
             write(*,*) ' '
             write(*,*) 'Calculation was completed in previous invocation: Passing to next stage'
             write(*,*) ' '
                
             ! Note: No cache retrieval here: In order to get to this position, the
             ! cache would already have been retrieved
          
          else

             call rsp_energy_calculate(D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, cache_next)
          
             call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
             
          end if
          
          call prog_incr(prog_info, 2)
          
          if (cache_next%last) then
             traverse_end = .TRUE.
          end if
          
          cache_next => cache_next%next
          
       end do

       call cpu_time(time_end)

       write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
       write(id_outp,*) 'Finished calculating HF energy-type contributions'
       write(id_outp,*) ' '
       
       call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
    
    end if
    
    call prog_incr(prog_info, 1)
   
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'HF energy-type contribution assembly was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(props_retrieved)) then
       
          call mat_scal_retrieve(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
          props_retrieved = .TRUE.
          
       end if
           
    else
    
       ! For each property: Recurse to identify HF energy-type contributions and 
       ! add to the property under consideration

       k = 1
    
       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)

             write(id_outp,*) ' '
             write(id_outp,*) 'Assembling HF-energy type contributions'
             write(id_outp,*) ' '

             call cpu_time(time_start)
             call rsp_energy_recurse(p_tuples(k), p_tuples(k)%npert, kn_rule(k,:), 1, (/emptypert/), &
                  0, D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, .FALSE., &
                  contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
             call cpu_time(time_end)

             write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
             write(id_outp,*) 'Finished assembling HF energy-type contributions'
             write(id_outp,*) ' '

             write(*,*) 'Property sample', props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1: &
             min(sum(prop_sizes(1:k)) - prop_sizes(k) + 100, sum(prop_sizes(1:k))))
          
             k = k + 1
       
          end do
        
       end do
       
       call mat_scal_store(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)

    end if
    
    deallocate(contribution_cache)
    
    contrib_retrieved = .FALSE.
    
    call prog_incr(prog_info, 1)
    
    ! For each property: Recurse to identify two-factor contributions and store in cache
    
    
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Two-factor type contribution identification was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(contrib_retrieved)) then
       
          allocate(contribution_cache)
    
          call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
           
    else
    
       k  = 1
    
       call contrib_cache_allocate(contribution_cache)
    
       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)
          
             write(*,*) 'This kn rule', kn_rule(k,:)
       
             write(id_outp,*) ' '
             write(id_outp,*) 'Identifying two-factor contributions'
             write(id_outp,*) ' '

             call cpu_time(time_start)
          
             call rsp_twofact_recurse(p_tuples(k), &
                  kn_rule(k,:), (/emptypert, emptypert/), &
                  .TRUE., contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
             call cpu_time(time_end)

             write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
             write(id_outp,*) 'Finished identifying two-factor contributions'
             write(id_outp,*) ' '
          
             k = k + 1
       
          end do
       
       end do
    
       call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
    
    end if
    
    call prog_incr(prog_info, 1)
    
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Two-factor type contribution calculation was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(contrib_retrieved)) then
       
          allocate(contribution_cache)
    
          call contrib_cache_retrieve(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          contrib_retrieved = .TRUE.
          
       end if
       
       
    else
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Calculating two-factor contributions'
       write(id_outp,*) ' '
    
       call cpu_time(time_start)
    
       traverse_end = .FALSE.

       cache_next => contribution_cache

       ! Cycle cache
       do while (cache_next%last .eqv. .FALSE.)
          cache_next => cache_next%next
       end do

       cache_next => cache_next%next
       
       if (cache_next%p_inner%npert == 0) then
!           write(*,*) 'cycling dummy'
          cache_next => cache_next%next
       end if
       
       ! Traverse linked list and calculate
       do while (traverse_end .eqv. .FALSE.)
       
          write(*,*) 'Calculating contribution for factor 1 tuple'
          write(*,*) cache_next%p_inner%plab
          write(*,*) ' '
          
          if (rs_check(prog_info, rs_info, lvl=2)) then
          
             write(*,*) ' '
             write(*,*) 'Calculation was completed in previous invocation: Passing to next stage'
             write(*,*) ' '
                
             ! Note: No cache retrieval here: In order to get to this position, the
             ! cache would already have been retrieved
          
          else

             call rsp_twofact_calculate(S, D, F, get_ovl_exp, cache_next)
             
             call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
          
          end if
          
          call prog_incr(prog_info, 2)
          
          if (cache_next%last) then
             traverse_end = .TRUE.
          end if
          
          cache_next => cache_next%next
          
       end do

       call cpu_time(time_end)

       write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
       write(id_outp,*) 'Finished calculating two-factor contributions'
       write(id_outp,*) ' '
       
       call contrib_cache_store(contribution_cache, 'OPENRSP_CONTRIB_CACHE')
       
    end if

    call prog_incr(prog_info, 1)
    
    
    
    if (rs_check(prog_info, rs_info, lvl=1)) then
    
       write(id_outp,*) ' '
       write(id_outp,*) 'Two-factor type contribution assembly was completed'
       write(id_outp,*) 'in previous invocation: Passing to next stage of calculation'
       write(id_outp,*) ' '
       
       if (.NOT.(props_retrieved)) then
       
          call mat_scal_retrieve(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
          props_retrieved = .TRUE.
       
       end if
       
    else
    
       ! For each property: Recurse to identify two-factor contributions and 
       ! add to the property under consideration
       
       k = 1
    
       do i = 1, n_props
    
          do j = 1, n_freq_cfgs(i)

             write(id_outp,*) ' '
             write(id_outp,*) 'Assembling two-factor contributions'
             write(id_outp,*) ' '

             call cpu_time(time_start)
          
             call rsp_twofact_recurse(p_tuples(k), &
                  kn_rule(k,:), (/emptypert, emptypert/), &
                  .FALSE., contribution_cache, prop_sizes(k), &
                  props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
             call cpu_time(time_end)

             write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
             write(id_outp,*) 'Finished assembling two-factor contributions'
             write(id_outp,*) ' '

             k = k + 1
          
          end do
       
       end do
       
       call mat_scal_store(sum(prop_sizes), 'OPENRSP_PROP_CACHE', scal=props)
  
    end if
  
    deallocate(contribution_cache)
    
    k = 1
    
    do i = 1, n_props
    
       do j = 1, n_freq_cfgs(i)
          
          write(*,*) 'Property', i, ', freq. config', j
          write(*,*) props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1: &
                                      sum(prop_sizes(1:k)))
          write(*,*) ' '
          
          
          k = k + 1
          
       end do
       
    end do
    
    call prog_incr(prog_info, 1)
   
    
  end subroutine
   

   recursive subroutine rsp_energy_recurse(pert, total_num_perturbations, kn, num_p_tuples, &
                                  p_tuples, density_order, D, get_nucpot, get_1el_exp, &
                                  get_t_exp, get_2el_exp, dryrun, cache, p_size, prop)

    implicit none

    logical :: e_knskip, dryrun, traverse_end
    type(p_tuple) :: pert
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, id_outp
    integer :: p_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(contrib_cache_outer) :: D
    type(contrib_cache), target :: cache
    type(contrib_cache), pointer :: cache_next
    complex(8), dimension(p_size), optional :: prop
    external :: get_nucpot, get_1el_exp, get_t_exp, get_2el_exp


    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'


    if (p_tuples(1)%npert == 0) then

       call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples, (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
       density_order, D, get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, &
       dryrun, cache, p_size=p_size, prop=prop)

    else

       call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations,  &
       kn, num_p_tuples, (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
       p_tuples(2:size(p_tuples))/), density_order, D,  &
       get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, &
       dryrun, cache, p_size=p_size, prop=prop)

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
          kn, num_p_tuples, t_new, density_order + 1, D, &
          get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, dryrun, cache, p_size=p_size, prop=prop)

       end do

       ! Since we are only calculating Hartree-Fock type energy terms here,
       ! we don't need to go beyond to perturbed contraction density matrices
       ! (but that is in general needed for XC contribs)
       if (num_p_tuples < 3) then

          ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
          ! a(nother) pert D contraction)

          call rsp_energy_recurse(p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
          density_order + 1, D, get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, &
          dryrun, cache, p_size=p_size, prop=prop)

       end if

    ! At the final recursion level: Calculate the contrib (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allows it)

    else

       p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)
    
       e_knskip = .FALSE.

       do i = 1, num_p_tuples
 
          if (i > 1) then

             if(kn_skip(p_tuples(i)%npert, p_tuples(i)%pid, kn)) then

                e_knskip = .TRUE.

             end if
          
          elseif (i == 1) then
          
          

          end if

       end do


       if (e_knskip .EQV. .FALSE.) then
       
          if (contrib_cache_already(cache, num_p_tuples, p_tuples)) then
          
             if (.NOT.(dryrun)) then

                write(*,*) 'Cache retrieval: getting contribution'
                
                do i = 1, num_p_tuples

                   if (i == 1) then
    
                      write(*,*) 'E', p_tuples(i)%pid
    
                   else 
    
                      write(*,*) 'D', p_tuples(i)%pid
                       
                   end if
    
                end do

                    write(*,*) ''
             
                ! NOTE (MaR): EVERYTHING MUST BE STANDARD ORDER IN 
                ! THIS CALL (LIKE property_cache_getdata ASSUMES)
                call contrib_cache_getdata(cache, num_p_tuples, &
                   p_tuples_standardorder(num_p_tuples, p_tuples), p_size, 0, scal=prop)

             end if

          else

             if (dryrun) then
             
                write(*,*) 'Adding cache element'
!                 write(*,*) 'num p tuples when adding', num_p_tuples
             
                do i = 1, num_p_tuples

                   if (i == 1) then
    
                      write(*,*) 'E', p_tuples(i)%pid
    
                   else 
    
                      write(*,*) 'D', p_tuples(i)%pid
                       
                   end if
    
                end do
             
             
                call contrib_cache_add_element(cache, num_p_tuples, &
                     p_tuples_standardorder(num_p_tuples, p_tuples))

             end if

          end if

       end if

    end if

  end subroutine


  subroutine rsp_energy_calculate(D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, cache)

    implicit none

    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n, p, offset
    integer :: id_outp, c1_ctr, c2_ctr, lhs_ctr_1, lhs_ctr_2, rhs_ctr_2
    integer :: total_outer_size_1, total_outer_size_2
    integer :: num_0, num_1, num_pert, tot_num_pert
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
    character(30) :: mat_str, fmt_str
    type(contrib_cache_outer) :: D
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: outer_next
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, LHS_dmat_2, RHS_dmat_2
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    integer, allocatable, dimension(:,:) :: outer_contract_sizes_2, blk_sizes
    complex(8), allocatable, dimension(:) :: contrib_0, contrib_1, contrib_2, data_tmp
    external :: get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp
    
    ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)
    
    outer_next => cache%contribs_outer
    
    write(*,*) 'Number of outer contributions: ', cache%num_outer

    
    allocate(outer_contract_sizes_1(cache%num_outer))
    allocate(outer_contract_sizes_2(cache%num_outer,2))
        
   
    ! Traversal: Find number of density matrices for contraction for nuc-nuc, 1-el, 2-el cases
    
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
       
      
    total_outer_size_1 = 0
    total_outer_size_2 = 0
    num_0 = 0
    num_1 = 0
        
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
       write(*,*) 'Outer contribution:'
    
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'D', outer_next%p_tuples(i)%pid
       
       end do
       
       write(*,*) ' '
  
  
       if (outer_next%num_dmat == 0) then

          num_0 = 1
          num_1 = num_1 + 1
          outer_contract_sizes_1(k) = 1
          outer_contract_sizes_2(k, :) = (/1,1/)
          
          total_outer_size_1 = total_outer_size_1 + 1
          total_outer_size_2 = total_outer_size_2 + 1
          
       else if (outer_next%num_dmat == 1) then
       
          num_1 = num_1 + 1
       
          outer_contract_sizes_1(k) = outer_next%blks_tuple_triang_size(1)
          outer_contract_sizes_2(k, :) = (/outer_next%blks_tuple_triang_size(1),1/)
          
          total_outer_size_1 = total_outer_size_1 + outer_next%blks_tuple_triang_size(1)
          total_outer_size_2 = total_outer_size_2 + outer_next%blks_tuple_triang_size(1)
       
       else if (outer_next%num_dmat == 2) then
       
          outer_contract_sizes_1(k) = 0
          outer_contract_sizes_2(k, :) = (/outer_next%blks_tuple_triang_size(1), &
                                         outer_next%blks_tuple_triang_size(2)/)
          
          total_outer_size_2 = total_outer_size_2 + outer_next%blks_tuple_triang_size(1)*outer_next%blks_tuple_triang_size(2)
       
       end if
   

    
       
    
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do
 
    ! Make collapsed contraction sizes array for 1-el call
 
    allocate(outer_contract_sizes_1_coll(num_1))
    
!     write(*,*) 'num_1', num_1
!     write(*,*) 'outer contract sizes 1', outer_contract_sizes_1
    
    k = 1 
     do i = 1, cache%num_outer
        if (outer_contract_sizes_1(i) > 0) then
           outer_contract_sizes_1_coll(k) = outer_contract_sizes_1(i)
          k = k + 1
        end if
    end do
    
    ! Allocate and set up outer
    
    allocate(LHS_dmat_1(sum(outer_contract_sizes_1(:))))
    allocate(LHS_dmat_2(sum(outer_contract_sizes_2(:, 1))))
    allocate(RHS_dmat_2(sum(outer_contract_sizes_2(:, 2))))
    
    do i = 1, size(LHS_dmat_1)
    
       call QcMatInit(LHS_dmat_1(i))
    
    end do
    
    do i = 1, size(LHS_dmat_2)
    
       call QcMatInit(LHS_dmat_2(i))
    
    end do
    
    do i = 1, size(RHS_dmat_2)
    
       call QcMatInit(RHS_dmat_2(i))
    
    end do
    
    
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
       
    k = 1
    lhs_ctr_1 = 1
    lhs_ctr_2 = 1
    rhs_ctr_2 = 1
    
    do while (traverse_end .EQV. .FALSE.)
    
       ! No chain rule applications
       if (outer_next%num_dmat == 0) then
       
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_1(lhs_ctr_1))
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_2(lhs_ctr_2))
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))

       ! One chain rule application
       else if (outer_next%num_dmat == 1) then
       
       
          do m = 1, outer_contract_sizes_1(k) 
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                  1, ind_len=size(outer_next%indices, 2), ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_1(lhs_ctr_1 + m  - 1))
          end do
       
          do m = 1, outer_contract_sizes_2(k, 1) 
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                  1, ind_len=size(outer_next%indices, 2), ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
          end do
          
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))
       
       ! Two chain rule applications
       else if (outer_next%num_dmat == 2) then
       
       
          do m = 1, outer_contract_sizes_2(k, 1) 
          
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                  1, ind_len=outer_next%p_tuples(1)%npert, &
                  ind_unsorted=outer_next%indices(1 + &
                  (m - 1) * outer_next%blks_tuple_triang_size(2), &
                  1:outer_next%p_tuples(1)%npert), &
                  mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
          end do
          
          do n = 1, outer_contract_sizes_2(k, 2) 
          
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(2)/), .FALSE., &
                  1, ind_len=outer_next%p_tuples(2)%npert, &
                  ind_unsorted=outer_next%indices(n, outer_next%p_tuples(1)%npert + 1: &
                  outer_next%p_tuples(1)%npert + outer_next%p_tuples(2)%npert), &
                  mat_sing=RHS_dmat_2(rhs_ctr_2 + n  - 1))
          end do          
       
       end if
   
       if (outer_next%last) then
          traverse_end = .TRUE.
       end if

       lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
       lhs_ctr_2 = lhs_ctr_2 + outer_contract_sizes_2(k, 1)
       rhs_ctr_2 = rhs_ctr_2 + outer_contract_sizes_2(k, 2)
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    
    allocate(contrib_0(cache%blks_triang_size))
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    allocate(contrib_2(cache%blks_triang_size*total_outer_size_2))
    
    ! Calculate contributions
    
    ! Calculate nuclear-nuclear repulsion contribution
    if (num_0 > 0) then
    
       contrib_0 = 0.0
       call get_nucpot(num_pert, pert_ext, size(contrib_0), contrib_0)
       
       write(*,*) 'nucpot contribution: ', contrib_0(1:min(12, size(contrib_0)))
    
    end if
    
    ! Calculate one-electron contributions
    if (num_1 > 0) then
    
       contrib_1 = 0.0
       call get_1el_exp(num_pert, pert_ext, total_outer_size_1, &
                        LHS_dmat_1, size(contrib_1), contrib_1)
      
       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()
      
!        call rsp_ovlave_t_matrix_2014(get_ovl_exp, cache%p_inner, cache%p_inner%npert, &
!                                 t_matrix_bra, t_matrix_ket, outer_contract_sizes_1_coll, &
!                                 LHS_dmat_1, size(contrib_1), contrib_1)
    
    write(*,*) '1-el contribution: ', contrib_1(1:min(12, size(contrib_1)))
    
    end if
    
!     write(*,*) '2-el'
    
    ! Calculate two-electron contributions
    contrib_2 = 0.0
    call get_2el_exp(num_pert, pert_ext, cache%num_outer, outer_contract_sizes_2(:, 1), LHS_dmat_2, & 
                     outer_contract_sizes_2(:, 2), RHS_dmat_2, size(contrib_2), contrib_2)
                       
    
    write(*,*) '2-el contribution: ', contrib_2(1:min(12, size(contrib_2)))
    
    ! Traversal: Add nuc-nuc, 1-el and two-el contributions together (put in contrib_2)
    
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
      
    k = 1
    
    c1_ctr = 1
    c2_ctr = 1
    
    
    do while (traverse_end .EQV. .FALSE.)
  
       ! Nuc-nuc, one-el and two-el contribution
       if (outer_next%num_dmat == 0) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size))
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
       else if (outer_next%num_dmat == 1) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
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
       else if (outer_next%num_dmat == 2) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)))
          allocate(data_tmp(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)))
                   
          data_tmp = contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size * &
                            outer_contract_sizes_2(k, 1) * outer_contract_sizes_2(k, 2) - 1)
       
          c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)
       end if
       
       if (outer_next%num_dmat == 0) then
       
       
          outer_next%data_scal = data_tmp
          
       else
       
          
       
          if (cache%p_inner%npert > 0) then
          
             tot_num_pert = cache%p_inner%npert + &
             sum((/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/))
                   
             allocate(blks_tuple_info(outer_next%num_dmat + 1,tot_num_pert, 3))
             allocate(blk_sizes(outer_next%num_dmat + 1, tot_num_pert))
             
             blks_tuple_info = 0
             blk_sizes = 0
                
             do j = 1, outer_next%num_dmat + 1
                
                if (j == 1) then
                
                   do m = 1, cache%nblks
                   
                      blks_tuple_info(j, m, :) = cache%blk_info(m, :)
                      
                   end do
                   
                   blk_sizes(j, 1:cache%nblks) = cache%blk_sizes
                
                
                else
                
                   do m = 1, outer_next%nblks_tuple(j - 1)
                
                      do p = 1, 3
                   
                         blks_tuple_info(j, m, :) = outer_next%blks_tuple_info(j - 1, m, :)
                   
                      end do
                
                   end do
                   
                   blk_sizes(j, 1:outer_next%nblks_tuple(j-1)) = &
                   outer_next%blk_sizes(j-1, 1:outer_next%nblks_tuple(j-1))
                   
                end if
                
             end do
       
             do i = 1, size(outer_next%indices, 1)
          
                do j = 1, size(cache%indices, 1)
                
!                    write(*,*) 'getting size', blk_sizes
!                    write(*,*) 'cache blk sizes', cache%blk_sizes
!                    write(*,*) 'outer blk sizes', (/(outer_next%blk_sizes(m,:), m = 1, outer_next%num_dmat)/)
                   
             
                   offset = get_triang_blks_tuple_offset(outer_next%num_dmat + 1, &
                   cache%p_inner%npert + sum((/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/)), &
                   (/cache%nblks, (/(outer_next%nblks_tuple(m), m = 1, outer_next%num_dmat) /) /), &
                   (/cache%p_inner%npert, (/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/)/), &
                   blks_tuple_info, &
                   blk_sizes, &
                   (/cache%blks_triang_size, &
                   (/(outer_next%blks_tuple_triang_size(m), m = 1, outer_next%num_dmat)/)/), &
                   (/cache%indices(j, :), outer_next%indices(i, :)/))
                
                   outer_next%data_scal(offset) = data_tmp(j + size(cache%indices, 1) * (i - 1))
                            
             
                end do
          
             end do
             
             deallocate(blk_sizes)
             deallocate(blks_tuple_info)
          
          else
          
             write(*,*) 'ERROR: UNEXPECTED: NO INNER PERTURBATIONS'
          
          
          end if
                 
       end if
       
       deallocate(data_tmp)
   
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do

        
    deallocate(outer_contract_sizes_2)
    
  end subroutine

  
   recursive subroutine rsp_twofact_recurse(pert, kn, p12, dryrun, cache, p_size, prop)

    implicit none

    logical :: dryrun, lag_eligible
    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: curr_outer
    integer ::  i, j
    integer :: hard_offset
    integer, dimension(2) :: kn
    integer :: nblks, block_size, p_size
    integer, allocatable, dimension(:,:,:) :: blk_info
    complex(8), dimension(p_size) :: prop
    
    if (pert%npert > 0) then

       call rsp_twofact_recurse(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), dryrun, &
       cache, p_size, prop)
       
       call rsp_twofact_recurse(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), dryrun, &
       cache, p_size, prop)

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

       if (kn_skip(p12(2)%npert, p12(2)%pid, kn) .EQV. .FALSE.) then

          if (contrib_cache_already(cache, 2, p12, n_rule=kn(2))) then

             if (.NOT.(dryrun)) then

             ! FIXME: CHANGE THIS TO WORK WITH CONTRIBUTION TYPE
             
                write(*,*) 'Retrieving Pulay n contribution:'
                write(*,*) 'S', p12(1)%pid
                write(*,*) 'W', p12(2)%pid
             
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
             else
             
                call contrib_cache_cycle_outer(cache, 2, p12, curr_outer, n_rule=kn(2))
                if (curr_outer%contrib_type == 3) then
                
                   curr_outer%contrib_type = 4
                
                end if
             
             end if
       
          else
          
             if (dryrun) then
             
                write(*,*) 'Identified Pulay n contribution:'
                write(*,*) 'S', p12(1)%pid
                write(*,*) 'W', p12(2)%pid
                
!                 write(*,*) 'kn rule for cache creation', kn
                
                call contrib_cache_add_element(cache, 2, p12, n_rule=kn(2))
                call contrib_cache_cycle_outer(cache, 2, p12, curr_outer, n_rule=kn(2))
                curr_outer%contrib_type = 1
             
             else
             
                write(*,*) 'ERROR: Expected to find Pulay contribution but it was not present'
                write(*,*) 'S', p12(1)%pid
                write(*,*) 'W', p12(2)%pid
             
             end if

          end if

       end if
       
       lag_eligible = .FALSE.
       
       do j = 1, p12(1)%npert
       
          if (p12(1)%pid(j) == 1) then
          
             lag_eligible = .TRUE.
          
          end if
       
       end do
       
       if ((kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) .AND. &
           (p12(1)%npert > 0) .AND. lag_eligible) then

          if (contrib_cache_already(cache, 2, p12, n_rule=kn(2))) then

             if (.NOT.(dryrun)) then
             
             ! FIXME: CHANGE THIS TO WORK WITH CONTRIBUTION TYPE
             
                write(*,*) 'Retrieving Lagrange contributions:'
                write(*,*) 'A', p12(1)%pid
                write(*,*) 'B', p12(2)%pid
                
                ! Pulay Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'Prop after Pulay Lagrange retrieval', prop(1:min(12, size(prop)))
                
!                 Idempotency Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'Prop after idempotency Lagrange retrieval', prop(1:min(12, size(prop)))

!                 SCFE Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                
!                 write(*,*) 'Prop after SCFE Lagrange retrieval', prop(1:min(12, size(prop)))
            
             else
             
                call contrib_cache_cycle_outer(cache, 2, p12, curr_outer, n_rule=kn(2))
                if (curr_outer%contrib_type == 1) then
                
                   curr_outer%contrib_type = 4
                
                end if
            
             end if
       
          else
          
             if (dryrun) then
             
                write(*,*) 'Identified Lagrange contributions:'
                write(*,*) 'A', p12(1)%pid
                write(*,*) 'B', p12(2)%pid
                
                call contrib_cache_add_element(cache, 2, p12, n_rule=kn(2))
                
                call contrib_cache_cycle_outer(cache, 2, p12, curr_outer, n_rule=kn(2))
                curr_outer%contrib_type = 3
             
             else
             
                write(*,*) 'ERROR: Expected to find Lagrange contributions but they were not present'
                write(*,*) 'A', p12(1)%pid
                write(*,*) 'B', p12(2)%pid
             
             end if
          


          end if

       end if 

    end if

  end subroutine
  
  
  subroutine rsp_twofact_calculate(S, D, F, get_ovl_exp, cache)

    implicit none

    logical :: traverse_end, any_lagrange
    integer :: cache_offset, i, j, k, m, n, p, c_ctr, c_snap, lagrange_max_n
    integer :: id_outp, i_supsize, o_triang_size, offset, tot_num_pert
    integer :: ctr_lagrange, ctr_pulay_n, o_ctr, size_lagrange, size_pulay_n
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
    type(contrib_cache_outer) :: S, D, F
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: outer_next
    type(QcMat), allocatable, dimension(:) :: Lambda, Zeta, W
    type(QcMat) :: Y, Z, D_unp
           
    integer, allocatable, dimension(:) :: pert_ext
    
    external :: get_ovl_exp
    
    ! Getting unperturbed D for template
    
    call QCMatInit(D_unp)

    call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
         contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)
    
    ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    ! Unsure about npert argument
    call p_tuple_to_external_tuple(cache%p_inner, cache%p_inner%npert, pert_ext)

    
    
    outer_next => cache%contribs_outer

    
!     write(*,*) 'num outer', cache%num_outer

    i_supsize = 3**cache%p_inner%npert
    
    
    allocate(o_supsize(cache%num_outer))
    allocate(o_supsize_prime(cache%num_outer))
    allocate(o_size(cache%num_outer))

    ! Prepare matrices for inner tuple
    
    i_supsize = derivative_superstructure_getsize(p_tuple_remove_first(cache%p_inner), &
                (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(d_struct_inner(i_supsize, 3))

    sstr_incr = 0
    
    call derivative_superstructure(p_tuple_remove_first(cache%p_inner), &
          (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
          (/get_emptypert(), get_emptypert(), get_emptypert()/), &
          i_supsize, sstr_incr, d_struct_inner)  
          
    sstr_incr = 0
          
    
    ! Traverse outer elements to determine which contributions are Pulay n and/or Lagrange
    ! Allocate accordingly: One array for Pulay n/Pulay Lagrange combined, one each for last two
    ! Alternatively, put directly in cache as they are calculated
    ! Do in turn (later introduce memory management for memory-intensive jobs): 
    !
    ! - Make all W matrices
    ! - Do S * W contraction and store
    ! - Make all Zeta matrices
    ! - Make Z in turn, contract and store
    ! - Make all Lambda matrices
    ! - Make Y in turn, contract and store
   
    ! Traversal: Find number of matrices/terms for contraction
    
    any_lagrange = .FALSE.
    lagrange_max_n = 0    
    
    traverse_end = .FALSE.
    
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
        
    size_pulay_n = 0
    size_lagrange = 0
    
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
       if (outer_next%p_tuples(1)%npert == 0) then
       
          o_triang_size = 1
       
       else
       
          o_triang_size = outer_next%blks_tuple_triang_size(1)
       
       end if
  
  
       write(*,*) 'Outer contribution, type', outer_next%contrib_type
    
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'B', outer_next%p_tuples(i)%pid
       
       end do
    
       write(*,*) ' '
    
       ! Here and elsewhere: k in kn rule does not matter as long as it 
       ! is > n; it will always > n so used like this
       
       o_supsize(k) = derivative_superstructure_getsize(outer_next%p_tuples(1), &
                   (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/))
                   
       o_supsize_prime(k) = derivative_superstructure_getsize(outer_next%p_tuples(1), &
                   (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                   (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
       
       o_size(k) = outer_next%contrib_type * o_triang_size
       
       if ((outer_next%contrib_type == 1) .OR. (outer_next%contrib_type == 4)) then
       
          size_pulay_n = size_pulay_n + o_triang_size
          
          
       
       end if
       
       if (outer_next%contrib_type >= 3) then
       
          size_lagrange = size_lagrange + 3 * o_triang_size
          
          any_lagrange = .TRUE.
          
          if (outer_next%n_rule > lagrange_max_n) then
          
             lagrange_max_n = outer_next%n_rule            
          
          end if
       
       end if
       
       allocate(outer_next%data_scal(o_size(k) * cache%blks_triang_size))
       
       outer_next%data_scal = 0.0
         
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
     
    end do
    
    
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
    
    allocate(which_index_is_pid(20))
    
    ! Get zeta and lambda matrices if applicable
    if (any_lagrange) then
    
       allocate(Lambda(cache%blks_triang_size))
       allocate(Zeta(cache%blks_triang_size))
    
       ! FIXME: Change to max order of all properties
    
       which_index_is_pid = 0
          
       do i = 1, cache%p_inner%npert
          
          which_index_is_pid(cache%p_inner%pid(i)) = i
          
       end do
          
       do i = 1, cache%blks_triang_size
    
          ! CONTINUE HERE: CREATE INNER INDICES WHEN MAKING CACHE ENTRY (IT IS MISSING NOW)
          call QcMatInit(Lambda(i), D_unp)
          call QcMatInit(Zeta(i), D_unp)
          call QcMatZero(Lambda(i))
          call QcMatZero(Zeta(i))
       
          call rsp_get_matrix_zeta(p_tuple_getone(cache%p_inner, 1), (/lagrange_max_n, &
               lagrange_max_n/), i_supsize, d_struct_inner, maxval(cache%p_inner%pid), &
               which_index_is_pid(1:maxval(cache%p_inner%pid)), size(cache%indices(i,:)), &
               cache%indices(i,:), F, D, S, Zeta(i))
                  
          call rsp_get_matrix_lambda(p_tuple_getone(cache%p_inner, 1), i_supsize, &
               d_struct_inner, maxval(cache%p_inner%pid), &
               which_index_is_pid(1:maxval(cache%p_inner%pid)), &
               size(cache%indices(i,:)), cache%indices(i,:), D, S, Lambda(i))
    
       end do
    
    
    end if    
    
    ! Traversal: Make W matrices and store
    
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if

    ctr_pulay_n = 0
    ctr_lagrange = 0
    o_ctr = 1
    
    allocate(W(size_pulay_n + size_lagrange/3))
    
    do i = 1, size(W)
    
       call QcMatInit(W(i), D_unp)
       call QcMatZero(W(i))
       
    end do
    
    allocate(contrib_pulay(cache%blks_triang_size* ( size_pulay_n + size_lagrange/3)))
       
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
!        write(*,*) 'Outer contribution:'!, outer_next%num_dmat, outer_next%dummy_entry
!     
!        do i = 1, outer_next%num_dmat
!           
!           write(*,*) 'B', outer_next%p_tuples(i)%pid
!        
!        end do
    
          allocate(d_struct_o(o_supsize(k), 3))
          allocate(d_struct_o_prime(o_supsize_prime(k), 3))

          
          which_index_is_pid = 0
          
          do j = 1, outer_next%p_tuples(1)%npert
          
             which_index_is_pid(outer_next%p_tuples(1)%pid(j)) = j
          
          end do
          
          sstr_incr = 0
      
          if(outer_next%p_tuples(1)%npert ==0) then
          
             if (outer_next%contrib_type == 1 .OR. outer_next%contrib_type == 4) then
                    
                call derivative_superstructure(get_emptypert(), &
                (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                o_supsize(k), sstr_incr, d_struct_o)
                
                sstr_incr = 0
             
             end if
            
             if (outer_next%contrib_type == 3 .OR. outer_next%contrib_type == 4) then
             
                call derivative_superstructure(get_emptypert(), &
                (/outer_next%n_rule, outer_next%n_rule/), .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                o_supsize_prime(k), sstr_incr, d_struct_o)
                
                sstr_incr = 0

             end if
             
             o_triang_size = 1
          
          else
      
             if (outer_next%contrib_type == 1 .OR. outer_next%contrib_type == 4) then             
      
                call derivative_superstructure(outer_next%p_tuples(1), &
                (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                o_supsize(k), sstr_incr, d_struct_o)
                
                sstr_incr = 0
             
             end if

             if (outer_next%contrib_type == 3 .OR. outer_next%contrib_type == 4) then
             
                call derivative_superstructure(outer_next%p_tuples(1), &
                (/outer_next%n_rule, outer_next%n_rule/), .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                o_supsize_prime(k), sstr_incr, d_struct_o_prime)
                
                sstr_incr = 0
                
             end if

             o_triang_size = size(outer_next%indices, 1)
             
          end if
            
          do j = 1, o_triang_size
          
          
             select case (outer_next%contrib_type)
             
             case (1)
             
!                 write(*,*) 'supsize', o_supsize(k)
!                 write(*,*) 'size outer ind', size(outer_next%indices(j,:))
!                 write(*,*) 'outer_next%indices(j,:)', outer_next%indices(j,:)
!                 write(*,*) 'inner npert', cache%p_inner%npert
!                 write(*,*) 'outer npert', outer_next%p_tuples(1)%npert
!                 write(*,*) 'which ind is pid', which_index_is_pid(1:cache%p_inner%npert + &
!                      outer_next%p_tuples(1)%npert)
!                 write(*,*) 'o ctr', o_ctr
!                 write(*,*) 'size of W', size(W)
             
                call rsp_get_matrix_w(o_supsize(k), d_struct_o, cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert, &
                     which_index_is_pid(1:cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                     F, D, S, W(o_ctr))
             
             case (3)
             
                call rsp_get_matrix_w(o_supsize_prime(k), d_struct_o_prime, &
                     cache%p_inner%npert + outer_next%p_tuples(1)%npert, &
                     which_index_is_pid(1:cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                     F, D, S, W(o_ctr))
             
             
             case (4)
             
                call rsp_get_matrix_w(o_supsize(k), d_struct_o, cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert, &
                     which_index_is_pid(1:cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                     F, D, S, W(o_ctr))
                     
                call rsp_get_matrix_w(o_supsize_prime(k), d_struct_o_prime, &
                     cache%p_inner%npert + outer_next%p_tuples(1)%npert, &
                     which_index_is_pid(1:cache%p_inner%npert + &
                     outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                     F, D, S, W(o_ctr + o_triang_size))
             
             end select
                  
             o_ctr = o_ctr + 1
       
          end do
       
       
          deallocate(d_struct_o)
          deallocate(d_struct_o_prime)
       
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do
    
    ! Outside traversal: Calculate contributions
    
    contrib_pulay = 0.0
    
    call get_ovl_exp(0, nof, 0, nof, size(pert_ext), pert_ext, size(W), W, &
                     size(contrib_pulay), contrib_pulay)
    
    
    contrib_pulay = -2.0 * contrib_pulay
    
    ! Traversal: Store Pulay contributions, calculate/store idempotency/SCFE contributions
    
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if

    call QcMatInit(Z, D_unp)
    call QcMatZero(Z)
    call QcMatInit(Y, D_unp)
    call QcMatZero(Y)
    
    k = 1
    o_ctr = 0
    
    do while (traverse_end .EQV. .FALSE.)

!        write(*,*) 'Outer contribution:'!, outer_next%num_dmat, outer_next%dummy_entry
!     
!        do i = 1, outer_next%num_dmat
!           
!           write(*,*) 'B', outer_next%p_tuples(i)%pid
!        
!        end do
    
    
    
       c_ctr = 0
    
       if (outer_next%p_tuples(1)%npert ==0) then
          
          o_triang_size = 1
          
       else
      
          o_triang_size = size(outer_next%indices, 1)
             
       end if

       
       
       tot_num_pert = cache%p_inner%npert + &
       sum((/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/))
                   
       allocate(blks_tuple_info(outer_next%num_dmat + 1,tot_num_pert, 3))
       allocate(blk_sizes(outer_next%num_dmat + 1, tot_num_pert))
               
       do j = 1, outer_next%num_dmat + 1
          
          if (j == 1) then
             
             do m = 1, cache%nblks
                
                blks_tuple_info(j, m, :) = cache%blk_info(m, :)
             
             end do
             
             blk_sizes(j, 1:cache%nblks) = cache%blk_sizes
          
          else
             
             if (size(outer_next%p_tuples) > 0) then
                
                if (outer_next%p_tuples(1)%npert > 0) then
                   
                   do m = 1, outer_next%nblks_tuple(j - 1)
                      
                      do p = 1, 3
                         
                         blks_tuple_info(j, m, :) = outer_next%blks_tuple_info(j - 1, m, :)
                      
                      end do
                   
                   end do
                   
                   blk_sizes(j, 1:outer_next%nblks_tuple(j-1)) = &
                   outer_next%blk_sizes(j-1, 1:outer_next%nblks_tuple(j-1))
                
                end if
             
             end if
          
          end if
       
       end do
       
        
       ! Store any Pulay n and Lagrange terms
       
          do i = 1, o_triang_size
          
             do j = 1, size(cache%indices, 1)
                
                if (size(outer_next%p_tuples) > 0) then
                   if (outer_next%p_tuples(1)%npert > 0) then

                      offset = get_triang_blks_tuple_offset(2, tot_num_pert, &
                      (/cache%nblks, outer_next%nblks_tuple(1)/), &
                      (/cache%p_inner%npert, outer_next%p_tuples(1)%npert/), &
                      blks_tuple_info, &
                      blk_sizes, &
                      (/cache%blks_triang_size, outer_next%blks_tuple_triang_size(1)/), &
                      (/cache%indices(j, :), outer_next%indices(i, :)/))
                      
                   else
                   
                      offset = j
                   
                   end if
                   
                end if
                
                if (outer_next%contrib_type == 1) then
                
                   outer_next%data_scal(offset) = contrib_pulay(j + &
                   size(cache%indices, 1) * (i - 1) + o_ctr)
                   c_ctr = c_ctr + 1
                   
                else if (outer_next%contrib_type == 3) then

                   outer_next%data_scal(offset) = contrib_pulay(j + &
                   size(cache%indices, 1) * (i - 1) + o_ctr)
                   c_ctr = c_ctr + 1
                   
                else if (outer_next%contrib_type == 4) then
                
                   outer_next%data_scal(offset) = contrib_pulay(j + &
                   size(cache%indices, 1) * (i - 1) + o_ctr)
                   
                   outer_next%data_scal(offset + size(cache%indices, 1) * o_triang_size) = &
                   contrib_pulay(j + size(cache%indices, 1) * (i - 1) + &
                   size(cache%indices, 1) * o_triang_size  + o_ctr)
                   c_ctr = c_ctr + 2
                
                
                end if
                   
                
             end do
             
          end do
          
       c_snap = c_ctr
       o_ctr = o_ctr + c_snap
       
       ! Calculate and store idempotency and SCFE terms
       if ((outer_next%contrib_type == 3) .OR. (outer_next%contrib_type == 4)) then
      
          if (outer_next%contrib_type == 3) then

             cache_offset = outer_next%blks_tuple_triang_size(1)
             
          elseif (outer_next%contrib_type == 4) then
             
             cache_offset = 2 * outer_next%blks_tuple_triang_size(1)
          
          end if
      
          allocate(d_struct_o(o_supsize_prime(k), 3))

          
          which_index_is_pid = 0
          
          do j = 1, outer_next%p_tuples(1)%npert
          
             which_index_is_pid(outer_next%p_tuples(1)%pid(j)) = j
          
          end do
          
          sstr_incr = 0
          
          if(outer_next%p_tuples(1)%npert ==0) then
          
             call derivative_superstructure(get_emptypert(), &
             (/outer_next%n_rule, outer_next%n_rule/), .TRUE., &
             (/get_emptypert(), get_emptypert(), get_emptypert()/), &
             o_supsize_prime(k), sstr_incr, d_struct_o)
             
           else
      
             call derivative_superstructure(outer_next%p_tuples(1), &
             (/outer_next%n_rule, outer_next%n_rule/), .TRUE., &
             (/get_emptypert(), get_emptypert(), get_emptypert()/), &
             o_supsize_prime(k), sstr_incr, d_struct_o)
            
          end if
      

          do i = 1, o_triang_size

             call QcMatZero(Y)
             call rsp_get_matrix_y(o_supsize_prime(k), d_struct_o, cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert, &
                  which_index_is_pid(1:cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert), size(outer_next%indices(i,:)), outer_next%indices(i,:), &
                  F, D, S, Y)
             
             ! NOTE: Rule choice very likely to give correct exclusion but
             ! have another look if something goes wrong
             call QcMatZero(Z)
             call rsp_get_matrix_z(o_supsize_prime(k), d_struct_o, &
                  (/outer_next%n_rule, outer_next%n_rule/), cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert, &
                  which_index_is_pid(1:cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert), size(outer_next%indices(i,:)), outer_next%indices(i,:), &
                  F, D, S, Z)
                  
             do j = 1, cache%blks_triang_size
           
                if (size(outer_next%p_tuples) > 0) then
                   if (outer_next%p_tuples(1)%npert > 0) then

                      offset = get_triang_blks_tuple_offset(2, tot_num_pert, &
                      (/cache%nblks, outer_next%nblks_tuple(1)/), &
                      (/cache%p_inner%npert, outer_next%p_tuples(1)%npert/), &
                      blks_tuple_info, &
                      blk_sizes, &
                      (/cache%blks_triang_size, outer_next%blks_tuple_triang_size(1)/), &
                      (/cache%indices(j, :), outer_next%indices(i, :)/))
                      
                   else
                   
                      offset = j
                   
                   end if
                   
                end if
                
!                 write(*,*) 'offset, c_snap', offset, c_snap
                
                call QcMatTraceAB(Zeta(j), Z, outer_next%data_scal(c_snap + offset))
                call QcMatTraceAB(Lambda(j), Y, outer_next%data_scal(c_snap + &
                cache%blks_triang_size*o_triang_size + offset))
                
                outer_next%data_scal(c_snap + offset) = &
                -2.0 * outer_next%data_scal(c_snap + offset)
                
                outer_next%data_scal(c_snap + &
                cache%blks_triang_size*o_triang_size + offset) = &
                -2.0 * outer_next%data_scal(c_snap + &
                cache%blks_triang_size*o_triang_size + offset)
             
             end do
                  
          end do
       
          deallocate(d_struct_o)
          
       end if
       
       deallocate(blk_sizes)
       deallocate(blks_tuple_info)
       
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
       
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    
  end subroutine

  
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
