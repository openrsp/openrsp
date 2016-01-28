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
                                    rsp_get_matrix_zeta_2014,          &
                                    rsp_get_matrix_lambda_2014,        &
                                    rsp_get_matrix_z_2014,             &
                                    rsp_get_matrix_w_2014,             &
                                    rsp_get_matrix_y_2014,             &
                                    rsp_get_matrix_zeta,          &
                                    rsp_get_matrix_lambda,        &
                                    rsp_get_matrix_z,             &
                                    rsp_get_matrix_w,             &
                                    rsp_get_matrix_y
  use rsp_perturbed_sdf, only: rsp_fds
  use rsp_property_caching, only: property_cache,                      &
                                  contrib_cache_outer,                 &
                                  contrib_cache,                       &
                                  property_cache_initialize,           &
                                  property_cache_next_element,         &
                                  property_cache_add_element,          &
                                  property_cache_already,              &
                                  property_cache_getdata,              &
                                  property_cache_allocate,             &
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
                                  contrib_cache_allocate
  
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
          
!           write(*,*) 'first, last', dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + &
!           1 + (j - 1)*np(i), dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + (j)*np(i)
          
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
                  id_outp, prop_sizes, rsp_tensor)

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
                  id_outp, prop_sizes, props)

    implicit none

    
    logical :: traverse_end
    integer :: n_props, id_outp, i, j, k
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
    type(property_cache), pointer :: contribs_cache
    type(contrib_cache_outer) :: F, D, S
    
    call empty_p_tuple(emptypert)
    emptyp_tuples = (/emptypert, emptypert/)
  
  
    ! Get all necessary F, D, S derivatives
    
    write(id_outp,*) ' '
    write(id_outp,*) 'Calculating perturbed overlap/density/Fock matrices'
    write(id_outp,*) ' '

    call cpu_time(time_start)
    call contrib_cache_allocate(contribution_cache)
        
    call rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, &
                 get_rsp_sol, get_ovl_mat, get_1el_mat, &
                 get_2el_mat, get_xc_mat, .TRUE., contribution_cache, id_outp)
                     
    deallocate(contribution_cache)
    call cpu_time(time_end)

    write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
    write(id_outp,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
    write(id_outp,*) ' '
    

    
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
    cache_next => cache_next%next
       
    ! Traverse linked list and calculate
    do while (traverse_end .eqv. .FALSE.)
       
       write(*,*) 'Calculating contribution for inner perturbation tuple'
       write(*,*) cache_next%p_inner%plab

       call rsp_energy_calculate(D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, cache_next)
       write(*,*) ' '
          
       if (cache_next%last) then
          traverse_end = .TRUE.
       end if
          
       cache_next => cache_next%next
          
    end do

    call cpu_time(time_end)

    write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
    write(id_outp,*) 'Finished calculating HF energy-type contributions'
    write(id_outp,*) ' '

    ! For each property: Recurse to identify HF energy-type contributions and 
    ! add to the property under consideration

    k = 1
    
    do i = 1, n_props
    
       do j = 1, n_freq_cfgs(i)

          call property_cache_allocate(contribs_cache)
          call contrib_cache_allocate(contribution_cache)
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
       
          write(*,*) 'Property is now', &
          props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k)))

          k = k + 1
       
       end do
       
    end do
    
    
    deallocate(contribution_cache)
    
! Notes:
! 
! The primed contributions are a subset of the nonprimed contributions
! 
! If the nonprimed recursion is done, then the primed recursion can be 
! a point-by-point test of the resulting elements upon traversal
! 
! Scheme:
!
! - Do the nonprimed recursion for identification
!
! - Traverse for calculation:
! - Check each element for primed-validity
! - If no: Calculate only Pulay n contribution (allocate accordingly)
! - If yes: Calculate both Pulay n and all Lagrange contributions (allocate accordingly)
! 
! - Do the nonprimed recursion again for retrieval: 
! - Check each nonprimed "yes" element for primed-validity
! - If no: Retrieve only Pulay n contribution
! - If yes: Retrieve both Pulay n and all Lagrange contributions
! 
! Comments:
!
! - Need to add optional "hard offset" argument for addressing proper contribution
! - Ensures sparing recursion, groups Pulay n and Pulay Lagrange contributions together
!   for increased efficiency
! - Fits well with existing framework

    
    
    
    ! For each property: Recurse to two-factor contributions and store in cache
    
    k  = 1
    
    call contrib_cache_allocate(contribution_cache)
    
    do i = 1, n_props
    
       do j = 1, n_freq_cfgs(i)
       
          write(id_outp,*) ' '
          write(id_outp,*) 'Identifying two-factor contributions'
          write(id_outp,*) ' '

          call cpu_time(time_start)
          
          call rsp_twofact_recurse(p_tuple_remove_first(p_tuples(k)), &
               kn_rule(k,:), (/p_tuple_getone(p_tuples(k), 1), emptypert/), &
               .TRUE., contribution_cache, prop_sizes(k), &
               props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
          call cpu_time(time_end)

          write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
          write(id_outp,*) 'Finished identifying two-factor contributions'
          write(id_outp,*) ' '
          
          k = k + 1
       
       end do
       
    end do
    
       
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
    cache_next => cache_next%next
       
    ! Traverse linked list and calculate
    do while (traverse_end .eqv. .FALSE.)
       
       write(*,*) 'Calculating contribution for factor 1 tuple'
       write(*,*) cache_next%p_inner%plab

       call rsp_twofact_calculate(S, D, F, get_ovl_exp, cache_next)
          
       write(*,*) ' '
          
       if (cache_next%last) then
          traverse_end = .TRUE.
       end if
          
       cache_next => cache_next%next
          
    end do

    call cpu_time(time_end)

    write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
    write(id_outp,*) 'Finished calculating two-factor contributions'
    write(id_outp,*) ' '

!     stop
    
    ! For each property: Recurse to identify two-factor contributions and 
    ! add to the property under consideration
       
    k = 1
    
    do i = 1, n_props
    
       do j = 1, n_freq_cfgs(i)

          write(id_outp,*) ' '
          write(id_outp,*) 'Assembling two-factor contributions'
          write(id_outp,*) ' '

          call cpu_time(time_start)
          
          call rsp_twofact_recurse(p_tuple_remove_first(p_tuples(k)), &
               kn_rule(k,:), (/p_tuple_getone(p_tuples(k), 1), emptypert/), &
               .FALSE., contribution_cache, prop_sizes(k), &
               props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k))))
          
          call cpu_time(time_end)

          write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
          write(id_outp,*) 'Finished assembling two-factor contributions'
          write(id_outp,*) ' '
       
          write(*,*) 'Property is now', &
          props(sum(prop_sizes(1:k)) - prop_sizes(k) + 1:sum(prop_sizes(1:k)))

          k = k + 1
          
       end do
       
    end do
  
    deallocate(contribution_cache)
  
   
    
  end subroutine
  
  
  

!   subroutine get_prop_2014(pert, kn, num_blks, blk_sizes, blk_info, F, D, S, get_rsp_sol, &
!                   get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
!                   get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, &
!                   id_outp, prop_size, prop)
! 
! !     use rsp_property_caching
!                   
!     implicit none
! 
!     type(p_tuple) :: pert, emptypert
!     type(p_tuple), dimension(2) :: emptyp_tuples
!     integer :: num_blks, id_outp, prop_size
!     integer, dimension(2) :: kn
!     integer, dimension(num_blks) :: blk_sizes
!     integer, dimension(num_blks,3) :: blk_info
!     type(sdf_2014) :: F, D, S
!     external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
!     external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
!     complex(8), dimension(*) :: prop
!     type(contrib_cache), pointer :: contribution_cache
!     type(property_cache), pointer :: contribs_cache
!     
! 
!     call empty_p_tuple(emptypert)
!     emptyp_tuples = (/emptypert, emptypert/)

    !prop = 0.0

    ! Get all necessary F, D, S derivatives as dictated by
    ! number of perturbations and kn
! 
!     write(id_outp,*) ' '
!     write(id_outp,*) 'Calculating perturbed overlap/density/Fock matrices'
!     write(id_outp,*) ' '
! 
!     call cpu_time(time_start)
!     call rsp_fds_2014(pert, kn, F, D, S, get_rsp_sol, get_ovl_mat, get_1el_mat, &
!                  get_2el_mat, get_xc_mat, id_outp)
!     call cpu_time(time_end)
! 
!     write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!     write(id_outp,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
!     write(id_outp,*) ' '

! 
! write(*,*) 'num perts', pert%npert
! 
! 
!     call property_cache_allocate(contribs_cache)
!     call contrib_cache_allocate(contribution_cache)
!     write(id_outp,*) ' '
!     write(id_outp,*) 'Calculating HF-energy type contribs'
!     write(id_outp,*) ' '
! 
!     call cpu_time(time_start)
! !     call rsp_energy_2014(pert, pert%npert, kn, 1, (/emptypert/), 0, D, get_nucpot, &
! !                     get_1el_exp, get_ovl_exp, get_2el_exp, contrib_cache, prop_size, prop)
!     call rsp_energy_recurse(pert, pert%npert, kn, 1, (/emptypert/), 0, D, get_nucpot, &
!                     get_1el_exp, get_ovl_exp, get_2el_exp, .TRUE., contribution_cache, prop_size, prop)
!     call cpu_time(time_end)
! 
!     write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!     write(id_outp,*) 'Finished calculating HF energy-type contribs'
!     write(id_outp,*) ' '
!     deallocate(contribs_cache)
!     
!     write(*,*) 'Property is now', prop(1:prop_size)

! 
!     write(*,*) ' '
!     write(*,*) 'Calculating exchange/correlation contribs'
!     write(*,*) ' '
!     call cpu_time(time_start)
!     ! CHANGE TO USE CALLBACK FUNCTIONALITY
! !     call rsp_xcave_interface(pert, kn, num_blks, blk_sizes, blk_info, D, &
! !                              get_xc_exp, prop)
!     call cpu_time(time_end)
!     write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!     write(*,*) 'Finished calculating exchange/correlation contribs'
!     write(*,*) ' '
! 
!     
!     write(*,*) 'Property is now', prop(1:prop_size)
! 
!     call property_cache_allocate(contribs_cache)
!     write(*,*) ' '
!     write(*,*) 'Calculating Pulay n type contribs'
!     write(*,*) ' '
!     call cpu_time(time_start)
!     call rsp_pulay_n_2014(pert, kn, (/emptypert, emptypert/), S, D, F, &
!                      get_ovl_exp, contribs_cache, prop_size, prop)
!     call cpu_time(time_end)
!     write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!     write(*,*) ' '
!     write(*,*) 'Finished calculating Pulay n type contribs'
!     write(*,*) ' '
!     deallocate(contribs_cache)
!     
!         write(*,*) 'Property is now', prop(1:prop_size)
! 
!     ! There are Lagrangian type contribs only when not using n + 1 rule
!     if (kn(2) < pert%npert) then
! 
!        call property_cache_allocate(contribs_cache)
!        write(*,*) ' '
!        write(*,*) 'Calculating Pulay Lagrangian type contribs'
!        write(*,*) ' '
!        call cpu_time(time_start)
!        call rsp_pulay_lag_2014(p_tuple_remove_first(pert), kn, &
!                           (/p_tuple_getone(pert,1), emptypert/), &
!                           S, D, F, get_ovl_exp, contribs_cache, prop_size, prop)
!        call cpu_time(time_end)
!        write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!        write(*,*) ' '
!        write(*,*) 'Finished calculating Pulay Lagrangian type contribs' 
!        write(*,*) ' '
!    
!        deallocate(contribs_cache)
!    
!        write(*,*) 'Property is now', prop(1:prop_size)
!    
!        call property_cache_allocate(contribs_cache)
!        write(*,*) ' '
!        write(*,*) 'Calculating idempotency Lagrangian type contribs'
!        write(*,*) ' '
!        call cpu_time(time_start)
!        ! MaR: Unchanged by introduction of callback functionality
!        call rsp_idem_lag_2014(p_tuple_remove_first(pert), kn, &
!                          (/p_tuple_getone(pert,1), emptypert/), &
!                          S, D, F, contribs_cache, prop_size, prop)
!        call cpu_time(time_end)
!        write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!        write(*,*) ' '
!        write(*,*) 'Finished calculating idempotency Lagrangian type contribs'
!        write(*,*) ' '
!    
!        deallocate(contribs_cache)
!    
!        write(*,*) 'Property is now', prop(1:prop_size)
!        
!        call property_cache_allocate(contribs_cache)
!        write(*,*) ' '
!        write(*,*) 'Calculating SCF Lagrangian type contribs'
!        write(*,*) ' '
!        call cpu_time(time_start)
!        ! MaR: Unchanged by introduction of callback functionality
!        call rsp_scfe_lag_2014(p_tuple_remove_first(pert), kn, &
!                          (/p_tuple_getone(pert,1), emptypert/), &
!                          S, D, F, contribs_cache, prop_size, prop)
!        call cpu_time(time_end)
!    
!    
!        write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
!        write(*,*) ' '
!        write(*,*) 'Finished calculating SCF Lagrangian type contribs'
!        write(*,*) ' '
! 
!        deallocate(contribs_cache)
! 
!     write(*,*) 'Property is now', prop(1:prop_size)       
!        
!     end if
! 
!   end subroutine
  
  
  ! END NEW 2014
  


! BEGIN NEW 2014

  recursive subroutine rsp_energy_2014(pert, total_num_perturbations, kn, num_p_tuples, &
                                p_tuples, density_order, D, get_nucpot, get_1el_exp, &
                                get_t_exp, get_2el_exp, cache, p_size, prop)

    implicit none

    logical :: e_knskip
    type(p_tuple) :: pert, p_rf1, p_rf2, p_rf3
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new, p_new1, p_stord
    type(p_tuple), dimension(num_p_tuples + 1) :: p_new3
    external :: get_nucpot, get_1el_exp, get_t_exp, get_2el_exp
    type(SDF_2014) :: D
    type(property_cache) :: cache
    complex(8), dimension(*) :: prop
    integer :: num_blks_full, p_size
    integer, allocatable, dimension(:,:) :: blk_info_full


    
    
!     allocate(blk_info_full(num_blks_full, 3))
        
!     num_blks_full = get_num_blks(pert)
    
!     blk_info_full = get_blk_info(num_blks_full, pert)
!     p_size = get_triangulated_size(num_blks_full, blk_info_full)
    
!     deallocate(blk_info_full)
        
    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'

    if (p_tuples(1)%npert == 0) then

       p_rf1 = p_tuple_remove_first(pert)
       p_new1 = (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/)
    
       call rsp_energy_2014(p_rf1, total_num_perturbations, &
       kn, num_p_tuples, p_new1, &
       density_order, D, get_nucpot, get_1el_exp, &
       get_t_exp, get_2el_exp, cache, p_size, prop)

    else

       p_rf1 = p_tuple_remove_first(pert)
       p_new1 = (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
       p_tuples(2:size(p_tuples))/)

       call rsp_energy_2014(p_rf1, total_num_perturbations,  &
       kn, num_p_tuples, p_new1, density_order, D, get_nucpot, get_1el_exp, &
       get_t_exp, get_2el_exp, cache, p_size, prop)

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

          p_rf2 = p_tuple_remove_first(pert)
          
          call rsp_energy_2014(p_rf2, total_num_perturbations, &
          kn, num_p_tuples, t_new, density_order + 1, D, get_nucpot, get_1el_exp, &
          get_t_exp, get_2el_exp, cache, p_size, prop)

       end do

       ! MaR: Since we are only calculating Hartree-Fock type energy terms here,
       ! we don't need to go beyond to perturbed contraction density matrices
       ! (but that is in general needed for XC contributions)
       if (num_p_tuples < 3) then

          ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
          ! a(nother) pert D contraction)

          p_rf3 = p_tuple_remove_first(pert)
          p_new3 = (/p_tuples(:), p_tuple_getone(pert, 1)/)
          
          call rsp_energy_2014(p_rf3, total_num_perturbations, &
          kn, num_p_tuples + 1, p_new3, &
          density_order + 1, D, get_nucpot, get_1el_exp, &
          get_t_exp, get_2el_exp, cache, p_size, prop)

       end if

    ! At the final recursion level: Calculate the contribution (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allows it)

    else

       e_knskip = .FALSE.


!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)

!        write(*,*) 'Getting energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

!              write(*,*) 'D ', p_tuples(i)%pid

             if(kn_skip(p_tuples(i)%npert, p_tuples(i)%pid, kn) .EQV. .TRUE.) then

                e_knskip = .TRUE.

             end if
          
          elseif (i == 1) then

!              write(*,*) 'E ', p_tuples(i)%pid

          end if

       end do


       if (e_knskip .EQV. .FALSE.) then

          p_stord = p_tuples_standardorder(num_p_tuples, p_tuples)
       
!          open(unit=257, file='totterms', status='old', action='write', &
!               position='append') 
!          write(257,*) 'T'
!          close(257)
          
!           write(*,*) 'Evaluating property_cache_already'

          if (property_cache_already(cache, num_p_tuples, p_stord)) then

!             open(unit=257, file='cachehit', status='old', action='write', &
!                  position='append') 
!             write(257,*) 'T'
!             close(257)

!              write(*,*) 'Getting values from cache'

             ! NOTE (MaR): EVERYTHING MUST BE STANDARD ORDER IN 
             ! THIS CALL (LIKE property_cache_getdata ASSUMES)
             call property_cache_getdata(cache, num_p_tuples, p_stord, p_size, prop)

!              write(*,*) ' '
       
          else


       write(*,*) 'Calculating energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

             write(*,*) 'D ', p_tuples(i)%pid

!              if(kn_skip(p_tuples(i)%npert, p_tuples(i)%pid, kn) .EQV. .TRUE.) then
! 
!                 e_knskip = .TRUE.
! 
!              end if
          
          elseif (i == 1) then

             write(*,*) 'E ', p_tuples(i)%pid

          end if

       end do

        write(*,*) 'Property before get energy is', prop(1:p_size)   
             call get_energy_2014(num_p_tuples, total_num_perturbations, & 
                  p_stord, density_order, D, get_nucpot, get_1el_exp, &
                  get_t_exp, get_2el_exp, cache, prop)

                  write(*,*) 'Calculated energy contribution'
                  write(*,*) ' '
        write(*,*) 'Property after get energy is', prop(1:p_size)   
                  
                  
          end if

       else

!           write(*,*) 'Energy contribution was k-n skipped'
!           write(*,*) ' '

       end if

    end if

  end subroutine



  recursive subroutine rsp_energy_recurse(pert, total_num_perturbations, kn, num_p_tuples, &
                                  p_tuples, density_order, D, get_nucpot, get_1el_exp, &
                                  get_t_exp, get_2el_exp, dryrun, cache, p_size, prop)

    implicit none

    logical :: e_knskip, dryrun, traverse_end
    type(p_tuple) :: pert
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, id_outp
    integer, optional :: p_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(contrib_cache_outer) :: D
    type(contrib_cache), target :: cache
    type(contrib_cache), pointer :: cache_next
    complex(8), dimension(*), optional :: prop
    external :: get_nucpot, get_1el_exp, get_t_exp, get_2el_exp

!     write(*,*) 'Dryrun?', dryrun
!     write(*,*) 'for pert', pert%npert
    

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
       
!        write(*,*) 'about to check already'

          if (contrib_cache_already(cache, num_p_tuples, p_tuples)) then
          
!           write(*,*) 'was already'

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
                ! call contrib_cache_getdata(cache, num_p_tuples, &
                !    p_tuples_standardorder(num_p_tuples, p_tuples), p_size, .FALSE., prop=prop)

             end if

          else
!           write(*,*) 'was not already'
             

             if (dryrun) then
             
!              write (*,*) 'not found, adding element'

                call contrib_cache_add_element(cache, num_p_tuples, p_tuples)

!              else
! 
!                 write(*,*) 'ERROR: Contribution should be in cache but was not found', kn
!                 
!                     do i = 1, num_p_tuples
! 
!     
!     
!     
!     write(*,*) 'p tuples', i
!     write(*,*) p_tuples(i)%npert
!     write(*,*) p_tuples(i)%pid
!    write(*,*) ' '
!     write(*,*) ' '
! 
!     
!                     end do

             end if

          end if

       end if

    end if



  end subroutine

    subroutine get_energy_2014(num_p_tuples, total_num_perturbations, &
                        p_tuples, density_order, D, get_nucpot, get_1el_exp, &
                        get_ovl_exp, get_2el_exp, cache, prop)

    implicit none

    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(p_tuple) :: merged_p_tuple, t_matrix_bra, t_matrix_ket, t_matrix_newpid
    type(SDF_2014) :: D
    type(property_cache) :: cache
    type(QcMat) :: D_unp
    type(QcMat), allocatable, dimension(:) :: dens_tuple
!    type(rsp_field), allocatable, dimension(:) :: nucpot_pert
    external :: get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp
    integer :: i, j, k, m, n, num_p_tuples, total_num_perturbations, density_order, &
                offset, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks, npert_ext
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged, pert_ext
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter, &
                                              pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, pidoutersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, contrib, prop_forcache
    complex(8), dimension(*) :: prop

!    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
!    ncouter = nc_only(total_num_perturbations, total_num_perturbations - &
!              p_tuples(1)%npert, num_p_tuples - 1, &
!              p_tuples(2:num_p_tuples), ncarray)
!    ncinner = nc_only(total_num_perturbations, p_tuples(1)%npert, 1, &
!                      p_tuples(1), ncarray)

    allocate(dens_tuple(num_p_tuples))
!    allocate(nucpot_pert(p_tuples(1)%npert))
 !   allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%npert))
 !   allocate(ncinnersmall(p_tuples(1)%npert))
 !   allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%npert))

 !   ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
 !                  p_tuples(1)%npert, num_p_tuples - 1, &
 !                  p_tuples(2:num_p_tuples), ncarray)
 !   ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%npert, &
 !                  1, p_tuples(1), ncarray)
 !   pidoutersmall = get_pidoutersmall(total_num_perturbations - &
 !                   p_tuples(1)%npert, num_p_tuples - 1, &
 !                   p_tuples(2:num_p_tuples))

!    call QcMatInit(D_unp)

    call p_tuple_to_external_tuple(p_tuples(1), npert_ext, pert_ext)
 
    call p1_cloneto_p2(p_tuples(1), t_matrix_newpid)
    t_matrix_newpid%pid = (/(i, i = 1, t_matrix_newpid%npert)/)

    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(nfields(num_p_tuples))
    allocate(nblks_tuple(num_p_tuples))

    do i = 1, num_p_tuples

       nfields(i) = p_tuples(i)%npert
       nblks_tuple(i) = get_num_blks(p_tuples(i))

    end do

    allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(num_p_tuples))
    allocate(blk_sizes(num_p_tuples, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, num_p_tuples

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))

       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    outer_indices_size = product(blks_tuple_triang_size(2:num_p_tuples))

    if (p_tuples(1)%npert == 0) then

       inner_indices_size = 1

    else

       inner_indices_size = blks_tuple_triang_size(1)

    end if
    
    allocate(tmp(inner_indices_size))
    allocate(contrib(inner_indices_size))
    allocate(prop_forcache(inner_indices_size*outer_indices_size))

    prop_forcache = 0.0
    contrib = 0.0

!    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
!                      p_tuples(1)%npert, pidoutersmall, &
!                      ncarray, ncoutersmall, o_whichpert)

    call sdf_getdata_s_2014(D, get_emptypert(), (/1/), D_unp)


    if (total_num_perturbations > p_tuples(1)%npert) then

       allocate(outer_indices(outer_indices_size,total_num_perturbations - &
                p_tuples(1)%npert))
       allocate(inner_indices(inner_indices_size,p_tuples(1)%npert))

       k = 1
    
       do i = 2, num_p_tuples
          do j = 1, p_tuples(i)%npert
    
             o_wh_forave(p_tuples(i)%pid(j)) = k
             k = k + 1
    
          end do
       end do
    
       k = 1
   
!       do i = 2, num_p_tuples
!          do j = 1, p_tuples(i)%npert
!   
!             ncoutersmall(k) =  p_tuples(i)%pdim(j)
!             k = k + 1
!   
!          end do
!       end do
   
       do i = 2, num_p_tuples
   
          call QcMatInit(dens_tuple(i), D_unp)
   
       end do
   
       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)
   
   
       if (p_tuples(1)%npert > 0) then
   
          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)
   
       end if
   
       do i = 1, size(outer_indices, 1)
   
          dtup_ind = 0
   
          do j = 2, num_p_tuples
   
             call sdf_getdata_s_2014(D, p_tuples(j), outer_indices(i, &
                  dtup_ind+1:dtup_ind + p_tuples(j)%npert), dens_tuple(j))
   
             dtup_ind = dtup_ind + p_tuples(j)%npert
   
          end do
   
          tmp = 0.0
          contrib = 0.0
   
          if (num_p_tuples == 2) then
          
             call get_1el_exp(npert_ext, pert_ext, 1, (/dens_tuple(2)/), &
                              size(contrib), contrib)
          
             
!              call rsp_oneave(p_tuples(1)%npert, p_tuples(1)%plab, &
!                             (/ (1, j = 1, p_tuples(1)%npert) /), &
!                             p_tuples(1)%pdim, dens_tuple(2), &
!                             nblks_tuple(1),  blks_tuple_info(1, 1:nblks_tuple(1), :), &
!                             blk_sizes(1, 1:nblks_tuple(1)), inner_indices_size, contrib)
   
          end if
          
          write(*,*) ' '
          write(*,*) 'oneave contrib',  contrib
   
          tmp = tmp + contrib
          contrib = 0.0
   
          if (num_p_tuples == 2) then

             t_matrix_bra = get_emptypert()
             t_matrix_ket = get_emptypert()

             call rsp_ovlave_t_matrix_2014(t_matrix_newpid%npert, t_matrix_newpid, &
                                      t_matrix_bra, t_matrix_ket, &
                                      dens_tuple(2), get_ovl_exp, inner_indices_size, contrib)
   
          end if
   
   write(*,*) ' '
          write(*,*) 'ovlave t contrib', contrib
   
   
          tmp = tmp - contrib
          contrib = 0.0
   
          if (num_p_tuples == 2) then
   
             call get_2el_exp(npert_ext, pert_ext, 1, (/1/), (/dens_tuple(2)/), &
                              (/1/), (/D_unp/), size(contrib), contrib)
   
!              call rsp_twoave(p_tuples(1)%npert, p_tuples(1)%plab, &
!                              (/ (1, j = 1, p_tuples(1)%npert) /), &
!                              p_tuples(1)%pdim, dens_tuple(2), &
!                              D_unp, inner_indices_size, contrib)
    
          elseif (num_p_tuples == 3) then
          
             call get_2el_exp(npert_ext, pert_ext, 1, (/1/), (/dens_tuple(2)/), &
                              (/1/), (/dens_tuple(3)/), size(contrib), contrib)
          
!              call rsp_twoave(p_tuples(1)%npert, p_tuples(1)%plab, &
!                              (/ (1, j = 1, p_tuples(1)%npert) /), &
!                              p_tuples(1)%pdim, dens_tuple(2), dens_tuple(3), &
!                              inner_indices_size, contrib)
    
          end if

          write(*,*) ' '
          write(*,*) 'twoave contrib', contrib
   
              
          tmp = tmp + contrib
    
          if (p_tuples(1)%npert > 0) then
    
             do j = 1, size(inner_indices, 1)
    
                offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, &
                         nblks_tuple, (/ (p_tuples(k)%npert, k = 1, num_p_tuples) /), &
                         blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/inner_indices(j, :), outer_indices(i, :) /)) 
!     write(*,*) 'indices', (/inner_indices(j, :), outer_indices(i, :) /)
!     write(*,*) 'offset in cache', offset, 'is j', j
    
                prop_forcache(offset) = prop_forcache(offset) + tmp(j)
    
             end do
    
          else
    
             offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations,  &
                      nblks_tuple(2:num_p_tuples), &
                      (/ (p_tuples(k)%npert, k = 2, num_p_tuples) /), &
                      blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
                      blks_tuple_triang_size(2:num_p_tuples), (/outer_indices(i, :) /)) 
    
             prop_forcache(offset) = prop_forcache(offset) + tmp(1)
    
          end if
    
       end do
    
       if (p_tuples(1)%npert > 0) then
    
          call p1_cloneto_p2(p_tuples(1), merged_p_tuple)
    
          do i = 2, num_p_tuples
    
             ! MaR: This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))
   
          end do
   
       else
   
          call p1_cloneto_p2(p_tuples(2), merged_p_tuple)
   
          do i = 3, num_p_tuples
   
             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))
   
          end do
   
       end if
   
       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
   
       ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
       ! PIDS ARE IN STANDARD ORDER? FIND OUT
   
       k = 1
       do i = 1, num_p_tuples
          do j = 1, p_tuples(i)%npert
             pids_current_contribution(k) = p_tuples(i)%pid(j)
             k = k + 1
          end do
       end do
    
       merged_nblks = get_num_blks(merged_p_tuple)
       
       allocate(merged_blk_info(1, merged_nblks, 3))
       
       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
    
       allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))
    
       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
                                      merged_triang_size, triang_indices_pr)
    
       do i = 1, size(triang_indices_pr, 1)
    
          pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                      (/sum(nfields)/), &
                      (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                      (/triang_indices_pr(i, :) /))
    
          do j = 1, total_num_perturbations
    
             translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))
    
          end do
    
          if (p_tuples(1)%npert > 0) then
    
             ec_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                         total_num_perturbations, nblks_tuple, &
                         nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/ translated_index(:) /))
    
          else
    
             ec_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
                         total_num_perturbations, nblks_tuple(2:num_p_tuples), &
                         nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
                         blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
                         (/ translated_index(:) /))
    
          end if
   
!    write(*,*) 'pr indices are', triang_indices_pr(i,:)
!       write(*,*) 'translated index is', translated_index(:)
!    write(*,*) 'offset in prop', pr_offset, 'is cache offset', ec_offset
   
          prop(pr_offset) = prop(pr_offset) + prop_forcache(ec_offset)
   
       end do
   
       call property_cache_add_element(cache, num_p_tuples, p_tuples, &
                                       inner_indices_size * outer_indices_size, prop_forcache)   
   
       deallocate(triang_indices_pr)
       deallocate(outer_indices)
       deallocate(inner_indices)

    else

       ! MR: THIS IS THE CASE 'ALL INDICES ARE INNER INDICES'
       ! THIS TERM OCCURS ONLY ONCE AND DOES NOT NEED TO BE CACHED

!       do i = 1, p_tuples(1)%npert
!
!          nucpot_pert(i) = rsp_field(p_tuples(1)%plab(i), p_tuples(1)%freq(i), 1, &
!                                     p_tuples(1)%pdim(i))
!
!       end do

       tmp = 0.0
       contrib = 0.0

       call get_nucpot(p_tuples(1)%npert, p_tuples(1)%pdim, &
                       (/ (1, j = 1, p_tuples(1)%npert) /), &
                       p_tuples(1)%plab, contrib)
       
       
!       call rsp_nucpot(nucpot_pert, contrib) 

       tmp = tmp + contrib
       contrib = 0.0

       write(*,*) 'Getting a 1el contribution with npert pert ord ='
       write(*,*) npert_ext
       write(*,*) pert_ext
       
       
       call get_1el_exp(npert_ext, pert_ext, 1, D_unp, size(contrib), contrib)
       
!        call rsp_oneave(p_tuples(1)%npert, p_tuples(1)%plab, &
!                        (/ (1, j = 1, p_tuples(1)%npert) /), p_tuples(1)%pdim, &
!                        D_unp, nblks_tuple(1),  blks_tuple_info(1, 1:nblks_tuple(1), :), &
!                        blk_sizes(1, 1:nblks_tuple(1)), contrib)

                                 write(*,*) ' '
          write(*,*) 'oneave contrib', contrib
                       
       tmp = tmp + contrib
       contrib = 0.0

       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()

       call rsp_ovlave_t_matrix_2014(t_matrix_newpid%npert, t_matrix_newpid, &
                                t_matrix_bra, t_matrix_ket, &
                                D_unp, get_ovl_exp, inner_indices_size, contrib)

                                          write(*,*) ' '
          write(*,*) 'ovlave t mat contrib', contrib
                                
       tmp = tmp - contrib
       contrib = 0.0

       write(*,*) 'pert definition', npert_ext
       write(*,*) pert_ext

       call get_2el_exp(npert_ext, pert_ext, 1, (/1/), (/D_unp/), &
                        (/1/), (/D_unp/), size(contrib), contrib)
       
!        call rsp_twoave(p_tuples(1)%npert, p_tuples(1)%plab, &
!                        (/ (1, j = 1, p_tuples(1)%npert) /), p_tuples(1)%pdim, &
!                        D_unp, D_unp, contrib)

                                 write(*,*) ' '
          write(*,*) 'twoave contrib', contrib
                       
       tmp = tmp + 0.5*(contrib)

       do i = 1, inner_indices_size
           prop(i) =  prop(i) + tmp(i)
       end do

       call p1_cloneto_p2(p_tuples(1), merged_p_tuple)

       do i = 2, num_p_tuples

          ! This can be problematic - consider rewriting merge_p_tuple as subroutine
          merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

       end do

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
       merged_nblks = get_num_blks(merged_p_tuple)
       allocate(merged_blk_info(1, merged_nblks, 3))
       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

    end if

!  write(*,*) 'energy contribution'
!  call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!  (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
!  merged_blk_info, prop_forcache)

    call QcMatDst(D_unp)

    do i = 2, num_p_tuples
   
       call QcMatDst(dens_tuple(i))
   
    end do

    deallocate(pert_ext)
    
    deallocate(merged_blk_info)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
!    deallocate(nucpot_pert)
    deallocate(dens_tuple)
!    deallocate(ncoutersmall)
!    deallocate(ncinnersmall)
!    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(tmp)
    deallocate(contrib)
    deallocate(prop_forcache)

  end subroutine
  
  
  subroutine rsp_energy_calculate(D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, cache)

    implicit none

    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n 
    integer :: id_outp, c1_ctr, c2_ctr, lhs_ctr_1, lhs_ctr_2, rhs_ctr_2
    integer :: total_outer_size_1, total_outer_size_2
    integer :: num_0, num_1, num_pert
    character(30) :: mat_str, fmt_str
    type(contrib_cache_outer) :: D
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: outer_next
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, LHS_dmat_2, RHS_dmat_2
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    integer, allocatable, dimension(:,:) :: outer_contract_sizes_2
    complex(8), allocatable, dimension(:) :: contrib_0, contrib_1, contrib_2
    external :: get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp
    
    
    ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)
    
    outer_next => cache%contribs_outer
    
    write(*,*) 'num outer', cache%num_outer

    
    allocate(outer_contract_sizes_1(cache%num_outer))
    allocate(outer_contract_sizes_2(cache%num_outer,2))
        
   
    ! Traversal: Find number of density matrices for contraction for nuc-nuc, 1-el, 2-el cases
    
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
       
    total_outer_size_1 = 0
    total_outer_size_2 = 0
    num_0 = 0
    num_1 = 0
        
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
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
   
       write(*,*) 'Outer contribution:'!, outer_next%num_dmat, outer_next%dummy_entry
    
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'D', outer_next%p_tuples(i)%pid
       
       end do
    
       
    
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do
 
    ! Make collapsed contraction sizes array for 1-el call
 
    allocate(outer_contract_sizes_1_coll(num_1))
    
    k = 1 
     do i = 1, cache%num_outer
        if (outer_contract_sizes_1(i) > 0) then
           outer_contract_sizes_1_coll(k) = outer_contract_sizes_1(i)
          k = k + 1
        end if
    end do
    
!     write(*,*) 'outer 1', outer_contract_sizes_1
!     write(*,*) 'outer 2 1', outer_contract_sizes_2(:,1)
!     write(*,*) 'outer 2 2', outer_contract_sizes_2(:,2)
    

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
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
       
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
             call contrib_cache_getdata_outer(D, 1, outer_next%p_tuples, .FALSE., &
                  1, ind_len=outer_next%p_tuples(1)%npert, &
                  ind_unsorted=outer_next%indices(m, 1:outer_next%p_tuples(1)%npert), &
                  mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
          end do
          
          do n = 1, outer_contract_sizes_2(k, 2) 
             call contrib_cache_getdata_outer(D, 1, outer_next%p_tuples, .FALSE., &
                  1, ind_len=outer_next%p_tuples(1)%npert + 1, &
                  ind_unsorted=outer_next%indices(n, outer_next%p_tuples(1)%npert + 1: &
                  outer_next%p_tuples(1)%npert + outer_next%p_tuples(2)%npert), &
                  mat_sing=RHS_dmat_2(rhs_ctr_2 + n  - 1))
          end do          
       
       
       end if
   
       if (outer_next%next%dummy_entry) then
          traverse_end = .TRUE.
       end if

       lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
       lhs_ctr_2 = lhs_ctr_2 + outer_contract_sizes_2(k, 1)
       rhs_ctr_2 = rhs_ctr_2 + outer_contract_sizes_2(k, 2)
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    
!     write(*,*) 'cache%blks_triang_size', cache%blks_triang_size
!     write(*,*) 'total_outer_size_1', total_outer_size_1
!     write(*,*) 'total_outer_size_2', total_outer_size_2
    
    
    allocate(contrib_0(cache%blks_triang_size))
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    allocate(contrib_2(cache%blks_triang_size*total_outer_size_2))
    
    
!     do i = 1, size(LHS_dmat_1)
!     
!     
!           if (i < 10) then
!           
!              fmt_str = "(A11, I1)"
!           
!           else if (i < 100) then
!           
!              fmt_str = "(A11, I2)"
!           
!           else
!           
!              fmt_str = "(A11, I3)"
!           
!           end if
!           
! !           write(mat_str, fmt_str) 'LHS_dmaE_1_', i
!           
! !           write(*,*) 'i', i
! !           write(*,*) 'fname:', mat_str
!           
!           
!           
!           j = QcMatWrite_f(LHS_dmat_1(i), trim(mat_str), ASCII_VIEW)
!     
!     end do
    
    
    ! Calculate contributions
    
    ! Calculate nuclear-nuclear repulsion contribution
    if (num_0 > 0) then
    
    
       write(*,*) 'nucpot'
    
       contrib_0 = 0.0
       call get_nucpot(num_pert, pert_ext, size(contrib_0), contrib_0)
       
!        write(*,*) 'nucpot contribution: ', contrib_0
    
    end if
    
    ! Calculate one-electron contributions
    if (num_1 > 0) then
    
       write(*,*) '1-el'
       contrib_1 = 0.0
       call get_1el_exp(num_pert, pert_ext, total_outer_size_1, &
                        LHS_dmat_1, size(contrib_1), contrib_1)
      
       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()
      
!        call rsp_ovlave_t_matrix_2014(get_ovl_exp, cache%p_inner, cache%p_inner%npert, &
!                                 t_matrix_bra, t_matrix_ket, outer_contract_sizes_1_coll, &
!                                 LHS_dmat_1, size(contrib_1), contrib_1)
    
    write(*,*) '1-el contribution: ', contrib_1
    
    end if
    
!     write(*,*) '2-el'
    
    ! Calculate two-electron contributions
    contrib_2 = 0.0
    call get_2el_exp(num_pert, pert_ext, cache%num_outer, outer_contract_sizes_2(:, 1), LHS_dmat_2, & 
                     outer_contract_sizes_2(:, 2), RHS_dmat_2, size(contrib_2), contrib_2)
                       
    
!     write(*,*) '2-el contribution: ', contrib_2
    
    ! Add nuc-nuc, 1-el and two-el contributions together (put in contrib_2)
    
    ! Find out about "length of contributions" attribute in cache datatype: Seems like %contrib_size
    ! must be cache(inner)%blks_triang_size * product(cache_outer%blks_tuple_triang_size) and it is
    ! this attribute that should specify the contribution size: However, when outer is used as F, D, S
    ! storage, the inner size is not applicable and it should be the outer size only: Each situation
    ! must be handled appropriately
    
    
    ! Traversal: Add nuc-nuc, 1-el and two-el contributions together (put in contrib_2)
    
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
      
    k = 1
    
    c1_ctr = 1
    c2_ctr = 1
    
    
    do while (traverse_end .EQV. .FALSE.)
  
       ! Nuc-nuc, one-el and two-el contribution
       if (outer_next%num_dmat == 0) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size))

          outer_next%data_scal = contrib_0
          outer_next%data_scal = outer_next%data_scal + contrib_1(c1_ctr:c1_ctr + cache%blks_triang_size - 1)
          outer_next%data_scal = outer_next%data_scal + contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size - 1)

          c1_ctr = c1_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
          c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)

       ! One-el and two-el contribution
       else if (outer_next%num_dmat == 1) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)))

          outer_next%data_scal = contrib_1(c1_ctr:c1_ctr + cache%blks_triang_size * &
                            outer_contract_sizes_2(k, 1) - 1)
          
          outer_next%data_scal = outer_next%data_scal + contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size * &
                            outer_contract_sizes_2(k, 1) - 1)

          c1_ctr = c1_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
          c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1)
                   
       ! Only two-electron contribution
       else if (outer_next%num_dmat == 2) then
       
          allocate(outer_next%data_scal(cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)))
                   
          outer_next%data_scal = contrib_2(c2_ctr:c2_ctr + cache%blks_triang_size * &
                            outer_contract_sizes_2(k, 1) * outer_contract_sizes_2(k, 2) - 1)
       
          c2_ctr = c2_ctr + cache%blks_triang_size * outer_contract_sizes_2(k, 1) * &
                   outer_contract_sizes_2(k, 2)
       end if
   
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do

        
    deallocate(outer_contract_sizes_2)
    
  end subroutine

  
   recursive subroutine rsp_twofact_recurse(pert, kn, p12, dryrun, cache, p_size, prop)

    implicit none

    logical :: dryrun
    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: curr_outer
    integer ::  i
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
                
                call contrib_cache_add_element(cache, 2, p12, n_rule=kn(2))
!                 write(*,*) 'added element'
                call contrib_cache_cycle_outer(cache, 2, p12, curr_outer, n_rule=kn(2))
!                 write(*,*) 'cycled'
                curr_outer%contrib_type = 1
             
             else
             
                write(*,*) 'ERROR: Expected to find Pulay contribution but it was not present'
                write(*,*) 'S', p12(1)%pid
                write(*,*) 'W', p12(2)%pid
             
             end if

          end if

       end if
       
       
       if ((kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) .AND. (p12(1)%npert > 0)) then

          if (contrib_cache_already(cache, 2, p12, n_rule=kn(2))) then

             if (.NOT.(dryrun)) then
             
             ! FIXME: CHANGE THIS TO WORK WITH CONTRIBUTION TYPE
             
                write(*,*) 'Retrieving Lagrange contributions:'
                write(*,*) 'A', p12(1)%pid
                write(*,*) 'B', p12(2)%pid

                
!                 write(*,*) 'Pulay lagrange hard offset', hard_offset
                ! Pulay Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'Idempotency lagrange hard offset', hard_offset
                ! Idempotency Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
                hard_offset = hard_offset + block_size
                
!                 write(*,*) 'SCFE lagrange hard offset', hard_offset
                ! SCFE Lagrange contribution
                call contrib_cache_getdata(cache, 2, p12, p_size, 0, hard_offset=hard_offset, &
                scal=prop, n_rule=kn(2))
            
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
    integer :: cache_offset, i, j, k, m, n, c_ctr, c_snap, lagrange_max_n
    integer :: id_outp, i_supsize, o_triang_size
    integer :: ctr_lagrange, ctr_pulay_n, o_ctr, size_lagrange, size_pulay_n
    integer :: sstr_incr
    integer, dimension(0) :: nof
    integer, allocatable, dimension(:) :: o_supsize, o_supsize_prime, o_size
    integer, allocatable, dimension(:) :: which_index_is_pid
    complex(8), allocatable, dimension(:) :: contrib_pulay
    character(30) :: mat_str, fmt_str
    type(p_tuple) :: p_inner
    type(p_tuple), allocatable, dimension(:,:) :: d_struct_inner, d_struct_o, d_struct_o_prime
    type(contrib_cache_outer) :: S, D, F
    type(contrib_cache) :: cache
    type(contrib_cache_outer), pointer :: outer_next
    type(QcMat), allocatable, dimension(:) :: Lambda, Zeta, W
    type(QcMat) :: Y, Z
           
    integer, allocatable, dimension(:) :: pert_ext
    
    external :: get_ovl_exp
    
   ! Assume indices for inner, outer blocks are calculated earlier during the recursion
    
    ! Unsure about npert argument
    call p_tuple_to_external_tuple(cache%p_inner, cache%p_inner%npert, pert_ext)

    
    
    outer_next => cache%contribs_outer

    
    write(*,*) 'num outer', cache%num_outer

    i_supsize = 3**cache%p_inner%npert
    
    
    allocate(o_supsize(cache%num_outer))
    allocate(o_supsize_prime(cache%num_outer))
    allocate(o_size(cache%num_outer))

!     write(*,*) 'n rule', outer_next%n_rule
    
    
    ! Prepare matrices for inner tuple
    
    i_supsize = derivative_superstructure_getsize(p_tuple_remove_first(cache%p_inner), &
                (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

!     write(*,*) 'i supsize', i_supsize
                
    allocate(d_struct_inner(i_supsize, 3))

    sstr_incr = 0
    
    call derivative_superstructure(p_tuple_remove_first(cache%p_inner), &
          (/outer_next%n_rule, outer_next%n_rule/), .FALSE., &
          (/get_emptypert(), get_emptypert(), get_emptypert()/), &
          i_supsize, sstr_incr, d_struct_inner)  
          
    sstr_incr = 0
          
!     write(*,*) 'a'
          
    
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
    
!     write(*,*) 'b'
    
    traverse_end = .FALSE.
    
    
    do while (outer_next%dummy_entry .eqv. .FALSE.)
        outer_next => outer_next%next
    end do
    outer_next => outer_next%next
        
    size_pulay_n = 0
    size_lagrange = 0
    
    k = 1
    
!     write(*,*) 'c'
    
    do while (traverse_end .EQV. .FALSE.)
  
!         write(*,*) 'traversing'
  
      
       if (outer_next%p_tuples(1)%npert == 0) then
       
          o_triang_size = 1
       
       else
       
          o_triang_size = outer_next%blks_tuple_triang_size(1)
       
       end if
  
  
       write(*,*) 'Outer contribution, type', outer_next%contrib_type
!        write(*,*) 'triang size', o_triang_size
    
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'B', outer_next%p_tuples(i)%plab
       
       end do
    
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
!        write(*,*) 'b'
       
       if (outer_next%contrib_type >= 3) then
       
          size_lagrange = size_lagrange + 3 * o_triang_size
          
          any_lagrange = .TRUE.
          
          if (outer_next%n_rule > lagrange_max_n) then
          
             lagrange_max_n = outer_next%n_rule            
          
          end if
       
       end if
       
!        write(*,*) 'c'
!        write(*,*) 'allocation:', o_size(k) * cache%blks_triang_size
       
!        write(*,*) 'c2'
      
       allocate(outer_next%data_scal(o_size(k) * cache%blks_triang_size))
       
       outer_next%data_scal = 0.0
         
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
     
!        write(*,*) 'd'
    end do
    
    do while (outer_next%dummy_entry .eqv. .FALSE.)
        outer_next => outer_next%next
    end do
    outer_next => outer_next%next
!     write(*,*) 'e'
    
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
    
    
!           write(*,*) 'i', i, ' of', cache%blks_triang_size
!           write(*,*) 'supsize', i_supsize
       
!           write(*,*) 'indices size', size(cache%indices(i,:))
!           write(*,*) 'lagrange max n', lagrange_max_n
       
          ! CONTINUE HERE: CREATE INNER INDICES WHEN MAKING CACHE ENTRY (IT IS MISSING NOW)
          call QcMatInit(Lambda(i))
          call QcMatInit(Zeta(i))
       
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
    
    
!     write(*,*) 'f'
    
    ! Traversal: Make W matrices and store
    
    traverse_end = .FALSE.
    
    do while (outer_next%dummy_entry .eqv. .FALSE.)
        outer_next => outer_next%next
    end do
    outer_next => outer_next%next

!     write(*,*) 'g'
    
    ctr_pulay_n = 0
    ctr_lagrange = 0
    o_ctr = 1
    
    allocate(W(size_pulay_n + size_lagrange/3))
    
    do i = 1, size(w)
    
       call QcMatInit(W(i))
    
    end do
    
    allocate(contrib_pulay(cache%blks_triang_size* ( size_pulay_n + size_lagrange/3)))
       
    k = 1
    
!     write(*,*) 'h'    
    
    do while (traverse_end .EQV. .FALSE.)
  
       write(*,*) 'Outer contribution:'!, outer_next%num_dmat, outer_next%dummy_entry
    
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'B', outer_next%p_tuples(i)%pid
       
       end do
    
       
       
       
!        if ((outer_next%contrib_type == 1) .OR. (outer_next%contrib_type == 4)) then
      
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
          
          
          
          
          
!              write(*,*) 'calling w', j
!              write(*,*) 'wiip', which_index_is_pid
             
                  
!              write(*,*) 'back from w'
                  
             o_ctr = o_ctr + 1
       
              i = QcMatWrite_f(W(1), 'W1', ASCII_VIEW)
       
          end do
       
       
          deallocate(d_struct_o)
          deallocate(d_struct_o_prime)
       
!        end if
    
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do
    
!     if (o_triang_size == 78) then
!     
!     
!     i = QcMatWrite_f(W(1), 'W1', ASCII_VIEW)
!     i = QcMatWrite_f(W(12), 'W12', ASCII_VIEW)
!     i = QcMatWrite_f(W(43), 'W43', ASCII_VIEW)
!     i = QcMatWrite_f(W(52), 'W52', ASCII_VIEW)
! !     i = QcMatWrite_f(W(164), 'W164', ASCII_VIEW)
! !     i = QcMatWrite_f(W(343), 'W343', ASCII_VIEW)
! !     i = QcMatWrite_f(W(364), 'W364', ASCII_VIEW)
!     end if
    
    ! Outside traversal: Calculate contributions
    
    contrib_pulay = 0.0
    
    call get_ovl_exp(0, nof, 0, nof, size(pert_ext), pert_ext, size(W), W, &
                     size(contrib_pulay), contrib_pulay)
    
    
    
    write(*,*) 'Pulay contribution', contrib_pulay
    
    ! Traversal: Store Pulay contributions, calculate/store idempotency/SCFE contributions
    
    traverse_end = .FALSE.
    
    do while (outer_next%dummy_entry .eqv. .FALSE.)
        outer_next => outer_next%next
    end do
    outer_next => outer_next%next

    
    call QcMatInit(Y)
    call QcMatInit(Z)
    
    do i = 1, size(w)
    
       call QcMatInit(W(i))
    
    end do
    
    k = 1
    o_ctr = 1
    
    do while (traverse_end .EQV. .FALSE.)

!        write(*,*) 'contrib type', outer_next%contrib_type
    
       c_ctr = 0
    
       if (outer_next%p_tuples(1)%npert ==0) then
          
          o_triang_size = 1
          
       else
      
          o_triang_size = size(outer_next%indices, 1)
             
       end if
       
        
        
       ! Store Pulay n terms (if any)
       if ((outer_next%contrib_type == 1)) then
          
          do j = 1, o_triang_size * cache%blks_triang_size
          
             outer_next%data_scal(j) = contrib_pulay(o_ctr)
             o_ctr = o_ctr + 1
             c_ctr = c_ctr + 1
          
          end do
          
       end if
       
       c_snap = c_ctr
       
       ! Store Pulay Lagrange terms (if any)
       if ((outer_next%contrib_type >= 3)) then
          
          do j = 1, o_triang_size * cache%blks_triang_size
          
             outer_next%data_scal(j + c_snap) = &
             contrib_pulay(o_ctr)
             o_ctr = o_ctr + 1
             c_ctr = c_ctr + 1             
          
          end do
          
       end if
       
       c_snap = c_ctr
       
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
      
          
          
         
          
          do j = 1, o_triang_size
          
             call rsp_get_matrix_y(o_supsize_prime(k), d_struct_o, cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert, &
                  which_index_is_pid(1:cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                  F, D, S, Y)
             
             ! NOTE: Rule choice very likely to give correct exclusion but
             ! have another look if something goes wrong
             call rsp_get_matrix_z(o_supsize_prime(k), d_struct_o, &
                  (/outer_next%n_rule, outer_next%n_rule/), cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert, &
                  which_index_is_pid(1:cache%p_inner%npert + &
                  outer_next%p_tuples(1)%npert), size(outer_next%indices(j,:)), outer_next%indices(j,:), &
                  F, D, S, Z)
                  
             do i = 1, cache%blks_triang_size
           
!                 write(*,*) 'mat trace i j,', i, j
!                 write(*,*) 'c snap, j disp', c_snap, (j-1) * cache%blks_triang_size
!                 write(*,*) 'property displ', cache%blks_triang_size*o_triang_size
             
                call QcMatTraceAB(Lambda(i), Y, outer_next%data_scal(c_snap + &
                (j-1) * cache%blks_triang_size + i))
                call QcMatTraceAB(Zeta(i), Z, outer_next%data_scal(c_snap + &
                cache%blks_triang_size*o_triang_size + (j-1) * cache%blks_triang_size + i))
             
                write(*,*) 'Idempotency:', outer_next%data_scal(c_snap + &
                (j-1) * cache%blks_triang_size + i)
                
                write(*,*) 'SCFE:', outer_next%data_scal(c_snap + &
                cache%blks_triang_size*o_triang_size + (j-1) * cache%blks_triang_size + i)
             
             
             end do
                  
                  
             
       
       
          end do
       
       
          deallocate(d_struct_o)
      
      
       end if
       
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
       
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    

    
    
  end subroutine
  

  recursive subroutine rsp_pulay_n_2014(pert, kn, p12, S, D, F, get_ovl_exp, cache, p_size, prop)

    implicit none

    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    integer ::  i
    integer, dimension(2) :: kn
    external :: get_ovl_exp
    complex(8), dimension(*) :: prop
    integer :: num_blks_full, p_size
    integer, allocatable, dimension(:,:) :: blk_info_full

!     allocate(blk_info_full(num_blks_full, 3))
!         
!     num_blks_full = get_num_blks(pert)
!     
!     blk_info_full = get_blk_info(num_blks_full, pert)
!     p_size = get_triangulated_size(num_blks_full, blk_info_full)
!     
!     deallocate(blk_info_full)
    
    if (pert%npert > 0) then

       call rsp_pulay_n_2014(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, &
       get_ovl_exp, cache, p_size, prop)

       call rsp_pulay_n_2014(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, &
       get_ovl_exp, cache, p_size, prop)

    else

       if (kn_skip(p12(2)%npert, p12(2)%pid, kn) .EQV. .FALSE.) then



!          open(unit=257, file='totterms', status='old', action='write', &
!               position='append')
!          write(257,*) 'T'
!          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

!             open(unit=257, file='cachehit', status='old', action='write', &
!                  position='append') 
!             write(257,*) 'T'
!             close(257)

             call property_cache_getdata(cache, 2, p12, p_size, prop)
       
          else

          write(*,*) 'Calculating Pulay k-n contribution:'
          write(*,*) 'S', p12(1)%pid
          write(*,*) 'W', p12(2)%pid


             call get_pulay_n_2014((/ (p_tuple_standardorder(p12(i)) , i = 1, 2)  /), & 
                               kn, F, D, S, get_ovl_exp, cache, prop)

             write(*,*) 'Calculated Pulay k-n contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Pulay k-n contribution was k-n skipped:'
!           write(*,*) 'S ', p12(1)%pid 
!           write(*,*) 'W ', p12(2)%pid 
!           write(*,*) ' '

       end if 

    end if

  end subroutine


  subroutine get_pulay_n_2014(p12, kn, F, D, S, get_ovl_exp, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    type(QcMat) :: W
    external :: get_ovl_exp
    integer :: i, j, k, sstr_incr, offset, total_num_perturbations, &
                dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks, npert_ext
    integer, dimension(p12(1)%npert + p12(2)%npert) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged, pert_ext
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer :: d_supsize
    integer, dimension(0) :: noc
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, prop_forcache
    complex(8), dimension(*) :: prop

    call p_tuple_to_external_tuple(p12(1), npert_ext, pert_ext)

    
    d_supsize = derivative_superstructure_getsize(p12(2), kn, .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structb(d_supsize, 3))

    sstr_incr = 0

    call derivative_superstructure(p12(2), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, sstr_incr, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))
    
    
    
    do i = 1, 2
    
       nfields(i) = p12(i)%npert
       nblks_tuple(i) = get_num_blks(p12(i))
    
    end do
    
    total_num_perturbations = sum(nfields)
    
    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))
    
    do i = 1, 2
    
       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
    
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
    
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do
    
    if (p12(2)%npert == 0) then
    
       outer_indices_size = 1
    
    else
    
       outer_indices_size = blks_tuple_triang_size(2)
    
    end if
    
    inner_indices_size = blks_tuple_triang_size(1)
    
    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%npert + p12(2)%npert))
    allocate(ncinner(p12(1)%npert))
    allocate(tmp(inner_indices_size))
    allocate(inner_offsets(inner_indices_size))
    allocate(outer_indices(outer_indices_size, p12(2)%npert))
    allocate(inner_indices(inner_indices_size, p12(1)%npert))
    allocate(which_index_is_pid(p12(1)%npert + p12(2)%npert))

    prop_forcache = 0.0
    
    ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
    ncinner = nc_onlysmall(p12(1)%npert + p12(2)%npert, &
                           p12(1)%npert, 1, p12(1), ncarray)
    
    which_index_is_pid = 0
    
    do i = 1, p12(2)%npert
    
       which_index_is_pid(p12(2)%pid(i)) = i
    
    end do
    
    
    if (p12(2)%npert > 0) then
    
       call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
            1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices)
    
    end if
    
    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)
    
    call QcMatInit(W)
    
    do i = 1, size(outer_indices, 1)
    
       tmp = 0.0

       call QcMatZero(W)

       call rsp_get_matrix_w_2014(d_supsize, deriv_structb, p12(1)%npert + &
                             p12(2)%npert, which_index_is_pid, &
                             p12(2)%npert, outer_indices(i,:), F, D, S, W)
 
 
       call get_ovl_exp(0, noc, 0, noc, npert_ext, pert_ext, 1, (/W/), size(tmp), tmp)
       
 
!        call rsp_ovlave(p12(1)%npert, p12(1)%plab, &
!                           (/ (j/j, j = 1, p12(1)%npert) /), &
!                           p12(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
!                           1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), &
!                           size(tmp), W, tmp)
 
       do j = 1, size(inner_indices, 1)
    
          if (p12(2)%npert > 0) then
    
             offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                      (/ p12(1)%npert, p12(2)%npert /), &
                      blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/inner_indices(j, :), outer_indices(i, :) /)) 
          else
    
             offset = get_triang_blks_tuple_offset(1, total_num_perturbations, &
                      nblks_tuple(1), (/ p12(1)%npert /), &
                      blks_tuple_info(1,:,:), blk_sizes(1,:), blks_tuple_triang_size(1), &
                      (/inner_indices(j, :) /)) 

          end if
   
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do
    
    end do
    
    call p1_cloneto_p2(p12(1), merged_p_tuple)
    
    if (p12(2)%npert > 0) then
    
       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))
    
    end if
    
    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
    
    ! MaR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%npert
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do
    
    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
    
    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))
    
    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)
    
    do i = 1, size(triang_indices_pr, 1)
    
       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))
   
       do j = 1, total_num_perturbations
    
          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))
    
       end do
    
       if (p12(2)%npert > 0) then
    
          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))
    
       else
    
          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))
    
       end if
    
       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)
    
    end do
    

    write(*,*) 'after pulay kn contribution', prop(1:78)
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    
   
    call QcMatDst(W)

    deallocate(pert_ext)
    
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)

  end subroutine

  
  recursive subroutine rsp_pulay_lag_2014(pert, kn, p12, S, D, F, &
                       get_ovl_exp, cache, p_size, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    integer ::  i
    integer, dimension(2) :: kn
    external :: get_ovl_exp
    complex(8), dimension(*) :: prop
    integer :: num_blks_full, p_size
    integer, allocatable, dimension(:,:) :: blk_info_full

!     allocate(blk_info_full(num_blks_full, 3))
!         
!     num_blks_full = get_num_blks(pert)
!     
!     blk_info_full = get_blk_info(num_blks_full, pert)
!     p_size = get_triangulated_size(num_blks_full, blk_info_full)
!     
!     deallocate(blk_info_full)
    
    if (pert%npert > 0) then

       call rsp_pulay_lag_2014(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
       S, D, F, get_ovl_exp, cache, p_size, prop)
       call rsp_pulay_lag_2014(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
       S, D, F, get_ovl_exp, cache, p_size, prop)

    else

       ! At lowest level:
       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then



!       open(unit=257, file='totterms', status='old', action='write', position='append') 
!       write(257,*) 'T'
!       close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '
       
!             open(unit=257, file='cachehit', status='old', action='write', &
!             position='append') 
!             write(257,*) 'T'
!             close(257)

             call property_cache_getdata(cache, 2, p12, p_size, prop)

          else

       write(*,*) 'Calculating Pulay lagrange contribution:'
       write(*,*) 'S', p12(1)%pid
       write(*,*) 'W', p12(2)%pid, 'primed', kn(2)

             call get_pulay_lag_2014((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, get_ovl_exp, cache, prop)

             write(*,*) 'Calculated Pulay lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Pulay lagrange contribution was k-n skipped:'
!           write(*,*) 'S', p12(1)%pid 
!           write(*,*) 'W', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_pulay_lag_2014(p12, kn, F, D, S, get_ovl_exp, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    type(QcMat) :: W
    external :: get_ovl_exp
    integer :: i, j, k ,m, incr, npert_ext
    integer :: d_supsize, total_num_perturbations, &
              offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%npert + p12(2)%npert) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged, pert_ext, pert_ord_ext
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn
    integer, dimension(0) :: noc
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:) :: outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, prop_forcache
    complex(8), dimension(*) :: prop

    call p_tuple_to_external_tuple(p12(1), npert_ext, pert_ext)

    
    d_supsize = derivative_superstructure_getsize(p12(2), kn, .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
    allocate(deriv_structb(d_supsize, 3))

    incr = 0

    call derivative_superstructure(p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, incr, deriv_structb)


    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%npert
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%npert == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%npert + p12(2)%npert))
    allocate(ncinner(p12(1)%npert + p12(2)%npert))
    allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))
    allocate(tmp(inner_indices_size))
    allocate(inner_offsets(inner_indices_size))
    allocate(outer_indices(outer_indices_size, p12(2)%npert))
    allocate(inner_indices(inner_indices_size, p12(1)%npert))
    allocate(which_index_is_pid(p12(1)%npert + p12(2)%npert))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
    ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
              p12(1)%npert, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%npert

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices)


    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)

    call QcMatInit(W)

    do i = 1, size(outer_indices, 1)

       call QcMatZero(W)

       call rsp_get_matrix_w_2014(d_supsize, deriv_structb, p12(1)%npert + &
                            p12(2)%npert, which_index_is_pid, &
                            p12(2)%npert, outer_indices(i,:), F, D, S, W)

       tmp = 0.0
       call get_ovl_exp(0, noc, noc, 0, noc, noc, npert_ext, pert_ext, pert_ord_ext, 1, (/W/), size(tmp), tmp)
                            
!        call rsp_ovlave(p12(1)%npert, p12(1)%plab, &
!                        (/ (j/j, j = 1, p12(1)%npert) /), &
!                        p12(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
!                        1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), &
!                        size(tmp), W, tmp)

       do j = 1, size(inner_indices, 1)

          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%npert, p12(2)%npert /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/inner_indices(j, :), outer_indices(i, :) /)) 

          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do

    call p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%npert > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

    ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%npert
          pids_current_contribution(k) = p12(i)%pid(j)
       k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%npert > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do


!     write(*,*) 'pulay lag contribution'
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    


    deallocate(pert_ext)
    deallocate(pert_ord_ext)         
         
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_b_large)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(triang_indices_pr)
    deallocate(merged_blk_info)
    deallocate(prop_forcache)
    call QcMatDst(W)

  end subroutine
  
  
  ! NEW CODE
  
  
  recursive subroutine rsp_idem_lag_2014(pert, kn, p12, S, D, F, &
                                     cache, p_size, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    integer ::  i
    integer, dimension(2) :: kn
    complex(8), dimension(*) :: prop
    integer :: num_blks_full, p_size
    integer, allocatable, dimension(:,:) :: blk_info_full

!     allocate(blk_info_full(num_blks_full, 3))
!         
!     num_blks_full = get_num_blks(pert)
!     
!     blk_info_full = get_blk_info(num_blks_full, pert)
!     p_size = get_triangulated_size(num_blks_full, blk_info_full)
!     
!     deallocate(blk_info_full)
    
    if (pert%npert > 0) then

       call rsp_idem_lag_2014(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, &
       cache, p_size, prop)
       call rsp_idem_lag_2014(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, &
       cache, p_size, prop)

    else

       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then



!          open(unit=257, file='totterms', status='old', action='write', &
!               position='append') 
!          write(257,*) 'T'
!          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

!             open(unit=257, file='cachehit', status='old', action='write', &
!                  position='append')
!             write(257,*) 'T'
!             close(257)

             call property_cache_getdata(cache, 2, p12, p_size, prop)
      
          else

          write(*,*) 'Calculating idempotency lagrange contribution'
          write(*,*) 'Zeta', p12(1)%pid
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)

             ! At lowest level:
             call get_idem_lag_2014((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, cache, prop)

             write(*,*) 'Calculated idempotency lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Idempotency lagrange contribution was k-n skipped:'
!           write(*,*) 'Zeta', p12(1)%pid 
!           write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_idem_lag_2014(p12, kn, F, D, S, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    type(QcMat) :: Zeta, Z
    integer :: i, j, k, m, n, p, incr1, incr2, total_num_perturbations, &
              offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%npert + p12(2)%npert) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn
    integer, dimension(2) :: d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, which_index_is_pid1, &
                                          which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8) :: tmp_tr
    complex(8), dimension(*) :: prop
    complex(8), allocatable, dimension(:) :: prop_forcache

    d_supsize = 0
    d_supsize(1) = derivative_superstructure_getsize(p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(p_tuple_remove_first(p12(1)), kn, .FALSE., & 
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(2), incr2, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%npert
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%npert == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%npert + p12(2)%npert))
    allocate(ncinner(p12(1)%npert + p12(2)%npert))
    allocate(outer_ind_a_large(p12(1)%npert + p12(2)%npert))
    allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))
    allocate(outer_indices_a(inner_indices_size, p12(1)%npert))
    allocate(outer_indices_b(outer_indices_size, p12(2)%npert))
    allocate(which_index_is_pid1(p12(1)%npert + p12(2)%npert))
    allocate(which_index_is_pid2(p12(1)%npert + p12(2)%npert))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
    ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
                      p12(1)%npert, 1, p12(1), ncarray)

    which_index_is_pid1 = 0

    do i = 1, p12(1)%npert

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%npert

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), outer_indices_a)


    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices_b)

    offset = 0.0

    call QcMatInit(Z)
    call QcMatInit(Zeta)
    
    do i = 1, size(outer_indices_a, 1)

       call QcMatZero(Zeta)

       call rsp_get_matrix_zeta_2014(p_tuple_getone(p12(1), 1), kn, d_supsize(1), &
            deriv_structa, p12(1)%npert + p12(2)%npert, &
            which_index_is_pid1, p12(1)%npert, outer_indices_a(i,:), &
            F, D, S, Zeta)

       do j = 1, size(outer_indices_b, 1)

          call QcMatZero(Z)

          call rsp_get_matrix_z_2014(d_supsize(2), deriv_structb, kn, &
               p12(1)%npert + p12(2)%npert, which_index_is_pid2, &
               p12(2)%npert, outer_indices_b(j,:), F, D, S, Z)

          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%npert, p12(2)%npert /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/outer_indices_a(i, :), outer_indices_b(j, :) /)) 
                   
          call QcMatTraceAB(Zeta, Z, tmp_tr)
          prop_forcache(offset) = prop_forcache(offset) - tmp_tr

       end do

    end do

    call p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%npert > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%npert
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%npert > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do

!     write(*,*) 'idempotency contribution'
!  call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!  (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
!  merged_blk_info, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)
    
    call QcMatDst(Zeta)
    call QcMatDst(Z)

  end subroutine



  recursive subroutine rsp_scfe_lag_2014(pert, kn, p12, S, D, F, &
                                     cache, p_size, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    integer ::  i
    integer, dimension(2) :: kn
    complex(8), dimension(*) :: prop
    integer :: num_blks_full, p_size
    integer, allocatable, dimension(:,:) :: blk_info_full

!     allocate(blk_info_full(num_blks_full, 3))
!         
!     num_blks_full = get_num_blks(pert)
!     
!     blk_info_full = get_blk_info(num_blks_full, pert)
!     p_size = get_triangulated_size(num_blks_full, blk_info_full)
!     
!     deallocate(blk_info_full)
    
    if (pert%npert > 0) then

       call rsp_scfe_lag_2014(p_tuple_remove_first(pert), kn, &
            (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
            S, D, F, cache, p_size, prop)
       call rsp_scfe_lag_2014(p_tuple_remove_first(pert), kn, &
            (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
            S, D, F, cache, p_size, prop)

    else

       if (kn_skip(p12(1)%npert, p12(1)%pid, kn) .EQV. .FALSE.) then



!          open(unit=257, file='totterms', status='old', action='write', &
!               position='append') 
!          write(257,*) 'T'
!          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!             open(unit=257, file='cachehit', status='old', action='write', &
!                  position='append') 
!             write(257,*) 'T'
!             close(257)

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

             call property_cache_getdata(cache, 2, p12, p_size, prop)
       
          else

          write(*,*) 'Calculating scfe lagrange contribution'
          write(*,*) 'Lambda', p12(1)%pid
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)

             ! At lowest level:
             call get_scfe_lag_2014((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), &
             kn, F, D, S, cache, prop)

             write(*,*) 'Calculated scfe lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'scfe lagrange contribution was k-n skipped:'
!           write(*,*) 'Lambda', p12(1)%pid 
!           write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_scfe_lag_2014(p12, kn, F, D, S, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF_2014) :: S, D, F
    type(property_cache) :: cache
    type(QcMat) :: L, Y
    integer :: i, j, k, m, n, p, incr1, incr2, total_num_perturbations, &
              offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%npert + p12(2)%npert) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged, nfields, nblks_tuple
    integer, allocatable, dimension(:) :: blks_tuple_triang_size
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn
    integer, dimension(2) :: d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, which_index_is_pid1, which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8) :: tmp_tr
    complex(8), dimension(*) :: prop
    complex(8), allocatable, dimension(:) :: prop_forcache

    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(p_tuple_remove_first(p12(1)), kn, .FALSE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), & 
                    d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(p12(2), kn, .TRUE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                    d_supsize(2), incr2, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%npert
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%npert == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%npert + p12(2)%npert))
    allocate(ncinner(p12(1)%npert + p12(2)%npert))
    allocate(outer_indices_a(inner_indices_size, p12(1)%npert))
    allocate(outer_indices_b(outer_indices_size, p12(2)%npert))
    allocate(outer_ind_a_large(p12(1)%npert + p12(2)%npert))
    allocate(outer_ind_b_large(p12(1)%npert + p12(2)%npert))
    allocate(which_index_is_pid1(p12(1)%npert + p12(2)%npert))
    allocate(which_index_is_pid2(p12(1)%npert + p12(2)%npert))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%npert + p12(2)%npert, 2, p12)
    ncinner = nc_only(p12(1)%npert + p12(2)%npert, &
              p12(1)%npert, 1, p12(1), ncarray)

    which_index_is_pid1 = 0

    do i = 1, p12(1)%npert

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%npert

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), outer_indices_a)

    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices_b)

    offset = 0

    ! ASSUME CLOSED SHELL
    call QcMatInit(Y)
    
    ! ASSUME CLOSED SHELL
    call QcMatInit(L)

    do i = 1, size(outer_indices_a, 1)

       call QcMatZero(L)

       call rsp_get_matrix_lambda_2014(p_tuple_getone(p12(1), 1), d_supsize(1), &
            deriv_structa, p12(1)%npert + p12(2)%npert, &
            which_index_is_pid1, p12(1)%npert, outer_indices_a(i,:), D, S, L)

       do j = 1, size(outer_indices_b, 1)

          call QcMatZero(Y)

          call rsp_get_matrix_y_2014(d_supsize(2), deriv_structb, &
               p12(1)%npert + p12(2)%npert, which_index_is_pid2, &
               p12(2)%npert, outer_indices_b(j,:), F, D, S, Y)


          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%npert, p12(2)%npert /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/outer_indices_a(i, :), outer_indices_b(j, :) /)) 
          
          call QcMatTraceAB(L, Y, tmp_tr)
          prop_forcache(offset) = prop_forcache(offset) - tmp_tr

       end do

    end do

    call p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%npert > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

    ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%npert
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%npert > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do

!     write(*,*) 'scfe contribution'
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)

    call QcMatDst(L)
    call QcMatDst(Y)

  end subroutine
  
  
  ! END NEW CODE
  
  
  
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
