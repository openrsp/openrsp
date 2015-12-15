! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

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
                                    rsp_get_matrix_y_2014
  use rsp_perturbed_sdf, only: rsp_fds_2014, &
                               rsp_fds
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
                                  contrib_cache_add_element,           &
                                  contrib_cache_already,               &
                                  contrib_cache_getdata,               &
                                  contrib_cache_getdata_outer,               &
                                  contrib_cache_allocate
  
  use rsp_sdf_caching
  
  use qcmatrix_f

  implicit none

!   public get_prop_2014
  public rsp_energy_2014
  public get_energy_2014
  public rsp_pulay_n_2014
  public get_pulay_n_2014
  public rsp_pulay_lag_2014
  public get_pulay_lag_2014
  public rsp_idem_lag_2014
  public get_idem_lag_2014
  public rsp_scfe_lag_2014
  public get_scfe_lag_2014
  public print_rsp_tensor
  public print_rsp_tensor_stdout
  public print_rsp_tensor_stdout_tr
  public openrsp_get_property_2014
  
  private

  real(8) :: time_start
  real(8) :: time_end

  contains
  
!   subroutine openrsp_get_property_2014(nprops, np, pert_dims, pert_first_comp, pert_labels, num_freq_cfgs, pert_freqs, &
!                                    kn_rules, F_unpert, S_unpert, D_unpert, get_rsp_sol, get_nucpot, &
!                                    get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
!                                    get_2el_mat, get_2el_exp, get_xc_mat, & 
!                                    get_xc_exp, id_outp, rsp_tensor, file_id)
!     implicit none
! 
!     integer(kind=QINT), intent(in) :: nprops
!     integer(kind=QINT), dimension(nprops), intent(in) :: np, num_freq_cfgs
!     integer(kind=4), intent(in) :: id_outp
!     integer(kind=QINT), dimension(sum(np)), intent(in) :: pert_dims, pert_first_comp
!     character(4), dimension(sum(np)), intent(in) :: pert_labels
!     integer :: i, j, num_blks, prop_size
!     integer(kind=QINT), intent(in), dimension(nprops) :: kn_rules
!     integer, dimension(2) :: kn_rule
!     character, optional, dimension(20) :: file_id
!     integer, allocatable, dimension(:) :: blk_sizes
!     integer, allocatable, dimension(:,:) :: blk_info
!     complex(8), dimension(dot_product(np, num_freq_cfgs)), intent(in) :: pert_freqs
!     integer(kind=QINT) num_perts
!     real :: timing_start, timing_end
!     type(p_tuple) :: perturbations
!     external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
!     external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
!     complex(8), dimension(*) :: rsp_tensor
!     type(QcMat) :: S_unpert, D_unpert, F_unpert
!     type(SDF_2014), pointer :: S, D, F
!     integer kn(2)
! 
!     if (nprops/=1) then
!        write(*,*) 'ERROR: Only one property at a time supported for now'
! 
!     else
! 
!     if (num_freq_cfgs(1)/=1) then
! 
!        write(*,*) 'ERROR: Only one frequency configuration for each property supported for now'
! 
!     else
! 
!     call sdf_setup_datatype_2014(S, S_unpert)
!     call sdf_setup_datatype_2014(D, D_unpert)
!     call sdf_setup_datatype_2014(F, F_unpert)
!     
!     
!     ! Loop over properties should start here
!     
!     num_perts = np(1)
! 
!     kn_rule(1) = kn_rules(1)
!     kn_rule(2) = num_perts - 1 - kn_rules(1)
! 
!     
!     kn(1) = kn_rule(1)
!     kn(2) = kn_rule(2)
!     
!     perturbations%npert = num_perts
!     allocate(perturbations%pdim(num_perts))
! !     allocate(perturbations%pfcomp(num_perts))
!     allocate(perturbations%plab(num_perts))
!     allocate(perturbations%pid(num_perts))
!     allocate(perturbations%freq(num_perts))
!     perturbations%pdim = pert_dims
! 
!     do i = 1, num_perts
!          perturbations%plab(i) = pert_labels(i)
!     end do
!     
!     perturbations%pid = (/(i, i = 1, num_perts)/)
!     
!     write(*,*), 'num perts first', num_perts
!     write(*,*), 'in p tuple first', perturbations%npert
! 
!     
!     ! Loop over frequencies should start here
!     
!     
!     perturbations%freq = pert_freqs(1:num_perts)
!     
!     num_blks = get_num_blks(perturbations)
!     allocate(blk_info(num_blks, 3))
!     allocate(blk_sizes(num_blks))
!     blk_info = get_blk_info(num_blks, perturbations)
!     blk_sizes = get_triangular_sizes(num_blks, blk_info(1:num_blks, 2), &
!                                      blk_info(1:num_blks, 3))
! 
!     prop_size = get_triangulated_size(num_blks, blk_info)
!                                      
!                                      
!     write(id_outp,*) ' '
!     write(id_outp,*) 'OpenRSP lib called'
!     write(id_outp,*) ' '
!     write(id_outp,*) 'Calculating a property of order ', num_perts
!     write(id_outp,*) 'The choice of k, n is ', kn(1), ' and ', kn(2)
!     write(id_outp,*) ' '
!     write(id_outp,*) 'The number of components for each perturbation is:    ', pert_dims
! !     write(id_outp,*) 'The first component of each perturbation is:          ', pert_first_comp
!     write(id_outp,*) 'The perturbation labels are:                          ', pert_labels
!     write(id_outp,*) 'The frequencies of the perturbations (real part) are: ', (/(real(pert_freqs(i)), i = 1, num_perts)/)
!     write(id_outp,*) 'The frequencies of the perturbations (imag. part) are:', (/(aimag(pert_freqs(i)), i = 1, num_perts)/)
!     write(id_outp,*) ' '
!  
!     if ((kn(1) - kn(2) > 1) .OR. .NOT.(kn(1) + kn(2) == num_perts - 1)) then
! 
!        write(id_outp,*) 'ERROR: Invalid choice of (k,n)'
!        write(id_outp,*) 'Valid choices for k are integers between and including 0 and ', (num_perts - 1)/2
!        write(id_outp,*) 'Valid choices of n are such that k + n =', num_perts - 1
!        write(id_outp,*) 'Cannot proceed with calculation: Exiting OpenRSP lib'
!        write(id_outp,*) ' '
!        return
!  
!     end if
! 
!     
!     write(id_outp,*) 'Starting clock: About to call get_prop routine'
!     write(id_outp,*) ' '
! 
!     call cpu_time(timing_start)
! 
!     call get_prop_2014(perturbations, kn, num_blks, blk_sizes, blk_info, F, D, S, get_rsp_sol, &
!                   get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp, &
!                   get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp, &
!                   id_outp, prop_size, rsp_tensor)
! 
!     call cpu_time(timing_end)
! 
!     write(id_outp,*) 'Clock stopped: Property was calculated'
!     write(id_outp,*) 'Time spent in get_prop:',  timing_end - timing_start, ' seconds'
!     write(id_outp,*) ' '
! 
!     write(id_outp,*) 'Property was calculated'
!     write(id_outp,*) ' '
!     
! 
! 
!     if (present(file_id)) then
! !       open(unit=260, file='rsp_tensor_' // trim(adjustl(file_id)), &
! !             status='replace', action='write') 
! !        open(unit=261, file='rsp_tensor_human_' // trim(adjustl(file_id)), &
! !             status='replace', action='write') 
!     else
!        open(unit=260, file='rsp_tensor', &
!             status='replace', action='write') 
!        open(unit=261, file='rsp_tensor_human', &
!             status='replace', action='write') 
!             
!        call print_rsp_tensor_tr(1, perturbations%npert, perturbations%pdim, &
!        (/ (1, j = 1, (perturbations%npert - 1) ) /), num_blks, blk_sizes, &
!        blk_info, prop_size, rsp_tensor, 260, 261)
! 
!     close(260)
!     close(261)
!     end if
! 
! 
! 
!     if (present(file_id)) then
! !        write(*,*) 'Property was printed to rsp_tensor_' // trim(adjustl(file_id))
! !        write(*,*) 'Property (formatted print) was printed to rsp_tensor_human_' &
! !                    // trim(adjustl(file_id)) 
!     else
!        write(*,*) 'Property was printed to rsp_tensor'
!        write(*,*) 'Property (formatted print) was printed to rsp_tensor_human'
!     end if
! 
!     write(*,*) ' '
!     write(*,*) 'End of print'
! 
! !     call print_rsp_tensor(1, perturbations%npert, perturbations%pdim, &
! !     (/ (1, j = 1, (perturbations%npert - 1) ) /), num_blks, blk_sizes, &
! !     blk_info, prop, id_outp, id_outp)
! 
!     deallocate(blk_info)
! 
!   end if
! 
!   end if
! 
!   end subroutine
  
  
  
  
  
  
  
  
  
  subroutine openrsp_get_property_2014(n_props, np, pert_dims, pert_first_comp, pert_labels, n_freq_cfgs, pert_freqs, &
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
    integer kn(2)

    call contrib_cache_outer_allocate(S)
    call contrib_cache_outer_allocate(D)
    call contrib_cache_outer_allocate(F)
    
    call contrib_cache_outer_add_element(S, .FALSE., 1, (/get_emptypert()/), &
         data_size = 1, data_mat=(/S_unpert/))
    call contrib_cache_outer_add_element(D, .FALSE., 1, (/get_emptypert()/), &
         data_size = 1, data_mat=(/D_unpert/))
    call contrib_cache_outer_add_element(F, .FALSE., 1, (/get_emptypert()/), &
         data_size = 1, data_mat=(/F_unpert/))         

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
          ! Unsure of addressing in next line
          

          p_tuples(k)%freq = pert_freqs(dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + &
          1 + (j - 1)*np(i):dot_product(np(1:i), n_freq_cfgs(1:i)) - np(i)*n_freq_cfgs(i) + 1 + (j)*np(i))
       
          write(id_outp,*) 'Frequency configuration', j
          write(id_outp,*) ' '
          write(id_outp,*) 'Frequencies (real part):', (/(real(p_tuples(k)%freq(m)), m = 1, np(i))/)
          write(id_outp,*) 'Frequencies (imag. part):', (/(aimag(p_tuples(k)%freq(m)), m = 1, np(i))/)
          write(id_outp,*) ' '
       
              
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

    integer :: n_props, id_outp, i
    integer, dimension(n_props) :: n_freq_cfgs
    integer, dimension(sum(n_freq_cfgs)) :: prop_sizes
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    type(p_tuple) :: emptypert
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), dimension(2) :: emptyp_tuples
    external :: get_rsp_sol, get_nucpot, get_ovl_mat, get_ovl_exp, get_1el_mat, get_1el_exp
    external :: get_2el_mat, get_2el_exp, get_xc_mat, get_xc_exp
    complex(8), dimension(*) :: props
    type(contrib_cache), pointer :: contribution_cache
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
    
!     write(*,*) 'num perts', pert%npert

    do i = 1, n_props


       call property_cache_allocate(contribs_cache)
       call contrib_cache_allocate(contribution_cache)
       write(id_outp,*) ' '
       write(id_outp,*) 'Calculating HF-energy type contribs'
       write(id_outp,*) ' '

       call cpu_time(time_start)
       call rsp_energy_recurse(p_tuples(i), p_tuples(i)%npert, kn_rule(i,:), 1, (/emptypert/), &
            0, D, get_nucpot, get_1el_exp, get_ovl_exp, get_2el_exp, .TRUE., &
            contribution_cache, prop_sizes(i), &
            props(sum(prop_sizes(1:i)) - prop_sizes(i) + 1:sum(prop_sizes(1:i))))
       call cpu_time(time_end)

       write(id_outp,*) 'Time spent:', time_end - time_start, 'seconds'
       write(id_outp,*) 'Finished calculating HF energy-type contribs'
       write(id_outp,*) ' '
       deallocate(contribs_cache)
       
       write(*,*) 'Property is now', &
       props(sum(prop_sizes(1:i)) - prop_sizes(i) + 1:sum(prop_sizes(1:i)))

    end do
  
  
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

    ! After dryrun recursion is done, calculate all necessary cache elements 
    ! and recurse again to get cache retrieval situations
    if ((pert%npert == total_num_perturbations) .AND. dryrun) then
    
    
       traverse_end = .FALSE.
!        write(*,*) 'dryrun is done'

       cache_next => cache

       ! Cycle to last element of cache
       do while (cache_next%last .eqv. .FALSE.)
          cache_next => cache_next%next
       end do

       cache_next => cache_next%next
       cache_next => cache_next%next
       
       ! Traverse linked list while getting contribs until at last element again
       do while (traverse_end .eqv. .FALSE.)
       
       write(*,*) 'Calculating energy-type contribs for inner perturbation tuple'
       write(*,*) cache_next%p_inner%pid

          call rsp_energy_calculate(D, get_nucpot, get_1el_exp, get_t_exp, get_2el_exp, cache_next)
                 write(*,*) ' '
          
          if (cache_next%last) then
             traverse_end = .TRUE.
          end if
          
          cache_next => cache_next%next
          
       end do
       
!        write(*,*) 'finished calculation'
! 
!        ! Do new recursion for cache retrieval situations
!        call rsp_energy_recurse(pert, total_num_perturbations, kn, num_p_tuples, &
!                        p_tuples, density_order, D, get_nucpot, get_1el_exp,  &
!                        get_t_exp, get_2el_exp, .FALSE., cache_next, p_size=p_size, prop=prop)
!                        
! !        write(*,*) 'finished new recursion'

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
    
    ! Allocate and set up outer
    
    allocate(LHS_dmat_1(sum(outer_contract_sizes_1(:))))
    allocate(LHS_dmat_2(sum(outer_contract_sizes_2(:, 1))))
    allocate(RHS_dmat_2(sum(outer_contract_sizes_2(:, 2))))
    
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

          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .TRUE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=LHS_dmat_2(lhs_ctr_2))
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .TRUE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))

       ! One chain rule application
       else if (outer_next%num_dmat == 1) then
       
          do m = 1, outer_contract_sizes_1(k) 
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .TRUE., &
                  1, ind_len=size(outer_next%indices, 2), ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_1(lhs_ctr_1 + m  - 1))
          end do
       
          do m = 1, outer_contract_sizes_2(k, 1) 
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .TRUE., &
                  1, ind_len=size(outer_next%indices, 2), ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
          end do
          
          call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .TRUE., &
               1, ind_len=1, ind_unsorted=(/1/), mat_sing=RHS_dmat_2(rhs_ctr_2))
       
       ! Two chain rule applications
       else if (outer_next%num_dmat == 2) then
       
          do m = 1, outer_contract_sizes_2(k, 1) 
             call contrib_cache_getdata_outer(D, 1, outer_next%p_tuples, .TRUE., &
                  1, ind_len=outer_next%p_tuples(1)%npert, &
                  ind_unsorted=outer_next%indices(m, 1:outer_next%p_tuples(1)%npert), &
                  mat_sing=LHS_dmat_2(lhs_ctr_2 + m  - 1))
          end do
          
          do n = 1, outer_contract_sizes_2(k, 2) 
             call contrib_cache_getdata_outer(D, 1, outer_next%p_tuples, .TRUE., &
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
    
    allocate(contrib_0(cache%blks_triang_size))
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    allocate(contrib_2(cache%blks_triang_size*total_outer_size_2))
    
    ! Calculate contributions
    
    ! Calculate nuclear-nuclear repulsion contribution
    if (num_0 > 0) then
    
       contrib_0 = 0.0
       call get_nucpot(num_pert, pert_ext, size(contrib_0), contrib_0)
    
    end if
    
    ! Calculate one-electron contributions
    if (num_1 > 0) then
    
       contrib_1 = 0.0
       call get_1el_exp(num_pert, pert_ext, num_1, outer_contract_sizes_1_coll, &
                        LHS_dmat_1, size(contrib_1), contrib_1)
      
       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()
      
!        call rsp_ovlave_t_matrix_2014(get_ovl_exp, cache%p_inner, cache%p_inner%npert, &
!                                 t_matrix_bra, t_matrix_ket, outer_contract_sizes_1_coll, &
!                                 LHS_dmat_1, size(contrib_1), contrib_1)
    
    end if
    
    ! Calculate two-electron contributions
    contrib_2 = 0.0
    call get_2el_exp(num_pert, pert_ext, cache%num_outer, outer_contract_sizes_2(:, 1), LHS_dmat_2, & 
                     outer_contract_sizes_2(:, 2), RHS_dmat_2, size(contrib_2), contrib_2)
                       
    
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
!     total_contrib_size = 0
!     total_num_outer = total_num_perturbations - cache%p_inner%npert
!     
!     inner_nblks = get_num_blks(cache%p_inner)
!     inner_blks_tuple_info = get_blk_info(inner_nblks, cache%p_inner)
!     inner_triang_size = get_triangulated_size(inner_nblks, inner_blks_tuple_info(1:inner_nblks, :))
!     inner_blk_sizes = get_triangular_sizes(inner_nblks, inner_blks_tuple_info(1:inner_nblks, 2), &
!                       inner_blks_tuple_info(1:inner_nblks, 3))
!     
!     call make_triangulated_tuples_indices(1, cache%p_inner%npert, (/inner_nblks/), &
!          inner_blks_tuple_info, (/inner_triang_size/), inner_indices)
!     
!     if (total_num_outer > 0) then
!     
!        do i = 1, cache%num_outer
!     
!           outer_next => cache%contribs_outer%next
!        
!           allocate(outer_next%nblks_tuple(outer_next%num_dmat))
!           allocate(outer_next%blk_sizes(outer_next%num_dmat, total_num_outer))
!           allocate(outer_next%blks_tuple_info(outer_next%num_dmat, total_num_outer, 3))
!           allocate(outer_next%blks_tuple_triang_size((outer_next%num_dmat)))
!        
!           do j = 1, outer_next%num_dmat
!        
!              outer_next%nblks_tuple(j) = get_num_blks(outer_next%outer_p_tuples(j))
!              outer_next%blks_tuple_info(j, :, :) = get_blk_info(outer_next%nblks_tuple(j), & 
!                                                    outer_next%outer_p_tuples(j))
!              outer_next%blks_tuple_triang_size(j) = get_triangulated_size(outer_next%nblks_tuple(j), &
!                                                     outer_next%blks_tuple_info(j, 1:outer_next%nblks_tuple(j), :))
!              outer_next%blk_sizes(j, 1:outer_next%nblks_tuple(j)) = &
!              get_triangular_sizes(outer_next%nblks_tuple(j), &
!              outer_next%blks_tuple_info(i,1:outer_next%nblks_tuple(j),2), &
!              outer_next%blks_tuple_info(i,1:outer_next%nblks_tuple(j),3))
!              outer_next%outer_size = product(outer_next%blks_tuple_triang_size)
!        
!           end do
!        
!           allocate(outer_next%indices(outer_next%outer_size, total_num_outer))
!        
!           ! MaR: Unsure about 2nd argument (maybe total_num_perturbations)
!           call make_triangulated_tuples_indices(outer_next%num_dmat, total_num_outer, & 
!                outer_next%nblks_tuple, outer_next%blks_tuple_info, &
!                outer_next%blks_tuple_triang_size, outer_next%indices)
!        
!           total_outer_size = total_contrib_size + outer_next%outer_size
!        
!     
!        end do
!     
!        allocate(dens_tuples(total_outer_size,2))
!        allocate(contrib(inner_triang_size * total_outer_size))    
!     
!        call empty_p_tuple(emptypert)
!     
!        do i = 1, total_outer_size
!     
!           call sdf_getdata_s(D, emptypert, (/1/), dens_tuples(i, 2))
!     
!        end do
!     
!        m = 0
!     
!        do i = 1, cache%num_outer
!     
!           outer_next => cache%contribs_outer%next
!        
!           do j = 1, outer_next%outer_size
!        
!              m = m + 1
!        
!              istart = 1
!              iend = 0
!        
!              do k = 1, outer_next%num_dmat
!           
!                 iend = iend + outer_next%outer_p_tuples(k)%npert
!           
!                 call sdf_getdata_s(D, outer_next%outer_p_tuples(k), outer_next%indices(j,istart:iend), &
!                      dens_tuples(m, k))
!                   
!                 istart = iend + 1
!        
!              end do
!           
!           end do
!           
!        end do
!                 ! MaR: Restore this later   
! !        call get_1el_contrib(cache%p_inner%npert, cache%p_inner%pdim, cache%p_inner%pfcomp, &
! !                             cache%p_inner%plab, cache%p_inner%freq, 1, total_outer_size, &
! !                             dens_tuples(:,1), contrib)
! 
! !        call empty_p_tuples(t_matrix_bra)
! !        call empty_p_tuples(t_matrix_ket)
! 
!        t_mat_p_tuple = cache%p_inner
!        t_mat_p_tuple%pid = (/(i, i = 1, t_mat_p_tuple%npert)/)
! 
!        ! MaR: Restore this later
!            
! !        call rsp_ovlave_t_matrix(get_tmatrix_contrib, t_mat_p_tuple%npert, t_mat_p_tuple, &
! !                                 t_matrix_bra, t_matrix_ket, &
! !                                 dens_tuples(:,1), total_outer_size, contrib)                         
!                          
! !        call get_2el_contrib(cache%p_inner%npert, cache%p_inner%pdim, cache%p_inner%pfcomp, &
! !                             cache%p_inner%plab, cache%p_inner%freq, 1, total_outer_size, &
! !                             dens_tuples, contrib)
!  
!      ! Loop to get offset and put in rsp tensor and cache
!     
!        do i = 1, cache%num_outer
!     
!           outer_next => cache%contribs_outer%next
!        
!           do j = 1, outer_next%outer_size
! 
!              if (cache%p_inner%npert > 0) then
!           
!                 do k = 1, inner_triang_size
!              
!                    ! Merging of blks_tuple_info can give mem issues
!                    offset = get_triang_blks_tuple_offset(outer_next%num_dmat + 1, &
!                    total_num_perturbations, (/inner_nblks, outer_next%nblks_tuple(:)/), &
!                    (/cache%p_inner%npert, (outer_next%outer_p_tuples(m)%npert, m = 1, outer_next%num_dmat)/), &
!                    (/inner_blks_tuple_info, (outer_next%blks_tuple_info(m, :, :), m = 1, outer_next%num_dmat)/), &
!                    (/inner_blk_sizes, (outer_next%blk_sizes(m,:), m = 1, outer_next%num_dmat)/), &
!                    (/inner_triang_size, (outer_next%blks_tuple_triang_size(m), m = 1, outer_next%num_dmat)/), &
!                    (/inner_indices(k,:), outer_next%indices(j,:)/))
! 
!                    prop_forcache(offset) = prop_forcache(offset) + contrib(k)
! 
!                 end do
!              
!              else
!              
!                 offset = get_triang_blks_tuple_offset(outer_next%num_dmat, total_num_perturbations, &
!                 outer_next%nblks_tuple, (/(outer_next%outer_p_tuples(m)%npert, m = 1,outer_next% num_dmat)/), &
!                 outer_next%blks_tuple_info, outer_next%blk_sizes, outer_next%blks_tuple_triang_size, &
!                 outer_next%indices(j,:))
! 
!                 prop_forcache(offset) = prop_forcache(offset) + contrib(k)
! 
!              end if
!           
!           end do
!           
!                    
!           
!           if (cache%p_inner%npert > 0) then
!     
!              call p1_cloneto_p2(cache%p_inner, merged_p_tuple)
!     
!              do j = 1, outer_next%num_dmat
!     
!                 call p1_merge_p2(outer_next%outer_p_tuples(j), merged_p_tuple, merged_p_tuple)
!    
!              end do
!    
!           else
!    
!              call p1_cloneto_p2(outer_next%outer_p_tuples(1), merged_p_tuple)
!    
!              do j = 2, outer_next%num_dmat
!    
!                 call p1_merge_p2(outer_next%outer_p_tuples(j), merged_p_tuple, merged_p_tuple)
!    
!              end do
!    
!           end if
!    
!           call p_tuple_ordered(merged_p_tuple, merged_p_tuple)
!    
!           ! MaR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
!           ! PIDS ARE IN STANDARD ORDER? FIND OUT
!    
!           m = 1
!       
!           do j = 1, cache%p_inner%npert
!              pids_current_contrib(j) = cache%p_inner%pid(j)
!              m = m + 1
!           end do
!           
!           do j = 1, outer_next%num_dmat
!              do k = 1, outer_next%outer_p_tuples(j)%npert
!                 pids_current_contrib(k) = outer_next%outer_p_tuples(j)%pid(k)
!                 m = m + 1
!              end do
!           end do
!     
!           merged_nblks = get_num_blks(merged_p_tuple)
!        
!           allocate(merged_blk_info(1, merged_nblks, 3))
!        
!           merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
!           blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
!           merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
!           merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
!     
!           allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))
!     
!           call make_triangulated_indices(merged_nblks, merged_blk_info, & 
!                                          merged_triang_size, triang_indices_pr)
!     
!        
!           do j = 1, size(triang_indices_pr, 1)
!     
!              prop_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
!                          (/merged_p_tuple%npert/), &
!                          (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
!                          (/triang_indices_pr(j, :) /))
!     
!              do k = 1, total_num_perturbations
!     
!                 translated_index(k) = triang_indices_pr(j,pids_current_contrib(k))
!     
!              end do
! 
!              
!              if (cache%p_inner%npert > 0) then
! 
!                 ! Merging of blks_tuple_info can give mem issues
!                 cache_offset = get_triang_blks_tuple_offset(outer_next%num_dmat + 1, &
!                 total_num_perturbations, (/inner_nblks, outer_next%nblks_tuple(:)/), &
!                 (/cache%p_inner%npert, (outer_next%outer_p_tuples(m)%npert, m = 1, num_dmat)/), &
!                 (/inner_blks_tuple_info, (outer_next%blks_tuple_info(m, :, :), m = 1, num_dmat)/), &
!                 (/inner_blk_sizes, (outer_next%blk_sizes(m,:), m = 1, num_dmat)/), &
!                 (/inner_triang_size, (outer_next%blks_tuple_triang_size(m), m = 1, num_dmat)/), &
!                 (/ translated_index(:) /))
!     
!              else
! 
!                 cache_offset = get_triang_blks_tuple_offset(outer_next%num_dmat, total_num_perturbations, &
!                 outer_next%nblks_tuple, (/(outer_next%outer_p_tuples(m)%npert, m = 1, num_dmat)/), &
!                 outer_next%blks_tuple_info, outer_next%blk_sizes, outer_next%blks_tuple_triang_size, &
!                 (/ translated_index(:) /))
! 
!              end if
!    
!              prop(prop_offset) = prop(prop_offset) + prop_forcache(cache_offset)
!    
!           end do
!    
!           outer_next%contrib_size = inner_triang_size * outer_next%outer_size
!           allocate(outer_next%data_scal(outer_next%contrib_size))
!           outer_next%data_scal = prop_forcache
!            
!        end do
!  
!  
!     else
!                   ! MaR: Restore this later
! !        call get_nucpot_contrib(cache%p_inner%npert, cache%p_inner%pdim, cache%p_inner%pfcomp, &
! !                             cache%p_inner%plab, cache%p_inner%freq, prop)
! 
!        call empty_p_tuple(emptypert)
!        call sdf_getdata_s(D, emptypert, (/1/), D_unp)
!                ! MaR: Restore this later                           
! !        call get_1el_contrib(cache%p_inner%npert, cache%p_inner%pdim, cache%p_inner%pfcomp, &
! !                             cache%p_inner%plab, cache%p_inner%freq, 1, 1, &
! !                             (/D_unp/), prop)
!     
!        t_mat_p_tuple = cache%p_inner
!        t_mat_p_tuple%pid = (/(i, i = 1, t_mat_p_tuple%npert)/)
! 
!               ! MaR: Restore this later
!            
! !        call rsp_ovlave_t_matrix(get_tmatrix_contrib, t_mat_p_tuple%npert, t_mat_p_tuple, &
! !                                 t_matrix_bra, t_matrix_ket, &
! !                                 (/D_unp/), 1, prop)         
!     
!     
! !        call get_2el_contrib(cache%p_inner%npert, cache%p_inner%pdim, cache%p_inner%pfcomp, &
! !                              cache%p_inner%plab, cache%p_inner%freq, 1, 1, &
! !                              (/D_unp, D_unp/), prop)
!    
!     end if
! 
! !  write(*,*) 'energy contrib'
! !  call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
! !  (/ (1, j = 1, (merged_p_tuple%npert - 1) ) /), merged_nblks, blk_sizes_merged, &
! !  merged_blk_info, prop_forcache)
! 
! 
!   end subroutine

  ! END NEW 2014
  
  
  
  
  

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
