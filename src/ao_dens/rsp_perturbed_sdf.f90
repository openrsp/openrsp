! Copyright 2015 Magnus Ringholm
!
! Contains routines for the calculation of perturbed
! overlap, density and Fock matrices used throughout the
! rsp_general calculation.

module rsp_perturbed_sdf

  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_property_caching
  use qcmatrix_f
  
  implicit none

  public rsp_fds
  public rsp_xc_wrapper

  private

  contains
    
  ! Main routine for managing calculation of perturbed Fock, density and overlap matrices
  subroutine rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, len_sdf, F, D, S, get_rsp_sol, get_ovl_mat, &
                               get_1el_mat, get_2el_mat, get_xc_mat, out_print, dryrun, &
                               prog_info, rs_info, r_flag, sdf_retrieved, mem_mgr, Xf)

    implicit none

    type(mem_manager) :: mem_mgr
    integer :: n_props
    integer, dimension(n_props) :: n_freq_cfgs
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), allocatable, dimension(:) :: p_dummy_orders
    logical :: termination, dryrun, lof_retrieved, sdf_retrieved, rsp_eqn_retrieved, traverse_end
    logical :: residue_select
    integer, dimension(3) :: prog_info, rs_info
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    integer :: i, j, k, m, max_order, max_npert, o_size, lof_mem_total
    integer :: r_flag
    integer :: c_ord
    integer :: len_curr_outer, len_d, len_lof_cache
    integer :: len_cache, len_sdf
    integer, allocatable, dimension(:) :: size_i
    type(QcMat) :: Fp_dum
    type(contrib_cache_outer), allocatable, dimension(:) :: F, D, S
    type(contrib_cache_outer), allocatable, dimension(:) , optional :: Xf
    type(contrib_cache), allocatable, dimension(:) :: cache, lof_cache
    
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat,  get_2el_mat, get_xc_mat
    
    external :: out_print
    character(len=2047) :: out_str

    max_order = max(maxval(kn_rule(:, 1)), maxval(kn_rule(:, 2)))
    max_npert = maxval((/(p_tuples(i)%npert, i = 1, sum(n_freq_cfgs))/))

    ! Make dummy perturbation tuples to head cache to identify by perturbation order: Identifier is frequencies
    allocate(p_dummy_orders(max_npert))
    do i = 1, max_npert
       p_dummy_orders(i)%npert = 1
       allocate(p_dummy_orders(i)%pdim(1))
       allocate(p_dummy_orders(i)%plab(1))
       allocate(p_dummy_orders(i)%pid(1))
       allocate(p_dummy_orders(i)%freq(1))
       p_dummy_orders(i)%pdim = (/1/)
       p_dummy_orders(i)%plab = (/'DMMY'/)
       p_dummy_orders(i)%pid = (/1/)
       p_dummy_orders(i)%freq = (/1.0*i/)
    end do


    ! Check if this stage passed previously and if so, then retrieve and skip execution
    call prog_incr(prog_info, r_flag, 2)
        
    if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
          
       write(out_str, *) ' '
       call out_print(out_str, 1)
       write(out_str, *) 'SDF identification was completed in previous'
       call out_print(out_str, 1)
       write(out_str, *) 'invocation: Passing to next stage of calculation'
       call out_print(out_str, 1)
       write(out_str, *) ' '
       call out_print(out_str, 1)
          
!        allocate(cache)
!        call contrib_cache_retrieve(cache, 'OPENRSP_FDS_ID')
       lof_retrieved = .TRUE.
             
    else

       call contrib_cache_allocate(cache)
    
       ! Recurse to identify all necessary perturbed F, D, S
       k = 1 
       do i = 1, n_props
          do j = 1, n_freq_cfgs(i)
          
             len_cache = size(cache)
          
             call rsp_fds_recurse(p_tuples(k), kn_rule(k, :), max_npert, p_dummy_orders, len_cache, cache, out_print)
             k = k + 1

          end do
       end do

!        call contrib_cache_store(cache, r_flag, 'OPENRSP_FDS_ID')
       
    end if

    ! NOTE: Something may be wrong about this progress increase, revisit if problems
    call prog_incr(prog_info, r_flag, 2)
    
    lof_retrieved = .FALSE.
    rsp_eqn_retrieved = .FALSE.
    
    ! For each order of perturbation identified (lowest to highest):
    do i = 1, max_order
    
       lof_mem_total = 0
       
       ! Will be determined below, this to avoir maybe-uninitialized warning
       c_ord = 0
    
       ! Find entry for this order
       do m = 1, size(cache)
       
          if (cache(m)%p_inner%freq(1) == 1.0*i) then
          
             c_ord = m
          
             exit
             
          end if
       
       end do
       
       ! Contains number of components of perturbed matrices for each perturbation
       allocate(size_i(cache(c_ord)%num_outer))
       
       ! Traverse to set up size information
       k = 1
       do m = 1, size(cache(c_ord)%contribs_outer)
       
          ! Skip the dummy entry, who cares
          if (cache(c_ord)%contribs_outer(m)%dummy_entry) then
             cycle
          end if
          
          size_i(k) = cache(c_ord)%contribs_outer(m)%blks_tuple_triang_size(1)
          k = k + 1
       
       end do
       
       ! Check if this stage passed previously and if so, then retrieve and skip execution       
       if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
          
          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'Lower-order Fock matrix contribution identification at order', i
          call out_print(out_str, 1)
          write(out_str, *) 'was completed in previous invocation: Passing to next stage of calculation'
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
       
          if (.NOT.(lof_retrieved)) then
          
!              call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
       
          if (.NOT.(lof_retrieved)) then
    
             call contrib_cache_allocate(lof_cache)
             
          end if
          
          
          ! Traverse all elements of outer cache of present cache element
          k = 1
          do m = 1, size(cache(c_ord)%contribs_outer)
       
             ! Skip the dummy entry
             if (cache(c_ord)%contribs_outer(m)%dummy_entry) then
                cycle
             end if
          
             residue_select = .FALSE.
             residue_select = .NOT.find_complete_residualization( & 
                      cache(c_ord)%contribs_outer(m)%p_tuples(1)) &
                     .AND.find_residue_info(cache(c_ord)%contribs_outer(m)%p_tuples(1))

             len_lof_cache = size(lof_cache)
             
             ! Recurse to identify lower-order Fock matrix contributions
             ! The p_tuples attribute should always be length 1 here, so OK to take the first element
             call rsp_lof_recurse(cache(c_ord)%contribs_outer(m)%p_tuples(1), &
                                  cache(c_ord)%contribs_outer(m)%p_tuples(1)%npert, &
                                  1, (/get_emptypert()/), .TRUE., len_lof_cache, lof_cache, &
                                  1, (/Fp_dum/), out_print, residue_select)
             
             k = k + 1
       
          end do
          
!           call contrib_cache_store(lof_cache, r_flag, 'OPENRSP_LOF_CACHE')
       
       end if
       
       ! Check if this stage passed previously and if so, then retrieve and skip execution
       call prog_incr(prog_info, r_flag, 2)
       
       
       if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
          
          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'Lower-order Fock matrix contribution calculation at order', i
          call out_print(out_str, 1)
          write(out_str, *) 'was completed in previous invocation: Passing to next stage of calculation'
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
                    
          if (.NOT.(lof_retrieved)) then
          
!              call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
       
          ! Traverse lower-order Fock cache and precalculate elements
          k = 1
          len_lof_cache = size(lof_cache)
       
          do m = 1, len_lof_cache
       
             ! Check if this stage passed previously and if so, then retrieve and skip execution
             if (rs_check(prog_info, rs_info, r_flag, lvl=3)) then
             
                write(out_str, *) ' '
                call out_print(out_str, 1)
                write(out_str, *) 'Lower-order Fock matrix contribution calculation at order', i, 'inside traversal'
                call out_print(out_str, 1)
                write(out_str, *) 'was completed in previous invocation: Passing to next stage of calculation'
                call out_print(out_str, 1)
                write(out_str, *) ' '
                call out_print(out_str, 1)
           
                ! Note: No cache retrieval here: In order to get to this position, the
                ! cache would already have been retrieved
          
             else
       
                o_size = 0
                call rsp_lof_calculate(size(D), D, get_1el_mat, get_ovl_mat, get_2el_mat, out_print, &
                                       lof_cache(k), o_size, mem_mgr)

                if (mem_exceed(mem_mgr)) then
                
                   return
                
                end if
                
                lof_mem_total = lof_mem_total + lof_cache(k)%blks_triang_size * o_size
                
                if (.NOT.(mem_mgr%calibrate)) then
                
!                    call contrib_cache_store(lof_cache, r_flag, 'OPENRSP_LOF_CACHE')
                   
                end if
                
                 
                
             end if
                      
             call prog_incr(prog_info, r_flag, 3)
             k = k + 1
       
          end do
       
       end if
       
       ! Set 'retrieved' flag to false: If retrieved at current order, state should now 
       ! be as if not restarted (i.e. as if calculated in present run)
       lof_retrieved = .FALSE.

      
       ! Check if this stage passed previously and if so, then retrieve and skip execution
       call prog_incr(prog_info, r_flag, 2)
       
       
       if (rs_check(prog_info, rs_info, r_flag, lvl=2)) then
          
          write(out_str, *) ' '
          call out_print(out_str, 1)
          write(out_str, *) 'Perturbed S, D, F calculation at order', i, 'was completed'
          call out_print(out_str, 1)
          write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
          
          if (.NOT.(sdf_retrieved)) then
          
!              call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE.)
!              call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE.)
!              call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE.)
             
             if (present(Xf)) then
       
!              call contrib_cache_outer_retrieve(Xf, 'OPENRSP_Xf_CACHE', .FALSE.)
          
             end if
                       
          end if
          
          sdf_retrieved = .TRUE.
       
       else
       
          write(out_str, *) 'Calculating perturbed S, D, F at order', i
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
       
          len_lof_cache = size(lof_cache)
          len_d = size(D)
          len_curr_outer = size(cache(c_ord)%contribs_outer)
       
          ! Calculate all perturbed S, D, F at this order
       
          if (present(Xf)) then
          
             call rsp_sdf_calculate(len_curr_outer, cache(c_ord)%contribs_outer, cache(c_ord)%num_outer, size_i,&
                  get_rsp_sol, get_ovl_mat, get_2el_mat, get_xc_mat, out_print, len_d, F, D, S, len_lof_cache, lof_cache, &
                  rsp_eqn_retrieved, prog_info, rs_info, r_flag, mem_mgr, Xf=Xf)
               
               
          else
          
             call rsp_sdf_calculate(len_curr_outer, cache(c_ord)%contribs_outer, cache(c_ord)%num_outer, size_i,&
                  get_rsp_sol, get_ovl_mat, get_2el_mat, get_xc_mat, out_print, len_d, F, D, S, len_lof_cache, lof_cache, &
                  rsp_eqn_retrieved, prog_info, rs_info, r_flag, mem_mgr)
          
          
          end if
               
          if (mem_exceed(mem_mgr)) then
                
                   return
                
          end if
          
          if (.NOT.(mem_mgr%calibrate)) then
          
             len_d = size(D)
          
!              call contrib_cache_outer_store(len_d, S, 'OPENRSP_S_CACHE', r_flag)
!              call contrib_cache_outer_store(len_d, D, 'OPENRSP_D_CACHE', r_flag)
!              call contrib_cache_outer_store(len_d, F, 'OPENRSP_F_CACHE', r_flag)
             
          end if
       
       end if
       
       deallocate(size_i)
       
       call mem_decr(mem_mgr, lof_mem_total)
       
       call prog_incr(prog_info, r_flag, 2)
       
          
    end do
    
    deallocate(cache)
    
    if (max_order > 0) then
    
       deallocate(lof_cache)
       
    end if   
    
    deallocate(p_dummy_orders)
    
  end subroutine

  ! Recursive routine to identify necessary perturbed F, D, S
  recursive subroutine rsp_fds_recurse(pert, kn, max_npert, p_dummy_orders, &
                       len_cache, contribution_cache, out_print)

    implicit none

    integer :: n_props, max_npert, len_cache
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%npert) :: psub
    type(p_tuple), dimension(max_npert) :: p_dummy_orders
    integer, dimension(2) :: kn
    integer :: i, j, k
    type(contrib_cache), allocatable, dimension(:) :: contribution_cache
    external :: out_print
    character(len=2047) :: out_str
    
    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    if (pert%npert > 1) then

       call make_p_tuple_subset(pert, psub)
 
       do i = 1, size(psub)
       
          if (contrib_cache_already(len_cache, contribution_cache, 2, (/p_dummy_orders(psub(i)%npert), &
                              p_tuple_standardorder(psub(i))/)) .eqv. .FALSE.) then

             call rsp_fds_recurse(psub(i), kn, max_npert, p_dummy_orders, len_cache, contribution_cache, out_print)

          end if

       end do       

    end if

    ! See if already identified, if not and if keeping, then store element
    if (contrib_cache_already(len_cache, contribution_cache, 2, (/p_dummy_orders(pert%npert), &
                              p_tuple_standardorder(pert)/)) .eqv. .FALSE.) then
         
       if (kn_skip(pert%npert, pert%pid, kn) .eqv. .FALSE.) then

          write(out_str, *) 'Identified perturbed S, D, F for calculation with labels ', pert%plab, &
                            'and perturbation id ', pert%pid, ' with frequencies (real part)', &
                             real(pert%freq)
          call out_print(out_str, 2)
          write(out_str, *) ' '
          call out_print(out_str, 2)
                 
          k = 1

          do j = 1, pert%npert
             pert%pid(j) = k
             k = k + 1
          end do
          
          call contrib_cache_add_element(len_cache, contribution_cache, 2, (/p_dummy_orders(pert%npert), &
                                                                  p_tuple_standardorder(pert)/))

       end if

    end if

  end subroutine
  
  ! Recursive routine to identify (dryrun == .TRUE.) or assemble (dryrun == .FALSE.) 
  ! lower-order perturbed Fock matrix terms
  recursive subroutine rsp_lof_recurse(pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, dryrun, len_lof_cache, &
                       fock_lowerorder_cache, &
                       fp_size, Fp, out_print, residue_select)

    implicit none

    ! fp_size and Fp are dummy if this is a dryrun
    
    logical :: density_order_skip, residue_skip, dryrun, residue_select
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations
    integer :: fp_size
    integer :: len_lof_cache
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(contrib_cache), allocatable, dimension(:) :: fock_lowerorder_cache
    type(QcMat), dimension(fp_size) :: Fp
    
    external :: out_print
    character(len=2047) :: out_str

    
    ! If perturbation tuple not empty, recurse further
    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%npert == 0) then

          call rsp_lof_recurse(p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               dryrun, len_lof_cache, fock_lowerorder_cache, fp_size, &
               Fp, out_print,residue_select)

       else

          call rsp_lof_recurse(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
               p_tuples(2:size(p_tuples))/), dryrun, len_lof_cache, &
               fock_lowerorder_cache, fp_size, Fp, out_print,residue_select)

       end if
    
       ! 2. Differentiate all of the contraction densities in turn

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%npert == 0) then
             t_new(i) = p_tuple_getone(pert, 1)
          else
             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))
          end if

          call rsp_lof_recurse(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               t_new, dryrun, len_lof_cache, fock_lowerorder_cache, &
               fp_size, Fp, out_print,residue_select)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       if (num_p_tuples < 2) then
       
          call rsp_lof_recurse(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples + 1, &
               (/p_tuples(:), p_tuple_getone(pert, 1)/), dryrun, &
               len_lof_cache, fock_lowerorder_cache, fp_size, Fp, &
               out_print, residue_select)
               
       end if

    ! If recursion is done, proceed to cache manipulation stage
    else

       p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)
       density_order_skip = .FALSE.
       residue_skip = .FALSE.

       if (residue_select.and.find_residue_info(p_tuples(1))) residue_skip = .TRUE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%npert >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if ( (density_order_skip .EQV. .FALSE.).AND. &
            (residue_skip .EQV. .FALSE.) ) then
       
          if (contrib_cache_already(len_lof_cache, fock_lowerorder_cache, &
          num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples))) then

             ! If in cache and not dryrun, then put precalculated data in answer array
             if (.NOT.(dryrun)) then
             
                write(out_str, *) 'Getting lower-order Fock matrix contribution from cache:'
                call out_print(out_str, 2)
             
                do i = 1, num_p_tuples
                
                   if (i == 1) then

                      write(out_str, *) 'F', p_tuples(i)%plab
                      call out_print(out_str, 2)
          
                   else
                   
                      write(out_str, *) 'D', p_tuples(i)%plab
                      call out_print(out_str, 2)
                      
                   end if
                   
                end do
             
                call contrib_cache_getdata(len_lof_cache, fock_lowerorder_cache, &
                num_p_tuples, p_tuples, contrib_size=fp_size, ind_len=0, mat=Fp)
             
                write(out_str, *) ' '
                call out_print(out_str, 2)

             end if   
          
          else

             ! If dryrun and not in cache, then add identified element to cache
             ! for later calculation upon traversal
             if (dryrun) then
          
                write(out_str, *) 'Adding newly identified lower-order perturbed Fock matrix contribution to cache:'
                call out_print(out_str, 2)
                
                do i = 1, num_p_tuples
                   
                   if (i == 1) then
                   
                      write(out_str, *) 'F', p_tuples(i)%plab
                      call out_print(out_str, 2)
                   
                   else
                   
                      write(out_str, *) 'D', p_tuples(i)%plab
                      call out_print(out_str, 2)
                  
                   end if
                end do
                
                write(out_str, *) ' '
                call out_print(out_str, 2)
                
                call contrib_cache_add_element(len_lof_cache, fock_lowerorder_cache, &
                     num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples))
                
             else
             
                write(out_str, *) 'ERROR: Wanted lower-order perturbed Fock matrix contribution but it was not in cache'
                call out_print(out_str, -1)
             
             end if
          
          end if

       end if

    end if

  end subroutine
  
  ! Calculate lower-order Fock contributions for a given inner perturbation tuple
  subroutine rsp_lof_calculate(len_d, D, get_1el_mat, get_ovl_mat, get_2el_mat, out_print, &
                               cache, total_outer_size_1, mem_mgr)

    implicit none

    type(mem_manager) :: mem_mgr
    integer :: mctr, mcurr, miter, msize, octr, mem_track
    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n, s, p
    integer :: total_outer_size_1, c1_ctr, lhs_ctr_1, num_pert
    integer :: num_0, num_1, tot_num_pert, offset
    integer :: len_d
    character(30) :: mat_str, fmt_str
    type(contrib_cache) :: cache
    type(contrib_cache_outer), dimension(len_d) :: D
    
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, contrib_0, contrib_1
    type(QcMat) :: D_unp
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    integer, allocatable, dimension(:,:) :: blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
    external :: get_1el_mat, get_ovl_mat, get_2el_mat
    
    external :: out_print
    character(len=2047) :: out_str
    
    write(out_str, *) 'Calculating lower-order Fock matrix contribution for inner tuple', cache%p_inner%plab
    call out_print(out_str, 1)
    write(out_str, *) ' '
    call out_print(out_str, 1)
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)

    allocate(outer_contract_sizes_1(cache%num_outer))    
   
    call mem_incr(mem_mgr, 1)
   
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QCMatInit(D_unp)

       call contrib_cache_getdata_outer(len_d, D, 1, (/get_emptypert()/), .FALSE., &
            contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)

    end if
   
    total_outer_size_1 = 0
    num_1 = 0
    k = 1
   
    ! Traversal: Find number of density matrices for contraction for 1-el, 2-el cases
    do m = 1, size(cache%contribs_outer)
    
       if (cache%contribs_outer(m)%dummy_entry) then
       
          cycle
    
       end if
       
       ! No chain rule application: Both 1-el and 2-el contributions
       if (cache%contribs_outer(m)%num_dmat == 0) then

          outer_contract_sizes_1(k) = 1
          total_outer_size_1 = total_outer_size_1 + 1
       
       ! One chain rule application: Only 2-el contributions
       else if (cache%contribs_outer(m)%num_dmat == 1) then
       
          num_1 = num_1 + 1
          outer_contract_sizes_1(k) = cache%contribs_outer(m)%blks_tuple_triang_size(1)
          total_outer_size_1 = total_outer_size_1 + cache%contribs_outer(m)%blks_tuple_triang_size(1)

       end if
       
       ! Check if enough memory
       if (.NOT.(mem_enough(mem_mgr, cache%blks_triang_size*outer_contract_sizes_1(k)))) then
       
          ! Not possible to run; flag it and return
          call mem_set_status(mem_mgr, 2)
          return
          
       end if
       
       call mem_incr(mem_mgr, cache%blks_triang_size*outer_contract_sizes_1(k))
       allocate(cache%contribs_outer(m)%data_mat(cache%blks_triang_size*outer_contract_sizes_1(k)))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          do i = 1, cache%blks_triang_size*outer_contract_sizes_1(k)
    
             call QcMatInit(cache%contribs_outer(m)%data_mat(i), D_unp)
             call QcMatZero(cache%contribs_outer(m)%data_mat(i))
                
          end do
          
       end if
       
   
       if (cache%contribs_outer(m)%num_dmat == 0) then
       
          write(out_str, *) 'All inner contribution'
          call out_print(out_str, 1)
          
       else
          
          write(out_str, *) 'Outer contribution:'
          call out_print(out_str, 1)
          
       end if
       
       do i = 1, cache%contribs_outer(m)%num_dmat
          
          write(out_str, *) 'D', cache%contribs_outer(m)%p_tuples(i)%plab
          call out_print(out_str, 1)
       
       end do
       
       write(out_str, *) ' '
       call out_print(out_str, 1)
    
       k = k + 1
    
    end do
    

  
    ! Memory savings loop start
    ! Additional needed for non-savings: size_inner * (2 + total_outer_size_1) (+ additions from callback)
    ! Savings: if > size_inner*2: Do "all inner" first, then appropriate number of outer
    ! If < size_inner*2, do subset of inner for all inner first, then in the same way for outer
    ! Minimum needed: 2 or 3 (deallocate/not deallocate contraction mat before data storage)
    
    mem_track = total_outer_size_1 + 1 + &
              cache%blks_triang_size*(1 + total_outer_size_1)
    
    if (mem_enough(mem_mgr, mem_track)) then
    
       ! Possible to run in non-savings mode
       mem_track = total_outer_size_1

    else if (mem_enough(mem_mgr, 2 + 2 * cache%blks_triang_size)) then

       ! Possible to run in savings mode
       call mem_set_status(mem_mgr, 1)
       
       ! Find maximum number of components to do at the same time
       mem_track = (mem_mgr%remain - 2)/(2 * cache%blks_triang_size)
              
    else
    
       ! Not possible to run; flag it and return
       call mem_set_status(mem_mgr, 2)
       return
    
    end if
    

    
    mcurr = 1
    
    do while (mcurr <= total_outer_size_1)

       msize = min(mem_track, total_outer_size_1 - mcurr + 1)
       mctr = 0
       num_0 = 0
    
       ! Initializing data and arrays for external calls
       
       call mem_incr(mem_mgr, msize)
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(LHS_dmat_1(msize))
       
          do i = 1, size(LHS_dmat_1)
       
            call QCMatInit(LHS_dmat_1(i), D_unp)
       
          end do
          
       end if
       
       call mem_incr(mem_mgr, cache%blks_triang_size)
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(contrib_0(cache%blks_triang_size))
          
       end if
       
       call mem_incr(mem_mgr, cache%blks_triang_size*msize)
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(contrib_1(cache%blks_triang_size*msize))
          
       end if
       
       if (.NOT.(mem_mgr%calibrate)) then       
          do i = 1, size(contrib_0)
            call QCMatInit(contrib_0(i), D_unp)
            call QCMatZero(contrib_0(i))
          end do
       
          do i = 1, size(contrib_1)
            call QCMatInit(contrib_1(i), D_unp)
            call QCMatZero(contrib_1(i))
          end do
      
       end if
       
       k = 1
       lhs_ctr_1 = 1
       
       ! Traversal: Get matrices for contraction from cache
       do m = 1, size(cache%contribs_outer)
    
          if (cache%contribs_outer(m)%dummy_entry) then
       
             cycle
    
          end if
       

       
          ! One chain rule application
          if (cache%contribs_outer(m)%num_dmat == 1) then
          
             do n = 1, outer_contract_sizes_1(k) 
             
                if ((lhs_ctr_1 + n - 1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                   if (.NOT.(mem_mgr%calibrate)) then
                
                      call contrib_cache_getdata_outer(size(D), D, 1, &
                           (/cache%contribs_outer(m)%p_tuples(1)/), .FALSE., &
                           contrib_size=1, ind_len=1, &
                           ind_unsorted=cache%contribs_outer(m)%indices(n, :), &
                           mat_sing=LHS_dmat_1(mctr + 1))
                           
                   end if
                        
                        mctr = mctr + 1
                
                end if
       
             end do
             
          ! No chain rule applications: Contraction matrix is only unperturbed D
          elseif(cache%contribs_outer(m)%num_dmat == 0) then
             
             if ((lhs_ctr_1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                num_0 = 1
                
                if (.NOT.(mem_mgr%calibrate)) then
                
                   call contrib_cache_getdata_outer(size(D), D, 1, (/get_emptypert()/), .FALSE., &
                        contrib_size=1, ind_len=1, ind_unsorted=(/1/), &
                        mat_sing=LHS_dmat_1(mctr + 1))
                        ! Was: mat_sing=LHS_dmat_1(1)), change if problems
                        
                end if
                     
                mctr = mctr + 1
                
             end if
             
          end if
          
          
          lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
          k = k + 1
       
       
       end do
       
       ! Calculate contributions
       
       ! Calculate one-electron contributions
       if (num_0 > 0) then
       
          if (.NOT.(mem_mgr%calibrate)) then
          
             write(out_str, *) 'Calculating density matrix-independent contribution'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
       
             call get_1el_mat(num_pert, pert_ext, size(contrib_0), contrib_0)
             
          end if
         
          t_matrix_bra = get_emptypert()
          t_matrix_ket = get_emptypert()
         
   ! T matrix terms are not reintroduced yet, skipping for now      
         
   !        call rsp_ovlave_t_matrix(get_ovl_mat, cache%p_inner, cache%p_inner%npert, &
   !                                 t_matrix_bra, t_matrix_ket, 1, &
   !                                 D_unp, size(contrib_0), contrib_0)
       
       end if
       
       ! Calculate two-electron contributions
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          write(out_str, *) 'Calculating first-order density matrix-dependent contribution'
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
       
          call get_2el_mat(num_pert, pert_ext, size(LHS_dmat_1), LHS_dmat_1, &
                           size(contrib_1), contrib_1)
                           
       end if
       
       call mem_decr(mem_mgr, msize)
       
       if (.NOT.(mem_mgr%calibrate)) then

          do i = 1, msize
          
             call QcMatDst(LHS_dmat_1(i))
             
          end do
        
          deallocate(LHS_dmat_1)                        
          
       end if
                          
       
       ! Traversal: Add 1-el and two-el contributions together
       
       k = 1
       c1_ctr = 1
       mctr = 0
       octr = 1
       
       do m = 1, size(cache%contribs_outer)
    
          if (cache%contribs_outer(m)%dummy_entry) then
       
             cycle
    
          end if

          ! One-el and two-el contributions ("all inner contribution")
          if (cache%contribs_outer(m)%num_dmat == 0) then
          
             ! For "all inner contribution", the data is already ordered correctly
             if ((octr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                if (.NOT.(mem_mgr%calibrate)) then
             
                   do i = 1, cache%blks_triang_size
             
                      call QcMatkAB(1.0d0, contrib_0(i), contrib_1(c1_ctr + i - 1), &
                      cache%contribs_outer(m)%data_mat(i))
          
                   end do
                   
                end if
                
                c1_ctr = c1_ctr + cache%blks_triang_size
                mctr = mctr + 1
                
             end if
             
             octr = octr + 1
     
          ! Only two-el contribution
          else if (cache%contribs_outer(m)%num_dmat == 1) then
                     
             ! Initialize block information for cache indexing
             tot_num_pert = cache%p_inner%npert + &
             sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1, &
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
          
             do i = 1, size(cache%contribs_outer(m)%indices, 1)
             
                if ((octr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                   do j = 1, size(cache%indices, 1)
                   
                      offset = get_triang_blks_tuple_offset(cache%contribs_outer(m)%num_dmat + 1, &
                      cache%p_inner%npert + sum((/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1,&
                      cache%contribs_outer(m)%num_dmat)/)), &
                      (/cache%nblks, (/(cache%contribs_outer(m)%nblks_tuple(n), n = 1, &
                      cache%contribs_outer(m)%num_dmat) /) /), &
                      (/cache%p_inner%npert, (/(cache%contribs_outer(m)%p_tuples(n)%npert, n = 1, &
                      cache%contribs_outer(m)%num_dmat)/)/), &
                      blks_tuple_info, &
                      blk_sizes, &
                      (/cache%blks_triang_size, &
                      (/(cache%contribs_outer(m)%blks_tuple_triang_size(n), &
                      n = 1, cache%contribs_outer(m)%num_dmat)/)/), &
                      (/cache%indices(j, :), cache%contribs_outer(m)%indices(i, :)/))

                      if (.NOT.(mem_mgr%calibrate)) then
                      
                         ! Store result in cache
                         call QcMatRAXPY(1.0d0, contrib_1(c1_ctr + j - 1), &
                         cache%contribs_outer(m)%data_mat(offset))
                         
                      end if
                   
                   end do
                
                   c1_ctr = c1_ctr + cache%blks_triang_size
                   mctr = mctr + 1
                
                end if
                
                octr = octr + 1
           
             end do
                
             deallocate(blk_sizes)
             deallocate(blks_tuple_info)
                      
          end if
          
          k = k + 1
          
       end do
       
       mcurr = mcurr + msize
    
       if (.NOT.(mem_mgr%calibrate)) then
    
          do i = 1, size(contrib_0)
             call QcMatDst(contrib_0(i))
          end do
       
          deallocate(contrib_0)
          
       end if
       
       call mem_decr(mem_mgr, cache%blks_triang_size)

       
       if (.NOT.(mem_mgr%calibrate)) then
       
          do i = 1, size(contrib_1)
             call QcMatDst(contrib_1(i))
          end do
       
          deallocate(contrib_1)
          
       end if
       
       call mem_decr(mem_mgr, cache%blks_triang_size*msize)
       
    end do
       
    deallocate(outer_contract_sizes_1)

    if (.NOT.(mem_mgr%calibrate)) then
    
       call QCMatDst(D_unp)

    end if
    
    call mem_decr(mem_mgr, 1)
   
  end subroutine

  
  ! Do main part of perturbed S, D, F calculation at one order
  subroutine rsp_sdf_calculate(len_cache_outer, cache_outer, num_outer, size_i, &
  get_rsp_sol, get_ovl_mat, get_2el_mat, get_xc_mat, out_print, len_d, F, D, S, &
  len_lof_cache, lof_cache, rsp_eqn_retrieved, prog_info, rs_info, r_flag, mem_mgr, Xf)
  
    implicit none
    
    type(mem_manager) :: mem_mgr
    integer :: mctr, mcurr, miter, msize, octr, r_flag, mem_track
    logical :: termination, rsp_eqn_retrieved, residue_select, residualization
    logical :: init_xx
    integer :: num_outer, ind_ctr, npert_ext, sstr_incr, superstructure_size
    integer :: i, j, k, m, n, w, nblks
    integer :: first, last, ierr
    integer, dimension(0) :: noc
    integer, dimension(3) :: prog_info, rs_info
    integer, allocatable, dimension(:) :: pert_ext, blk_sizes, ind
    integer, allocatable, dimension(:,:) :: blk_info
    integer, allocatable, dimension(:,:) :: indices
    integer, dimension(num_outer) :: size_i
    character(30) :: mat_str, fmt_str
    complex(8), dimension(num_outer) :: freq_sums
    complex(8) :: xrtm
    real(8), parameter :: xtiny=1.0d-8
    type(p_tuple) :: pert, pert_xc_null
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    
    integer :: len_cache_outer, len_d, len_lof_cache
    
    type(contrib_cache), allocatable, dimension(:) :: lof_cache
    type(contrib_cache_outer), dimension(len_cache_outer) :: cache_outer
    type(contrib_cache_outer), allocatable, dimension(:) :: F, D, S
    type(contrib_cache_outer), optional, allocatable, dimension(:) :: Xf

    type(Qcmat), allocatable, dimension(:) :: Dh, Dp, Fp, Sp, RHS, X, Xx
    type(Qcmat) :: A, B, C, T, U
    external :: get_rsp_sol, get_ovl_mat,  get_2el_mat, get_xc_mat
    external :: out_print
    character(len=2047) :: out_str

    write(out_str, *) 'Number of different perturbation tuples at this order:', num_outer
    call out_print(out_str, 1)
    
    write(out_str, *) ' '
    call out_print(out_str, 1)
    

    ! Set up null perturbation for XC use
    
    pert_xc_null%npert = 1
    allocate(pert_xc_null%pdim(1))
    allocate(pert_xc_null%plab(1))
    allocate(pert_xc_null%freq(1))
    allocate(pert_xc_null%pid(1))
    
    pert_xc_null%plab(1) = 'NULL'
    pert_xc_null%freq(1) = 0.0
    pert_xc_null%pid(1) = 1
    
    ! May be inaccurate, revisit if problems
    if (.NOT.(mem_enough(mem_mgr, 8 * sum(size_i) + 9))) then
    
       ! Not possible to run; flag it and return
       call mem_set_status(mem_mgr, 2)
       return
    
    end if
    
    
    ! Initialize matrices

    call mem_incr(mem_mgr, 5)
   
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QcMatInit(A)
       call QcMatInit(B)
       call QcMatInit(C)
       call QcMatInit(T)
       call QcMatInit(U)
    
       call contrib_cache_getdata_outer(len_d, D, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=A)
       call contrib_cache_getdata_outer(len_d, S, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=B)
       call contrib_cache_getdata_outer(len_d, F, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=C)
    
    end if
    
    init_xx = .FALSE.
    
    do m = 1, size(cache_outer)
       
          if (.NOT.(cache_outer(m)%dummy_entry)) then
       
             if (cache_outer(m)%p_tuples(1)%do_residues.gt.0) then
       
                if (.NOT.(init_xx)) then
       
                   call mem_incr(mem_mgr, 1)

                   allocate(Xx(1))
                   call QcMatInit(Xx(1))
                   call contrib_cache_getdata_outer(len_d, Xf, 1, (/get_emptypert()/), &
                   .FALSE., contrib_size = 1, &
                   ind_len = 1, ind_unsorted = (/1/), mat_sing = Xx(1))
                   
                   init_xx = .TRUE.
                
                end if
          
             end if
             
          end if   
          
    end do   
    
    
     
    ! Traverse cache elements and calculate pertubed S
    ind_ctr = 1
    k = 1
        
    do m = 1, size(cache_outer)
    
       if (cache_outer(m)%dummy_entry) then
       
          cycle
    
       end if
       
       pert = cache_outer(m)%p_tuples(1)

       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(Sp(size_i(k)))
          
          do i = 1, size_i(k)
          
             call QcMatInit(Sp(i), A)
             call QcMatZero(Sp(i))
     
          end do
       

          write(out_str, *) 'Calculating perturbed overlap matrix for perturbation tuple with labels'
          call out_print(out_str, 1)
          write(out_str, *) pert%plab
          call out_print(out_str, 1)
          write(out_str, *) ''
          call out_print(out_str, 1)
       
          ! For each cache element:
          ! Calculate Sp
          call p_tuple_to_external_tuple(pert, npert_ext, pert_ext)
          call get_ovl_mat(0, noc, 0, noc, npert_ext, pert_ext, &
                        size_i(k), Sp)
          
          
          deallocate(pert_ext)
          
       end if
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
          
          len_d = size(S)
          
          call contrib_cache_outer_add_element(len_d, S, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Sp )
          
          do i = 1, size_i(k)
          
             call QcMatDst(Sp(i))
     
          end do
          
          deallocate(Sp)
          
       end if
       
       call mem_decr(mem_mgr, size_i(k))
       
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
               
    end do
    
    
    call mem_incr(mem_mgr, 5 * sum(size_i))

    
    
    
    if (.NOT.(mem_mgr%calibrate)) then    
    
       allocate(Dh(sum(size_i)))
       allocate(Dp(sum(size_i)))
       allocate(Fp(sum(size_i)))
       allocate(RHS(sum(size_i)))
       allocate(X(sum(size_i)))
    
       do i = 1, sum(size_i)
    
          call QcMatInit(Dh(i), A)
          call QcMatZero(Dh(i))
          call QcMatInit(Dp(i), A)
          call QcMatZero(Dp(i))
          call QcMatInit(Fp(i), A)
          call QcMatZero(Fp(i))
          call QcMatInit(RHS(i), A)
          call QcMatZero(RHS(i))
          call QcMatInit(X(i), A)
          call QcMatZero(X(i))
    
       end do
    
    end if
    
    
 
 
    ! Traverse cache elements and do various stages before collective 2-el
    ! particular Fock matrix contribution call
    ind_ctr = 1
    k = 1

    do m = 1, size(cache_outer)
    
       if (cache_outer(m)%dummy_entry) then
       
          cycle
    
       end if
       
       pert = cache_outer(m)%p_tuples(1)
       
       ! Get frequency sum
       freq_sums(k) = sum(real(pert%freq(:))) 
       
       ! Set up block info
       nblks = get_num_blks(pert)

       allocate(blk_info(nblks, 3))
       allocate(blk_sizes(pert%npert))
       blk_info = get_blk_info(nblks, pert)
       blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))
       
       ! Do we treat a perturbation which requires residue selection of terms?
       residue_select = .not.find_complete_residualization(pert).and.find_residue_info(pert)
 
       ! Add the initialized Dp to cache
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          write(out_str, *) 'Getting lower-order Fock matrix terms for perturbation tuple with labels'
          call out_print(out_str, 1)
          write(out_str, *) pert%plab
          call out_print(out_str, 1)
          write(out_str, *) ''
          call out_print(out_str, 1)
       
          len_d = size(D)
       
          call contrib_cache_outer_add_element(len_d, D, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )

          ! Assemble Fp (lower-order) for all components and add to cache
 
          call rsp_lof_recurse(pert, pert%npert, &
                               1, (/get_emptypert()/), .FALSE., size(lof_cache), lof_cache, size_i(k), &
                               Fp(ind_ctr:ind_ctr + size_i(k) - 1), out_print, &
                               residue_select = residue_select)
                               
          ! XC call should go here
          
          write(out_str, *) 'Calculating intermediate exchange-correlation contributions'
          call out_print(out_str, 1)
          write(out_str, *) 'to perturbed Fock matrix for perturbation tuple with labels'
          call out_print(out_str, 1)
          write(out_str, *) pert%plab
          call out_print(out_str, 1)
          write(out_str, *) ' '
          call out_print(out_str, 1)
          
          
          len_d = size(D)
          
          ! Currently only one freq. configuration
          ! The (k,n) rule argument is adapted for the Fock contribution case
          call rsp_xc_wrapper(1, (/pert/), (/pert%npert, pert%npert/), len_d, D, get_xc_mat, out_print, &
                                size_i(k), mem_mgr, fock=Fp(ind_ctr:ind_ctr + size_i(k) - 1))
          
       
       end if
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          len_d = size(F)
       
          call contrib_cache_outer_add_element(len_d, F, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
       end if
       
       
       ! Calculate Dp for all components and add to cache
       
       ! Set up Z matrix calculation
       superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%npert, pert%npert/), &
                          .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

       sstr_incr = 0
       
       allocate(derivative_structure(superstructure_size, 3))
       allocate(indices(size_i(k), pert%npert))
       allocate(ind(pert%npert))
       
       call derivative_superstructure(pert, &
            (/pert%npert, pert%npert/), .FALSE., &
            (/get_emptypert(pert), get_emptypert(pert), get_emptypert(pert)/), &
            superstructure_size, sstr_incr, derivative_structure)
      
       call make_triangulated_indices(nblks, blk_info, size_i(k), indices)
       
       ! Calculate Z matrices and process to get Dp
       do j = 1, size(indices, 1)
       
          ind = indices(j, :)
          
          call mem_incr(mem_mgr, 4)
          
          if (.NOT.(mem_mgr%calibrate)) then

             call rsp_get_matrix_z(superstructure_size, derivative_structure, &
                  (/pert%npert,pert%npert/), pert%npert, &
                  (/ (m, m = 1, pert%npert) /), pert%npert, &
                  ind, size(F), F, size(D), D, size(S), S, Dp(ind_ctr + j - 1), &
                  select_terms_arg=residue_select)

          end if
          
          call mem_decr(mem_mgr, 4)               
       
          if (.NOT.(mem_mgr%calibrate)) then
       
             call QcMatkABC(-1.0d0, Dp(ind_ctr + j - 1), B, A, T)
             call QcMatkABC(-1.0d0, A, B, Dp(ind_ctr + j - 1), U)
             call QcMatRAXPY(1.0d0, T, U)
             call QcMatRAXPY(1.0d0, U, Dp(ind_ctr + j - 1))
             
          end if
       
       end do
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          len_d = size(D)
       
          call contrib_cache_outer_add_element(len_d, D, .FALSE., 1, & 
               (/pert/), data_size = size_i(k),  data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
               
       end if
       
       ! Clean up and set up next iteration
       deallocate(indices)
       deallocate(ind)
       deallocate(derivative_structure)
       deallocate(blk_sizes)
       deallocate(blk_info)
              
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
    end do
    
       
    
! ! Debug printing kept for later use    
!     do i = 1, size(Dp)
!     
!           if (i < 10) then
!              fmt_str = "(A3, I1)"
!           else if (i < 100) then
!              fmt_str = "(A3, I2)"
!           else
!              fmt_str = "(A3, I3)"
!           end if
!           
!           write(mat_str, fmt_str) 'Dp_', i
!           write(*,*) 'i', i
!           write(*,*) 'fname:', mat_str
!           
!           j = QcMatWrite_f(Dp(i), trim(mat_str), ASCII_VIEW)
!     
!     end do
    
    

    
!   Kept for possible reintroduction: Calculate for subset of all contractions
!   if necessary because of memory considerations (now always all contractions)


!     k = 129
!     
!     if (size(Dp) > k) then
!     
!        do i = 1, size(Dp)/k + 1
!     
!           if (.NOT.((i - 1) *k >= size(Dp))) then
!     
!              first = (i - 1) * k + 1
!              last = min(i * k, size(Dp))
!     
!              call get_2el_mat(0, noc, last - first + 1, Dp(first:last), last - first + 1, Fp(first:last))
!     
!           end if
!        
!        end do
!   
!     else

    write(out_str, *) 'Calculating particular contribution remainders'
    call out_print(out_str, 1)
    write(out_str, *) ''
    call out_print(out_str, 1)
   
    if (.NOT.(mem_mgr%calibrate)) then
    
       ! Outside traversal:
       ! Complete Fp using Dp
       call get_2el_mat(0, noc, sum(size_i), Dp, sum(size_i), Fp)
       
     
       ! Set up null perturbation dimensionality
       pert_xc_null%pdim(1) = sum(size_i)
    
       call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), size(D), D, get_xc_mat, out_print, &
                                sum(size_i), mem_mgr, null_dmat=Dp, fock=Fp)
    
    
    end if
    
!     end if

    ! Traversal: Adding Fp to cache, constructing RHS
    
    ind_ctr = 1
    k = 1
    
    do m = 1, size(cache_outer)
    
       if (cache_outer(m)%dummy_entry) then
       
          cycle
    
       end if
       
       pert = cache_outer(m)%p_tuples(1)
       
       ! Do we treat a perturbation which requires residue selection of terms?
       residue_select = .not.find_complete_residualization(pert).and.find_residue_info(pert)

       ! Set up block info
       nblks = get_num_blks(pert)

       allocate(blk_info(nblks, 3))
       allocate(blk_sizes(pert%npert))
       blk_info = get_blk_info(nblks, pert)
       blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

       if (.NOT.(mem_mgr%calibrate)) then
       
         len_d = size(F)
        
         ! Add the completed Fp to cache
           call contrib_cache_outer_add_element(len_d, F, .FALSE., 1, & 
               (/pert/), data_size = size_i(k),  data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
      
       end if
       
       ! Construct right-hand side for all components
       
       ! Set up Y matrix calculation
       superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%npert, pert%npert/), &
                          .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

       sstr_incr = 0
       
       allocate(derivative_structure(superstructure_size, 3))
       allocate(indices(size_i(k), pert%npert))
       allocate(ind(pert%npert))
       
       call derivative_superstructure(pert, &
            (/pert%npert, pert%npert/), .FALSE., &
            (/get_emptypert(pert), get_emptypert(pert), get_emptypert(pert)/), &
            superstructure_size, sstr_incr, derivative_structure)
      
       call make_triangulated_indices(nblks, blk_info, size_i(k), indices)
       
       ! Calculate RHS matrices
       do j = 1, size(indices, 1)
       
          ind = indices(j, :)
          
          call mem_incr(mem_mgr, 4)
          if (.NOT.(mem_mgr%calibrate)) then
    
             call rsp_get_matrix_y(superstructure_size, derivative_structure, &
                  pert%npert, (/ (m, m = 1, pert%npert) /), &
                  pert%npert, ind, size(F), F, size(D), D, size(S), S, RHS(ind_ctr + j - 1), &
                  select_terms_arg = residue_select)
                  
          end if        
                  
          call mem_decr(mem_mgr, 4)          
                
       end do

       ! Clean up and set up next iteration
       deallocate(indices)
       deallocate(ind)
       deallocate(derivative_structure)
       deallocate(blk_sizes)
       deallocate(blk_info)
       
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
    end do
    
    ! Outside of traversal:
    ! Solve all response equations (opportunities for optimization for e.g.
    ! first-order EL, but that can be introduced later)
    
    ! IN NEW FORMAT, AWAITING CHANGES TO CALLBACK STRUCTURE
!     call get_rsp_sol(num_outer, size_i, (/(i - i + 1, i = 1, num_outer)/), freq_sums, RHS, X)

    m  = 1
    ind_ctr = 1
    k = 1
    
    ! Traverse to solve response equations
    ! May extend to do collectively all cache elements for which the
    ! frequency sum is the same
    
    do n = 1, size(cache_outer)
    
       if (cache_outer(n)%dummy_entry) then
       
          cycle
    
       end if
       
       ! MaR: I added this line because otherwise, looks like pert from Daniel's
       ! residualization flag determination below will not refer to the correct perturbation
       ! revisit this if problems
       pert = cache_outer(n)%p_tuples(1)
       
       if (pert%do_residues.gt.0) then       
        
           residualization = dabs(dabs(dble(freq_sums(k)))-dabs(dble(pert%exenerg(1)))).lt.xtiny
        
        else
        
           residualization = .false.
        
        end if

        if (size_i(k) > m) then
    
          do i = 1, size_i(k)/m + 1
    
             if (.NOT.((i - 1) * m >= size_i(k))) then
    
                first = (i - 1) * m + 1
                last = min(i * m, size_i(k))
                
                if (rs_check(prog_info, rs_info, r_flag, lvl=3)) then
                
                   write(out_str, *) 'Response equation solution batch was completed'
                   call out_print(out_str, 1)
                   write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
                   call out_print(out_str, 1)
                   write(out_str, *) ' '
                   call out_print(out_str, 1)
                   
                   if (.NOT.(rsp_eqn_retrieved)) then
                   
                      write(out_str, *) 'Retrieving response equation solutions from storage'
                      call out_print(out_str, 2)
                      write(out_str, *) ' '
                      call out_print(out_str, 2)
                   
!                       call mat_scal_retrieve(rs_info(3), 'OPENRSP_MAT_RSP', mat=X(1:rs_info(3)))
                      rsp_eqn_retrieved = .TRUE.
                      
             
                   end if

                else

                   if (.NOT.(mem_mgr%calibrate)) then
                   
                      if (.not.residualization) then
                      
                         write(out_str, *) 'Solving response equations'
                         call out_print(out_str, 1)
                         write(out_str, *) 'Frequency sum:', freq_sums(k)
                         call out_print(out_str, 2)
                         write(out_str, *) ' '
                         call out_print(out_str, 1)
                   
                         call get_rsp_sol(1,                                    &
                                          (/last-first+1/),                     &
                                          (/1/),                                &
                                          dcmplx(real((/freq_sums(k)/)),0.0d0), &
                                          RHS(ind_ctr+first-1:ind_ctr+last-1),  &
                                          X(ind_ctr+first-1:ind_ctr+last-1))

                      else

                         ! DaF: Replace LES solution by contraction for residualized perturbations
                         ! DaF: This is only for debugging! In principle we replace X by Xx.  
                         do j = first, last

                            ierr = QcMatDuplicate_f(Xx(1),COPY_PATTERN_AND_VALUE,X(ind_ctr+j-1))                                                
                            call QcMatTraceATrB(RHS(ind_ctr+j-1),Xx(1),xrtm)                                                                                        
                            write(*,*)'xrtm=',xrtm
                            write(*,*)'Nullifying Fp!'
                            call QcMatZero(Fp(ind_ctr+j-1))

                         end do
  
                      end if
                  
!                       call mat_scal_store(last - first + 1, 'OPENRSP_MAT_RSP', r_flag, &
!                            mat=X(ind_ctr+first-1:ind_ctr+last-1), start_pos = ind_ctr+first-1)
                    
                   end if
                   
                end if
                
                call prog_incr(prog_info, r_flag, 3)
                
             end if
       
          end do
       
        else
       
          ! Check if this stage passed previously and if so, then retrieve and skip execution
          if (rs_check(prog_info, rs_info, r_flag, lvl=3)) then
          
             write(out_str, *) 'Response equation solution batch was completed'
             call out_print(out_str, 1)
             write(out_str, *) 'in previous invocation: Passing to next stage of calculation'
             call out_print(out_str, 1)
             write(out_str, *) ' '
             call out_print(out_str, 1)
         
             if (.NOT.(rsp_eqn_retrieved)) then
          
!                 call mat_scal_retrieve(rs_info(3), 'OPENRSP_MAT_RSP', mat=X(1:rs_info(3)))
                rsp_eqn_retrieved = .TRUE.
             
             end if

          else
       
             if (.NOT.(mem_mgr%calibrate)) then
             
                write(out_str, *) 'Solving response equations'
                call out_print(out_str, 1)
                write(out_str, *) ''
                call out_print(out_str, 1)
       
                if (residualization) then 
                
                   write(out_str, *) 'ERROR: No residualization yet for 2nd solver call in rsp_sdf_calc'
                   call out_print(out_str, 0)
                   write(out_str, *) 'Cannot proceed with calculation, halting'
                   call out_print(out_str, -1)
                
                   stop 
                
                end if    
            
                call get_rsp_sol(1,                                    &
                                 (/size_i(k)/),                        &
                                 (/1/),                                &
                                 dcmplx(real((/freq_sums(k)/)),0.0d0), &
                                 RHS(ind_ctr:ind_ctr+size_i(k)-1),     &
                                 X(ind_ctr:ind_ctr+size_i(k)-1))
          
!                 call mat_scal_store(size_i(k), 'OPENRSP_MAT_RSP', r_flag, &
!                            mat=X(ind_ctr:ind_ctr+size_i(k)-1), start_pos = ind_ctr)
                           
             end if
                   
          end if
          
          call prog_incr(prog_info, r_flag, 3)
    
        end if
  
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
    end do
    
        
    if (.NOT.(mem_mgr%calibrate)) then
       
       do i = 1, size(RHS)
          call QcmatDst(RHS(i))
       end do
       deallocate(RHS)
       
    end if
    
    call mem_decr(mem_mgr, sum(size_i))
        
    ind_ctr = 1
    k = 1
            
    ! Traverse: Make homogeneous contribution to perturbed D
    do m = 1, size(cache_outer)
    
       if (cache_outer(m)%dummy_entry) then
       
          cycle
    
       end if
       
       ! For each cache element:
       ! Construct Dh for all components
    
       if (.NOT.(mem_mgr%calibrate)) then
    
          do j = 1, size_i(k)
       
             call QcMatkABC(-1.0d0, X(ind_ctr + j - 1), B, A, T)
             call QcMatkABC(1.0d0, A, B, X(ind_ctr + j - 1), U)
             call QcMatRAXPY(1.0d0, T, U)
             call QcMatAEqB(Dh(ind_ctr + j - 1), U)
     
          end do
          
       end if
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
    end do
       
        
    if (.NOT.(mem_mgr%calibrate)) then
    
       do i = 1, size(X)
          call QcmatDst(X(i))
       end do
       deallocate(X)
    
    end if
    
    call mem_decr(mem_mgr, sum(size_i))
    
    
    ! Outside traversal
    ! Calculate Fh for all components using Dh
    
    ! Kept for possible reintroduction: Call for subset of all contributions
    
!     if (.NOT.(mem_mgr%calibrate)) then
!     
!        ! Currently limiting number of simultaneous calls
!        ! Will later be tuned according to memory requirements
!        k = 65536
!        
!        if (size(Dh) > k) then
! 
!           do i = 1, size(Dh)/k + 1
!           
!              if (.NOT.((i - 1) *k >= size(Dh))) then
!        
!                 first = (i - 1) * k + 1
!                 last = min(i * k, size(Dh))
!        
!                 call get_2el_mat(0, noc, last - first + 1, Dh(first:last), &
!                      last - first + 1, Fp(first:last))
!                      
!                 ! Set up null perturbation dimensionality
!                 pert_xc_null%pdim(1) = last - first + 1
!                 
!                 call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), size(D), D, get_xc_mat, out_print, &
!                                     last - first + 1, mem_mgr, &
!                                     null_dmat=Dh(first:last), fock=Fp(first:last))
!     
!                 
!                         
!              end if
!              
!           end do
!        
!        else

    write(out_str, *) 'Calculating homogeneous contribution remainders'
    call out_print(out_str, 1)
    write(out_str, *) ''
    call out_print(out_str, 1)

    call get_2el_mat(0, noc, sum(size_i), Dh, sum(size_i), Fp)
          
    ! Set up null perturbation dimensionality
    pert_xc_null%pdim(1) = sum(size_i)
                
    call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), size(D), D, get_xc_mat, out_print, &
                             sum(size_i), mem_mgr, null_dmat=Dh, fock=Fp)
       
!        end if
! 
!     end if
    
    ind_ctr = 1
    k = 1
    

    if (.NOT.(mem_mgr%calibrate)) then
    
       ! Traverse: Add together Dp and Dh and store, and store perturbed F
       do m = 1, size(cache_outer)
    
          if (cache_outer(m)%dummy_entry) then
       
             cycle
    
          end if
          
          pert = cache_outer(m)%p_tuples(1)
          
          if (pert%do_residues.gt.0) then 
          
             residualization = dabs(dabs(dble(freq_sums(k)))-dabs(dble(pert%exenerg(1)))).lt.xtiny
             
          else
          
             residualization = .false.
             
          end if
         
          do j = 1, size_i(k)
          
             ! DaF: Residue case: Leave alone the Dp part
             if (residualization) then
             
                ierr = QcMatDuplicate_f(Dh(ind_ctr + j - 1),COPY_PATTERN_AND_VALUE,Dp(ind_ctr + j -1))
                
             else
             
                call QcMatRAXPY(1.0d0, Dh(ind_ctr + j - 1), Dp(ind_ctr + j - 1))
                
             end if
             
          end do
       
          len_d = size(F)
       
          call contrib_cache_outer_add_element(len_d, F, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
          len_d = size(D)
       
          call contrib_cache_outer_add_element(len_d, D, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
          
          ind_ctr = ind_ctr + size_i(k)
          k = k + 1
          
       end do
       
       do i = 1, sum(size_i)
       
          call QcMatDst(Dh(i))
          call QcMatDst(Dp(i))
          call QcMatDst(Fp(i))
       
       end do

       deallocate(Dh)
       deallocate(Dp)
       deallocate(Fp)
       
       call QcMatDst(A)
       call QcMatDst(B)
       call QcMatDst(C)
       call QcMatDst(T)
       call QcMatDst(U)
       
       ! MaR: Potential memory issue here, trying a fix
       do m = 1, size(cache_outer)
       
          if (.NOT.(cache_outer(m)%dummy_entry)) then
       
             if (cache_outer(m)%p_tuples(1)%do_residues.gt.0) then
       
                if (allocated(Xx)) then
       
                   call QcMatDst(Xx(1))
                        
                   deallocate(Xx)   
        
                   call mem_decr(mem_mgr, 1)
                
                end if
          
             end if
             
          end if   
          
       end do   
 
    end if
    
    call mem_decr(mem_mgr, 3 * sum(size_i) + 5)
    
    
  end subroutine
  
  
  ! NEW: XC handling
  
  recursive subroutine rsp_xcave_setup_dmat_perts(pert,           &
                                                   sofar,          &
                                                   kn,             &
                                                   rec_prog,       &
                                                   enc_len,        &
                                                   dmat_tuple_len, &
                                                   pert_ids,       &
                                                   enc_perts,      &
                                                   dmat_perts)
  
  
    implicit none
  
    type(p_tuple), intent(in) :: pert
    integer,       intent(in) :: kn(2)
    integer,       intent(in) :: enc_len
    integer,       intent(in) :: dmat_tuple_len
    type(p_tuple), intent(in) :: enc_perts(enc_len)
    type(p_tuple), intent(in) :: dmat_perts(dmat_tuple_len)
    integer, intent(inout) :: sofar
    integer, intent(inout) :: rec_prog
    integer, intent(inout) :: pert_ids(dmat_tuple_len)
    integer       :: i
    integer       :: j
    logical       :: dmat_already
    type(p_tuple) :: psub(pert%npert)

      dmat_already = .false.

      ! unless at final recursion level, recurse further
      ! make all size (n - 1) subsets of the perturbations and recurse
      ! then (at final recursion level) get perturbed F, D, S
      if (pert%npert > 1) then

         call make_p_tuple_subset(pert, psub)

         do i = size(psub), 1, -1

            dmat_already = .false.
            do j = 1, sofar
               if (psub(i)%npert == dmat_perts(j)%npert) then
                  if (pid_compare(psub(i)%npert, psub(i)%pid, dmat_perts(j)%pid)) then
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
         if (pert%npert == enc_perts(j)%npert) then
            if (pid_compare(pert%npert, pert%pid, enc_perts(j)%pid)) then
               dmat_already = .true.
            end if
         end if
      end do

      if (.not. dmat_already) then
         rec_prog = rec_prog + 1
         call p1_cloneto_p2(pert, enc_perts(rec_prog))

         if (.not. kn_skip(pert%npert, pert%pid, kn)) then
            sofar = sofar + 1
            call p1_cloneto_p2(pert, dmat_perts(sofar))
            pert_ids(sofar) = rec_prog + 1
         end if
      end if

   end subroutine
  
  
  
  
  subroutine rsp_xc_wrapper(n_freq_cfgs, pert, kn, len_d, D, get_xc, out_print, &
                                prop_size_total, mem_mgr, null_dmat, fock, prop)
  
    implicit none
    
    type(mem_manager) :: mem_mgr
    integer :: n_freq_cfgs, i, j, k, m, dmat_array_ctr, rec_prog, triang_size
    integer :: dmat_length, dmat_total_size, enc_length, ind_dmat_perts, num_blks
    type(p_tuple) :: emptypert
    type(p_tuple), dimension(n_freq_cfgs), intent(in) :: pert
    type(p_tuple), dimension(:), allocatable :: dmat_perts, enc_perts
    integer, dimension(:), allocatable :: pert_ids, dmat_tuple_sizes, blk_sizes
    integer, dimension(:), allocatable :: pert_freq_category, pert_ext
    integer, dimension(:,:), allocatable ::  blk_info, curr_dmat_indices
    integer,       intent(in) :: kn(2)
    type(QcMat), dimension(:), allocatable :: dmat_total_array
    integer :: len_d
    type(contrib_cache_outer), dimension(len_d), intent(in) :: D
    integer,       intent(in) :: prop_size_total
    type(QcMat), optional, dimension(prop_size_total) :: fock
    type(QcMat), optional, dimension(prop_size_total) :: null_dmat
    complex(8), optional, dimension(prop_size_total) :: prop
    real(8), dimension(prop_size_total) :: prop_dummy
    logical :: srch_fin
    external :: get_xc
    
    external :: out_print
    character(len=2047) :: out_str
    

    ! Special handling for null perturbation case
    
    if (pert(1)%npert == 1) then
    
       if (pert(1)%plab(1) == 'NULL') then
       
          write(out_str, *) 'XC wrapper: Null treatment'
          call out_print(out_str, 3)
          
          

          call p_tuple_to_external_tuple(pert(1), pert(1)%npert, pert_ext)
          
          if (.NOT.(mem_mgr%calibrate)) then
    
             ! Allocate dmat holder array to counter size
             allocate(dmat_total_array(prop_size_total + 1))
       
             call mem_incr(mem_mgr, prop_size_total + 1)

          end if
    
          ! Make the unperturbed D the first elements of the dmat holder array
    
          if (.NOT.(mem_mgr%calibrate)) then

             if (present(null_dmat)) then
          
                call QCMatInit(dmat_total_array(1))

                call contrib_cache_getdata_outer(len_d, D, 1, (/get_emptypert()/), .FALSE., &
                     contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=dmat_total_array(1))
                  
                do i = 1, prop_size_total
             
                   call QCMatInit(dmat_total_array(i + 1))
                   call QcMatAEqB(dmat_total_array(i + 1), null_dmat(i))
             
                end do

             else
             
                write(out_str, *) 'ERROR: Null density matrices required but not present'
                call out_print(out_str, -1)
                
             end if

          end if
          
          
          if (.NOT.(n_freq_cfgs == 1)) then
          
             write(out_str, *) 'ERROR: Null treatment needs only one freq cfg, currently', n_freq_cfgs
             call out_print(out_str, -1)
          
          else
          
             call get_xc(pert(1)%npert, pert_ext, 1, (/1/), &
                 2, (/1, 2/), prop_size_total + 1, dmat_total_array, prop_size_total, fock)
                  
          end if
          
          
          if (.NOT.(mem_mgr%calibrate)) then
    
             do i = 1, prop_size_total + 1
    
                call QcMatDst(dmat_total_array(i))
    
             end do
    
             deallocate(dmat_total_array)
       
          end if
    
          call mem_decr(mem_mgr, prop_size_total + 1)
          
          return
       
       end if
    
    end if
    
    
    

    
    if (present(fock)) then
    
       dmat_length = 2**pert(1)%npert
    
    elseif (present(prop)) then
    
       ! MaR: It seems to hold generally in this case that dmat_length is 2**npert - 1 
       ! even though I can't prove it
    
       dmat_length = 2**(pert(1)%npert - 1)
    
    end if
    
    
    enc_length  = 2**pert(1)%npert
    
    allocate(dmat_perts(dmat_length))
    allocate(enc_perts(enc_length))
    allocate(pert_ids(dmat_length))
    allocate(dmat_tuple_sizes(dmat_length))
    
    emptypert = get_emptypert()
    call p1_cloneto_p2(emptypert, dmat_perts(1))
    

    
    
    pert_ids(1) = 1
    ind_dmat_perts = 1
    rec_prog = 0

    call rsp_xcave_setup_dmat_perts(pert(1),        &
                                    ind_dmat_perts, &
                                    kn,             &
                                    rec_prog,       &
                                    enc_length,     &
                                    dmat_length,    &
                                    pert_ids,       &
                                    enc_perts,      &
                                    dmat_perts)
    
    
    allocate(pert_freq_category(pert(1)%npert * n_freq_cfgs))
    
    ! For each freq config in the original perturbation tuple:
    ! Make "frequency equality array"
    
    do i = 1, n_freq_cfgs
    
       m = 1
    
       do j = 1, pert(1)%npert
       
          srch_fin = .FALSE.
       
          do k = 1, j - 1
          
             if (.NOT.(srch_fin)) then
          
                if (pert(i)%freq(k) == pert(i)%freq(j)) then
             
                   pert_freq_category((i - 1) * pert(1)%npert + j) = &
                   pert_freq_category((i - 1) * pert(1)%npert + k)
                   
                   srch_fin = .TRUE.
                   
                end if
             
             end if
             
          end do
          
          if (.NOT.(srch_fin)) then
          
             pert_freq_category((i - 1) * pert(1)%npert + j) = m
             m = m + 1
          
          end if
       
       
       end do
    
    end do
    
    dmat_total_size = 1
    
    ! For each perturbation subset in dmat_perts:
    
    do i = 2, dmat_length
    
       ! For each freq config:
              
       do j = 1, n_freq_cfgs
       
          ! Dress dmat pert tuple with freqs from present freq config
          do k = 1, dmat_perts(i)%npert
          
             dmat_perts(i)%freq(k) = pert(j)%freq(dmat_perts(i)%pid(k))
          
          end do
          
          num_blks = get_num_blks(dmat_perts(i))
       
          allocate(blk_info(num_blks, 3))
          allocate(blk_sizes(num_blks))
          blk_info = get_blk_info(num_blks, dmat_perts(i))

          ! Get size and add to counter                                           
          dmat_total_size = dmat_total_size + get_triangulated_size(num_blks, blk_info)
          
          deallocate(blk_info)
          deallocate(blk_sizes)
       
       
       
    
       end do
    
    
    
    end do
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       ! Allocate dmat holder array to counter size
       allocate(dmat_total_array(dmat_total_size))
       
       call mem_incr(mem_mgr, dmat_total_size)

    end if
    
    
    ! Make the unperturbed D the first elements of the dmat holder array
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QCMatInit(dmat_total_array(1))

       call contrib_cache_getdata_outer(len_d, D, 1, (/get_emptypert()/), .FALSE., &
            contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=dmat_total_array(1))

    end if
    
        
    dmat_array_ctr = 1
        
    
    ! For each perturbation subset in dmat_perts:
   
    do i = 2, dmat_length
    
       ! For each freq config:
              
       do j = 1, n_freq_cfgs
       
          ! Dress dmat pert tuple with freqs from present freq config
          do k = 1, dmat_perts(i)%npert
          
             dmat_perts(i)%freq(k) = pert(j)%freq(dmat_perts(i)%pid(k))
          
          end do
          
          num_blks = get_num_blks(dmat_perts(i))
       
          allocate(blk_info(num_blks, 3))
          allocate(blk_sizes(num_blks))
          blk_info = get_blk_info(num_blks, dmat_perts(i))
          triang_size = get_triangulated_size(num_blks, blk_info)


          ! Allocate and make indices
          allocate(curr_dmat_indices(triang_size, dmat_perts(i)%npert))
          
          call make_triangulated_indices(num_blks, blk_info, triang_size, curr_dmat_indices)
          
          ! Fill dmat holder array with appropriate matrices
          do k = 1, size(curr_dmat_indices,1)
          
             dmat_array_ctr = dmat_array_ctr + 1
             
             if (.NOT.(mem_mgr%calibrate)) then
             
                call QCMatInit(dmat_total_array(dmat_array_ctr))
           
                call contrib_cache_getdata_outer(len_d, D, 1, (/dmat_perts(i)/), .FALSE., &
                              1, ind_len=size(curr_dmat_indices, 2), &
                              ind_unsorted=curr_dmat_indices(k, :), &
                              mat_sing=dmat_total_array(dmat_array_ctr))
                              
             end if
          
          
          end do
          
          deallocate(curr_dmat_indices)
          
          deallocate(blk_info)
          deallocate(blk_sizes)
    
       end do
    
    end do


    ! Get perturbation tuple in external representation
    call p_tuple_to_external_tuple(pert(1), pert(1)%npert, pert_ext)
    
       write(out_str, *) 'XC wrapper argument summary:'
       call out_print(out_str, 3)
       write(out_str, *) 'pert(1)%npert', pert(1)%npert
       call out_print(out_str, 3)
       write(out_str, *) 'pert_ext', pert_ext
       call out_print(out_str, 3)
       write(out_str, *) 'n_freq_cfgs', n_freq_cfgs
       call out_print(out_str, 3)
       write(out_str, *) 'pert freq category', pert_freq_category
       call out_print(out_str, 3)
       write(out_str, *) 'dmat_length', dmat_length
       call out_print(out_str, 3)
       write(out_str, *) 'pert_ids', pert_ids
       call out_print(out_str, 3)
       write(out_str, *) 'dmat_total_size', dmat_total_size
       call out_print(out_str, 3)
       write(out_str, *) 'prop_size_total', prop_size_total
       call out_print(out_str, 3)

    ! Invoke callback routine
    
    if (present(fock)) then
    
       call get_xc(pert(1)%npert, pert_ext, n_freq_cfgs, pert_freq_category, &
            dmat_length, pert_ids, dmat_total_size, dmat_total_array, prop_size_total, fock)


    elseif (present(prop)) then
    
       call get_xc(pert(1)%npert, pert_ext, n_freq_cfgs, pert_freq_category, &
            dmat_length, pert_ids, dmat_total_size, dmat_total_array, prop_size_total, prop)

    end if
    
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       do i = 1, dmat_total_size

          call QcMatDst(dmat_total_array(i))
    
       end do
    
       deallocate(dmat_total_array)
       
    end if
    
    call mem_decr(mem_mgr, dmat_total_size)
        
    deallocate(pert_freq_category)
    
  
  end subroutine
  
  
  
  
  ! END NEW: XC handling
  
  end module
