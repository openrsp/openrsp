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
  subroutine rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, get_rsp_sol, get_ovl_mat, &
                               get_1el_mat, get_2el_mat, get_xc_mat, dryrun, id_outp, &
                               prog_info, rs_info, sdf_retrieved, mem_mgr,Xf)

    implicit none

    type(mem_manager) :: mem_mgr
    integer :: n_props
    integer, dimension(n_props) :: n_freq_cfgs
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), allocatable, dimension(:) :: p_dummy_orders
    logical :: termination, dryrun, lof_retrieved, sdf_retrieved, rsp_eqn_retrieved, traverse_end
    integer, dimension(3) :: prog_info, rs_info
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    integer :: i, j, k, id_outp, max_order, max_npert, o_size, lof_mem_total
    integer, allocatable, dimension(:) :: size_i
    type(QcMat) :: Fp_dum
    type(contrib_cache_outer) :: F, D, S
    type(contrib_cache_outer), optional :: Xf
    type(contrib_cache), pointer :: cache
    type(contrib_cache), pointer :: cache_next, lof_cache, lof_next
    type(contrib_cache_outer), pointer :: cache_outer_next, lof_outer_next
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat,  get_2el_mat, get_xc_mat

    max_order = max(maxval(kn_rule(:, 1)), maxval(kn_rule(:, 2)))
    max_npert = maxval((/(p_tuples(i)%npert, i = 1, sum(n_freq_cfgs))/))

    ! Make sure that Xf vectors are present in residue calculations
    if (.not.present(Xf).and.p_tuples(1)%do_residues.eq.0) then
      stop 'Residue calculation needs Xf vectors in rsp_fds!'
    end if

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
    call prog_incr(prog_info, 2)
        
    if (rs_check(prog_info, rs_info, lvl=2)) then
          
       write(*,*) ' '
       write(*,*) 'SDF identification was completed in previous'
       write(*,*) 'invocation: Passing to next stage of calculation'
       write(*,*) ' '
       
       allocate(cache)
       call contrib_cache_retrieve(cache, 'OPENRSP_FDS_ID')
       lof_retrieved = .TRUE.
             
    else

       call contrib_cache_allocate(cache)
    
       ! Recurse to identify all necessary perturbed F, D, S
       k = 1 
       do i = 1, n_props
          do j = 1, n_freq_cfgs(i)
       
             call rsp_fds_recurse(p_tuples(k), kn_rule(k, :), max_npert, p_dummy_orders, cache, id_outp)
             k = k + 1
       
          end do
       end do

       call contrib_cache_store(cache, 'OPENRSP_FDS_ID')
       
    end if

    ! NOTE: Something may be wrong about this progress increase, revisit if problems
    call prog_incr(prog_info, 2)
    
    cache_next => cache
    
    lof_retrieved = .FALSE.
    rsp_eqn_retrieved = .FALSE.
    
    ! For each order of perturbation identified (lowest to highest):
    do i = 1, max_order
    
       lof_mem_total = 0
    
       ! Cycle until order reached
       do while(.NOT.(cache_next%p_inner%freq(1) == 1.0*i))
          cache_next => cache_next%next
       end do
      
       ! Contains number of components of perturbed matrices for each perturbation
       allocate(size_i(cache_next%num_outer))       
       k = 1
       
       ! Cycle until at start of outer cache
       cache_outer_next => contrib_cache_outer_cycle_first(cache_next%contribs_outer)
       if (cache_outer_next%dummy_entry) then
          cache_outer_next => cache_outer_next%next
       end if

       ! Traverse to set up size information
       termination = .FALSE.

       do while(.NOT.(termination))
          ! Get number of perturbed matrices for this tuple
          size_i(k) = cache_outer_next%blks_tuple_triang_size(1)
          k = k + 1

          termination = cache_outer_next%last
          cache_outer_next => cache_outer_next%next
          
       end do
       
       ! Cycle until at start of outer cache
       cache_outer_next => contrib_cache_outer_cycle_first(cache_next%contribs_outer)
       if (cache_outer_next%dummy_entry) then
          cache_outer_next => cache_outer_next%next
       end if
             
       ! Check if this stage passed previously and if so, then retrieve and skip execution       
       if (rs_check(prog_info, rs_info, lvl=2)) then
          
          write(*,*) ' '
          write(*,*) 'LOF identification at order', i, 'was completed'
          write(*,*) 'in previous invocation: Passing to next stage of calculation'
          write(*,*) ' '
          
          if (.NOT.(lof_retrieved)) then
          
             allocate(lof_cache)
          
             call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_next => lof_cache
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
       
          if (.NOT.(lof_retrieved)) then
    
             call contrib_cache_allocate(lof_cache)
             
          end if

       
          ! Traverse all elements of outer cache of present cache element
          termination = .FALSE.
          do while(.NOT.(termination))
       
             ! Recurse to identify lower-order Fock matrix contributions
             ! The p_tuples attribute should always be length 1 here, so OK to take the first element
             call rsp_lof_recurse(cache_outer_next%p_tuples(1), cache_outer_next%p_tuples(1)%npert, &
                                         1, (/get_emptypert()/), .TRUE., lof_cache, 1, (/Fp_dum/))
       
             
             termination = cache_outer_next%last
             cache_outer_next => cache_outer_next%next
          
          end do
          
          call contrib_cache_store(lof_cache, 'OPENRSP_LOF_CACHE')
       
       end if
       
       ! Check if this stage passed previously and if so, then retrieve and skip execution
       call prog_incr(prog_info, 2)
       
       if (rs_check(prog_info, rs_info, lvl=2)) then
          
          write(*,*) ' '
          write(*,*) 'LOF calculation at order', i, 'was completed'
          write(*,*) 'in previous invocation: Passing to next stage of calculation'
          write(*,*) ' '
          
          if (.NOT.(lof_retrieved)) then
          
             allocate(lof_cache)
          
             call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_next => lof_cache
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
       
          lof_next => lof_cache
       
          ! Cycle lower-order Fock cache until at start
          do while(.NOT.(lof_next%last))
             lof_next => lof_next%next
          end do
          lof_next => lof_next%next
          if (lof_next%p_inner%npert == 0) then
!             write(*,*) 'cycling dummy'
             lof_next => lof_next%next
          end if
          
       
          ! Traverse lower-order Fock cache and precalculate elements
          termination = .FALSE.
          do while (.NOT.(termination))
       
             ! Check if this stage passed previously and if so, then retrieve and skip execution
             if (rs_check(prog_info, rs_info, lvl=3)) then
          
                write(*,*) ' '
                write(*,*) 'LOF calculation at order', i, ' inside traversal was completed'
                write(*,*) 'in previous invocation: Passing to next stage of calculation'
                write(*,*) ' '
                
                ! Note: No cache retrieval here: In order to get to this position, the
                ! cache would already have been retrieved
          
             else
       
       
                o_size = 0
                
                call rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, &
                                       lof_next, o_size, mem_mgr)
                
                if (mem_exceed(mem_mgr)) then
                
                   return
                
                end if
                
                lof_mem_total = lof_mem_total + lof_next%blks_triang_size * o_size
                
                if (.NOT.(mem_mgr%calibrate)) then
                
                   call contrib_cache_store(lof_next, 'OPENRSP_LOF_CACHE')
                   
                end if
                
             end if
             
             call prog_incr(prog_info, 3)
                       
             termination = (lof_next%last)
             lof_next => lof_next%next
          
          end do
          
          
       
       end if
       
       ! Check if this stage passed previously and if so, then retrieve and skip execution
       call prog_incr(prog_info, 2)
       
       if (rs_check(prog_info, rs_info, lvl=2)) then
          
          write(*,*) ' '
          write(*,*) 'SDF calculation at order', i, 'was completed'
          write(*,*) 'in previous invocation: Passing to next stage of calculation'
          write(*,*) ' '
          
          if (.NOT.(sdf_retrieved)) then
          
             call contrib_cache_outer_retrieve(S, 'OPENRSP_S_CACHE', .FALSE.)
             call contrib_cache_outer_retrieve(D, 'OPENRSP_D_CACHE', .FALSE.)
             call contrib_cache_outer_retrieve(F, 'OPENRSP_F_CACHE', .FALSE.)
             sdf_retrieved = .TRUE.
             
          end if
       
       
       else
       
          ! Calculate all perturbed S, D, F at this order
          call rsp_sdf_calculate(cache_outer_next, cache_next%num_outer, size_i,&
               get_rsp_sol, get_ovl_mat, get_2el_mat, get_xc_mat, F, D, S, lof_next, &
               rsp_eqn_retrieved, prog_info, rs_info, mem_mgr)
               
          if (mem_exceed(mem_mgr)) then
                
                   return
                
          end if
          
          if (.NOT.(mem_mgr%calibrate)) then
          
             call contrib_cache_outer_store(S, 'OPENRSP_S_CACHE')
             call contrib_cache_outer_store(D, 'OPENRSP_D_CACHE')
             call contrib_cache_outer_store(F, 'OPENRSP_F_CACHE')
             
          end if
       
       end if
       
       deallocate(size_i)
       deallocate(lof_cache)
       call mem_decr(mem_mgr, lof_mem_total)
          
    end do
    
    deallocate(p_dummy_orders)
    
  end subroutine

  ! Recursive routine to identify necessary perturbed F, D, S
  recursive subroutine rsp_fds_recurse(pert, kn, max_npert, p_dummy_orders, contribution_cache, id_outp)

    implicit none

    integer :: n_props, max_npert
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%npert) :: psub
    type(p_tuple), dimension(max_npert) :: p_dummy_orders
    integer, dimension(2) :: kn
    integer :: i, j, k, id_outp
    type(contrib_cache) :: contribution_cache
    
    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    if (pert%npert > 1) then

       call make_p_tuple_subset(pert, psub)
 
       do i = 1, size(psub)
       
          if (contrib_cache_already(contribution_cache, 2, (/p_dummy_orders(psub(i)%npert), &
                              p_tuple_standardorder(psub(i))/)) .eqv. .FALSE.) then

             call rsp_fds_recurse(psub(i), kn, max_npert, p_dummy_orders, contribution_cache, id_outp)

          end if

       end do       

    end if

    ! See if already identified, if not and if keeping, then store element
    if (contrib_cache_already(contribution_cache, 2, (/p_dummy_orders(pert%npert), &
                              p_tuple_standardorder(pert)/)) .eqv. .FALSE.) then
         
       if (kn_skip(pert%npert, pert%pid, kn) .eqv. .FALSE.) then

          write(id_outp,*) 'Identified necessary ovlint/fock/density with labels ', pert%plab, &
                     ' and perturbation id ', pert%pid, ' with frequencies (real part)', &
                     real(pert%freq)
          write(id_outp,*) ' '
                 
          k = 1

          do j = 1, pert%npert
             pert%pid(j) = k
             k = k + 1
          end do

          call contrib_cache_add_element(contribution_cache, 2, (/p_dummy_orders(pert%npert), &
                                                                  p_tuple_standardorder(pert)/))

       end if

    end if

  end subroutine
  
  ! Recursive routine to identify (dryrun == .TRUE.) or assemble (dryrun == .FALSE.) 
  ! lower-order perturbed Fock matrix terms
  recursive subroutine rsp_lof_recurse(pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, dryrun, fock_lowerorder_cache, &
                       fp_size, Fp)

    implicit none

    ! fp_size and Fp are dummy if this is a dryrun
    
    logical :: density_order_skip, residue_skip, dryrun
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations
    integer :: fp_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(contrib_cache) :: fock_lowerorder_cache
    type(QcMat), dimension(fp_size) :: Fp

    
    ! If perturbation tuple not empty, recurse further
    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%npert == 0) then

          call rsp_lof_recurse(p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               dryrun, fock_lowerorder_cache, fp_size, Fp)

       else

          call rsp_lof_recurse(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
               p_tuples(2:size(p_tuples))/), dryrun, fock_lowerorder_cache, &
               fp_size, Fp)

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
               t_new, dryrun, fock_lowerorder_cache, fp_size, Fp)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       if (num_p_tuples < 2) then
       
          call rsp_lof_recurse(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples + 1, &
               (/p_tuples(:), p_tuple_getone(pert, 1)/), dryrun, &
               fock_lowerorder_cache, fp_size, Fp)
               
       end if

    ! If recursion is done, proceed to cache manipulation stage
    else

       p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)
       density_order_skip = .FALSE.
       residue_skip = .FALSE.

       if (find_residue_info(p_tuples(1))) residue_skip = .TRUE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%npert >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if ( (density_order_skip .EQV. .FALSE.).AND. &
            (residue_skip .EQV. .FALSE.) ) then
       
          if (contrib_cache_already(fock_lowerorder_cache, &
          num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples))) then

             ! If in cache and not dryrun, then put precalculated data in answer array
             if (.NOT.(dryrun)) then
             
                do i = 1, num_p_tuples
                   if (i == 1) then
                      write(*,*) 'F ', p_tuples(i)%plab
                   else
                      write(*,*) 'D ', p_tuples(i)%plab
                   end if
                end do
             
                write(*,*) 'Getting lower-order perturbed Fock matrix contribution from cache'
                call contrib_cache_getdata(fock_lowerorder_cache, num_p_tuples, p_tuples, &
                contrib_size=fp_size, ind_len=0, mat=Fp)
             
!                 write(*,*) 'Got contribution'
             
             end if   
          
          else

             ! If dryrun and not in cache, then add identified element to cache
             ! for later calculation upon traversal
             if (dryrun) then
          
                write(*,*) 'Identified lower-order perturbed Fock matrix contribution'

                do i = 1, num_p_tuples
                   if (i == 1) then
                      write(*,*) 'F', p_tuples(i)%pid
                   else
                      write(*,*) 'D', p_tuples(i)%pid
                   end if
                end do
                
                call contrib_cache_add_element(fock_lowerorder_cache, num_p_tuples, &
                     p_tuples_standardorder(num_p_tuples, p_tuples))
                
             else
             
                write(*,*) 'ERROR: Wanted lower-order perturbed Fock matrix contribution but it was not in cache'
                
             end if
          
          end if

       end if

    end if

  end subroutine
  
  ! Calculate lower-order Fock contributions for a given inner perturbation tuple
  subroutine rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, cache, &
                               total_outer_size_1, mem_mgr)

    implicit none

    type(mem_manager) :: mem_mgr
    integer :: mctr, mcurr, miter, msize, octr, mem_track
    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n, s, p
    integer :: id_outp
    integer :: total_outer_size_1, c1_ctr, lhs_ctr_1, num_pert
    integer :: num_0, num_1, tot_num_pert, offset
    character(30) :: mat_str, fmt_str
    type(contrib_cache) :: cache
    type(contrib_cache_outer) :: D
    type(contrib_cache_outer), pointer :: outer_next
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, contrib_0, contrib_1
    type(QcMat) :: D_unp
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    integer, allocatable, dimension(:,:) :: blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
    external :: get_1el_mat, get_ovl_mat, get_2el_mat
    
    write(*,*) 'Calculating lower-order Fock matrix contribution for inner tuple', cache%p_inner%plab
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)
    outer_next => cache%contribs_outer

    allocate(outer_contract_sizes_1(cache%num_outer))    
   
    call mem_incr(mem_mgr, 1)
   
    if (.NOT.(mem_mgr%calibrate)) then
    
       call QCMatInit(D_unp)

       call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
            contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)

    end if
   
   
    ! Traversal: Find number of density matrices for contraction for nuc-nuc, 1-el, 2-el cases
    traverse_end = .FALSE.
    
    outer_next => contrib_cache_outer_cycle_first(outer_next)
    if (outer_next%dummy_entry) then
       outer_next => outer_next%next
    end if
       
    total_outer_size_1 = 0
    num_1 = 0
        
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
       ! No chain rule application: Both 1-el and 2-el contributions
       if (outer_next%num_dmat == 0) then

          outer_contract_sizes_1(k) = 1
          total_outer_size_1 = total_outer_size_1 + 1
       
       ! One chain rule application: Only 2-el contributions
       else if (outer_next%num_dmat == 1) then
       
          num_1 = num_1 + 1
          outer_contract_sizes_1(k) = outer_next%blks_tuple_triang_size(1)
          total_outer_size_1 = total_outer_size_1 + outer_next%blks_tuple_triang_size(1)

       end if
       
       ! Check if enough memory
       if (.NOT.(mem_enough(mem_mgr, cache%blks_triang_size*outer_contract_sizes_1(k)))) then
       
          ! Not possible to run; flag it and return
          call mem_set_status(mem_mgr, 2)
          return
          
       end if
       
       call mem_incr(mem_mgr, cache%blks_triang_size*outer_contract_sizes_1(k))
       allocate(outer_next%data_mat(cache%blks_triang_size*outer_contract_sizes_1(k)))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          do i = 1, cache%blks_triang_size*outer_contract_sizes_1(k)
    
             call QcMatInit(outer_next%data_mat(i), D_unp)
             call QcMatZero(outer_next%data_mat(i))
                
          end do
          
       end if
       
   
       if (outer_next%num_dmat == 0) then
       
          write(*,*) 'All inner contribution'
 
       else
       
          write(*,*) 'Outer contribution:'
          
       end if
       
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'D ', outer_next%p_tuples(i)%plab
       
       end do
    
       if (outer_next%last) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
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
       
       
       ! Traversal: Get matrices for contraction from cache
       
       traverse_end = .FALSE.
       
       outer_next => contrib_cache_outer_cycle_first(outer_next)
       if (outer_next%dummy_entry) then
          outer_next => outer_next%next
       end if
          
       k = 1
       lhs_ctr_1 = 1
       
       do while (traverse_end .EQV. .FALSE.)
       
          ! One chain rule application
          if (outer_next%num_dmat == 1) then
          
             do m = 1, outer_contract_sizes_1(k) 
             
                if ((lhs_ctr_1 + m - 1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                   if (.NOT.(mem_mgr%calibrate)) then
                
                      call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                           contrib_size=1, ind_len=1, ind_unsorted=outer_next%indices(m, :), &
                           mat_sing=LHS_dmat_1(mctr + 1))
                           
                   end if
                        
                        mctr = mctr + 1
                
                end if
       
             end do
             
          ! No chain rule applications: Contraction matrix is only unperturbed D
          elseif(outer_next%num_dmat == 0) then
             
             if ((lhs_ctr_1 == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                num_0 = 1
                
                if (.NOT.(mem_mgr%calibrate)) then
                
                   call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
                        contrib_size=1, ind_len=1, ind_unsorted=(/1/), &
                        mat_sing=LHS_dmat_1(mctr + 1))
                        ! Was: mat_sing=LHS_dmat_1(1)), change if problems
                        
                end if
                     
                mctr = mctr + 1
                
             end if
             
          end if
          
          if (outer_next%last) then
             traverse_end = .TRUE.
          end if
          
          lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
          k = k + 1
       
          outer_next => outer_next%next
       
       end do
       
       ! Calculate contributions
       
       ! Calculate one-electron contributions
       if (num_0 > 0) then
       
          if (.NOT.(mem_mgr%calibrate)) then
       
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
       
       traverse_end = .FALSE.
       
       outer_next => contrib_cache_outer_cycle_first(outer_next)
       if (outer_next%dummy_entry) then
          outer_next => outer_next%next
       end if
         
       k = 1
       c1_ctr = 1
       mctr = 0
       octr = 1
       
       do while (traverse_end .EQV. .FALSE.)
     
          ! One-el and two-el contributions ("all inner contribution")
          if (outer_next%num_dmat == 0) then
          
             ! For "all inner contribution", the data is already ordered correctly
             if ((octr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                if (.NOT.(mem_mgr%calibrate)) then
             
                   do i = 1, cache%blks_triang_size
             
                      call QcMatkAB(1.0d0, contrib_0(i), contrib_1(c1_ctr + i - 1), outer_next%data_mat(i))
          
                   end do
                   
                end if
                
                c1_ctr = c1_ctr + cache%blks_triang_size
                mctr = mctr + 1
                
             end if
             
             octr = octr + 1
     
          ! Only two-el contribution
          else if (outer_next%num_dmat == 1) then
                     
             ! Initialize block information for cache indexing
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
             
                if ((octr == mcurr + mctr) .AND. .NOT.(msize <= mctr)) then
             
                   do j = 1, size(cache%indices, 1)
                   
                      offset = get_triang_blks_tuple_offset(outer_next%num_dmat + 1, &
                      cache%p_inner%npert + sum((/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/)), &
                      (/cache%nblks, (/(outer_next%nblks_tuple(m), m = 1, outer_next%num_dmat) /) /), &
                      (/cache%p_inner%npert, (/(outer_next%p_tuples(m)%npert, m = 1, outer_next%num_dmat)/)/), &
                      blks_tuple_info, &
                      blk_sizes, &
                      (/cache%blks_triang_size, &
                      (/(outer_next%blks_tuple_triang_size(m), m = 1, outer_next%num_dmat)/)/), &
                      (/cache%indices(j, :), outer_next%indices(i, :)/))

                      if (.NOT.(mem_mgr%calibrate)) then
                      
                         ! Store result in cache
                         call QcMatRAXPY(1.0d0, contrib_1(c1_ctr + j - 1), &
                         outer_next%data_mat(offset))
                         
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
          
          if (outer_next%last) then
             traverse_end = .TRUE.
          end if
       
          k = k + 1
          outer_next => outer_next%next
       
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
  subroutine rsp_sdf_calculate(cache_outer, num_outer, size_i, &
  get_rsp_sol, get_ovl_mat, get_2el_mat, get_xc_mat, F, D, S, lof_cache, &
  rsp_eqn_retrieved, prog_info, rs_info, mem_mgr)
  
    implicit none
    
    type(mem_manager) :: mem_mgr
    integer :: mctr, mcurr, miter, msize, octr, mem_track
    logical :: termination, rsp_eqn_retrieved
    integer :: num_outer, ind_ctr, npert_ext, sstr_incr, superstructure_size
    integer :: i, j, k, m, w, nblks
    integer :: first, last
    integer, dimension(0) :: noc
    integer, dimension(3) :: prog_info, rs_info
    integer, allocatable, dimension(:) :: pert_ext, blk_sizes, ind
    integer, allocatable, dimension(:,:) :: blk_info
    integer, allocatable, dimension(:,:) :: indices
    integer, dimension(num_outer) :: size_i
    character(30) :: mat_str, fmt_str
    complex(8), dimension(num_outer) :: freq_sums
    type(p_tuple) :: pert, pert_xc_null
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(contrib_cache) :: lof_cache
    type(contrib_cache_outer), target :: cache_outer
    type(contrib_cache_outer) :: F, D, S
    type(contrib_cache_outer), pointer :: cache_outer_next
    type(Qcmat), allocatable, dimension(:) :: Dh, Dp, Fp, Sp, RHS, X
    type(Qcmat) :: A, B, C, T, U
    external :: get_rsp_sol, get_ovl_mat,  get_2el_mat, get_xc_mat
    

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
    
       call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=A)
       call contrib_cache_getdata_outer(S, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=B)
       call contrib_cache_getdata_outer(F, 1, (/get_emptypert()/), .FALSE., contrib_size=1, & 
       ind_len=1, ind_unsorted=(/1/), mat_sing=C)
    
    end if
    
    cache_outer_next => cache_outer
    
    
    
    
    
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
             cache_outer_next => cache_outer_next%next
    end if
 
 
 
    ! Traverse cache elements and calculate pertubed S
    ind_ctr = 1
    k = 1
    termination = .FALSE.
    
    do while(.NOT.(termination))

       pert = cache_outer_next%p_tuples(1)

       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          allocate(Sp(size_i(k)))
          
          do i = 1, size_i(k)
          
             call QcMatInit(Sp(i), A)
             call QcMatZero(Sp(i))
     
          end do
       
     
          ! For each cache element:
          ! Calculate Sp
          call p_tuple_to_external_tuple(pert, npert_ext, pert_ext)
          call get_ovl_mat(0, noc, 0, noc, npert_ext, pert_ext, &
                        size_i(k), Sp)
          
          
          deallocate(pert_ext)
          
       end if
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
          
          
          call contrib_cache_outer_add_element(S, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Sp )
          
          do i = 1, size_i(k)
          
             call QcMatDst(Sp(i))
     
          end do
          
          deallocate(Sp)
          
       end if
       
       call mem_decr(mem_mgr, size_i(k))
       
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
       
       ! Set up next iteration
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
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
    
    
    
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
             cache_outer_next => cache_outer_next%next
    end if
 
 
 
    ! Traverse cache elements and do various stages before collective 2-el
    ! particular Fock matrix contribution call
    ind_ctr = 1
    k = 1
    termination = .FALSE.
    
    do while(.NOT.(termination))

       pert = cache_outer_next%p_tuples(1)
       
       ! Get frequency sum
       freq_sums(k) = sum(real(pert%freq(:))) 
       
       ! Set up block info
       nblks = get_num_blks(pert)

       allocate(blk_info(nblks, 3))
       allocate(blk_sizes(pert%npert))
       blk_info = get_blk_info(nblks, pert)
       blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))
    
       ! Add the initialized Dp to cache
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          call contrib_cache_outer_add_element(D, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
      
          ! Assemble Fp (lower-order) for all components and add to cache
            
          call rsp_lof_recurse(pert, pert%npert, &
                               1, (/get_emptypert()/), .FALSE., lof_cache, size_i(k), &
                               Fp=Fp(ind_ctr:ind_ctr + size_i(k) - 1))
                               
          ! XC call should go here
          
          ! Currently only one freq. configuration
          ! The (k,n) rule argument is adapted for the Fock contribution case
          call rsp_xc_wrapper(1, (/pert/), (/pert%npert, pert%npert/), D, get_xc_mat, &
                                size_i(k), mem_mgr, fock=Fp(ind_ctr:ind_ctr + size_i(k) - 1))
          
       
       end if
       
       call mem_incr(mem_mgr, size_i(k))
       
       if (.NOT.(mem_mgr%calibrate)) then
       
          call contrib_cache_outer_add_element(F, .FALSE., 1, & 
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
            (/get_emptypert(), get_emptypert(), get_emptypert()/), &
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
                  ind, F, D, S, Dp(ind_ctr + j - 1))
                  
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
       
          call contrib_cache_outer_add_element(D, .FALSE., 1, & 
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
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
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
    
    

    
!   Kept for possible reintroduction: Calculate for subset of all contraction
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
    
    if (.NOT.(mem_mgr%calibrate)) then
    
       ! Outside traversal:
       ! Complete Fp using Dp
       call get_2el_mat(0, noc, sum(size_i), Dp, sum(size_i), Fp)
       
     
       ! Set up null perturbation dimensionality
       pert_xc_null%pdim(1) = sum(size_i)
    
       call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), D, get_xc_mat, &
                                sum(size_i), mem_mgr, null_dmat=Dp, fock=Fp)
    
    
    end if
    
!     end if

    
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
       cache_outer_next => cache_outer_next%next
    end if
       
    ! Traversal: Adding Fp to cache, constructing RHS
    termination = .FALSE.
    do while(.NOT.(termination))

       pert = cache_outer_next%p_tuples(1)
    
       ! Set up block info
       nblks = get_num_blks(pert)

       allocate(blk_info(nblks, 3))
       allocate(blk_sizes(pert%npert))
       blk_info = get_blk_info(nblks, pert)
       blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

       if (.NOT.(mem_mgr%calibrate)) then
       
         ! Add the completed Fp to cache
           call contrib_cache_outer_add_element(F, .FALSE., 1, & 
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
            (/get_emptypert(), get_emptypert(), get_emptypert()/), &
            superstructure_size, sstr_incr, derivative_structure)
      
       call make_triangulated_indices(nblks, blk_info, size_i(k), indices)
       
       ! Calculate RHS matrices
       do j = 1, size(indices, 1)
       
          ind = indices(j, :)
          
          call mem_incr(mem_mgr, 4)
          if (.NOT.(mem_mgr%calibrate)) then
          
             call rsp_get_matrix_y(superstructure_size, derivative_structure, &
                   pert%npert, (/ (m, m = 1, pert%npert) /), &
                  pert%npert, ind, F, D, S, RHS(ind_ctr + j - 1))
                  
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
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do

    ! Outside traversal:
    ! Solve all response equations (opportunities for optimization for e.g.
    ! first-order EL, but that can be introduced later)
    
    ! IN NEW FORMAT, AWAITING CHANGES TO CALLBACK STRUCTURE
!     call get_rsp_sol(num_outer, size_i, (/(i - i + 1, i = 1, num_outer)/), freq_sums, RHS, X)

    m  = 1
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
       cache_outer_next => cache_outer_next%next
    end if

    
    
    ! Traverse to solve response equations
    ! May extend to do collectively all cache elements for which the
    ! frequency sum is the same
    termination = .FALSE.
    do while(.NOT.(termination))

       write(*,*) 'Frequency sum:', freq_sums(k)
    
       if (size_i(k) > m) then
    
          do i = 1, size_i(k)/m + 1
    
             if (.NOT.((i - 1) *m >= size_i(k))) then
    
                first = (i - 1) * m + 1
                last = min(i * m, size_i(k))
                
                if (rs_check(prog_info, rs_info, lvl=3)) then
                
                   write(*,*) ' '
                   write(*,*) 'RSP eqn solution batch was completed'
                   write(*,*) 'in previous invocation: Passing to next stage of calculation'
                   write(*,*) ' '
          
                   if (.NOT.(rsp_eqn_retrieved)) then
                   
                      write(*,*) 'Retrieving RSP eqn solutions from disk'
                      write(*,*) ' '
          
                      call mat_scal_retrieve(rs_info(3), 'OPENRSP_MAT_RSP', mat=X(1:rs_info(3)))
                      rsp_eqn_retrieved = .TRUE.
                      
                      write(*,*) 'Finished retrieval'
             
                   end if

                else

                   if (.NOT.(mem_mgr%calibrate)) then
                   
                      call get_rsp_sol(1,                                    &
                                       (/last-first+1/),                     &
                                       (/1/),                                &
                                       dcmplx(real((/freq_sums(k)/)),0.0d0), &
                                       RHS(ind_ctr+first-1:ind_ctr+last-1),  &
                                       X(ind_ctr+first-1:ind_ctr+last-1))
                   
                      call mat_scal_store(last - first + 1, 'OPENRSP_MAT_RSP', &
                           mat=X(ind_ctr+first-1:ind_ctr+last-1), start_pos = ind_ctr+first-1)
                    
                    end if
                   
                end if
                
                call prog_incr(prog_info, 3)
                
   
             end if
       
          end do
       
       else
       
          ! Check if this stage passed previously and if so, then retrieve and skip execution
          if (rs_check(prog_info, rs_info, lvl=3)) then
                
             write(*,*) ' '
             write(*,*) 'RSP eqn solution batch was completed'
             write(*,*) 'in previous invocation: Passing to next stage of calculation'
             write(*,*) ' '
         
             if (.NOT.(rsp_eqn_retrieved)) then
          
                call mat_scal_retrieve(rs_info(3), 'OPENRSP_MAT_RSP', mat=X(1:rs_info(3)))
                rsp_eqn_retrieved = .TRUE.
             
             end if

          else
       
             if (.NOT.(mem_mgr%calibrate)) then
       
                call get_rsp_sol(1,                                    &
                                 (/size_i(k)/),                        &
                                 (/1/),                                &
                                 dcmplx(real((/freq_sums(k)/)),0.0d0), &
                                 RHS(ind_ctr:ind_ctr+size_i(k)-1),     &
                                 X(ind_ctr:ind_ctr+size_i(k)-1))
          
                call mat_scal_store(size_i(k), 'OPENRSP_MAT_RSP', &
                           mat=X(ind_ctr:ind_ctr+size_i(k)-1), start_pos = ind_ctr)
                           
             end if
                   
          end if
          
          call prog_incr(prog_info, 3)
          
          
    
       end if
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
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
            
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
             cache_outer_next => cache_outer_next%next
    end if
       
    ! Traverse: Make homogeneous contribution to perturbed D
    termination = .FALSE.
    do while(.NOT.(termination))

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
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
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
!                 call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), D, get_xc_mat, &
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
          
          call get_2el_mat(0, noc, sum(size_i), Dh, sum(size_i), Fp)
          
          ! Set up null perturbation dimensionality
          pert_xc_null%pdim(1) = sum(size_i)
                
          call rsp_xc_wrapper(1, (/pert_xc_null/), (/1, 1/), D, get_xc_mat, &
                                    sum(size_i), mem_mgr, null_dmat=Dh, fock=Fp)
       
!        end if
! 
!     end if
    
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next => contrib_cache_outer_cycle_first(cache_outer_next)
    if (cache_outer_next%dummy_entry) then
       cache_outer_next => cache_outer_next%next
    end if
    
    if (.NOT.(mem_mgr%calibrate)) then
       
       ! Traverse: Add together Dp and Dh and store, and store perturbed F
       termination = .FALSE.
       do while(.NOT.(termination))

          pert = cache_outer_next%p_tuples(1)
         
          do j = 1, size_i(k)
             call QcMatRAXPY(1.0d0, Dh(ind_ctr + j - 1), Dp(ind_ctr + j - 1))
          end do
       
          call contrib_cache_outer_add_element(F, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
          call contrib_cache_outer_add_element(D, .FALSE., 1, & 
               (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
          
          ind_ctr = ind_ctr + size_i(k)
          k = k + 1
       
          termination = cache_outer_next%last
          cache_outer_next => cache_outer_next%next
          
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
  
  
  
  
  subroutine rsp_xc_wrapper(n_freq_cfgs, pert, kn, D, get_xc, &
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
    type(contrib_cache_outer), intent(in) :: D
    integer,       intent(in) :: prop_size_total
    type(QcMat), optional, dimension(prop_size_total) :: fock
    type(QcMat), optional, dimension(prop_size_total) :: null_dmat
    complex(8), optional, dimension(prop_size_total) :: prop
    logical :: srch_fin
    external :: get_xc
    

    ! Special handling for null perturbation case
    
    if (pert(1)%npert == 1) then
    
       if (pert(1)%plab(1) == 'NULL') then
       
          write(*,*) 'Special case: null treatment'
          
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

                call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
                     contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=dmat_total_array(1))
                  
                do i = 1, prop_size_total
             
                   call QCMatInit(dmat_total_array(i + 1))
                   call QcMatAEqB(dmat_total_array(i + 1), null_dmat(i))
             
                end do

             else
             
                write(*,*) 'ERROR: Null density matrices required but not present'
                
             end if

          end if
          
          
          if (.NOT.(n_freq_cfgs == 1)) then
          
             write(*,*) 'ERROR: Null treatment needs only one freq cfg, currently', n_freq_cfgs
          
          else
          
!              call get_xc(pert(1)%npert, pert_ext, 1, (/1/), &
!                   2, (/1, 2/), prop_size_total + 1, dmat_total_array, prop_size_total, fock)
                  
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
    
    
    write(*,*) 'xc dmatlen', dmat_length
    
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
    
    do i = 1, dmat_length
    
       ! For each freq config:
              
       do j = 1, n_freq_cfgs
       
          ! Dress dmat pert tuple with freqs from present freq config
          do k = 1, dmat_perts(i)%npert
          
             write(*,*) 'dmat freq', dmat_perts(i)%freq
          
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

       call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
            contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=dmat_total_array(1))

    end if
    
        
    dmat_array_ctr = 1
        
    
    ! For each perturbation subset in dmat_perts:
   
    do i = 1, dmat_length
    
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
           
                call contrib_cache_getdata_outer(D, 1, (/dmat_perts(i)/), .FALSE., &
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
    
    write(*,*) 'XC wrapper argument summary:'
    
    write(*,*) 'pert(1)%npert', pert(1)%npert
    write(*,*) 'pert_ext', pert_ext
    write(*,*) 'n_freq_cfgs', n_freq_cfgs
    write(*,*) 'pert freq category', pert_freq_category
    
    write(*,*) 'dmat_length', dmat_length
    write(*,*) 'pert_ids', pert_ids
    write(*,*) 'dmat_total_size', dmat_total_size
    write(*,*) 'prop_size_total', prop_size_total
    
    ! Invoke callback routine
    
    if (present(fock)) then
    
!        call get_xc(pert(1)%npert, pert_ext, n_freq_cfgs, pert_freq_category, &
!             dmat_length, pert_ids, dmat_total_size, dmat_total_array, prop_size_total, fock)

    elseif (present(prop)) then
    
!     call get_xc(pert(1)%npert, pert_ext, n_freq_cfgs, pert_freq_category, &
!          dmat_length, pert_ids, dmat_total_size, dmat_total_array, prop_size_total, prop)
    
    
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
