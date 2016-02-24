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

  private

  contains
  
  recursive subroutine rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, get_rsp_sol, get_ovl_mat, &
                               get_1el_mat, get_2el_mat, get_xc_mat, dryrun, cache, id_outp, &
                               prog_info, rs_info, sdf_retrieved)

    implicit none

    integer :: n_props
    integer, dimension(n_props) :: n_freq_cfgs
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), allocatable, dimension(:) :: p_dummy_orders
    logical :: termination, dryrun, lof_retrieved, sdf_retrieved, rsp_eqn_retrieved
    integer, dimension(3) :: prog_info, rs_info
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    integer :: i, j, k, id_outp, max_order, max_npert
    integer, allocatable, dimension(:) :: size_i
    type(QcMat) :: Fp_dum
    type(contrib_cache_outer) :: F, D, S
    type(contrib_cache), target :: cache
    type(contrib_cache), pointer :: cache_next, lof_cache, lof_next
    type(contrib_cache_outer), pointer :: cache_outer_next
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat,  get_2el_mat, get_xc_mat

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

    call prog_incr(prog_info, 2)
    
    if (rs_check(prog_info, rs_info, lvl=2)) then
          
       write(*,*) ' '
       write(*,*) 'SDF identification was completed in previous'
       write(*,*) 'invocation: Passing to next stage of calculation'
       write(*,*) ' '
          
       call contrib_cache_retrieve(cache, 'OPENRSP_FDS_ID')
       lof_retrieved = .TRUE.
             
    else
    
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
           
    cache_next => cache
    
    lof_retrieved = .FALSE.
    rsp_eqn_retrieved = .FALSE.
    ! CONTINUE HERE: INTRODUCE RESTARTING SCHEME
    
    !   For each order of perturbation identified (lowest to highest):
    do i = 1, max_order
    
       call prog_incr(prog_info, 2)
       
       if (rs_check(prog_info, rs_info, lvl=2)) then
          
          write(*,*) ' '
          write(*,*) 'LOF identification at order', i, 'was completed'
          write(*,*) 'in previous invocation: Passing to next stage of calculation'
          write(*,*) ' '
          
          if (.NOT.(lof_retrieved)) then
          
             call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
    
          call contrib_cache_allocate(lof_cache)

          ! Cycle until order reached
          do while(.NOT.(cache_next%p_inner%freq(1) == 1.0*i))
             cache_next => cache_next%next
          end do
      
          ! Contains number of components of perturbed matrices for each perturbation
          allocate(size_i(cache_next%num_outer))       
          k = 1
       
          ! Cycle until at start of outer cache
          cache_outer_next => contrib_cache_outer_cycle_first(cache_next%contribs_outer)
          cache_outer_next => cache_outer_next%next
       
          ! Traverse all elements of outer cache of present cache element
          termination = .FALSE.
          do while(.NOT.(termination))
       
             ! Recurse to identify lower-order Fock matrix contributions
             ! The p_tuples attribute should always be length 1 here, so OK to take the first element
             call rsp_lof_recurse(cache_outer_next%p_tuples(1), cache_outer_next%p_tuples(1)%npert, &
                                         1, (/get_emptypert()/), .TRUE., lof_cache, 1, (/Fp_dum/))
       
             ! Get number of perturbed matrices for this tuple
             size_i(k) = cache_outer_next%blks_tuple_triang_size(1)
             k = k + 1

             termination = cache_outer_next%last
             cache_outer_next => cache_outer_next%next
          
          end do
          
          call contrib_cache_store(lof_cache, 'OPENRSP_LOF_CACHE')
       
       end if
       
       call prog_incr(prog_info, 2)
       
       if (rs_check(prog_info, rs_info, lvl=2)) then
          
          write(*,*) ' '
          write(*,*) 'LOF calculation at order', i, 'was completed'
          write(*,*) 'in previous invocation: Passing to next stage of calculation'
          write(*,*) ' '
          
          if (.NOT.(lof_retrieved)) then
          
             call contrib_cache_retrieve(lof_cache, 'OPENRSP_LOF_CACHE')
             lof_retrieved = .TRUE.
             
          end if
       
       
       else
       
          lof_next => lof_cache
       
          ! Cycle lower-order Fock cache until at start
          do while(.NOT.(lof_next%last))
             lof_next => lof_next%next
          end do
          lof_next => lof_next%next
          lof_next => lof_next%next
       
          ! Traverse lower-order Fock cache and precalculate elements
          termination = .FALSE.
          do while (.NOT.(termination))
       
             call rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, lof_next)
          
             termination = (lof_next%last)
             lof_next => lof_next%next
          
          end do
          
          call contrib_cache_store(lof_next, 'OPENRSP_LOF_CACHE')
       
       end if
       
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
               get_rsp_sol, get_ovl_mat, get_2el_mat, F, D, S, lof_next, &
               rsp_eqn_retrieved, prog_info, rs_info)
            
            
          call contrib_cache_outer_store(S, 'OPENRSP_S_CACHE')
          call contrib_cache_outer_store(D, 'OPENRSP_D_CACHE')
          call contrib_cache_outer_store(F, 'OPENRSP_F_CACHE')
       
       end if
       
       deallocate(size_i)
       deallocate(lof_cache)
          
    end do
    
    deallocate(p_dummy_orders)
    
  end subroutine


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
  
  recursive subroutine rsp_lof_recurse(pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, dryrun, fock_lowerorder_cache, &
                       fp_size, Fp)

    implicit none

    ! fp_size and Fp are dummy if this is a dryrun
    
    logical :: density_order_skip, dryrun
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

!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)
       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%npert >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then
       
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
                
                call contrib_cache_add_element(fock_lowerorder_cache, num_p_tuples, p_tuples)
                
             else
             
                write(*,*) 'ERROR: Wanted lower-order perturbed Fock matrix contribution but it was not in cache'
                
             end if
          
          end if

       end if

    end if

  end subroutine
  
  subroutine rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, cache)

    implicit none

    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n, s
    integer :: id_outp
    integer :: total_outer_size_1, c1_ctr, lhs_ctr_1, num_pert
    integer :: num_0, num_1
    character(30) :: mat_str, fmt_str
    type(contrib_cache) :: cache
    type(contrib_cache_outer) :: D
    type(contrib_cache_outer), pointer :: outer_next
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, contrib_0, contrib_1
    type(QcMat) :: D_unp
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    external :: get_1el_mat, get_ovl_mat, get_2el_mat
    
    write(*,*) 'Calculating lower-order Fock matrix contribution for perturbation', cache%p_inner%plab
    
    call p_tuple_to_external_tuple(cache%p_inner, num_pert, pert_ext)
    outer_next => cache%contribs_outer

    allocate(outer_contract_sizes_1(cache%num_outer))    
   
    ! Traversal: Find number of density matrices for contraction for nuc-nuc, 1-el, 2-el cases
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
       
    total_outer_size_1 = 0
    num_0 = 0
    num_1 = 0
        
    k = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
       if (outer_next%num_dmat == 0) then

          num_0 = 1
          outer_contract_sizes_1(k) = 1
          total_outer_size_1 = total_outer_size_1 + 1
          
       else if (outer_next%num_dmat == 1) then
       
          num_1 = num_1 + 1
          outer_contract_sizes_1(k) = outer_next%blks_tuple_triang_size(1)
          total_outer_size_1 = total_outer_size_1 + outer_next%blks_tuple_triang_size(1)

       end if
   
       if (outer_next%num_dmat == 0) then
       
          write(*,*) 'All inner contribution', total_outer_size_1
 
       else
       
          write(*,*) 'Outer contribution:', total_outer_size_1
          
       end if
       
       do i = 1, outer_next%num_dmat
          
          write(*,*) 'D ', outer_next%p_tuples(i)%plab
       
       end do
    
       if (outer_next%next%dummy_entry) then
    
          traverse_end = .TRUE.
    
       end if
    
       k = k + 1
    
       outer_next => outer_next%next
    
    end do

    ! Initializing data and arrays for external calls
    
    allocate(LHS_dmat_1(sum(outer_contract_sizes_1)))
    
    call QCMatInit(D_unp)

    call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
         contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)

    do i = 1, size(LHS_dmat_1)
    
      call QCMatInit(LHS_dmat_1(i), D_unp)
    
    end do
    
    allocate(contrib_0(cache%blks_triang_size))
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    
    
    do i = 1, size(contrib_0)
      call QCMatInit(contrib_0(i), D_unp)
      call QCMatZero(contrib_0(i))
    end do
    
    do i = 1, size(contrib_1)
      call QCMatInit(contrib_1(i), D_unp)
      call QCMatZero(contrib_1(i))
    end do

    ! Traversal: Get matrices for contraction from cache

    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
       
    k = 1
    lhs_ctr_1 = 1
    
!     write(*,*) 'outer ctr sizes 1', outer_contract_sizes_1
    
    do while (traverse_end .EQV. .FALSE.)

    
!        write(*,*) 'assembling contraction D'
       ! One chain rule application
       if (outer_next%num_dmat == 1) then
       
          do m = 1, outer_contract_sizes_1(k) 
          
!           write(*,*) 'Data ctr', lhs_ctr_1 + m  - 1
!           stop

             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                  contrib_size=1, ind_len=1, ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_1(lhs_ctr_1 + m  - 1))

          end do
          
          
       elseif(outer_next%num_dmat == 0) then
           
!            write(*,*) 'all inner'
       
           call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
                  contrib_size=1, ind_len=1, ind_unsorted=(/1/), &
                  mat_sing=LHS_dmat_1(1))
          
       
       end if
   
       if (outer_next%next%dummy_entry) then
          traverse_end = .TRUE.
       end if

       lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    
    ! Calculate contributions
    
    ! Calculate one-electron contributions
    if (num_0 > 0) then
    
       call get_1el_mat(num_pert, pert_ext, size(contrib_0), contrib_0)
      
!       s = QcMatWrite_f(contrib_0(1), 'contrib_0', ASCII_VIEW)
      
       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()
      
!        call rsp_ovlave_t_matrix(get_ovl_mat, cache%p_inner, cache%p_inner%npert, &
!                                 t_matrix_bra, t_matrix_ket, 1, &
!                                 D_unp, size(contrib_0), contrib_0)
    
    end if
    
    ! Calculate two-electron contributions
    call get_2el_mat(num_pert, pert_ext, size(LHS_dmat_1), LHS_dmat_1, &
                     size(contrib_1), contrib_1)
                       

                       
!          s = QcMatWrite_f(LHS_dmat_1(1), 'lof_lhs_dmat', ASCII_VIEW)   
!         s = QcMatWrite_f(contrib_1(1), 'contrib_1', ASCII_VIEW)                       
    ! Traversal: Add 1-el and two-el contributions together
    
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
      
    k = 1
    c1_ctr = 1
    
    do while (traverse_end .EQV. .FALSE.)
  
       ! One-el and two-el contributions
       if (outer_next%num_dmat == 0) then
       
          allocate(outer_next%data_mat(cache%blks_triang_size*outer_contract_sizes_1(k)))
          
          do i = 1, cache%blks_triang_size*outer_contract_sizes_1(k)

!              write(*,*) 'i', i
             call QcMatInit(outer_next%data_mat(i), D_unp)
             call QcMatZero(outer_next%data_mat(i))
             call QcMatkAB(1.0d0, contrib_0(i), contrib_1(c1_ctr + i - 1), outer_next%data_mat(i))
          
          end do
          
!           s = QcMatWrite_f( outer_next%data_mat(1), 'lof_inside', ASCII_VIEW)     
          
          c1_ctr = c1_ctr + cache%blks_triang_size

       ! Only two-el contribution
       else if (outer_next%num_dmat == 1) then
       
          allocate(outer_next%data_mat(cache%blks_triang_size*outer_contract_sizes_1(k)))

          do i = 1, cache%blks_triang_size*outer_contract_sizes_1(k)

             call QcMatInit(outer_next%data_mat(i), D_unp)
             call QcMatZero(outer_next%data_mat(i))
             call QcMatRAXPY(1.0d0, contrib_1(c1_ctr + i - 1), outer_next%data_mat(i))
          
          end do

          c1_ctr = c1_ctr + cache%blks_triang_size*outer_contract_sizes_1(k)
                   
       end if
       
       if (outer_next%next%dummy_entry) then
          traverse_end = .TRUE.
       end if
    
       k = k + 1
       outer_next => outer_next%next
    
    end do
       
    deallocate(outer_contract_sizes_1)
    
  end subroutine

  ! Do main part of perturbed S, D, F calculation at one order
  subroutine rsp_sdf_calculate(cache_outer, num_outer, size_i, &
  get_rsp_sol, get_ovl_mat, get_2el_mat, F, D, S, lof_cache, &
  rsp_eqn_retrieved, prog_info, rs_info)
  
    implicit none
    
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
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(contrib_cache) :: lof_cache
    type(contrib_cache_outer), target :: cache_outer
    type(contrib_cache_outer) :: F, D, S
    type(contrib_cache_outer), pointer :: cache_outer_next
    type(Qcmat), dimension(sum(size_i)) :: Dh, Dp, Fp, Sp, RHS, X
    type(Qcmat) :: A, B, C, T, U
    external :: get_rsp_sol, get_ovl_mat,  get_2el_mat
    
    ! Initialize all matrices
    
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
    
    do i = 1, sum(size_i)
    
       call QcMatInit(Dh(i), A)
       call QcMatZero(Dh(i))
       call QcMatInit(Dp(i), A)
       call QcMatZero(Dp(i))
       call QcMatInit(Fp(i), A)
       call QcMatZero(Fp(i))
       call QcMatInit(Sp(i), A)
       call QcMatZero(Sp(i))
       call QcMatInit(RHS(i), A)
       call QcMatZero(RHS(i))
       call QcMatInit(X(i), A)
       call QcMatZero(X(i))
    
    end do
    
    cache_outer_next => cache_outer
    
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
 
    ! Traverse
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
    
    
       ! For each cache element:
       ! Calculate Sp
       call p_tuple_to_external_tuple(pert, npert_ext, pert_ext)
       call get_ovl_mat(0, noc, 0, noc, npert_ext, pert_ext, &
                     size_i(k), Sp(ind_ctr:ind_ctr + size_i(k) - 1))
       
!        write(*,*) 'stage 1'
       
       
       deallocate(pert_ext)
       
       call contrib_cache_outer_add_element(S, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Sp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
       ! Add the initialized Dp to cache
       
       call contrib_cache_outer_add_element(D, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
      
       ! Assemble Fp (lower-order) for all components and add to cache
            
       call rsp_lof_recurse(pert, pert%npert, &
                            1, (/get_emptypert()/), .FALSE., lof_cache, size_i(k), &
                            Fp=Fp(ind_ctr:ind_ctr + size_i(k) - 1))
        
!        write(*,*) 'stage 2'
        
        
       call contrib_cache_outer_add_element(F, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
!         j = QcMatWrite_f(Fp(1), 'Fp_1_lo', ASCII_VIEW)
       
!        write(*,*) 'stage a'
       
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
          
          call rsp_get_matrix_z(superstructure_size, derivative_structure, &
               (/pert%npert,pert%npert/), pert%npert, &
               (/ (m, m = 1, pert%npert) /), pert%npert, &
               ind, F, D, S, Dp(ind_ctr + j - 1))
       
          call QcMatkABC(-1.0d0, Dp(ind_ctr + j - 1), B, A, T)
          call QcMatkABC(-1.0d0, A, B, Dp(ind_ctr + j - 1), U)
          call QcMatRAXPY(1.0d0, T, U)
          call QcMatRAXPY(1.0d0, U, Dp(ind_ctr + j - 1))
       
       end do
       
       call contrib_cache_outer_add_element(D, .FALSE., 1, & 
            (/pert/), data_size = size_i(k),  data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
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
    
!      j = QcMatWrite_f(Dp(1), 'Dp_after_Z', ASCII_VIEW)
    
! Debug printing kept for later use    
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
    
    
!     j = QcMatWrite_f(Fp(1), 'Fp_before_partic', ASCII_VIEW)
    
    k = 129
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
    
       ! Outside traversal:
       ! Complete Fp using Dp
       call get_2el_mat(0, noc, sum(size_i), Dp, sum(size_i), Fp)
    
!     end if

!     j = QcMatWrite_f(Fp(1), 'Fp_after_partic', ASCII_VIEW)
    
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
       
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
    
       ! Add the completed Fp to cache
       call contrib_cache_outer_add_element(F, .FALSE., 1, & 
            (/pert/), data_size = size_i(k),  data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
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
          
!           w = QcMatWrite_f(RHS(1), 'RHS_before', ASCII_VIEW)
          
          call rsp_get_matrix_y(superstructure_size, derivative_structure, &
                pert%npert, (/ (m, m = 1, pert%npert) /), &
                pert%npert, ind, F, D, S, RHS(ind_ctr + j - 1))
                
!           w = QcMatWrite_f(RHS(1), 'RHS_1', ASCII_VIEW)
          
!           stop
                
       
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
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next

    
! CONTINUE HERE: Retrieve and store rsp eqn solution vectors according to level 3 cache increment
    
    ! Traverse
    termination = .FALSE.
    do while(.NOT.(termination))

       if (size_i(k) > m) then
    
          do i = 1, size_i(k)/m + 1
    
             if (.NOT.((i - 1) *m >= size_i(k))) then
    
                first = (i - 1) * m + 1
                last = min(i * m, size_i(k))
                
                write(*,*) 'first, last, len', first, last, last - first + 1
!                 write(*,*) 'total ind first last', ind_ctr + first - 1, &
!                 ind_ctr + last - 1
                
!                 write(*,*) 'frequency passed:', real((/freq_sums(k)/))

                call prog_incr(prog_info, 3)
                
                if (rs_check(prog_info, rs_info, lvl=3)) then
                
                   write(*,*) ' '
                   write(*,*) 'RSP eqn solution batches were completed'
                   write(*,*) 'in previous invocation: Passing to next stage of calculation'
                   write(*,*) ' '
          
                   if (.NOT.(rsp_eqn_retrieved)) then
          
                      call mat_scal_retrieve(ind_ctr+last-1, 'OPENRSP_MAT_RSP', mat=X(1:ind_ctr+last-1))
                      rsp_eqn_retrieved = .TRUE.
             
                   end if

                else

          
                   !To Magnus: could you please check if the new callback subroutine works?
                   call get_rsp_sol(1,                                    &
                                    (/last-first+1/),                     &
                                    (/1/),                                &
                                    dcmplx(real((/freq_sums(k)/)),0.0d0), &
                                    RHS(ind_ctr+first-1:ind_ctr+last-1),  &
                                    X(ind_ctr+first-1:ind_ctr+last-1))
                   !call get_rsp_sol(1, dcmplx(real((/freq_sums(k)/)), 0.0d0), last - first + 1, &
                   !     RHS(ind_ctr + first - 1:ind_ctr + last - 1), &
                   !     X(ind_ctr + first - 1:ind_ctr + last - 1))
                   
                   call mat_scal_store(last - first + 1, 'OPENRSP_MAT_RSP', &
                        mat=X(ind_ctr+first-1:ind_ctr+last-1), start_pos = ind_ctr+first-1)
                   
                end if
   
             end if
       
          end do
       
       else
       
          if (rs_check(prog_info, rs_info, lvl=3)) then
                
             write(*,*) ' '
             write(*,*) 'RSP eqn solution batches were completed'
             write(*,*) 'in previous invocation: Passing to next stage of calculation'
             write(*,*) ' '
         
             if (.NOT.(rsp_eqn_retrieved)) then
          
                call mat_scal_retrieve(ind_ctr+last-1, 'OPENRSP_MAT_RSP', mat=X(1:ind_ctr+last-1))
                rsp_eqn_retrieved = .TRUE.
             
             end if

          else
       
             !To Magnus: could you please check if the new callback subroutine works?
             call get_rsp_sol(1,                                    &
                              (/size_i(k)/),                        &
                              (/1/),                                &
                              dcmplx(real((/freq_sums(k)/)),0.0d0), &
                              RHS(ind_ctr:ind_ctr+size_i(k)-1),     &
                              X(ind_ctr:ind_ctr+size_i(k)-1))
             ! "OLD NEW FORMAT": NOT SUPPOSED TO BE IN TRAVERSAL WHEN CHANGING TO NEW FORMAT
             !call get_rsp_sol(1, real((/freq_sums(k)/)), size_i(k), RHS(ind_ctr:ind_ctr + size_i(k) - 1), &
             !     X(ind_ctr:ind_ctr + size_i(k) - 1))
          
             call mat_scal_store(size_i(k), 'OPENRSP_MAT_RSP', &
                        mat=X(ind_ctr:ind_ctr+size_i(k)-1), start_pos = ind_ctr)
                   
          end if
          
          
    
       end if
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do
    
    do i = 1, size(RHS)
       call QcmatDst(RHS(i))
    end do
        
    ind_ctr = 1
    k = 1
            
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
       
    ! Traverse
    termination = .FALSE.
    do while(.NOT.(termination))

       ! For each cache element:
       ! Construct Dh for all components
    
       do j = 1, size_i(k)
    
          call QcMatkABC(-1.0d0, X(ind_ctr + j - 1), B, A, T)
          call QcMatkABC(1.0d0, A, B, X(ind_ctr + j - 1), U)
          call QcMatRAXPY(1.0d0, T, U)
          call QcMatAEqB(Dh(ind_ctr + j - 1), U)
    
       end do
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do
    
    k = 36
    
    if (size(Dh) > k) then

       do i = 1, size(Dh)/k + 1
       
          if (.NOT.((i - 1) *k >= size(Dh))) then
    
             first = (i - 1) * k + 1
             last = min(i * k, size(Dh))
    
             call get_2el_mat(0, noc, last - first + 1, Dh(first:last), last - first + 1, Fp(first:last))
    
          end if
          
       end do
    
    else
    
       ! Outside traversal
       ! Calculate Fh for all components using Dh
       call get_2el_mat(0, noc, sum(size_i), Dh, sum(size_i), Fp)
    
    end if
    
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
       
    ! Traverse
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
       
       ! Add together Dp and Dh and store
       ! Store Fp
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do
    
  end subroutine
  
  end module
