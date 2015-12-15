! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines for the calculation of perturbed
! overlap, density and Fock matrices used throughout the
! rsp_general calculation.

module rsp_perturbed_sdf

!  use matrix_defop, matrix => openrsp_matrix
!  use matrix_lowlevel, only: mat_init, mat_zero_like
  use rsp_contribs
  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_sdf_caching
  use rsp_lof_caching
  use rsp_property_caching
!  use interface_2el
  
  use qcmatrix_f
  
  implicit none


  public rsp_fds_2014
  public rsp_fds
  public get_fds_2014
  public get_fds_2014_mc
  public rsp_fock_lowerorder_2014
  public get_fock_lowerorder_2014

  private

  real(8) :: time_start
  real(8) :: time_end

  contains

  

  recursive subroutine rsp_fds_2014(pert, kn, F, D, S, get_rsp_sol, get_ovl_mat, &
                               get_1el_mat, get_2el_mat, get_xc_mat, id_outp)

    implicit none

    
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%npert) :: psub
    integer, dimension(2) :: kn
    integer :: i, j, k, id_outp
    type(SDF_2014) :: F, D, S
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat,  get_2el_mat, get_xc_mat



    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then (at final recursion level) get perturbed F, D, S 
    if (pert%npert > 1) then

       call make_p_tuple_subset(pert, psub)

       do i = 1, size(psub)

          if (sdf_already_2014(D, psub(i)) .eqv. .FALSE.) then

             call rsp_fds_2014(psub(i), kn, F, D, S, get_rsp_sol, get_ovl_mat, &
                          get_1el_mat, get_2el_mat, get_xc_mat, id_outp)

          end if

       end do       

    end if

    if (sdf_already_2014(D, pert) .eqv. .FALSE.) then
         
       if (kn_skip(pert%npert, pert%pid, kn) .eqv. .FALSE.) then

          write(id_outp,*) 'Calling ovlint/fock/density with labels ', pert%plab, &
                     ' and perturbation id ', pert%pid, ' with frequencies (real part)', &
                     real(pert%freq)
          write(id_outp,*) ' '
                 
          k = 1

          do j = 1, pert%npert

             pert%pid(j) = k
             k = k + 1

          end do

          call get_fds_2014(p_tuple_standardorder(pert), F, D, S, get_rsp_sol, &
                       get_ovl_mat, get_1el_mat, get_2el_mat, get_xc_mat, id_outp)

       else

!           write(*,*) 'Would have called ovlint/fock/density with labels ', &
!                      pert%plab, ' and perturbation id ', pert%pid, &
!                      ' but it was k-n forbidden'
!           write(*,*) ' '

       end if

    else

!        write(*,*) 'FDS for labels ', pert%plab, &
!                   'and perturbation id ', pert%pid, ' was found in cache'
!        write(*,*) ' '

    end if

  end subroutine
  
  

  
  
  
  
  
  recursive subroutine rsp_fds(n_props, n_freq_cfgs, p_tuples, kn_rule, F, D, S, get_rsp_sol, get_ovl_mat, &
                               get_1el_mat, get_2el_mat, get_xc_mat, dryrun, cache, id_outp)

    implicit none

    integer :: n_props
    integer, dimension(n_props) :: n_freq_cfgs
    type(p_tuple) :: p_test
    type(p_tuple), dimension(sum(n_freq_cfgs)) :: p_tuples
    type(p_tuple), allocatable, dimension(:) :: p_dummy_orders
    logical :: termination, dryrun
    integer, dimension(0) :: noc
    character(4), dimension(0) :: nof
    integer, dimension(sum(n_freq_cfgs), 2) :: kn_rule
    integer :: i, j, k, id_outp, max_order, max_npert
    integer, dimension(1) :: dfs
    integer, allocatable, dimension(:) :: size_i
    type(QcMat) :: Fp_dum, M_test
    type(contrib_cache_outer) :: F, D, S
    type(contrib_cache), target :: cache
    type(contrib_cache), pointer :: cache_next, lof_cache, lof_next
    type(contrib_cache_outer), pointer :: cache_outer_next
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat,  get_2el_mat, get_xc_mat

    max_order = max(maxval(kn_rule(:, 1)), maxval(kn_rule(:, 2)))
    
    max_npert = maxval((/(p_tuples(i)%npert, i = 1, sum(n_freq_cfgs))/))
    
!     write(*,*) 'max order', max_order
!     write(*,*) 'max npert', max_npert

       p_test%npert = 1
       allocate(p_test%pdim(1))
       allocate(p_test%plab(1))
       allocate(p_test%pid(1))
       allocate(p_test%freq(1))
       p_test%pdim = (/12/)
       p_test%plab = (/'GEO '/)
       p_test%pid = (/1/)
       p_test%freq = (/0.0/)

    
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


   k = 1 
    ! Recurse to identify all necessary perturbed F, D, S
    do i = 1, n_props
       do j = 1, n_freq_cfgs(i)
!             write(*,*) 'kn rules', kn_rule(k, :)
       
          call rsp_fds_recurse(p_tuples(k), kn_rule(k, :), max_npert, p_dummy_orders, cache, id_outp)
       
          k = k + 1
       
       end do
    end do
           
    cache_next => cache
    

    
    !   For each order of perturbation identified (lowest to highest):
    do i = 1, max_order
    
       call contrib_cache_allocate(lof_cache)

       ! Cycle until order reached
       do while(.NOT.(cache_next%p_inner%freq(1) == 1.0*i))
          cache_next => cache_next%next
!           write(*,*) 'inner pdim, num outer cycle', cache_next%p_inner%freq(1), cache_next%num_outer
       end do
       
!        write(*,*) 'num outer', cache_next%num_outer
!        write(*,*) 'outer example', cache_next%contribs_outer%last
       
       ! Contains number of components of perturbed matrices for each perturbation
       
       allocate(size_i(cache_next%num_outer))       
       k = 1
       
       ! Cycle until at start of outer cache
       cache_outer_next => contrib_cache_outer_cycle_first(cache_next%contribs_outer)
       cache_outer_next => cache_outer_next%next
       
       write(*,*) 'outer entry', cache_outer_next%p_tuples(1)%plab
             
       ! Traverse all elements of outer cache of present cache element
       termination = .FALSE.
       do while(.NOT.(termination))
       
          ! Recurse to identify lower-order Fock matrix contributions
          ! The p_tuples attribute should always be length 1 here, so OK to take the first element
          call rsp_lof_recurse(cache_outer_next%p_tuples(1), cache_outer_next%p_tuples(1)%npert, &
                                      1, (/get_emptypert()/), .TRUE., lof_cache, 1, (/Fp_dum/))
       
       write(*,*) 'lof recursion done'
       
          ! Get number of perturbed matrices for this tuple
          size_i(k) = cache_outer_next%blks_tuple_triang_size(1)
!           write(*,*) 'size i', size_i(k)
          k = k + 1

          termination = cache_outer_next%last
          cache_outer_next => cache_outer_next%next
          
       end do

       lof_next => lof_cache
       
       ! Cycle lower-order Fock cache until at start
       do while(.NOT.(lof_next%last))
          lof_next => lof_next%next
       end do
       lof_next => lof_next%next
       lof_next => lof_next%next
       
              write(*,*) 'lof calculation starting'
       
       ! Traverse cache and precalculate elements
       termination = .FALSE.
       do while (.NOT.(termination))
       
          call rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, lof_next)
          
          termination = (lof_next%last)
          lof_next => lof_next%next
          
       end do
       
       call rsp_sdf_calculate(cache_outer_next, cache_next%num_outer, size_i,&
            get_rsp_sol, get_ovl_mat,  get_2el_mat, F, D, S, lof_next)
       
!        call QcMatInit(Fp_dum)

!        write(*,*) 'getting sample matrix'
!        dfs = (/4/)
!        
!        call contrib_cache_getdata_outer(D, 1, (/p_test/), .FALSE., &
!                   contrib_size=1, ind_len=1, ind_unsorted=dfs, &
!                   mat_sing=Fp_dum)
       
!        j = QcMatWrite_f(Fp_dum, 'Fp_dum', ASCII_VIEW)
       
!        write(*,*) 'got sample matrix'
!        stop
       
       deallocate(size_i)
       
       deallocate(lof_cache)
          
    end do
    
!     deallocate(lof_cache)
    
    
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
    ! Then (at final recursion level) get perturbed F, D, S 
    if (pert%npert > 1) then

       call make_p_tuple_subset(pert, psub)
       
       

       do i = 1, size(psub)
       
!           write(*,*) 'psub npert', psub(i)%npert
!           write(*,*) 'psub dummy npert', p_dummy_orders(psub(i)%npert)%npert
       

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

       else

!           write(*,*) 'Would have identified necessary ovlint/fock/density with labels ', &
!                      pert%plab, ' and perturbation id ', pert%pid, &
!                      ' but it was k-n forbidden'
!           write(*,*) ' '

       end if

    else

!        write(*,*) 'FDS identifier for labels ', pert%plab, &
!                   'and perturbation id ', pert%pid, ' was found in cache'
!        write(*,*) ' '

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

    else

!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)

       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%npert >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then
       
!           write(*,*) 'num p tuples', num_p_tuples
!           write(*,*) 'plab', p_tuples(1)%plab

          if (contrib_cache_already(fock_lowerorder_cache, &
          num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples))) then

             if (.NOT.(dryrun)) then
             
                do i = 1, num_p_tuples
                   if (i == 1) then
                      write(*,*) 'F ', p_tuples(i)%plab
                   else
                      write(*,*) 'D ', p_tuples(i)%plab
                   end if
                end do
             
                write(*,*) 'Getting lower-order perturbed Fock order contribution from cache'
                call contrib_cache_getdata(fock_lowerorder_cache, num_p_tuples, p_tuples, &
                contrib_size=fp_size, ind_len=0, mat=Fp)
             
             else
             
                do i = 1, num_p_tuples
                   if (i == 1) then
                      write(*,*) 'F', p_tuples(i)%pid
                   else
                      write(*,*) 'D', p_tuples(i)%pid
                   end if
                end do
             
                write(*,*) 'Identified lower-order perturbed Fock contribution already in cache'
                
             end if   
          
          else

             if (dryrun) then
          
                write(*,*) 'Identified perturbed Fock matrix lower order contribution'

                do i = 1, num_p_tuples
                   if (i == 1) then
                      write(*,*) 'F', p_tuples(i)%pid
                   else
                      write(*,*) 'D', p_tuples(i)%pid
                   end if
                end do
                call contrib_cache_add_element(fock_lowerorder_cache, num_p_tuples, p_tuples)
                
                
             else
             
                write(*,*) 'ERROR: Wanted cache element that was not calculated'
                
             end if
          
          end if

       else

!           write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
!           write(*,*) 'at order for which perturbed D is to be found '
!           write(*,*) ' '

       end if

    end if

  end subroutine
  
  subroutine rsp_lof_calculate(D, get_1el_mat, get_ovl_mat, get_2el_mat, cache)

    implicit none

    logical :: traverse_end
    integer :: cache_offset, i, j, k, m, n 
    integer :: id_outp
    integer :: total_outer_size_1, c1_ctr, lhs_ctr_1, num_pert
    integer :: num_0, num_1
    character(30) :: mat_str, fmt_str, fmt_str2
    type(contrib_cache) :: cache
    type(contrib_cache_outer) :: D
    type(contrib_cache_outer), pointer :: outer_next
    type(p_tuple) :: t_mat_p_tuple, t_matrix_bra, t_matrix_ket
    type(QcMat), allocatable, dimension(:) :: LHS_dmat_1, contrib_0, contrib_1
    type(QcMat) :: D_unp
    integer, allocatable, dimension(:) :: outer_contract_sizes_1, outer_contract_sizes_1_coll
    integer, allocatable, dimension(:) :: pert_ext
    external :: get_1el_mat, get_ovl_mat, get_2el_mat
    
    write(*,*) 'calculating for perturbation', cache%p_inner%plab
    
    ! Assume indices for inner, outer blocks are calculated earlier
    
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
 
!     ! Make collapsed contraction sizes array for 1-el call
!  
!     allocate(outer_contract_sizes_1_coll(num_1))
!     
!     k = 1 
!      do i = 1, cache%num_outer
!         if (outer_contract_sizes_1(i) > 0 then
!            outer_contract_sizes_1_coll(k) = outer_contract_sizes_1(i)
!           k = k + 1
!         end if
!     end do
    
    ! Allocate and set up outer
    
!     write(*,*) 'outer contract sizes_1', outer_contract_sizes_1
!     if (sum(outer_contract_sizes_1) > 10000) then
!     stop
!     
!     end if
!     
    
    allocate(LHS_dmat_1(sum(outer_contract_sizes_1)))
    
    call QCMatInit(D_unp)
! write(*,*) 'A'
    call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .FALSE., &
         contrib_size=1, ind_len=1, ind_unsorted=(/1/), mat_sing=D_unp)

! write(*,*) 'b'


         
    do i = 1, size(LHS_dmat_1)
    
      call QCMatInit(LHS_dmat_1(i), D_unp)
    
    end do

!     write(*,*) 'C'
    
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
       
    k = 1
    lhs_ctr_1 = 1
    
    do while (traverse_end .EQV. .FALSE.)

!        ! No chain rule applications
!        if (outer_next%num_dmat == 0) then
! 
!           call contrib_cache_getdata_outer(D, 1, (/get_emptypert()/), .TRUE., &
!                1, ind_unsorted=(/1/), mat=D_unp)

       ! One chain rule application
       if (outer_next%num_dmat == 1) then
       
!              write(*,*) 'D', outer_next%blks_tuple_triang_size 
!              write(*,*) 'outer_next%indices', outer_next%indices
!              write(*,*) 'outer_next%indices(m, :)', outer_next%indices(m, :)
       
          do m = 1, outer_contract_sizes_1(k) 
             call contrib_cache_getdata_outer(D, 1, (/outer_next%p_tuples(1)/), .FALSE., &
                  contrib_size=1, ind_len=1, ind_unsorted=outer_next%indices(m, :), &
                  mat_sing=LHS_dmat_1(lhs_ctr_1 + m  - 1))
                  write(*,*) 'indices', outer_next%indices(m, :)
                  write(*,*) 'got matrix', lhs_ctr_1 + m  - 1
    j = QcMatWrite_f(LHS_dmat_1(lhs_ctr_1 + m  - 1), 'LHS_dmat_1', ASCII_VIEW)
          end do
       
!              write(*,*) 'D2' 
       end if
   
       if (outer_next%next%dummy_entry) then
          traverse_end = .TRUE.
       end if

       lhs_ctr_1 = lhs_ctr_1 + outer_contract_sizes_1(k)
       k = k + 1
       
       outer_next => outer_next%next
    
    end do
    
!     write(*,*) 'blks triang size', cache%blks_triang_size
    
       if (cache%blks_triang_size > 10000) then
    stop
    
    end if
    
    
    allocate(contrib_0(cache%blks_triang_size))
    allocate(contrib_1(cache%blks_triang_size*total_outer_size_1))
    

    
    do i = 1, size(contrib_0)
!       write(*,*) 'i is', i
      call QCMatInit(contrib_0(i), D_unp)
      call QCMatZero(contrib_0(i))
!       j = QcMatWrite_f(contrib_0(i), 'contrib 0', ASCII_VIEW)
    end do
    
    do i = 1, size(contrib_1)
      call QCMatInit(contrib_1(i), D_unp)
      call QCMatZero(contrib_1(i))
    end do

    
    ! Calculate contributions
    
    
    ! Calculate one-electron contributions
    if (num_0 > 0) then
    

!         write(*,*) 'num pert', num_pert
!         write(*,*) 'pert ext', pert_ext
!         write(*,*) 'num 0', num_0
    
       call get_1el_mat(num_pert, pert_ext, size(contrib_0), contrib_0)
      
       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()
      
!        call rsp_ovlave_t_matrix(get_ovl_mat, cache%p_inner, cache%p_inner%npert, &
!                                 t_matrix_bra, t_matrix_ket, 1, &
!                                 D_unp, size(contrib_0), contrib_0)
    
    end if
    
!            j = QcMatWrite_f(contrib_0(1), 'contrib 0_1', ASCII_VIEW)
    
    write(*,*) 'num pert', num_pert
    write(*,*) 'pert_ext', pert_ext
    write(*,*) 'LHS dmat 1 size', size(LHS_dmat_1)
    write(*,*) 'size contrib 1', size(contrib_1)
    
    
    if ((size(LHS_dmat_1) == 15) .AND. (size(contrib_1) == 180)) then
    
!        do i = 1, 15
!        
!           
!        
!           if (i < 10) then
!           
!              fmt_str = "(A8, I1)"
!           
!           else if (i < 100) then
!           
!              fmt_str = "(A8, I2)"
!           
!           else
!           
!              fmt_str = "(A8, I3)"
!           
!           end if
!           
!           write(mat_str, fmt_str) 'LHS_dmat_', i
!           
!           write(*,*) 'i', i
!           write(*,*) 'fname:', mat_str
!           
!           
!           
!           j = QcMatWrite_f(LHS_dmat_1(i), trim(mat_str), ASCII_VIEW)
!           
!           do k = 1, 12
!           
!           if ((i-1)* 12 + k < 10) then
!           
!              fmt_str = "(A8, I1)"
!           
!           else if ((i-1)* 12 + k < 100) then
!           
!              fmt_str = "(A8, I2)"
!           
!           else
!           
!              fmt_str = "(A8, I3)"
!           
!           end if
!           
!           write(mat_str, fmt_str) 'contrib_', (i-1) * 12 + k
!           write(*,*) 'k', k
!           write(*,*) 'fname:', mat_str
!           
!           j = QcMatWrite_f(contrib_1((i-1)*12 + k), trim(mat_str), ASCII_VIEW)
!           
!           end do
!           
!           
!           
!           
!           
!                  
!        end do
    
    

    
    end if
    
    ! Calculate two-electron contributions
    call get_2el_mat(num_pert, pert_ext, size(LHS_dmat_1), LHS_dmat_1, &
    size(contrib_1), contrib_1)
                       
                       
    ! Traversal: Add nuc-nuc, 1-el and two-el contributions together
    
    traverse_end = .FALSE.
    
    outer_next = contrib_cache_outer_cycle_first(outer_next)
    outer_next => outer_next%next
      
    k = 1
    
    c1_ctr = 1
    
    
    write(*,*) 'Inner perturbation ', cache%p_inner%plab
    write(*,*) 'Inner perturbation size ', cache%blks_triang_size
    write(*,*) ' '
    
    do while (traverse_end .EQV. .FALSE.)
    
       do i = 1, outer_next%num_dmat
       write(*,*) 'Outer perturbation(s) ', i, ' are ', outer_next%p_tuples(i)%plab
       
       end do
    
       
       write(*,*) 'Outer size ', outer_contract_sizes_1(k)

  
       ! One-el and two-el contributions
       if (outer_next%num_dmat == 0) then
       
          allocate(outer_next%data_mat(cache%blks_triang_size*outer_contract_sizes_1(k)))
          
          do i = 1, cache%blks_triang_size*outer_contract_sizes_1(k)
!           write(*,*) 'i is', i

             call QcMatInit(outer_next%data_mat(i), D_unp)
             call QcMatZero(outer_next%data_mat(i))
!              write(*,*) 'A'
             call QcMatkAB(1.0d0, contrib_0(i), contrib_1(c1_ctr + i - 1), outer_next%data_mat(i))
!              write(*,*) 'B'
!            j = QcMatWrite_f(contrib_0(1), 'contrib_0', ASCII_VIEW)
!            j = QcMatWrite_f(contrib_1(1), 'contrib_1', ASCII_VIEW)
!            j = QcMatWrite_f(outer_next%data_mat(1), 'lof_1', ASCII_VIEW)
          
          end do
          
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
       
       write(*,*) 'Size of data', size(outer_next%data_mat)
       write(*,*) ' '
   
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
  get_rsp_sol, get_ovl_mat,  get_2el_mat, F, D, S, lof_cache)
  
    implicit none
    
    logical :: termination
    integer :: num_outer, ind_ctr, npert_ext, sstr_incr, superstructure_size
    integer :: i, j, k, m, nblks
    integer :: first, last
    integer, dimension(0) :: noc
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
    
    ! Get perturbed overlap matrices
!     call rsp_ovl_mat(pert, size_i(k), A, S)
    
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
       
       deallocate(pert_ext)
       
       call contrib_cache_outer_add_element(S, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Sp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
       ! Add the initialized Dp to cache
       
       call contrib_cache_outer_add_element(D, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
      
       ! Assemble Fp (lower-order) for all components and add to cache
             
       write(*,*) 'About to get lof data'
             
       call rsp_lof_recurse(pert, pert%npert, &
                            1, (/get_emptypert()/), .FALSE., lof_cache, size_i(k), &
                            Fp=Fp(ind_ctr:ind_ctr + size_i(k) - 1))
                            
       write(*,*) 'Returned from lof'
                            
!               j = QcMatWrite_f(Fp(1), 'Fp 1 a', ASCII_VIEW)
!        j = QcMatWrite_f(Fp(13), 'Fp 13 a', ASCII_VIEW)                            
       call contrib_cache_outer_add_element(F, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
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
       
       

       
       call contrib_cache_outer_add_element(F, .FALSE., 1, & 
            (/pert/), data_size = size_i(k),  data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
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
    
!            j = QcMatWrite_f(Dp(1), 'Dp 1', ASCII_VIEW)
!        j = QcMatWrite_f(Dp(13), 'Dp 13', ASCII_VIEW)
!               j = QcMatWrite_f(Fp(1), 'Fp 1 before', ASCII_VIEW)
!        j = QcMatWrite_f(Fp(13), 'Fp 13 before', ASCII_VIEW)
    
    write(*,*) 'size of Dp, Fp', size(Dp), size(Fp)
    
    
!     do i = 1, size(Dp)
!     
!     
!           if (i < 10) then
!           
!              fmt_str = "(A3, I1)"
!           
!           else if (i < 100) then
!           
!              fmt_str = "(A3, I2)"
!           
!           else
!           
!              fmt_str = "(A3, I3)"
!           
!           end if
!           
!           write(mat_str, fmt_str) 'Dp_', i
!           
!           write(*,*) 'i', i
!           write(*,*) 'fname:', mat_str
!           
!           
!           
!           j = QcMatWrite_f(Dp(i), trim(mat_str), ASCII_VIEW)
!     
!     end do
    
    k = 36
    
    if (size(Dp) > k) then
    
    write(*,*) 'Dp multiples', size(Dp)/k + 1
    
    do i = 1, size(Dp)/k + 1
    
       if (.NOT.((i - 1) *k >= size(Dp))) then
    
          first = (i - 1) * k + 1
          last = min(i * k, size(Dp))
    
          call get_2el_mat(0, noc, last - first + 1, Dp(first:last), last - first + 1, Fp(first:last))
    
       end if
       
    end do
    
    
    
    
    else
    
    ! Outside traversal:
    ! Complete Fp using Dp
    call get_2el_mat(0, noc, sum(size_i), Dp, sum(size_i), Fp)
    
    end if

    
!            j = QcMatWrite_f(Fp(1), 'Fp 1', ASCII_VIEW)
!        j = QcMatWrite_f(Fp(13), 'Fp 13', ASCII_VIEW)
    
    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
       
    ! Traverse
    termination = .FALSE.
    do while(.NOT.(termination))

       pert = cache_outer_next%p_tuples(1)
    
    
       ! Set up block info
       nblks = get_num_blks(pert)

       allocate(blk_info(nblks, 3))
       allocate(blk_sizes(pert%npert))
       blk_info = get_blk_info(nblks, pert)
       blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))
    
       ! For each cache element:
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
          
          call rsp_get_matrix_y(superstructure_size, derivative_structure, &
                pert%npert, (/ (m, m = 1, pert%npert) /), &
                pert%npert, ind, F, D, S, RHS(ind_ctr + j - 1))
       
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
    
!           j = QcMatWrite_f(RHS(1), 'RHS 1', ASCII_VIEW)
!           j = QcMatWrite_f(X(1), 'X 1 before', ASCII_VIEW)
          
!           write(*,*) 'num outer', num_outer
!           write(*,*) 'freq sums', freq_sums
!           write(*,*) 'size i', size_i
!           write(*,*) 'sum size i', sum(size_i)
!           write(*,*) 'size of rhs', size(RHS)
!           write(*,*) 'size of X', size(RHS)
          

          m  = 1

    ind_ctr = 1
    k = 1
    
    ! Cycle to start
    cache_outer_next = contrib_cache_outer_cycle_first(cache_outer_next)
    cache_outer_next => cache_outer_next%next
       
    ! Traverse
    termination = .FALSE.
    do while(.NOT.(termination))

!        write(*,*) 'first component', ind_ctr
!     
!        write(*,*) 'freq sums', freq_sums(k)
!        write(*,*) 'size i', size_i(k)
! j = QcMatWrite_f(RHS(1), 'RHS_1', ASCII_VIEW)
!     j = QcMatWrite_f(RHS(13), 'RHS_13', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(15), 'RHS_15', ASCII_VIEW)
!         
!         if (size(RHS) > 15) then
!         
!         
!         j = QcMatWrite_f(RHS(16), 'RHS_16', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(17), 'RHS_17', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(18), 'RHS_18', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(19), 'RHS_19', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(20), 'RHS_20', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(21), 'RHS_21', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(22), 'RHS_22', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(23), 'RHS_23', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(24), 'RHS_24', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(25), 'RHS_25', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(26), 'RHS_26', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(27), 'RHS_27', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(28), 'RHS_28', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(29), 'RHS_29', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(30), 'RHS_30', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(31), 'RHS_31', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(32), 'RHS_32', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(33), 'RHS_33', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(34), 'RHS_34', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(35), 'RHS_35', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(36), 'RHS_36', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(37), 'RHS_37', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(38), 'RHS_38', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(39), 'RHS_39', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(40), 'RHS_40', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(41), 'RHS_41', ASCII_VIEW)
!         j = QcMatWrite_f(RHS(42), 'RHS_42', ASCII_VIEW)
!         
!         
!         end if
!         j = QcMatWrite_f(X(13), 'X i1 before', ASCII_VIEW)

       
       
    if (size_i(k) > m) then
    
    write(*,*) 'Rsp multiples', size_i(k)/m + 1
!     
    do i = 1, size_i(k)/m + 1
    
       if (.NOT.((i - 1) *m >= size_i(k))) then
    
          first = (i - 1) * m + 1
          last = min(i * m, size_i(k))
          
          write(*,*) 'first, last: ', first, last

          
           call get_rsp_sol(1, real((/freq_sums(k)/)), last - first + 1, RHS(ind_ctr + first - 1:ind_ctr + last - 1), &
                X(ind_ctr + first - 1:ind_ctr + last - 1))
          
!           call get_2el_mat(0, noc, last - first + 1, Dp(first:last), last - first + 1, Fp(first:last))
    
       end if
       
    end do
       
    else
          
       ! "OLD NEW FORMAT": NOT SUPPOSED TO BE IN TRAVERSAL WHEN CHANGING TO NEW FORMAT
       call get_rsp_sol(1, real((/freq_sums(k)/)), size_i(k), RHS(ind_ctr:ind_ctr + size_i(k) - 1), &
            X(ind_ctr:ind_ctr + size_i(k) - 1))
    
    end if
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do
    

    
    j = QcMatWrite_f(X(1), 'X_1_after', ASCII_VIEW)
    j = QcMatWrite_f(X(13), 'X_13_after', ASCII_VIEW)
     j = QcMatWrite_f(X(15), 'X_15_after', ASCII_VIEW)
       
    do i = 1, size(RHS)
!     write(*,*) 'i dst', i
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
    
    write(*,*) 'Dh multiples', size(Dh)/k + 1
    
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
       
       write(*,*) 'Pert: ', pert%plab, pert%freq
       
       
       do j = 1, size_i(k)
    
!             m = QcMatWrite_f(Dp(ind_ctr + j - 1), 'Dp_before', ASCII_VIEW)
            
          call QcMatRAXPY(1.0d0, Dh(ind_ctr + j - 1), Dp(ind_ctr + j - 1))
          
!           m = QcMatWrite_f(Dp(ind_ctr + j - 1), 'Dp_after', ASCII_VIEW)
    
       end do
    
           call contrib_cache_outer_add_element(F, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Fp(ind_ctr:ind_ctr + size_i(k) - 1) )
    
       call contrib_cache_outer_add_element(D, .FALSE., 1, & 
            (/pert/), data_size = size_i(k), data_mat = Dp(ind_ctr:ind_ctr + size_i(k) - 1) )
       
!        if(pert%npert ==1) then
!         if(pert%plab(1) .eq. 'GEO ') then
!               stop
!          end if
!        end if
       ! For each cache element:
       ! Add together Dp and Dh and store
       ! Store Fp
    
       ind_ctr = ind_ctr + size_i(k)
       k = k + 1
    
       termination = cache_outer_next%last
       cache_outer_next => cache_outer_next%next
       
    end do
    
!     do i = 1, size(Dp)
!     
!     
!           if (i < 10) then
!           
!              fmt_str = "(A3, I1)"
!           
!           else if (i < 100) then
!           
!              fmt_str = "(A3, I2)"
!           
!           else
!           
!              fmt_str = "(A3, I3)"
!           
!           end if
!           
!           write(mat_str, fmt_str) 'Dp_', i
!           
!           write(*,*) 'i', i
!           write(*,*) 'fname:', mat_str
!           
!           
!           
!           j = QcMatWrite_f(Dp(i), trim(mat_str), ASCII_VIEW)
!     
!     end do
    
!     stop
    
    
    
  end subroutine
  
  
  
  
  
  
  
  
  

    subroutine get_fds_2014(pert, F, D, S, get_rsp_sol, get_ovl_mat, &
                       get_1el_mat, get_2el_mat, get_xc_mat, id_outp)

use qcmatrix_f
!    use interface_rsp_solver, only: rsp_solver_exec
    implicit none

    
    integer :: sstr_incr, i, j, superstructure_size, nblks, perturbed_matrix_size, id_outp
    integer :: ierr, npert_ext
    integer, allocatable, dimension(:) :: ind, blk_sizes, pert_ext
    integer, allocatable, dimension(:,:) :: blk_info, indices
    integer, dimension(0) :: noc
    character(4), dimension(0) :: nof
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(SDF_2014) :: F, D, S
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat, get_2el_mat, get_xc_mat
    type(QcMat) :: X(1), RHS(1), A, B, C, zeromat, T, U
    type(QcMat), allocatable, dimension(:) :: Fp, Dp, Sp, Dh
    type(f_l_cache_2014), pointer :: fock_lowerorder_cache


    ! ASSUME CLOSED SHELL
    
    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)
    call QcMatInit(T)
    call QcMatInit(U)
    

    call sdf_getdata_s_2014(D, get_emptypert(), (/1/), A)
    call sdf_getdata_s_2014(S, get_emptypert(), (/1/), B)
    call sdf_getdata_s_2014(F, get_emptypert(), (/1/), C)
    
    

    nblks = get_num_blks(pert)

    allocate(blk_info(nblks, 3))
    allocate(blk_sizes(pert%npert))

    blk_info = get_blk_info(nblks, pert)
    perturbed_matrix_size = get_triangulated_size(nblks, blk_info)
    blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

    allocate(Fp(perturbed_matrix_size))
    allocate(Dp(perturbed_matrix_size))
    allocate(Sp(perturbed_matrix_size))
    allocate(Dh(perturbed_matrix_size))

    ! Process perturbation tuple for external call
    
    call p_tuple_to_external_tuple(pert, npert_ext, pert_ext)
    
    
    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix


    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Sp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Sp(i), A)

    end do
! write(*,*) 'Sp a', Sp(1)%elms

    call get_ovl_mat(0, noc, 0, noc, npert_ext, pert_ext, &
                     perturbed_matrix_size, Sp)

!     call rsp_ovlint(zeromat%nrow, pert%npert, pert%plab, &
!                        (/ (1, j = 1, pert%npert) /), pert%pdim, &
!                        nblks, blk_info, blk_sizes, &
!                        perturbed_matrix_size, Sp)

! write(*,*) 'Sp b', Sp(1)%elms

    call sdf_add_2014(S, pert, perturbed_matrix_size, Sp)

    deallocate(blk_sizes)

    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp

    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Dp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Dp(i), Sp(1))

!        call mat_init(Dh(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Dh(i), Sp(1))

!        call mat_init(Fp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Fp(i), Sp(1))

    end do

    call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

! write(*,*) 'Fp a', Fp(1)%elms


    call f_l_cache_allocate_2014(fock_lowerorder_cache)
    call rsp_fock_lowerorder_2014(pert, pert%npert, 1, (/get_emptypert()/), &
                         get_1el_mat, get_ovl_mat, get_2el_mat, get_xc_mat, 0, D, &
                         perturbed_matrix_size, Fp, fock_lowerorder_cache)

! write(*,*) 'Fp b', Fp(1)%elms
    deallocate(fock_lowerorder_cache)

    call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)

    ! b) For Dp: Create differentiation superstructure: First dryrun for size, and
    ! then the actual superstructure call

    superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%npert, pert%npert/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

    allocate(derivative_structure(superstructure_size, 3))
    allocate(indices(perturbed_matrix_size, pert%npert))
    allocate(ind(pert%npert))

    call derivative_superstructure(pert, (/pert%npert, &
         pert%npert/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)
    call make_triangulated_indices(nblks, blk_info, perturbed_matrix_size, indices)

    do i = 1, size(indices, 1)
ierr = QcMatWrite_f(Fp(i), "Fp_a", ASCII_VIEW)

       ind = indices(i, :)

! write(*,*) 'Dp 0', Dp(1)%elms

       call rsp_get_matrix_z_2014(superstructure_size, derivative_structure, &
               (/pert%npert,pert%npert/), pert%npert, &
               (/ (j, j = 1, pert%npert) /), pert%npert, &
               ind, F, D, S, Dp(i))
ierr = QcMatWrite_f(Dp(i), "Dp_a", ASCII_VIEW)

! write(*,*) 'Dp 1', Dp(1)%elms


!      Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A
       call QcMatkABC(-1.0d0, Dp(i), B, A, T)
       call QcMatkABC(-1.0d0, A, B, Dp(i), U)
       call QcMatRAXPY(1.0d0, T, U)
       call QcMatRAXPY(1.0d0, U, Dp(i))

ierr = QcMatWrite_f(Dp(i), "Dp_b", ASCII_VIEW)
       
       
! write(*,*) 'Dp 2', Dp(1)%elms


       call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

       ! 3. Complete the particular contribution to Fp
! write(*,*) 'Fp b2', Fp(1)%elms
       call cpu_time(time_start)
!get_2el_mat(npert_ext, pert_ext, 1, (/D_unp/), size(Fp), Fp)

       call get_2el_mat(0, noc, 1, (/Dp(i)/), 1, Fp(i:i))
!        call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i), &
!                           1, Fp(i:i))
ierr = QcMatWrite_f(Fp(i), "Fp_b", ASCII_VIEW)

       call cpu_time(time_end)
!        print *, 'seconds spent in 2-el particular contribution', time_end - time_start
! write(*,*) 'Fp b3', Fp(1)%elms

       call cpu_time(time_start)
       ! MaR: Reintroduce once callback functionality is complete
!        call get_xc_mat(0, pert%pdim, noc, nof, 2, (/A, Dp(i)/), Fp(i:i))
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dp(i) /) , 1, Fp(i:i))
       call cpu_time(time_end)
!       print *, 'seconds spent in XC particular contribution', time_end - time_start

       ! MaR: Considered to be part of the 2el call
!        call cpu_time(time_start)
!        call rsp_pe(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i) , 1, Fp(i))
!        call cpu_time(time_end)
! !       print *, 'seconds spent in PE particular contribution', time_end - time_start

! write(*,*) 'Fp b4', Fp(1)%elms

       call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)
! write(*,*) 'Fp c', Fp(1)%elms


       ! 4. Make right-hand side using Dp

    call QcMatInit(RHS(1))
    call QcMatInit(X(1))

       call rsp_get_matrix_y_2014(superstructure_size, derivative_structure, &
                pert%npert, (/ (j, j = 1, pert%npert) /), &
                pert%npert, ind, F, D, S, RHS(1))

ierr = QcMatWrite_f(Dp(i), "Dp", ASCII_VIEW)
ierr = QcMatWrite_f(RHS(1), "RHS", ASCII_VIEW)
ierr = QcMatWrite_f(Fp(i), "Fp", ASCII_VIEW)

       ! Note (MaR): Passing only real part of freq. Is this OK?
       ! MaR: May need to vectorize RHS and X
       
       call get_rsp_sol(1, (/sum(real(pert%freq(:)))/), 1, RHS, X)
       
!        call get_rsp_sol(RHS(1), 1, (/sum(real(pert%freq(:)))/), X)
!        call rsp_solver_exec(RHS(1), (/sum(real(pert%freq(:)))/), X)

       call QcMatDst(RHS(1))

       ! 5. Get Dh using the rsp equation solution X
       
!        Dh(i) = A*B*X(1) - X(1)*B*A
       call QcMatkABC(-1.0d0, X(1), B, A, T)
       call QcMatkABC(1.0d0, A, B, X(1), U)
       call QcMatRAXPY(1.0d0, T, U)
       call QcMatAEqB(Dh(i), U)
ierr = QcMatWrite_f(A, "A", ASCII_VIEW)
ierr = QcMatWrite_f(B, "B", ASCII_VIEW)
ierr = QcMatWrite_f(X(1), "X", ASCII_VIEW)

       ! 6. Make homogeneous contribution to Fock matrix

       call cpu_time(time_start)
       call get_2el_mat(0, noc, 1, (/Dh(i)/), 1, Fp(i:i))
ierr = QcMatWrite_f(Dh(i), "Dh_c", ASCII_VIEW)
ierr = QcMatWrite_f(Fp(i), "Fp_c", ASCII_VIEW)

       call cpu_time(time_end)
!        print *, 'seconds spent in 2-el homogeneous contribution', time_end - time_start
! write(*,*) 'Fp b3', Fp(1)%elms

       call cpu_time(time_start)
       ! MaR: Reintroduce once callback functionality is complete
!        call get_xc_mat(0, pert%pdim, noc, nof, 2, (/A, Dp(i)/), Fp(i:i))
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dh(i) /) , 1, Fp(i:i))
       call cpu_time(time_end)
!       print *, 'seconds spent in XC homogeneous contribution', time_end - time_start
       
       
       
!        call cpu_time(time_start)
!        call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), &
!                           1, Fp(i:i))
!        call cpu_time(time_end)
!        print *, 'seconds spent in 2-el homogeneous contribution', time_end - time_start

!        call cpu_time(time_start)
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dh(i) /) , 1, Fp(i:i))
!        call cpu_time(time_end)
!        print *, 'seconds spent in XC homogeneous contribution', time_end - time_start

!        call cpu_time(time_start)
!        call rsp_pe(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), 1, Fp(i))
!        call cpu_time(time_end)
!       print *, 'seconds spent in PE homogeneous contribution', time_end - time_start

       ! 7. Complete perturbed D with homogeneous part

!        Dp(i) = Dp(i) + Dh(i)
       call QcMatRAXPY(1.0d0, Dh(i), Dp(i))


if (perturbed_matrix_size < 10) then

write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'
  

else

if (mod(i, perturbed_matrix_size/10) == 1) then

write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'

end if

end if


! ierr = QcMatWrite_f(Dp(i), 'Dp_drop', ASCII_VIEW)
! ierr = QcMatWrite_f(Fp(i), 'Fp_drop', ASCII_VIEW)
! ierr = QcMatWrite_f(Sp(i), 'Sp_drop', ASCII_VIEW)
! if (i == 12) then
! stop
! 
! end if
!        write(*,*) ' '
!        write(*,*) 'Finally, Dp is:'
!        write(*,*) Dp(i)%elms
!        write(*,*) ' '
!        write(*,*) 'Finally, Fp is:'
!        write(*,*) Fp(i)%elms
!        write(*,*) ' '
!        write(*,*) 'Finally, Sp is:'
!        write(*,*) Sp(i)%elms_alpha
!        write(*,*) ' '

    end do

    ! Add the final values to cache

    call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)
    call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

    do i = 1, size(indices, 1)

       call QcMatDst(Dh(i))
       call QcMatDst(Dp(i))
       call QcMatDst(Fp(i))
       call QcMatDst(Sp(i))
       
    end do

    
    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)
    call QcMatDst(U)

    deallocate(derivative_structure)
    deallocate(ind)
    deallocate(pert_ext)
    deallocate(Fp)
    deallocate(Dp)
    deallocate(Sp)
    deallocate(Dh)
    deallocate(blk_info)

  end subroutine
  
  
  ! Routine with support for multiple right-hand sides
    subroutine get_fds_2014_mc(pert, F, D, S, get_rsp_sol, get_ovl_mat, &
                       get_1el_mat, get_2el_mat, get_xc_mat, id_outp)

!    use interface_rsp_solver, only: rsp_solver_exec
    implicit none

    
    integer :: sstr_incr, i, j, superstructure_size, nblks, perturbed_matrix_size, id_outp
    integer :: ierr, npert_ext
    integer, allocatable, dimension(:) :: ind, blk_sizes, pert_ext
    integer, allocatable, dimension(:,:) :: blk_info, indices
    integer, dimension(0) :: noc
    character(4), dimension(0) :: nof
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(SDF_2014) :: F, D, S
    external :: get_rsp_sol, get_ovl_mat, get_1el_mat, get_2el_mat, get_xc_mat
    type(QcMat) :: A, B, C, zeromat, T, U
    type(QcMat), allocatable, dimension(:) :: Fp, Dp, Sp, Dh, X, RHS
    type(f_l_cache_2014), pointer :: fock_lowerorder_cache


    ! ASSUME CLOSED SHELL
    
    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)
    call QcMatInit(T)
    call QcMatInit(U)
    
    
    call sdf_getdata_s_2014(D, get_emptypert(), (/1/), A)
    call sdf_getdata_s_2014(S, get_emptypert(), (/1/), B)
    call sdf_getdata_s_2014(F, get_emptypert(), (/1/), C)
    
    

    nblks = get_num_blks(pert)

    allocate(blk_info(nblks, 3))
    allocate(blk_sizes(pert%npert))

    blk_info = get_blk_info(nblks, pert)
    perturbed_matrix_size = get_triangulated_size(nblks, blk_info)
    blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

    allocate(Fp(perturbed_matrix_size))
    allocate(Dp(perturbed_matrix_size))
    allocate(Sp(perturbed_matrix_size))
    allocate(Dh(perturbed_matrix_size))

    ! Process perturbation tuple for external call
    
    call p_tuple_to_external_tuple(pert, npert_ext, pert_ext)
    
    
    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix


    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Sp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Sp(i))

    end do
! write(*,*) 'Sp a', Sp(1)%elms

    call get_ovl_mat(0, noc, 0, noc, npert_ext, pert_ext, &
                     perturbed_matrix_size, Sp)

!     call rsp_ovlint(zeromat%nrow, pert%npert, pert%plab, &
!                        (/ (1, j = 1, pert%npert) /), pert%pdim, &
!                        nblks, blk_info, blk_sizes, &
!                        perturbed_matrix_size, Sp)

! write(*,*) 'Sp b', Sp(1)%elms

    call sdf_add_2014(S, pert, perturbed_matrix_size, Sp)

    deallocate(blk_sizes)

    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp

    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Dp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Dp(i), Sp(1))

!        call mat_init(Dh(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Dh(i), Sp(1))

!        call mat_init(Fp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call QcMatInit(Fp(i), Sp(1))

    end do

    call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

! write(*,*) 'Fp a', Fp(1)%elms

    call f_l_cache_allocate_2014(fock_lowerorder_cache)
    call rsp_fock_lowerorder_2014(pert, pert%npert, 1, (/get_emptypert()/), &
                         get_1el_mat, get_ovl_mat, get_2el_mat, get_xc_mat, 0, D, &
                         perturbed_matrix_size, Fp, fock_lowerorder_cache)

! write(*,*) 'Fp b', Fp(1)%elms

    deallocate(fock_lowerorder_cache)

    call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)

    ! b) For Dp: Create differentiation superstructure: First dryrun for size, and
    ! then the actual superstructure call

    superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%npert, pert%npert/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

    allocate(derivative_structure(superstructure_size, 3))
    allocate(indices(perturbed_matrix_size, pert%npert))
    allocate(ind(pert%npert))

    call derivative_superstructure(pert, (/pert%npert, &
         pert%npert/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)
    call make_triangulated_indices(nblks, blk_info, perturbed_matrix_size, indices)

    do i = 1, size(indices, 1)

       ind = indices(i, :)

       call rsp_get_matrix_z_2014(superstructure_size, derivative_structure, &
               (/pert%npert,pert%npert/), pert%npert, &
               (/ (j, j = 1, pert%npert) /), pert%npert, &
               ind, F, D, S, Dp(i))




!      Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A
       call QcMatkABC(-1.0d0, Dp(i), B, A, T)
       call QcMatkABC(-1.0d0, A, B, Dp(i), U)
       call QcMatRAXPY(1.0d0, T, U)
       call QcMatRAXPY(1.0d0, U, Dp(i))


    
       

       call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

       ! 3. Complete the particular contribution to Fp
       call cpu_time(time_start)
       call get_2el_mat(0, noc, noc, 1, (/Dp(i)/), 1, Fp(i:i))
!        call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i), &
!                           1, Fp(i:i))
       call cpu_time(time_end)
!        print *, 'seconds spent in 2-el particular contribution', time_end - time_start

       call cpu_time(time_start)
       ! MaR: Reintroduce once callback functionality is complete
!        call get_xc_mat(0, pert%pdim, noc, nof, 2, (/A, Dp(i)/), Fp(i:i))
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dp(i) /) , 1, Fp(i:i))
       call cpu_time(time_end)
!       print *, 'seconds spent in XC particular contribution', time_end - time_start

       ! MaR: Considered to be part of the 2el call
!        call cpu_time(time_start)
!        call rsp_pe(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i) , 1, Fp(i))
!        call cpu_time(time_end)
! !       print *, 'seconds spent in PE particular contribution', time_end - time_start

       call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)


       call QcMatInit(RHS(i))
       call QcMatInit(X(i))

       ! 4. Make right-hand side using Dp

       call rsp_get_matrix_y_2014(superstructure_size, derivative_structure, &
                pert%npert, (/ (j, j = 1, pert%npert) /), &
                pert%npert, ind, F, D, S, RHS(i))
       
    end do

    
    ! Note (MaR): Passing only real part of freq. Is this OK?
    call get_rsp_sol(1, (/sum(real(pert%freq(:)))/), size(indices, 1), RHS, X)

    
    do i = 1, size(indices, 1)   

       ! 5. Get Dh using the rsp equation solution X
       
!        Dh(i) = A*B*X(1) - X(1)*B*A
       call QcMatkABC(-1.0d0, X(i), B, A, T)
       call QcMatkABC(1.0d0, A, B, X(i), U)
       call QcMatRAXPY(1.0d0, T, U)
       call QcMatAEqB(Dh(i), U)

       ! 6. Make homogeneous contribution to Fock matrix

       call cpu_time(time_start)
       call get_2el_mat(0, noc, noc, 1, (/Dh(i)/), 1, Fp(i:i))

       call cpu_time(time_end)
!        print *, 'seconds spent in 2-el homogeneous contribution', time_end - time_start

       call cpu_time(time_start)
       ! MaR: Reintroduce once callback functionality is complete
!        call get_xc_mat(0, pert%pdim, noc, nof, 2, (/A, Dp(i)/), Fp(i:i))
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dh(i) /) , 1, Fp(i:i))
       call cpu_time(time_end)
!       print *, 'seconds spent in XC homogeneous contribution', time_end - time_start
       
       
!        call cpu_time(time_start)
!        call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), &
!                           1, Fp(i:i))
!        call cpu_time(time_end)
!        print *, 'seconds spent in 2-el homogeneous contribution', time_end - time_start

!        call cpu_time(time_start)
!        call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
!             (/ A, Dh(i) /) , 1, Fp(i:i))
!        call cpu_time(time_end)
!        print *, 'seconds spent in XC homogeneous contribution', time_end - time_start

!        call cpu_time(time_start)
!        call rsp_pe(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), 1, Fp(i))
!        call cpu_time(time_end)
!       print *, 'seconds spent in PE homogeneous contribution', time_end - time_start

       ! 7. Complete perturbed D with homogeneous part

!        Dp(i) = Dp(i) + Dh(i)
       call QcMatRAXPY(1.0d0, Dh(i), Dp(i))


       if (perturbed_matrix_size < 10) then

          write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'
  
       else

          if (mod(i, perturbed_matrix_size/10) == 1) then

             write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'

          end if

       end if

    end do

    ! Add the final values to cache

    call sdf_add_2014(F, pert, perturbed_matrix_size, Fp)
    call sdf_add_2014(D, pert, perturbed_matrix_size, Dp)

    do i = 1, size(indices, 1)

       call QcMatDst(Dh(i))
       call QcMatDst(Dp(i))
       call QcMatDst(Fp(i))
       call QcMatDst(Sp(i))
       ! Maybe no need to destroy X(i)
       call QcMatDst(X(i))
       call QcMatDst(RHS(i))
       
    end do

    
    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)
    call QcMatDst(U)

    deallocate(derivative_structure)
    deallocate(ind)
    deallocate(pert_ext)
    deallocate(Fp)
    deallocate(Dp)
    deallocate(Sp)
    deallocate(Dh)
    deallocate(blk_info)

  end subroutine
  
  

  
  
    recursive subroutine rsp_fock_lowerorder_2014(pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
                       density_order, D, property_size, Fp, fock_lowerorder_cache)

    implicit none

    logical :: density_order_skip
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF_2014) :: D
    external :: get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat
    type(QcMat), dimension(property_size) :: Fp
    type(f_l_cache_2014) :: fock_lowerorder_cache

    if (pert%npert >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%npert == 0) then

          call rsp_fock_lowerorder_2014(p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
               density_order, D, property_size, Fp, fock_lowerorder_cache)

       else

          call rsp_fock_lowerorder_2014(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
               p_tuples(2:size(p_tuples))/), &
               get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
               density_order, D, property_size, Fp, fock_lowerorder_cache)

       end if
    
       ! 2. Differentiate all of the contraction densities in turn

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%npert == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call rsp_fock_lowerorder_2014(p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               t_new, get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
               density_order + 1, D, property_size, Fp, fock_lowerorder_cache)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_fock_lowerorder_2014(p_tuple_remove_first(pert), &
            total_num_perturbations, num_p_tuples + 1, &
            (/p_tuples(:), p_tuple_getone(pert, 1)/), &
            get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
            density_order + 1, D, property_size, Fp, fock_lowerorder_cache)

    else

!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)

       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%npert >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then

          if (f_l_cache_already_2014(fock_lowerorder_cache, &
          num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples)) .EQV. .FALSE.) then

       write(*,*) 'Calculating perturbed Fock matrix lower order contribution'

       do i = 1, num_p_tuples
 
          if (i == 1) then

             write(*,*) 'F', p_tuples(i)%pid

          else

             write(*,*) 'D', p_tuples(i)%pid

          end if

       end do

             call get_fock_lowerorder_2014(num_p_tuples, total_num_perturbations, &
                                      p_tuples_standardorder(num_p_tuples, p_tuples), &
                                      density_order, get_1el_mat, get_t_mat, get_2el_mat, get_xc_mat, &
                                      D, property_size, Fp, fock_lowerorder_cache)

             write(*,*) 'Calculated perturbed Fock matrix lower order contribution'
             write(*,*) ' '

          else

             call f_l_cache_getdata_2014(fock_lowerorder_cache, num_p_tuples, &
                                    p_tuples_standardorder(num_p_tuples, p_tuples), &
                                    property_size, Fp)

!              write(*,*) ' '

          end if

       else

!           write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
!           write(*,*) 'at order for which perturbed D is to be found '
!           write(*,*) ' '

       end if

    end if

  end subroutine




  subroutine get_fock_lowerorder_2014(num_p_tuples, total_num_perturbations, p_tuples, &
                                 density_order, get_1el_mat, get_ovl_mat, get_2el_mat, get_xc_mat, &
                                 D, property_size, Fp, fock_lowerorder_cache)

    implicit none
    
    type(p_tuple) :: merged_p_tuple, t_matrix_bra, t_matrix_ket, t_matrix_newpid
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF_2014) :: D
    type(QcMat), allocatable, dimension(:) :: dens_tuple
    integer :: i, j, k, m, num_p_tuples, total_num_perturbations, merged_nblks, &
               density_order, property_size, fp_offset, lo_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, offset, npert_ext
    integer, dimension(0) :: noc
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter, &
                                                pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: ncoutersmall, pidoutersmall, ncinnersmall
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: blk_sizes_merged, pert_ext
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    integer, allocatable, dimension(:,:) :: triang_indices_fp, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    external :: get_1el_mat, get_ovl_mat, get_2el_mat, get_xc_mat
    type(QcMat) :: zeromat, D_unp
    type(QcMat), allocatable, dimension(:) :: tmp, lower_order_contribution
    type(QcMat), dimension(property_size) :: Fp
    type(f_l_cache_2014) :: fock_lowerorder_cache

!    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
!    ncouter = nc_only(total_num_perturbations, total_num_perturbations - & 
!                      p_tuples(1)%npert, num_p_tuples - 1, &
!                      p_tuples(2:num_p_tuples), ncarray)
!    ncinner = nc_only(total_num_perturbations, p_tuples(1)%npert, 1, &
!                      p_tuples(1), ncarray)

    allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%npert))
    allocate(ncinnersmall(p_tuples(1)%npert))
    allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%npert))

!    ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
!                                p_tuples(1)%npert, num_p_tuples - 1, &
!                                p_tuples(2:num_p_tuples), ncarray)
!    ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%npert, &
!                   1, p_tuples(1), ncarray)
!    pidoutersmall = get_pidoutersmall(total_num_perturbations - &
!                    p_tuples(1)%npert, num_p_tuples - 1, &
 !                   p_tuples(2:num_p_tuples))

    ! MaR: Second way of blks_tuple_info can in the general case be larger than
    ! needed, but is allocated this way to get a prismic data structure
    allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(num_p_tuples))
    allocate(blk_sizes(num_p_tuples, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))
    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    !FIXME Gao: we do not need dens_tuple(1)?
    allocate(dens_tuple(2:num_p_tuples))
    allocate(nfields(num_p_tuples))
    allocate(nblks_tuple(num_p_tuples))

    
    call p_tuple_to_external_tuple(p_tuples(1), npert_ext, pert_ext)
    
    
    call p1_cloneto_p2(p_tuples(1), t_matrix_newpid)
    t_matrix_newpid%pid = (/(i, i = 1, t_matrix_newpid%npert)/)


    do i = 1, num_p_tuples

       nfields(i) = p_tuples(i)%npert
       nblks_tuple(i) = get_num_blks(p_tuples(i))

    end do

    do i = 1, num_p_tuples

       call get_blk_info_s(nblks_tuple(i), p_tuples(i), blks_tuple_info(i, 1:nblks_tuple(i), :))

! write(*,*) blks_tuple_info(i, :, :)
! write(*,*) 'sanitized'
! write(*,*) blks_tuple_info(i, 1:nblks_tuple(i), :)

       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))


! write(*,*)  blks_tuple_triang_size(i) 

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
    allocate(lower_order_contribution(inner_indices_size * outer_indices_size))

    o_whichpert = make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)
!    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
!                      p_tuples(1)%npert, pidoutersmall, &
!                      ncarray, ncoutersmall, o_whichpert)

    call sdf_getdata_s_2014(D, get_emptypert(), (/1/), D_unp)

    !FIXME Gao: we should allocate and initialize lower_order_contribution and
    !           dens_tuple in the if statement where they are used, also their
    !           deallocation should be moved
    do j = 1, size(lower_order_contribution)
       call QcMatInit(lower_order_contribution(j), Fp(1))
    end do

    do j = 1, size(tmp)
       call QcMatInit(tmp(j), Fp(1))
    end do

    do i = 2, num_p_tuples
       call QcMatInit(dens_tuple(j), Fp(1))
    end do


    if (total_num_perturbations > p_tuples(1)%npert) then

       k = 1

       do i = 2, num_p_tuples
          do j = 1, p_tuples(i)%npert

             o_wh_forave(p_tuples(i)%pid(j)) = k
             k = k + 1

          end do
       end do

       allocate(outer_indices(outer_indices_size,total_num_perturbations - &
                p_tuples(1)%npert))
       allocate(inner_indices(inner_indices_size,p_tuples(1)%npert))

       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)

       if (p_tuples(1)%npert > 0) then

          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)

       end if

       do i = 1, size(outer_indices, 1)

          do j = 1, size(lower_order_contribution)
             call QcMatZero(lower_order_contribution(j))
          end do
     
          do j = 1, size(tmp)
             call QcMatZero(tmp(j))
          end do


          do j = 2, num_p_tuples

             call sdf_getdata_s_2014(D, p_tuples(j), (/ &
                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
                             k = 1, p_tuples(j)%npert) /), dens_tuple(j))

          end do

          if (num_p_tuples <= 2) then

             call cpu_time(time_start)
             
             call get_2el_mat(npert_ext, pert_ext, 1, (/dens_tuple(2)/), size(tmp), tmp)
             
!              call rsp_twoint(zeromat%nrow, p_tuples(1)%npert, p_tuples(1)%plab, &
!                              (/ (1, j = 1, p_tuples(1)%npert) /), &
!                              p_tuples(1)%pdim, dens_tuple(2), size(tmp), tmp)
             call cpu_time(time_end)
!              print *, 'seconds spent in 2-el contribution', time_end - time_start
          end if

          ! MaR: Reintroduce after minimal working version is complete
          call cpu_time(time_start)
          
!           call get_xc_mat(p_tuples(1)%npert, p_tuples(1)%pdim, &
!                           (/ (1, j = 1, p_tuples(1)%npert) /), p_tuples(1)%plab, &
!                           num_p_tuples, (/ D_unp, (dens_tuple(k), k = 2, num_p_tuples) /), tmp)
          
!           call rsp_xcint_adapt(zeromat%nrow, p_tuples(1)%npert, &
!                p_tuples(1)%plab, (/ (1, j = 1, p_tuples(1)%npert) /), &
!                p_tuples(1)%pdim, (/ D_unp, &
!                (dens_tuple(k), k = 2, num_p_tuples) /), property_size, tmp)
          call cpu_time(time_end)
!           print *, 'seconds spent in XC contribution', time_end - time_start

          ! MaR: Remove and consider as part of 2el contribution
!           if (num_p_tuples <= 2) then
!              call cpu_time(time_start)
!              call rsp_pe(zeromat%nrow, p_tuples(1)%npert, p_tuples(1)%plab, &
!                              (/ (1, j = 1, p_tuples(1)%npert) /), &
!                              p_tuples(1)%pdim, dens_tuple(2), size(tmp), tmp)
!              call cpu_time(time_end)
! !             print *, 'seconds spent in PE contribution', time_end - time_start
!           end if

          if (p_tuples(1)%npert > 0) then

             do j = 1, size(inner_indices,1)

                offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, &
                nblks_tuple, (/ (p_tuples(k)%npert, k = 1, num_p_tuples) /), &
                blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                (/inner_indices(j, :), outer_indices(i, :) /)) 

                call QcMatAEqB(lower_order_contribution(offset),tmp(j))

             end do

          else

             ! MaR: There might be problems with this call (since the first p_tuple is empty)

             offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations, &
             nblks_tuple(2:num_p_tuples), &
             (/ (p_tuples(k)%npert, k = 2, num_p_tuples) /), &
             blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
             blks_tuple_triang_size(2:num_p_tuples), (/outer_indices(i, :) /)) 

             call QcMatAEqB(lower_order_contribution(offset),tmp(1))

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

             ! MaR: This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       end if

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

       k = 1
       do i = 1, num_p_tuples
          do j = 1, p_tuples(i)%npert
             pids_current_contribution(k) = p_tuples(i)%pid(j)
             k = k + 1
          end do
       end do

! write(*,*) 'merged plab', merged_p_tuple%plab

       merged_nblks = get_num_blks(merged_p_tuple)

! write(*,*) 'merged plab 2', merged_p_tuple%plab

       allocate(merged_blk_info(1, merged_nblks, 3))

! write(*,*) 'allocate OK', merged_p_tuple%plab

       call get_blk_info_s(merged_nblks, merged_p_tuple, merged_blk_info(1, :, :))

! write(*,*) 'merged plab 3', merged_p_tuple%plab
! 
! do i = 1, merged_nblks
! 
! write(*,*) 'merged block info', merged_blk_info(1,i,:)
! 
! end do

       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(triang_indices_fp(merged_triang_size, sum(merged_blk_info(1, :,2))))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, triang_indices_fp)

! do i = 1, size(triang_indices_fp,1)
! 
! write(*,*) 'triang indices', triang_indices_fp(i,:)
! 
! end do
! 
! write(*,*) 'size Fp', size(Fp)
! write(*,*) 'size loc', size(lower_order_contribution)


       do i = 1, size(triang_indices_fp, 1)

          fp_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                      (/sum(nfields)/), &
                      (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                      (/triang_indices_fp(i, :) /))

          do j = 1, total_num_perturbations
    
             translated_index(j) = triang_indices_fp(i,pids_current_contribution(j))
    
          end do

          if (p_tuples(1)%npert > 0) then

             lo_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                         total_num_perturbations, nblks_tuple, &
                         nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/translated_index(:)/))

          else

             lo_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
                         total_num_perturbations, nblks_tuple(2:num_p_tuples), &
                         nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
                         blk_sizes(2:num_p_tuples,:), &
                         blks_tuple_triang_size(2:num_p_tuples), & 
                         (/translated_index(:)/))

          end if

          call QcMatRAXPY(1.0d0, lower_order_contribution(lo_offset), Fp(fp_offset))

       end do

       call f_l_cache_add_element_2014(fock_lowerorder_cache, num_p_tuples, p_tuples, &
            inner_indices_size * outer_indices_size, lower_order_contribution)

       deallocate(merged_blk_info)
       deallocate(triang_indices_fp)
       deallocate(outer_indices)
       deallocate(inner_indices)

    else


       if (num_p_tuples <= 1) then

          call get_1el_mat(npert_ext, pert_ext, size(Fp), Fp)
       
!           call rsp_oneint(zeromat%nrow, p_tuples(1)%npert, p_tuples(1)%plab, &
!                           (/ (1, j = 1, p_tuples(1)%npert) /), &
!                           p_tuples(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
!                    1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), property_size, Fp)

! NOTE: Find out if necessary ovlint/oneint in "outer indices case" above
! NOTE (Oct 12): Probably not unless some hidden density matrix dependence


          t_matrix_bra = get_emptypert()
          t_matrix_ket = get_emptypert()
          call rsp_ovlint_t_matrix_2014(t_matrix_newpid%npert, t_matrix_newpid, &
                                   t_matrix_bra, t_matrix_ket, get_ovl_mat, property_size, Fp)
       end if

       if (num_p_tuples <= 2) then


          call cpu_time(time_start)
          
          call get_2el_mat(npert_ext, pert_ext, 1, (/D_unp/), size(Fp), Fp)
          
!           call rsp_twoint(zeromat%nrow, p_tuples(1)%npert, p_tuples(1)%plab, &
!                (/ (1, j = 1, p_tuples(1)%npert) /), &
!                p_tuples(1)%pdim, D_unp, &
!                property_size, Fp)
          call cpu_time(time_end)
!           print *, 'seconds spent in 2-el contribution', time_end - time_start

       end if

       ! MaR: Reintroduce when minimal working version is complete
!        call cpu_time(time_start)
!        call rsp_xcint_adapt(zeromat%nrow, p_tuples(1)%npert, p_tuples(1)%plab, &
!                       (/ (1, j = 1, p_tuples(1)%npert) /), &
!                       p_tuples(1)%pdim, &
!                       (/ D_unp /), &
!                       property_size, Fp)
!        call cpu_time(time_end)
! !        print *, 'seconds spent in XC contribution', time_end - time_start

       ! MaR: Remove and consider as part of 2el contribution
!        if (num_p_tuples <= 2) then
!           call cpu_time(time_start)
!           call rsp_pe(zeromat%nrow,                                 &
!                       p_tuples(1)%npert,                  &
!                       p_tuples(1)%plab,                             &
!                       (/ (1, j = 1, p_tuples(1)%npert) /),&
!                       p_tuples(1)%pdim,                             &
!                       D_unp,                                        &
!                       property_size,                                &
!                       Fp)
!           call cpu_time(time_end)
! !          print *, 'seconds spent in PE contribution', time_end - time_start
!        end if

       ! MaR: THERE IS NO NEED TO CACHE THE "ALL INNER" CONTRIBUTION
       ! It should be possible to just add it to Fp like already done above
       ! even with the extra complexity from the triangularization 

    end if

    call QcMatDst(D_unp)

    do i = 2, num_p_tuples
   
       call QcMatDst(dens_tuple(i))
   
    end do

    do i = 1, size(tmp)

       call QcMatDst(tmp(i))

    end do

    do i = 1, size(lower_order_contribution)

       call QcMatDst(lower_order_contribution(i))

    end do

    deallocate(dens_tuple)


    deallocate(pert_ext)
    
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(tmp)
    deallocate(lower_order_contribution)

    ! MaR: Why is the next line commented? Find out
!     deallocate(dens_tuple)

  end subroutine
  
  
  end module
