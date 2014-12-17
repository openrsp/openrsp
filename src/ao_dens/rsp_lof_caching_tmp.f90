! Copyright 2012      Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines/functions and definitions related to fock_lowerorder datatype
! in which (perturbed) overlap, density and Fock matrices are stored.

module rsp_lof_caching

  use rsp_field_tuple, only: p_tuple,                &
                             p_tuple_standardorder,  &
                             p_tuples_compare,       &
                             p_tuples_standardorder, &
                             merge_p_tuple,          &
                             p_tuple_p1_cloneto_p2
  use rsp_indices_and_addressing
  use matrix_defop, matrix => openrsp_matrix
  use matrix_lowlevel, only: mat_init
  use qmatrix

  implicit none

  public f_l_cache_init
  public f_l_cache_getdata
  public f_l_cache_next_element
  public f_l_cache_add_element
  public f_l_cache_already
  public f_l_cache_allocate

  public f_l_cache_init_2014
  public f_l_cache_getdata_2014
  public f_l_cache_next_element_2014
  public f_l_cache_add_element_2014
  public f_l_cache_already_2014
  public f_l_cache_allocate_2014
  
  ! Define Fock lower order contribution cache datatype

  type f_l_cache

     type(f_l_cache), pointer :: next
     logical :: last
     integer :: num_p_tuples
     ! MaR: Should all of the data attributes be pointers too?
     type(p_tuple), allocatable, dimension(:) :: p_tuples
     type(matrix), allocatable, dimension(:) :: data ! (Perturbed) matrix data

  end type

  type f_l_cache_2014

     type(f_l_cache_2014), pointer :: next
     logical :: last
     integer :: num_p_tuples
     ! MaR: Should all of the data attributes be pointers too?
     type(p_tuple), allocatable, dimension(:) :: p_tuples
     type(qmat), allocatable, dimension(:) :: data ! (Perturbed) matrix data

  end type
  
  contains

  ! Begin SDF linked list manipulation/data retrieval routines

  subroutine f_l_cache_init(new_element, num_p_tuples, p_tuples, property_size, data)

    implicit none

    integer :: i, num_p_tuples, property_size
    type(f_l_cache) :: new_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(matrix), dimension(property_size) :: data

    new_element%last = .TRUE.
    new_element%num_p_tuples = num_p_tuples
    allocate(new_element%p_tuples(num_p_tuples))

    do i = 1, num_p_tuples
       call p_tuple_p1_cloneto_p2(p_tuples(i), new_element%p_tuples(i))
    end do

    allocate(new_element%data(property_size))

    do i = 1, property_size
       ! ASSUME CLOSED SHELL
       !call mat_init(new_element%data(i), data(i)%nrow, data(i)%ncol)
       !call mat_init_like_and_zero(data(i), new_element%data(i))

       new_element%data(i) = data(i)

    end do
 
  end subroutine

  subroutine f_l_cache_getdata(cache, num_p_tuples, p_tuples, property_size, Fp)

    implicit none

    logical :: found
    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
               property_size, total_num_perturbations, fp_offset, cache_offset, &
               merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contribution, & 
                                          p_tuples_dimensions, &
                                          p_tuples_dimensions_cacheorder, &
                                          pids_merged_pert, translated_index
    integer, allocatable, dimension(:) :: blk_sizes_merged
    integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:,:) :: indices, blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
    type(f_l_cache), target :: cache
    type(f_l_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    type(p_tuple) :: merged_p_tuple
    type(matrix), dimension(property_size) :: Fp

    next_element => cache
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    found = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

       next_element => f_l_cache_next_element(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          found = p_tuples_compare(num_p_tuples, p_tuples_standardorder( &
                                   next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples_st_order)

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (found .eqv. .TRUE.) then

!        write(*,*) 'Getting f_l_cache data'

       total_num_perturbations = 0

       do i = 1, num_p_tuples
          total_num_perturbations = total_num_perturbations + p_tuples(i)%n_perturbations
       end do

       next_element%p_tuples = p_tuples_standardorder(next_element%num_p_tuples, &
                                                      next_element%p_tuples)

       ! get triangulated indices for Fp

       if (p_tuples(1)%n_perturbations > 0) then

          call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

          do i = 2, num_p_tuples

             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       else

          call p_tuple_p1_cloneto_p2(p_tuples(2), merged_p_tuple)

          do i = 3, num_p_tuples

             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       end if

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
       merged_nblks = get_num_blks(merged_p_tuple)

       allocate(merged_blk_info(1,merged_nblks, 3))

       merged_blk_info(1,:,:) = get_blk_info(merged_nblks, merged_p_tuple)
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(blk_sizes(num_p_tuples, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))

       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

       allocate(indices(property_size, total_num_perturbations))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, indices)


       allocate(translated_index(total_num_perturbations))
       allocate(pids_in_cache(total_num_perturbations))
       allocate(pids_current_contribution(total_num_perturbations))
       allocate(pids_merged_pert(total_num_perturbations))
       allocate(p_tuples_dimensions(total_num_perturbations))
       allocate(p_tuples_dimensions_cacheorder(total_num_perturbations))

       k = 1

       do i = 1, num_p_tuples
          do j = 1, p_tuples_st_order(i)%n_perturbations
             pids_current_contribution(k) = p_tuples_st_order(i)%pid(j)
          k = k + 1
          end do
       end do

       do i = 1, total_num_perturbations

          pids_merged_pert(i) = merged_p_tuple%pid(i)

       end do

       p_tuples_dimensions = get_ncarray(total_num_perturbations, num_p_tuples, &
                                         p_tuples_st_order)

       allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))

       do i = 1, num_p_tuples

          nfields(i) = p_tuples_st_order(i)%n_perturbations
          nblks_tuple(i) = get_num_blks(p_tuples_st_order(i))
          blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples_st_order(i))
          blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
          blks_tuple_info(i,1:nblks_tuple(i),:))
          blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
          blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

       end do

       do i = 1, size(indices, 1)

          ! To which element in the cached data does that 
          ! cache p_tuples index tuple correspond?
          ! Get that element
          ! To which element in the property tensor does that correspond? 
          ! Put the element there

          fp_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
          (/merged_nblks/), & 
          (/total_num_perturbations/), (/merged_blk_info/), blk_sizes_merged, &
          (/merged_triang_size/), &
          (/indices(i, : )/) )

          do j = 1, total_num_perturbations

             translated_index(j) = indices(i,pids_current_contribution(j))

          end do

          if (p_tuples(1)%n_perturbations > 0) then

             cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
             total_num_perturbations, nblks_tuple, & 
             nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)

          else

             cache_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
             total_num_perturbations, nblks_tuple(2:num_p_tuples), & 
             nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :),&
             blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
             translated_index)

          end if

          Fp(fp_offset) = &
          Fp(fp_offset) + &
          next_element%data(cache_offset) 

       end do

!        write(*,*) 'Got Fock lower order contribution cache data'

       deallocate(indices)
       deallocate(blks_tuple_info)
       deallocate(merged_blk_info)

    else

       write(*,*) 'Failed to retrieve data in property_cache_getdata: Element not found'

    end if

    deallocate(translated_index)
    deallocate(pids_in_cache)
    deallocate(pids_current_contribution)
    deallocate(pids_merged_pert)
    deallocate(p_tuples_dimensions)
    deallocate(p_tuples_dimensions_cacheorder)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)

  end subroutine

  ! BEGIN NEW CODE


  subroutine f_l_cache_init_2014(new_element, num_p_tuples, p_tuples, property_size, data)

    implicit none

    integer :: i, num_p_tuples, property_size
    type(f_l_cache_2014) :: new_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(qmat), dimension(property_size) :: data

    new_element%last = .TRUE.
    new_element%num_p_tuples = num_p_tuples
    allocate(new_element%p_tuples(num_p_tuples))

    do i = 1, num_p_tuples
       call p_tuple_p1_cloneto_p2(p_tuples(i), new_element%p_tuples(i))
    end do

    allocate(new_element%data(property_size))

    do i = 1, property_size

       call QMatInit(new_element%data(i))
       call QMatAEqB(new_element%data(i), data(i))

    end do
 
  end subroutine

  subroutine f_l_cache_getdata_2014(cache, num_p_tuples, p_tuples, property_size, Fp)

    implicit none

    logical :: found
    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
               property_size, total_num_perturbations, fp_offset, cache_offset, &
               merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contribution, & 
                                          p_tuples_dimensions, &
                                          p_tuples_dimensions_cacheorder, &
                                          pids_merged_pert, translated_index
    integer, allocatable, dimension(:) :: blk_sizes_merged
    integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:,:) :: indices, blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
    type(f_l_cache_2014), target :: cache
    type(f_l_cache_2014), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    type(p_tuple) :: merged_p_tuple
    type(qmat), dimension(property_size) :: Fp

    next_element => cache
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    found = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

       next_element => f_l_cache_next_element_2014(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          found = p_tuples_compare(num_p_tuples, p_tuples_standardorder( &
                                   next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples_st_order)

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (found .eqv. .TRUE.) then

!        write(*,*) 'Getting f_l_cache data'

       total_num_perturbations = 0

       do i = 1, num_p_tuples
          total_num_perturbations = total_num_perturbations + p_tuples(i)%n_perturbations
       end do

       next_element%p_tuples = p_tuples_standardorder(next_element%num_p_tuples, &
                                                      next_element%p_tuples)

       ! get triangulated indices for Fp

       if (p_tuples(1)%n_perturbations > 0) then

          call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

          do i = 2, num_p_tuples

             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       else

          call p_tuple_p1_cloneto_p2(p_tuples(2), merged_p_tuple)

          do i = 3, num_p_tuples

             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       end if

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
       merged_nblks = get_num_blks(merged_p_tuple)

       allocate(merged_blk_info(1,merged_nblks, 3))

       merged_blk_info(1,:,:) = get_blk_info(merged_nblks, merged_p_tuple)
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(blk_sizes(num_p_tuples, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))

       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

       allocate(indices(property_size, total_num_perturbations))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, indices)


       allocate(translated_index(total_num_perturbations))
       allocate(pids_in_cache(total_num_perturbations))
       allocate(pids_current_contribution(total_num_perturbations))
       allocate(pids_merged_pert(total_num_perturbations))
       allocate(p_tuples_dimensions(total_num_perturbations))
       allocate(p_tuples_dimensions_cacheorder(total_num_perturbations))

       k = 1

       do i = 1, num_p_tuples
          do j = 1, p_tuples_st_order(i)%n_perturbations
             pids_current_contribution(k) = p_tuples_st_order(i)%pid(j)
          k = k + 1
          end do
       end do

       do i = 1, total_num_perturbations

          pids_merged_pert(i) = merged_p_tuple%pid(i)

       end do

       p_tuples_dimensions = get_ncarray(total_num_perturbations, num_p_tuples, &
                                         p_tuples_st_order)

       allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))

       do i = 1, num_p_tuples

          nfields(i) = p_tuples_st_order(i)%n_perturbations
          nblks_tuple(i) = get_num_blks(p_tuples_st_order(i))
          blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples_st_order(i))
          blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
          blks_tuple_info(i,1:nblks_tuple(i),:))
          blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
          blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

       end do

       do i = 1, size(indices, 1)

          ! To which element in the cached data does that 
          ! cache p_tuples index tuple correspond?
          ! Get that element
          ! To which element in the property tensor does that correspond? 
          ! Put the element there

          fp_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
          (/merged_nblks/), & 
          (/total_num_perturbations/), (/merged_blk_info/), blk_sizes_merged, &
          (/merged_triang_size/), &
          (/indices(i, : )/) )

          do j = 1, total_num_perturbations

             translated_index(j) = indices(i,pids_current_contribution(j))

          end do

          if (p_tuples(1)%n_perturbations > 0) then

             cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
             total_num_perturbations, nblks_tuple, & 
             nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)

          else

             cache_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
             total_num_perturbations, nblks_tuple(2:num_p_tuples), & 
             nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :),&
             blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
             translated_index)

          end if

          call QMatRAXPY(1.0d0, next_element%data(cache_offset), Fp(fp_offset))           

       end do

!        write(*,*) 'Got Fock lower order contribution cache data'

       deallocate(indices)
       deallocate(blks_tuple_info)
       deallocate(merged_blk_info)

    else

       write(*,*) 'Failed to retrieve data in property_cache_getdata: Element not found'

    end if

    deallocate(translated_index)
    deallocate(pids_in_cache)
    deallocate(pids_current_contribution)
    deallocate(pids_merged_pert)
    deallocate(p_tuples_dimensions)
    deallocate(p_tuples_dimensions_cacheorder)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)

  end subroutine
  
  
  ! Add element routine
  ! NOTE(MaR): This routine assumes that the pert_tuple and data
  ! is already in standard order

  subroutine f_l_cache_add_element_2014(current_element, num_p_tuples, p_tuples, &
                                   property_size, data)

    implicit none

    integer :: num_p_tuples, property_size !, i
    type(f_l_cache_2014), target :: current_element
    type(f_l_cache_2014), pointer :: new_element
    type(f_l_cache_2014), pointer :: new_element_ptr
    type(f_l_cache_2014), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(qmat), dimension(property_size) :: data

    next_element => current_element

    allocate(new_element)

    call f_l_cache_init_2014(new_element, num_p_tuples, p_tuples, &
                              property_size, data)

    new_element_ptr => new_element

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while (next_element%last .eqv. .FALSE.)
       next_element => f_l_cache_next_element_2014(next_element)
    end do

    next_element%last = .FALSE.
    new_element%next => next_element%next
    next_element%next => new_element

  end subroutine

  function f_l_cache_next_element_2014(current_element) result(next_element)

    implicit none

    type(f_l_cache_2014), target :: current_element
    type(f_l_cache_2014), pointer :: next_element

    next_element => current_element%next

  end function
  
    function f_l_cache_already_2014(current_element, num_p_tuples, p_tuples)

    implicit none

    logical :: f_l_cache_already_2014
    integer :: passedlast, num_p_tuples
    type(f_l_cache_2014), target :: current_element
    type(f_l_cache_2014), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    
    next_element => current_element
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    f_l_cache_already_2014 = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (f_l_cache_already_2014 .eqv. .FALSE.))

       next_element => f_l_cache_next_element_2014(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          f_l_cache_already_2014 = p_tuples_compare(num_p_tuples, &
                                   p_tuples_standardorder(next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples_st_order)

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (f_l_cache_already_2014 .EQV. .TRUE.) then

!         write(*,*) 'f_l_cache_already: Found element in cache'

    else

!         write(*,*) 'f_l_cache_already: Element not in cache'

    end if

  end function
  
  
    subroutine f_l_cache_allocate_2014(current_element)

    implicit none

    type(f_l_cache_2014), pointer :: current_element

    allocate(current_element)

    current_element%next => current_element
    current_element%num_p_tuples = 1
    current_element%last = .TRUE.

    allocate(current_element%p_tuples(1))

    current_element%p_tuples(1)%n_perturbations = 0

    allocate(current_element%p_tuples(1)%pdim(1))
    allocate(current_element%p_tuples(1)%plab(1))
    allocate(current_element%p_tuples(1)%pid(1))
    allocate(current_element%p_tuples(1)%freq(1))
    allocate(current_element%data(1))

  end subroutine
  
  
  
  ! END NEW CODE
  
  
  function f_l_cache_next_element(current_element) result(next_element)

    implicit none

    type(f_l_cache), target :: current_element
    type(f_l_cache), pointer :: next_element

    next_element => current_element%next

  end function

  ! Add element routine
  ! NOTE(MaR): This routine assumes that the pert_tuple and data
  ! is already in standard order

  subroutine f_l_cache_add_element(current_element, num_p_tuples, p_tuples, &
                                   property_size, data)

    implicit none

    integer :: num_p_tuples, property_size !, i
    type(f_l_cache), target :: current_element
    type(f_l_cache), pointer :: new_element
    type(f_l_cache), pointer :: new_element_ptr
    type(f_l_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(matrix), dimension(property_size) :: data

    next_element => current_element

    allocate(new_element)

    call f_l_cache_init(new_element, num_p_tuples, p_tuples, &
                              property_size, data)

    new_element_ptr => new_element

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while (next_element%last .eqv. .FALSE.)
       next_element => f_l_cache_next_element(next_element)
    end do

    next_element%last = .FALSE.
    new_element%next => next_element%next
    next_element%next => new_element

  end subroutine

  function f_l_cache_already(current_element, num_p_tuples, p_tuples)

    implicit none

    logical :: f_l_cache_already
    integer :: passedlast, num_p_tuples
    type(f_l_cache), target :: current_element
    type(f_l_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    
    next_element => current_element
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    f_l_cache_already = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (f_l_cache_already .eqv. .FALSE.))

       next_element => f_l_cache_next_element(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          f_l_cache_already = p_tuples_compare(num_p_tuples, &
                                   p_tuples_standardorder(next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples_st_order)

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (f_l_cache_already .EQV. .TRUE.) then

!         write(*,*) 'f_l_cache_already: Found element in cache'

    else

!         write(*,*) 'f_l_cache_already: Element not in cache'

    end if

  end function


  subroutine f_l_cache_allocate(current_element)

    implicit none

    type(f_l_cache), pointer :: current_element

    allocate(current_element)

    current_element%next => current_element
    current_element%num_p_tuples = 1
    current_element%last = .TRUE.

    allocate(current_element%p_tuples(1))

    current_element%p_tuples(1)%n_perturbations = 0

    allocate(current_element%p_tuples(1)%pdim(1))
    allocate(current_element%p_tuples(1)%plab(1))
    allocate(current_element%p_tuples(1)%pid(1))
    allocate(current_element%p_tuples(1)%freq(1))
    allocate(current_element%data(1))

  end subroutine

end module
