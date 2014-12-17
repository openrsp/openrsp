! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and functions related to caching of contributions to the
! response property tensor.

module rsp_property_caching

  use rsp_field_tuple
  use rsp_indices_and_addressing
!  use matrix_backend

  implicit none

  public property_cache_initialize
  public property_cache_next_element
  public property_cache_add_element
  public property_cache_already
  public property_cache_getdata
  public property_cache_allocate

  ! NEW 2014
  
!  public contrib_cache_initialize
!  public contrib_cache_next_element
!  public contrib_cache_outer_next_element
!  public contrib_cache_add_element
!  public contrib_cache_already
!  public contrib_cache_getdata
!  public contrib_cache_allocate
  
  ! END NEW 2014
  
  
  !  Define property contribution cache datatype

 type property_cache

     type(property_cache), pointer :: next
     logical :: last
     integer :: num_p_tuples
     type(p_tuple), allocatable, dimension(:) :: p_tuples
     complex(8), allocatable, dimension(:) :: data ! Property data    

  end type
  
  ! NEW 2014
  
  ! Define contrib cache datatype

!  type contrib_cache_outer
!
!    type(contrib_cache_outer), pointer :: next
!    logical :: last
!    integer :: num_dmat
!    type(p_tuple), allocatable, dimension(:) :: outer_p_tuples
!    integer :: outer_size
!    integer :: contrib_size
!    integer, allocatable, dimension(:) :: nblks_tuple
!    integer, allocatable, dimension(:,:) :: blk_sizes
!    integer, allocatable, dimension(:,:) :: indices
!    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
!    integer, allocatable, dimension(:) :: blks_tuple_triang_size
!    type(openrsp_matrix), allocatable, dimension(:) :: data_int ! Fock matrix contribution data
!    complex(8), allocatable, dimension(:) :: data_ave ! Property data    
!
!  end type 
!  
!  type contrib_cache
!
!    type(contrib_cache), pointer :: next
!    logical :: last
!    type(p_tuple) :: p_inner
!    integer :: num_outer
!    type(contrib_cache_outer) :: contribs_outer
!
!  end type 
  
  
  ! END NEW 2014
  
  
  

  contains

 
  
  
  recursive subroutine printt_rsp_tensor_stdout(npert, lvl, pdim, prop, offset)

    implicit none

    integer :: npert, i, j, offset, lvl, new_offset
    integer, dimension(npert) :: pdim
    complex(8), dimension(product(pdim)) :: prop

    if (lvl > 1) then

    do i = 1, pdim(npert - lvl + 1)

       new_offset = offset + (i - 1)*product(pdim(npert - lvl + 1:npert))/ &
                                             pdim(npert - lvl + 1)

       call printt_rsp_tensor_stdout(npert, lvl - 1, pdim, prop, new_offset)

    end do

    write(*,*) ' '

    else

    write(*,*) real(prop(offset:offset+pdim(npert) - 1))

    end if

  end subroutine


  subroutine property_cache_initialize(new_element, num_p_tuples, p_tuples, &
                                       property_size, data)

    implicit none

    integer :: i, num_p_tuples, property_size
    type(property_cache) :: new_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    complex(8), dimension(property_size) :: data


    new_element%last = .TRUE.
    new_element%num_p_tuples = num_p_tuples
    allocate(new_element%p_tuples(num_p_tuples))

    do i = 1, num_p_tuples
       call p_tuple_p1_cloneto_p2(p_tuples(i), new_element%p_tuples(i))
    end do

    allocate(new_element%data(property_size))
    new_element%data = data

  end subroutine


  function property_cache_next_element(current_element) result(next_element)

    implicit none

    type(property_cache), target :: current_element
    type(property_cache), pointer :: next_element

    next_element => current_element%next

  end function


  subroutine property_cache_add_element(current_element, num_p_tuples, p_tuples, &
                                        property_size, data)

    implicit none

    integer :: num_p_tuples, property_size !, i
    type(property_cache), target :: current_element
    type(property_cache), pointer :: new_element
    type(property_cache), pointer :: new_element_ptr
    type(property_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    complex(8), dimension(property_size) :: data

    next_element => current_element

    allocate(new_element)

    call property_cache_initialize(new_element, num_p_tuples, p_tuples, &
                                   property_size, data)

    new_element_ptr => new_element

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while (next_element%last .eqv. .FALSE.)
       next_element => property_cache_next_element(next_element)
    end do

    next_element%last = .FALSE.
    new_element%next => next_element%next
    next_element%next => new_element

  end subroutine


  function property_cache_already(current_element, num_p_tuples, p_tuples)

    implicit none

    logical :: property_cache_already
    integer :: passedlast, num_p_tuples, i
    type(property_cache), target :: current_element
    type(property_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    
    next_element => current_element
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    property_cache_already = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (property_cache_already .eqv. .FALSE.))

       next_element => property_cache_next_element(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          property_cache_already = p_tuples_compare(num_p_tuples, &
                                   p_tuples_standardorder(next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples_st_order)

! if (property_cache_already .eqv. .true.) then
! 
! if (num_p_tuples == 3) then
! 
! if (p_tuples(2)%n_perturbations == 2 .AND. p_tuples(3)%n_perturbations == 2) then
! 
! if ((p_tuples(2)%pid(1) == 6) .AND. (p_tuples(2)%pid(1) == 6) &
!      .AND. (p_tuples(3)%pid(1) == 1)) then
! 
! write(*,*) ' '
! write(*,*) 'pr_cache already match: pristine:'
! 
! do i = 1, num_p_tuples
! 
! write(*,*) p_tuples(i)%pid
! write(*,*) p_tuples(i)%pdim
! write(*,*) p_tuples(i)%freq
! 
! end do
! 
! write(*,*) ' '
! write(*,*) 'st_order:'
! 
! do i = 1, num_p_tuples
! 
! write(*,*) p_tuples(i)%pid
! write(*,*) p_tuples(i)%pdim
! write(*,*) p_tuples(i)%freq
! 
! end do
! 
! 
! write(*,*) ' '
! 
! end if
! 
! end if
! 
! end if
! 
! end if

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (property_cache_already .EQV. .TRUE.) then

!         write(*,*) 'property_cache_already: Found element in cache'

    else

!         write(*,*) 'property_cache_already: Element not in cache'

    end if

  end function

  ! Assumes that p_tuples is in standard order
  subroutine property_cache_getdata(cache, num_p_tuples, p_tuples, property_size, prop)

    implicit none

    logical :: found
    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
               property_size, total_num_perturbations, pr_offset, cache_offset, &
               merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contribution, & 
                                          p_tuples_dimensions, &
                                          p_tuples_dimensions_cacheorder, &
                                          pids_merged_pert, translated_index
    integer, allocatable, dimension(:) :: blk_sizes_merged
    integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:,:) :: indices, blk_sizes
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
    type(property_cache), target :: cache
    type(property_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    type(p_tuple) :: merged_p_tuple
    complex(8), dimension(property_size) :: prop


    next_element => cache
    passedlast = 0
    p_tuples_st_order = p_tuples
!     p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    found = .FALSE.



! if (num_p_tuples == 3) then
! 
! if (p_tuples(2)%n_perturbations == 2 .AND. p_tuples(3)%n_perturbations == 2) then
! 
! if ((p_tuples(2)%pid(1) == 6) .AND. (p_tuples(2)%pid(1) == 6) &
!      .AND. (p_tuples(3)%pid(1) == 1)) then
! 
! write(*,*) ' '
! write(*,*) 'pr_cache getdata:'
! 
! do i = 1, num_p_tuples
! 
! write(*,*) p_tuples(i)%pid
! write(*,*) p_tuples(i)%pdim
! write(*,*) p_tuples(i)%freq
! 
! end do
! 
! 
! 
! end if
! 
! end if
! 
! end if


    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

       next_element => property_cache_next_element(next_element)

       if (next_element%num_p_tuples == num_p_tuples) then

          found = p_tuples_compare(num_p_tuples, p_tuples_standardorder( &
                                   next_element%num_p_tuples, &
                                   next_element%p_tuples), p_tuples)

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (found .eqv. .TRUE.) then

!        write(*,*) 'Getting property_cache data'

       total_num_perturbations = 0

       do i = 1, num_p_tuples
          total_num_perturbations = total_num_perturbations + p_tuples(i)%n_perturbations
       end do

       next_element%p_tuples = p_tuples_standardorder(next_element%num_p_tuples, &
                                                      next_element%p_tuples)


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
          do j = 1, p_tuples(i)%n_perturbations
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

          pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
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

          prop(pr_offset) = &
          prop(pr_offset) + &
          next_element%data(cache_offset) 

       end do

       deallocate(indices)
       deallocate(blks_tuple_info)
       deallocate(merged_blk_info)

    else

       write(*,*) 'Failed to retrieve data in property_cache_getdata: Element not found'

    end if

       deallocate(pids_in_cache)
       deallocate(pids_current_contribution)
       deallocate(pids_merged_pert)
       deallocate(p_tuples_dimensions)
       deallocate(p_tuples_dimensions_cacheorder)
       deallocate(blk_sizes)
       deallocate(blk_sizes_merged)

  end subroutine


  subroutine property_cache_allocate(current_element)

    implicit none

    type(property_cache), pointer :: current_element

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

    ! NEW 2014
  
  
!    function contrib_cache_next_element(current_element) result(next_element)
!
!    implicit none
!
!    type(contrib_cache), target :: current_element
!    type(contrib_cache), pointer :: next_element
!
!    next_element => current_element%next
!
!  end function
!
!  function contrib_cache_outer_next_element(current_element) result(next_element)
!
!    implicit none
!
!    type(contrib_cache_outer), target :: current_element
!    type(contrib_cache_outer), pointer :: next_element
!
!    next_element => current_element%next
!
!  end function
!  
!  
!  subroutine contrib_cache_initialize(new_element, num_p_tuples, p_tuples)
!
!    implicit none
!
!    integer :: num_p_tuples
!    type(contrib_cache) :: new_element
!    type(p_tuple), dimension(num_p_tuples) :: p_tuples
!
!    new_element%last = .TRUE.
!    call p1_cloneto_p2(p_tuples(1), new_element%p_inner)
!    
!    if (num_p_tuples > 1) then
!    
!       call contrib_cache_outer_initialize(new_element%contribs_outer, num_p_tuples - 1, &
!                                           p_tuples(2:num_p_tuples))
!    
!    end if
!
!  end subroutine
!  
!  
!  subroutine contrib_cache_outer_initialize(new_element, num_dmat, outer_p_tuples)
!
!    implicit none
!
!    integer :: num_dmat, i
!    type(contrib_cache_outer) :: new_element
!    type(p_tuple), dimension(num_dmat) :: outer_p_tuples
!
!    new_element%last = .TRUE.
!    
!    do i = 1, num_dmat
!       call p1_cloneto_p2(outer_p_tuples(i), new_element%outer_p_tuples(i))
!    end do
!       
!  end subroutine
!  
!  
!  subroutine contrib_cache_add_element(current_element, num_p_tuples, p_tuples)
!
!    implicit none
!
!    integer :: num_p_tuples
!    type(contrib_cache), target :: current_element
!    type(contrib_cache), pointer :: new_element
!    type(contrib_cache), pointer :: new_element_ptr
!    type(contrib_cache), pointer :: next_element
!    type(p_tuple), dimension(num_p_tuples) :: p_tuples
!    type(p_tuple) :: emptypert
!
!    ! If cache element for inner perturbations already exists, just add outer
!    if (contrib_cache_already_inner(current_element, p_tuples(1))) then
!    
!       next_element => current_element
!    
!       ! Skip to cache element for this inner
!       do while (p_tuple_compare(next_element%p_inner, p_tuples(1)) .EQV. .FALSE.)
!
!          next_element => next_element%next
!           
!       end do
!       
!       if (num_p_tuples > 1) then
!      
!          call contrib_cache_outer_initialize(next_element%contribs_outer, &
!               num_p_tuples - 1, p_tuples(2:num_p_tuples))
!             
!       else
!          
!          call empty_p_tuple(emptypert)
!!           call contrib_cache_outer_add_element(next_element%contribs_outer, 1, emptypert)
!        
!       end if
!
!
!    ! Otherwise, add both inner and outer    
!    else
!    
!       next_element => current_element
!
!       allocate(new_element)
!
!       call contrib_cache_initialize(new_element, num_p_tuples, p_tuples)
!
!       new_element_ptr => new_element
!
!       do while (next_element%last .EQV. .FALSE.)
!          next_element => contrib_cache_next_element(next_element)
!       end do
!
!       next_element%last = .FALSE.
!       new_element%next => next_element%next
!       next_element%next => new_element
!       
!    end if
!
!  end subroutine
!  
!  
!  ! MAYBE SOME WORK REMAINING ON THIS FUNCTION
!  
!  function contrib_cache_already(current_element, num_p_tuples, p_tuples)
!
!    implicit none
!
!    logical :: contrib_cache_already
!    integer :: passedlast, passedlast_outer, num_p_tuples, i
!    type(contrib_cache), target :: current_element
!    type(contrib_cache), pointer :: next_element
!    type(p_tuple), dimension(num_p_tuples) :: p_tuples
!    type(p_tuple) :: emptypert, p_tuple_ord, p_tmp_ord
!    
!    next_element => current_element
!    passedlast = 0
!    call p_tuple_ordered(p_tuples(1), p_tuple_ord)
!    contrib_cache_already = .FALSE.
!
!    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
!    ! COULD THIS BE DONE IN ANOTHER WAY?
!    do while ((passedlast < 2) .AND. (contrib_cache_already .eqv. .FALSE.))
!
!       next_element => contrib_cache_next_element(next_element)
!
!       call p_tuple_ordered(next_element%p_inner, p_tmp_ord)
!       
!       contrib_cache_already = p_tuple_compare(p_tmp_ord, p_tuple_ord)
!
!       if (contrib_cache_already) then
!          
!          if (num_p_tuples > 1) then
!          
!             contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
!                                     num_p_tuples - 1, p_tuples(2:num_p_tuples))
!             
!          else
!          
!             call empty_p_tuple(emptypert)          
!             contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
!                                     1, (/emptypert/))
!                       
!          end if
!       
!       end if
!
!       if (next_element%last .eqv. .TRUE.) then
!          passedlast = passedlast + 1
!       end if
!
!    end do
!
!  end function
!  
!
!  ! NOT DONE WITH THIS FUNCTION
!  
!  function contrib_cache_already_outer(current_element, num_dmat, p_tuples_outer)
!
!    implicit none
!
!    logical :: contrib_cache_already_outer
!    integer :: passedlast, passedlast_outer, num_dmat, i
!    type(contrib_cache_outer), target :: current_element
!    type(contrib_cache_outer), pointer :: next_element
!    type(p_tuple), dimension(num_dmat) :: p_tuples_outer, p_tuples_ord
!    
!    next_element => current_element
!    passedlast = 0
!    call p_tuples_ordered(num_dmat, p_tuples_outer, p_tuples_ord)
!    contrib_cache_already_outer = .FALSE.
!
!    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
!    ! COULD THIS BE DONE IN ANOTHER WAY?
!    do while ((passedlast < 2) .AND. (contrib_cache_already_outer .eqv. .FALSE.))
!
!    ! Make sure the tuples that need to be ordered are ordered
!    
!       next_element => contrib_cache_outer_next_element(next_element)
!
!       if (next_element%num_dmat == num_dmat) then
!       
!       contrib_cache_already_outer = p_tuples_compare(num_dmat, &
!                                     next_element%outer_p_tuples, p_tuples_ord)
!
!                                  
!
!                                  
!       end if
!
!       if (next_element%last .eqv. .TRUE.) then
!          passedlast = passedlast + 1
!       end if
!
!    end do
!
!
!  end function
!
!    ! NOT DONE WITH THIS FUNCTION
!  
!  function contrib_cache_already_inner(current_element, p_inner)
!
!    implicit none
!
!    logical :: contrib_cache_already_inner
!    integer :: passedlast,  i
!    type(contrib_cache), target :: current_element
!    type(contrib_cache), pointer :: next_element
!    type(p_tuple) :: p_inner, p_tuple_ord
!    
!    next_element => current_element
!    passedlast = 0
!    call p_tuple_ordered(p_inner, p_tuple_ord)
!    contrib_cache_already_inner = .FALSE.
!
!    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
!    ! COULD THIS BE DONE IN ANOTHER WAY?
!    do while ((passedlast < 2) .AND. (contrib_cache_already_inner .eqv. .FALSE.))
!
!    ! Make sure the tuples that need to be ordered are ordered
!    
!       next_element => contrib_cache_next_element(next_element)
!      
!       contrib_cache_already_inner = p_tuple_compare(next_element%p_inner, p_tuple_ord)
!
!       if (next_element%last .eqv. .TRUE.) then
!          passedlast = passedlast + 1
!       end if
!
!    end do
!
!
!  end function
!  
!     
!  ! Missing contents (to be taken from non-outer routine below)
!  
!  subroutine contrib_cache_getdata_outer(cache, num_p_tuples, p_tuples, &
!             contrib_size, fock_or_prop, fock, prop)
!
!    implicit none
!    logical :: found, fock_or_prop
!    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
!               contrib_size, total_num_perturbations, pr_offset, cache_offset, &
!               merged_triang_size, merged_nblks
!    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contrib, & 
!                                          p_tuples_dimensions, &
!                                          p_tuples_dimensions_cacheorder, &
!                                          pids_merged_pert, translated_index
!    integer, allocatable, dimension(:) :: blk_sizes_merged
!    integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
!    integer, allocatable, dimension(:,:) :: indices, blk_sizes
!    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
!    type(contrib_cache), target :: cache
!    type(contrib_cache), pointer :: next_element
!    type(contrib_cache_outer), pointer :: next_element_outer
!    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_ord
!    type(p_tuple) :: merged_p_tuple
!    type(openrsp_matrix), optional, dimension(contrib_size) :: fock
!    complex(8), optional, dimension(contrib_size) :: prop
!    
!  end subroutine
!  
!! Missing variable declaration and allocations, otherwise OK
!
!  ! Assumes that p_tuples is in standard order
!  subroutine contrib_cache_getdata(cache, num_p_tuples, p_tuples, contrib_size, fock_or_prop, fock, prop)
!
!    implicit none
!
!    logical :: found, fock_or_prop
!    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
!               contrib_size, total_num_perturbations, pr_offset, cache_offset, &
!               merged_triang_size, merged_nblks, res_offset
!    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contrib, & 
!                                          p_tuples_dimensions, &
!                                          p_tuples_dimensions_cacheorder, &
!                                          pids_merged_pert, translated_index
!    integer, allocatable, dimension(:) :: blk_sizes_merged
!    integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
!    integer, allocatable, dimension(:,:) :: indices, blk_sizes
!    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
!    type(contrib_cache), target :: cache
!    type(contrib_cache), pointer :: next_element
!    type(contrib_cache_outer), pointer :: next_element_outer
!    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_ord
!    type(p_tuple) :: merged_p_tuple
!    type(openrsp_matrix), optional, dimension(contrib_size) :: fock
!    complex(8), optional, dimension(contrib_size) :: prop
!
!
!    next_element => cache
!    passedlast = 0
!    found = .FALSE.
!
!
!    call p_tuples_ordered(num_p_tuples, p_tuples, p_tuples_ord)
!
!    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
!    ! COULD THIS BE DONE IN ANOTHER WAY?
!    do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))
!
!       next_element => contrib_cache_next_element(next_element)
!
!       found = p_tuple_compare(next_element%p_inner, p_tuples_ord(1))
!
!       if (next_element%last) then
!          passedlast = passedlast + 1
!       end if
!
!    end do
!
!    
!    if (found) then
!    
!       next_element_outer => next_element%contribs_outer
!       passedlast = 0
!       found = .FALSE.
!
!       ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
!       ! COULD THIS BE DONE IN ANOTHER WAY?
!       do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))
!
!          next_element_outer => contrib_cache_outer_next_element(next_element_outer)
!
!          if (next_element_outer%num_dmat == num_p_tuples - 1) then
!
!             found = p_tuples_compare(num_p_tuples - 1, next_element_outer%outer_p_tuples, &
!                                      p_tuples_ord(2:num_p_tuples))
!
!          end if
!
!          if (next_element%last) then
!             passedlast = passedlast + 1
!          end if
!
!       end do
!       
!       if (found) then
!
!          total_num_perturbations = 0
!
!          do i = 1, num_p_tuples
!             total_num_perturbations = total_num_perturbations + p_tuples(i)%n_perturbations
!          end do
!
!
!          if (p_tuples(1)%n_perturbations > 0) then
!
!             call p1_cloneto_p2(p_tuples(1), merged_p_tuple)
!
!             do i = 2, num_p_tuples
!
!                ! This can be problematic - consider rewriting merge_p_tuple as subroutine
!                call p1_merge_p2(merged_p_tuple, p_tuples(i), merged_p_tuple)
!
!             end do
!
!          else
!
!             call p1_cloneto_p2(p_tuples(2), merged_p_tuple)
!
!             do i = 3, num_p_tuples
!
!                ! This can be problematic - consider rewriting merge_p_tuple as subroutine
!                call p1_merge_p2(merged_p_tuple, p_tuples(i), merged_p_tuple)
!                
!             end do
!
!          end if
!
!          call p_tuple_ordered(merged_p_tuple, merged_p_tuple)
!          merged_nblks = get_num_blks(merged_p_tuple)
!
!          allocate(merged_blk_info(1,merged_nblks, 3))
!
!          merged_blk_info(1,:,:) = get_blk_info(merged_nblks, merged_p_tuple)
!          merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
!
!          allocate(blk_sizes(num_p_tuples, total_num_perturbations))
!          allocate(blk_sizes_merged(total_num_perturbations))
!
!          blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
!          merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
!
!          allocate(indices(contrib_size, total_num_perturbations))
!
!          call make_triangulated_indices(merged_nblks, merged_blk_info, & 
!               merged_triang_size, indices)
!
!          allocate(translated_index(total_num_perturbations))
!          allocate(pids_in_cache(total_num_perturbations))
!          allocate(pids_current_contrib(total_num_perturbations))
!          allocate(pids_merged_pert(total_num_perturbations))
!          allocate(p_tuples_dimensions(total_num_perturbations))
!          allocate(p_tuples_dimensions_cacheorder(total_num_perturbations))
!
!          k = 1
!
!          do i = 1, num_p_tuples
!             do j = 1, p_tuples(i)%n_perturbations
!!                 pids_current_contrib(k) = p_tuples_ord(i)%perts(j)%pid
!                k = k + 1
!             end do
!          end do
!
!          do i = 1, total_num_perturbations
!!              pids_merged_pert(i) = merged_p_tuple%perts(i)%pid
!          end do
!
!          p_tuples_dimensions = get_ncarray(total_num_perturbations, num_p_tuples, &
!                                            p_tuples_ord)
!
!          allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
!
!          do i = 1, num_p_tuples
!
!             nfields(i) = p_tuples_ord(i)%n_perturbations
!             nblks_tuple(i) = get_num_blks(p_tuples_ord(i))
!             blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples_ord(i))
!             blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
!                                         blks_tuple_info(i,1:nblks_tuple(i),:))
!             blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
!             blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))
!
!          end do
!
!          do i = 1, size(indices, 1)
!
!             res_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
!             (/merged_nblks/), & 
!             (/total_num_perturbations/), (/merged_blk_info/), blk_sizes_merged, &
!             (/merged_triang_size/), &
!             (/indices(i, : )/) )
!
!             do j = 1, total_num_perturbations
!                translated_index(j) = indices(i,pids_current_contrib(j))
!             end do
!
!             if (p_tuples(1)%n_perturbations > 0) then
!
!                cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
!                total_num_perturbations, nblks_tuple, & 
!                nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)
!
!             else
!
!                cache_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
!                total_num_perturbations, nblks_tuple(2:num_p_tuples), & 
!                nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :),&
!                blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
!                translated_index)
!
!             end if
!
!             if (fock_or_prop) then
!             
!             write(*,*) 'dummy matrix in prop caching'
!             
!                fock(res_offset) = &
!!                 fock(res_offset) + &
!                next_element_outer%data_int(cache_offset) 
!                
!             else
!
!                prop(res_offset) = &
!                prop(res_offset) + &
!                next_element_outer%data_ave(cache_offset)              
!             
!             end if
!
!          end do
!
!      
!       else
!
!          write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Element not found'
!
!       end if
!       
!       
!       
!    else
!
!       write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Element not found'
!
!    end if
!
!
!  end subroutine
!
!
!  
!  
!  
!  
!  subroutine contrib_cache_allocate(current_element)
!
!    implicit none
!
!    type(contrib_cache), pointer :: current_element
!
!    allocate(current_element)
!
!    current_element%next => current_element
!    current_element%last = .TRUE.
!
!    current_element%p_inner%n_perturbations = 0
!!     allocate(current_element%p_inner%perts(1))
!    ! Also set values of perturbation?
!    
!  end subroutine
  
  ! END NEW 2014
  
  
end module
