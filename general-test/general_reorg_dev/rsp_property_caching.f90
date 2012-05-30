! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and functions related to caching of contributions to the
! response property tensor.

module rsp_property_caching

  use rsp_field_tuple

  implicit none

  public property_cache_initialize
  public property_cache_next_element
  public property_cache_add_element
  public property_cache_already
  public property_cache_getdata
  public property_cache_allocate

!  Define property contribution cache datatype

 type property_cache

     type(property_cache), pointer :: next
     logical :: last
     integer :: num_p_tuples
     type(p_tuple), allocatable, dimension(:) :: p_tuples
     complex(8), allocatable, dimension(:) :: data ! Property data    

  end type 

  contains

  ! Begin property_cache linked list manipulation/data retrieval routines

  subroutine property_cache_initialize(new_element, num_p_tuples, p_tuples, &
                                       property_size, data)

    implicit none

    integer :: i, num_p_tuples
    type(property_cache) :: new_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    integer :: property_size
    complex(8), dimension(property_size) :: data


    new_element%last = .TRUE.
    new_element%num_p_tuples = num_p_tuples
    allocate(new_element%p_tuples(num_p_tuples))

    do i = 1, num_p_tuples
       call p_tuple_p1_cloneto_p2(p_tuples(i), new_element%p_tuples(i))
    end do

    allocate(new_element%data(property_size))
    new_element%data = data
!     write(*,*) new_element%data

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
    integer :: passedlast, num_p_tuples
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

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (property_cache_already .EQV. .TRUE.) then

        write(*,*) 'property_cache_already: Found element in cache'

    else

        write(*,*) 'property_cache_already: Element not in cache'

    end if

  end function


  subroutine property_cache_getdata(cache, num_p_tuples, p_tuples, property_size, prop)

    implicit none

    logical :: found
    integer :: i, j, k, first, last, passedlast, num_p_tuples, &
               property_size, total_num_perturbations
    integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contribution, & 
                                          p_tuples_dimensions, &
                                          p_tuples_dimensions_cacheorder
    integer, allocatable, dimension(:,:) :: indices
    type(property_cache), target :: cache
    type(property_cache), pointer :: next_element
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order
    complex(8), dimension(property_size) :: prop
complex(8), dimension(property_size) :: p_debug
integer :: ind_debug, ind2_debug


    next_element => cache
    passedlast = 0
    p_tuples_st_order = p_tuples_standardorder(num_p_tuples, p_tuples)
    found = .FALSE.

p_debug = 0.0

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

       next_element => property_cache_next_element(next_element)

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

       write(*,*) 'Getting property_cache data'

       total_num_perturbations = 0

       do i = 1, num_p_tuples
          total_num_perturbations = total_num_perturbations + p_tuples(i)%n_perturbations
       end do

       next_element%p_tuples = p_tuples_standardorder(next_element%num_p_tuples, &
                                                      next_element%p_tuples)


       allocate(pids_in_cache(total_num_perturbations))
       allocate(pids_current_contribution(total_num_perturbations))
       allocate(p_tuples_dimensions(total_num_perturbations))
       allocate(p_tuples_dimensions_cacheorder(total_num_perturbations))

       k = 1

       do i = 1, num_p_tuples

          do j = 1, p_tuples(i)%n_perturbations

             pids_in_cache(k) = next_element%p_tuples(i)%pid(j)
             pids_current_contribution(k) = p_tuples_st_order(i)%pid(j)
             p_tuples_dimensions_cacheorder(k) = p_tuples_st_order(i)%pdim(j)


          k = k + 1

          end do

       end do

       p_tuples_dimensions = get_ncarray(total_num_perturbations, num_p_tuples, &
                                         p_tuples_st_order)

       ! Making indices
       ! Note (MaR): This can take a lot of memory: 
       ! Consider splitting index generation into several
       ! steps in order to save memory - e.g. "every n ranks"
       ! starts a new loop (note that this probably means that 
       ! this procedure needs to be recursive)

       allocate(indices(product(p_tuples_dimensions), total_num_perturbations))

       call make_indices(total_num_perturbations, 1, p_tuples_dimensions_cacheorder, &
                         0, indices)

       do i = 1, size(indices, 1)



! ind_debug = get_one_tensor_offset(total_num_perturbations, indices(i, :), &
! pids_current_contribution, p_tuples_dimensions)
! 
! ind2_debug = get_one_tensor_offset(total_num_perturbations, indices(i, :), &
! pids_in_cache, p_tuples_dimensions)
! 
! write(*,*) indices(i,:), ': ', ind_debug, ' and ', ind2_debug
! 
! p_debug(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
! pids_current_contribution, p_tuples_dimensions)) = &
! p_debug(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
! pids_current_contribution, p_tuples_dimensions)) + &
! next_element%data(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
! pids_in_cache, p_tuples_dimensions)) 


          ! To which element in the cached data does that 
          ! cache p_tuples index tuple correspond?
          ! Get that element
          ! To which element in the property tensor does that correspond? 
          ! Put the element there

          prop(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
               pids_current_contribution, p_tuples_dimensions)) = &
          prop(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
               pids_current_contribution, p_tuples_dimensions)) + &
          next_element%data(get_one_tensor_offset(total_num_perturbations, indices(i, :), &
               pids_in_cache, p_tuples_dimensions)) 

       end do

! write(*,*) 'ind_debug'
! write(*,*) ind_debug
! 
! write(*,*) 'ind2_debug'
! write(*,*) ind2_debug
! 
 write(*,*) 'Got property cache data'
! write(*,*) real(p_debug)
!  call print_rsp_tensor_stdout(total_num_perturbations,total_num_perturbations, &
!                               p_tuples_dimensions, p_debug, 1)

       deallocate(indices)

    else

       write(*,*) 'Failed to retrieve data in property_cache_getdata: Element not found'

    end if

       deallocate(pids_in_cache)
       deallocate(pids_current_contribution)
       deallocate(p_tuples_dimensions)
       deallocate(p_tuples_dimensions_cacheorder)


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

  ! End property_cache linked list manipulation/data retrieval routines


end module