! Copyright 2012      Magnus Ringholm
!           2012      Dan Jonsson
!           2009-2011 Radovan Bast
!           2009-2011 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_general

!> General response routines. This module organizes, computes and prints
!> response function tensors.
module rsp_general

  use matrix_defop
  use matrix_backend
  use rsp_contribs
  use rsp_equations


  implicit none

  ! Note(MaR): Make sure that this private statement is placed correctly and then delete
  ! this comment. Should it maybe come after e.g. the type definitions?

  ! private

  !unit number for IO in this file
  !radovan: a unit nr this high may not work on all compilers and/or processors
  !         we should be careful going beyond 100
  integer, parameter :: iounit = 345645
  

  ! Zero matrix

  type(matrix) :: zromtrx

  ! Define perturbation tuple datatype

  type p_tuple

     integer :: n_perturbations ! Number of perturbations
     integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
     character(4), allocatable, dimension(:) :: plab ! Perturbation labels
     integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
     complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
     ! Add other perturbation identification info as needed

  end type

  ! Define perturbed S, D, or F linked list datatype

  type SDF

     type(SDF), pointer :: next
     logical :: last
     ! Should all of the data attributes be pointers too?
     type(p_tuple) :: perturb
     type(matrix), allocatable, dimension(:) :: data ! (Perturbed) matrix data

  ! Note(MaR): Like for property_cache, a good extension here would be 
  ! to let the data attribute point to a function/routine that manages retrieval 
  ! (and storage) of values. This will allow better management of very large pieces 
  ! of data, for instance by letting the function manage whether some (large) piece 
  ! of data is stored on disk (and retrieved from there) or if, for smaller pieces 
  ! of data or, in case of parallel implementations, the data is stored in memory.

  end type

  ! Define property contribution cache datatype

  type property_cache

     type(property_cache), pointer :: next
     logical :: last
     integer :: num_p_tuples
     type(p_tuple), allocatable, dimension(:) :: p_tuples
     complex(8), allocatable, dimension(:) :: data ! Property data    

  end type 

  ! MR: Like for SDF, a future extension could be made here in which the data
  ! attribute points to a function that manages retrieval (and storage) of values.

  contains

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

  function p_tuple_extend(pert, ext)

    implicit none

    type(p_tuple) :: pert, ext, p_tuple_extend
    integer :: i

    allocate(p_tuple_extend%pdim(pert%n_perturbations + 1))
    allocate(p_tuple_extend%plab(pert%n_perturbations + 1))
    allocate(p_tuple_extend%pid(pert%n_perturbations + 1))
    allocate(p_tuple_extend%freq(pert%n_perturbations + 1))

    if (pert%n_perturbations == 0) then

       p_tuple_extend%n_perturbations = pert%n_perturbations + 1
       p_tuple_extend%pdim = (/ext%pdim(:)/)
       p_tuple_extend%plab = (/ext%plab(:)/)
       p_tuple_extend%pid = (/ext%pid(:)/)
       p_tuple_extend%freq = (/ext%freq(:)/)

    else

       p_tuple_extend%n_perturbations = pert%n_perturbations + 1
       p_tuple_extend%pdim = (/(pert%pdim(i), i = 1, pert%n_perturbations), ext%pdim(:)/)
       p_tuple_extend%plab = (/(pert%plab(i), i = 1, pert%n_perturbations), ext%plab(:)/)
       p_tuple_extend%pid = (/(pert%pid(i), i = 1, pert%n_perturbations), ext%pid(:)/)
       p_tuple_extend%freq = (/(pert%freq(i), i = 1, pert%n_perturbations), ext%freq(:)/)
 
end if

  end function


  function p_tuple_getone(pert, which)

    implicit none

    type(p_tuple) :: pert, p_tuple_getone
    integer :: which

    allocate(p_tuple_getone%pdim(1))
    allocate(p_tuple_getone%plab(1))
    allocate(p_tuple_getone%pid(1))
    allocate(p_tuple_getone%freq(1))

    p_tuple_getone%n_perturbations = 1
    p_tuple_getone%pdim = (/pert%pdim(which)/)
    p_tuple_getone%plab = (/pert%plab(which)/)
    p_tuple_getone%pid = (/pert%pid(which)/)
    p_tuple_getone%freq = (/pert%freq(which)/)

  end function


  function p_tuple_remove_first(pert)

    implicit none

    type(p_tuple) :: pert, p_tuple_remove_first

    allocate(p_tuple_remove_first%pdim(pert%n_perturbations - 1))
    allocate(p_tuple_remove_first%plab(pert%n_perturbations - 1))
    allocate(p_tuple_remove_first%pid(pert%n_perturbations - 1))
    allocate(p_tuple_remove_first%freq(pert%n_perturbations - 1))

    if (pert%n_perturbations > 1) then

       p_tuple_remove_first%n_perturbations = pert%n_perturbations - 1
       p_tuple_remove_first%pdim = (/pert%pdim(2:pert%n_perturbations)/)
       p_tuple_remove_first%plab = (/pert%plab(2:pert%n_perturbations)/)
       p_tuple_remove_first%pid = (/pert%pid(2:pert%n_perturbations)/)
       p_tuple_remove_first%freq = (/pert%freq(2:pert%n_perturbations)/)

    else

       p_tuple_remove_first%n_perturbations = 0
       p_tuple_remove_first%pdim = (/0/)
       p_tuple_remove_first%plab = (/'NUTN'/)
       p_tuple_remove_first%pid = (/0/)
       p_tuple_remove_first%freq = (/0.0/)

    end if



  end function

  function merge_p_tuple(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2, merge_p_tuple

    allocate(merge_p_tuple%pdim(p1%n_perturbations + p2%n_perturbations))
    allocate(merge_p_tuple%plab(p1%n_perturbations + p2%n_perturbations))
    allocate(merge_p_tuple%pid(p1%n_perturbations + p2%n_perturbations))
    allocate(merge_p_tuple%freq(p1%n_perturbations + p2%n_perturbations))

    ! NOTE (MaR): ARE THESE CASE DISTINCTIONS UNNECESSARY? CONSIDER REWRITE.

    if (p1%n_perturbations > 0 .AND. p2%n_perturbations > 0) then

       merge_p_tuple%n_perturbations = p1%n_perturbations + p2%n_perturbations
       merge_p_tuple%pdim = (/p1%pdim(:), p2%pdim(:)/)
       merge_p_tuple%plab = (/p1%plab(:), p2%plab(:)/)
       merge_p_tuple%pid = (/p1%pid(:), p2%pid(:)/)
       merge_p_tuple%freq = (/p1%freq(:), p2%freq(:)/)

    elseif (p1%n_perturbations > 0 .AND. p2%n_perturbations == 0) then

       merge_p_tuple%n_perturbations = p1%n_perturbations
       merge_p_tuple%pdim = p1%pdim(:)
       merge_p_tuple%plab = p1%plab(:)
       merge_p_tuple%pid = p1%pid(:)
       merge_p_tuple%freq = p1%freq(:)

    elseif (p1%n_perturbations == 0 .AND. p2%n_perturbations > 0) then

       merge_p_tuple%n_perturbations = p2%n_perturbations
       merge_p_tuple%pdim = p2%pdim(:)
       merge_p_tuple%plab = p2%plab(:)
       merge_p_tuple%pid = p2%pid(:)
       merge_p_tuple%freq = p2%freq(:)

    elseif (p1%n_perturbations == 0 .AND. p2%n_perturbations == 0) then

       merge_p_tuple = get_emptypert()
       ! MaR: KEEP NEXT LINES FOR REVERSION IN CASE get_emptypert() DOESN'T WORK
       ! merge_p_tuple%n_perturbations = 0
       ! merge_p_tuple%pdim = (/0/)
       ! merge_p_tuple%plab = (/'NUTN'/)
       ! merge_p_tuple%pid = (/0/)
       ! merge_p_tuple%freq = (/0.0/)

    else

       write(*,*) 'Error in merge_p_tuple: Unrecognized size of p1 or p2 or both:', &
                   p1%n_perturbations, p2%n_perturbations

    end if

  end function

  subroutine p_tuple_p1_cloneto_p2(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2

    p2%n_perturbations = p1%n_perturbations

    allocate(p2%pdim(p1%n_perturbations)) 
    allocate(p2%plab(p1%n_perturbations))
    allocate(p2%pid(p1%n_perturbations))
    allocate(p2%freq(p1%n_perturbations))

    p2%pdim = p1%pdim
    p2%plab = p1%plab
    p2%pid = p1%pid
    p2%freq = p1%freq

  end subroutine

  subroutine p_tuple_deallocate(p1)

    type(p_tuple) :: p1

    p1%n_perturbations = 0

    deallocate(p1%pdim) 
    deallocate(p1%plab)
    deallocate(p1%pid)
    deallocate(p1%freq)

  end subroutine


  function get_emptypert() result(emptypert)

    implicit none

      type(p_tuple) :: emptypert

      emptypert%n_perturbations = 0
      allocate(emptypert%pdim(0))    
      allocate(emptypert%plab(0))
      allocate(emptypert%pid(0))
      allocate(emptypert%freq(0))

  end function


! Compare two perturbation tuples to each other 
! Assumes that input is already in standard order

  function p_tuple_p1_lt_p2(p1, p2)

    implicit none

    logical :: p_tuple_p1_lt_p2
    integer :: i
    type(p_tuple) :: p1, p2

    ! NOTE (MaR): COULD THERE BE FALSE NEGATIVES HERE? 
    p_tuple_p1_lt_p2 = .FALSE.


    ! Compare number of perturbations
    ! REMEMBER: INCREASING ORDER OF DIFFERENTIATION
    if (p1%n_perturbations < p2%n_perturbations) then

       p_tuple_p1_lt_p2 = .TRUE.

    elseif (p1%n_perturbations == p2%n_perturbations) then

       do i = 1, p1%n_perturbations

          ! Compare dimensionality
          ! REMEMBER: DECREASING ORDER OF DIMENSIONALITY
          if (p1%pdim(i) > p2%pdim(i)) then

             p_tuple_p1_lt_p2 = .TRUE.
             exit

          elseif (p1%pdim(i) == p2%pdim(i)) then

             if (llt(p1%plab(i), p2%plab(i)) .eqv. .TRUE.) then

                p_tuple_p1_lt_p2 = .TRUE.  
                exit

             elseif (p1%plab(i) == p2%plab(i)) then

                ! NOTE (MaR): IS IT SUFFICIENTLY GENERAL TO COMPARE ONLY THE REAL PART OF
                ! THE FREQS.? WHICH CASES WILL INCLUDE COMPLEX FREQS. IN THE PERTURBATIONS?
                if (real(p1%freq(i)) < real(p2%freq(i))) then

                   p_tuple_p1_lt_p2 = .TRUE.  
                   exit

                end if

                if (p_tuple_p1_lt_p2 .eqv. .TRUE.) exit

             end if

          end if

       end do

    end if

  end function


  ! FIXME (MaR): THIS FUNCTION IS POORLY WRITTEN AND, ALTHOUGH IT
  ! SEEMS TO BE WORKING, IS BULKY AND MAY CONTAIN
  ! DUPLICATION OF WORK. CONSIDER REWRITING.

  function p_tuple_standardorder(pert) result(p_tuple_st)

    implicit none

    type(p_tuple) :: pert, p_tuple_st
    integer :: i, j, new_minimum, position_in_result
    integer :: temporary_pdim, temporary_pid, current_first_position, current_last_position
    integer :: current_minimum_pdim
    character(4) :: temporary_plab, current_minimum_plab
    complex(8) :: temporary_freq, current_minimum_freq

    call p_tuple_p1_cloneto_p2(pert, p_tuple_st)

    position_in_result = 1

    do i = position_in_result, p_tuple_st%n_perturbations

       current_minimum_pdim = p_tuple_st%pdim(i)
       current_minimum_plab = p_tuple_st%plab(i)
       current_minimum_freq = p_tuple_st%freq(i)
       new_minimum = i

       do j = i + 1, p_tuple_st%n_perturbations

          if (p_tuple_st%pdim(j) > current_minimum_pdim) then

             if (p_tuple_st%pdim(j) == current_minimum_pdim) then

                if (lle(p_tuple_st%plab(j), current_minimum_plab) .EQV. .TRUE.) then

                   if (p_tuple_st%plab(j) == current_minimum_plab) then

                      ! NOTE (MaR): IS IT SUFFICIENTLY GENERAL TO COMPARE
                      ! ONLY THE REAL PART OF THE FREQS.?
                      if (abs(p_tuple_st%freq(j)) < abs(current_minimum_freq)) then

                         current_minimum_freq = p_tuple_st%freq(j)
                         new_minimum = j

                      end if

                   else

                      current_minimum_plab = p_tuple_st%plab(j)
                      new_minimum = j

                   end if

                end if

             else

                current_minimum_pdim = p_tuple_st%pdim(j)
                new_minimum = j

             end if

          end if

       end do

       temporary_pdim = p_tuple_st%pdim(new_minimum)
       temporary_plab = p_tuple_st%plab(new_minimum)
       temporary_pid = p_tuple_st%pid(new_minimum)
       temporary_freq = p_tuple_st%freq(new_minimum)

       p_tuple_st%pdim(new_minimum) = p_tuple_st%pdim(i)
       p_tuple_st%plab(new_minimum) = p_tuple_st%plab(i)
       p_tuple_st%pid(new_minimum) = p_tuple_st%pid(i)
       p_tuple_st%freq(new_minimum) = p_tuple_st%freq(i)

       p_tuple_st%pdim(i) = temporary_pdim
       p_tuple_st%plab(i) = temporary_plab
       p_tuple_st%pid(i) = temporary_pid
       p_tuple_st%freq(i) = temporary_freq

       position_in_result = position_in_result + 1

    end do

    current_first_position = 1
    current_last_position = 1

    do while (current_last_position <= p_tuple_st%n_perturbations)

       if (current_last_position < p_tuple_st%n_perturbations) then

          do while ((p_tuple_st%plab(current_last_position) == &
                     p_tuple_st%plab(current_first_position)))

             current_last_position = current_last_position + 1
             if (current_last_position > p_tuple_st%n_perturbations) exit

          end do

          current_last_position = current_last_position - 1

       end if

       do i = current_first_position, current_last_position, 1

          current_minimum_freq = p_tuple_st%freq(i)
          new_minimum = i

          do j = i + 1, current_last_position

             if (abs(p_tuple_st%freq(j)) < abs(current_minimum_freq)) then

                current_minimum_freq = p_tuple_st%freq(j)
                new_minimum = j

             end if

          end do

          temporary_pdim = p_tuple_st%pdim(new_minimum)
          temporary_plab = p_tuple_st%plab(new_minimum)
          temporary_pid = p_tuple_st%pid(new_minimum)
          temporary_freq = p_tuple_st%freq(new_minimum)

          p_tuple_st%pdim(new_minimum) = p_tuple_st%pdim(i)
          p_tuple_st%plab(new_minimum) = p_tuple_st%plab(i)
          p_tuple_st%pid(new_minimum) = p_tuple_st%pid(i)
          p_tuple_st%freq(new_minimum) = p_tuple_st%freq(i)

          p_tuple_st%pdim(i) = temporary_pdim
          p_tuple_st%plab(i) = temporary_plab
          p_tuple_st%pid(i) = temporary_pid
          p_tuple_st%freq(i) = temporary_freq

       end do

       current_last_position = current_last_position + 1
       current_first_position = current_last_position

    end do

  end function


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


  function p_tuples_standardorder(num_p_tuples, p_tuples) result(p_tuples_st)

    implicit none

    integer :: num_p_tuples
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st
    type(p_tuple) :: temporary_pert
    integer :: i, j, k, new_minimum, max_order_curr
    integer :: temporary_pdim, temporary_pid, current_first_position, &
               current_last_position
    character(4) :: temporary_plab
    character(4), dimension(:), allocatable :: min_plab_curr
    complex(8) :: temporary_freq
    complex(8), dimension(:), allocatable :: current_minimum_freq

    do i = 1, num_p_tuples

       call p_tuple_p1_cloneto_p2(p_tuples(i), p_tuples_st(i))

    end do

    do i = 2, num_p_tuples

       new_minimum = i

       do j = i + 1, num_p_tuples

          if (p_tuple_p1_lt_p2(p_tuple_standardorder(p_tuples_st(j)), &
              p_tuple_standardorder(p_tuples_st(new_minimum)))) then

             new_minimum = j

          end if

       end do

       call p_tuple_p1_cloneto_p2(p_tuples_st(new_minimum),temporary_pert)
       call p_tuple_deallocate(p_tuples_st(new_minimum))

       call p_tuple_p1_cloneto_p2(p_tuples_st(i),p_tuples_st(new_minimum))
       call p_tuple_deallocate(p_tuples_st(i))

       call p_tuple_p1_cloneto_p2(temporary_pert,p_tuples_st(i))
       call p_tuple_deallocate(temporary_pert)

    end do

  end function


  function property_cache_next_element(current_element) result(next_element)

    implicit none

    type(property_cache), target :: current_element
    type(property_cache), pointer :: next_element

    next_element => current_element%next

  end function


  function p_tuples_compare(num_p_tuples, p_tuples, p_tuples_st_order)

    implicit none

    logical :: p_tuples_compare, elem_by_elem_isequivalent
    integer ::  num_p_tuples, i
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order

    p_tuples_compare = .FALSE.
    elem_by_elem_isequivalent = .TRUE.

    do i = 1, num_p_tuples

       elem_by_elem_isequivalent = elem_by_elem_isequivalent .AND. &
                                   p_tuple_compare(p_tuple_standardorder(p_tuples(i)), &
                                   p_tuple_standardorder(p_tuples_st_order(i)))

    end do

    if (elem_by_elem_isequivalent .eqv. .TRUE.) then

       p_tuples_compare = .TRUE.

    end if

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


  function get_one_tensor_offset(total_num_perturbations, indices, pids, dims)

    implicit none

    integer :: total_num_perturbations, i, get_one_tensor_offset
    integer, dimension(total_num_perturbations) ::  indices, pids, dims

    get_one_tensor_offset = 1

    do i = 1, total_num_perturbations

       get_one_tensor_offset = get_one_tensor_offset + &
       (indices(i)- 1)*product(dims(pids(i):total_num_perturbations))/dims(pids(i))

    end do

  end function


  function get_ncarray(total_order, num_p_tuples, p_tuples)

    implicit none

    integer :: total_order, num_p_tuples, i, j, k
    integer, dimension(total_order) :: get_ncarray
    type(p_tuple), dimension(num_p_tuples) :: p_tuples

    do i = 1, total_order
       do j = 1, num_p_tuples
          do k = 1, p_tuples(j)%n_perturbations

             if (p_tuples(j)%pid(k) == i) then
                get_ncarray(i) = p_tuples(j)%pdim(k)
             end if

          end do
       end do

    end do

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

  ! End property_cache linked list manipulation/data retrieval routines


  ! Begin SDF linked list manipulation/data retrieval routines

  subroutine sdf_init(new_element, pert, data)

    implicit none

    type(SDF) :: new_element
    type(p_tuple) :: pert
    type(matrix), dimension(product(pert%pdim)) :: data
    integer :: i

    new_element%last = .TRUE.

    call p_tuple_p1_cloneto_p2(pert, new_element%perturb)
    allocate(new_element%data(product(pert%pdim)))

    do i = 1, product(pert%pdim)
    
       call mat_nullify(new_element%data(i))
       new_element%data(i)%nrow = data(i)%nrow
       new_element%data(i)%ncol = data(i)%nrow
       new_element%data(i)%closed_shell = .true. !*2 on tr(A,B) and dot(A,B)
       new_element%data(i)%magic_tag = 825169837 !mark as set-up
       call mat_alloc(new_element%data(i))
       call mat_axpy((0.0d0, 0.0d0), new_element%data(i), &
                      .false., .true., new_element%data(i))

       new_element%data(i) = data(i)

    end do
 
  end subroutine


  subroutine sdf_standardorder(pert_tuple)

    implicit none

    type(p_tuple) :: pert_tuple
    integer :: i, j,  new_minimum, position_in_result
    integer :: temporary_pdim, temporary_pid
    character(4) :: temporary_plab, current_minimum
    complex(8) :: temporary_freq

    position_in_result = 1

    do i = position_in_result, pert_tuple%n_perturbations

       new_minimum = i

       do j = i + 1, pert_tuple%n_perturbations

          if (llt(pert_tuple%plab(j), current_minimum) .EQV. .TRUE.) then

             new_minimum = j

          end if

       end do

       ! NOTE (MaR): CONSIDER REWRITING TO DO THIS SWITCH MORE ELEGANTLY

       temporary_pdim = pert_tuple%pdim(new_minimum)
       temporary_plab = pert_tuple%plab(new_minimum)
       temporary_pid = pert_tuple%pid(new_minimum)
       temporary_freq = pert_tuple%freq(new_minimum)

       pert_tuple%pdim(new_minimum) = pert_tuple%pdim(i)
       pert_tuple%plab(new_minimum) = pert_tuple%plab(i)
       pert_tuple%pid(new_minimum) = pert_tuple%pid(i)
       pert_tuple%freq(new_minimum) = pert_tuple%freq(i)

       pert_tuple%pdim(i) = temporary_pdim
       pert_tuple%plab(i) = temporary_plab
       pert_tuple%pid(i) = temporary_pid
       pert_tuple%freq(i) = temporary_freq

       position_in_result = position_in_result + 1

    end do

  end subroutine


  function sdf_next_element(current_element) result(next_element)

    implicit none

    type(SDF), target :: current_element
    type(SDF), pointer :: next_element

    next_element => current_element%next

  end function


  ! Add element routine
  ! NOTE(MaR): This routine assumes that the pert_tuple and data
  ! is already in standard order

  subroutine sdf_add(current_element, pert_tuple, data) 

    implicit none

    logical :: found_element
    integer :: passedlast, i
    type(SDF), target :: current_element
    type(SDF), pointer :: new_element
    type(SDF), pointer :: new_element_ptr
    type(SDF), pointer :: next_element
    type(p_tuple) :: pert_tuple, p_tuple_st_order
    type(matrix), dimension(product(pert_tuple%pdim)) :: data

    if (sdf_already(current_element, pert_tuple) .eqv. .TRUE.) then

       next_element => current_element
       passedlast = 0
       p_tuple_st_order = p_tuple_standardorder(pert_tuple)
       found_element = .FALSE.

       ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
       ! COULD THIS BE DONE IN ANOTHER WAY?
       do while ((passedlast < 2) .AND. (found_element .eqv. .FALSE.))

          next_element => sdf_next_element(next_element)
          found_element = p_tuple_compare(next_element%perturb, p_tuple_st_order)

          if (next_element%last .eqv. .TRUE.) then
             passedlast = passedlast + 1
          end if

       end do
! 
! write(*,*) 'entering data'
! do i = 1, size(data)
! 
! write(*,*) 'i is', i
! write(*,*) data(i)%elms
! 
! end do
! if (size(data) >= 5) then

! write(*,*) 'component 5 data'
! write(*,*) data(5)%elms
       next_element%data = data
! write(*,*) 'component 5 n elem'
! write(*,*) next_element%data(5)%elms

! end if

    else

       next_element => current_element
       allocate(new_element)
       call sdf_init(new_element, pert_tuple, data)
       new_element_ptr => new_element

       ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
       ! COULD THIS BE DONE IN ANOTHER WAY?
       do while (next_element%last .eqv. .FALSE.)
          next_element => sdf_next_element(next_element)
       end do

       next_element%last = .FALSE.
       new_element%next => next_element%next
       next_element%next => new_element

    end if

  end subroutine


  function plab_compare(n_perturbations, plab1, plab2)

    implicit none

    logical :: plab_compare
    integer :: n_perturbations, i, j
    character(4) :: plab1(n_perturbations), plab2(n_perturbations)

    plab_compare = .TRUE.

    do i = 1, n_perturbations
       if (.NOT.(plab1(i) == plab2(i))) then
          plab_compare = .FALSE.
       end if
    end do

  end function


  function pfreq_compare(n, p1, p2)

    implicit none

    logical :: pfreq_compare
    integer :: i, n
    complex(8), dimension(n) :: p1, p2
    
    pfreq_compare = .TRUE.

    do i = 1, n
       if ((p1(i) == p2(i)) .EQV. .FALSE.) then
          pfreq_compare = .FALSE.
       end if
    end do

  end function


  function p_tuple_compare(p1, p2)

    implicit none

    logical :: p_tuple_compare
    type(p_tuple) :: p1, p2

    if (p1%n_perturbations == p2%n_perturbations) then

       if (plab_compare(p1%n_perturbations, p1%plab, p2%plab) .eqv. .TRUE.) then

          if (pfreq_compare(p1%n_perturbations, p1%freq, p2%freq) .eqv. .TRUE.) then
             ! NOTE (MaR): Commented code is pseudo for comparing other perturbation info
             ! if (p1%other == p2%other) then
             p_tuple_compare = .TRUE.
             ! else
             ! p_tuple_compare = .FALSE.
             ! end if

          else 

             p_tuple_compare = .FALSE.

          end if

       else 

          p_tuple_compare = .FALSE.

       end if

    else 

       p_tuple_compare = .FALSE.

    end if

  end function


  function sdf_already(current_element, pert_tuple)

    implicit none

    logical :: sdf_already
    type(SDF), target :: current_element
    type(SDF), pointer :: next_element
    type(p_tuple) :: pert_tuple, p_tuple_st_order
    integer :: passedlast

    next_element => current_element

    passedlast = 0
    
    p_tuple_st_order = p_tuple_standardorder(pert_tuple)

    sdf_already = .FALSE.

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. (sdf_already .eqv. .FALSE.))

       next_element => sdf_next_element(next_element)
       sdf_already = p_tuple_compare(next_element%perturb, p_tuple_st_order)

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

  end function


  ! Get SDF element
  ! Assumes that pert_tuple and the p_tuple in current_element is in standard order

  subroutine sdf_getdata_s(current_element, pert_tuple, ind, M)

    implicit none

    logical :: found
    type(SDF), target :: current_element
    type(SDF), pointer :: next_element
    type(p_tuple) :: pert_tuple
    type(matrix) :: M
    integer, dimension(pert_tuple%n_perturbations) :: ind
    integer :: i, offset, passedlast

    next_element => current_element


    if (pert_tuple%n_perturbations > 0) then

       offset = 1

       do i = 1, pert_tuple%n_perturbations

          offset = offset + (ind(i) - 1)*product(pert_tuple%pdim &
                            (i:pert_tuple%n_perturbations))/pert_tuple%pdim(i)

       end do
!  write(*,*) 'offset', offset, ' at ind ', ind
    else

       offset = 1

    end if

    found = .FALSE.
    passedlast = 0

    ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
    ! COULD THIS BE DONE IN ANOTHER WAY?
    do while ((passedlast < 2) .AND. .NOT.(found) .eqv. .TRUE.)

       next_element => sdf_next_element(next_element)

       if (next_element%perturb%n_perturbations == pert_tuple%n_perturbations) then

          if (pert_tuple%n_perturbations == 0) then

             found = .TRUE.

          else

             found = p_tuple_compare(next_element%perturb, pert_tuple)

          end if

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (found .eqv. .TRUE.) then

! if (size(next_element%data) == 36) then
! 
!  write(*,*) 'found element in s', pert_tuple%plab
! write(*,*) 'indices are', ind
! write(*,*) 'offset is', offset
! write(*,*) next_element%data(offset)%elms
! end if

! write(*,*)'found element', pert_tuple%plab, offset
! do i =  1, size(next_element%data)
! 
! write(*,*) next_element%data(i)%elms
! 
! end do
! write(*,*) 'data size', size(next_element%data)
       M = next_element%data(offset)

    else

       write(*,*) 'Failed to retrieve data in sdf_getdata: Element not found'

    end if

  end subroutine


  ! Get SDF element
  ! Assumes that pert_tuple and the p_tuple in current_element is in standard order

  function sdf_getdata(current_element, pert_tuple, ind)

    implicit none

    logical :: found
    type(SDF), target :: current_element
    type(SDF), pointer :: next_element
    type(p_tuple) :: pert_tuple
    type(matrix) :: sdf_getdata
    integer, dimension(pert_tuple%n_perturbations) :: ind
    integer :: i, offset, passedlast

    sdf_getdata = 1.0d0*zromtrx
    next_element => current_element

    if (pert_tuple%n_perturbations > 0) then

       offset = 1
       
       do i = 1, pert_tuple%n_perturbations

          offset = offset + (ind(i) - 1)*product(pert_tuple%pdim(i: &
                   pert_tuple%n_perturbations))/pert_tuple%pdim(i)
 
       end do

    else

       offset = 1

    end if

    found = .FALSE.
    passedlast = 0

    do while ((passedlast < 2) .AND. .NOT.(found) .eqv. .TRUE.)

       next_element => sdf_next_element(next_element)

       if (next_element%perturb%n_perturbations == pert_tuple%n_perturbations) then

          if (pert_tuple%n_perturbations == 0) then

             found = .TRUE.

          else

             found = p_tuple_compare(next_element%perturb, pert_tuple)

          end if

       end if

       if (next_element%last .eqv. .TRUE.) then
          passedlast = passedlast + 1
       end if

    end do

    if (found .eqv. .TRUE.) then

!  write(*,*) 'found element in fn', pert_tuple%plab
! write(*,*) 'indices are', ind
! write(*,*) 'offset is', offset
! write(*,*) next_element%data(offset)%elms

       sdf_getdata = next_element%data(offset)

    else

       write(*,*) 'Failed to retrieve data in sdf_getdata: Element not found'

    end if

  end function





  ! Find out if kn rules say that this term should be skipped
  function kn_skip(n_perturbations, pertid, kn)

    implicit none

    logical :: kn_skip, p_tuple_hasfirst
    integer :: n_perturbations, i
    integer, dimension(n_perturbations) :: pertid
    integer, dimension(2) :: kn

    kn_skip = .FALSE.
    p_tuple_hasfirst = .FALSE.

    do i = 1, size(pertid)
       if (pertid(i) == 1) then
          p_tuple_hasfirst = .TRUE.
       end if
    end do

   
    if (p_tuple_hasfirst .eqv. .TRUE.) then

       if (kn(1) < size(pertid)) then

          kn_skip = .TRUE.

       end if

    else

       if (kn(2) < size(pertid)) then

          kn_skip = .TRUE.

       end if

    end if

  end function


  function nc_only(total_order, thisorder, num_p_tuples, p_tuples, ncarray)

    implicit none

    integer :: i, j, total_order, thisorder, num_p_tuples
    integer, dimension(total_order) :: ncarray
    integer, dimension(total_order) :: nc_only
    type(p_tuple), dimension(num_p_tuples) :: p_tuples

    do i = 1, size(ncarray)
       nc_only(i) = 1
    end do

    do i = 1, num_p_tuples
       do j = 1, p_tuples(i)%n_perturbations
          nc_only(p_tuples(i)%pid(j)) = ncarray(p_tuples(i)%pid(j))
       end do
    end do

  end function


  function nc_onlysmall(total_order, thisorder, num_p_tuples, p_tuples, ncarray)

    implicit none

    integer :: i, j, k, total_order, thisorder, num_p_tuples
    integer, dimension(total_order) :: ncarray
    integer, dimension(thisorder) :: nc_onlysmall
    type(p_tuple), dimension(num_p_tuples) :: p_tuples

    k = 1

    do i = 1, num_p_tuples
       do j = 1, p_tuples(i)%n_perturbations

          nc_onlysmall(k) = ncarray(p_tuples(i)%pid(j))
          k = k + 1

       end do
    end do

  end function


  recursive subroutine make_indices(tot_outer, lvl, ncarray, offset, outer_indices)

    implicit none

    integer :: i, j, k, tot_outer, lvl, offset
    integer, dimension(tot_outer) :: ncarray
    integer, dimension(product(ncarray), tot_outer) :: outer_indices

    k = 1

    if (tot_outer > 0) then
       do i = 1, ncarray(lvl)

          if (lvl < tot_outer) then

             call make_indices(tot_outer, lvl + 1, ncarray, &
             k + offset - 1, outer_indices)

          end if

          if (lvl <= tot_outer) then

             do j = 1, product(ncarray(lvl:size(ncarray)))/ncarray(lvl)

                outer_indices(k + offset, lvl) = i
                k = k + 1

             end do

          end if

       end do

    else

    end if

  end subroutine


  function make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)

    implicit none

    integer :: i, j, k, total_num_perturbations, num_p_tuples
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    integer, dimension(total_num_perturbations) :: make_outerwhichpert

    do i = 1, total_num_perturbations

       make_outerwhichpert(i) = 0

    end do

    k = 1

    do i = 2, num_p_tuples
       do j = 1, p_tuples(i)%n_perturbations

          make_outerwhichpert(p_tuples(i)%pid(j)) = k
          k = k + 1

       end do
    end do

  end function


  function get_pidoutersmall(totouter, len_outer, o_orders)

    implicit none

    integer :: totouter, len_outer, i, j, k
    integer, dimension(totouter) :: get_pidoutersmall
    type(p_tuple), dimension(len_outer) :: o_orders

    k = 1

    do i = 1, len_outer
       do j = 1, o_orders(i)%n_perturbations

          get_pidoutersmall(k) = o_orders(i)%pid(j)
          k = k + 1

       end do
    end do

  end function


  subroutine sortdimbypid(total_num_perturbations, totouter, pids, dims, dimsouter, whichs)

    implicit none

    integer :: totouter, total_num_perturbations, s, i, j, whichmax, whatmax
    integer, dimension(totouter) :: b, d, pids, dimsouter
    integer, dimension(total_num_perturbations) :: whichs, dims

    do i = 1, total_num_perturbations

       whichs(i) = 0

    end do

    s = totouter
    j = totouter
    d = pids

    do while (j > 0)

       whatmax = 0

       ! At which index is the pid largest?

       do i = 1, s
          if (d(i) > whatmax) then

             ! It is currently largest at index i
             whatmax = d(i)
             whichmax = i

          end if
       end do

       ! Then, put the dimension of that pid at the current end of the array to be returned

       b(j) = dims(whatmax)

       ! j is the (current) highest outer index

       whichs(j) = whatmax
       j = j - 1
       d(whichmax) = 0

    end do

    dimsouter = b

  end subroutine


  subroutine get_energy(mol, num_p_tuples, total_num_perturbations, p_tuples, density_order, &
                        D, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF) :: D
    type(property_cache) :: cache
    type(matrix), allocatable, dimension(:) :: dens_tuple
    type(rsp_field), allocatable, dimension(:) :: nucpot_pert
    integer :: i, j, k, m, n, num_p_tuples, total_num_perturbations, density_order, &
             property_size, offset, dtup_ind
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, pidoutersmall
    integer, allocatable, dimension(:) :: ncinnersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, contrib
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
    ncouter = nc_only(total_num_perturbations, total_num_perturbations - &
              p_tuples(1)%n_perturbations, num_p_tuples - 1, &
              p_tuples(2:num_p_tuples), ncarray)
    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
                      p_tuples(1), ncarray)

    allocate(dens_tuple(num_p_tuples))
    allocate(nucpot_pert(p_tuples(1)%n_perturbations))
    allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))
    allocate(ncinnersmall(p_tuples(1)%n_perturbations))
    allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))

    ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
                   p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                   p_tuples(2:num_p_tuples), ncarray)
    ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%n_perturbations, &
                   1, p_tuples(1), ncarray)
    pidoutersmall = get_pidoutersmall(total_num_perturbations - &
                    p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                    p_tuples(2:num_p_tuples))

    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(outer_indices(product(ncoutersmall),size(ncoutersmall)))
    allocate(inner_indices(product(ncinnersmall),size(ncinnersmall)))

    if (p_tuples(1)%n_perturbations > 0) then

       allocate(tmp(product(p_tuples(1)%pdim)))
       allocate(contrib(product(p_tuples(1)%pdim)))

    else

       allocate(tmp(1))
       allocate(contrib(1))

    end if

    contrib = 0.0

    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
                      p_tuples(1)%n_perturbations, pidoutersmall, &
                      ncarray, ncoutersmall, o_whichpert)


    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

    do i = 1, size(o_whichpert)

       if (.NOT.(o_whichpert(i) == 0)) then

       o_wh_forave(o_whichpert(i)) = i

       end if
  
    end do


k = 1

do i = 2, num_p_tuples

do j = 1, p_tuples(i)%n_perturbations

ncoutersmall(k) =  p_tuples(i)%pdim(j)

k = k + 1

end do

end do



    do i = 1, num_p_tuples

       call mat_nullify(dens_tuple(i))
       dens_tuple(i)%nrow = mol%zeromat%nrow
       dens_tuple(i)%ncol = mol%zeromat%ncol
       dens_tuple(i)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
       dens_tuple(i)%magic_tag = 825169837 !mark as set-up
       call mat_alloc(dens_tuple(i))
       call mat_axpy((0.0d0, 0.0d0), dens_tuple(i), .false., .true., dens_tuple(i))

    end do

! if (num_p_tuples >= 3) then
! 	  write(*,*) 'dt2 before', dens_tuple(2)%elms
!           call sdf_getdata_s(D, p_tuples(2), (/ 5 /), dens_tuple(2))
! 	  write(*,*) 'dt2', dens_tuple(2)%elms
! 
! 
! end if
    call make_indices(total_num_perturbations - p_tuples(1)%n_perturbations, &
                      1, ncoutersmall, 0, outer_indices)

    if (p_tuples(1)%n_perturbations > 0) then

       call make_indices(p_tuples(1)%n_perturbations, 1, &
                         p_tuples(1)%pdim, 0, inner_indices)

    end if

    do i = 1, size(outer_indices, 1)

       dtup_ind = 0

       do j = 2, num_p_tuples

! write(*,*) 'getting tuple for', p_tuples(j)%plab
! write(*,*) 'bad indices are', (/ &
!                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
!                              k = 1, p_tuples(j)%n_perturbations) /)
! 
! write(*,*) 'new indices are', outer_indices(i, &
! dtup_ind+1:dtup_ind + p_tuples(j)%n_perturbations)

          call sdf_getdata_s(D, p_tuples(j), outer_indices(i, &
dtup_ind+1:dtup_ind + p_tuples(j)%n_perturbations), dens_tuple(j))

dtup_ind = dtup_ind + p_tuples(j)%n_perturbations

! if (p_tuples(j)%n_perturbations == 2) then
! 
! if ( outer_indices(i,o_wh_forave(p_tuples(j)%pid(1))) == 1 )then
! 
! if ( outer_indices(i,o_wh_forave(p_tuples(j)%pid(2))) == 2 )then
! 
!           call sdf_getdata_s(D, p_tuples(j), (/ 2,1 /), dens_tuple(j))
! write(*,*) 'switcheroo'
! 
! end if
! end if
! 
! end if

       end do

       tmp = 0.0
       contrib = 0.0

       if (num_p_tuples == 1) then

          call rsp_oneave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                         (/ (1, j = 1, p_tuples(1)%n_perturbations) /), & 
                         p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), (/1/)), &
                         contrib)

       elseif (num_p_tuples == 2) then

          call rsp_oneave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                         (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                         p_tuples(1)%pdim, dens_tuple(2), contrib)

       end if

       tmp = tmp + contrib

! if (i == 1) then
! write(*,*) 'AFTER ONEAVE'
! write(*,*) real(contrib)
! write(*,*) ' '
! end if

       contrib = 0.0

       if (num_p_tuples == 1) then

          call rsp_twoave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), &
                          (/1/)), sdf_getdata(D, get_emptypert(), (/1/)) , contrib)

       elseif (num_p_tuples == 2) then

          call rsp_twoave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, dens_tuple(2), &
                          sdf_getdata(D, get_emptypert(), (/1/)) , contrib)

       elseif (num_p_tuples == 3) then

! write(*,*) 'dens tuple 2', dens_tuple(2)%elms
! write(*,*) 'dens tuple 2 plabs', p_tuples(2)%plab
! write(*,*) 'dens tuple 3', dens_tuple(3)%elms
! write(*,*) 'dens tuple 3 plabs', p_tuples(3)%plab

          call rsp_twoave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, dens_tuple(2), dens_tuple(3), contrib)

       end if

       tmp = tmp + contrib

! if (i == 1) then
! write(*,*) 'AFTER TWOAVE'
! write(*,*) real(contrib)
! write(*,*) ' '
! end if

! NOTE (MaR): XCAVE CALL REMOVED FOR NOW
! 
!        contrib = 0.0
! 
!        call rsp_xcave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
!                      num_p_tuples, (/ sdf_getdata(D, get_emptypert(), (/1/)), &
!                      (dens_tuple(k), k = 2, num_p_tuples) /), contrib)
! 
! 
!        tmp = tmp + contrib

       if (p_tuples(1)%n_perturbations > 0) then

          do j = 1, size(inner_indices, 1)

             offset = get_one_tensor_offset( &
                      sum( (/ (p_tuples(k)%n_perturbations, k=1, num_p_tuples ) /) ), &
                      (/ inner_indices(j,:), outer_indices(i,:) /), &
                      (/ (p_tuples(k)%pid, k=1, num_p_tuples ) /), ncarray)

! write(*,*) 'indices inner', inner_indices(j,:), 'and outer', outer_indices(i,:)
! write(*,*) 'offset', offset
! write(*,*) 'tmp at j is', tmp(j)

             prop(offset) = prop(offset) + tmp(j)
             prop_forcache(offset) = prop_forcache(offset) + tmp(j)

          end do

       else

          offset = get_one_tensor_offset( &
                   sum( (/ (p_tuples(k)%n_perturbations, k=2, num_p_tuples ) /) ), &
                   (/ outer_indices(i,:) /), &
                   (/ (p_tuples(k)%pid, k=2, num_p_tuples ) /), ncarray)

          prop(offset) = prop(offset) + tmp(1)
          prop_forcache(offset) = prop_forcache(offset) + tmp(1)

       end if

       end do

    else

       do i = 1, p_tuples(1)%n_perturbations

          nucpot_pert(i) = rsp_field(p_tuples(1)%plab(i), p_tuples(1)%freq(i), 1, &
                                     p_tuples(1)%pdim(i))

       end do

!        write(*,*) 'all indices inner'

       tmp = 0.0
       contrib = 0.0

       call rsp_nucpot(nucpot_pert, contrib) 
       tmp = tmp + contrib

! write(*,*) 'AFTER NUCPOT'
! write(*,*) real(contrib(1))
! write(*,*) ' '

       contrib = 0.0

       call rsp_oneave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , contrib)

       tmp = tmp + contrib

!  write(*,*) 'AFTER ONEAVE'
! 
! call print_rsp_tensor_stdout(total_num_perturbations,total_num_perturbations, &
!                               ncarray, contrib, 1)
! write(*,*) ' '

       contrib = 0.0

       call rsp_twoave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , &
                       sdf_getdata(D, get_emptypert(), (/1/)) , contrib)

       tmp = tmp + 0.5*(contrib)

! write(*,*) 'AFTER TWOAVE'
! write(*,*) 0.5*real(contrib(1))
! write(*,*) ' '

! NOTE (MaR): XCAVE CALL REMOVED FOR NOW
 
!        contrib = 0.0
!
!        call rsp_xcave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
!                      1, (/ sdf_getdata(D, get_emptypert(), (/1/)) /), contrib)
! 
!        tmp = tmp + contrib

       prop =  prop + tmp
       prop_forcache = prop_forcache + tmp

    end if


    call property_cache_add_element(cache, num_p_tuples, p_tuples, &
                                    property_size, prop_forcache)    

!  write(*,*) 'energy contribution'
!  call print_rsp_tensor_stdout(total_num_perturbations,total_num_perturbations, &
!                               ncarray, prop_forcache, 1)

    deallocate(nucpot_pert)
    deallocate(dens_tuple)
    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(tmp)
    deallocate(contrib)

  end subroutine

  ! Calculate and add all the energy contributions

  recursive subroutine rsp_ener(mol, pert, total_num_perturbations, kn, num_p_tuples, &
                                p_tuples, density_order, D, property_size, cache, prop)

    implicit none

    logical :: e_knskip
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF) :: D
    type(property_cache) :: cache
    complex(8), dimension(property_size) :: prop

    if (pert%n_perturbations >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'

    if (p_tuples(1)%n_perturbations == 0) then

       call rsp_ener(mol, p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples, (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
       density_order, D, property_size, cache, prop)

    else


       call rsp_ener(mol, p_tuple_remove_first(pert), total_num_perturbations,  &
       kn, num_p_tuples, (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
       p_tuples(2:size(p_tuples))/), density_order, D, property_size, cache, prop)

    end if
    
       ! 2. Differentiate all of the contraction densities in turn

       ! Find the number of terms

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%n_perturbations == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call rsp_ener(mol, p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples, t_new, density_order + 1, D, property_size, cache, prop)

       end do


       ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_ener(mol, p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
       density_order + 1, D, property_size, cache, prop)


    ! At the final recursion level: Calculate the contribution (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allow it)

    else


    e_knskip = .FALSE.

       write(*,*) 'Getting energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

             write(*,*) 'D ', p_tuples(i)%pid

             if(kn_skip(p_tuples(i)%n_perturbations, p_tuples(i)%pid, kn) .EQV. .TRUE.) then

                e_knskip = .TRUE.

             end if
          
          elseif (i == 1) then

             write(*,*) 'E ', p_tuples(i)%pid

          end if


       end do

       if (e_knskip .EQV. .FALSE.) then

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)
          
          write(*,*) 'Evaluating property_cache_already'

          if (property_cache_already(cache, num_p_tuples, p_tuples) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             write(*,*) 'Getting values from cache'


             ! NOTE (MaR): EVERYTHING IS IN STANDARD ORDER IN 
             ! THIS CALL (LIKE property_cache_getdata ASSUMES)
             call property_cache_getdata(cache, num_p_tuples, &
                  p_tuples_standardorder(num_p_tuples, p_tuples), property_size, prop)

             write(*,*) ' '
       
          else

             call get_energy(mol, num_p_tuples, total_num_perturbations, & 
                  (/ (p_tuple_standardorder(p_tuples(i)) , i = 1, num_p_tuples ) /), &
                  density_order, D, property_size, cache, prop)

                  write(*,*) 'Calculated energy contribution'
                  write(*,*) ' '

          end if

       else

          write(*,*) 'Energy contribution was k-n skipped'
          write(*,*) ' '

       end if

    end if

  end subroutine


  recursive function derivative_superstructure_getsize(mol, pert, kn, &
                     primed, current_derivative_term) result(superstructure_size)

    implicit none

    logical :: primed
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term
    integer, dimension(2) :: kn
    integer :: i, superstructure_size

    superstructure_size = 0

    if (pert%n_perturbations > 0) then

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol,p_tuple_remove_first(pert), kn, primed, &
                             (/p_tuple_extend(current_derivative_term(1), &
                             p_tuple_getone(pert, 1)), current_derivative_term(2:3)/))

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol, p_tuple_remove_first(pert), kn, primed, &
                             (/current_derivative_term(1), &
                             p_tuple_extend(current_derivative_term(2), &
                             p_tuple_getone(pert,1)), current_derivative_term(3)/))

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol, p_tuple_remove_first(pert), kn, primed, &
                             (/current_derivative_term(1:2), &
                             p_tuple_extend(current_derivative_term(3), &
                             p_tuple_getone(pert, 1))/))

    else

       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%n_perturbations <= kn(2)) .AND.&
              current_derivative_term(2)%n_perturbations <= kn(2) ) .AND. &
              current_derivative_term(3)%n_perturbations <= kn(2) ) .eqv. .TRUE.) then

             superstructure_size = 1

          else

             superstructure_size = 0

          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%n_perturbations, current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%n_perturbations, current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%n_perturbations, current_derivative_term(3)%pid, kn) ) .eqv. &
                  .FALSE.) then

             superstructure_size = 1

          else

             superstructure_size = 0

          end if

       end if

    end if

  end function


  recursive subroutine derivative_superstructure(mol, pert, kn, primed, &
                       current_derivative_term, superstructure_size, & 
                       new_element_position, derivative_structure)

    implicit none

    logical :: primed
    integer :: i, superstructure_size, new_element_position
    integer, dimension(2) :: kn    
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term
    type(p_tuple), dimension(superstructure_size, 3) :: derivative_structure

    if (pert%n_perturbations > 0) then

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/p_tuple_extend(current_derivative_term(1), p_tuple_getone(pert, 1)), &
            current_derivative_term(2:3)/), superstructure_size, new_element_position, &
            derivative_structure)

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/current_derivative_term(1), p_tuple_extend(current_derivative_term(2), &
            p_tuple_getone(pert, 1)), current_derivative_term(3)/), &
            superstructure_size, new_element_position, derivative_structure)

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/current_derivative_term(1:2), p_tuple_extend(current_derivative_term(3), &
            p_tuple_getone(pert, 1))/), superstructure_size, new_element_position, &
            derivative_structure)

    else


       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%n_perturbations <= kn(2)) .AND.&
              current_derivative_term(2)%n_perturbations <= kn(2) ) .AND. &
              current_derivative_term(3)%n_perturbations <= kn(2) ) .eqv. .TRUE.) then

             new_element_position = new_element_position + 1
             derivative_structure(new_element_position, :) = current_derivative_term
 
          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%n_perturbations,  &
                          current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%n_perturbations,  &
                          current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%n_perturbations, &
                          current_derivative_term(3)%pid, kn) ) .eqv. .FALSE.) then

             new_element_position = new_element_position + 1 
             derivative_structure(new_element_position, :) = current_derivative_term(:)

          end if

       end if

    end if

  end subroutine


  function get_fds_data_index(pert_tuple, total_num_perturbations, which_index_is_pid, &
                              indices_len, indices)

    implicit none

    type(p_tuple) :: pert_tuple
    integer :: i, total_num_perturbations, indices_len
    integer, allocatable, dimension(:) :: get_fds_data_index
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: indices

    if (pert_tuple%n_perturbations == 0) then

       allocate(get_fds_data_index(1))
       get_fds_data_index(1) = 1

    else

       allocate(get_fds_data_index(pert_tuple%n_perturbations))

       do i = 1, pert_tuple%n_perturbations

          get_fds_data_index(i) = indices(which_index_is_pid(pert_tuple%pid(i)))

       end do

    end if

  end function


  function frequency_zero_or_sum(pert_tuple)

    implicit none

    type(p_tuple) :: pert_tuple
    complex(8) :: frequency_zero_or_sum
    integer :: i

    frequency_zero_or_sum = 0.0

    if (pert_tuple%n_perturbations > 0) then

       do i = 1, pert_tuple%n_perturbations

          frequency_zero_or_sum = frequency_zero_or_sum + pert_tuple%freq(i)

       end do

    end if

  end function


  subroutine mat_manual_attribute_inherit(A, B)

    implicit none

    type(matrix) :: A, B

    A%complex = B%complex
    A%closed_shell = B%closed_shell
    A%open_shell = B%open_shell
    A%ih_sym = B%ih_sym
    A%pg_sym = B%pg_sym

  end subroutine


  function rsp_get_matrix_w(mol, superstructure_size, &
           deriv_struct, total_num_perturbations, which_index_is_pid, &
           indices_len, ind, F, D, S) result(W)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: W, A, B, C

    call mat_nullify(W)
    W%nrow = mol%zeromat%nrow
    W%ncol = mol%zeromat%ncol
    W%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
    W%magic_tag = 825169837 !mark as set-up
    call mat_alloc(W)
    call mat_axpy((0.0d0, 0.0d0), W, .false., .true., W)

    A = tiny(0.0d0)*mol%zeromat
    call mat_alloc(A)

    B = tiny(0.0d0)*mol%zeromat
    call mat_alloc(B)

    C = tiny(0.0d0)*mol%zeromat
    call mat_alloc(C)

    do i = 1, superstructure_size

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(F, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + A * B * C

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) *  &
               A * B * C

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + ((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3))  * &
               A * B * C

    end do

    A = 0
    B = 0
    C = 0

  end function


  function rsp_get_matrix_y(mol, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S) result(Y)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Y, A, B, C

    call mat_manual_attribute_inherit(Y, mol%zeromat)
    call mat_manual_attribute_inherit(A, mol%zeromat)
    call mat_manual_attribute_inherit(B, mol%zeromat)
    call mat_manual_attribute_inherit(C, mol%zeromat)

    call mat_nullify(Y)
    Y%nrow = mol%zeromat%nrow
    Y%ncol = mol%zeromat%ncol
    Y%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
    Y%magic_tag = 825169837 !mark as set-up
    call mat_alloc(Y)
    call mat_axpy((0.0d0, 0.0d0), Y, .false., .true., Y)

    A = mol%zeromat
    call mat_alloc(A)
    B = mol%zeromat
    call mat_alloc(B)
    C = mol%zeromat
    call mat_alloc(C)

    ! NOTE (MaR): THIS CAN PROBABLY BE DONE MORE ECONOMICALLY

    Y%elms= 0.0

    do i = 1, superstructure_size

       call sdf_getdata_s(F, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y + A*B*C

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(F, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y - A*B*C

       ! Note (MaR): Can do frequency_zero_or_sum in 
       ! if evaluation to save matrix fetching
        
       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y - ((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3)) * A*B*C
        
       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y - ((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) * A*B*C

    end do

    A = 0
    B = 0
    C = 0

  end function


  function rsp_get_matrix_z(mol, superstructure_size, deriv_struct, kn, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S) result(Z)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    type(p_tuple) :: merged_p_tuple
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Z, A, B, C

    call mat_manual_attribute_inherit(Z, mol%zeromat)
    call mat_manual_attribute_inherit(A, mol%zeromat)
    call mat_manual_attribute_inherit(B, mol%zeromat)
    call mat_manual_attribute_inherit(C, mol%zeromat)

    call mat_nullify(Z)
    Z%nrow = mol%zeromat%nrow
    Z%ncol = mol%zeromat%ncol
    Z%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
    Z%magic_tag = 825169837 !mark as set-up
    call mat_alloc(Z)
    call mat_axpy((0.0d0, 0.0d0), Z, .false., .true., Z)

    A = mol%zeromat
    call mat_alloc(A)
    B = mol%zeromat
    call mat_alloc(B)
    C = mol%zeromat
    call mat_alloc(C)

    Z%elms = 0.0

    do i = 1, superstructure_size

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Z = Z + A*B*C

!     write(*,*) 'i is',i
!     write(*,*) 'Z', Z%elms
!     write(*,*) 'A', A%elms
!     write(*,*) 'B', B%elms
!     write(*,*) 'C', C%elms
    end do

    merged_p_tuple = merge_p_tuple(deriv_struct(1,1), merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3)))

    if (kn_skip(total_num_perturbations, merged_p_tuple%pid, kn) .eqv. .FALSE.) then

        call sdf_getdata_s(D, merged_p_tuple, get_fds_data_index(merged_p_tuple, &
        total_num_perturbations, which_index_is_pid, indices_len, ind), A)

        Z = Z - A

    end if


    A = 0
    B = 0
    C = 0

  end function


  function rsp_get_matrix_lambda(mol, p_tuple_a, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, ind, D, S) result(L)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple) :: p_tuple_a, merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: D, S
    type(matrix) :: L, A, B, C

    call mat_manual_attribute_inherit(L, mol%zeromat)
    call mat_manual_attribute_inherit(A, mol%zeromat)
    call mat_manual_attribute_inherit(B, mol%zeromat)
    call mat_manual_attribute_inherit(C, mol%zeromat)

    call mat_nullify(L)
    L%nrow = mol%zeromat%nrow
    L%ncol = mol%zeromat%ncol
    L%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
    L%magic_tag = 825169837 !mark as set-up
    call mat_alloc(L)
    call mat_axpy((0.0d0, 0.0d0), L, .false., .true., L)

    A = mol%zeromat
    call mat_alloc(A)
    B = mol%zeromat
    call mat_alloc(B)
    C = mol%zeromat
    call mat_alloc(C)

    L%elms = 0.0

    do i = 1, superstructure_size

       merged_A = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       call sdf_getdata_s(D, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       L = L + A * B * C
        
       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       L = L - A * B * C
       
    end do

    A = 0
    B = 0
    C = 0

  end function


  function rsp_get_matrix_zeta(mol, p_tuple_a, kn, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S) result(Zeta)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple) :: p_tuple_a, merged_p_tuple, merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Zeta, A, B, C

    call mat_manual_attribute_inherit(Zeta, mol%zeromat)
    call mat_manual_attribute_inherit(A, mol%zeromat)
    call mat_manual_attribute_inherit(B, mol%zeromat)
    call mat_manual_attribute_inherit(C, mol%zeromat)

    call mat_nullify(Zeta)
    Zeta%nrow = mol%zeromat%nrow
    Zeta%ncol = mol%zeromat%ncol
    Zeta%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
    Zeta%magic_tag = 825169837 !mark as set-up
    call mat_alloc(Zeta)
    call mat_axpy((0.0d0, 0.0d0), Zeta, .false., .true., Zeta)

    A = mol%zeromat
    call mat_alloc(A)
    B = mol%zeromat
    call mat_alloc(B)
    C = mol%zeromat
    call mat_alloc(C)

    do i = 1, superstructure_size

       merged_A = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       call sdf_getdata_s(F, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + A * B * C

       call sdf_getdata_s(F, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta - A * B * C

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,1)) * A * B * C
         
       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(F, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta +  A * B * C
         
       call sdf_getdata_s(S, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(F, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta - A * B * C
         
       call sdf_getdata_s(S, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,3)) * A * B * C
         
       call sdf_getdata_s(S, merged_A, get_fds_data_index(merged_A,  &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C
         
    end do

    merged_p_tuple = merge_p_tuple(p_tuple_a, merge_p_tuple(deriv_struct(1,1), &
                     merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3))))

    if (kn_skip(merged_p_tuple%n_perturbations, &
        merged_p_tuple%pid, kn) .eqv. .FALSE.) then

       call sdf_getdata_s(F, merged_p_tuple, get_fds_data_index(merged_p_tuple, & 
       total_num_perturbations, which_index_is_pid, indices_len, ind), A)

       Zeta = Zeta - A

    end if

    A = 0
    B = 0
    C = 0

  end function


  subroutine get_pulay_kn(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, sstr_incr, offset
    integer :: property_size, d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = derivative_superstructure_getsize(mol, p12(2), kn, .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structb(d_supsize, 3))

    sstr_incr = 0

    call derivative_superstructure(mol, p12(2), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, sstr_incr, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations))
    allocate(tmp(product(p12(1)%pdim)))
    allocate(inner_offsets(product(p12(1)%pdim)))
    allocate(outer_indices(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(inner_indices(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_onlysmall(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                      p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, inner_indices)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices)

    do i = 1, size(outer_indices, 1)

       tmp = 0.0

       W = mol%zeromat
       call mat_alloc(W)

       W = rsp_get_matrix_w(mol, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                            p12(2)%n_perturbations, which_index_is_pid, &
                            p12(2)%n_perturbations, outer_indices(i,:), F, D, S)

! write(*,*) 'got W', W%elms
! write(*,*) 'ncinner', ncinner

       call rsp_ovlave(mol, p12(1)%n_perturbations, p12(1)%plab, &
                      (/ (j/j, j = 1, p12(1)%n_perturbations) /), p12(1)%pdim, W, tmp)

! write(*,*) 'tmp tensor', tmp

       do j = 1, size(inner_indices, 1)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/inner_indices(j,:), &
                   outer_indices(i,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

! write(*,*) 'indices', (/inner_indices(j,:), outer_indices(i,:) /)
! write(*,*) 'offset', offset
! write(*,*) 'tmp(j)', tmp(j)

          prop(offset) = prop(offset) + tmp(j)
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do

!     write(*,*) 'pulay kn contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)


    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)    

    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    W = 0

  end subroutine


  recursive subroutine rsp_pulay_kn(mol, pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_kn(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)

       call rsp_pulay_kn(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(2)%n_perturbations, p12(2)%pid, kn) .EQV. .FALSE.) then


          write(*,*) 'Getting Pulay k-n contribution:'
          write(*,*) 'S', p12(1)%pid
          write(*,*) 'W', p12(2)%pid

          open(unit=257, file='totterms', status='old', action='write', &
               position='append')
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

             call get_pulay_kn(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2)  /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay k-n contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Pulay k-n contribution was k-n skipped:'
          write(*,*) 'S ', p12(1)%pid 
          write(*,*) 'W ', p12(2)%pid 
          write(*,*) ' '

       end if 

    end if

  end subroutine


  subroutine get_pulaylag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, k ,m, incr, offset
    integer :: property_size, d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:) :: outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = derivative_superstructure_getsize(mol, p12(2), kn, .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
    allocate(deriv_structb(d_supsize, 3))

    incr = 0

    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, incr, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(tmp(product(p12(1)%pdim)))
    allocate(inner_offsets(product(p12(1)%pdim)))
    allocate(outer_indices(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(inner_indices(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
              p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices)
    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, inner_indices)

    W = mol%zeromat
    call mat_alloc(W)

    do i = 1, size(outer_indices, 1)

       tmp = 0.0

       W = rsp_get_matrix_w(mol, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                            p12(2)%n_perturbations, which_index_is_pid, &
                            p12(2)%n_perturbations, outer_indices(i,:), F, D, S)

       call rsp_ovlave(mol, p12(1)%n_perturbations, p12(1)%plab, &
                       (/ (j/j, j = 1, p12(1)%n_perturbations) /), &
                       p12(1)%pdim, W, tmp)

       do j = 1, size(inner_indices, 1)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/inner_indices(j,:), &
                   outer_indices(i,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) + tmp(j)
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do
!     write(*,*) 'pulay lag contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)

    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_b_large)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    W = 0

  end subroutine


  recursive subroutine rsp_pulay_lag(mol, pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
       S, D, F, property_size, cache, prop)
       call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
       S, D, F, property_size, cache, prop)

    else

       ! At lowest level:
       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

       write(*,*) 'Getting Pulay lagrange contribution:'
       write(*,*) 'S', p12(1)%pid
       write(*,*) 'W', p12(2)%pid, 'primed', kn(2)

       open(unit=257, file='totterms', status='old', action='write', position='append') 
       write(257,*) 'T'
       close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '
       
             open(unit=257, file='cachehit', status='old', action='write', &
             position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)

          else

             call get_pulaylag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Pulay lagrange contribution was k-n skipped:'
          write(*,*) 'S', p12(1)%pid 
          write(*,*) 'W', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_idem_lag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: Zeta, Z
    integer :: i, j, k, m, n, p, incr1, incr2
    integer :: property_size, offset
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, ncprod, which_index_is_pid1, &
                                          which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(mol, p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(mol, p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(mol, p_tuple_remove_first(p12(1)), kn, .FALSE., & 
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(2), incr2, deriv_structb)


    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncprod(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                      p12(1)%n_perturbations, 1, p12(1), ncarray)

    do i = 1, size(ncarray)

       ncprod(i) = product(ncarray(i:size(ncarray)))/ncarray(i)

    end do

    which_index_is_pid1 = 0

    do i = 1, p12(1)%n_perturbations

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, outer_indices_a)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices_b)

    offset = 0.0

    do i = 1, size(outer_indices_a, 1)

       Zeta = rsp_get_matrix_zeta(mol, p_tuple_getone(p12(1), 1), kn, d_supsize(1), &
           deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
           which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), F, D, S)

       do j = 1, size(outer_indices_b, 1)

          Z = rsp_get_matrix_z(mol, d_supsize(2), deriv_structb, kn, &
              p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
              p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/outer_indices_a(i,:), &
                   outer_indices_b(j,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) -tr(Zeta, Z)
          prop_forcache(offset) = prop_forcache(offset) -tr(Zeta, Z)

       end do

    end do

!     write(*,*) 'idempotency contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache) 

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(ncprod)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    Zeta = 0
    Z = 0

  end subroutine


  recursive subroutine rsp_idem_lag(mol, pert, kn, p12, S, D, F, &
                                    property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)
       call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

          write(*,*) 'Getting idempotency lagrange contribution'
          write(*,*) 'Zeta', p12(1)%pid
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append')
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
      
          else

             ! At lowest level:
             call get_idem_lag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated idempotency lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Idempotency lagrange contribution was k-n skipped:'
          write(*,*) 'Zeta', p12(1)%pid 
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_scfe_lag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: L, Y
    integer :: i, j, k, m, n, p, incr1, incr2
    integer :: property_size, offset
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, ncprod, which_index_is_pid1, which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0
    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(mol, p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(mol, p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(mol, p_tuple_remove_first(p12(1)), kn, .FALSE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), & 
                    d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                    d_supsize(2), incr2, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
              p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid1 = 0

    do i = 1, p12(1)%n_perturbations

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, outer_indices_a)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices_b)

    offset = 0

    do i = 1, size(outer_indices_a, 1)

       L = rsp_get_matrix_lambda(mol, p_tuple_getone(p12(1), 1), d_supsize(1), &
           deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
           which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), D, S)

       do j = 1, size(outer_indices_b, 1)

          Y = rsp_get_matrix_y(mol, d_supsize(2), deriv_structb, &
              p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
              p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/outer_indices_a(i,:), &
          outer_indices_b(j,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) - tr(L, Y)
          prop_forcache(offset) = prop_forcache(offset) - tr(L, Y)

       end do

    end do

!     write(*,*) 'scfe contribution'
!  
! call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)

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
    L = 0
    Y = 0

  end subroutine


  recursive subroutine rsp_scfe_lag(mol, pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
            (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
            S, D, F, property_size, cache, prop)
       call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
            (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
            S, D, F, property_size, cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

          write(*,*) 'Getting scfe lagrange contribution'
          write(*,*) 'Lambda', p12(1)%pid
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

             ! At lowest level:
             call get_scfe_lag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), &
             kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated scfe lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'scfe lagrange contribution was k-n skipped:'
          write(*,*) 'Lambda', p12(1)%pid 
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine make_p_tuple_subset(pert, psub)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer :: i, j, k, m

    do i = 1, pert%n_perturbations

       psub(i)%n_perturbations = pert%n_perturbations - 1

       allocate(psub(i)%pdim(pert%n_perturbations - 1))
       allocate(psub(i)%plab(pert%n_perturbations - 1))
       allocate(psub(i)%pid(pert%n_perturbations - 1))
       allocate(psub(i)%freq(pert%n_perturbations - 1))

       k = i
       m = 1

       do j = 1, (pert%n_perturbations)

          if (.NOT.(j == k)) then

             psub(i)%pdim(m) = pert%pdim(j)
             psub(i)%plab(m) = pert%plab(j)
             psub(i)%pid(m) = pert%pid(j)
             psub(i)%freq(m) = pert%freq(j)

             m = m + 1

          end if

       end do
    end do
  
  end subroutine


  subroutine get_fock_lowerorder(mol, num_p_tuples, total_num_perturbations, p_tuples, &
                                 density_order, D, property_size, Fp)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF) :: D
    type(property_cache) :: energy_cache
    type(matrix), allocatable, dimension(:) :: dens_tuple
    integer :: i, j, k, m, num_p_tuples, total_num_perturbations, &
               density_order, property_size, offset
    integer, dimension(0) :: noc
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, &
                                          pidoutersmall, ncinnersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    type(matrix), allocatable, dimension(:) :: tmp
    type(matrix), dimension(property_size) :: Fp
    
    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
    ncouter = nc_only(total_num_perturbations, total_num_perturbations - & 
                      p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                      p_tuples(2:num_p_tuples), ncarray)
    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
                      p_tuples(1), ncarray)

    allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))
    allocate(ncinnersmall(p_tuples(1)%n_perturbations))
    allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))

    ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
                                p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                                p_tuples(2:num_p_tuples), ncarray)
    ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%n_perturbations, &
                   1, p_tuples(1), ncarray)
    pidoutersmall = get_pidoutersmall(total_num_perturbations - &
                    p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                    p_tuples(2:num_p_tuples))

    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(dens_tuple(num_p_tuples))
    allocate(outer_indices(product(ncoutersmall),size(ncoutersmall)))
    allocate(inner_indices(product(ncinnersmall),size(ncinnersmall)))
    allocate(tmp(product(ncinnersmall)))

    o_whichpert = make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)

    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
                      p_tuples(1)%n_perturbations, pidoutersmall, &
                      ncarray, ncoutersmall, o_whichpert)

    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

       do i = 1, size(o_whichpert)

          if (.NOT.(o_whichpert(i) == 0)) then

             o_wh_forave(o_whichpert(i)) = i

          end if
  
       end do


       do i = 1, num_p_tuples

          call mat_nullify(dens_tuple(i))
          dens_tuple(i)%nrow = mol%zeromat%nrow
          dens_tuple(i)%ncol = mol%zeromat%ncol
          dens_tuple(i)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
          dens_tuple(i)%magic_tag = 825169837 !mark as set-up
          call mat_alloc(dens_tuple(i))
          call mat_axpy((0.0d0, 0.0d0), dens_tuple(i), .false., .true., dens_tuple(i))

       end do

       do j = 1, size(tmp)

          call mat_nullify(tmp(j))
          tmp(j)%nrow = mol%zeromat%nrow
          tmp(j)%ncol = mol%zeromat%ncol
          tmp(j)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
          tmp(j)%magic_tag = 825169837 !mark as set-up
          call mat_alloc(tmp(j))
          call mat_axpy((0.0d0, 0.0d0), tmp(j), .false., .true., tmp(j))

       end do

       call make_indices(total_num_perturbations - p_tuples(1)%n_perturbations, 1, & 
                         ncoutersmall, 0, outer_indices)

       if (p_tuples(1)%n_perturbations > 0) then

          call make_indices(p_tuples(1)%n_perturbations, 1, p_tuples(1)%pdim, 0, &
          inner_indices)

       end if

       allocate(inner_offsets(product(ncinner)))

       do i = 1, size(outer_indices, 1)

          do j = 2, num_p_tuples

             call sdf_getdata_s(D, p_tuples(j), (/ &
                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
                             k = 1, p_tuples(j)%n_perturbations) /), dens_tuple(j))

          end do

          do j = 1, size(tmp)

             call mat_nullify(tmp(j))
             tmp(j)%nrow = mol%zeromat%nrow
             tmp(j)%ncol = mol%zeromat%ncol
             tmp(j)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
             tmp(j)%magic_tag = 825169837 !mark as set-up
             call mat_alloc(tmp(j))
             call mat_axpy((0.0d0, 0.0d0), tmp(j), .false., .true., tmp(j))

          end do


          if (num_p_tuples <= 1) then

             call rsp_oneint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, tmp)

          end if


          if (num_p_tuples <= 2) then

             call rsp_twoint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), tmp)

          end if

!              write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'
!
!              call rsp_xcint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                            (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                            p_tuples(1)%pdim, density_order, &
!                            (/ sdf_getdata(D, get_emptypert(), (/1/)), &
!                            (dens_tuple(k), k = 2, num_p_tuples) /), &
!                            tmp)

          if (p_tuples(1)%n_perturbations > 0) then

! write(*,*) 'np1 > 0'
             do j = 1, size(inner_indices, 1)

                offset = get_one_tensor_offset( sum( (/ (p_tuples(k)%n_perturbations, & 
                         k=1, num_p_tuples ) /) ), &
                         (/ inner_indices(j,:), outer_indices(i,:) /), &
                         (/ (p_tuples(k)%pid, k=1, num_p_tuples ) /), ncarray)
! write(*,*) 'indices', (/ inner_indices(j,:), outer_indices(i,:) /)
! write(*,*) 'offset', offset
! write(*,*) 'Fp tmp contribution', tmp(j)%elms
                Fp(offset) = Fp(offset) + tmp(j)

             end do

          else
! write(*,*) 'np1 = 0'
             offset = get_one_tensor_offset( sum( (/ (p_tuples(k)%n_perturbations, &
                      k=2, num_p_tuples ) /) ), &
                      (/ outer_indices(i,:) /), &
                      (/ (p_tuples(k)%pid, k=2, num_p_tuples ) /), ncarray)
! write(*,*) 'indices', outer_indices(i,:)
! write(*,*) 'offset', offset
! write(*,*) 'Fp tmp contribution', tmp(1)%elms
                Fp(offset) = Fp(offset) + tmp(1)

          end if

       end do

    deallocate(inner_offsets)

    else

       if (num_p_tuples <= 1) then

          call rsp_oneint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, &
                          Fp)

       end if

       if (num_p_tuples <= 2) then

          call rsp_twoint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
               (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
               p_tuples(1)%pdim, &
               sdf_getdata(D, get_emptypert(), (/1/) ), &
               Fp)

       end if

!        write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'

!        call rsp_xcint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                       (/product(p_tuples(1)%pdim)/), 1, &
!                       (/ sdf_getdata(D, get_emptypert(), (/1/)) /), &
!                       Fp)


    end if


    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(dens_tuple)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(tmp)

  end subroutine


  recursive subroutine fock_lowerorder(mol, pert, total_num_perturbations, num_p_tuples, p_tuples, &
                       density_order, D, property_size, Fp)

    implicit none

    logical :: density_order_skip
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF) :: D
    type(property_cache) :: energy_cache
    type(matrix), dimension(property_size) :: Fp

    if (pert%n_perturbations >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%n_perturbations == 0) then

          call fock_lowerorder(mol, p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               density_order, D, property_size, Fp)

       else

          call fock_lowerorder(mol, p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
               p_tuples(2:size(p_tuples))/), &
               density_order, D, property_size, Fp)

       end if
    
       ! 2. Differentiate all of the contraction densities in turn

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%n_perturbations == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call fock_lowerorder(mol, p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               t_new, density_order + 1, D, property_size, Fp)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call fock_lowerorder(mol, p_tuple_remove_first(pert), total_num_perturbations, &
            num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
            density_order + 1, D, property_size, Fp)

    else

       write(*,*) 'Getting perturbed Fock matrix lower order contribution'

       do i = 1, num_p_tuples
 
          if (i == 1) then

             write(*,*) 'F', p_tuples(i)%pid

          else

             write(*,*) 'D', p_tuples(i)%pid

          end if

       end do

       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%n_perturbations >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then

          call get_fock_lowerorder(mol, num_p_tuples, total_num_perturbations, &
                                   p_tuples, density_order, D, property_size, Fp)

          write(*,*) 'Calculated perturbed Fock matrix lower order contribution'
          write(*,*) ' '

       else

          write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
          write(*,*) 'at order for which perturbed D is to be found '
          write(*,*) ' '

       end if


    end if

  end subroutine


  function make_one_index_tuple(n_perturbations, pdim, icomp_in)

    implicit none

    integer :: n_perturbations, icomp, i, icomp_in
    integer, dimension(n_perturbations) :: pdim, make_one_index_tuple

    icomp = icomp_in - 1

    do i = 1, n_perturbations

       ! The integer division (rounding to nearest lower integer) should make this work
       make_one_index_tuple(i) = icomp/(product(pdim(i:n_perturbations))/pdim(i)) + 1
       icomp = icomp - (make_one_index_tuple(i) - 1) * &
                       (product(pdim(i:n_perturbations))/pdim(i))

    end do

  end function


  ! ASSUMES THAT PERTURBATION TUPLE IS IN STANDARD ORDER
  subroutine get_fds(mol, pert, F, D, S)

    use interface_rsp_solver, only: rsp_mosolver_exec

    implicit none

    type(rsp_cfg) :: mol
    integer :: sstr_incr, i, j, superstructure_size
    integer, allocatable, dimension(:) :: ind
    integer, dimension(0) :: noc
    character(4), dimension(0) :: nof
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(SDF) :: F, D, S
    type(matrix) :: X(1), RHS(1), A, B, TEST
    type(matrix), dimension(product(pert%pdim)) :: Fp, Dp, Sp, Dh

! write(*,*) 'Starting'

    A = mol%zeromat
    call mat_alloc(A)
! write(*,*) 'Starting'
    B = mol%zeromat
    call mat_alloc(B)

    call sdf_getdata_s(D, get_emptypert(), (/1/), A)
! write(*,*) 'Starting'
    call sdf_getdata_s(S, get_emptypert(), (/1/), B)

! write(*,*) 'Starting'
    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix

    call rsp_ovlint(mol, pert%n_perturbations, pert%plab, &
                    (/ (1, j = 1, pert%n_perturbations) /), pert%pdim, Sp)
    call sdf_add(S, pert, Sp)

! write(*,*) 'Got Sp'
! do i = 1, product(pert%pdim)
! 
! write(*,*) 'Sp at', i
!        write(*,*) Sp(i)%elms
! 
! end do


    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp

    do i = 1, product(pert%pdim)

       call mat_nullify(Dp(i))
          Dp(i)%nrow = mol%zeromat%nrow
          Dp(i)%ncol = mol%zeromat%ncol
          Dp(i)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
          Dp(i)%magic_tag = 825169837 !mark as set-up
          call mat_alloc(Dp(i))
          call mat_axpy((0.0d0, 0.0d0), Dp(i), .false., .true., Dp(i))
! write(*,*) 'Set Dp'
!        Dp(i) = mol%zeromat
!        call mat_alloc(Dp(i))
!        Dp(i)%elms = 0.0

       call mat_nullify(Dh(i))
          Dh(i)%nrow = mol%zeromat%nrow
          Dh(i)%ncol = mol%zeromat%ncol
          Dh(i)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
          Dh(i)%magic_tag = 825169837 !mark as set-up
          call mat_alloc(Dh(i))
          call mat_axpy((0.0d0, 0.0d0), Dh(i), .false., .true., Dh(i))
! write(*,*) 'Set Dh'

       call mat_nullify(Fp(i))
          Fp(i)%nrow = mol%zeromat%nrow
          Fp(i)%ncol = mol%zeromat%ncol
          Fp(i)%closed_shell = .true.       !*2 on tr(A,B) and dot(A,B)
          Fp(i)%magic_tag = 825169837 !mark as set-up
          call mat_alloc(Fp(i))
          call mat_axpy((0.0d0, 0.0d0), Fp(i), .false., .true., Fp(i))
!        Fp(i) = mol%zeromat
!        call mat_alloc(Fp(i))
!        Fp(i)%elms = 0.0
! write(*,*) 'Set Fp'
    end do

    call sdf_add(D, pert, Dp)

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

    call fock_lowerorder(mol, pert, pert%n_perturbations, 1, (/get_emptypert()/), &
                         0, D, product(pert%pdim), Fp)


    call sdf_add(F, pert, Fp)

! write(*,*) 'added Fp'
! 
!        write(*,*) Dp(1)%elms
!        write(*,*) Fp(1)%elms

    ! b) For Dp: Create differentiation superstructure: First dryrun, then actual call

    superstructure_size = derivative_superstructure_getsize(mol, pert, &
                          (/pert%n_perturbations, pert%n_perturbations/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

    allocate(derivative_structure(superstructure_size,3))
    allocate(ind(pert%n_perturbations))

    call derivative_superstructure(mol, pert, (/pert%n_perturbations, &
         pert%n_perturbations/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)

    do i = 1, product(pert%pdim)


! write(*,*) 'i', i
       ind = make_one_index_tuple(pert%n_perturbations, pert%pdim, i)

!         write(*,*) 'Dp before z', Dp(i)%elms

       Dp(i) = 1.0d0 * rsp_get_matrix_z(mol, superstructure_size, derivative_structure, &
               (/pert%n_perturbations,pert%n_perturbations/), pert%n_perturbations, &
               (/ (j, j = 1, pert%n_perturbations) /), pert%n_perturbations, &
               ind, F, D, S)

! write(*,*) 'got z'

!         write(*,*) 'Dp at z', Dp(i)%elms

       Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A


!        Dp(i) = Dp(i) - sdf_getdata(D, get_emptypert(), (/1/)) * &
!                        sdf_getdata(S, get_emptypert(), (/1/)) * Dp(i) - &
!                Dp(i) * sdf_getdata(S, get_emptypert(), (/1/)) * &
!                        sdf_getdata(D, get_emptypert(), (/1/)) 

!        write(*,*) 'Dp at projection', Dp(i)%elms

       call sdf_add(D, pert, Dp)

       ! 3. Complete the particular contribution to Fp
       ! NOTE (MaR): THERE SHOULD BE A CALL TO rsp_xcint HERE TOO
       ! THE IF CRITERION HERE NEEDS ANOTHER LOOK

!        if (pert%n_perturbations <=2) then
! 
          call rsp_twoint(mol, 0, nof, noc, pert%pdim, Dp(i), Fp(i:i))
! 
!        end if


       call sdf_add(F, pert, Fp)

! write(*,*) 'added Fp'
! 
!        write(*,*) Dp(i)%elms
!        write(*,*) Fp(i)%elms


       ! 4. Make right-hand side using Dp

       RHS(1) = mol%zeromat
       call mat_alloc(RHS(1))

       RHS(1) = rsp_get_matrix_y(mol, superstructure_size, derivative_structure, &
                pert%n_perturbations, (/ (j, j = 1, pert%n_perturbations) /), &
                pert%n_perturbations, ind, F, D, S)
     
       X(1) = 0*RHS(1)
       call mat_alloc(X(1))
       X(1)%elms= 0.0

! write(*,*) 'RHS', RHS(1)%elms
! write(*,*) 'set up solver'

       ! Note (MaR): What does the second argument in rsp_mosolver_exec mean?
       call rsp_mosolver_exec(RHS(1), (/0d0/), X)
       ! Note (MaR): Why multiply by -2 like below?
       X(1) = -2d0*X(1)
       RHS(1) = 0

! write(*,*) 'solved', X(1)%elms
       ! 5. Get Dh using the rsp equation solution X

       Dh(i) = X(1) * B * A - A * B * X(1)
! 
! sdf_getdata(S, get_emptypert(), (/1/)) * &
!                       sdf_getdata(D, get_emptypert(), (/1/)) - &
!                       sdf_getdata(D, get_emptypert(), (/1/)) * &
!                       sdf_getdata(S, get_emptypert(), (/1/)) * X(1)

       ! 6. Make homogeneous contribution to Fock matrix

! write(*,*) 'got Dh'

!        if (pert%n_perturbations <=2) then

          call rsp_twoint(mol, 0, nof, noc, pert%pdim, Dh(i), Fp(i:i))

!        end if

! write(*,*) 'added last'

       ! 'NOTE (MaR): XCINT CALL SKIPPED FOR NOW'

       ! call rsp_xcint(mol, 0, nof, noc, pert%pdim, 2, &
       ! (/sdf_getdata(D, get_emptypert(), (/1/)), Dh(i)/), Fp(i:i))

       ! 7. Complete perturbed D with homogeneous part

       Dp(i) = Dp(i) + Dh(i)

! write(*,*) 'Finally, Dp(i) at i', i, 'is'
! write(*,*) Dp(i)%elms
! write(*,*) 'Fp(i) at i', i, 'is'
! write(*,*) Fp(i)%elms


! write(*,*) 'added homogeneous'

    end do

    deallocate(derivative_structure)
    deallocate(ind)

    ! Add the final values to cache

! write(*,*) 'adding final values to cache'
    call sdf_add(F, pert, Fp)
    call sdf_add(D, pert, Dp)
! write(*,*) 'added final values to cache'

! NOTE (MaR): FOR SOME REASON, IF I DON'T PRINT HERE THE VALUES ARE ZEROED LATER

! if (pert%pdim(1) == 12) then

! TEST = mol%zeromat
! call mat_alloc(TEST)
! do i = 1, pert%pdim(1)
! 
!           call sdf_getdata_s(D, pert, (/ i /), TEST)
! 	  write(*,*) 'dt2 3 in get fds', TEST%elms
! 
! 
! 
! end do
! end if
  end subroutine


  recursive subroutine rsp_fds(mol, pert, kn, F, D, S)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer, dimension(2) :: kn
    integer :: i, j, k
    type(SDF) :: F, D, S

    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then (at final recursion level) get perturbed F, D, S 

    if (pert%n_perturbations > 1) then

       call make_p_tuple_subset(pert, psub)

       do i = 1, size(psub)

          call rsp_fds(mol, psub(i), kn, F, D, S)

       end do       

    end if

    if (sdf_already(D, pert) .eqv. .FALSE.) then
         
       if (kn_skip(pert%n_perturbations, pert%pid, kn) .eqv. .FALSE.) then

          write(*,*) 'Calling ovlint/fock/density with labels ', pert%plab, &
                     ' and perturbation id ', pert%pid
          write(*,*) ' '
                 
          ! FIXME (MaR) Quick fix: Reenumerate pids from 1 and up so that 
          ! get_fds doesn't stumble. It seems to work, but consider rewrite.
         
          k = 1

          do j = 1, pert%n_perturbations

             pert%pid(j) = k
             k = k + 1

          end do

          call get_fds(mol, p_tuple_standardorder(pert), F, D, S)

       else

          write(*,*) 'Would have called ovlint/fock/density with labels ', &
                     pert%plab, ' and perturbation id ', pert%pid, &
                     ' but it was k-n forbidden'
          write(*,*) ' '

       end if

    else

       write(*,*) 'FDS for labels ', pert%plab, &
                  'and perturbation id ', pert%pid, ' was found in cache'
       write(*,*) ' '

    end if

  end subroutine

  function rm_dim(np, dims, skip)

    implicit none

    integer :: skip, i, j, np
    integer, dimension(np) :: dims
    integer, dimension(np - 1) :: rm_dim

    j = 1

    do i = 1, np

       if (i == skip) then

       else

          rm_dim(j) = dims(i)
  
          j = j + 1   

       end if

    end do

  end function


  recursive function kn_prod(np, lvl, dims, k_or_n, sofar) result (kn_p)

    implicit none

    integer :: np, lvl, k_or_n, sofar, kn_p, i
    integer, dimension(np) :: dims

    kn_p = 0

    if (lvl < k_or_n) then

       do i = 1, np

          kn_p = kn_p + kn_prod(np - 1, lvl + 1, &
                                rm_dim(np, dims, i), k_or_n, dims(i)*sofar)

       end do

    else

       kn_p = sofar

    end if

  end function


  function get_bestkn(pert)

    implicit none

    type(p_tuple) :: pert
    integer, dimension(2) :: get_bestkn
    integer :: min_n, n, csize, minsize

    ! Do the case n = 1 first to establish one value

    minsize = kn_prod(pert%n_perturbations - 1, 1, pert%pdim(2:pert%n_perturbations), &
              pert%n_perturbations - 1 - 1, pert%pdim(1)) 


    minsize = minsize + kn_prod(pert%n_perturbations - 1, &
              0, pert%pdim(2:pert%n_perturbations), 1, 1)

    min_n = 1

    ! Get the products for the pert dimensions as dictated by k and n
    ! Assume that integer division rounds down to nearest lower integer

    do n = (pert%n_perturbations/2), pert%n_perturbations - 1

       csize = kn_prod(pert%n_perturbations - 1, 1, pert%pdim(2:pert%n_perturbations), &
                 pert%n_perturbations - n - 1, pert%pdim(1))


       csize = csize + kn_prod(pert%n_perturbations - 1, &
                 0, pert%pdim(2:pert%n_perturbations), n, 1)

       ! If the products are smaller than the previous min, update 
       ! index identifer and minsize

       if (csize < minsize) then

          min_n = n
          minsize = csize

       end if

    end do

    get_bestkn(2) = min_n
    get_bestkn(1) = pert%n_perturbations - min_n - 1

  end function


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


  subroutine get_prop(mol, pert, kn, prop, F, D, S)

    implicit none

    type(matrix) :: TEST
    type(SDF) :: F, D, S
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: emptyp_tuples
    integer, dimension(2) :: kn
    complex(8), dimension(product(pert%pdim)) :: prop
    type(property_cache), pointer :: energy_cache, pulay_kn_cache, &
                       pulay_lag_cache, idem_cache, scfe_cache

    emptypert%n_perturbations = 0
    allocate(emptypert%pdim(0))    
    allocate(emptypert%plab(0))
    allocate(emptypert%pid(0))
    allocate(emptypert%freq(0))

    emptyp_tuples = (/emptypert, emptypert/)

    ! Get all necessary F, D, S derivatives as dictated by
    ! number of perturbations and kn
    ! Maybe also extend to get W for kn part of rsp_pulay

    call rsp_fds(mol, pert, kn, F, D, S)
 
    write(*,*) ' '
    write(*,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
    write(*,*) ' '
! TEST = mol%zeromat
! call mat_alloc(TEST)
! 
!           call sdf_getdata_s(D, p_tuple_getone(pert,2), (/ 3 /), TEST)
! 	  write(*,*) 'dt2 3', TEST%elms
!           call sdf_getdata_s(D, p_tuple_getone(pert,2), (/ 5 /), TEST)
! 	  write(*,*) 'dt2 5', TEST%elms



    call property_cache_allocate(energy_cache)
    call rsp_ener(mol, pert, pert%n_perturbations, kn, 1, (/emptypert/), 0, D, &
                  product(pert%pdim), energy_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating energy-type contributions'
    write(*,*) ' '

    deallocate(energy_cache)

    call property_cache_allocate(pulay_kn_cache)
    call rsp_pulay_kn(mol, pert, kn, (/emptypert, emptypert/), S, D, F, &
                      product(pert%pdim), pulay_kn_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating Pulay k-n type contributions'
    write(*,*) ' '

    deallocate(pulay_kn_cache)

    call property_cache_allocate(pulay_lag_cache)
    call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
                       (/p_tuple_getone(pert,1), emptypert/), &
                       S, D, F, product(pert%pdim), pulay_lag_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating Pulay lagrangian type contributions' 
    write(*,*) ' '

    deallocate(pulay_lag_cache)

    call property_cache_allocate(idem_cache)
    call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
                      (/p_tuple_getone(pert,1), emptypert/), &
                      S, D, F, product(pert%pdim), idem_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating idempotency lagrangian type contributions'
    write(*,*) ' '

    deallocate(idem_cache)

    call property_cache_allocate(scfe_cache)
    call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
                      (/p_tuple_getone(pert,1), emptypert/), &
                      S, D, F, product(pert%pdim), scfe_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating SCF lagrangian type contributions'
    write(*,*) ' '

    deallocate(scfe_cache)

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

  subroutine rsp_prop(mol, pert_unordered, kn, F_unperturbed, D_unperturbed, S_unperturbed)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, pert_unordered
    type(matrix) :: F_unperturbed, D_unperturbed, S_unperturbed
    type(SDF), pointer :: F, D, S
    integer, dimension(2) :: kn
    integer :: i, j
    complex(8), allocatable, dimension(:) :: prop

    open(unit=257, file='totterms', status='replace', action='write') 
    write(257,*) 'START'
    close(257)

    open(unit=257, file='cachehit', status='replace', action='write') 
    write(257,*) 'START'
    close(257)


    pert = p_tuple_standardorder(pert_unordered)


! The get_bestkn function is taken out of use for now until it has been improved
!    kn = get_bestkn(pert)
   write(*,*) ' '
   write(*,*) 'Choice of k, n is ', kn(1), ' and ', kn(2)
   write(*,*) ' '

! kn(1) = 0
! kn(2) = 2
! 
! write(*,*) 'NOTE: (k, n) forced for debug purposes:', kn(1), kn(2)

    allocate(S)
    S%next => S
    S%last = .TRUE.
    S%perturb%n_perturbations = 0
    allocate(S%perturb%pdim(0))
    allocate(S%perturb%plab(0))
    allocate(S%perturb%pid(0))
    allocate(S%perturb%freq(0))
    allocate(S%data(1))
    S%data = S_unperturbed

    allocate(D)
    D%next => D
    D%last = .TRUE.
    D%perturb%n_perturbations = 0
    allocate(D%perturb%pdim(0))
    allocate(D%perturb%plab(0))
    allocate(D%perturb%pid(0))
    allocate(D%perturb%freq(0))
    allocate(D%data(1))
    D%data = D_unperturbed

    allocate(F)
    F%next => F
    F%last = .TRUE.
    F%perturb%n_perturbations = 0
    allocate(F%perturb%pdim(0))
    allocate(F%perturb%plab(0))
    allocate(F%perturb%pid(0))
    allocate(F%perturb%freq(0))
    allocate(F%data(1))
    F%data = F_unperturbed

    zromtrx = mol%zeromat
    call mat_alloc(zromtrx)

    allocate(prop(product(pert%pdim)))
    prop = 0.0

    call get_prop(mol, pert, kn, prop, F, D, S)

    write(*,*) 'Property was calculated and printed to rsp_tensor'
    write(*,*) ' '
    open(unit=260, file='rsp_tensor', status='replace', action='write') 
    write(260,*) ' '
    close(260)
    call print_rsp_tensor(size(pert%pdim),size(pert%pdim),pert%pdim, prop, 1)


    write(*,*) 'End of print'

    open(unit=257, file='totterms', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

    open(unit=257, file='cachehit', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

  end subroutine

end module
