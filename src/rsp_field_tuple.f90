! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains the definition and code for manipulation 
! of the field tuple datatype.

module rsp_field_tuple

  implicit none

  public p_tuple_extend
  public p_tuple_getone
  public p_tuple_remove_first
  public merge_p_tuple
  public p_tuple_p1_cloneto_p2
  public p_tuple_deallocate
  public get_emptypert
  public p_tuple_p1_lt_p2
  public p_tuple_standardorder
  public p_tuples_standardorder
  public p_tuples_compare
  public plab_compare
  public pfreq_compare
  public p_tuple_compare
  public make_p_tuple_subset

  ! Define perturbation tuple datatype

  type p_tuple

     integer :: n_perturbations ! Number of perturbations
     integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
     character(4), allocatable, dimension(:) :: plab ! Perturbation labels
     integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
     complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
     ! Add other perturbation identification info as needed

  end type

  contains

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

end module
