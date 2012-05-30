! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and associated functions/definitions for the
! creation of indices and addressing of elements of the (collapsed)
! arrays used throughout the rsp_general calculation

module rsp_indices_and_addressing

  use rsp_field_tuple

  implicit none

  public get_one_tensor_offset
  public test_making_triangulated_indices
  public get_num_blks
  public get_blk_info
  public make_triangulated_indices
  public index_blks_direct_product
  public make_one_triang_index_blk
  public get_triangulated_size
  public get_one_triangular_size
  public fact_terminate_lower
  public make_one_index_tuple

  ! Define triangulated index block datatype

  type triangulated_index_block

     integer, allocatable, dimension(:,:) :: t_ind

  end type

  contains

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


  subroutine test_making_triangulated_indices(fields)

    implicit none

    type(p_tuple) :: fields
    integer :: i, num_blks, triangulated_size
    integer, allocatable, dimension(:,:) :: blk_info, indices


    num_blks = get_num_blks(fields)

    write(*,*) 'num_blks is', num_blks    

    allocate(blk_info(num_blks, 3))

    blk_info = get_blk_info(num_blks, fields)

    write(*,*) 'blk_info is', blk_info

    triangulated_size = get_triangulated_size(num_blks, blk_info)

    write(*,*) 'triang_size is', triangulated_size

    allocate(indices(triangulated_size, sum(blk_info(:,2))))

    call make_triangulated_indices(num_blks, blk_info, triangulated_size, indices)

    do i = 1, size(indices, 1)

       write(*,*) 'indices at i =', i, ' :', indices(i,:)
       
    end do

    deallocate(blk_info)
    deallocate(indices)

  end subroutine




 ! Assumes fields is sorted
  function get_num_blks(fields)

    implicit none

    integer :: i, curr_blk_start, get_num_blks
    type(p_tuple) :: fields
    type(p_tuple), allocatable, dimension(:) :: each_field

    allocate(each_field(fields%n_perturbations))

    do i = 1, fields%n_perturbations

       each_field(i) = p_tuple_getone(fields, i)

    end do

    curr_blk_start = 1
    get_num_blks = 1

    do i = 1, fields%n_perturbations

       if (p_tuple_p1_lt_p2(each_field(curr_blk_start),each_field(i)) .EQV. .TRUE.) then

          curr_blk_start = i
          get_num_blks = get_num_blks + 1

       end if

    end do

    deallocate(each_field)

  end function

! Return array of size (nblks, 3) where row 1 is curr_blk_start, row 2 is block_len
! row 3 is pdim

  function get_blk_info(nblks, fields)

    implicit none

    integer :: nblks, i, this_blk, curr_blk_start
    integer, dimension(nblks,3) :: get_blk_info   
    type(p_tuple) :: fields
    type(p_tuple), allocatable, dimension(:) :: each_field

    allocate(each_field(fields%n_perturbations))

    if (fields%n_perturbations > 0) then

       do i = 1, fields%n_perturbations

          each_field(i) = p_tuple_getone(fields, i)

       end do

       this_blk = 1
       curr_blk_start = 1
       get_blk_info(1, 1) = 1
       get_blk_info(1, 3) = each_field(1)%pdim(1)


       do i = 1, fields%n_perturbations

          if (p_tuple_p1_lt_p2(each_field(curr_blk_start),each_field(i)) &
             .EQV. .TRUE.) then

             curr_blk_start = i
             this_blk = this_blk + 1
             get_blk_info(this_blk, 1) = i
             get_blk_info(this_blk - 1, 2) = i - get_blk_info(this_blk - 1, 1)
             get_blk_info(this_blk, 3) = each_field(i)%pdim(1)

          end if

       end do

       if (nblks > 1) then

          get_blk_info(nblks, 2) = fields%n_perturbations - get_blk_info(nblks, 1) + 1

       elseif (nblks == 1) then

          get_blk_info(1, 2) = fields%n_perturbations

       end if

    else

       get_blk_info(1, 1) = 1
       get_blk_info(1, 2) = 1
       get_blk_info(1, 3) = 1

    end if

    deallocate(each_field)

  end function

  function get_triang_blks_tuple_offset(ntuple, nblks_tuple, nfields, blks_info, &
                                        blks_sizes, blk_sizes, inds) result(offset)

    implicit none

    offset = 0
    k = 1

    do i = 1, ntuple

          offset = offset + get_triang_blks_offset(nblks_tuple(i, :), nfields(i), &
                            blks_info(i,:,:), blks_sizes(i, :),  &
                            inds(i, k:k + nfields(i) - 1)) * &
                            blks_sizes(i:ntuple)/blks_sizes(i)

          k = k + nfields(i)

    end do


  end function

  function get_triang_blks_offset(nblks, nfield, blk_info, blk_sizes, ind) &
           result(offset)

    implicit none

    offset = 0

    do i = 1, nblks

       offset = offset + get_triang_offset(blk_info(i,2), ind(blk_info(i,1): &
                            blk_info(i,1) + blk_info(i,2) - 1), blk_info(i,3)) * &
                            blk_sizes(i:nblks)/blk_sizes(i)

    end do


  end function


  function get_triang_offset(nfield, ind, pdim) result(offset)

    implicit none

    offset = 1

    do i = 1, nfield

       offset = offset + get_one_triangular_size(nfield - i, pdim - ind(i) + 1)

    end do


  end function

  subroutine sort_triangulated_indices(nfield, nblks, blk_info, indices)

    implicit none

    integer :: nfield, nblks, i, j, current_way
    integer, dimension(nblks, 3) :: blk_info
    integer, dimension(nfield) :: indices, sorted_indices

    current_blk_start = 1

    do i = 1, nblks

       do j = 1, blk_info(i, 2)

          current_minimum_index_position = current_blk_start
          current_minimum_index = indices(current_blk_start)

          do k = j + 1, blk_info(i, 2)
   
             if (indices(current_blk_start + k) < current_minimum_index) then

                current_minimum_index = indices(current_blk_start + k)
                current_minimum_index_position = current_blk_start + k

             end if

                ind_tmp = indices(current_blk_start + j - 1)
                indices(current_blk_start + j - 1) = current_minimum_index
                indices(current_blk_start + k - 1) = ind_tmp

          end do
      
       end do

    end do


  end function

subroutine make_triangulated_tuples_indices(ntuples, field_tuples, &
           triangulated_size, indices)

implicit none

allocate(individual_block_tuple_indices)

do i = 1, ntuples

call make_triangulated indices(..., individual_block_tuple_indices(i))

end do

call index_blks_direct_product(..., indices, ...)



deallocate(individual_block_tuple_indices)



end subroutine



  subroutine make_triangulated_indices(nblks, blk_info, triangulated_size, indices)

    implicit none

    integer :: nblks, i, triangulated_size, j
    integer, dimension(nblks) :: triang_sizes
    integer, dimension(nblks, 3) :: blk_info
    integer, dimension(triangulated_size, sum(blk_info(:,2))) :: indices
    type(triangulated_index_block), allocatable, dimension(:) :: blks

allocate(blks(nblks))

    do i = 1, nblks

       triang_sizes(i) = get_one_triangular_size(blk_info(i, 2), blk_info(i,3))

       allocate(blks(i)%t_ind(triang_sizes(i), blk_info(i, 2)))

       call make_one_triang_index_blk(blk_info(i, 2), blk_info(i, 3), 1, 1, 1, &
                                      triang_sizes(i), blks(i)%t_ind)

    end do

    call index_blks_direct_product(nblks, triang_sizes, blks, indices, &
                                   sum(blk_info(:,2)), 1, 1, 1)

    do i = 1, nblks

       deallocate(blks(i)%t_ind)

    end do

deallocate(blks)

  end subroutine

  recursive subroutine index_blks_direct_product(nblks, blk_sizes, blks, indices, &
                       nways, current_way, lvl, offset)

    implicit none

    integer :: nblks, current_way, lvl, i, j, offset, increment, nways, new_offset
    integer, dimension(nblks) :: blk_sizes
    type(triangulated_index_block), dimension(nblks) :: blks
    integer, dimension(product(blk_sizes), nways) :: indices

    if (lvl < nblks) then

       increment = product(blk_sizes(lvl:nblks))/blk_sizes(lvl)

       do i = 0, blk_sizes(lvl) - 1

          do j = 0, increment - 1

             indices(offset + i * increment + j, &
                     current_way:current_way + size(blks(lvl)%t_ind, 2) - 1) = &
                     blks(lvl)%t_ind(i + 1, :)

          end do

          new_offset = offset + i * increment
          call index_blks_direct_product(nblks, blk_sizes, blks, indices, &
               nways, current_way + size(blks(lvl)%t_ind, 2), lvl + 1, new_offset)

       end do


    elseif (lvl == nblks) then

       do i = 0, blk_sizes(lvl) - 1

          indices(offset + i, current_way:current_way + size(blks(lvl)%t_ind, 2) - 1) = &
          blks(lvl)%t_ind(i + 1, :)

       end do

    end if


  end subroutine


  recursive subroutine make_one_triang_index_blk(blk_size, pdim, st_ind, lvl, offset, &
                                                   triang_size, index_blk)

    implicit none

    integer :: blk_size, pdim, st_ind, lvl, i, j, offset, &
               triang_size, new_offset, increment
    integer, dimension(triang_size, blk_size) :: index_blk

    if (lvl < blk_size) then

       new_offset = offset

       do i = 0, pdim - st_ind

          increment = get_one_triangular_size(blk_size - lvl, pdim - st_ind - i + 1)
          index_blk(new_offset:new_offset + increment, lvl) = & 
          (i + st_ind) * (/ (j/j, j = 1, increment)/)

          call make_one_triang_index_blk(blk_size, pdim, st_ind + i, lvl + 1, &
                                         new_offset, triang_size, index_blk)

          new_offset = new_offset + increment

       end do

    elseif (lvl == blk_size) then

       do i = 0, pdim - st_ind

          index_blk(offset + i, lvl) = i + st_ind

       end do


    end if


  end subroutine

  function get_triangular_sizes(nblks, blk_nfield, pdims) result(blk_sizes)

    implicit none

    do i = 1, nblks

       blk_sizes(i) = get_one_triangular_size(blk_nfield(i), pdims(i))

    end do 

  end function


  function get_triangulated_size(nblks, blk_info)

    implicit none

    integer :: nblks, i, get_triangulated_size
    integer, dimension(nblks, 3) :: blk_info

    get_triangulated_size = 1

    do i = 1, nblks

       get_triangulated_size = get_triangulated_size * & 
                               get_one_triangular_size(blk_info(i, 2), blk_info(i, 3))

    end do

  end function


  function get_one_triangular_size(blk_size, pdim)

    implicit none

    integer :: get_one_triangular_size, blk_size, pdim

    get_one_triangular_size = fact_terminate_lower(pdim + blk_size - 1, pdim) / &
                              fact_terminate_lower(blk_size,  1)

  end function


  recursive function fact_terminate_lower(highest, lowest) result(ftl)

    implicit none

    integer :: ftl, highest, lowest

    if (highest == lowest) then

       ftl = highest   

    else

       ftl = highest * fact_terminate_lower(highest - 1, lowest)

    end if

  end function



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

end module