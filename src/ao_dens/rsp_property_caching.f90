! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and functions related to caching of contributions to the
! response property tensor.

module rsp_property_caching

  use rsp_field_tuple
  use rsp_indices_and_addressing
  use qcmatrix_f
!  use matrix_backend

  implicit none

  public property_cache_initialize
  public property_cache_next_element
  public property_cache_add_element
  public property_cache_already
  public property_cache_getdata
  public property_cache_allocate

  ! NEW 2014
  
 public contrib_cache_initialize
 public contrib_cache_next_element
 public contrib_cache_outer_next_element
 public contrib_cache_outer_cycle_first
 public contrib_cache_cycle_outer
 public contrib_cache_outer_add_element
 public contrib_cache_add_element
 public contrib_cache_already
 public contrib_cache_getdata
 public contrib_cache_allocate
 public contrib_cache_outer_allocate
 public contrib_cache_store
 public contrib_cache_retrieve
 public contrib_cache_outer_store
 public contrib_cache_outer_retrieve
 public mat_scal_store
 public mat_scal_retrieve
 public rs_check
 public prog_incr
 public prog_init
  
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

  ! NOTE FEB 2016: SET DEFAULT VALUES FOR ALL NON-ALLOCATE ATTRIBUTES OF CACHE TYPES

  
 type contrib_cache_outer

   type(contrib_cache_outer), pointer :: next
   logical :: last
   logical :: dummy_entry
   integer :: num_dmat
   type(p_tuple), allocatable, dimension(:) :: p_tuples

   ! Contribution type: 1: Only Pulay n, 3: Only Lagrange, 4: Both Pulay and Lagrange
   integer :: contrib_type = 0
   
   integer :: n_rule = 0
   integer :: contrib_size = 0
   integer, allocatable, dimension(:) :: nblks_tuple
   integer, allocatable, dimension(:,:) :: blk_sizes
   integer, allocatable, dimension(:,:) :: indices
   integer, allocatable, dimension(:,:,:) :: blks_tuple_info
   integer, allocatable, dimension(:) :: blks_tuple_triang_size
   type(QcMat), allocatable, dimension(:) :: data_mat ! Fock matrix contribution data
   complex(8), allocatable, dimension(:) :: data_scal ! Property data    

 end type 

 type contrib_cache

   type(contrib_cache), pointer :: next
   logical :: last
   type(p_tuple) :: p_inner

   integer :: num_outer
   integer :: nblks
   integer, allocatable, dimension(:) :: blk_sizes
   integer :: blks_triang_size
   integer, allocatable, dimension(:,:) :: blk_info
   integer, allocatable, dimension(:,:) :: indices
         
   type(contrib_cache_outer), pointer :: contribs_outer

 end type 

  
  
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
       call p1_cloneto_p2(p_tuples(i), new_element%p_tuples(i))
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
          total_num_perturbations = total_num_perturbations + p_tuples(i)%npert
       end do

       next_element%p_tuples = p_tuples_standardorder(next_element%num_p_tuples, &
                                                      next_element%p_tuples)


       if (p_tuples(1)%npert > 0) then

          call p1_cloneto_p2(p_tuples(1), merged_p_tuple)

          do i = 2, num_p_tuples

             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       else

          call p1_cloneto_p2(p_tuples(2), merged_p_tuple)

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
          do j = 1, p_tuples(i)%npert
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

          nfields(i) = p_tuples_st_order(i)%npert
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

          if (p_tuples(1)%npert > 0) then

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

    current_element%p_tuples(1)%npert = 0

    allocate(current_element%p_tuples(1)%pdim(1))
    allocate(current_element%p_tuples(1)%plab(1))
    allocate(current_element%p_tuples(1)%pid(1))
    allocate(current_element%p_tuples(1)%freq(1))
    allocate(current_element%data(1))

  end subroutine

    ! NEW 2014
    
  ! Initialize progress/restarting framework
  subroutine prog_init(rs_info)
  
    implicit none
    
    integer, dimension(3) :: rs_info
    logical :: r_exist
    
    inquire(file='OPENRSP_RESTART', exist=r_exist)
    
    if (r_exist) then
    
       open(unit=260, file='OPENRSP_RESTART', action='read') 
       read(260,*) rs_info(1)
       read(260,*) rs_info(2)
       read(260,*) rs_info(3)
       close(260)
    
    else
    
       rs_info = (/0, 0, 0/)
    
    end if
    
  end subroutine
     
  ! Check progress counter against restart checkpoint
  ! Return true if checkpoint not passed, false if passed
  ! If optional argument 'lvl' specified, only check progress at that level
  function rs_check(prog_info, rs_info, lvl)
  
    implicit none
    
    logical :: rs_check
    integer, dimension(3) :: prog_info, rs_info
    integer :: i
    integer, optional :: lvl
    
    rs_check = .TRUE.
    
    if (present(lvl)) then
    
       if (prog_info(lvl) > rs_info(lvl)) then
       
          rs_check = .FALSE.
       
       end if    
    
    else
    
       if (.NOT.(prog_info(1) < rs_info(1))) then
       
          if (prog_info(1) == rs_info(1)) then
       
             if (.NOT.(prog_info(2) < rs_info(2))) then
             
                if (prog_info(2) == rs_info(2)) then
          
                   if (.NOT.(prog_info(3) <= rs_info(3))) then
                   
                      rs_check = .FALSE.
                                
                   end if
                   
                else
                
                   rs_check = .FALSE.
                 
                end if
             
             end if
          
          else 
          
             rs_check = .FALSE.
          
          end if
                    
       end if
       
       
    
    end if
    
    
  end function
  
  
  ! Increment calculation progress counter and store in file
  subroutine prog_incr(prog_info, lvl)
    
    implicit none
    
    integer, dimension(3) :: prog_info
    integer :: lvl, i
        
    prog_info(lvl) = prog_info(lvl) + 1
    
    do i = lvl + 1, 3
    
       prog_info(i) = 0
    
    end do
    
    write(*,*) 'Increased progress at level', lvl
    write(*,*) 'Progress is now', prog_info
    
    
    open(unit=260, file='OPENRSP_RESTART', status='replace', action='write') 
       write(260,*) prog_info(1)
       write(260,*) prog_info(2)
       write(260,*) prog_info(3)
    close(260)
    
  end subroutine
 
 
 subroutine mat_scal_store(array_size, fname, mat, scal, start_pos)
 
   implicit none
   
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
   integer :: array_size, funit, i, k, first_i
   integer, optional :: start_pos
   complex(8), dimension(array_size), optional :: scal
   type(QcMat), dimension(array_size), optional :: mat
 
   funit = 260
 
   if (present(scal)) then
 
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='replace', action='write')
   
      write(funit) array_size
      write(funit) scal(1:array_size)
        
      close(funit)
      
    elseif (present(mat)) then
    
       fid_fmt = '(I10.10)'
          
       
       if (present(start_pos)) then
       
          first_i = start_pos
          
       else 
       
          first_i = 1
          
       end if
       
    
       do i = first_i, array_size
       
          write(str_fid, fid_fmt) i
       
          k = QcMatWrite_f(mat(i), trim(adjustl(fname)) // '_MAT_' // &
                           trim(str_fid) // '.DAT', BINARY_VIEW)
       
       end do
    
    end if
 
 end subroutine
 
  
 subroutine mat_scal_retrieve(array_size, fname, mat, scal)
 
   implicit none
   
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
   integer :: array_size, funit, i, k
   complex(8), dimension(array_size), optional :: scal
   type(QcMat), dimension(array_size), optional :: mat
 
   funit = 260
 
   if (present(scal)) then
 
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='old', action='read')
   
      read(funit) array_size
      read(funit) scal(1:array_size)
        
      close(funit)
      
    elseif (present(mat)) then
    
       fid_fmt = '(I10.10)'
          
       do i = 1, array_size
       
          write(str_fid, fid_fmt) i
       
          k = QcMatRead_f(mat(i), trim(adjustl(fname)) // '_MAT_' // &
                           trim(str_fid) // '.DAT', BINARY_VIEW)
       
       end do
    
    end if
 
 end subroutine
 
 
 
 
 

 
    
 subroutine contrib_cache_store(cache, fname)
 
   implicit none
   
   logical :: termination
   integer :: num_entries, i, j, funit, mat_acc
   character(*) :: fname
  
   type(contrib_cache), target :: cache
   type(contrib_cache), pointer :: cache_next
   type(contrib_cache_outer), pointer :: outer_next
 
 
   mat_acc = 0
 
   funit = 260
 
   open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='replace', action='write')
   
   ! Cycle to start
   cache_next => cache
   do while(.NOT.(cache_next%last))
      cache_next => cache_next%next
   end do
   cache_next => cache_next%next
   cache_next => cache_next%next
   
   num_entries = 0
   
   ! Traverse to get number of entries
   termination = .FALSE.
   do while(.NOT.(termination))
   
      num_entries = num_entries + 1
   
      termination = cache_next%last
      cache_next => cache_next%next
   
   end do
   
   write(funit) num_entries
   
   ! Cycle to start
   do while(.NOT.(cache_next%last))
      cache_next => cache_next%next
   end do
   cache_next => cache_next%next
   cache_next => cache_next%next
   
      
   ! Traverse and store
   termination = .FALSE.
   do i = 1, num_entries
   
      write(funit) cache_next%last
      write(funit) cache_next%p_inner%npert
      do j = 1, cache_next%p_inner%npert
      
         write(funit) cache_next%p_inner%pdim(j)
         write(funit) cache_next%p_inner%plab(j)
         write(funit) cache_next%p_inner%pid(j)
         write(funit) cache_next%p_inner%freq(j)
      
      end do
      
      write(funit) cache_next%num_outer
      write(funit) cache_next%nblks
      
      write(funit) size(cache_next%blk_sizes)
      write(funit) cache_next%blk_sizes
      
      write(funit) cache_next%blks_triang_size

      write(funit) size(cache_next%blk_info, 1)
      write(funit) size(cache_next%blk_info, 2)
      write(funit) cache_next%blk_info
     
      write(funit) size(cache_next%indices, 1)
      write(funit) size(cache_next%indices, 2)
      write(funit) cache_next%indices
      
         
      
      outer_next => cache_next%contribs_outer
      
      ! num_outer may not be well-defined; consider 
      ! switching to traversal with termination instead
         
      call contrib_cache_outer_put(outer_next, fname, funit, mat_acc_in=mat_acc)
   
      cache_next => cache_next%next
   
   end do
   
   
   close(funit)
 
 end subroutine
 
 
 subroutine contrib_cache_retrieve(cache, fname)
 
   implicit none
   
   logical :: termination, cache_ext
   integer :: num_entries, i, j, funit
   character(*) :: fname
  
   type(contrib_cache), target :: cache
   type(contrib_cache), pointer :: cache_next, cache_new, cache_orig
   type(contrib_cache_outer), pointer :: outer_next
 
   inquire(file=fname // '.DAT', exist=cache_ext)
   
   if (.NOT.(cache_ext)) then
   
      write(*,*) 'ERROR: The expected cache file ', trim(adjustl(fname)), ' does not exist'
      stop
   
   end if
 
   funit = 260
 
   open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='old', action='read')
   
   read(funit) num_entries

   ! Cache allocation may need own subroutine because of target/pointer issues
   ! Set up first element of cache
   
   ! NOTE: DEALLOCATION/ALLOCATION MAY BE PROBLEMATIC   
    
!    cache_next => cache
   
!    if (allocated(cache)) then
!    
!       nullify(cache)
!       
!      
!    end if
!    
!    allocate(cache)
   cache_next => cache
   
   call contrib_cache_read_one_inner(cache_next, fname, funit)
   
   allocate(cache_next%contribs_outer)
   outer_next => cache_next%contribs_outer
   call contrib_cache_outer_retrieve(outer_next, fname, .TRUE., funit)
   
   cache_orig => cache_next
   
   write(*,*) 'outer read 1'
   
   ! Get all remaining elements
   do i = 1, num_entries - 1
   
      call contrib_cache_retrieve_one_inner(cache_next%next, fname, funit)
      write(*,*) 'inner read', i + 1
      
      cache_next => cache_next%next
      allocate(cache_next%contribs_outer)
      outer_next => cache_next%contribs_outer
         
      call contrib_cache_outer_retrieve(outer_next, fname, .TRUE., funit)
      
      write(*,*) 'outer read', i + 1
  
   end do
   
   ! Bite the tail
   cache_next%next => cache_orig
   
   
   close(funit)
 
 end subroutine
 
 subroutine contrib_cache_retrieve_one_inner(cache_new, fname, funit)
 
    implicit none
 
    integer :: funit 
    character(*) :: fname
!     type(contrib_cache), target :: cache
    type(contrib_cache), pointer :: cache_new
 
    allocate(cache_new)
 
    call contrib_cache_read_one_inner(cache_new, fname, funit)
    
!     cache%next => cache_new
 
 end subroutine
 
 
 
 
 subroutine contrib_cache_read_one_inner(cache, fname, funit)
 
   implicit none
   
   integer :: num_entries, i, j, funit
   integer :: size_i, size_j
   character(*) :: fname
  
   type(contrib_cache), pointer :: cache
   
   read(funit) cache%last
   read(funit) cache%p_inner%npert
   
   write(*,*) 'cache last', cache%last
   write(*,*) 'cache%p_inner%npert', cache%p_inner%npert
   
   if (cache%p_inner%npert == 0) then

      allocate(cache%p_inner%pdim(1))
      allocate(cache%p_inner%plab(1))
      allocate(cache%p_inner%pid(1))
      allocate(cache%p_inner%freq(1))
      
   else
    
      allocate(cache%p_inner%pdim(cache%p_inner%npert))
      allocate(cache%p_inner%plab(cache%p_inner%npert))
      allocate(cache%p_inner%pid(cache%p_inner%npert))
      allocate(cache%p_inner%freq(cache%p_inner%npert))
   
   end if
   
   do j = 1, cache%p_inner%npert
      
         read(funit) cache%p_inner%pdim(j)
         read(funit) cache%p_inner%plab(j)
         read(funit) cache%p_inner%pid(j)
         read(funit) cache%p_inner%freq(j)
      
   end do
      
      read(funit) cache%num_outer
      read(funit) cache%nblks
      
      read(funit) size_i
      allocate(cache%blk_sizes(size_i))
      read(funit) cache%blk_sizes
      
      read(funit) cache%blks_triang_size

      read(funit) size_i
      read(funit) size_j
      allocate(cache%blk_info(size_i, size_j))
      read(funit) cache%blk_info
     
      read(funit) size_i
      read(funit) size_j
      allocate(cache%indices(size_i, size_j))
      read(funit) cache%indices
      
           
         
      
   
 end subroutine
 
 

 
 
 
 subroutine contrib_cache_outer_store(cache, fname, mat_acc_in)
 
   implicit none
   
   type(contrib_cache_outer) :: cache
   character(*) :: fname
   integer, optional :: mat_acc_in
   integer :: funit, mat_acc
   
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if   
   
   funit = 260
 
   open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='replace', action='write')
        
   call contrib_cache_outer_put(cache, fname, funit, mat_acc_in=mat_acc)
        
        
   close(funit)
 
 
 end subroutine
 
 
 ! Assumes file 'fname' is already opened with I/O unit 'funit'
 subroutine contrib_cache_outer_put(cache, fname, funit, mat_acc_in)
 
   implicit none
   
   logical :: termination
   integer :: num_entries, i, j, k, mat_scal_none, funit, mat_acc
   integer, optional :: mat_acc_in
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
  
   type(contrib_cache_outer), target :: cache
   type(contrib_cache_outer), pointer :: cache_next
 
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
   
 
   ! Cycle to start
   cache_next => cache
   do while(.NOT.(cache_next%last))
      cache_next => cache_next%next
   end do
   cache_next => cache_next%next
   cache_next => cache_next%next
   
   num_entries = 0
   
   ! Traverse to get number of entries
   termination = .FALSE.
   do while(.NOT.(termination))
   
      num_entries = num_entries + 1
   
      termination = cache_next%last
      cache_next => cache_next%next
   
   end do
   
   write(funit) num_entries
   
   ! Cycle to start
   do while(.NOT.(cache_next%last))
      cache_next => cache_next%next
   end do
   cache_next => cache_next%next
   cache_next => cache_next%next
  
   ! Traverse and store
   do i = 1, num_entries
   
      write(funit) size(cache_next%p_tuples)
      do j = 1, size(cache_next%p_tuples)
      
         write(funit) cache_next%p_tuples(j)%npert
         write(*,*) 'wrote npert', cache_next%p_tuples(j)%npert
         
         do k = 1, cache_next%p_tuples(j)%npert
         
            write(funit) cache_next%p_tuples(j)%pdim(k)
            write(funit) cache_next%p_tuples(j)%plab(k)
            write(funit) cache_next%p_tuples(j)%pid(k)
            write(funit) cache_next%p_tuples(j)%freq(k)
         
         end do
         
         write(*,*) 'stored plab', cache_next%p_tuples(j)%plab
      
      end do
      
   
      write(funit) cache_next%last
      write(funit) cache_next%dummy_entry
      write(funit) cache_next%num_dmat
      write(funit) cache_next%contrib_type
      write(funit) cache_next%n_rule
      write(funit) cache_next%contrib_size
      
      write(*,*) 'Storage: Dummy entry?', cache_next%dummy_entry
      
      if (cache_next%p_tuples(1)%npert > 0) then
      
         write(funit) size(cache_next%nblks_tuple)
         write(funit) cache_next%nblks_tuple
      
         write(funit) size(cache_next%blk_sizes, 1)
         write(funit) size(cache_next%blk_sizes, 2)
         do j = 1, size(cache_next%blk_sizes, 1)
            write(funit) cache_next%blk_sizes(j, :)
         end do
      
         write(funit) size(cache_next%indices, 1)
         write(funit) size(cache_next%indices, 2)
         
         do j = 1, size(cache_next%indices, 1)
            write(funit) cache_next%indices(j,:)
         end do
      
         write(funit) size(cache_next%blks_tuple_info, 1)
         write(funit) size(cache_next%blks_tuple_info, 2)
         write(funit) size(cache_next%blks_tuple_info, 3)
      
         do j = 1, size(cache_next%blks_tuple_info, 1)
            do k = 1, size(cache_next%blks_tuple_info, 2)
               write(funit) cache_next%blks_tuple_info(j, k, :)
            end do
         end do
      
         write(funit) size(cache_next%blks_tuple_triang_size)
         write(funit) cache_next%blks_tuple_triang_size
      
      end if
      
      if (allocated(cache_next%data_mat)) then
      
         mat_scal_none = 0
         write(funit) mat_scal_none
         write(funit) size(cache_next%data_mat)
         
         fid_fmt = '(I10.10)'
         
         do j = 1, size(cache_next%data_mat)
            
            write(str_fid, fid_fmt) j + mat_acc
         
            k = QcMatWrite_f(cache_next%data_mat(j), trim(adjustl(fname)) // '_MAT_' // &
                 trim(str_fid) // '.DAT', BINARY_VIEW)
         
         end do
         
         mat_acc = mat_acc + size(cache_next%data_mat)
         
      
      elseif (allocated(cache_next%data_scal)) then
      
         mat_scal_none = 1
         write(funit) mat_scal_none
         write(funit) size(cache_next%data_scal)
         write(funit) cache_next%data_scal
         
      
      else
      
         mat_scal_none = 2 
         write(funit) mat_scal_none
      
      end if
            
      cache_next => cache_next%next
   
   
   end do
 
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
 
 
 end subroutine
 
 
 ! Add condition to handle "not from inner"
 subroutine contrib_cache_outer_retrieve(cache, fname, from_inner, funit_in, mat_acc_in)
 
   implicit none
   
   logical :: from_inner, cache_ext
   integer :: funit, i, j, k, mat_acc
   integer, optional :: funit_in, mat_acc_in
   integer :: num_entries
   character(*) :: fname
   type(contrib_cache_outer), target :: cache
   type(contrib_cache_outer), pointer :: cache_next, cache_orig
   
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
   

   funit = 260
   if (present(funit_in)) then
      funit = funit_in
   end if
   
   
   ! If not as part of inner cache, then open file
   if (.NOT.(from_inner)) then
   
      write(*,*) 'Opening ', trim(adjustl(fname))
   
      inquire(file=fname // '.DAT', exist=cache_ext)
      
      if (.NOT.(cache_ext)) then
   
         write(*,*) 'ERROR: The expected cache file ', trim(adjustl(fname)), ' does not exist'
         stop
   
      end if
      
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='old', action='read')
           
   end if
   
   
   
   read(funit) num_entries
   
   write(*,*) 'num entries', num_entries
 
   call contrib_cache_read_one_outer(cache, fname, funit, mat_acc_in=mat_acc)

   write(*,*) 'read one outer'
   
   cache_next => cache
   cache_orig => cache
   
   do i = 1, num_entries - 1
   
      write(*,*) 'reading another outer'
   
      call contrib_cache_retrieve_one_outer(cache_next, fname, funit, mat_acc_in=mat_acc)
      cache_next => cache_next%next
    
   end do

   cache_next%next => cache_orig
   
   
   if (.NOT.(from_inner)) then
      
      close(260)
      
   end if
   
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
 
 end subroutine
 
 
 
 subroutine contrib_cache_retrieve_one_outer(cache, fname, funit, mat_acc_in)
 
   implicit none
   
   integer :: funit, mat_acc
   integer, optional :: mat_acc_in
   character(*) :: fname
   type(contrib_cache_outer), target :: cache
   type(contrib_cache_outer), pointer :: cache_new
   
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
   
   allocate(cache_new)
   
   call contrib_cache_read_one_outer(cache_new, fname, funit, mat_acc_in=mat_acc)
   
   cache%next => cache_new
  
 
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
 
 end subroutine
 
 
 subroutine contrib_cache_read_one_outer(cache, fname, funit, mat_acc_in)
 
   implicit none
   
   logical :: from_inner
   integer :: funit, i, j, k
   integer :: size_i, size_j, size_k
   integer :: mat_scal_none
   integer, optional :: mat_acc_in
   integer :: mat_acc
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
   type(contrib_cache_outer) :: cache
   
   
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
  
   
   read(funit) size_i
   write(*,*) 'num p tuples', size_i
   allocate(cache%p_tuples(size_i))
   
      do j = 1, size(cache%p_tuples)
      
         read(funit) cache%p_tuples(j)%npert
         
         write(*,*) 'read npert', cache%p_tuples(j)%npert
         
         allocate(cache%p_tuples(j)%pdim(cache%p_tuples(j)%npert))
         allocate(cache%p_tuples(j)%plab(cache%p_tuples(j)%npert))
         allocate(cache%p_tuples(j)%pid(cache%p_tuples(j)%npert))
         allocate(cache%p_tuples(j)%freq(cache%p_tuples(j)%npert))
         
         do k = 1, cache%p_tuples(j)%npert
         
            read(funit) cache%p_tuples(j)%pdim(k)
            read(funit) cache%p_tuples(j)%plab(k)
            read(funit) cache%p_tuples(j)%pid(k)
            read(funit) cache%p_tuples(j)%freq(k)
         
         end do
         
         write(*,*) 'cache plab at j = ', j, ' is', cache%p_tuples(j)%plab
      
      end do
   
  
   read(funit) cache%last
   read(funit) cache%dummy_entry
   read(funit) cache%num_dmat
   read(funit) cache%contrib_type
   read(funit) cache%n_rule
   read(funit) cache%contrib_size
   
   write(*,*) 'Dummy entry?', cache%dummy_entry

   if (cache%p_tuples(1)%npert > 0) then
   
      read(funit) size_i
      allocate(cache%nblks_tuple(size_i))
      read(funit) cache%nblks_tuple
      
      read(funit) size_i
      read(funit) size_j
      allocate(cache%blk_sizes(size_i, size_j))
      do j = 1, size_i
         read(funit) cache%blk_sizes(j, :)
      end do   
     
      read(funit) size_i
      read(funit) size_j
      allocate(cache%indices(size_i, size_j))
      do j = 1, size_i
         read(funit) cache%indices(j,:)
      end do
      
      read(funit) size_i
      read(funit) size_j
      read(funit) size_k
      allocate(cache%blks_tuple_info(size_i, size_j, size_k))
      do j = 1, size_i
         do k = 1, size_j
            read(funit) cache%blks_tuple_info(j, k, :)
         end do
      end do
      
      read(funit) size_i
      allocate(cache%blks_tuple_triang_size(size_i))
      read(funit) cache%blks_tuple_triang_size
   
   end if
   
   read(funit) mat_scal_none
      
      
   if (mat_scal_none == 0) then
   
      fid_fmt = '(I10.10)'
   
      read(funit) size_i
      
      allocate(cache%data_mat(size_i))
      
      do i = 1, size_i
      
         write(str_fid, fid_fmt) i + mat_acc
         
         write(*,*) 'Reading matrix data from', trim(adjustl(fname)) // '_MAT_' // &
                          trim(str_fid) // '.DAT'
      
         call QcMatInit(cache%data_mat(i))
         k = QcMatRead_f(cache%data_mat(i), trim(adjustl(fname)) // '_MAT_' // &
                          trim(str_fid) // '.DAT', BINARY_VIEW)
      
      end do
      
      mat_acc = mat_acc + size_i
   
   
   elseif (mat_scal_none == 1) then
   
      read(funit) size_i
      allocate(cache%data_scal(size_i))
      read(funit) cache%data_scal
   
   end if
      
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
   
 end subroutine
 
 
    
  
 subroutine contrib_cache_cycle_outer(current_element, num_p_tuples, p_tuples, &
            next_outer, n_rule)

   implicit none

   integer :: num_p_tuples, i, passedlast
   logical :: found
   integer, optional :: n_rule
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: new_element
   type(contrib_cache), pointer :: next_element
   type(contrib_cache_outer), pointer :: next_outer
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert

!    write(*,*) 'plab 1', p_tuples(1)%plab
!       write(*,*) 'curr elem plab', current_element%p_inner%plab
!       write(*,*) 'next', current_element%next%p_inner%plab
!       write(*,*) 'next', current_element%next%next%p_inner%plab
!       write(*,*) 'next', current_element%next%next%next%p_inner%plab
!       write(*,*) 'next', current_element%next%next%next%next%p_inner%plab
!       write(*,*) 'next', current_element%next%next%next%next%next%p_inner%plab
      
   
   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(current_element, p_tuples(1))) then
   
      next_element => current_element
      
!       write(*,*) 'inner', next_element%p_inner%plab

   
      ! Skip to cache element for this inner
      do while (p_tuple_compare(next_element%p_inner, p_tuples(1)) .EQV. .FALSE.)

!        write(*,*) 'inner compare', next_element%p_inner%plab
      
        next_element => next_element%next
         
      end do
      
!       write(*,*) 'inner skip OK'
      next_outer => next_element%contribs_outer
      
      passedlast = 0
      found = .FALSE.

      do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

         next_outer => contrib_cache_outer_next_element(next_outer)

!          write(*,*) 'next outer', next_outer%p_tuples(1)%plab
         
         if (num_p_tuples > 1) then
         
            found = p_tuples_compare(num_p_tuples - 1, next_outer%p_tuples, &
                                     p_tuples(2:))
                                  
         else
         
            found = p_tuples_compare(num_p_tuples - 1, next_outer%p_tuples, &
                                     (/get_emptypert()/))
         
         end if
                                  
         if (present(n_rule)) then
      
            found = found .AND. (n_rule == next_outer%n_rule)
      
         end if

         if (next_outer%last) then
            passedlast = passedlast + 1
         end if

      end do
      
      if (.NOT.(found)) then
      
         write(*,*) 'ERROR: Did not find expected outer cache element'
         
      end if
      
   else
   
      write(*,*) 'ERROR: Did not find expected cache element'
      
   end if
   

!    write(*,*) 'end of function', found
   
 end subroutine
 
 
 function contrib_cache_cycle_first(current_element) result(next_element)

   implicit none

   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element
   
   next_element => current_element
   
   do while (next_element%last .eqv. .FALSE.)
          next_element => next_element%next
   end do

   next_element => next_element%next   

 end function

 function contrib_cache_outer_cycle_first(current_element) result(next_element)

   implicit none

   type(contrib_cache_outer), target :: current_element
   type(contrib_cache_outer), pointer :: next_element
   
   next_element => current_element
!    write(*,*) 'out cycling', next_element%dummy_entry, next_element%last
   do while (next_element%last .eqv. .FALSE.)
!           write(*,*) 'cycling', next_element%dummy_entry, next_element%last
   
          next_element => next_element%next
   end do

   next_element => next_element%next   

 end function
 
 function contrib_cache_next_element(current_element) result(next_element)

   implicit none

   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element

   next_element => current_element%next

 end function

 function contrib_cache_outer_next_element(current_element) result(next_element)

   implicit none

   type(contrib_cache_outer), target :: current_element
   type(contrib_cache_outer), pointer :: next_element

   next_element => current_element%next

 end function
 
 
  subroutine contrib_cache_allocate(current_element)

   implicit none

   type(contrib_cache), pointer :: current_element
   !type(contrib_cache), pointer :: current_element
   allocate(current_element)

   current_element%next => current_element
   current_element%last = .TRUE.

   current_element%p_inner%npert = 0
   
   allocate(current_element%p_inner%pdim(1))
   allocate(current_element%p_inner%plab(1))
   allocate(current_element%p_inner%pid(1))
   allocate(current_element%p_inner%freq(1))
   
   current_element%p_inner%pdim = (/0/)
   current_element%p_inner%plab = (/'NUTN'/)
   current_element%p_inner%pid = (/0/)
   current_element%p_inner%freq = (/0.0/)
   
   call contrib_cache_outer_allocate(current_element%contribs_outer)
   
 end subroutine

 subroutine contrib_cache_outer_allocate(current_element)

    implicit none

    type(contrib_cache_outer), pointer :: current_element

    allocate(current_element)

    current_element%next => current_element
    current_element%num_dmat = 0
    current_element%last = .TRUE.
    current_element%dummy_entry = .TRUE.

    allocate(current_element%p_tuples(1))

    current_element%p_tuples(1)%npert = 0

    allocate(current_element%p_tuples(1)%pdim(1))
    allocate(current_element%p_tuples(1)%plab(1))
    allocate(current_element%p_tuples(1)%pid(1))
    allocate(current_element%p_tuples(1)%freq(1))
        
   current_element%p_tuples(1)%pdim = (/0/)
   current_element%p_tuples(1)%plab = (/'NUTN'/)
   current_element%p_tuples(1)%pid = (/0/)
   current_element%p_tuples(1)%freq = (/0.0/)

  end subroutine
 
 subroutine contrib_cache_initialize(new_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   integer :: num_p_tuples
   integer, optional :: n_rule
   type(contrib_cache) :: new_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples

!    write(*,*) 'new contrib cache inner: ', p_tuples(1)%plab
   
   new_element%last = .TRUE.
   call p1_cloneto_p2(p_tuples(1), new_element%p_inner)
   
   new_element%nblks = get_num_blks(new_element%p_inner)
  
   allocate(new_element%blk_sizes(new_element%nblks))
   allocate(new_element%blk_info(new_element%nblks, 3))
   new_element%blk_info = get_blk_info(new_element%nblks, new_element%p_inner)
   new_element%blk_sizes = get_triangular_sizes(new_element%nblks, &
                           new_element%blk_info(:,2),new_element%blk_info(:,3))
                   
   new_element%blks_triang_size = product(new_element%blk_sizes)
  
   allocate(new_element%indices(new_element%blks_triang_size, new_element%p_inner%npert))
   
   call make_triangulated_indices(new_element%nblks, new_element%blk_info, &
        new_element%blks_triang_size, new_element%indices)
  
!    write(*,*) 'nblks', new_element%nblks
!    write(*,*) 'blk sizes', new_element%blk_sizes
!    write(*,*) 'blks triang size', new_element%blks_triang_size
  
  
   if (num_p_tuples > 1) then
   
   

   
   
!    write(*,*) 'outer case a'
      call contrib_cache_outer_allocate(new_element%contribs_outer)

      if (present(n_rule)) then
      
         call contrib_cache_outer_add_element(new_element%contribs_outer, &
              .NOT.( p_tuples(2)%npert > 0), num_p_tuples - 1, &
               p_tuples(2:num_p_tuples), n_rule=n_rule)   
         
      else
         
         call contrib_cache_outer_add_element(new_element%contribs_outer, .FALSE., num_p_tuples - 1, &
                                              p_tuples(2:num_p_tuples))   
         
      end if
      
   
   
   else
   
!    write(*,*) 'outer case b'
      call contrib_cache_outer_allocate(new_element%contribs_outer)
      
      if (present(n_rule)) then
         
         call contrib_cache_outer_add_element(new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                              p_tuples(2:num_p_tuples), n_rule=n_rule)
         
      else
      
         call contrib_cache_outer_add_element(new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                              p_tuples(2:num_p_tuples))
         
      end if
      

      
   end if
      
 end subroutine
 
 
 subroutine contrib_cache_outer_initialize(new_element, unperturbed, num_dmat, outer_p_tuples)

   implicit none

   logical :: unperturbed
   integer :: num_dmat, i, total_npert
!    integer, allocatable, dimension(:) :: nblks_tuple, blks_tuple_triang_size
!    integer, allocatable, dimension(:,:) :: blk_sizes
!    integer, allocatable, dimension(:,:,:) :: blks_tuple_info
   type(contrib_cache_outer) :: new_element
   type(p_tuple), dimension(num_dmat) :: outer_p_tuples

   new_element%num_dmat = num_dmat
   new_element%last = .TRUE.
   new_element%dummy_entry = .FALSE.
   
! 
!       write(*,*) 'New outer cache'
!       write(*,*) 'num dmats', num_dmat

   if (.NOT.(unperturbed)) then
!       do i = 1, num_dmat
!         write(*,*) 'i', i
!         write(*,*) outer_p_tuples(i)%npert
!         write(*,*) outer_p_tuples(i)%pdim
!         write(*,*) outer_p_tuples(i)%plab
!         write(*,*) 'D', outer_p_tuples(i)%pid
!       end do
   
   
     allocate(new_element%p_tuples(num_dmat))
   
     do i = 1, num_dmat
   
        call p1_cloneto_p2(outer_p_tuples(i), new_element%p_tuples(i))
     end do

     total_npert = sum((/(outer_p_tuples(i)%npert, i = 1, num_dmat)/))
     
!      write(*,*) 'total_npert', total_npert
     
     if (total_npert > 0) then
     
        allocate(new_element%nblks_tuple(num_dmat))
        allocate(new_element%blk_sizes(num_dmat, total_npert))
        allocate(new_element%blks_tuple_info(num_dmat, total_npert, 3))
        allocate(new_element%blks_tuple_triang_size(num_dmat))
     
     
     new_element%nblks_tuple = (/(get_num_blks(outer_p_tuples(i)), i = 1, num_dmat)/)
   
!      write(*,*) 'new_element%nblks_tuple', new_element%nblks_tuple   
   
     do i = 1, num_dmat
     
        new_element%blks_tuple_info(i, :, :) = get_blk_info(new_element%nblks_tuple(i), outer_p_tuples(i))
        new_element%blk_sizes(i, 1:new_element%nblks_tuple(i)) = &
        get_triangular_sizes(new_element%nblks_tuple(i), &
        new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),2), &
        new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),3))
        new_element%blks_tuple_triang_size(i) = get_triangulated_size(new_element%nblks_tuple(i), &
                                   new_element%blks_tuple_info(i, 1:new_element%nblks_tuple(i), :))
                                   
     end do
     
!      write(*,*) 'new_element%blks_tuple_triang_size', new_element%blks_tuple_triang_size
     allocate(new_element%indices(product(new_element%blks_tuple_triang_size), total_npert))
     call make_triangulated_tuples_indices(num_dmat, total_npert, new_element%nblks_tuple, &
          new_element%blks_tuple_info, new_element%blks_tuple_triang_size, new_element%indices)
     
     end if
     
   
   
   else
   
     allocate(new_element%p_tuples(1)) 
   
      new_element%p_tuples%npert = 0
   
   allocate(new_element%p_tuples(1)%pdim(1))
   allocate(new_element%p_tuples(1)%plab(1))
   allocate(new_element%p_tuples(1)%pid(1))
   allocate(new_element%p_tuples(1)%freq(1))
   
   new_element%p_tuples(1)%pdim = (/0/)
   new_element%p_tuples(1)%plab = (/'NUTN'/)
   new_element%p_tuples(1)%pid = (/0/)
   new_element%p_tuples(1)%freq = (/0.0/)
   
!    write(*,*) 'initialized empty outer'
   end if
   
      
 end subroutine
   

   

 subroutine contrib_cache_outer_add_element(curr_element, unperturbed, num_dmat, &
            outer_p_tuples, data_size, data_mat, data_scal, n_rule)

   implicit none

   logical :: unperturbed, found_element, already
   integer :: num_dmat, i, j, passedlast
   integer, optional :: data_size, n_rule
   type(contrib_cache_outer), target :: curr_element
   type(contrib_cache_outer), pointer :: new_element
   type(contrib_cache_outer), pointer :: new_element_ptr
   type(contrib_cache_outer), pointer :: next_element
   type(p_tuple), dimension(num_dmat) :: outer_p_tuples, p_tuples_st_order
   
   type(Qcmat), optional, dimension(*) :: data_mat
   complex(8), optional, dimension(*) :: data_scal

   if (present(n_rule)) then
   
      already = contrib_cache_already_outer(curr_element, num_dmat, outer_p_tuples, n_rule=n_rule)
   
   else
   
      already = contrib_cache_already_outer(curr_element, num_dmat, outer_p_tuples)
   
   end if
   
!    write(*,*) 'already?', already
   
   if (already) then
   
      next_element => curr_element
      passedlast = 0
      p_tuples_st_order = p_tuples_standardorder(num_dmat, outer_p_tuples)
      found_element = .FALSE.

      do while ((passedlast < 2) .AND. (found_element .eqv. .FALSE.))
      
         

         next_element => next_element%next
         
!          write(*,*) 'size data mat a', size(next_element%data_mat)
         
         if (next_element%num_dmat == num_dmat .AND. .NOT.(next_element%dummy_entry)) then
            found_element = p_tuples_compare(next_element%num_dmat, &
            p_tuples_standardorder(next_element%num_dmat, &
            next_element%p_tuples), p_tuples_st_order)
         end if
            
         if (next_element%last .eqv. .TRUE.) then
            passedlast = passedlast + 1
         end if

      end do
      
!       write(*,*) 'size data mat b', size(next_element%data_mat)
! 
!       write(*,*) 'Found?', found_element
!       write(*,*) 'lab 1', p_tuples_st_order(1)%plab
!       write(*,*) 'lab 2', next_element%p_tuples(1)%plab
!       
!       write(*,*) 'size data mat c', size(next_element%data_mat)
      
      if(present(data_mat)) then
      
!       write(*,*) 'lab, freq', outer_p_tuples(1)%plab, outer_p_tuples(1)%freq
!       write(*,*) 'size data mat d', size(next_element%data_mat)
!          write(*,*) 'data size', data_size
      
         do i = 1, data_size
            call QcMatAEqB(next_element%data_mat(i),  data_mat(i))
           
!            if (size(outer_p_tuples(1)%plab) > 0) then
!            if (outer_p_tuples(1)%plab(1) .eq. 'GEO ') then
!            
!            write(*,*) 'adding element'
!            
!            j = QcMatWrite_f(next_element%data_mat(i), 'data_mat_i_cache', ASCII_VIEW)
!            j = QcMatWrite_f(data_mat(i), 'data_mat_i', ASCII_VIEW)
!            
!            end if
!            end if
           
         end do
   
   
   
      end if
   
      if (present(data_scal)) then
   
      end if
      
   
   else
   
      next_element => curr_element
      allocate(new_element)
      call contrib_cache_outer_initialize(new_element, unperturbed, num_dmat, outer_p_tuples)
      new_element_ptr => new_element
   
      do while (next_element%last .eqv. .FALSE.)
         next_element => contrib_cache_outer_next_element(next_element)
      end do
   
      if(present(data_mat)) then
      
         if (.NOT.(allocated(new_element%data_mat))) then

            allocate(new_element%data_mat(data_size))
   
            do i = 1, data_size
               call QcMatInit(new_element%data_mat(i), data_mat(i))
            end do
      
         end if
   
         do i = 1, data_size
            call QcMatAEqB(new_element%data_mat(i),  data_mat(i))
!             j = QcMatWrite_f(data_mat(i), 'dmi', ASCII_VIEW)
!             j = QcMatWrite_f(new_element%data_mat(i), 'dmj', ASCII_VIEW)
         end do
   
      end if
   
      if (present(data_scal)) then
   
      end if
      
      if (present(n_rule)) then
      
!          write(*,*) 'assigned n rule', n_rule
         new_element%n_rule = n_rule
      
      end if

      next_element%last = .FALSE.
      new_element%next => next_element%next
      next_element%next => new_element
   
   end if
      
 end subroutine
 
 subroutine contrib_cache_add_element(current_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   integer :: num_p_tuples, i
   integer, optional :: n_rule
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: new_element
   type(contrib_cache), pointer :: new_element_ptr
   type(contrib_cache), pointer :: next_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert

   

   
!     write(*,*) 'adding cache element, p tuples = ', num_p_tuples
!     
!     do i = 1, num_p_tuples
!     
!     if (i == 1) then
!     
!     write(*,*) 'E', p_tuples(i)%pid
!     
!     else
!     
!     write(*,*) 'D',  p_tuples(i)%pid
!     
!     end if
!     
!     
!     end do
   
   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(current_element, p_tuples(1))) then
   
      next_element => current_element
   
      ! Skip to cache element for this inner
      do while (p_tuple_compare(next_element%p_inner, p_tuples(1)) .EQV. .FALSE.)

!       write(*,*) 'skiparoo'
      
         next_element => next_element%next

          
      end do
!          write(*,*) 'pdim of skip', next_element%p_inner%freq
      ! NOTE: AND CONDITION MAY ADD WRONG KIND OF CACHE ELEMENT, REVISIT IF ERROR
      if (num_p_tuples > 1 .AND. p_tuples(2)%npert > 0) then
     
!       write(*,*) 'add case a'
     
         next_element%num_outer = next_element%num_outer + 1
     
         if (present(n_rule)) then
         
            call contrib_cache_outer_add_element(next_element%contribs_outer, .FALSE., &
                 num_p_tuples - 1, p_tuples(2:num_p_tuples), n_rule)
         
         else
         
            call contrib_cache_outer_add_element(next_element%contribs_outer, .FALSE., &
                 num_p_tuples - 1, p_tuples(2:num_p_tuples))
         
         end if
     
         
            
      else
!       write(*,*) 'add case b'
     
         next_element%num_outer = next_element%num_outer + 1

          if (present(n_rule)) then
         
             call contrib_cache_outer_add_element(next_element%contribs_outer, .TRUE., 1, (/emptypert/), n_rule)
         
          else
             ! call empty_p_tuple(emptypert)
             call contrib_cache_outer_add_element(next_element%contribs_outer, .TRUE., 1, (/emptypert/))
         
          end if
         

       
      end if


   ! Otherwise, add both inner and outer    
   else
   
      next_element => current_element


      
      allocate(new_element)
      
!       write(*,*) 'adding inner and outer'
      
      
      if (present(n_rule)) then
         
         call contrib_cache_initialize(new_element, num_p_tuples, p_tuples, n_rule)   
         
      else
         
         call contrib_cache_initialize(new_element, num_p_tuples, p_tuples)   
         
      end if
         

      

!       write(*,*) 'pdim of non-skip', new_element%p_inner%freq
      
      new_element%num_outer = 1
      new_element_ptr => new_element

      do while (next_element%last .EQV. .FALSE.)
         next_element => contrib_cache_next_element(next_element)
      end do

      next_element%last = .FALSE.
      new_element%next => next_element%next
      next_element%next => new_element
      
   end if

 end subroutine
 
 
 ! MAYBE SOME WORK REMAINING ON THIS FUNCTION
 
 function contrib_cache_already(current_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   logical :: contrib_cache_already
   integer :: passedlast, passedlast_outer, num_p_tuples, i
   integer, optional :: n_rule
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert, p_tuple_ord, p_tmp_ord

!       write(*,*) 'n p tuples for already', num_p_tuples
   
   next_element => current_element
   passedlast = 0
   
   

   
!    do i = 1, num_p_tuples
!    
!    write(*,*) 'p tuple', i
!    
!    write(*,*) 'npert', p_tuples(i)%npert
!    write(*,*) 'pid', p_tuples(i)%pid
!    write(*,*) 'plab', p_tuples(i)%plab
!    write(*,*) 'pfreq', p_tuples(i)%freq
!    write(*,*) 'pdim', p_tuples(i)%pdim
!    
!    write(*,*) ' '
!    
!    end do

   p_tuple_ord = p_tuple_standardorder(p_tuples(1))
   contrib_cache_already = .FALSE.

   call empty_p_tuple(emptypert)
   
   ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
   ! COULD THIS BE DONE IN ANOTHER WAY?
   do while ((passedlast < 2) .AND. (contrib_cache_already .eqv. .FALSE.))

      next_element => contrib_cache_next_element(next_element)

      p_tmp_ord = p_tuple_standardorder(next_element%p_inner)
      
!       write(*,*) 'p tmp ord', p_tmp_ord%plab
!       write(*,*) 'p tuple ord', p_tuple_ord%plab
      
      contrib_cache_already = p_tuple_compare(p_tmp_ord, p_tuple_ord)
      
!       write(*,*) 'Match?', contrib_cache_already

      if (contrib_cache_already) then
         
         if (num_p_tuples > 1) then
         
! write(*,*) 'match ,comparing with np > 1'
! write(*,*) next_element%contribs_outer%num_dmat


            if (present(n_rule)) then
            
               contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
                                       num_p_tuples - 1, p_tuples(2:num_p_tuples), n_rule=n_rule)
            
            else

               contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
                                       num_p_tuples - 1, p_tuples(2:num_p_tuples))
            
            end if

            
         else
         
            if (present(n_rule)) then
            
               contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
                                       0, (/emptypert/), n_rule=n_rule)
            
            else
            
               contrib_cache_already = contrib_cache_already_outer(next_element%contribs_outer, &
                                       0, (/emptypert/))
            
            end if            
            
                      
         end if
      
      end if

      if (next_element%last .eqv. .TRUE.) then
         passedlast = passedlast + 1
      end if

   end do

!    write(*,*) 'End result', contrib_cache_already
   
 end function
 

 ! NOT DONE WITH THIS FUNCTION
 
 function contrib_cache_already_outer(current_element, num_dmat, p_tuples_outer, n_rule)

   implicit none

   logical :: contrib_cache_already_outer
   integer :: passedlast, passedlast_outer, num_dmat, i
   integer, optional :: n_rule
   type(contrib_cache_outer), target :: current_element
   type(contrib_cache_outer), pointer :: next_element
   type(p_tuple), dimension(num_dmat) :: p_tuples_outer, p_tuples_ord
   
   next_element => current_element
   passedlast = 0
   p_tuples_ord = p_tuples_standardorder(num_dmat, p_tuples_outer)
   contrib_cache_already_outer = .FALSE.

   ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
   ! COULD THIS BE DONE IN ANOTHER WAY?
   do while ((passedlast < 2) .AND. (contrib_cache_already_outer .eqv. .FALSE.))

   ! Make sure the tuples that need to be ordered are ordered
   
!    write(*,*) 'curr element is last', next_element%last
   
      next_element => contrib_cache_outer_next_element(next_element)
      
!       write(*,*) 'next element is last', next_element%last

      if (next_element%num_dmat == num_dmat .AND. .NOT.(next_element%dummy_entry)) then
      
         contrib_cache_already_outer = p_tuples_compare(num_dmat, &
                                       next_element%p_tuples, p_tuples_ord)

                                    
         if (present(n_rule)) then
      
            contrib_cache_already_outer = contrib_cache_already_outer .AND. (n_rule == next_element%n_rule)

         end if
         
                                 
      end if

      if (next_element%last .eqv. .TRUE.) then
         passedlast = passedlast + 1
      end if

   end do
!  write(*,*) 'Was it found?', contrib_cache_already_outer

 end function

   ! NOT DONE WITH THIS FUNCTION
 
 function contrib_cache_already_inner(current_element, p_inner)

   implicit none

   logical :: contrib_cache_already_inner
   integer :: passedlast,  i
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element
   type(p_tuple) :: p_inner, p_tuple_ord
   
   next_element => current_element
   passedlast = 0
   p_tuple_ord = p_tuple_standardorder(p_inner)
   contrib_cache_already_inner = .FALSE.

   ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
   ! COULD THIS BE DONE IN ANOTHER WAY?
   do while ((passedlast < 2) .AND. (contrib_cache_already_inner .eqv. .FALSE.))

   ! Make sure the tuples that need to be ordered are ordered
   
      next_element => contrib_cache_next_element(next_element)
     
      contrib_cache_already_inner = p_tuple_compare(next_element%p_inner, p_tuple_ord)

      if (next_element%last .eqv. .TRUE.) then
         passedlast = passedlast + 1
      end if

   end do


 end function
 
    
 ! Missing contents (to be taken from non-outer routine below)
 
 subroutine contrib_cache_getdata_outer(cache, num_p_tuples, p_tuples, &
            from_inner, contrib_size, ind_len, ind_unsorted, hard_offset, mat, mat_sing, &
            scal, n_rule)

! Proposed: Puts all data into return array if prop or specified matrix if mat
            
! ind_unsorted, mat, prop are optional            
            
   implicit none
   logical :: found, from_inner
   integer :: i, j, k, first, last, passedlast, num_p_tuples, &
              total_num_perturbations, pr_offset, cache_offset, &
              merged_triang_size, merged_nblks, inner_rm, res_offset, &
              ind_len, nblks, offset, cache_hard_offset
   integer, optional :: hard_offset, n_rule
   integer :: contrib_size
   integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contrib, & 
                                         p_tuples_dimensions, &
                                         p_tuples_dimensions_cacheorder, &
                                         pids_merged_pert, translated_index
   integer, allocatable, dimension(:) :: blk_sizes_merged
   integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
   integer, allocatable, dimension(:,:) :: indices, blk_sizes
   integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
   integer, allocatable, dimension(:,:) :: blk_info
   integer, allocatable, dimension(:) :: blk_sizes_sing
   type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_ord
   type(p_tuple), allocatable, dimension(:) :: p_tuples_srch, p_tuples_srch_ord
   type(p_tuple) :: merged_p_tuple
   integer, dimension(ind_len), optional ::  ind_unsorted
   type(contrib_cache_outer), target :: cache
   type(contrib_cache_outer), pointer :: next_element_outer
   type(QcMat), optional, dimension(contrib_size) :: mat
   type(QcMat), optional :: mat_sing
   complex(8), optional, dimension(contrib_size) :: scal

 
   
   
   
   
   if (present(hard_offset)) then
   
      cache_hard_offset = hard_offset
   
   else
   
      cache_hard_offset = 0
   
   end if
   
   
!       write(*,*) 'b0'
   if (from_inner) then
      
      allocate(p_tuples_srch(num_p_tuples - 1))
      allocate(p_tuples_srch_ord(num_p_tuples - 1))
      inner_rm = 1

      do i = 1, num_p_tuples - 1
         call p1_cloneto_p2(p_tuples(i + 1), p_tuples_srch(i))
      end do
  
   else
   
      allocate(p_tuples_srch(num_p_tuples))
      allocate(p_tuples_srch_ord(num_p_tuples))
      inner_rm = 0

      do i = 1, num_p_tuples
         call p1_cloneto_p2(p_tuples(i), p_tuples_srch(i))
      end do
   
   end if

   
   p_tuples_srch_ord = p_tuples_standardorder(size(p_tuples_srch), p_tuples_srch)
   
   p_tuples_ord = p_tuples_standardorder(num_p_tuples, p_tuples)
  
      next_element_outer => cache
      passedlast = 0
      found = .FALSE.

   do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

!       write(*,*) 'searching'
   
      next_element_outer => contrib_cache_outer_next_element(next_element_outer)

!       write(*,*) 'num  - inner', num_p_tuples - inner_rm
!       write(*,*) 'len outer', size(next_element_outer%p_tuples)
!        if (size(next_element_outer%p_tuples) == 1) then
!       write(*,*) 'outer plab', next_element_outer%p_tuples(1)%plab
!       end if
!       write(*,*) 'len srch', size(p_tuples_srch_ord)
      
!             if (size(next_element_outer%p_tuples) == 1) then
!          write(*,*) next_element_outer%p_tuples(1)%plab
         if ((size(p_tuples_srch_ord) == 0) .AND. &
             (size(next_element_outer%p_tuples) == 1) ) then
             
!             if ((next_element_outer%p_tuples(1)%plab(1) == 'NUTN')) then
            if ((next_element_outer%p_tuples(1)%npert == 0)) then
            
                found = .TRUE.
            
            end if
            
         end if
!          write(*,*) 'outer case'
      
      
!       end if
      
      
      if (size(next_element_outer%p_tuples) == size(p_tuples_srch_ord)) then
         found = p_tuples_compare(num_p_tuples - inner_rm, next_element_outer%p_tuples, &
                                  p_tuples_srch_ord)
      end if   
                               
                                  
      if (present(n_rule)) then
      
         found = found .AND. (n_rule == next_element_outer%n_rule)
      
      end if

      if (next_element_outer%last) then
         passedlast = passedlast + 1
      end if

   end do
   
   if (.NOT.(found)) then
   
      write(*,*) 'ERROR: Element not found'
      stop
   
   end if
   
   if (present(mat) .OR. present(scal)) then
   
      if (found) then
      
!         write(*,*) 'THIS IS THE CACHE DATA'
!     write(*,*) next_element_outer%data_scal
!      write(*,*)  ' ' 
   
         

         total_num_perturbations = 0

         do i = 1, num_p_tuples
            total_num_perturbations = total_num_perturbations + p_tuples(i)%npert
         end do


         if (p_tuples(1)%npert == 0) then

            call p1_cloneto_p2(p_tuples(2), merged_p_tuple)

            do i = 3, num_p_tuples

               call p1_merge_p2(merged_p_tuple, p_tuples(i), merged_p_tuple)
               
            end do
            

         else
            
            call p1_cloneto_p2(p_tuples(1), merged_p_tuple) 
            
            if (num_p_tuples > 1) then
            
               if (p_tuples(2)%npert > 0) then
            
                  do i = 2, num_p_tuples

                     call p1_merge_p2(merged_p_tuple, p_tuples(i), merged_p_tuple)

                  end do
               
               end if
            
            end if
     
         end if
         
!       write(*,*) 'b1'

         merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
         merged_nblks = get_num_blks(merged_p_tuple)

         allocate(merged_blk_info(1,merged_nblks, 3))

         merged_blk_info(1,:,:) = get_blk_info(merged_nblks, merged_p_tuple)
         merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

         allocate(blk_sizes(num_p_tuples, total_num_perturbations))
         allocate(blk_sizes_merged(total_num_perturbations))

         blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
         merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

         allocate(indices(contrib_size, total_num_perturbations))

         call make_triangulated_indices(merged_nblks, merged_blk_info, & 
              merged_triang_size, indices)
!       write(*,*) 'b2'
         allocate(translated_index(total_num_perturbations))
         allocate(pids_in_cache(total_num_perturbations))
         allocate(pids_current_contrib(total_num_perturbations))
         allocate(pids_merged_pert(total_num_perturbations))
         allocate(p_tuples_dimensions(total_num_perturbations))
         allocate(p_tuples_dimensions_cacheorder(total_num_perturbations))

         k = 1

         do i = 1, num_p_tuples
            do j = 1, p_tuples(i)%npert
                pids_current_contrib(k) = p_tuples_ord(i)%pid(j)
               k = k + 1
            end do
         end do

         do i = 1, total_num_perturbations
             pids_merged_pert(i) = merged_p_tuple%pid(i)
         end do

         p_tuples_dimensions = get_ncarray(total_num_perturbations, num_p_tuples, &
                                           p_tuples_ord)

         allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
!       write(*,*) 'b3'
         do i = 1, num_p_tuples

            nfields(i) = p_tuples_ord(i)%npert
            nblks_tuple(i) = get_num_blks(p_tuples_ord(i))
            blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples_ord(i))
            blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                        blks_tuple_info(i,1:nblks_tuple(i),:))
            blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
            blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

         end do
!       write(*,*) 'b4'
         do i = 1, size(indices, 1)

            res_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
            (/merged_nblks/), & 
            (/total_num_perturbations/), (/merged_blk_info/), blk_sizes_merged, &
            (/merged_triang_size/), &
            (/indices(i, : )/) )

            do j = 1, total_num_perturbations
               translated_index(j) = indices(i,pids_current_contrib(j))
            end do

         if (p_tuples(1)%npert == 0) then

            cache_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
               total_num_perturbations, nblks_tuple(2:num_p_tuples), & 
               nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :),&
               blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
               translated_index)

         else
            
            if (num_p_tuples > 1) then
            
               if (p_tuples(2)%npert > 0) then
            
                  cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                  total_num_perturbations, nblks_tuple, & 
                  nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)

               else
               
!                   write(*,*) 'getting cache offset'
!                   write(*,*) 'total_num_perturbations', total_num_perturbations
!                   write(*,*) 'nblks_tuple', nblks_tuple(1)
!                   write(*,*) 'nfields', nfields(1)
!                   write(*,*) 'blks_tuple_info', blks_tuple_info(1,:,:)
!                   write(*,*) 'blk_sizes', blk_sizes(1,:)
!                   write(*,*) 'blks_tuple_triang_size', blks_tuple_triang_size(1)
!                   write(*,*) 'translated_index', translated_index
            
                  cache_offset = get_triang_blks_tuple_offset(1, &
                  total_num_perturbations, nblks_tuple(1), & 
                  nfields(1), blks_tuple_info(1,:,:), blk_sizes(1,:), &
                  blks_tuple_triang_size(1), translated_index)
               
!                   write(*,*) 'got cache offset'
               
               end if
            
            else
            
               cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
               total_num_perturbations, nblks_tuple, & 
               nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)
            
               
            end if
     
         end if
            

!       write(*,*) 'b5 i', i
!       write(*,*) 'data scal len', size(next_element_outer%data_scal)
            if (present(mat)) then
            
!                write(*,*) 'res offset', res_offset
!                write(*,*) 'res size', size(mat)
!                write(*,*) 'cache offset', cache_offset
!                write(*,*) 'cache size', size(next_element_outer%data_mat)
               
            
               call QcMatAEqB(mat(res_offset), next_element_outer%data_mat(cache_offset + cache_hard_offset))

               
            else if (present(scal)) then

            
!                if (i < 24) then
!                
!                   write(*,*) 'Indices', indices(i,:)
!                   write(*,*) 'Results offset', res_offset
!                   write(*,*) 'Cache offset', cache_offset
!              
!                
!                end if
            
               scal(res_offset) = &
               scal(res_offset) + &
               next_element_outer%data_scal(cache_offset + cache_hard_offset)              
            
            end if

         end do


         
     
      else

         write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Element not found'

      end if

   else if (present(mat_sing)) then
   
    if (p_tuples_srch(1)%npert > 0) then

       nblks = get_num_blks(p_tuples_srch(1))
       
       allocate(blk_sizes_sing(nblks))
       allocate(blk_info(nblks, 3))
       blk_info = get_blk_info(nblks, p_tuples_srch(1))
       blk_sizes_sing = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

!        write(*,*) ' p_tuples_srch(1)%npert',p_tuples_srch(1)%npert
!        write(*,*) ' nblks',nblks
!        write(*,*) ' blk_info ', blk_info
!        write(*,*) ' ind_unsorted ', ind_unsorted
       
       
       call sort_triangulated_indices(p_tuples_srch(1)%npert, nblks, &
                                         blk_info, ind_unsorted)

       offset = get_triang_blks_offset(nblks, p_tuples_srch(1)%npert, &
                                       blk_info, blk_sizes_sing, ind_unsorted)

       deallocate(blk_sizes_sing)
       deallocate(blk_info)

    else

       offset = 1

    end if
   
      if (found) then

!          write(*,*) 'cache element retrieval, offset, hard', offset, cache_hard_offset
!          write(*,*) 'data mat len', size(next_element_outer%data_mat)
!          write(*,*) 'length of mat', size(next_element_outer%data_mat), next_element_outer%last, &
!          next_element_outer%p_tuples(1)%plab, next_element_outer%p_tuples(1)%freq
!          write(*,*) 'length of mat', size(next_element_outer%next%data_mat), &
!          next_element_outer%next%last, next_element_outer%next%p_tuples(1)%plab
!          write(*,*) 'length of mat', size(next_element_outer%next%next%data_mat), &
!          next_element_outer%next%next%last, next_element_outer%next%next%p_tuples(1)%plab
         
         
         
      
         call QcMatAEqB(mat_sing, next_element_outer%data_mat(offset + cache_hard_offset))
         
!          if (offset == 4) then
!            write(*,*) 'printing matrix'
!            j = QcMatWrite_f(mat_sing, 'mat_sing', ASCII_VIEW)
!            
!          end if

      else

         write(*,*) 'Failed to retrieve data in sdf_getdata: Element not found'

      end if
   
    end if
!          write(*,*) 'b6'
 end subroutine
 
! Missing variable declaration and allocations, otherwise OK

 ! Assumes that p_tuples is in standard order
 subroutine contrib_cache_getdata(cache, num_p_tuples, p_tuples, contrib_size, &
                                  ind_len, ind_unsorted, hard_offset, mat, mat_sing, scal, n_rule)

   implicit none

   logical :: found
   integer :: i, j, k, first, last, passedlast, num_p_tuples, &
              contrib_size, total_num_perturbations, pr_offset, cache_offset, &
              merged_triang_size, merged_nblks, res_offset, ind_len, &
              cache_hard_offset
   integer, optional :: hard_offset, n_rule
   integer, optional, dimension(ind_len) :: ind_unsorted
   type(contrib_cache), target :: cache
   type(contrib_cache), pointer :: next_element
   type(contrib_cache_outer), pointer :: next_element_outer
   type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_ord
   type(p_tuple) :: merged_p_tuple
   type(QcMat), optional, dimension(contrib_size) :: mat
   type(QcMat), optional :: mat_sing
   complex(8), optional, dimension(contrib_size) :: scal

   if (present(hard_offset)) then
   
      cache_hard_offset = hard_offset
   
   else
   
      cache_hard_offset = 0
   
   end if
   

   next_element => cache
   passedlast = 0
   found = .FALSE.


   p_tuples_ord = p_tuples_standardorder(num_p_tuples, p_tuples)

   do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

      next_element => contrib_cache_next_element(next_element)

      found = p_tuple_compare(next_element%p_inner, p_tuples_ord(1))

      if (next_element%last) then
         passedlast = passedlast + 1
      end if

   end do

   
   if (found) then
   
      if (present(mat)) then

         if (present(n_rule)) then
             
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, hard_offset, &
                 mat=mat, n_rule=n_rule)
       
         else
         
!             write(*,*) 'outer mat case b'
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, hard_offset, mat=mat)
         
         end if
      
         

              
      else if (present(mat_sing)) then
      
         if (present(n_rule)) then
       
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 mat_sing=mat_sing, n_rule=n_rule)
       
         else
         
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 mat_sing=mat_sing)
                 
         end if
            
         
              
      else if (present(scal)) then
      
         if (present(n_rule)) then
         
!             write(*,*) 'getting scal nrule'
         
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 scal=scal, n_rule=n_rule)       
       
         else
         
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 scal=scal)
         
         end if
      

         
         
      end if
      
   else

      write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Inner element not found'

   end if


 end subroutine


 
 
 
 

  
  ! END NEW 2014
  
  
end module
