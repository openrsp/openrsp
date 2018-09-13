! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and functions related to caching of
! various contributions obtained during the calculation

module rsp_property_caching

 use rsp_field_tuple
 use rsp_indices_and_addressing
 use qcmatrix_f

 implicit none
 
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
  
 ! Define contrib cache datatypes
 
 ! "Outer" type: can be "independent" or attached to inner
 type contrib_cache_outer

   type(contrib_cache_outer), pointer :: next
   logical :: last
   logical :: dummy_entry
   ! Number of chain rule applications
   integer :: num_dmat
   ! Perturbation tuples
   type(p_tuple), allocatable, dimension(:) :: p_tuples

   ! Contribution type for two-factor terms: 1: Only Pulay n, 
   ! 3: Only Lagrange, 4: Both Pulay and Lagrange
   integer :: contrib_type = 0
   
   ! Choice of n rule for contribution (for use in two-factor terms)
   integer :: n_rule = 0
   ! Size of contribution data
   integer :: contrib_size = 0
   ! Perturbation block information
   integer, allocatable, dimension(:) :: nblks_tuple
   integer, allocatable, dimension(:,:) :: blk_sizes
   integer, allocatable, dimension(:,:) :: indices
   integer, allocatable, dimension(:,:,:) :: blks_tuple_info
   integer, allocatable, dimension(:) :: blks_tuple_triang_size
   ! Matrix data (if used)
   type(QcMat), allocatable, dimension(:) :: data_mat 
   ! Scalar data (if used)
   complex(8), allocatable, dimension(:) :: data_scal

 end type 

 ! "Inner" type: For use in e.g. perturbed Fock and energy-type terms
 ! Attaches to one or more outer cache instances
 type contrib_cache

   type(contrib_cache), pointer :: next
   logical :: last
   ! Perturbation tuple
   type(p_tuple) :: p_inner

   ! Number of outer cache instances attached to this inner
   integer :: num_outer
   ! Perturbation block/indices information
   integer :: nblks
   integer, allocatable, dimension(:) :: blk_sizes
   integer :: blks_triang_size
   integer, allocatable, dimension(:,:) :: blk_info
   integer, allocatable, dimension(:,:) :: indices
         
   ! Pointer to attached outer cache instances
   type(contrib_cache_outer), pointer :: contribs_outer

 end type 

 
 contains

 
  
  
  
  
 
 
    
  ! Initialize progress/restarting framework if dictated
  ! by restart flag r_flag
  subroutine prog_init(rs_info, r_flag)
  
    implicit none
    
    integer, dimension(3) :: rs_info
    integer :: r_flag
    logical :: r_exist
    
    if (r_flag == 3) then
    
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
       
    end if
    
  end subroutine
     
  ! Check progress counter against restart checkpoint
  ! Return true if checkpoint not passed, false if passed
  ! If optional argument 'lvl' specified, only check progress at that level
  ! Currently, only 'lvl' style check implemented
  function rs_check(prog_info, rs_info, r_flag, lvl)
  
    implicit none
    
    logical :: rs_check, rs_past
    integer, dimension(3) :: prog_info, rs_info
    integer :: i, r_flag
    integer, optional :: lvl
    
    rs_check = .FALSE.
    rs_past = .FALSE.
    
    if (r_flag == 0) then
    
       rs_check = .FALSE.
    
    else if (r_flag == 3) then
    
       if (present(lvl)) then
    
          do i = 1, lvl
       
             if (prog_info(i) > rs_info(i)) then
          
                rs_past = .TRUE.
             
             end if
               
             if (.NOT.(rs_past)) then
            
                if (rs_info(i) > prog_info(i)) then
          
                   rs_check = .TRUE.
                   
                end if
             
             end if
        
          end do   
    
       else
    
          write(*,*) 'ERROR: lvl keyword must be present for restart check'
           
       end if
   
    end if
    
  end function
  
  
  ! Increment calculation progress counter prog_info at level 'lvl'
  ! and store in file if dictated by restarting setup flag r_flag 
  subroutine prog_incr(prog_info, r_flag, lvl)
    
    implicit none
    
    integer, dimension(3) :: prog_info
    integer :: lvl, i, r_flag
        
    prog_info(lvl) = prog_info(lvl) + 1
    
    do i = lvl + 1, 3
    
       prog_info(i) = 0
    
    end do
    
    if (r_flag == 3) then
    
       open(unit=260, file='OPENRSP_RESTART', status='replace', action='write') 
       write(260,*) prog_info(1)
       write(260,*) prog_info(2)
       write(260,*) prog_info(3)
       close(260)
       
    end if
    
  end subroutine
 
 ! Store array of matrices or scalars in file fname
 subroutine mat_scal_store(array_size, fname, r_flag, mat, scal, start_pos)
 
   implicit none
   
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
   integer :: array_size, funit, i, k, first_i, r_flag
   integer, optional :: start_pos
   complex(8), dimension(array_size), optional :: scal
   type(QcMat), dimension(array_size), optional :: mat
 
   if (r_flag == 0) then
    
      ! "No restarting or storing" case
    
   elseif (r_flag == 3) then
 
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
       
    
         do i = 1, array_size
       
            write(str_fid, fid_fmt) i + first_i - 1
       
            k = QcMatWrite_f(mat(i), trim(adjustl(fname)) // '_MAT_' // &
                             trim(str_fid) // '.DAT', BINARY_VIEW)
       
         end do
    
      end if
    
   end if
 
 end subroutine
 
 ! Retrieve array of matrices or scalars from file fname
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
 
 ! Store contribution cache structure (including outer attached instances) in file fname
 subroutine contrib_cache_store(cache, r_flag, fname)
 
   implicit none
   
   logical :: termination
   integer :: num_entries, i, j, k, funit, mat_acc, r_flag
   character(*) :: fname
  
   type(contrib_cache), target :: cache
   type(contrib_cache), pointer :: cache_next
   type(contrib_cache_outer), pointer :: outer_next
 
   if (r_flag == 0) then
   
      ! "No restarting or storing" case
      
   else if (r_flag == 3) then
 
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
      
         ! MaR: New code for residue handling
      
         write(funit) cache_next%p_inner%do_residues
         write(funit) cache_next%p_inner%n_pert_res_max
         write(funit) cache_next%p_inner%n_states
      
         if (cache_next%p_inner%do_residues > 0) then
      
            if (allocated(cache_next%p_inner%states)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_inner%states)
         
               do j = 1, size(cache_next%p_inner%states)
         
                  write(funit) cache_next%p_inner%states(j)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
            if (allocated(cache_next%p_inner%exenerg)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_inner%exenerg)
         
               do j = 1, size(cache_next%p_inner%exenerg)
         
                  write(funit) cache_next%p_inner%exenerg(j)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         
            if (allocated(cache_next%p_inner%part_of_residue)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_inner%part_of_residue, 1)
               write(funit) size(cache_next%p_inner%part_of_residue, 2)
         
               do j = 1, size(cache_next%p_inner%part_of_residue, 1)
            
                  do k = 1, size(cache_next%p_inner%part_of_residue, 2)
         
                     write(funit) cache_next%p_inner%part_of_residue(j, k)
                  
                  end do
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         end if
      
      
         ! End new
      
      
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
 
   end if
 
 end subroutine
 
 ! Retrieve contribution cache structure (including outer attached instances) from file fname
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
    
   cache_next => cache
   
   call contrib_cache_read_one_inner(cache_next, fname, funit)
   
   allocate(cache_next%contribs_outer)
   outer_next => cache_next%contribs_outer
   call contrib_cache_outer_retrieve(outer_next, fname, .TRUE., funit)
   
   cache_orig => cache_next
   
   ! Get all remaining elements
   do i = 1, num_entries - 1
   
      call contrib_cache_retrieve_one_inner(cache_next%next, fname, funit)
      
      cache_next => cache_next%next
      allocate(cache_next%contribs_outer)
      outer_next => cache_next%contribs_outer
         
      call contrib_cache_outer_retrieve(outer_next, fname, .TRUE., funit)
      
   end do
   
   ! Bite the tail
   cache_next%next => cache_orig
   
   
   close(funit)
 
 end subroutine
 
 ! Wrapper: Get one inner contribution cache instance from file fname
 subroutine contrib_cache_retrieve_one_inner(cache_new, fname, funit)
 
    implicit none
 
    integer :: funit 
    character(*) :: fname
    type(contrib_cache), pointer :: cache_new
 
    allocate(cache_new)
 
    call contrib_cache_read_one_inner(cache_new, fname, funit)

 end subroutine
 
 ! Read one inner contribution cache instance from file fname
 subroutine contrib_cache_read_one_inner(cache, fname, funit)
 
   implicit none
   
   integer :: num_entries, i, j, funit
   integer :: size_i, size_j
   character(*) :: fname
   logical :: alloc_indic
  
   type(contrib_cache), pointer :: cache
   
   read(funit) cache%last
   read(funit) cache%p_inner%npert
   
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
   
   ! MaR: New code for residue handling
      
      read(funit) cache%p_inner%do_residues
      read(funit) cache%p_inner%n_pert_res_max
      read(funit) cache%p_inner%n_states
      
      if (cache%p_inner%do_residues > 0) then
      
         read(funit) alloc_indic
       
         if (alloc_indic) then
         
            read(funit) size_i
            allocate(cache%p_inner%states(size_i))
         
            do i = 1, size_i
         
               read(funit) cache%p_inner%states(i)
               
            end do
            
         end if
         
         read(funit) alloc_indic
       
         if (alloc_indic) then
         
            read(funit) size_i
            allocate(cache%p_inner%exenerg(size_i))
         
            do i = 1, size_i
         
               read(funit) cache%p_inner%exenerg(i)
               
            end do
            
         end if
         
         read(funit) alloc_indic
       
         if (alloc_indic) then
         
            read(funit) size_i
            read(funit) size_j
            allocate(cache%p_inner%part_of_residue(size_i, size_j))
         
            do i = 1, size_i
            
               do j = 1, size_j
         
                  read(funit) cache%p_inner%part_of_residue(i, j)
               
               end do
               
            end do
            
         end if
         
      end if
      
      
      ! End new
   
   
   
      
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
 
 ! Wrapper: Store outer type contribution cache element in file fname
 subroutine contrib_cache_outer_store(cache, fname, r_flag, mat_acc_in)
 
   implicit none
   
   type(contrib_cache_outer) :: cache
   character(*) :: fname
   integer, optional :: mat_acc_in
   integer :: funit, mat_acc, r_flag
   
   if (r_flag == 0) then
   
      ! "No restarting or storing" case
      
   else if (r_flag == 3) then
   
      mat_acc = 0
      if (present(mat_acc_in)) then
         mat_acc = mat_acc_in
      end if   
   
      funit = 260
 
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='replace', action='write')
        
      call contrib_cache_outer_put(cache, fname, funit, mat_acc_in=mat_acc)
        
        
      close(funit)
 
   end if
 
 
 end subroutine
 
 ! Write one outer contribution cache element to file fname
 ! Assumes file 'fname' is already opened with I/O unit 'funit'
 subroutine contrib_cache_outer_put(cache, fname, funit, mat_acc_in)
 
   implicit none
   
   logical :: termination
   integer :: num_entries, i, j, k, m, mat_scal_none, funit, mat_acc
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
         
         do k = 1, cache_next%p_tuples(j)%npert
         
            write(funit) cache_next%p_tuples(j)%pdim(k)
            write(funit) cache_next%p_tuples(j)%plab(k)
            write(funit) cache_next%p_tuples(j)%pid(k)
            write(funit) cache_next%p_tuples(j)%freq(k)
         
         end do
         
      
      
         ! MaR: New code for residue handling
      
         write(funit) cache_next%p_tuples(j)%do_residues
         write(funit) cache_next%p_tuples(j)%n_pert_res_max
         write(funit) cache_next%p_tuples(j)%n_states
      
         if (cache_next%p_tuples(j)%do_residues > 0) then
      
            if (allocated(cache_next%p_tuples(j)%states)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_tuples(j)%states)
         
               do k = 1, size(cache_next%p_tuples(j)%states)
         
                  write(funit) cache_next%p_tuples(j)%states(k)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
            if (allocated(cache_next%p_tuples(j)%exenerg)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_tuples(j)%exenerg)
         
               do k = 1, size(cache_next%p_tuples(j)%exenerg)
         
                  write(funit) cache_next%p_tuples(j)%exenerg(k)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         
            if (allocated(cache_next%p_tuples(j)%part_of_residue)) then
         
               write(funit) .TRUE.
               write(funit) size(cache_next%p_tuples(j)%part_of_residue, 1)
               write(funit) size(cache_next%p_tuples(j)%part_of_residue, 2)
         
               do k = 1, size(cache_next%p_tuples(j)%part_of_residue, 1)
            
                  do m = 1, size(cache_next%p_tuples(j)%part_of_residue, 2)
         
                     write(funit) cache_next%p_tuples(j)%part_of_residue(k, m)
                  
                  end do
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         end if
      
      
      ! End new
      
      end do
      
   
      write(funit) cache_next%last
      write(funit) cache_next%dummy_entry
      write(funit) cache_next%num_dmat
      write(funit) cache_next%contrib_type
      write(funit) cache_next%n_rule
      write(funit) cache_next%contrib_size
      
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
 
 ! Wrapper 1: Get outer contribution cache instance from file fname
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
   
      inquire(file=fname // '.DAT', exist=cache_ext)
      
      if (.NOT.(cache_ext)) then
   
         write(*,*) 'ERROR: The expected cache file ', trim(adjustl(fname)), ' does not exist'
         stop
   
      end if
      
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='old', action='read')
           
   end if
   
   
   
   read(funit) num_entries
   
   call contrib_cache_read_one_outer(cache, fname, funit, mat_acc_in=mat_acc)

   cache_next => cache
   cache_orig => cache
   
   do i = 1, num_entries - 1
   
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
 
 
 ! Wrapper 2: Get one outer cache element from file fname
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
 
 ! Read one outer contribution cache element from file fname
 subroutine contrib_cache_read_one_outer(cache, fname, funit, mat_acc_in)
 
   implicit none
   
   logical :: from_inner
   logical :: alloc_indic
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
   allocate(cache%p_tuples(size_i))
   
   do j = 1, size(cache%p_tuples)
      
      read(funit) cache%p_tuples(j)%npert
         
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
         
   
      
      ! MaR: New code for residue handling
      ! NOTE: Potential size_i conflict - likely no problem 
      ! but revisit if there are issues
      
      read(funit) cache%p_tuples(j)%do_residues
      read(funit) cache%p_tuples(j)%n_pert_res_max
      read(funit) cache%p_tuples(j)%n_states
      
      if (cache%p_tuples(j)%do_residues > 0) then
      
         read(funit) alloc_indic
       
         if (alloc_indic) then
         
            read(funit) size_i
            allocate(cache%p_tuples(j)%states(size_i))
         
            do i = 1, size_i
         
               read(funit) cache%p_tuples(j)%states(i)
               
            end do
            
         end if
         
         read(funit) alloc_indic
      
         if (alloc_indic) then
         
            read(funit) size_i
            allocate(cache%p_tuples(j)%exenerg(size_i))
         
            do i = 1, size_i
         
               read(funit) cache%p_tuples(j)%exenerg(i)
               
            end do
            
         end if
         
         read(funit) alloc_indic
       
         if (alloc_indic) then
         
            read(funit) size_i
            read(funit) size_j
            allocate(cache%p_tuples(j)%part_of_residue(size_i, size_j))
         
            do i = 1, size_i
            
               do k = 1, size_j
         
                  read(funit) cache%p_tuples(j)%part_of_residue(i, k)
               
               end do
               
            end do
            
         end if
         
      end if
      
      
      ! End new
      
   end do   
   
  
   read(funit) cache%last
   read(funit) cache%dummy_entry
   read(funit) cache%num_dmat
   read(funit) cache%contrib_type
   read(funit) cache%n_rule
   read(funit) cache%contrib_size
   
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
      
   else
   
      allocate(cache%indices(1,1))
      cache%indices(1,1) = 1
      
   
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
 
 
 ! Cycle outer instances attached to 'current_element' until specified
 ! element is reached
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

   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(current_element, p_tuples(1))) then
   
      next_element => current_element
      
      ! Skip to cache element for this inner
      do while (p_tuple_compare(next_element%p_inner, p_tuples(1)) .EQV. .FALSE.)

        next_element => next_element%next
         
      end do
      
      next_outer => next_element%contribs_outer
      
      passedlast = 0
      found = .FALSE.

      do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

         next_outer => contrib_cache_outer_next_element(next_outer)
         
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
   

 end subroutine
 
 ! Cycle contribution cache until first element reached
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

 ! Cycle outer contribution cache element until first element reached
 function contrib_cache_outer_cycle_first(current_element) result(next_element)

   implicit none

   type(contrib_cache_outer), target :: current_element
   type(contrib_cache_outer), pointer :: next_element
   
   next_element => current_element
   do while (next_element%last .eqv. .FALSE.)
   
          next_element => next_element%next
   end do

   next_element => next_element%next   

 end function
 
 ! Return next element in contribution cache linked list
 function contrib_cache_next_element(current_element) result(next_element)

   implicit none

   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element

   next_element => current_element%next

 end function

 ! Return next element in outer contribution cache linked list
 function contrib_cache_outer_next_element(current_element) result(next_element)

   implicit none

   type(contrib_cache_outer), target :: current_element
   type(contrib_cache_outer), pointer :: next_element

   next_element => current_element%next

 end function
 
  ! Allocate and set up new contribution cache element
  subroutine contrib_cache_allocate(current_element)

   implicit none

   type(contrib_cache), pointer :: current_element
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
 
 ! Allocate and set up new outer contribution cache element
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
 
 ! Initialize contribution cache (and associated outer contribution cache)
 ! as specified by perturbation tuples 'p_tuples'
 subroutine contrib_cache_initialize(new_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   integer :: num_p_tuples
   integer, optional :: n_rule
   type(contrib_cache) :: new_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples

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
  
   if (num_p_tuples > 1) then
   
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
   
      call contrib_cache_outer_allocate(new_element%contribs_outer)

      ! MaR: Passing empty perturbation tuple as dummy argument since there should
      ! be no outer perturbation tuples in this case; revisit if problems.
      
      if (present(n_rule)) then
         
         call contrib_cache_outer_add_element(new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                              (/get_emptypert()/), n_rule=n_rule)
         
      else
      
         call contrib_cache_outer_add_element(new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                              (/get_emptypert()/))
         
      end if
      
   end if
      
 end subroutine
 
 
 ! Initialize outer contribution cache as specified by perturbation tuples 'outer_p_tuples'
 ! 'unperturbed' flag signifies empty perturbation tuple (typically used for 'all inner' case)
 subroutine contrib_cache_outer_initialize(new_element, unperturbed, num_dmat, outer_p_tuples)

   implicit none

   logical :: unperturbed
   integer :: num_dmat, i, total_npert
   type(contrib_cache_outer) :: new_element
   type(p_tuple), dimension(num_dmat) :: outer_p_tuples

   new_element%num_dmat = num_dmat
   new_element%last = .TRUE.
   new_element%dummy_entry = .FALSE.
   
   if (.NOT.(unperturbed)) then
   
     allocate(new_element%p_tuples(num_dmat))
   
     do i = 1, num_dmat
   
        call p1_cloneto_p2(outer_p_tuples(i), new_element%p_tuples(i))
        
     end do

     total_npert = sum((/(outer_p_tuples(i)%npert, i = 1, num_dmat)/))
     
     if (total_npert > 0) then
     
        allocate(new_element%nblks_tuple(num_dmat))
        allocate(new_element%blk_sizes(num_dmat, total_npert))
        allocate(new_element%blks_tuple_info(num_dmat, total_npert, 3))
        allocate(new_element%blks_tuple_triang_size(num_dmat))
     
     
        new_element%nblks_tuple = (/(get_num_blks(outer_p_tuples(i)), i = 1, num_dmat)/)
   
        do i = 1, num_dmat
     
           new_element%blks_tuple_info(i, :, :) = get_blk_info(new_element%nblks_tuple(i), outer_p_tuples(i))
           new_element%blk_sizes(i, 1:new_element%nblks_tuple(i)) = &
           get_triangular_sizes(new_element%nblks_tuple(i), &
           new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),2), &
           new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),3))
           new_element%blks_tuple_triang_size(i) = get_triangulated_size(new_element%nblks_tuple(i), &
                                      new_element%blks_tuple_info(i, 1:new_element%nblks_tuple(i), :))
                                   
        end do
     
        allocate(new_element%indices(product(new_element%blks_tuple_triang_size), total_npert))
        call make_triangulated_tuples_indices(num_dmat, total_npert, new_element%nblks_tuple, &
             new_element%blks_tuple_info, new_element%blks_tuple_triang_size, new_element%indices)
     
     
     else
     
     ! Reserved if other treatment necessary for "no external perturbations" case
     
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
      
      ! MaR: Creating indices for unperturbed entry: Return here if problems arise
      
      allocate(new_element%indices(1,1))
      new_element%indices(1,1) = 1 
   
   end if
      
 end subroutine
   
   
 ! Add outer contribution cache element to existing linked list
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
   
   if (already) then
   
      next_element => curr_element
      passedlast = 0
      p_tuples_st_order = p_tuples_standardorder(num_dmat, outer_p_tuples)
      found_element = .FALSE.

      do while ((passedlast < 2) .AND. (found_element .eqv. .FALSE.))
      
         next_element => next_element%next
         
         if (next_element%num_dmat == num_dmat .AND. .NOT.(next_element%dummy_entry)) then
            found_element = p_tuples_compare(next_element%num_dmat, &
            p_tuples_standardorder(next_element%num_dmat, &
            next_element%p_tuples), p_tuples_st_order)
         end if
            
         if (next_element%last .eqv. .TRUE.) then
            passedlast = passedlast + 1
         end if

      end do
      
      if(present(data_mat)) then
      
         do i = 1, data_size
         
            call QcMatAEqB(next_element%data_mat(i),  data_mat(i))
           
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
         end do
   
      end if
   
      if (present(data_scal)) then
   
      end if
      
      if (present(n_rule)) then
      
         new_element%n_rule = n_rule
      
      end if

      next_element%last = .FALSE.
      new_element%next => next_element%next
      next_element%next => new_element
   
   end if
      
 end subroutine
 
 ! Add contribution cache element to existing linked list
 subroutine contrib_cache_add_element(current_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   logical :: empty_outer
   integer :: num_p_tuples, i
   integer, optional :: n_rule
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: new_element
   type(contrib_cache), pointer :: new_element_ptr
   type(contrib_cache), pointer :: next_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert

   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(current_element, p_tuples(1))) then
   
      next_element => current_element
   
      ! Skip to cache element for this inner
      do while (p_tuple_compare(next_element%p_inner, p_tuples(1)) .EQV. .FALSE.)

         next_element => next_element%next
          
      end do
      
      ! NOTE: AND CONDITION MAY ADD WRONG KIND OF CACHE ELEMENT, REVISIT IF ERROR
      
      empty_outer = .TRUE.
      
      if (num_p_tuples > 1) then
     
         if (p_tuples(2)%npert > 0) then
         
            empty_outer = .FALSE.
            next_element%num_outer = next_element%num_outer + 1
     
            if (present(n_rule)) then
               
               call contrib_cache_outer_add_element(next_element%contribs_outer, .FALSE., &
                    num_p_tuples - 1, p_tuples(2:num_p_tuples), n_rule=n_rule)
         
            else
         
               call contrib_cache_outer_add_element(next_element%contribs_outer, .FALSE., &
                    num_p_tuples - 1, p_tuples(2:num_p_tuples))
         
            end if
         
         end if
         
      end if
            
      if (empty_outer) then
     
         next_element%num_outer = next_element%num_outer + 1

         call empty_p_tuple(emptypert)
         
         ! MaR: WARNING: CHANGED num_p_tuples ARGUMENT TO 0 BELOW; MAY PRODUCE ERRORS ELSEWHERE
         ! UPDATE: Looks like it didn't make errors, return here if there are problems anyway

         if (present(n_rule)) then
         
            call contrib_cache_outer_add_element(next_element%contribs_outer, .TRUE., 1, (/emptypert/), n_rule=n_rule)
         
         else
            
            call contrib_cache_outer_add_element(next_element%contribs_outer, .TRUE., 0, (/emptypert/))
         
         end if
       
      end if


   ! Otherwise, add both inner and outer    
   else
   
      next_element => current_element
      allocate(new_element)
      
      if (present(n_rule)) then
         
         call contrib_cache_initialize(new_element, num_p_tuples, p_tuples, n_rule)   
         
      else
         
         call contrib_cache_initialize(new_element, num_p_tuples, p_tuples)   
         
      end if
         
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
 
 ! Check if element (inner + outer combination) as specified by 'p_tuples' (and possible 'n_rule')
 ! already exists in cache
 function contrib_cache_already(current_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   logical :: contrib_cache_already
   integer :: passedlast, passedlast_outer, num_p_tuples, i
   integer, optional :: n_rule
   type(contrib_cache), target :: current_element
   type(contrib_cache), pointer :: next_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert, p_tuple_ord, p_tmp_ord

   next_element => current_element
   passedlast = 0
   
   p_tuple_ord = p_tuple_standardorder(p_tuples(1))
   contrib_cache_already = .FALSE.

   call empty_p_tuple(emptypert)
   
   ! NOTE (MaR): WHILE LOOP POTENTIALLY NON-TERMINATING
   ! COULD THIS BE DONE IN ANOTHER WAY?
   do while ((passedlast < 2) .AND. (contrib_cache_already .eqv. .FALSE.))

      next_element => contrib_cache_next_element(next_element)
      p_tmp_ord = p_tuple_standardorder(next_element%p_inner)
      
      contrib_cache_already = p_tuple_compare(p_tmp_ord, p_tuple_ord)
      
      if (contrib_cache_already) then
         
         if (num_p_tuples > 1) then
         
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

 end function
 
 
 ! Check if element of outer contribution cache as specified by 'p_tuples_outer'
 ! (and possible 'n_rule') exists in linked list
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
   
      next_element => contrib_cache_outer_next_element(next_element)
      
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

 end function

 ! Check if contribution cache (inner only) element (as specified by 'p_inner') 
 ! already exists in linked list 
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
 
 ! Get data from outer contribution cache element    
 subroutine contrib_cache_getdata_outer(cache, num_p_tuples, p_tuples, &
            from_inner, contrib_size, ind_len, ind_unsorted, hard_offset, mat, mat_sing, &
            scal, n_rule)

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

   ! Offset for two-factor type terms
   if (present(hard_offset)) then
   
      cache_hard_offset = hard_offset
   
   else
   
      cache_hard_offset = 0
   
   end if
   
   ! Remove first tuple from tuple of perturbation tuples if this
   ! outer cache was attached to an inner cache
   if (from_inner) then
      
      allocate(p_tuples_srch(num_p_tuples - 1))
      allocate(p_tuples_srch_ord(num_p_tuples - 1))
      inner_rm = 1

      do i = 1, num_p_tuples - 1
         call p1_cloneto_p2(p_tuples(i + 1), p_tuples_srch(i))
      end do
  
   ! Otherwise, keep entire tuple of tuples
   else
   
      allocate(p_tuples_srch(num_p_tuples))
      allocate(p_tuples_srch_ord(num_p_tuples))
      inner_rm = 0

      do i = 1, num_p_tuples
         call p1_cloneto_p2(p_tuples(i), p_tuples_srch(i))
      end do
   
   end if

   
   ! Search to find cache element
   p_tuples_srch_ord = p_tuples_standardorder(size(p_tuples_srch), p_tuples_srch)
   
   p_tuples_ord = p_tuples_standardorder(num_p_tuples, p_tuples)
  
      next_element_outer => cache
      passedlast = 0
      found = .FALSE.

   do while ((passedlast < 2) .AND. (found .eqv. .FALSE.))

   
      next_element_outer => contrib_cache_outer_next_element(next_element_outer)

         if ((size(p_tuples_srch_ord) == 0) .AND. &
             (size(next_element_outer%p_tuples) == 1) ) then
             
            if ((next_element_outer%p_tuples(1)%npert == 0)) then
            
                found = .TRUE.
            
            end if
            
         end if
      
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
   
   ! If matrix (array) or scalar (array) data requested
   ! Typical for resp. lower-order Fock and response tensor calls
   if (present(mat) .OR. present(scal)) then
   
      if (found) then
      
         ! Merge perturbation tuples for indexing
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
         
         ! Get block information for indexing
         
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

         do i = 1, num_p_tuples

            nfields(i) = p_tuples_ord(i)%npert
            nblks_tuple(i) = get_num_blks(p_tuples_ord(i))
            blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples_ord(i))
            blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                        blks_tuple_info(i,1:nblks_tuple(i),:))
            blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
            blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

         end do

         ! Loop over indices, get appropriate offsets and put data in return array
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
               
                     cache_offset = get_triang_blks_tuple_offset(1, &
                     total_num_perturbations, nblks_tuple(1), & 
                     nfields(1), blks_tuple_info(1,:,:), blk_sizes(1,:), &
                     blks_tuple_triang_size(1), translated_index)
               
                  end if
            
               else
            
                  cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                  total_num_perturbations, nblks_tuple, & 
                  nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)
                        
               end if
     
            end if
            
            if (present(mat)) then
            
               call QcMatRAXPY(1.0d0, next_element_outer%data_mat(cache_offset + &
                               cache_hard_offset), mat(res_offset))
               
            else if (present(scal)) then
            
               !write(*,*) 'res, cache offset', res_offset, cache_offset + cache_hard_offset
               !write(*,*) 'val', next_element_outer%data_scal(cache_offset + cache_hard_offset)              

               scal(res_offset) = &
               scal(res_offset) + &
               next_element_outer%data_scal(cache_offset + cache_hard_offset)              
            
            end if

         end do

         deallocate(merged_blk_info)
         deallocate(blk_sizes)
         deallocate(blk_sizes_merged)
         deallocate(indices)
         deallocate(translated_index)
         deallocate(pids_in_cache)
         deallocate(pids_current_contrib)
         deallocate(pids_merged_pert)
         deallocate(p_tuples_dimensions)
         deallocate(p_tuples_dimensions_cacheorder)
         
         
      else

         write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Element not found'

      end if
      
      deallocate(p_tuples_srch)
      deallocate(p_tuples_srch_ord)
      

   ! If single matrix requested
   ! Typical for single 'get perturbed S, D, or F' calls
   else if (present(mat_sing)) then
   
    if (p_tuples_srch(1)%npert > 0) then

       ! Get block information for indexing
       nblks = get_num_blks(p_tuples_srch(1))
       
       allocate(blk_sizes_sing(nblks))
       allocate(blk_info(nblks, 3))
       blk_info = get_blk_info(nblks, p_tuples_srch(1))
       blk_sizes_sing = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

       call sort_triangulated_indices(p_tuples_srch(1)%npert, nblks, &
                                         blk_info, ind_unsorted)

       offset = get_triang_blks_offset(nblks, p_tuples_srch(1)%npert, &
                                       blk_info, blk_sizes_sing, ind_unsorted)

       deallocate(blk_sizes_sing)
       deallocate(blk_info)

    else

       ! If unperturbed, offset is 1
       offset = 1

    end if
   
      if (found) then
      
         call QcMatAEqB(mat_sing, next_element_outer%data_mat(offset + cache_hard_offset))
         
      else

         write(*,*) 'Failed to retrieve data in sdf_getdata: Element not found'

      end if
   
    end if

 end subroutine
 
 ! Get data from contribution cache
 ! Will search linked list for appropriate inner entry, then
 ! search associated outer contribution cache linked list
 ! and retrieve with contrib_cache_getdata_outer
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

   ! Hard offset typical for two-factor terms
   if (present(hard_offset)) then
   
      cache_hard_offset = hard_offset
   
   else
   
      cache_hard_offset = 0
   
   end if
   

   ! Search inner cache for entry
   
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

   ! If found, move to search (and retrieve) from outer cache
   ! Several cases depending on which type of data is requested,
   ! and whether or not (k,n) rule choice (specified as n in 'n_rule')
   ! is a discriminating factor in search (the latter typical for two-factor terms)
   if (found) then
   
      if (present(mat)) then

         if (present(n_rule)) then
             
            call contrib_cache_getdata_outer(next_element%contribs_outer, num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, hard_offset, &
                 mat=mat, n_rule=n_rule)
       
         else
         
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

end module
