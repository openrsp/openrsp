! Copyright 2012 Magnus Ringholm
!
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

! Contains routines and functions related to caching of
! various contributions obtained during the calculation

module rsp_property_caching

 use rsp_field_tuple
 use rsp_indices_and_addressing
 use qcmatrix_f

 implicit none
 
 public contrib_cache_initialize
 public contrib_cache_locate
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
 public contrib_cache_outer_dealloc_mat
 public mat_scal_store
 public mat_scal_retrieve
 public rs_check
 public prog_incr
 public prog_init
  
 ! Define contrib cache datatypes
 
 ! New "outer" type
 ! An instance of this is one cache entry to be used in array of such
 type contrib_cache_outer
 
   logical :: dummy_entry
   ! Number of chain rule applications
   integer :: num_dmat
   ! Perturbation tuples
   type(p_tuple), allocatable, dimension(:) :: p_tuples
   
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
 
 ! New "inner" type
 
 ! "Inner" type: For use in e.g. perturbed Fock and energy-type terms
 ! An instance of this is one cache entry to be used in array of such
 ! Attaches to one or more outer cache instances
 type contrib_cache

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
         
   ! Associated outer cache instances
   type(contrib_cache_outer), allocatable, dimension(:) :: contribs_outer

 end type 
 
 contains
     
  ! Initialize progress/restarting framework if dictated
  ! by restart flag r_flag
  subroutine prog_init(rs_info, r_flag)
  
    implicit none
    
    integer, dimension(3) :: rs_info
    integer :: r_flag, funit
    logical :: r_exist
    
    if (r_flag == 3) then
    
       inquire(file='OPENRSP_RESTART', exist=r_exist)
    
       if (r_exist) then
       
          open(newunit=funit, file='OPENRSP_RESTART', action='read') 
          read(funit,*) rs_info(1)
          read(funit,*) rs_info(2)
          read(funit,*) rs_info(3)
          close(funit)
    
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
    integer :: lvl, i, r_flag, funit
        
    prog_info(lvl) = prog_info(lvl) + 1
    
    do i = lvl + 1, 3
    
       prog_info(i) = 0
    
    end do
    
    if (r_flag == 3) then
    
       open(newunit=funit, file='OPENRSP_RESTART', status='replace', action='write') 
       write(funit,*) prog_info(1)
       write(funit,*) prog_info(2)
       write(funit,*) prog_info(3)
       close(funit)
       
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
 
      if (present(scal)) then
 
         open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
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
 
   if (present(scal)) then
 
      open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
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
 
 ! FIXME: Definitely ll-to-array chgs needed for this routine UPD: FIRST NEW VERSION
 
 ! Store contribution cache structure (including outer attached instances) in file fname
 subroutine contrib_cache_store(len_cache, cache, r_flag, fname)
 
   implicit none
   
   integer :: num_entries, i, j, k, funit, mat_acc, r_flag
   integer :: len_cache, n_mpt
   character(*) :: fname
  
   type(contrib_cache), dimension(len_cache) :: cache
 
   if (r_flag == 0) then
   
      ! "No restarting or storing" case
      
   else if (r_flag == 3) then

      ! FIXME: I may have rmd something around mat_acc which I shouldn't have,
      ! don't know but return to this later if problems

      mat_acc = 0
 
      open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='replace', action='write')
   
      ! Find how many unperturbed elements: They should be skippable when writing
      ! since they were just "initialization" elements
      n_mpt = 0
      
      num_entries = len_cache
      
      do i = 1, num_entries
   
         if (cache(i)%p_inner%npert == 0) then
         
            n_mpt = n_mpt + 1
           
         end if
   
      end do

      write(funit) num_entries - n_mpt
   
      ! Traverse and store
      do i = 1, num_entries
      
        if (cache(i)%p_inner%npert == 0) then
        
            cycle
        
        end if
   
         write(funit) cache(i)%p_inner%npert
         do j = 1, cache(i)%p_inner%npert
      
            write(funit) cache(i)%p_inner%pdim(j)
            write(funit) cache(i)%p_inner%plab(j)
            write(funit) cache(i)%p_inner%pid(j)
            write(funit) cache(i)%p_inner%freq(j)
      
         end do
      
         ! MaR: New code for residue handling

         write(funit) cache(i)%p_inner%do_residues
         write(funit) cache(i)%p_inner%n_pert_res_max
         write(funit) cache(i)%p_inner%n_states
      
         if (cache(i)%p_inner%do_residues > 0) then
      
            if (allocated(cache(i)%p_inner%states)) then

               write(funit) .TRUE.
               write(funit) size(cache(i)%p_inner%states)
         
               do j = 1, size(cache(i)%p_inner%states)
         
                  write(funit) cache(i)%p_inner%states(j)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
            if (allocated(cache(i)%p_inner%exenerg)) then

               write(funit) .TRUE.
               write(funit) size(cache(i)%p_inner%exenerg)
         
               do j = 1, size(cache(i)%p_inner%exenerg)
         
                  write(funit) cache(i)%p_inner%exenerg(j)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         
            if (allocated(cache(i)%p_inner%part_of_residue)) then
         
               write(funit) .TRUE.
               write(funit) size(cache(i)%p_inner%part_of_residue, 1)
               write(funit) size(cache(i)%p_inner%part_of_residue, 2)
         
               do j = 1, size(cache(i)%p_inner%part_of_residue, 1)
            
                  do k = 1, size(cache(i)%p_inner%part_of_residue, 2)

                     write(funit) cache(i)%p_inner%part_of_residue(j, k)
                  
                  end do
               
               end do
            
            else

               write(funit) .FALSE.
         
            end if
         
         end if
      
      
         ! End new
      
         write(funit) size(cache(i)%contribs_outer)
         
         write(funit) cache(i)%nblks
      
         write(funit) size(cache(i)%blk_sizes)
         write(funit) cache(i)%blk_sizes
      
         write(funit) cache(i)%blks_triang_size

         write(funit) size(cache(i)%blk_info, 1)
         write(funit) size(cache(i)%blk_info, 2)
         write(funit) cache(i)%blk_info
     
         write(funit) size(cache(i)%indices, 1)
         write(funit) size(cache(i)%indices, 2)
         write(funit) cache(i)%indices
      
         
         call contrib_cache_outer_put(size(cache(i)%contribs_outer), &
              cache(i)%contribs_outer, fname, funit, mat_acc_in=mat_acc)
   
   
      end do

      close(funit)

   end if
 
 end subroutine
 
 ! FIXME: Definitely ll-to-array chgs needed for this routine
 
 ! Retrieve contribution cache structure (including outer attached instances) from file fname
 subroutine contrib_cache_retrieve(cache, fname)
 
   implicit none
   
   logical :: cache_ext
   integer :: num_entries, i, j, funit
   integer :: dum, mat_acc
   character(*) :: fname
   
   type(contrib_cache), dimension(:), allocatable :: cache, cache_head_elem

   ! Old LL functionality: Rewrite together with store routine into array form   
   
   inquire(file=fname // '.DAT', exist=cache_ext)

   if (.NOT.(cache_ext)) then
   
      write(*,*) 'ERROR: The expected cache file ', trim(adjustl(fname)), ' does not exist'
      stop
   
   end if
 
   open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
        form='unformatted', status='old', action='read')
   
   read(funit) num_entries
   
   ! NOTE: I choose to always deallocate the cache if already allocated since anything else
   ! would tend to confuse things
   
   if (allocated(cache)) then
   
      deallocate(cache)
   
   end if 
   
   allocate(cache(num_entries + 1))
   
   call contrib_cache_allocate(cache_head_elem)
   
   cache(1) = cache_head_elem(1)
   
   mat_acc = 0
   
   do i = 1, num_entries
   
      call contrib_cache_read_one_inner(cache(i + 1), fname, funit)
      
      allocate(cache(i + 1)%contribs_outer(cache(i + 1)%num_outer))
    
      do j = 1, cache(i + 1)%num_outer
     
         call contrib_cache_read_one_outer(cache(i + 1)%contribs_outer(j), fname, &
              funit, mat_acc_in=mat_acc)
      
      end do
      
      ! Since I wrote it with one more outer element (maybe because of dummy)
      ! I now need to decrease the number of outer by one to get the "acceptable"
      ! number of outer entries which avoids indexing issues, due to idiosyncrasies 
      ! of the current setup
      cache(i + 1)%num_outer = cache(i + 1)%num_outer - 1
   
   end do

   close(funit)
 
 end subroutine
 
 ! FIXME: Definitely ll-to-array chgs needed for this routine
 
 ! Wrapper: Get one inner contribution cache instance from file fname  FIXME: UPD BUT INCOMPLETE
 subroutine contrib_cache_retrieve_one_inner(cache_new, fname, funit)
 
    implicit none
 
    integer :: funit 
    character(*) :: fname
    type(contrib_cache) :: cache_new

! Old code to be revisited
 
!     allocate(cache_new)
!  
!     call contrib_cache_read_one_inner(cache_new, fname, funit)

 end subroutine

 ! FIXME: Definitely ll-to-array chgs needed for this routine FIXME: UPD BUT MAY BE INCOMPLETE
 
 ! Read one inner contribution cache instance from file fname
 subroutine contrib_cache_read_one_inner(cache, fname, funit)
 
   implicit none
   
   integer :: num_entries, i, j, funit
   integer :: size_i, size_j
   character(*) :: fname
   logical :: alloc_indic
  
   type(contrib_cache) :: cache
  
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
 subroutine contrib_cache_outer_store(len_cache, cache, fname, r_flag, mat_acc_in)
 
   implicit none
   
   integer :: len_cache
   type(contrib_cache_outer), dimension(len_cache) :: cache
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
   
      open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='replace', action='write')
           
      write(funit) size(cache)
        
      call contrib_cache_outer_put(size(cache), cache, fname, funit, mat_acc_in=mat_acc)
        
        
      close(funit)
 
   end if
 
 
 end subroutine

 ! FIXME: Definitely ll-to-array chgs needed for this routine UPD: FIRST NEW VERSION
 
 ! Write one outer contribution cache element to file fname
 ! Assumes file 'fname' is already opened with I/O unit 'funit'
 subroutine contrib_cache_outer_put(len_cache, cache, fname, funit, mat_acc_in)
 
   implicit none
   
   logical :: termination
   integer :: num_entries, i, j, k, m, mat_scal_none, funit, mat_acc
   integer, optional :: mat_acc_in
   integer :: len_cache
   character(*) :: fname
   character(10) :: str_fid
   character(8) :: fid_fmt
  
   type(contrib_cache_outer), dimension(len_cache) :: cache

   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
   
   num_entries = len_cache
   

   
  
   ! Traverse and store
   do i = 1, num_entries
   
!       write(funit) len_cache
   
      write(funit) size(cache(i)%p_tuples)
      do j = 1, size(cache(i)%p_tuples)
      
         write(funit) cache(i)%p_tuples(j)%npert
         
         do k = 1, cache(i)%p_tuples(j)%npert
         
            write(funit) cache(i)%p_tuples(j)%pdim(k)
            write(funit) cache(i)%p_tuples(j)%plab(k)
            write(funit) cache(i)%p_tuples(j)%pid(k)
            write(funit) cache(i)%p_tuples(j)%freq(k)
         
         end do
         
         ! MaR: New code for residue handling
      
         write(funit) cache(i)%p_tuples(j)%do_residues
         write(funit) cache(i)%p_tuples(j)%n_pert_res_max
         write(funit) cache(i)%p_tuples(j)%n_states
      
         if (cache(i)%p_tuples(j)%do_residues > 0) then
      
            if (allocated(cache(i)%p_tuples(j)%states)) then
         
               write(funit) .TRUE.
               write(funit) size(cache(i)%p_tuples(j)%states)
         
               do k = 1, size(cache(i)%p_tuples(j)%states)
         
                  write(funit) cache(i)%p_tuples(j)%states(k)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
            if (allocated(cache(i)%p_tuples(j)%exenerg)) then
         
               write(funit) .TRUE.
               write(funit) size(cache(i)%p_tuples(j)%exenerg)
         
               do k = 1, size(cache(i)%p_tuples(j)%exenerg)
         
                  write(funit) cache(i)%p_tuples(j)%exenerg(k)
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         
            if (allocated(cache(i)%p_tuples(j)%part_of_residue)) then
         
               write(funit) .TRUE.
               write(funit) size(cache(i)%p_tuples(j)%part_of_residue, 1)
               write(funit) size(cache(i)%p_tuples(j)%part_of_residue, 2)
         
               do k = 1, size(cache(i)%p_tuples(j)%part_of_residue, 1)
            
                  do m = 1, size(cache(i)%p_tuples(j)%part_of_residue, 2)
         
                     write(funit) cache(i)%p_tuples(j)%part_of_residue(k, m)
                  
                  end do
               
               end do
            
            else
         
               write(funit) .FALSE.
         
            end if
         
         end if
      
      
      ! End new
      
      end do
      
      write(funit) cache(i)%dummy_entry
      write(funit) cache(i)%num_dmat
      write(funit) cache(i)%contrib_type
      write(funit) cache(i)%n_rule
      write(funit) cache(i)%contrib_size
      
      if (cache(i)%p_tuples(1)%npert > 0) then
      
         write(funit) size(cache(i)%nblks_tuple)
         write(funit) cache(i)%nblks_tuple
      
         write(funit) size(cache(i)%blk_sizes, 1)
         write(funit) size(cache(i)%blk_sizes, 2)
         do j = 1, size(cache(i)%blk_sizes, 1)
            write(funit) cache(i)%blk_sizes(j, :)
         end do
      
         write(funit) size(cache(i)%indices, 1)
         write(funit) size(cache(i)%indices, 2)
         
         do j = 1, size(cache(i)%indices, 1)
            write(funit) cache(i)%indices(j,:)
         end do
      
         write(funit) size(cache(i)%blks_tuple_info, 1)
         write(funit) size(cache(i)%blks_tuple_info, 2)
         write(funit) size(cache(i)%blks_tuple_info, 3)
      
         do j = 1, size(cache(i)%blks_tuple_info, 1)
            do k = 1, size(cache(i)%blks_tuple_info, 2)
               write(funit) cache(i)%blks_tuple_info(j, k, :)
            end do
         end do
      
         write(funit) size(cache(i)%blks_tuple_triang_size)
         write(funit) cache(i)%blks_tuple_triang_size
      
      end if
      
      if (allocated(cache(i)%data_mat)) then
      
         mat_scal_none = 0
         write(funit) mat_scal_none
         write(funit) size(cache(i)%data_mat)
         
         fid_fmt = '(I10.10)'
         
         do j = 1, size(cache(i)%data_mat)
            
            write(str_fid, fid_fmt) j + mat_acc
         
            k = QcMatWrite_f(cache(i)%data_mat(j), trim(adjustl(fname)) // '_MAT_' // &
                 trim(str_fid) // '.DAT', BINARY_VIEW)
         
         end do
         
         mat_acc = mat_acc + size(cache(i)%data_mat)
         
      
      elseif (allocated(cache(i)%data_scal)) then
      
         mat_scal_none = 1
         write(funit) mat_scal_none
         write(funit) size(cache(i)%data_scal)
         write(funit) cache(i)%data_scal
         
      
      else
      
         mat_scal_none = 2 
         write(funit) mat_scal_none
      
      end if
   
   end do
 
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
   

 end subroutine

 ! FIXME: Definitely ll-to-array chgs needed for this routine 
 ! Don't know what to do about mem mgr update yet
 
 ! Wrapper: Get outer contribution cache instance from file fname
 ! This is the routine to call if one wants to retrieve an entire
 ! outer cache array from a file and that outer cache is "standalone", i.e.
 ! not part of an inner cache instance
 subroutine contrib_cache_outer_retrieve(cache, fname, funit_in, mat_acc_in)
 
   implicit none
   
   logical :: cache_ext
   integer :: funit, i, j, k, mat_acc
   integer, optional :: funit_in, mat_acc_in
   integer :: num_entries
   character(*) :: fname
   type(contrib_cache_outer), dimension(:), allocatable :: cache
   
   mat_acc = 0
   if (present(mat_acc_in)) then
      mat_acc = mat_acc_in
   end if
   

   if (present(funit_in)) then
      funit = funit_in
   end if
   
  
   inquire(file=fname // '.DAT', exist=cache_ext)
      
   if (.NOT.(cache_ext)) then
   
      write(*,*) 'ERROR: The expected cache file ', trim(adjustl(fname)), ' does not exist'
      stop
   
   end if
   
   if present(funit_in) then
   
      open(unit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='old', action='read')
           
   else

      open(newunit=funit, file=trim(adjustl(fname)) // '.DAT', &
           form='unformatted', status='old', action='read')
   
   end if
  
   read(funit) num_entries
   
   if (allocated(cache)) then
   
      deallocate(cache)
      
   end if
   
   allocate(cache(num_entries))
   
   do i = 1, num_entries

      call contrib_cache_read_one_outer(cache(i), fname, funit, mat_acc_in=mat_acc)
   
   end do
   
   close(funit)
   
   if (present(mat_acc_in)) then
      mat_acc_in = mat_acc
   end if
 
 end subroutine

 ! FIXME: Definitely ll-to-array chgs needed for this routine 
 
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
 
 ! FIXME: Definitely ll-to-array chgs needed for this routine, maybe RM altogether
 
 ! Locate which inner and outer cache element matches the perturbation tuples
 function contrib_cache_locate(len_cache, cache, num_p_tuples, p_tuples, n_rule)

   implicit none

   integer, dimension(2) :: contrib_cache_locate
   integer :: num_p_tuples, i, j, k
   logical :: found
   integer, optional :: n_rule
   integer :: len_cache
   type(contrib_cache), dimension(len_cache) :: cache
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert

   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(len_cache, cache, p_tuples(1))) then
     
      do i = 1, len_cache
      
         if (p_tuple_compare(cache(i)%p_inner, p_tuples(1))) then
         
            contrib_cache_locate(1) = i
            exit
         
         end if
      
      end do
      
      found = .false.
      
      do i = 1, size(cache(contrib_cache_locate(1))%contribs_outer)
      
         if (cache(contrib_cache_locate(1))%contribs_outer(i)%dummy_entry) then
        
            cycle
        
         end if
         
         if (num_p_tuples > 1) then
         
            found = p_tuples_compare(num_p_tuples - 1, &
            cache(contrib_cache_locate(1))%contribs_outer(i)%p_tuples, &
            p_tuples(2:))
            
            if (present(n_rule)) then
      
               found = found .AND. (n_rule == &
               cache(contrib_cache_locate(1))%contribs_outer(i)%n_rule)
      
            end if
                                  
         else
         
            found = p_tuples_compare(num_p_tuples - 1, &
            cache(contrib_cache_locate(1))%contribs_outer(i)%p_tuples, &
            (/get_emptypert()/))

            if (present(n_rule)) then
      
               found = found .AND. (n_rule == &
               cache(contrib_cache_locate(1))%contribs_outer(i)%n_rule)
      
            end if            
            
         end if
         
         if (found) then
         
            contrib_cache_locate(2) = i
            exit
         
         end if
      
      end do
     
      if (.NOT.(found)) then
      
         write(*,*) 'ERROR: Did not find expected outer cache element'
         
      end if
      
   else
   
      write(*,*) 'ERROR: Did not find expected cache element'
      
   end if
   

 end function

  ! Allocate and set up new contribution cache element
  subroutine contrib_cache_allocate(current_element)

   implicit none

   type(contrib_cache), allocatable, dimension(:) :: current_element

   allocate(current_element(1))

   current_element%p_inner%npert = 0
   
   allocate(current_element(1)%p_inner%pdim(1))
   allocate(current_element(1)%p_inner%plab(1))
   allocate(current_element(1)%p_inner%pid(1))
   allocate(current_element(1)%p_inner%freq(1))
   
   current_element(1)%p_inner%pdim = (/0/)
   current_element(1)%p_inner%plab = (/'NUTN'/)
   current_element(1)%p_inner%pid = (/0/)
   current_element(1)%p_inner%freq = (/0.0/)
   
  call contrib_cache_outer_allocate(current_element(1)%contribs_outer)
   
 end subroutine

 ! Set up new outer contribution cache
 ! This creates a one-entry cache with a dummy entry 
 ! The cache can then later be extended
 subroutine contrib_cache_outer_allocate(current_element)

    implicit none

    type(contrib_cache_outer), allocatable, dimension(:) :: current_element
    
    allocate(current_element(1))
    
    current_element(1)%dummy_entry = .TRUE.
    current_element(1)%num_dmat = 0
        
    allocate(current_element(1)%p_tuples(1))

    current_element(1)%p_tuples(1)%npert = 0

    allocate(current_element(1)%p_tuples(1)%pdim(1))
    allocate(current_element(1)%p_tuples(1)%plab(1))
    allocate(current_element(1)%p_tuples(1)%pid(1))
    allocate(current_element(1)%p_tuples(1)%freq(1))
        
    current_element(1)%p_tuples(1)%pdim = (/0/)
    current_element(1)%p_tuples(1)%plab = (/'NUTN'/)
    current_element(1)%p_tuples(1)%pid = (/0/)
    current_element(1)%p_tuples(1)%freq = (/0.0/)

  end subroutine

 ! Initialize contribution cache (and associated outer contribution cache)
 ! as specified by perturbation tuples 'p_tuples'
 subroutine contrib_cache_initialize(new_element, num_p_tuples, p_tuples, n_rule)

   implicit none

   integer :: num_p_tuples
   integer, optional :: n_rule
   type(contrib_cache) :: new_element
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple), dimension(1) :: empty_tuple

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
      
         call contrib_cache_outer_add_element(size(new_element%contribs_outer), &
               new_element%contribs_outer, &
              .NOT.( p_tuples(2)%npert > 0), num_p_tuples - 1, &
               p_tuples(2:num_p_tuples), n_rule=n_rule)   
         
      else
         
         call contrib_cache_outer_add_element(size(new_element%contribs_outer), &
               new_element%contribs_outer, .FALSE., num_p_tuples - 1, &
                                              p_tuples(2:num_p_tuples))   
         
      end if
   
   else
   
      call contrib_cache_outer_allocate(new_element%contribs_outer)

      ! MaR: Passing empty perturbation tuple as dummy argument since there should
      ! be no outer perturbation tuples in this case; revisit if problems.
       
        call empty_p_tuple(empty_tuple(1))

      if (present(n_rule)) then
         
         call contrib_cache_outer_add_element(size(new_element%contribs_outer), &
              new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                             empty_tuple, n_rule=n_rule)
         
      else
      
         call contrib_cache_outer_add_element(size(new_element%contribs_outer), &
              new_element%contribs_outer, .TRUE., num_p_tuples - 1, &
                                              empty_tuple)
         
      end if
      
   end if
      
 end subroutine
 
 ! Initialize outer contribution cache element as specified by perturbation tuples 'outer_p_tuples'
 ! 'unperturbed' flag signifies empty perturbation tuple (typically used for 'all inner' case)
 subroutine contrib_cache_outer_initialize(new_element, unperturbed, num_dmat, outer_p_tuples)

   implicit none

   logical :: unperturbed
   integer :: num_dmat, i, total_npert
   integer, dimension(:), allocatable :: blki_a, blki_b
   integer, dimension(:,:), allocatable :: blksi_a
   type(contrib_cache_outer) :: new_element
   type(p_tuple), dimension(num_dmat) :: outer_p_tuples

   new_element%num_dmat = num_dmat
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

           allocate(blki_a(new_element%nblks_tuple(i)))
           allocate(blki_b(new_element%nblks_tuple(i)))
           allocate(blksi_a(new_element%nblks_tuple(i), &
                    size(new_element%blks_tuple_info, 3) ))

           new_element%blks_tuple_info(i, :, :) = get_blk_info(new_element%nblks_tuple(i), outer_p_tuples(i))


           blki_a = new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),2)
           blki_b = new_element%blks_tuple_info(i,1:new_element%nblks_tuple(i),3)
           blksi_a = new_element%blks_tuple_info(i, 1:new_element%nblks_tuple(i), :)

           new_element%blk_sizes(i, 1:new_element%nblks_tuple(i)) = &
           get_triangular_sizes(new_element%nblks_tuple(i), &
           blki_a, blki_b)

           new_element%blks_tuple_triang_size(i) = get_triangulated_size(new_element%nblks_tuple(i), &
           blksi_a)
                    
           deallocate(blki_a)
           deallocate(blki_b)
           deallocate(blksi_a)
               
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
   
      allocate(new_element%p_tuples(1)%pdim(0))
      allocate(new_element%p_tuples(1)%plab(0))
      allocate(new_element%p_tuples(1)%pid(0))
      allocate(new_element%p_tuples(1)%freq(0))
      
      ! MaR: Creating indices for unperturbed entry: Return here if problems arise
      
      allocate(new_element%indices(1,1))
      new_element%indices(1,1) = 1 
   
   end if
      
 end subroutine
   
 ! Add new element to outer contribution cache
 subroutine contrib_cache_outer_add_element(len_cache, cache, unperturbed, num_dmat, &
            outer_p_tuples, data_size, data_mat, data_scal, n_rule)

   implicit none

   logical :: unperturbed, already, found_element
   integer :: num_dmat, i, j, passedlast
   integer, optional :: data_size, n_rule
   integer :: len_cache, found_loc
   type(contrib_cache_outer), allocatable, dimension(:) :: cache
   type(contrib_cache_outer), allocatable, dimension(:) :: cache_store
   type(contrib_cache_outer) :: new_element
  

   type(p_tuple), dimension(num_dmat) :: outer_p_tuples, p_tuples_st_order, pstcmp
   
   type(Qcmat), optional, dimension(*) :: data_mat
   complex(8), optional, dimension(*) :: data_scal


   if (present(n_rule)) then
   
      already = contrib_cache_already_outer(len_cache, cache, num_dmat, outer_p_tuples, n_rule=n_rule)
   
   else
   
      already = contrib_cache_already_outer(len_cache, cache, num_dmat, outer_p_tuples)
   
   end if
   
   ! If an entry for these pert. tuples already exists in cache, then locate it
   ! and update the data there
   if (already) then
   
      found_element = .FALSE.
      p_tuples_st_order = p_tuples_standardorder(num_dmat, outer_p_tuples)
   
      do i = 1, len_cache

      
         if (cache(i)%num_dmat == num_dmat .AND. .NOT.(cache(i)%dummy_entry)) then

         pstcmp = p_tuples_standardorder(cache(i)%num_dmat, &
            cache(i)%p_tuples)


            found_element = p_tuples_compare(cache(i)%num_dmat, &
            pstcmp, p_tuples_st_order)
            
            if (found_element) then
               found_loc = i
               exit
            end if
            
         end if
            
      end do
   
      ! Update matrix data
      if(present(data_mat)) then
      
         do i = 1, data_size
         
            call QcMatAEqB(cache(found_loc)%data_mat(i), data_mat(i))
           
         end do
   
      end if
   
      ! Update scalar data: Currently not used
      if (present(data_scal)) then
   
      end if
   
   ! Otherwise, if an entry doesn't exist, then make a new entry in the cache
   else
   
      ! Create the new entry
      call contrib_cache_outer_initialize(new_element, unperturbed, num_dmat, outer_p_tuples)
         
      ! Fill it with data if provided  
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
   
      ! Scalar data not presently used in this part of the code
      if (present(data_scal)) then
   
      end if
      
      if (present(n_rule)) then
      
         new_element%n_rule = n_rule
      
      end if

      allocate(cache_store(len_cache))
      
      ! FIXME: May be necessary to copy element by element
      ! But I think that the default assignment method will work here
      cache_store = cache
      
      deallocate(cache)
      
      allocate(cache(len_cache + 1))
      
      ! FIXME: May be necessary to copy element by element
      ! But I think that the default assignment method will work here
      cache(1:len_cache) = cache_store
      cache(len_cache + 1) = new_element
      
      deallocate(cache_store)
   
   end if
      
 end subroutine
 
 ! Add contribution cache element to existing linked list
 subroutine contrib_cache_add_element(len_cache, cache, num_p_tuples, p_tuples, n_rule)

   implicit none

   logical :: empty_outer
   integer :: num_p_tuples, i, len_cache
   integer :: found_loc
   integer, optional :: n_rule
   type(contrib_cache), allocatable, dimension(:) :: cache
   type(contrib_cache), allocatable, dimension(:) :: cache_store
   type(contrib_cache), allocatable, dimension(:) :: new_element
   
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) :: emptypert

   ! If cache element for inner perturbations already exists, just add outer
   if (contrib_cache_already_inner(len_cache, cache, p_tuples(1))) then
   
      found_loc = 1
   
      ! Locate cache element for this inner
      do while (p_tuple_compare(cache(found_loc)%p_inner, p_tuples(1)) .EQV. .FALSE.)

         found_loc = found_loc + 1
          
      end do
      
      ! NOTE: AND CONDITION MAY ADD WRONG KIND OF CACHE ELEMENT, REVISIT IF ERROR
      
      empty_outer = .TRUE.
      
      if (num_p_tuples > 1) then
     
         if (p_tuples(2)%npert > 0) then
         
            empty_outer = .FALSE.
            cache(found_loc)%num_outer = cache(found_loc)%num_outer + 1
     
            if (present(n_rule)) then
               
               call contrib_cache_outer_add_element(size(cache(found_loc)%contribs_outer), &
                    cache(found_loc)%contribs_outer, .FALSE., &
                    num_p_tuples - 1, p_tuples(2:num_p_tuples), n_rule=n_rule)
         
            else
         
               call contrib_cache_outer_add_element(size(cache(found_loc)%contribs_outer), &
                    cache(found_loc)%contribs_outer, .FALSE., &
                    num_p_tuples - 1, p_tuples(2:num_p_tuples))
         
            end if
         
         end if
         
      end if
            
      if (empty_outer) then
     
         cache(found_loc)%num_outer = cache(found_loc)%num_outer + 1

         call empty_p_tuple(emptypert)
         
         ! MaR: WARNING: CHANGED num_p_tuples ARGUMENT TO 0 BELOW; MAY PRODUCE ERRORS ELSEWHERE
         ! UPDATE: Looks like it didn't make errors, return here if there are problems anyway

         if (present(n_rule)) then
         
            call contrib_cache_outer_add_element(size(cache(found_loc)%contribs_outer), &
                    cache(found_loc)%contribs_outer, .TRUE., 1, (/emptypert/), n_rule=n_rule)
         
         else
            
            call contrib_cache_outer_add_element(size(cache(found_loc)%contribs_outer), &
                    cache(found_loc)%contribs_outer, .TRUE., 0, (/emptypert/))
         
         end if
       
      end if


   ! Otherwise, add both inner and outer    
   else
   

      allocate(new_element(1))
      
      if (present(n_rule)) then
         
         call contrib_cache_initialize(new_element(1), num_p_tuples, p_tuples, n_rule)   
         
      else
         
         call contrib_cache_initialize(new_element(1), num_p_tuples, p_tuples)   
         
      end if
         
      new_element%num_outer = 1
      
      allocate(cache_store(len_cache))
      
      ! FIXME: May be necessary to copy element by element
      ! But I think that the default assignment method will work here
      cache_store = cache
      
      deallocate(cache)
      
      allocate(cache(len_cache + 1))
      
      ! FIXME: May be necessary to copy element by element
      ! But I think that the default assignment method will work here
      cache(1:len_cache) = cache_store
      cache(len_cache + 1) = new_element(1)
      
      deallocate(cache_store)
      
      len_cache = len_cache + 1
      
      
   end if

 end subroutine
 
 ! Check if element (inner + outer combination) as specified by 'p_tuples' (and possible 'n_rule')
 ! already exists in cache
 function contrib_cache_already(len_cache, cache, num_p_tuples, p_tuples, n_rule)

   implicit none

   logical :: contrib_cache_already
   integer :: num_p_tuples, i, len_cache
   integer, optional :: n_rule
   
   type(contrib_cache), dimension(len_cache) :: cache
   
   type(p_tuple), dimension(num_p_tuples) :: p_tuples
   type(p_tuple) ::  p_tuple_ord, p_tmp_ord
   type(p_tuple), dimension(1) :: emptypert


   p_tuple_ord = p_tuple_standardorder(p_tuples(1))
   contrib_cache_already = .FALSE.

   call empty_p_tuple(emptypert(1))
   
   do i = 1, len_cache
   
      p_tmp_ord = p_tuple_standardorder(cache(i)%p_inner)
      contrib_cache_already = p_tuple_compare(p_tmp_ord, p_tuple_ord)
      
      if (contrib_cache_already) then
         
         if (num_p_tuples > 1) then
         
            if (present(n_rule)) then
            
               contrib_cache_already = contrib_cache_already_outer(size(cache(i)%contribs_outer), &
                                       cache(i)%contribs_outer, &
                                       num_p_tuples - 1, p_tuples(2:num_p_tuples), n_rule=n_rule)
            
            else

               contrib_cache_already = contrib_cache_already_outer(size(cache(i)%contribs_outer), &
                                       cache(i)%contribs_outer, &
                                       num_p_tuples - 1, p_tuples(2:num_p_tuples))
            
            end if
            
            if (contrib_cache_already) then
            
               exit
               
            end if

            
         else
         
            if (present(n_rule)) then
            
               contrib_cache_already = contrib_cache_already_outer(size(cache(i)%contribs_outer), &
                                       cache(i)%contribs_outer, 0, emptypert, n_rule=n_rule)
            
            else
            
               contrib_cache_already = contrib_cache_already_outer(size(cache(i)%contribs_outer), &
                                       cache(i)%contribs_outer, 0, emptypert)
            
            end if    
            
            if (contrib_cache_already) then
            
               exit
               
            end if
                      
         end if
      
      end if
      
   end do
   
 end function
 
 ! Check if element of outer contribution cache as specified by 'p_tuples_outer'
 ! (and possible 'n_rule') exists in linked list
 function contrib_cache_already_outer(len_cache, cache, num_dmat, p_tuples_outer, n_rule)

   implicit none

   logical :: contrib_cache_already_outer
   integer :: num_dmat, i
   integer :: len_cache
   integer, optional :: n_rule
   type(contrib_cache_outer), dimension(len_cache) :: cache
   
   type(p_tuple), dimension(num_dmat) :: p_tuples_outer, p_tuples_ord
   
   p_tuples_ord = p_tuples_standardorder(num_dmat, p_tuples_outer)

   contrib_cache_already_outer = .FALSE.

   do i = 1, len_cache
   
      if (cache(i)%num_dmat == num_dmat .AND. .NOT.(cache(i)%dummy_entry)) then
      
         contrib_cache_already_outer = p_tuples_compare(num_dmat, &
                                       cache(i)%p_tuples, p_tuples_ord)
                                    
         if (present(n_rule)) then
      
            contrib_cache_already_outer = contrib_cache_already_outer .AND. (n_rule == cache(i)%n_rule)

         end if
         
                                 
      end if
      
      if (contrib_cache_already_outer) then
      
         exit
         
      end if
   
   end do

end function
 
 ! Check if contribution cache (inner only) element (as specified by 'p_inner') 
 ! already exists in linked list 
 function contrib_cache_already_inner(len_cache, cache, p_inner)

   implicit none

   logical :: contrib_cache_already_inner
   integer :: passedlast, i, len_cache
   type(contrib_cache), dimension(len_cache) :: cache
   
   type(p_tuple) :: p_inner, p_tuple_ord
   
   p_tuple_ord = p_tuple_standardorder(p_inner)
   contrib_cache_already_inner = .FALSE.

   do i = 1, len_cache
   
      contrib_cache_already_inner = p_tuple_compare(cache(i)%p_inner, p_tuple_ord)
      
      if (contrib_cache_already_inner) then
      
         exit
         
      end if
   
   end do
   
 end function
 
 ! Get data from outer contribution cache element    
 subroutine contrib_cache_getdata_outer(len_cache, cache, num_p_tuples, p_tuples, &
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
   integer :: len_cache, loc_found
   integer, allocatable, dimension(:) :: pids_in_cache, pids_current_contrib, & 
                                         p_tuples_dimensions, &
                                         p_tuples_dimensions_cacheorder, &
                                         pids_merged_pert, translated_index
   integer, allocatable, dimension(:) :: blk_sizes_merged
   integer, dimension(num_p_tuples) :: nfields, nblks_tuple, blks_tuple_triang_size
   integer, allocatable, dimension(:,:) :: indices, blk_sizes
   integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info
   integer, allocatable, dimension(:,:,:,:) :: blk_arg_iii
   integer, allocatable, dimension(:,:) :: blk_info, blk_arg_ii
   integer, allocatable, dimension(:) :: blk_sizes_sing, blk_arg_a, blk_arg_b
   integer, allocatable, dimension(:) :: blk_arg_c, blk_arg_d
   type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_ord
   type(p_tuple), allocatable, dimension(:) :: p_tuples_srch, p_tuples_srch_ord
   type(p_tuple) :: merged_p_tuple
   integer, dimension(ind_len), optional ::  ind_unsorted
   type(contrib_cache_outer), dimension(len_cache) :: cache
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
  
   found = .FALSE.
      
   do i = 1, len_cache
   
      if (cache(i)%dummy_entry) then
      
         cycle
         
      end if
   
      if ((size(p_tuples_srch_ord) == 0) .AND. &
          (size(cache(i)%p_tuples) == 1) ) then
             
         if ((cache(i)%p_tuples(1)%npert == 0)) then
            
            found = .TRUE.
            
         end if
            
      end if
      
      if (size(cache(i)%p_tuples) == size(p_tuples_srch_ord)) then
         found = p_tuples_compare(num_p_tuples - inner_rm, cache(i)%p_tuples, &
                                  p_tuples_srch_ord)
      end if   
                               
                                  
      if (present(n_rule)) then
      
         found = found .AND. (n_rule == cache(i)%n_rule)
      
      end if

      if (found) then
      
         loc_found = i
      
      exit
      
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

         allocate(blk_arg_a(merged_nblks))
         allocate(blk_arg_b(merged_nblks))

         blk_arg_a=merged_blk_info(1, 1:merged_nblks, 2)
         blk_arg_b=merged_blk_info(1, 1:merged_nblks, 3)

         blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
         blk_arg_a, blk_arg_b)

         deallocate(blk_arg_a)
         deallocate(blk_arg_b)

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

            allocate(blk_arg_a(nblks_tuple(i)))
            allocate(blk_arg_b(nblks_tuple(i)))
            allocate(blk_arg_ii(nblks_tuple(i), size(blks_tuple_info, 3)))

            blk_arg_a = blks_tuple_info(i, 1:nblks_tuple(i), 2)
            blk_arg_b = blks_tuple_info(i, 1:nblks_tuple(i), 3)
            blk_arg_ii = blks_tuple_info(i, 1:nblks_tuple(i), :)


            blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                        blk_arg_ii)
            blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
            blk_arg_a, blk_arg_b)

            deallocate(blk_arg_a)
            deallocate(blk_arg_b)
            deallocate(blk_arg_ii)

         end do

         ! Loop over indices, get appropriate offsets and put data in return array
         do i = 1, size(indices, 1)

            allocate(blk_arg_a(1))
            allocate(blk_arg_b(1))
            allocate(blk_arg_c(1))
            allocate(blk_arg_d(size(indices, 2)))
         

            blk_arg_a(1) = merged_nblks
            blk_arg_b(1) = total_num_perturbations
            blk_arg_c(1) = merged_triang_size
            blk_arg_d = indices(i, :)

            res_offset = get_triang_blks_tuple_offset(1, merged_nblks, &
            blk_arg_a, & 
            blk_arg_b, merged_blk_info, blk_sizes_merged, &
            blk_arg_c, &
            blk_arg_d )

            deallocate(blk_arg_d)
            deallocate(blk_arg_c)
            deallocate(blk_arg_b)
            deallocate(blk_arg_a)

            do j = 1, total_num_perturbations

            ! FIXME: NEXT LINE WENT OUT OF BOUNDS IN EARLIER TESTS, SHOULD BE OK NOW BUT
            ! THIS COMMENT KEPT TO MARK POTENTIAL SOURCE OF LATER PROBLEMS
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

                     allocate(blk_arg_a(size(blk_sizes(1,:))))
                     
                     blk_arg_a = blk_sizes(1,:)
                     
                     allocate(blk_arg_ii(size(blks_tuple_info,2),size(blks_tuple_info,3)))

                     blk_arg_ii = blks_tuple_info(1,:,:)

               
                     cache_offset = get_triang_blks_tuple_offset(1, &
                     total_num_perturbations, nblks_tuple(1), & 
                     nfields(1), blk_arg_ii, blk_arg_a, &
                     blks_tuple_triang_size(1), translated_index)

                     deallocate(blk_arg_ii)

                     deallocate(blk_arg_a)
               
                  end if
            
               else
            
                  cache_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                  total_num_perturbations, nblks_tuple, & 
                  nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, translated_index)
                        
               end if
     
            end if
            
            if (present(mat)) then
            
               call QcMatRAXPY(1.0d0, cache(loc_found)%data_mat(cache_offset + &
                               cache_hard_offset), mat(res_offset))
               
            else if (present(scal)) then
            
               scal(res_offset) = &
               scal(res_offset) + &
               cache(loc_found)%data_scal(cache_offset + cache_hard_offset)              
            
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
      
!         write(*,*) 'mat sing cache offset', offset + cache_hard_offset
      
         call QcMatAEqB(mat_sing, cache(loc_found)%data_mat(offset + cache_hard_offset))
         
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
 subroutine contrib_cache_getdata(len_cache, cache, num_p_tuples, p_tuples, contrib_size, &
                                  ind_len, ind_unsorted, hard_offset, mat, mat_sing, scal, n_rule)

   implicit none

   logical :: found
   integer :: i, j, k, first, last, passedlast, num_p_tuples, &
              contrib_size, total_num_perturbations, pr_offset, cache_offset, &
              merged_triang_size, merged_nblks, res_offset, ind_len, &
              cache_hard_offset
   integer :: outer_size 
   integer :: loc_found
   integer :: len_cache
   integer, optional :: hard_offset, n_rule
   integer, optional, dimension(ind_len) :: ind_unsorted
   type(contrib_cache), dimension(len_cache) :: cache

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
   found = .FALSE.
   p_tuples_ord = p_tuples_standardorder(num_p_tuples, p_tuples)

   do i = 1, len_cache

! FIXME CMNTD OUT, LOOKS LIKE NOT NEEDED BUT KEPT FOR REVISIT IF LATER PROBLEMS
!       if (cache(i)%dummy_entry) then
!       
!          cycle
!          
!       end if
   
      found = p_tuple_compare(cache(i)%p_inner, p_tuples_ord(1))
      
      if (found) then
      
         loc_found = i
         exit
         
      end if   
      
   end do
   

   ! If found, move to search (and retrieve) from outer cache
   ! Several cases depending on which type of data is requested,
   ! and whether or not (k,n) rule choice (specified as n in 'n_rule')
   ! is a discriminating factor in search (the latter typical for two-factor terms)
   if (found) then
   
      outer_size = size(cache(loc_found)%contribs_outer)
   
   
      if (present(mat)) then

         if (present(n_rule)) then
             
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, hard_offset, &
                 mat=mat, n_rule=n_rule)
       
         else
         
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, hard_offset, mat=mat)
         
         end if
      
      else if (present(mat_sing)) then
      
         if (present(n_rule)) then
       
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 mat_sing=mat_sing, n_rule=n_rule)
       
         else
         
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 mat_sing=mat_sing)
                 
         end if
              
      else if (present(scal)) then
      
         if (present(n_rule)) then
         
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 scal=scal, n_rule=n_rule)       
       
         else
         
            call contrib_cache_getdata_outer(outer_size, cache(loc_found)%contribs_outer, &
                 num_p_tuples, &
                 p_tuples, .TRUE., contrib_size, ind_len, ind_unsorted, cache_hard_offset, &
                 scal=scal)
         
         end if
         
      end if
      
   else

      write(*,*) 'Failed to retrieve data in contrib_cache_getdata: Inner element not found'

   end if


 end subroutine
 
 subroutine contrib_cache_outer_dealloc_mat(len_cache, cache)
 
   implicit none
 
   integer :: len_cache, i, j
   type(contrib_cache_outer), dimension(len_cache) :: cache
   
   do i = 1, len_cache
   
      if (allocated(cache(i)%data_mat)) then
      
         do j = 1, size(cache(i)%data_mat)
         
            call QcMatDst(cache(i)%data_mat(j))
         
         end do
      
      end if
   
   end do
 
 
 end subroutine

end module
