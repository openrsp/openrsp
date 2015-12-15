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
  public p1_cloneto_p2
  public p_tuple_deallocate
  public get_emptypert
  public empty_p_tuple
  public p_tuple_p1_lt_p2
  public p_tuple_standardorder
  public p_tuples_standardorder
  public p_tuples_compare
  public plab_compare
  public pfreq_compare
  public pid_compare
  public p_tuple_compare
  public make_p_tuple_subset

  ! KEEP 2014
  
  
  
  
  
  ! END KEEP 2014
  
  ! NEW 2014
  
  public p_tuple_to_external_orders
  public p_tuple_to_external_tuple
  
  
  
  public p_tuple_dealloc
  public p1_lt_p2
!   public p_tuples_compare_2014
!   public p_tuple_compare_2014
  public get_p_subset
!   public p1_cloneto_p2
  public p1_switch_p2
  public p1_merge_p2
!   public p_tuple_ordered
!   public p_tuples_ordered
  public make_p_tuple_subsets
  
  
!   ! Perturbation datatype
!   type pert
! 
!      integer :: pdim ! Number of components
!      integer :: pfcomp ! First component is component #pfcomp
!      integer :: plab ! Perturbation label
!      integer :: pid ! Perturbation ID
!      complex(8) :: freq ! Frequency of perturbation
! 
!   end type
  
!   ! Perturbation tuple datatype
!   type p_tuple_2014
! 
! !      integer :: n_perturbations ! Number of perturbations
!      integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
!      integer, allocatable, dimension(:) :: plab ! Perturbation labels
!      integer, allocatable, dimension(:) :: pfcomp ! First component is component #pfcomp
!      integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
!      complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
!      ! Add other perturbation identification info as needed
! 
!   end type
  
  ! END NEW 2014
  
  
  ! Define perturbation tuple datatype

  type p_tuple

!      integer :: n_perturbations ! Number of perturbations
     integer :: npert
     integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
     character(4), allocatable, dimension(:) :: plab ! Perturbation labels
     integer, allocatable, dimension(:) :: pfcomp ! First component is component #pfcomp
     integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
     complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
     ! Add other perturbation identification info as needed

  end type

  contains

  ! NEW NOVEMBER 2014
  
  subroutine p_tuple_to_external_orders(pert, num_pert, perturbations, pert_orders)
  
    implicit none
    
    type(p_tuple) :: pert
    integer :: num_pert, i, mp, j
    integer, dimension(pert%npert) :: plab_int
    integer, dimension(:), allocatable :: perturbations, pert_orders
    
    
    do i = 1, pert%npert
    
       if (pert%plab(i) == 'GEO ') then
       
          plab_int(i) = 1
       
       else if (pert%plab(i) == 'EL  ') then

          plab_int(i) = 2
       
       else if (pert%plab(i) == 'ELGR') then
       
          plab_int(i) = 3
       
       else if (pert%plab(i) == 'MAG ') then
       
          plab_int(i) = 4
       
       else if (pert%plab(i) == 'MAG0 ') then
    
          plab_int(i) = 5
          
       end if
    
    end do
    
    mp = maxval(plab_int)
    
    j = 0
    
    do i = 1, mp
    
       if (count(plab_int == i) > 0) then
          
          j = j + 1
       
       end if
    
    end do
    
    num_pert = j
    allocate(perturbations(num_pert))
    allocate(pert_orders(num_pert))
    
    j = 0
    
    do i = 1, mp
    
       if (count(plab_int == i) > 0) then
          
!FIXME Gao: j should be increment of 1?
          j = j+1
          perturbations(j) = i
          pert_orders(j) = count(plab_int == i)
       
       end if
    
    end do
      
  end subroutine


  subroutine p_tuple_to_external_tuple(pert, num_pert, perturbations)
  
    implicit none
    
    type(p_tuple) :: pert
    integer :: num_pert, i
    integer, dimension(:), allocatable :: perturbations
    
    num_pert = pert%npert
    allocate(perturbations(num_pert))

    do i = 1, pert%npert
    
       if (pert%plab(i) == 'GEO ') then
       
          perturbations(i) = 1
       
       else if (pert%plab(i) == 'EL  ') then

          perturbations(i) = 2
       
       else if (pert%plab(i) == 'ELGR') then
       
          perturbations(i) = 3
       
       else if (pert%plab(i) == 'MAG ') then
       
          perturbations(i) = 4
       
       else if (pert%plab(i) == 'MAG0 ') then
    
          perturbations(i) = 5
          
       end if
    
    end do
    
  end subroutine
  

  ! NEW 2014
  
  
   ! Deallocate perturbations in perturbation tuple
    subroutine p_tuple_dealloc(p1)

      implicit none

      type(p_tuple) :: p1

      p1%npert = 0

      deallocate(p1%pdim) 
      deallocate(p1%plab)
      deallocate(p1%pid)
      if (allocated(p1%pfcomp)) then
         deallocate(p1%pfcomp)
      end if    
      deallocate(p1%freq)

     
    end subroutine

    ! Is perturbation tuple p1 "less than" tuple p2 by the canonical metric?
    function p1_lt_p2(p1, p2)

      implicit none

      type(p_tuple), intent(in) :: p1, p2
      type(p_tuple) :: p1_ord, p2_ord
      integer :: i
      logical :: p1_lt_p2

      p1_ord = p_tuple_standardorder(p1)
      p2_ord = p_tuple_standardorder(p2)
      
      p1_lt_p2 = .FALSE.

      if (p1%npert < p2%npert) then

         p1_lt_p2 = .TRUE.
         return

      elseif (p1%npert == p2%npert) then

         do i = 1, p1%npert

            if (p1%pdim(i) < p2%pdim(i)) then

               p1_lt_p2 = .TRUE.
               return

            elseif (p1%pdim(i) == p2%pdim(i)) then

               if (p1%pfcomp(i) < p2%pfcomp(i)) then

                  p1_lt_p2 = .TRUE.
                  return

               elseif (p1%pfcomp(i) == p2%pfcomp(i)) then

                  if (p1%plab(i)< p2%plab(i)) then

                     p1_lt_p2 = .TRUE.
                     return

                  elseif (p1%plab(i) == p2%plab(i)) then

                     if (real(p1%freq(i)) < real(p2%freq(i))) then

                        p1_lt_p2 = .TRUE.
                        return

                     end if

                  end if

               end if

            end if

         end do

      end if

      call p_tuple_dealloc(p1_ord)
      call p_tuple_dealloc(p2_ord)

    end function
    
    
!    ! Compare two tuples of perturbation tuples to see if they are equivalent
!     function p_tuples_compare_2014(num_p_tuples, p_tuples_1, p_tuples_2)
!     
!       implicit none
!       
!       logical :: p_tuples_compare_2014, none_unequal_yet
!       integer :: num_p_tuples
!       integer :: i
!       type(p_tuple), dimension(num_p_tuples) :: p_tuples_1, p_tuples_2
!       type(p_tuple), dimension(num_p_tuples) :: p_tuples_1_ord, p_tuples_2_ord
!       
!       call p_tuples_ordered(num_p_tuples, p_tuples_1, p_tuples_1_ord)
!       call p_tuples_ordered(num_p_tuples, p_tuples_2, p_tuples_2_ord)
!       
!       p_tuples_compare_2014 = .FALSE.
!       none_unequal_yet = .TRUE.
!       
!       do i = 1, num_p_tuples
!       
!          none_unequal_yet = none_unequal_yet .AND. &
!                             p_tuple_compare_2014(p_tuples_1_ord(i), p_tuples_2_ord(i))
!       
!       end do
!       
!       if (none_unequal_yet) p_tuples_compare_2014 = .TRUE.
!       
!     end function
  
!     ! Compare two perturbation tuples to see if they are equivalent
!     function p_tuple_compare_2014(p1, p2)
!     
!       implicit none
!       
!       logical :: p_tuple_compare_2014
!       integer :: i
!       type (p_tuple) :: p1, p2
!       
!       p_tuple_compare_2014 = .TRUE.
!       
!       if (p1%npert == p2%npert) then
!       
!          do i = 1, p1%npert
!       
!            if (.NOT. (p1%pdim(i) == p2%pdim(i))) p_tuple_compare_2014 = .FALSE.
!            if (.NOT. (p1%pfcomp(i) == p2%pfcomp(i))) p_tuple_compare_2014 = .FALSE.
!            if (.NOT. (p1%plab(i) == p2%plab(i))) p_tuple_compare_2014 = .FALSE.
!            if (.NOT. (p1%freq(i) == p2%freq(i))) p_tuple_compare_2014 = .FALSE.
!             
!          end do
!       
!       else
!       
!          p_tuple_compare_2014 = .FALSE.
!       
!       end if
!       
!     end function
  
    ! Get a general subset of perturbation tuple p1   
    subroutine get_p_subset(p1, sub_size, which_perts, p_subset)
    
      implicit none
      
      type(p_tuple), intent(in) :: p1
      type(p_tuple), intent(out) :: p_subset
      integer :: i, sub_size
      integer, dimension(sub_size) :: which_perts
      
!       p_subset%npert = sub_size
!       allocate(p_subset%perts(sub_size))
!       
!       p_subset%perts(1) = p1%perts(which_perts(1))
!       
!       do i = 2, sub_size
!       
!          p_subset%perts = (/p_subset%perts(:), p1%perts(which_perts(i))/)
!       
!       end do
      
    end subroutine
  
  ! Put all perturbations in p1 into p2 except perturbation 'except'
    subroutine p1_cloneto_p2_except(p1, p2, except)
    
      implicit none
      
      type(p_tuple), intent(in) :: p1
      type(p_tuple), intent(inout) :: p2
      integer :: except, startpos, i
     
      p2%npert = p1%npert - 1
!       allocate(p2%perts(p2%npert))
!       
!       startpos = 2
!       
!       if (.NOT.(except == 1)) then
!       
!          p2%perts(1) = p1%perts(1)
!          
!       else
!       
!          if (p1%npert == 1) then
!          
!             return
!             
!          else
!       
!             p2%perts(1) = p1%perts(2)
!             startpos = 3
!             
!          end if
!            
!       end if
!       
!       do i = startpos, p1%npert
!       
!          p2%perts = (/p2%perts(1:startpos - 1), p1%perts(i)/)
!          startpos = startpos + 1
!          
!       end do
      
    end subroutine
  
  
  
!    ! Clone perturbation tuple p1 into p2
!     subroutine p1_cloneto_p2(p1, p2)
! 
!       implicit none
! 
!       type(p_tuple), intent(in) :: p1
!       type(p_tuple), intent(out) :: p2
! 
!       p2%npert = p1%npert
! !       if (allocated(p2%perts)) deallocate(p2%perts)
! !       allocate(p2%perts(p2%npert))
! !       p2%perts = p1%perts
! 
!     end subroutine


    ! Switch perturbation tuples p1 and p2
    subroutine p1_switch_p2(p1, p2)

      implicit none

      type(p_tuple), intent(inout) :: p1
      type(p_tuple) :: p_tmp
      type(p_tuple), intent(out) :: p2

      call p1_cloneto_p2(p1, p_tmp)
      call p1_cloneto_p2(p2, p1)
      call p1_cloneto_p2(p_tmp, p2)
      call p_tuple_dealloc(p_tmp)

    end subroutine

    
    ! Merge two perturbation tuples p1 and p2 into p_merge
    subroutine p1_merge_p2(p1, p2, p_merge)
    
      implicit none
      
      type(p_tuple) :: p1, p2, p_merge
      
      p_merge = merge_p_tuple(p1, p2)
      
      
!       call p1_cloneto_p2(p1, p_merge)
!       p_merge%npert = p_merge%npert + p2%npert
! !       p_merge%perts = (/p_merge%perts, p2%perts/)
      
    end subroutine

!      ! Order perturbation tuple
!     subroutine p_tuple_ordered(p1, p_ord)
! 
!       implicit none
! 
!       integer :: i, j, min_pdim, min_pfcomp, min_plab, new_min
!       complex(8) :: min_freq
!       type(p_tuple), intent(in) :: p1
! !       type(pert) :: temp
!       type(p_tuple), intent(out) :: p_ord
! 
!       call p1_cloneto_p2(p1, p_ord)
! 
!          write(*,*) 'npert', p1%npert
!    write(*,*) 'pid', p1%pid
!    write(*,*) 'pfreq', p1%freq
!    write(*,*) 'pdim', p1%pdim
!       
!          write(*,*) 'p ord npert', p_ord%npert
!    write(*,*) 'pid', p_ord%pid
!    write(*,*) 'pfreq', p_ord%freq
!    write(*,*) 'pdim', p_ord%pdim
!       
!       
!       
!       do i = 1, p_ord%npert
! 
!          min_pdim = p_ord%pdim(i)
!          min_pfcomp = p_ord%pfcomp(i)
! !          min_plab = p_ord%plab(i)
!          min_freq = p_ord%freq(i)
!          new_min = i
! 
!          do j = i + 1, p_ord%npert
! 
!             if (p_ord%pdim(j) < min_pdim) then
! 
!                new_min = j
! 
!             elseif (p_ord%pdim(j) == min_pdim) then
! 
!                if (p_ord%pfcomp(j) < min_pfcomp) then
! 
!                   new_min = j
! 
!                elseif (p_ord%pfcomp(j) == min_pfcomp) then
! 
! !                   if (p_ord%plab(j) < min_plab) then
! ! 
! !                      new_min = j
! ! 
! !                   elseif (p_ord%plab(j) == min_plab) then
! ! 
! !                      if (real(p_ord%freq(j)) < real(min_freq)) then
! ! 
! !                         new_min = j
! ! 
! !                      end if
! ! 
! !                   end if
! 
!                end if
! 
!             end if
! 
!          end do
!          
! !          temp = p_ord%perts(new_min)
! !          p_ord%perts(new_min) = p_ord%perts(i)
! !          p_ord%perts(i) = temp
! 
!       end do
! 
!       ! Maybe some sorting remains: Find out during testing
! 
!     end subroutine
!     
!       ! Order tuples of perturbation tuples
!     ! Also order the tuples themselves
!     subroutine p_tuples_ordered(num_p_tuples, p_tuples, p_tuples_ord)
! 
!       implicit none
! 
!       integer, intent(in) :: num_p_tuples
!       integer :: i, j, new_min
!       type(p_tuple), dimension(num_p_tuples), intent(in) :: p_tuples
!       type(p_tuple), dimension(num_p_tuples), intent(out) :: p_tuples_ord
! 
!       do i = 1, num_p_tuples
!          
!          call p_tuple_ordered(p_tuples(i), p_tuples_ord(i))
! 
!       end do
! 
!       ! No reordering of first tuple since it contains "inner" perturbations
!       ! Only reordering of D contraction tuples: Therefore starting with i = 2
!       do i = 2, num_p_tuples
! 
!          new_min = i
! 
!          do j = i + 1, num_p_tuples 
! 
!             if (p1_lt_p2(p_tuples_ord(i), p_tuples_ord(j))) then
! 
!                new_min = j
!                call p_tuple_ordered(p_tuples_ord(new_min), p_tuples_ord(i))
! 
!             end if
! 
!          end do
! 
!          if (.NOT.(new_min == i)) then
! 
!             call p1_switch_p2(p_tuples_ord(i), p_tuples_ord(new_min))
! 
!          end if
! 
!       end do
! 
!     end subroutine
    
    
    
     ! Make all (p1%npert - 1) subsets of perturbation tuple p1 
    ! Assumes that p_sub is unallocated
    subroutine make_p_tuple_subsets(p1, p_sub)
    
      implicit none
      
      type(p_tuple) :: p1
      type(p_tuple), dimension(p1%npert) :: p_sub
      integer :: i, j, k, m
      
      do i = 1, p1%npert
      
         p_sub(i)%npert = p1%npert - 1
!          allocate(p_sub(i)%perts(p_sub(i)%npert))
         
         k = i
         m = 1
         
         do j = 1, p1%npert
            
            if (.NOT.(j == k)) then
            
!                p_sub(i)%perts(m) = p1%perts(j)
               
               m = m +1
            
            end if
            
         end do
      
      end do
      
    end subroutine
      
  ! END NEW 2014
  
  
  
  function p_tuple_extend(pert, ext)

    implicit none

    type(p_tuple) :: pert, ext, p_tuple_extend
    integer :: i

    allocate(p_tuple_extend%pdim(pert%npert + 1))
    allocate(p_tuple_extend%plab(pert%npert + 1))
    allocate(p_tuple_extend%pid(pert%npert + 1))
    allocate(p_tuple_extend%freq(pert%npert + 1))

    if (pert%npert == 0) then

       p_tuple_extend%npert = pert%npert + 1
       p_tuple_extend%pdim = (/ext%pdim(:)/)
       p_tuple_extend%plab = (/ext%plab(:)/)
       p_tuple_extend%pid = (/ext%pid(:)/)
       p_tuple_extend%freq = (/ext%freq(:)/)

    else

       p_tuple_extend%npert = pert%npert + 1
       p_tuple_extend%pdim = (/(pert%pdim(i), i = 1, pert%npert), ext%pdim(:)/)
       p_tuple_extend%plab = (/(pert%plab(i), i = 1, pert%npert), ext%plab(:)/)
       p_tuple_extend%pid = (/(pert%pid(i), i = 1, pert%npert), ext%pid(:)/)
       p_tuple_extend%freq = (/(pert%freq(i), i = 1, pert%npert), ext%freq(:)/)
 
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

    p_tuple_getone%npert = 1
    p_tuple_getone%pdim = (/pert%pdim(which)/)
    p_tuple_getone%plab = (/pert%plab(which)/)
    p_tuple_getone%pid = (/pert%pid(which)/)
    p_tuple_getone%freq = (/pert%freq(which)/)

  end function


  function p_tuple_remove_first(pert)

    implicit none

    type(p_tuple) :: pert, p_tuple_remove_first

    allocate(p_tuple_remove_first%pdim(pert%npert - 1))
    allocate(p_tuple_remove_first%plab(pert%npert - 1))
    allocate(p_tuple_remove_first%pid(pert%npert - 1))
    allocate(p_tuple_remove_first%freq(pert%npert - 1))

    if (pert%npert > 1) then

       p_tuple_remove_first%npert = pert%npert - 1
       p_tuple_remove_first%pdim = (/pert%pdim(2:pert%npert)/)
       p_tuple_remove_first%plab = (/pert%plab(2:pert%npert)/)
       p_tuple_remove_first%pid = (/pert%pid(2:pert%npert)/)
       p_tuple_remove_first%freq = (/pert%freq(2:pert%npert)/)

    else

       p_tuple_remove_first%npert = 0
       p_tuple_remove_first%pdim = (/0/)
       p_tuple_remove_first%plab = (/'NUTN'/)
       p_tuple_remove_first%pid = (/0/)
       p_tuple_remove_first%freq = (/0.0/)

    end if



  end function

  function merge_p_tuple(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2, merge_p_tuple

    allocate(merge_p_tuple%pdim(p1%npert + p2%npert))
    allocate(merge_p_tuple%plab(p1%npert + p2%npert))
    allocate(merge_p_tuple%pid(p1%npert + p2%npert))
    allocate(merge_p_tuple%freq(p1%npert + p2%npert))

    ! NOTE (MaR): ARE THESE CASE DISTINCTIONS UNNECESSARY? CONSIDER REWRITE.

    if (p1%npert > 0 .AND. p2%npert > 0) then

       merge_p_tuple%npert = p1%npert + p2%npert
       merge_p_tuple%pdim = (/p1%pdim(:), p2%pdim(:)/)
       merge_p_tuple%plab = (/p1%plab(:), p2%plab(:)/)
       merge_p_tuple%pid = (/p1%pid(:), p2%pid(:)/)
       merge_p_tuple%freq = (/p1%freq(:), p2%freq(:)/)

    elseif (p1%npert > 0 .AND. p2%npert == 0) then

       merge_p_tuple%npert = p1%npert
       merge_p_tuple%pdim = p1%pdim(:)
       merge_p_tuple%plab = p1%plab(:)
       merge_p_tuple%pid = p1%pid(:)
       merge_p_tuple%freq = p1%freq(:)

    elseif (p1%npert == 0 .AND. p2%npert > 0) then

       merge_p_tuple%npert = p2%npert
       merge_p_tuple%pdim = p2%pdim(:)
       merge_p_tuple%plab = p2%plab(:)
       merge_p_tuple%pid = p2%pid(:)
       merge_p_tuple%freq = p2%freq(:)

    elseif (p1%npert == 0 .AND. p2%npert == 0) then

       merge_p_tuple = get_emptypert()
       ! MaR: KEEP NEXT LINES FOR REVERSION IN CASE get_emptypert() DOESN'T WORK
       ! merge_p_tuple%npert = 0
       ! merge_p_tuple%pdim = (/0/)
       ! merge_p_tuple%plab = (/'NUTN'/)
       ! merge_p_tuple%pid = (/0/)
       ! merge_p_tuple%freq = (/0.0/)

    else

       write(*,*) 'Error in merge_p_tuple: Unrecognized size of p1 or p2 or both:', &
                   p1%npert, p2%npert

    end if

  end function

  subroutine p1_cloneto_p2(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2

!     write(*,*) 'presenting p1', p1%npert
!     write(*,*) 'dim', p1%pdim
!     write(*,*) 'lab', p1%plab
!     write(*,*) 'pid', p1%pid
!     write(*,*) 'freq', p1%freq
    
    p2%npert = p1%npert

    allocate(p2%pdim(p1%npert)) 
    allocate(p2%plab(p1%npert))
    allocate(p2%pid(p1%npert))
    allocate(p2%freq(p1%npert))

    p2%pdim = p1%pdim
    p2%plab = p1%plab
    p2%pid = p1%pid
    p2%freq = p1%freq

  end subroutine

  subroutine p_tuple_deallocate(p1)

    type(p_tuple) :: p1

    p1%npert = 0

    deallocate(p1%pdim) 
    deallocate(p1%plab)
    deallocate(p1%pid)
    if (allocated(p1%pfcomp)) then
       deallocate(p1%pfcomp)
    end if    
    deallocate(p1%freq)

    
  end subroutine


  function get_emptypert() result(emptypert)

    implicit none

      type(p_tuple) :: emptypert

      emptypert%npert = 0
      allocate(emptypert%pdim(0))    
      allocate(emptypert%plab(0))
      allocate(emptypert%pid(0))
      allocate(emptypert%freq(0))

  end function

  
  subroutine empty_p_tuple(emptypert)

    implicit none

      type(p_tuple) :: emptypert

      emptypert%npert = 0
      allocate(emptypert%pdim(0))    
      allocate(emptypert%plab(0))
      allocate(emptypert%pid(0))
      allocate(emptypert%freq(0))

  end subroutine


! Compare two perturbation tuples to each other 
! Assumes that input is already in standard order

  function p_tuple_p1_lt_p2(p1, p2)

    implicit none

    logical :: p_tuple_p1_lt_p2
    integer :: i, prn
    type(p_tuple) :: p1, p2

prn = 0

if (p1%npert == 2 .and. p2%npert == 2) then

if (p1%pid(1) == 6 .or. p2%pid(1) == 6) then

if (p1%pid(2) == 7 .or. p2%pid(2) == 7) then

if (p1%pid(1) == 1 .or. p2%pid(1) == 1) then

prn = 1

end if

end if

end if

end if

    ! NOTE (MaR): COULD THERE BE FALSE NEGATIVES HERE? 
    p_tuple_p1_lt_p2 = .FALSE.


    ! Compare number of perturbations
    ! REMEMBER: INCREASING ORDER OF DIFFERENTIATION
    if (p1%npert < p2%npert) then

       p_tuple_p1_lt_p2 = .TRUE.

    elseif (p1%npert == p2%npert) then

       do i = 1, p1%npert

          ! Compare dimensionality
          ! REMEMBER: DECREASING ORDER OF DIMENSIONALITY
          if (p1%pdim(i) > p2%pdim(i)) then

             p_tuple_p1_lt_p2 = .TRUE.
             exit

          elseif (p1%pdim(i) < p2%pdim(i)) then


             p_tuple_p1_lt_p2 = .FALSE.
             exit

          elseif (p1%pdim(i) == p2%pdim(i)) then

             if (llt(p1%plab(i), p2%plab(i)) .eqv. .TRUE.) then

                p_tuple_p1_lt_p2 = .TRUE.  
                exit


             elseif (lgt(p1%plab(i), p2%plab(i)) .eqv. .TRUE.) then


                p_tuple_p1_lt_p2 = .FALSE.  
                exit

             elseif (p1%plab(i) == p2%plab(i)) then


                ! NOTE (MaR): IS IT SUFFICIENTLY GENERAL TO COMPARE ONLY THE REAL PART OF
                ! THE FREQS.? WHICH CASES WILL INCLUDE COMPLEX FREQS. IN THE PERTURBATIONS?
                if (real(p1%freq(i)) < real(p2%freq(i))) then


                   p_tuple_p1_lt_p2 = .TRUE.  
                   exit


                elseif (real(p1%freq(i)) > real(p2%freq(i))) then

                   p_tuple_p1_lt_p2 = .FALSE.  
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


    call p1_cloneto_p2(pert, p_tuple_st)

    position_in_result = 1



    do i = position_in_result, p_tuple_st%npert

       current_minimum_pdim = p_tuple_st%pdim(i)
       current_minimum_plab = p_tuple_st%plab(i)
       current_minimum_freq = p_tuple_st%freq(i)
       new_minimum = i

       do j = i + 1, p_tuple_st%npert

! write(*,*) 'j', j

          if (p_tuple_st%pdim(j) >= current_minimum_pdim) then

             if (p_tuple_st%pdim(j) == current_minimum_pdim) then

                if (lge(p_tuple_st%plab(j), current_minimum_plab) .EQV. .TRUE.) then


! write(*,*) p_tuple_st%plab(j), 'lt', current_minimum_plab

                   if (p_tuple_st%plab(j) == current_minimum_plab) then

                      ! NOTE (MaR): IS IT SUFFICIENTLY GENERAL TO COMPARE
                      ! ONLY THE REAL PART OF THE FREQS.?
                      if (real(p_tuple_st%freq(j)) < real(current_minimum_freq)) then

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

    do while (current_last_position <= p_tuple_st%npert)

       if (current_last_position < p_tuple_st%npert) then

          do while ((p_tuple_st%plab(current_last_position) == &
                     p_tuple_st%plab(current_first_position)))

             current_last_position = current_last_position + 1
             if (current_last_position > p_tuple_st%npert) exit

          end do

          current_last_position = current_last_position - 1

       end if

       do i = current_first_position, current_last_position, 1

          current_minimum_freq = p_tuple_st%freq(i)
          new_minimum = i

          do j = i + 1, current_last_position

             if (real(p_tuple_st%freq(j)) < real(current_minimum_freq)) then

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
    integer :: i, j, k, new_minimum, max_order_curr, prn
    integer :: temporary_pdim, temporary_pid, current_first_position, &
               current_last_position
    character(4) :: temporary_plab
    character(4), dimension(:), allocatable :: min_plab_curr
    complex(8) :: temporary_freq
    complex(8), dimension(:), allocatable :: current_minimum_freq

! prn = 0

! if (num_p_tuples == 3) then
! 
! if (p_tuples(2)%npert == 2 .AND. p_tuples(3)%npert == 2) then

! if ((p_tuples(2)%pid(1) == 6) .AND. (p_tuples(2)%pid(1) == 6) &
!      .AND. (p_tuples(3)%pid(1) == 1) .AND. (p_tuples(3)%pid(2) == 10)) then

! prn = 1

! write(*,*) ' '
! write(*,*) 'incoming'
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

! end if

! end if
! 
! end if

    do i = 1, num_p_tuples

       call p1_cloneto_p2(p_tuples(i), p_tuples_st(i))

    end do

    do i = 2, num_p_tuples

       new_minimum = i

       do j = i + 1, num_p_tuples

          if (p_tuple_p1_lt_p2(p_tuple_standardorder(p_tuples_st(j)), &
              p_tuple_standardorder(p_tuples_st(new_minimum)))) then

! if (prn == 1) then
! 
! write(*,*) 'changing places:', j, 'and', new_minimum
! 
! end if

             new_minimum = j

          end if

       end do


       if (.NOT.(new_minimum == i)) then

          call p1_cloneto_p2(p_tuples_st(new_minimum),temporary_pert)
          call p_tuple_deallocate(p_tuples_st(new_minimum))

          p_tuples_st(i) = p_tuple_standardorder(p_tuples_st(i))

          call p1_cloneto_p2(p_tuples_st(i),p_tuples_st(new_minimum))
          call p_tuple_deallocate(p_tuples_st(i))

          temporary_pert = p_tuple_standardorder(temporary_pert)

          call p1_cloneto_p2(temporary_pert,p_tuples_st(i))
          call p_tuple_deallocate(temporary_pert)

       end if

    end do



! if (prn == 1) then
! 
! write(*,*) ' '
! write(*,*) 'outgoing'
! 
! do i = 1, num_p_tuples
! 
! write(*,*) p_tuples_st(i)%pid
! write(*,*) p_tuples_st(i)%pdim
! write(*,*) p_tuples_st(i)%freq
! 
! end do
! 
! 
! write(*,*) ' '
! 
! end if

  end function

  function p_tuples_compare(num_p_tuples, p_tuples, p_tuples_st_order)

    implicit none

    logical :: p_tuples_compare, elem_by_elem_isequivalent
    integer ::  num_p_tuples, i
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, p_tuples_st_order

    p_tuples_compare = .FALSE.
    elem_by_elem_isequivalent = .TRUE.

    do i = 1, num_p_tuples

!        write(*,*) 'plab a', p_tuples(i)%plab
!        write(*,*) 'plab b', p_tuples_st_order(i)%plab
!     
!        write(*,*) 'Comparison:', p_tuple_compare(p_tuple_standardorder(p_tuples(i)), &
!                                    p_tuple_standardorder(p_tuples_st_order(i)))
    
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
  
  function pid_compare(npert, pid1, pid2)
  
    logical :: pid_compare
    integer :: npert, i
    integer, dimension(npert) :: pid1, pid2
    
    pid_compare = .TRUE.
    
    do i = 1, npert
    
       if (.NOT.(pid1(i) == pid2(i))) then
       
          pid_compare = .FALSE.
       
       end if
    
    end do
    
  end function
    
    
  


  function p_tuple_compare(p1, p2)

    implicit none

    logical :: p_tuple_compare
    type(p_tuple) :: p1, p2

    if (p1%npert == p2%npert) then

       if (plab_compare(p1%npert, p1%plab, p2%plab) .eqv. .TRUE.) then

          if (pfreq_compare(p1%npert, p1%freq, p2%freq) .eqv. .TRUE.) then
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
    type(p_tuple), dimension(pert%npert) :: psub
    integer :: i, j, k, m

    do i = 1, pert%npert

       psub(i)%npert = pert%npert - 1

       allocate(psub(i)%pdim(pert%npert - 1))
       allocate(psub(i)%plab(pert%npert - 1))
       allocate(psub(i)%pid(pert%npert - 1))
       allocate(psub(i)%freq(pert%npert - 1))

       k = i
       m = 1

       do j = 1, (pert%npert)

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
