! Copyright 2012 Magnus Ringholm
!
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

! Contains the definition and code for manipulation 
! of the field tuple datatype.

module rsp_field_tuple

  implicit none

  public p_tuple_extend
  public p_tuple_getone
  public p_tuple_remove_first
  public merge_p_tuple
  public p1_cloneto_p2
  public p_tuple_add_stateinfo
  public p_tuple_deallocate
  public get_emptypert
  public recognize_contribution
  public find_residue_info
  public find_complete_residualization
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

  
  public p_tuple_to_external_orders
  public p_tuple_to_external_tuple
  
  
  
  public p_tuple_dealloc
  public p1_lt_p2
  public get_p_subset
  public p1_switch_p2
  public p1_merge_p2
  public make_p_tuple_subsets
  
  ! Define perturbation tuple datatype

  type p_tuple

     integer :: npert

     integer, allocatable, dimension(:) :: pdim ! Dimensions of perturbations
     character(4), allocatable, dimension(:) :: plab ! Perturbation labels
     integer, allocatable, dimension(:) :: pfcomp ! First component is component #pfcomp
     integer, allocatable, dimension(:) :: pid ! Pert. ID - for k,n rule evaluations
     complex(8), allocatable, dimension(:) :: freq ! Frequencies of perturbations
     
     integer :: do_residues=0 ! Degree of residue to be calculated (0: no residue)
     integer :: n_pert_res_max=1  ! Max. num. of pert. in residues
     integer :: n_states=1 ! number of states (for residues only)
     integer, allocatable, dimension(:) :: states ! indices of states involved
     complex(8), allocatable, dimension(:) :: exenerg ! Excitation energies (for residues) 
     logical, allocatable, dimension(:,:) :: part_of_residue
     ! Add other perturbation identification info as needed

  end type

  contains

  
  ! Make array of integers to transform character labels to integer labelling of perturbations
  ! This is a temporary solution: The long-term scheme should inherit numbers from host
  ! to designate the perturbation labels and  do away with character-type labels
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
          j = j + 1
          perturbations(j) = i
          pert_orders(j) = count(plab_int == i)
       
       end if
    
    end do
      
  end subroutine

  ! Another transformation routine to be made redundant with long-term scheme
  subroutine p_tuple_to_external_tuple(pert, num_pert, perturbations)
  
    implicit none
    
    type(p_tuple) :: pert
    integer :: num_pert, i
    integer, dimension(:), allocatable :: perturbations
    
    num_pert = pert%npert
    allocate(perturbations(num_pert))

    do i = 1, pert%npert
    
       if (pert%plab(i) == 'NULL') then
       
          perturbations(i) = 0
    
       else if (pert%plab(i) == 'GEO ') then
       
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
      if (allocated(p1%exenerg)) then
         deallocate(p1%exenerg)
         deallocate(p1%part_of_residue)
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
    
    
    ! Get a general subset of perturbation tuple p1   
    ! Doesn't seem to be currently in use, maybe reintroduce later
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
  ! Doesn't seem to be currently in use, maybe reintroduce later
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
      
    end subroutine

    
    ! Make all (p1%npert - 1) subsets of perturbation tuple p1 
    ! Doesn't seem to be currently in use, maybe reintroduce later
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
      
  
  ! Extend perturbation tuple pert by small extension ext
  function p_tuple_extend(pert, ext)

    implicit none

    type(p_tuple) :: pert, ext, p_tuple_extend
    integer :: i, j

    allocate(p_tuple_extend%pdim(pert%npert + 1))
    allocate(p_tuple_extend%plab(pert%npert + 1))
    allocate(p_tuple_extend%pid(pert%npert + 1))
    allocate(p_tuple_extend%freq(pert%npert + 1))
    call p_tuple_add_stateinfo(p_tuple_extend,pert)
    if (pert%do_residues.gt.0) then
      allocate(p_tuple_extend%part_of_residue(pert%npert + 1, pert%do_residues))
    end if
    if (pert%npert == 0) then

       p_tuple_extend%npert = pert%npert + 1
       p_tuple_extend%pdim = (/ext%pdim(:)/)
       p_tuple_extend%plab = (/ext%plab(:)/)
       p_tuple_extend%pid = (/ext%pid(:)/)
       p_tuple_extend%freq = (/ext%freq(:)/)
       if(allocated(p_tuple_extend%part_of_residue)) then
         p_tuple_extend%part_of_residue(:,:)= ext%part_of_residue(:,:)
       end if


    else

       p_tuple_extend%npert = pert%npert + 1
       p_tuple_extend%pdim = (/(pert%pdim(i), i = 1, pert%npert), ext%pdim(:)/)
       p_tuple_extend%plab = (/(pert%plab(i), i = 1, pert%npert), ext%plab(:)/)
       p_tuple_extend%pid = (/(pert%pid(i), i = 1, pert%npert), ext%pid(:)/)
       p_tuple_extend%freq = (/(pert%freq(i), i = 1, pert%npert), ext%freq(:)/)
       if(allocated(p_tuple_extend%part_of_residue)) then
          do j = 1, pert%do_residues
            p_tuple_extend%part_of_residue(:,j) = (/(pert%part_of_residue(i,j), i = 1, pert%npert), &
                ext%part_of_residue(:,j)/)
          end do
       end if
 
end if

  end function


  ! Get perturbation 'which' from perturbation tuple pert
  function p_tuple_getone(pert, which, istate)

    implicit none

    type(p_tuple) :: pert, p_tuple_getone
    integer :: which
    integer, optional :: istate

    allocate(p_tuple_getone%pdim(1))
    allocate(p_tuple_getone%plab(1))
    allocate(p_tuple_getone%pid(1))
    allocate(p_tuple_getone%freq(1))
    call p_tuple_add_stateinfo(p_tuple_getone,pert)
    allocate(p_tuple_getone%part_of_residue(1,pert%do_residues))

    p_tuple_getone%npert = 1
    p_tuple_getone%pdim = (/pert%pdim(which)/)
    p_tuple_getone%plab = (/pert%plab(which)/)
    p_tuple_getone%pid = (/pert%pid(which)/)
    p_tuple_getone%freq = (/pert%freq(which)/)
    if(pert%do_residues.gt.0) then
      p_tuple_getone%part_of_residue(1,1:pert%do_residues) = pert%part_of_residue(which,1:pert%do_residues)
    end if

  end function

  ! Remove the first perturbation from tuple pert
  function p_tuple_remove_first(pert)

    implicit none

    type(p_tuple) :: pert, p_tuple_remove_first

    allocate(p_tuple_remove_first%pdim(pert%npert - 1))
    allocate(p_tuple_remove_first%plab(pert%npert - 1))
    allocate(p_tuple_remove_first%pid(pert%npert - 1))
    allocate(p_tuple_remove_first%freq(pert%npert - 1))
    if (pert%do_residues.gt.0) then

      allocate(p_tuple_remove_first%part_of_residue(pert%npert - 1,pert%do_residues))
      p_tuple_remove_first%part_of_residue = .false.

    end if

    if (pert%npert > 1) then

       p_tuple_remove_first%npert = pert%npert - 1
       p_tuple_remove_first%pdim = (/pert%pdim(2:pert%npert)/)
       p_tuple_remove_first%plab = (/pert%plab(2:pert%npert)/)
       p_tuple_remove_first%pid = (/pert%pid(2:pert%npert)/)
       p_tuple_remove_first%freq = (/pert%freq(2:pert%npert)/)
       if(pert%do_residues.gt.0) then 
         p_tuple_remove_first%part_of_residue(:,1) = (/pert%part_of_residue(2:pert%npert,1)/)
         p_tuple_remove_first%part_of_residue(:,pert%do_residues) = (/pert%part_of_residue(2:pert%npert,pert%do_residues)/)
       end if

    else

       p_tuple_remove_first%npert = 0
       p_tuple_remove_first%pdim = (/0/)
       p_tuple_remove_first%plab = (/'NUTN'/)
       p_tuple_remove_first%pid = (/0/)
       p_tuple_remove_first%freq = (/0.0/)

    end if

    if (pert%do_residues.gt.0) then

       call p_tuple_add_stateinfo(p_tuple_remove_first,pert)
  
    end if

  end function

  recursive subroutine recognize_exenerg(pert,n,val,idx_preserve,energies,found_count)

  ! Recognizes whether the term descibed on pert depends on the perturbations 
  ! whose frequency tends to the excitation energy in a residue calculation

    implicit none

    integer found_count,integerout
    type(p_tuple) pert
    integer n, idx, i, j, idx_preserve,idx_forward
    real(8) val, energies(pert%do_residues)
    real(8), parameter :: xtiny=1.0d-6

    idx = idx_preserve
    loop: do i= idx, pert%npert - n + 1

        val = val + dble(pert%freq(i))
        if (n.gt.1) then
          idx_forward = i+1
          call recognize_exenerg(pert,n-1,val,idx_forward,energies,found_count)
        else if (n.eq.1) then
          do j = 1, pert%do_residues
            if(energies(j).gt.0) then
              if (abs(abs(val)-dble(energies(j))).lt.xtiny) then
                 found_count = found_count + 1
                 val = 0.0d0
                 energies(j) = -energies(j)
              end if
            end if
          end do
        else
          write(*,*) 'Index error in recognize_exenerg'
        end if
        val = 0.0d0
   end do loop

  end subroutine

  logical function recognize_contribution(pert,n_exenerg)

    integer i,n_exenerg,idx,jdx,found_count
    type(p_tuple) pert
    logical found_exenerg(2)
    real*8 energies(2),val

    if (pert%do_residues.eq.0) then
       recognize_contribution = .true.
       return
    end if

    idx = 1
    recognize_contribution=.false.
    energies = -9999.9d9
    found_count = 0
    val = 0.0d0

    do jdx = 1, pert%n_states
       energies(jdx) = pert%exenerg(jdx)
    end do

    call recognize_exenerg(pert,pert%n_pert_res_max,val,idx,energies,found_count)

    if (found_count .eq. n_exenerg) then
       recognize_contribution = .true.
    else
       recognize_contribution = .false.
    end if
    return
  end function recognize_contribution

  logical function find_residue_info(pert)

    implicit none
    type(p_tuple) :: pert
    integer :: i

    find_residue_info=.false.
    if(pert%do_residues.eq.0) then
      return
    end if

    do i = 1, pert%npert
       if (pert%part_of_residue(i,1).or.pert%part_of_residue(i,pert%do_residues)) then
          find_residue_info = .true.
          return
       end if
    end do
    return

  end function find_residue_info

  logical function find_complete_residualization(pert)

    implicit none
    type(p_tuple) :: pert
    logical :: lresidue(2)
    integer :: i,j

    if(pert%do_residues.eq.0) then
      find_complete_residualization = .false.
      return
    end if

    lresidue = .true.
    do i = 1, pert%npert
      do j = 1, pert%do_residues
        lresidue(j) = lresidue(j).and.pert%part_of_residue(i,j)
      end do
    end do

    find_complete_residualization = lresidue(1).or.lresidue(pert%do_residues)
    return

  end function find_complete_residualization

  
  ! Return a perturbation tuple which is the merged tuple of tuples p1 and p2
  function merge_p_tuple(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2, merge_p_tuple
    integer :: i

    allocate(merge_p_tuple%pdim(p1%npert + p2%npert))
    allocate(merge_p_tuple%plab(p1%npert + p2%npert))
    allocate(merge_p_tuple%pid(p1%npert + p2%npert))
    allocate(merge_p_tuple%freq(p1%npert + p2%npert))

    call p_tuple_add_stateinfo(merge_p_tuple,p1)

    if (merge_p_tuple%do_residues.gt.0) then
      allocate(merge_p_tuple%part_of_residue(p1%npert+p2%npert,merge_p_tuple%do_residues))
      merge_p_tuple%part_of_residue = .false.
    end if

    ! NOTE (MaR): ARE THESE CASE DISTINCTIONS UNNECESSARY? CONSIDER REWRITE.

    if (p1%npert > 0 .AND. p2%npert > 0) then

       merge_p_tuple%npert = p1%npert + p2%npert
       merge_p_tuple%pdim = (/p1%pdim(:), p2%pdim(:)/)
       merge_p_tuple%plab = (/p1%plab(:), p2%plab(:)/)
       merge_p_tuple%pid = (/p1%pid(:), p2%pid(:)/)
       merge_p_tuple%freq = (/p1%freq(:), p2%freq(:)/)
       do i = 1, merge_p_tuple%do_residues
          merge_p_tuple%part_of_residue(:,i) = (/p1%part_of_residue(:,i),p2%part_of_residue(:,i)/)
       end do

    elseif (p1%npert > 0 .AND. p2%npert == 0) then

       merge_p_tuple%npert = p1%npert
       merge_p_tuple%pdim = p1%pdim(:)
       merge_p_tuple%plab = p1%plab(:)
       merge_p_tuple%pid = p1%pid(:)
       merge_p_tuple%freq = p1%freq(:)
       do i = 1, merge_p_tuple%do_residues
          merge_p_tuple%part_of_residue(:,i) = (/p1%part_of_residue(:,i)/)
       end do


    elseif (p1%npert == 0 .AND. p2%npert > 0) then

       merge_p_tuple%npert = p2%npert
       merge_p_tuple%pdim = p2%pdim(:)
       merge_p_tuple%plab = p2%plab(:)
       merge_p_tuple%pid = p2%pid(:)
       merge_p_tuple%freq = p2%freq(:)
       do i = 1, merge_p_tuple%do_residues
          merge_p_tuple%part_of_residue(:,i) = (/p2%part_of_residue(:,i)/)
       end do


    elseif (p1%npert == 0 .AND. p2%npert == 0) then

       merge_p_tuple = get_emptypert()

    else

       write(*,*) 'Error in merge_p_tuple: Unrecognized size of p1 or p2 or both:', &
                   p1%npert, p2%npert

    end if

  end function

  subroutine p_tuple_add_stateinfo(emptytuple,tuple)

    implicit none
    type(p_tuple) :: emptytuple,tuple
    integer i

    if (tuple%do_residues.eq.0) then
       return
    end if

    if (tuple%npert > 0) then

    emptytuple%n_states = tuple%n_states
    emptytuple%n_pert_res_max = tuple%n_pert_res_max
    emptytuple%do_residues = tuple%do_residues

    end if

    allocate(emptytuple%states(emptytuple%do_residues))
    allocate(emptytuple%exenerg(emptytuple%do_residues))

    do i = 1, emptytuple%do_residues

      if (tuple%npert > 0) then

         emptytuple%states(i)=tuple%states(i)
         emptytuple%exenerg(i)=tuple%exenerg(i)

      end if

    end do

  end subroutine

  ! Take perturbation tuple p1 and clone it into p2
  subroutine p1_cloneto_p2(p1, p2)

    implicit none

    type(p_tuple) :: p1, p2

    p2%npert = p1%npert

    allocate(p2%pdim(p1%npert)) 
    allocate(p2%plab(p1%npert))
    allocate(p2%pid(p1%npert))
    allocate(p2%freq(p1%npert))

    if (p1%npert > 0) then

       p2%pdim = p1%pdim
       p2%plab = p1%plab
       p2%pid = p1%pid
       p2%freq = p1%freq

       p2%do_residues = p1%do_residues

    end if

    call p_tuple_add_stateinfo(p2,p1)

    if (p2%do_residues.gt.0) then
      allocate(p2%part_of_residue(p2%npert,p2%do_residues))
 
      if (p1%npert > 0) then

         p2%part_of_residue = p1%part_of_residue

      end if      

    end if

  end subroutine

  ! Deallocate perturbation tuple p1
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
    
    if (allocated(p1%states)) then
    
       deallocate(p1%states)
       
    end if
    
    if (allocated(p1%exenerg)) then
       
       deallocate(p1%exenerg)
       
    end if
    
    
    if (allocated(p1%part_of_residue)) then
    
       deallocate(p1%part_of_residue)
       
    end if

    
  end subroutine

  ! Return an emtpy perturbation tuple
  function get_emptypert(template) result(emptypert)

    implicit none

      type(p_tuple), optional :: template
      type(p_tuple) :: emptypert

      if (present(template)) then
      
         call p_tuple_add_stateinfo(emptypert,template)
         
      end if
      
      emptypert%npert = 0
      allocate(emptypert%pdim(0))    
      allocate(emptypert%plab(0))
      allocate(emptypert%pid(0))
      allocate(emptypert%freq(0))
      
      if (emptypert%do_residues.gt.0) then
    
         allocate(emptypert%part_of_residue(0,emptypert%do_residues))
         
      end if


  end function

  ! Make emptypert into an empty perturbation tuple
  subroutine empty_p_tuple(emptypert,template)

    implicit none

      type(p_tuple) :: emptypert
      type(p_tuple), optional :: template

      if (present(template)) then 
      
         call p_tuple_add_stateinfo(emptypert,template)
         
      end if
      
      emptypert%npert = 0
      allocate(emptypert%pdim(0))    
      allocate(emptypert%plab(0))
      allocate(emptypert%pid(0))
      allocate(emptypert%freq(0))
      
      
      if (emptypert%do_residues.gt.0) then
      
         allocate(emptypert%part_of_residue(0,emptypert%do_residues))
         
      end if

  end subroutine


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


  
  ! Return a perturbation tuple which is the standard ordered version of tuple pert
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
    logical :: temporary_part_of_residue(2)
    integer :: do_residues


    call p1_cloneto_p2(pert, p_tuple_st)

    position_in_result = 1

    do i = position_in_result, p_tuple_st%npert

       current_minimum_pdim = p_tuple_st%pdim(i)
       current_minimum_plab = p_tuple_st%plab(i)
       current_minimum_freq = p_tuple_st%freq(i)
       new_minimum = i

       do j = i + 1, p_tuple_st%npert

          if (p_tuple_st%pdim(j) >= current_minimum_pdim) then

             if (p_tuple_st%pdim(j) == current_minimum_pdim) then

                if (lge(p_tuple_st%plab(j), current_minimum_plab) .EQV. .TRUE.) then

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

       ! Do the exchange also for the residualization information
       do_residues = p_tuple_st%do_residues
       if (do_residues.gt.0) then
         temporary_part_of_residue(1:do_residues) =  p_tuple_st%part_of_residue(new_minimum,1:do_residues)
         p_tuple_st%part_of_residue(new_minimum,1:do_residues) = p_tuple_st%part_of_residue(i,1:do_residues)
         p_tuple_st%part_of_residue(i,1:do_residues) = temporary_part_of_residue(1:do_residues)
       end if

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

          ! Do the exchange also for the residualization information
          do_residues = p_tuple_st%do_residues
          if (do_residues.gt.0) then
            temporary_part_of_residue(1:do_residues) = p_tuple_st%part_of_residue(new_minimum,1:do_residues)
            p_tuple_st%part_of_residue(new_minimum,1:do_residues) = p_tuple_st%part_of_residue(i,1:do_residues)
            p_tuple_st%part_of_residue(i,1:do_residues) = temporary_part_of_residue(1:do_residues)
          end if

       end do

       current_last_position = current_last_position + 1
       current_first_position = current_last_position

    end do

  end function

  ! Return a tuple of pert tuples which is the standard order of the pert tuple tuple p_tuples
  function p_tuples_standardorder(num_p_tuples, p_tuples) result(p_tuples_st)

    implicit none

    integer :: num_p_tuples
    type(p_tuple), dimension(:), allocatable :: p_tuples_st
    type(p_tuple), dimension(num_p_tuples), intent(inout) :: p_tuples
    type(p_tuple) :: temporary_pert
    integer :: i, j, new_minimum

write(*,*) 'npt', num_p_tuples

    allocate(p_tuples_st(num_p_tuples))

write(*,*) 'allocd'

    do i = 1, num_p_tuples

       call p1_cloneto_p2(p_tuples(i), p_tuples_st(i))

    end do

write(*,*) 'a'

    do i = 2, num_p_tuples

       new_minimum = i

       do j = i + 1, num_p_tuples

          if (p_tuple_p1_lt_p2(p_tuple_standardorder(p_tuples_st(j)), &
              p_tuple_standardorder(p_tuples_st(new_minimum)))) then

             new_minimum = j

          end if

       end do

write(*,*) 'b'

       if (.NOT.(new_minimum == i)) then

          call p1_cloneto_p2(p_tuples_st(new_minimum),temporary_pert)
          call p_tuple_deallocate(p_tuples_st(new_minimum))

write(*,*) 'c'
          p_tuples_st(i) = p_tuple_standardorder(p_tuples_st(i))

          call p1_cloneto_p2(p_tuples_st(i),p_tuples_st(new_minimum))
          call p_tuple_deallocate(p_tuples_st(i))

write(*,*) 'd'
          temporary_pert = p_tuple_standardorder(temporary_pert)

          call p1_cloneto_p2(temporary_pert,p_tuples_st(i))
          call p_tuple_deallocate(temporary_pert)

       end if

    end do

write(*,*) 'e'

  end function

  ! Determine if p_tuples and p_tuples_st_order are equivalent
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

  ! Determine if perturbation label tuples plab1 and plab2 are equivalent
  function plab_compare(npert, plab1, plab2)

    implicit none

    logical :: plab_compare
    integer :: npert, i, j
    character(4) :: plab1(npert), plab2(npert)

    plab_compare = .TRUE.

    do i = 1, npert
       if (.NOT.(plab1(i) == plab2(i))) then
          plab_compare = .FALSE.
       end if
    end do

  end function

  ! Determine if perturbation frequency tuples p1 and p2 are equivalent
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
  
  ! Determine if perturbation ID tuples pid1 and pid2 are equivalent
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

  
  ! Determine if perturbation tuples p1 and p2 are equivalent
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

  ! Make the pert%npert size (pert%npert - 1) subsets of pert and put them in psub
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

       call p_tuple_add_stateinfo(psub(i),pert)
       if (pert%do_residues.gt.0) then
         allocate(psub(i)%part_of_residue(pert%npert - 1,pert%do_residues))
       end if
       k = i
       m = 1

       do j = 1, (pert%npert)

          if (.NOT.(j == k)) then

             psub(i)%pdim(m) = pert%pdim(j)
             psub(i)%plab(m) = pert%plab(j)
             psub(i)%pid(m) = pert%pid(j)
             psub(i)%freq(m) = pert%freq(j)
             if(pert%do_residues.gt.0) then
              psub(i)%part_of_residue(m,1:pert%do_residues) = pert%part_of_residue(j,1:pert%do_residues)
             end if

             m = m + 1

          end if

       end do
    end do
  
  end subroutine

end module
