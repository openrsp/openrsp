
! (c) Radovan Bast and Stefan Knecht
! licensed under the GNU Lesser General Public License

module input_reader

   implicit none

   public lowercase
   public uppercase
   public word_contains
   public reset_available_kw_list
   public check_whether_kw_found
   public kw_matches
   public kw_read
   public set_file_unit
   public get_file_unit

   integer, parameter, public :: kw_length = 7

   private

   interface kw_read
      module procedure kw_read_c
      module procedure kw_read_i1
      module procedure kw_read_i3
      module procedure kw_read_r1
      module procedure kw_read_r2
      module procedure kw_read_r3
      module procedure kw_read_r4
      module procedure kw_read_ivec
   end interface

   integer, parameter :: max_nr_kw   = 200
   integer, parameter :: line_length = 80

   integer :: nr_available_kw
   integer :: file_unit
   logical :: kw_found

   character(kw_length) :: available_kw_list(max_nr_kw)

contains

   subroutine set_file_unit(u)
      integer, intent(in) :: u
      file_unit = u
   end subroutine

   integer function get_file_unit()
      get_file_unit = file_unit
   end function

   function lowercase(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len(s))        :: lowercase
!     --------------------------------------------------------------------------
      integer                  :: off, i, ia
!     --------------------------------------------------------------------------

      lowercase = s

      off = iachar('a') - iachar('A')

      do i = 1, len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('A') .and. ia <= iachar('Z')) then
            lowercase(i:i) = achar(ia + off)
         end if
      enddo

   end function

   function uppercase(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len(s))        :: uppercase
!     --------------------------------------------------------------------------
      integer                  :: off, i, ia
!     --------------------------------------------------------------------------

      uppercase = s

      off = iachar('A') - iachar('a')

      do i = 1, len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('a') .and. ia <= iachar('Z')) then
            uppercase(i:i) = achar(ia + off)
         end if
      enddo

   end function

   function prefix_zeros(i, n)

!     prefix_zeros(137, 6) returns '000137'

!     --------------------------------------------------------------------------
      integer,      intent(in) :: i
      integer,      intent(in) :: n
      character(n)             :: prefix_zeros
!     --------------------------------------------------------------------------
      integer                  :: k
      character(1)             :: c09(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
!     --------------------------------------------------------------------------

      do k = 1, n
         prefix_zeros(n-k+1:n-k+1) = c09(mod(i, 10**k)/10**(k-1))
      end do

   end function

   function word_count(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      integer                  :: word_count
!     --------------------------------------------------------------------------
      integer                  :: i
      logical                  :: is_blank
!     --------------------------------------------------------------------------

      word_count = 0

      if (len(s) <= 0) return

      is_blank = .true.

      do i = 1, len(s)
         if (s(i:i) == ' ') then
            is_blank = .true.
         else if (is_blank) then
            word_count = word_count + 1
            is_blank = .false.
         end if
      end do

   end function

   function word_contains(word, substring)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: word
      character(*), intent(in) :: substring
!     --------------------------------------------------------------------------
      logical                  :: word_contains
!     --------------------------------------------------------------------------

      word_contains = .false.
      if (index(word, substring) > 0) then
         word_contains = .true.
      end if

   end function

   subroutine kw_read_c(kw_input, c)

      character(kw_length), intent(in)  :: kw_input
      character(*),         intent(out) :: c

      read(file_unit, *, err=1) c
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_i1(kw_input, i)

      character(kw_length), intent(in)  :: kw_input
      integer,              intent(out) :: i

      read(file_unit, *, err=1) i
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_i3(kw_input, i1, i2, i3)

      character(kw_length), intent(in)  :: kw_input
      integer,              intent(out) :: i1, i2, i3

      read(file_unit, *, err=1) i1, i2, i3
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_r1(kw_input, r)

      character(kw_length), intent(in)  :: kw_input
      real(8),              intent(out) :: r

      read(file_unit, *, err=1) r
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_r2(kw_input, r1, r2)

      character(kw_length), intent(in)  :: kw_input
      real(8),              intent(out) :: r1, r2

      read(file_unit, *, err=1) r1, r2
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_r3(kw_input, r1, r2, r3)

      character(kw_length), intent(in)  :: kw_input
      real(8),              intent(out) :: r1, r2, r3

      read(file_unit, *, err=1) r1, r2, r3
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_r4(kw_input, r1, r2, r3, r4)

      character(kw_length), intent(in)  :: kw_input
      real(8),              intent(out) :: r1, r2, r3, r4

      read(file_unit, *, err=1) r1, r2, r3, r4
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_ivec(kw_input, v, k)

      character(kw_length), intent(in)  :: kw_input
      integer,              intent(out) :: v(:)
      integer,              intent(in)  :: k
      integer                           :: j

      read(file_unit, *,err=1) (v(j), j=1,k)
      return

1     call kw_read_error(kw_input)

   end subroutine

   subroutine kw_read_error(kw_input)

      character(kw_length), intent(in) :: kw_input
      character(line_length)           :: line

      backspace file_unit
      read(file_unit, *) line

      write(*, *) 'error in input line:'
      write(*, *) line
      write(*, *) 'following keyword '//kw_input

      call quit('error in line following keyword '//kw_input)

   end subroutine

   subroutine reset_available_kw_list()

      nr_available_kw = 0
      kw_found        = .false.

   end subroutine

   subroutine check_whether_kw_found(kw_input, kw_section)

      character(kw_length), intent(in) :: kw_input
      character(kw_length), intent(in) :: kw_section
      integer                          :: i

      if (.not. kw_found) then

         write(*, *) 'illegal keyword '//kw_input//' in section '//kw_section

         write(*, *) 'list of available keywords in section '//kw_section//':'
         do i = 1, nr_available_kw
            write(*, *) available_kw_list(i)
         end do

         call quit('illegal keyword '//kw_input//' in section '//kw_section)

      end if

   end subroutine

   function kw_matches(kw_input, kw_option)

      character(kw_length), intent(in) :: kw_input
      character(kw_length), intent(in) :: kw_option
      logical                          :: kw_matches

      if (lowercase(kw_input) == lowercase(kw_option)) then
         kw_matches = .true.
         kw_found   = .true.
      else
         kw_matches = .false.
         nr_available_kw = nr_available_kw + 1
         available_kw_list(nr_available_kw) = kw_option
      end if

   end function

end module
