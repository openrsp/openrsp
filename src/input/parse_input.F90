   subroutine read_openrsp_input(file_name, file_unit)

      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(*), intent(in)  :: file_name
      integer,      intent(in)  :: file_unit
!     --------------------------------------------------------------------------
      character(kw_length)      :: word
      character(kw_length)      :: kw_section
      logical                   :: file_is_open
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

      inquire(file=file_name, opened=file_is_open)
      if (.not. file_is_open) then
         open(file_unit, file=file_name)
      end if

      rewind(file_unit)
      call set_file_unit(file_unit)

      do while (.true.)
         read(file_unit, '(a7)', end=1) word
         if (word == '       ') then
!            blank line
         else
            select case (word(1:1))
               case ('!', '#')
!                 comment
               case ('*')
!                 section
                  kw_section = uppercase(word)
                  if (word == '*END OF' .or. word == '**END O') then
                     exit
                  end if
               case default
!                 keyword
                  if (kw_section == '**OPENR') then
                     call read_openrsp_section(word, kw_section)
                  end if
            end select
         end if
      end do

1     if (.not. file_is_open) then
         close(file_unit, status='keep')
      end if

   end subroutine

   subroutine read_openrsp_section(word, kw_section)

      use input_reader
      use openrsp_cfg

      implicit none

!     ----------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ----------------------------------------------------------------------------
      integer                          :: i
!     ----------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, openrsp_cfg_print_level)
      end if

      if (kw_matches(word, '.MAXITR')) then
         call kw_read(word, openrsp_cfg_solver_maxitr)
      end if

      if (kw_matches(word, '.MAXPHP')) then
         call kw_read(word, openrsp_cfg_solver_maxphp)
      end if

      if (kw_matches(word, '.MAXRED')) then
         call kw_read(word, openrsp_cfg_solver_maxred)
      end if

      if (kw_matches(word, '.THRESH')) then
         call kw_read(word, openrsp_cfg_solver_thresh)
      end if

      if (kw_matches(word, '.OPTORB')) then
         openrsp_cfg_solver_optorb = .true.
      end if

      if (kw_matches(word, '.FREQ  ')) then
         call kw_read(word, openrsp_cfg_nr_real_freqs)
         allocate(openrsp_cfg_real_freqs(openrsp_cfg_nr_real_freqs))
         do i = 1, openrsp_cfg_nr_real_freqs
            call kw_read(word, openrsp_cfg_real_freqs(i))
         end do
      end if

      if (kw_matches(word, '.IMFREQ')) then
         call kw_read(word, openrsp_cfg_nr_imag_freqs)
         allocate(openrsp_cfg_imag_freqs(openrsp_cfg_nr_imag_freqs))
         do i = 1, openrsp_cfg_nr_imag_freqs
            call kw_read(word, openrsp_cfg_imag_freqs(i))
         end do
      end if

      if (kw_matches(word, '.CVALUE')) then
         call kw_read(word, openrsp_cfg_speed_of_light)
      end if

      if (kw_matches(word, '.NOLLSS')) then
         openrsp_cfg_skip_llss = .true.
      end if

      if (kw_matches(word, '.GRADIE')) then
         openrsp_cfg_gradient = .true.
      end if

      if (kw_matches(word, '.HESSIA')) then
         openrsp_cfg_hessian = .true.
      end if

      if (kw_matches(word, '.DIPHES')) then
         openrsp_cfg_dipole_hessian = .true.
      end if

      if (kw_matches(word, '.CUBICF')) then
         openrsp_cfg_cubic_ff = .true.
      end if

      if (kw_matches(word, '.QUARTI')) then
         openrsp_cfg_quartic_ff = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.F     ')) then
         openrsp_cfg_general_f = .true.
      end if

      if (kw_matches(word, '.G     ')) then
         openrsp_cfg_general_g = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.FF    ')) then
         openrsp_cfg_general_ff = .true.
      end if

      if (kw_matches(word, '.GF    ')) then
         openrsp_cfg_general_gf = .true.
      end if

      if (kw_matches(word, '.GG    ')) then
         openrsp_cfg_general_gg = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.FFF   ')) then
         openrsp_cfg_general_fff = .true.
      end if

      if (kw_matches(word, '.GFF   ')) then
         openrsp_cfg_general_gff = .true.
      end if

      if (kw_matches(word, '.GGF   ')) then
         openrsp_cfg_general_ggf = .true.
      end if

      if (kw_matches(word, '.GGG   ')) then
         openrsp_cfg_general_ggg = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.FFFF  ')) then
         openrsp_cfg_general_ffff = .true.
      end if

      if (kw_matches(word, '.GFFF  ')) then
         openrsp_cfg_general_gfff = .true.
      end if

      if (kw_matches(word, '.GGFF  ')) then
         openrsp_cfg_general_ggff = .true.
      end if

      if (kw_matches(word, '.GGGF  ')) then
         openrsp_cfg_general_gggf = .true.
      end if

      if (kw_matches(word, '.GGGG  ')) then
         openrsp_cfg_general_gggg = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.FFFFF ')) then
         openrsp_cfg_general_fffff = .true.
      end if

      if (kw_matches(word, '.gFFFF ')) then
         openrsp_cfg_general_gffff = .true.
      end if

      if (kw_matches(word, '.GGFFF ')) then
         openrsp_cfg_general_ggfff = .true.
      end if

      if (kw_matches(word, '.GGGFF ')) then
         openrsp_cfg_general_gggff = .true.
      end if

      if (kw_matches(word, '.GGGGF ')) then
         openrsp_cfg_general_ggggf = .true.
      end if

      if (kw_matches(word, '.GGGGG ')) then
         openrsp_cfg_general_ggggg = .true.
      end if

!     --------------------------------------------------------------------------

      if (kw_matches(word, '.FFFFFF')) then
         openrsp_cfg_general_ffffff = .true.
      end if

      if (kw_matches(word, '.GFFFFF')) then
         openrsp_cfg_general_gfffff = .true.
      end if

      if (kw_matches(word, '.GGFFFF')) then
         openrsp_cfg_general_ggffff = .true.
      end if

      if (kw_matches(word, '.GGGFFF')) then
         openrsp_cfg_general_gggfff = .true.
      end if

      if (kw_matches(word, '.GGGGFF')) then
         openrsp_cfg_general_ggggff = .true.
      end if

      if (kw_matches(word, '.GGGGGF')) then
         openrsp_cfg_general_gggggf = .true.
      end if

      if (kw_matches(word, '.GGGGGG')) then
         openrsp_cfg_general_gggggg = .true.
      end if

!     --------------------------------------------------------------------------

      call check_whether_kw_found(word, kw_section)

   end subroutine
