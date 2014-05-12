   subroutine read_openrsp_input(file_name, file_unit)

      use openrsp_input_reader

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

      use openrsp_input_reader
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

      ! reads the broadening (damping) parameter, as such the complex
      ! polarization propagator solver will be used
      if (kw_matches(word, '.DAMPING')) then
         openrsp_cfg_cpp_used = .true.
         call kw_read(word, openrsp_cfg_cpp_damping)
      end if

      if (kw_matches(word, '.FREQ  ')) then
         call kw_read(word, openrsp_cfg_nr_real_freqs)
         allocate(openrsp_cfg_real_freqs(openrsp_cfg_nr_real_freqs))
         do i = 1, openrsp_cfg_nr_real_freqs
            call kw_read(word, openrsp_cfg_real_freqs(i))
         end do
      end if

      ! MaR: Support for multiple frequency tuples
      ! Extend this support to ordinary rsp property calculations

      if (kw_matches(word, '.MFREQS')) then
         call kw_read(word, openrsp_cfg_nr_freq_tuples)
         call kw_read(word, openrsp_cfg_nr_real_freqs) ! Here used as freqs. per tuple
         allocate(openrsp_cfg_real_freqs(openrsp_cfg_nr_freq_tuples * &
                                         openrsp_cfg_nr_real_freqs))
         do i = 1, openrsp_cfg_nr_freq_tuples
               read(get_file_unit(), *) openrsp_cfg_real_freqs( &
               (i - 1) * openrsp_cfg_nr_real_freqs + 1 : &
               i * openrsp_cfg_nr_real_freqs)
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

      if (kw_matches(word, '.DOSSSS')) then
         openrsp_cfg_skip_ssss = .false.
      end if

      if (kw_matches(word, '.MAGNET')) then
         openrsp_cfg_magnetizability = .true.
      end if

      if (kw_matches(word, '.VCD   ')) then
         openrsp_cfg_vcd = .true.
      end if

      ! MaR: For entering temperature in Kelvin

      if (kw_matches(word, '.TMPRTR')) then
         call kw_read(word, openrsp_cfg_temperature)
      end if


!     --------------------------------------------------------------------------

      ! MaR: Keywords for custom property specification

      if (kw_matches(word, '.CUSTOM')) then
         openrsp_cfg_general_specify = .true.
      end if

      if (kw_matches(word, '.SPORDR ')) then
         call kw_read(word, openrsp_cfg_specify_order)
      end if

      if (kw_matches(word, '.SPRULE ')) then
         call kw_read(word, openrsp_cfg_specify_kn(1))
         call kw_read(word, openrsp_cfg_specify_kn(2))
      end if

      if (kw_matches(word, '.SPPLAB ')) then
         allocate(openrsp_cfg_specify_plab(openrsp_cfg_specify_order))
         do i = 1, openrsp_cfg_specify_order
            call kw_read(word, openrsp_cfg_specify_plab(i))
         end do
      end if

!     --------------------------------------------------------------------------

      ! MaR: Keywords for pure vibrational contribution calculations

      ! PV contribution to alpha

      if (kw_matches(word, '.PV2FLD')) then
         openrsp_cfg_general_pv2f = .true.
      end if

      ! PV contribution to beta

      if (kw_matches(word, '.PV3FLD')) then
         openrsp_cfg_general_pv3f = .true.
      end if

      ! PV contribution to gamma

      if (kw_matches(word, '.PV4FLD')) then
         openrsp_cfg_general_pv4f = .true.
      end if

      ! Order of electrical anharmonicity (default: 0)

      if (kw_matches(word, '.PVEANH')) then
         call kw_read(word, openrsp_cfg_general_pv_el_anh)
      end if

      ! Order of mechanical anharmonicity (default: 0)

      if (kw_matches(word, '.PVMANH')) then
         call kw_read(word, openrsp_cfg_general_pv_mech_anh)
      end if

      ! Order of total (both electric and mechanical) anharmonicity (default: 0)
      ! Overrides choices of electrical and mechanical anharmonicity above

      if (kw_matches(word, '.PVTANH')) then
         call kw_read(word, openrsp_cfg_general_pv_total_anh)
         openrsp_cfg_general_pv_el_anh = openrsp_cfg_general_pv_total_anh
         openrsp_cfg_general_pv_mech_anh = openrsp_cfg_general_pv_total_anh
      end if

!     --------------------------------------------------------------------------

      ! MaR: Keywords for ZPVA contribution calculations

      ! ZPVA contribution to alpha

      if (kw_matches(word, '.ZPVA2F')) then
         openrsp_cfg_general_zpva2f = .true.
      end if

      ! ZPVA contribution to beta

      if (kw_matches(word, '.ZPVA3F')) then
         openrsp_cfg_general_zpva3f = .true.
      end if

      ! ZPVA contribution to gamma

      if (kw_matches(word, '.ZPVA4F')) then
         openrsp_cfg_general_zpva4f = .true.
      end if

!     --------------------------------------------------------------------------

      ! MaR: Keywords for SHG (and later possibly higher order) calculations

      if (kw_matches(word, '.RSPSHG')) then
         openrsp_cfg_general_shg = .true.
      end if

!     --------------------------------------------------------------------------

      ! MaR: Keywords for EFISHG-CID calculation

      if (kw_matches(word, '.EFSHCI')) then
         openrsp_cfg_general_efishgcid = .true.
      end if


!     --------------------------------------------------------------------------

      ! MaR: Keywords for SFG calculation

      if (kw_matches(word, '.RSPSFG')) then
         openrsp_cfg_general_sfg = .true.
      end if


!     --------------------------------------------------------------------------

      ! MaR: Keywords for hyper-Raman calculation

      if (kw_matches(word, '.HYPRAM')) then
         openrsp_cfg_general_hyper_raman = .true.
      end if

      ! MaR: Normal mode frequency scaling (e.g. to account for
      ! lack of correlation or anharmonicity)

      if (kw_matches(word, '.FSCALE')) then
         call kw_read(word, openrsp_cfg_general_hypram_freqscale)
      end if

!     --------------------------------------------------------------------------

      ! MaR: Keywords for Cartesian to normal mode transformation routines

      if (kw_matches(word, '.CARTNC')) then
         openrsp_cfg_general_trans_cartnc = .true.
         call kw_read(word, openrsp_cfg_general_cartnc_order_geo)
         call kw_read(word, openrsp_cfg_general_cartnc_order_field)
      end if

!     --------------------------------------------------------------------------

      ! keywords that control the XCint grid (not the Dalton grid)
      if (kw_matches(word, '.XCGRID')) then
         openrsp_cfg_use_xcint_grid = .true.
         call kw_read(word, openrsp_cfg_radint)
         call kw_read(word, openrsp_cfg_angmin)
         call kw_read(word, openrsp_cfg_angint)
      end if

!     --------------------------------------------------------------------------

      call check_whether_kw_found(word, kw_section)

   end subroutine
