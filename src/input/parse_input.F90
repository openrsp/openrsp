module openrsp_input

   implicit none

   private

   public read_openrsp_input

contains

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
      
      
      if (kw_matches(word, '.JUSTTR')) then
         openrsp_cfg_general_suppress_calc = .true.
      end if

!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------
      ! DaF: Keywords for residue calculations, mimicing the custom-section
      ! We use the kn-input from the custom property specification

      ! Determine whether single or double residues are to be calced
      ! Double residues not yet to be calculated in this framework!
      if (kw_matches(word, '.RESIDU')) then
         openrsp_cfg_residues = .true.
         call kw_read(word, openrsp_cfg_residue_order)
      end if

      ! Number of perturbations in the corresponding response function
      if (kw_matches(word, '.RENPER')) then
         call kw_read(word, openrsp_cfg_residue_n_pert)
      end if 

      if (kw_matches(word, '.RSPECP')) then
         allocate (openrsp_cfg_residue_spec_pert(openrsp_cfg_residue_order))
         do i = 1, openrsp_cfg_residue_order
           call kw_read(word, openrsp_cfg_residue_spec_pert(i))
         end do
      end if

      if (kw_matches(word, '.RSPCIX')) then
         allocate (openrsp_cfg_residue_spec_index(max(openrsp_cfg_residue_spec_pert(1),&
               openrsp_cfg_residue_spec_pert(openrsp_cfg_residue_order)),openrsp_cfg_residue_order))
         openrsp_cfg_residue_spec_index = -99999
         do i = 1, openrsp_cfg_residue_order
            call kw_read(word,openrsp_cfg_residue_spec_index(:,i),openrsp_cfg_residue_spec_pert(i))
         end do
      end if

      ! Shortcut for TPCD using London orbitals
      if (kw_matches(word, '.RSTPCD')) then
         openrsp_cfg_residue_tpcd = .true.
         openrsp_cfg_residue_n_pert = 3
         openrsp_cfg_residue_order = 1
         openrsp_cfg_specify_kn(1)=0
         openrsp_cfg_specify_kn(2)=2
         allocate(openrsp_cfg_residue_spec_pert(1))
         allocate(openrsp_cfg_residue_spec_index(1,1))
         openrsp_cfg_residue_spec_pert(1) = 1
         openrsp_cfg_residue_spec_index(1,1) = 3
      end if
      
      ! Shortcut for 3PCD using London orbitals
      if (kw_matches(word, '.RS3PCD')) then
         openrsp_cfg_residue_3pcd = .true.
         openrsp_cfg_residue_n_pert = 4
         openrsp_cfg_residue_order = 1
         openrsp_cfg_specify_kn(1)=0
         openrsp_cfg_specify_kn(2)=3
         allocate(openrsp_cfg_residue_spec_pert(1))
         allocate(openrsp_cfg_residue_spec_index(1,1))
         openrsp_cfg_residue_spec_pert(1) = 1
         openrsp_cfg_residue_spec_index(1,1) = 4
      end if

      if (kw_matches(word, '.RESRUL')) then
         call kw_read(word, openrsp_cfg_specify_kn(1))
         call kw_read(word, openrsp_cfg_specify_kn(2))
      end if


      ! Read the number of states for which residues are calcd.
      if (kw_matches(word, '.RSNSTA')) then
         call kw_read(word, openrsp_cfg_residue_nstates)
      end if

      ! Read for which states residues are calcd.
      if (kw_matches(word, '.RESSTA')) then
         allocate(openrsp_cfg_residue_states(openrsp_cfg_residue_order,openrsp_cfg_residue_nstates))
         openrsp_cfg_residue_states=-999999
         do i = 1, openrsp_cfg_residue_order 
           call kw_read(word, openrsp_cfg_residue_states(i,:))
           ! Additional modification for double residue of one state, not yet needed
           if( openrsp_cfg_residue_order.gt.0..and.(openrsp_cfg_residue_order.gt.openrsp_cfg_residue_nstates )) then
               openrsp_cfg_residue_states(i,openrsp_cfg_residue_nstates)=openrsp_cfg_residue_states(i,1)
           end if
         end do
      end if

      ! Read the corresponding frequencies
      ! Initialized, postprocessing must follow
      if (kw_matches(word, '.RESFRQ')) then
         allocate(openrsp_cfg_residue_freq(openrsp_cfg_residue_n_pert))
         openrsp_cfg_residue_freq = -9999999.9d9
          do j = 1, openrsp_cfg_residue_n_pert !- 1
            call kw_read(word, openrsp_cfg_residue_freq(j))
          end do
      end if

      ! Read the corresponding perturbation labels
      ! Initialized, postprocessing must follow
      if (kw_matches(word, '.SPPLBR')) then
         allocate(openrsp_cfg_residue_plab(openrsp_cfg_residue_n_pert))
         openrsp_cfg_residue_plab = 'XXXX'
           do j = 1, openrsp_cfg_residue_n_pert
            call kw_read(word, openrsp_cfg_residue_plab(j))
           end do
      end if

!     --------------------------------------------------------------------------

      ! OrL: Keyword for Coriolis coupling

      if (kw_matches(word, '.CORIOL')) then
         openrsp_cfg_general_coriolis = .true.
      endif
      
!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------

	! OrL: Keyword for cubic and quartic force fields
	if (kw_matches(word, '.CUBFF ')) then
		openrsp_cfg_general_cubic_force = .true.
	endif
		
	if (kw_matches(word, '.QUARFF')) then
		openrsp_cfg_general_quartic_force = .true.
	endif

!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------
	
	! OrL: Keyword for dipole geometrical derivatives
	if (kw_matches(word, '.DIPGRA ')) then
		openrsp_cfg_general_dipole_gradient = .true.
	end if
	
	if (kw_matches(word, '.DIPHES')) then
		openrsp_cfg_general_dipole_hessian = .true.
	end if
	
	if (kw_matches(word, '.DIPCUB')) then
		openrsp_cfg_general_dipole_cubic = .true.
	end if
	
!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------

	! OrL: Keyword for polarizability geometrical derivatives
	if (kw_matches(word, '.POLGRA')) then
		openrsp_cfg_general_polarizability_gradient = .true.
	end if
	
	if (kw_matches(word, '.POLHES')) then
		openrsp_cfg_general_polarizability_hessian = .true.
	end if
	
	if (kw_matches(word, '.POLCUB')) then
		openrsp_cfg_general_polarizability_cubic = .true.
	end if
	
!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------

	! OrL: Keyword for first hyper-polarizability geometrical derivatives
	if (kw_matches(word, '.FHPGRA')) then
		openrsp_cfg_general_hyper_polarizability_gradient = .true.
	end if
	
	if (kw_matches(word, '.FHPHES')) then
		openrsp_cfg_general_hyper_polarizability_hessian = .true.
	end if
	
	if (kw_matches(word, '.FHPCUB')) then
		openrsp_cfg_general_hyper_polarizability_cubic = .true.
	end if
	
!     --------------------------------------------------------------------------

!     --------------------------------------------------------------------------

 	if (kw_matches(word, '.FHPFRQ')) then
         call kw_read(word, openrsp_cfg_general_hyper_polarizability_geom_freq)
    end if
    
    if (kw_matches(word, '.POLFRQ')) then
         call kw_read(word, openrsp_cfg_general_polarizability_geom_freq)
    end if
    
!     --------------------------------------------------------------------------

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

end module
