module openrsp_cfg

   implicit none

   integer :: openrsp_cfg_print_level

   ! maximum number of micro iterations in the iterative solution of
   ! the frequency independent linear response functions
   integer :: openrsp_cfg_solver_maxitr = 100
   ! maximum dimension of the sub-block of the configuration Hessian
   integer :: openrsp_cfg_solver_maxphp = 0
   ! maximum dimension of the reduced space to which new basis vectors are added
   integer :: openrsp_cfg_solver_maxred = 400
   ! convergence threshold for the solution of the frequency-independent response equations
   real(8) :: openrsp_cfg_solver_thresh = 1.0d-10
   ! true for optimal orbital trial vectors in the iterative solution of
   ! the frequency-dependent linear response equations
   logical :: openrsp_cfg_solver_optorb = .false.

   ! broadening (damping) parameter for the complex polarization propagator solver
   real(8) :: openrsp_cfg_cpp_damping = 0.0
   ! if using the complex polarization propagator solver
   logical :: openrsp_cfg_cpp_used = .false.

   real(8) :: openrsp_cfg_speed_of_light = 137.0359998d0
   logical :: openrsp_cfg_skip_llss      = .false.
   logical :: openrsp_cfg_skip_ssss      = .true.

   integer              :: openrsp_cfg_nr_real_freqs = 1
   real(8), allocatable :: openrsp_cfg_real_freqs(:)
   integer              :: openrsp_cfg_nr_imag_freqs = 1
   real(8), allocatable :: openrsp_cfg_imag_freqs(:)

   logical :: openrsp_cfg_magnetizability = .false.
   logical :: openrsp_cfg_vcd = .false.

   ! MaR: Temperature (in Kelvin) for purposes where it might apply (default: 298 K)
   real(8) :: openrsp_cfg_temperature = 298.0d0

   ! MaR: Support for multiple frequency tuples
   integer              :: openrsp_cfg_nr_freq_tuples = 1

   ! MaR: Configuration parameters for custom property specification
   logical :: openrsp_cfg_general_specify = .false.
   integer :: openrsp_cfg_specify_order = 0
   integer, dimension(2) :: openrsp_cfg_specify_kn
   character(4), allocatable, dimension(:) :: openrsp_cfg_specify_plab

   ! MaR: Keywords for pure vibrational contribution calculation
   logical :: openrsp_cfg_general_pv2f = .false.
   logical :: openrsp_cfg_general_pv3f = .false.
   logical :: openrsp_cfg_general_pv4f = .false.
   ! Order of electrical, mechanical, or total anharmonicity (default for all: 0)
   integer :: openrsp_cfg_general_pv_el_anh = 0
   integer :: openrsp_cfg_general_pv_mech_anh = 0
   integer :: openrsp_cfg_general_pv_total_anh = 0

! MaR: Keywords for ZPVA contribution calculations
   logical :: openrsp_cfg_general_zpva2f = .true.
   logical :: openrsp_cfg_general_zpva3f = .true.
   logical :: openrsp_cfg_general_zpva4f = .true.

   ! MaR: Keywords for SHG (and possibly higher-orders later) calculation
   logical :: openrsp_cfg_general_shg = .false.

   ! MaR: Keywords for EFISHG-CID calculation
   logical :: openrsp_cfg_general_efishgcid = .false.
   
   ! MaR: Keywords for SFG calculation
   logical :: openrsp_cfg_general_sfg = .false.

   ! MaR: Keywords for hyper-Raman calculation
   logical :: openrsp_cfg_general_hyper_raman = .false.
   real(8) :: openrsp_cfg_general_hypram_freqscale = 1.0d0

   ! MaR: Keywords for Cartesian to normal mode transformation routines
   logical :: openrsp_cfg_general_trans_cartnc = .false.
   integer :: openrsp_cfg_general_cartnc_order_field = 0
   integer :: openrsp_cfg_general_cartnc_order_geo = 0
   
   ! MaR: Keyword to suppress calculation of Cartesian basis tensor in conversion routines
   ! (just perform transformation on tensor which was previously calculated and stored in file)
   
   logical :: openrsp_cfg_general_suppress_calc = .false.
   
	! OrL: Keyword for Coriolis coupling
	logical :: openrsp_cfg_general_coriolis = .false.
	
	! OrL: Keyword for cubic and quartic force fields
	logical :: openrsp_cfg_general_cubic_force = .false.
	logical :: openrsp_cfg_general_quartic_force = .false.
	
	! OrL: Keyword for dipole geometrical derivatives
	logical :: openrsp_cfg_general_dipole_gradient = .false.
	logical :: openrsp_cfg_general_dipole_hessian = .false.
	logical :: openrsp_cfg_general_dipole_cubic = .false.
	
	! OrL: Keyword for polarizability geometrical derivatives
	logical :: openrsp_cfg_general_polarizability_gradient = .false.
	logical :: openrsp_cfg_general_polarizability_hessian = .false.
	logical :: openrsp_cfg_general_polarizability_cubic = .false.
	
	! OrL: Keyword for first hyper-polarizability geometrical derivatives
	logical :: openrsp_cfg_general_hyper_polarizability_gradient = .false.
	logical :: openrsp_cfg_general_hyper_polarizability_hessian = .false.
	logical :: openrsp_cfg_general_hyper_polarizability_cubic = .false.
	
	real(8) :: openrsp_cfg_general_polarizability_geom_freq       = 0.0d0
	real(8) :: openrsp_cfg_general_hyper_polarizability_geom_freq = 0.0d0

   ! keywords that control the XCint grid (not the Dalton grid)
   logical :: openrsp_cfg_use_xcint_grid = .false.
   real(8) :: openrsp_cfg_radint = 1.0d-12
   integer :: openrsp_cfg_angmin = 86
   integer :: openrsp_cfg_angint = 302

end module
