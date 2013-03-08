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

   real(8) :: openrsp_cfg_speed_of_light = 137.0359998d0
   logical :: openrsp_cfg_skip_llss      = .false.
   logical :: openrsp_cfg_skip_ssss      = .true.

   integer              :: openrsp_cfg_nr_real_freqs = 1
   real(8), allocatable :: openrsp_cfg_real_freqs(:)
   integer              :: openrsp_cfg_nr_imag_freqs = 1
   real(8), allocatable :: openrsp_cfg_imag_freqs(:)

   logical :: openrsp_cfg_magnetizability = .false.

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

   ! MaR: Keywords for SHG (and possibly higher-orders later) calculation
   logical :: openrsp_cfg_general_shg = .false.

   ! MaR: Keywords for EFISHG-CID calculation
   logical :: openrsp_cfg_general_efishgcid = .false.
   
   ! MaR: Keywords for SFG calculation
   logical :: openrsp_cfg_general_sfg = .false.

   ! MaR: Keywords for hyper-Raman calculation
   logical :: openrsp_cfg_general_hyper_raman = .false.
   real(8) :: openrsp_cfg_general_hypram_freqscale = 1.0d0



end module
