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

   ! MaR: Support for multiple frequency tuples
   integer              :: openrsp_cfg_nr_freq_tuples = 1

   logical :: openrsp_cfg_gradient        = .false.
   logical :: openrsp_cfg_hessian         = .false.
   logical :: openrsp_cfg_dipole_gradient = .false.
   logical :: openrsp_cfg_dipole_hessian  = .false.
   logical :: openrsp_cfg_cubic_ff        = .false.
   logical :: openrsp_cfg_quartic_ff      = .false.

   logical :: openrsp_cfg_pnc_gradient   = .false.
   logical :: openrsp_cfg_pnc_hessian    = .false.
   integer :: openrsp_cfg_pnc_center     = 1

   logical :: openrsp_cfg_general_f      = .false.
   logical :: openrsp_cfg_general_g      = .false.
   logical :: openrsp_cfg_general_ff     = .false.
   logical :: openrsp_cfg_general_gf     = .false.
   logical :: openrsp_cfg_general_gg     = .false.
   logical :: openrsp_cfg_general_fff    = .false.
   logical :: openrsp_cfg_general_gff    = .false.
   logical :: openrsp_cfg_general_ggf    = .false.
   logical :: openrsp_cfg_general_ggg    = .false.
   logical :: openrsp_cfg_general_ffff   = .false.
   logical :: openrsp_cfg_general_gfff   = .false.
   logical :: openrsp_cfg_general_ggff   = .false.
   logical :: openrsp_cfg_general_gggf   = .false.
   logical :: openrsp_cfg_general_gggg   = .false.
   logical :: openrsp_cfg_general_fffff  = .false.
   logical :: openrsp_cfg_general_gffff  = .false.
   logical :: openrsp_cfg_general_ggfff  = .false.
   logical :: openrsp_cfg_general_gggff  = .false.
   logical :: openrsp_cfg_general_ggggf  = .false.
   logical :: openrsp_cfg_general_ggggg  = .false.
   logical :: openrsp_cfg_general_ffffff = .false.
   logical :: openrsp_cfg_general_gfffff = .false.
   logical :: openrsp_cfg_general_ggffff = .false.
   logical :: openrsp_cfg_general_gggfff = .false.
   logical :: openrsp_cfg_general_ggggff = .false.
   logical :: openrsp_cfg_general_gggggf = .false.
   logical :: openrsp_cfg_general_gggggg = .false.

   ! MaR: Configuration parameters for custom property specification

   logical :: openrsp_cfg_general_specify = .false.
   integer :: openrsp_cfg_specify_order = 0
   integer, dimension(2) :: openrsp_cfg_specify_kn
   character(4), allocatable, dimension(:) :: openrsp_cfg_specify_plab

   ! MaR: Keywords for pure vibrational contribution calculation

   logical :: openrsp_cfg_general_pv2f = .false.
   logical :: openrsp_cfg_general_pv3f = .false.
   logical :: openrsp_cfg_general_pv4f = .false.


end module
