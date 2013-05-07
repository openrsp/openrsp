module interface_1el

   use openrsp_cfg
   use matrix_defop
   use matrix_lowlevel, only: mat_init, mat_zero_like
   use interface_molecule
   use interface_basis
   use interface_f77_memory
   use interface_pcm
   use interface_io
   use interface_1el_dirac
   use rsp_indices_and_addressing, only: get_triang_blks_offset

   implicit none

   public interface_1el_init
   public interface_1el_ovlave_tr
   public interface_1el_ovlave_half_diff
   public interface_1el_oneave_tr
   public interface_1el_ovlint_tr
   public interface_1el_ovlint_half_diff
   public interface_1el_oneint_tr

   public oneint_ave
   public get1in_ave_ifc
   public onedrv_ave_ifc
   public legacy_read_integrals

   private

contains

   subroutine interface_1el_init()
   end subroutine

#ifdef VAR_LSDALTON
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ave, DFD)
      ! Gen1Int interface
      use gen1int_host
      !> number of fields
      integer,       intent(in)  :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> energy-weighted density matrix
      type(matrix),  intent(in) :: DFD
      !> output average
      complex(8),    intent(inout) :: ave(propsize)
      STOP 'interface_1el_ovlave_tr not implemented for VAR_LSDALTON'
   end subroutine
#else
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ave, DFD)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      !> number of fields
      integer,       intent(in)  :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      integer :: num_derv, j, k, ave_offset
      integer, allocatable, dimension(:) :: order_derv, tmp_index
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> energy-weighted density matrix
      type(matrix),  intent(in) :: DFD
      !> output average
      complex(8),    intent(inout) :: ave(propsize)
      !----------------------------------------------
      type(matrix) A(2)
      real(8), parameter   :: fdistep = 2d0**(-25)
      real(8), allocatable :: tmp(:,:,:)
      real(8), allocatable :: fdi(:,:)
      integer i
      integer order_geo                      !order of total geometric derivatives
      integer order_mag                      !order of total magnetic derivatives
      integer num_atom                       !number of atoms
      integer num_coord                      !number of atomic coordinates
      integer num_geom                       !number of total geometric derivatives
      integer num_expt                       !number of all expectation values
      integer :: num_addr, num_order
      integer, allocatable :: address_list(:,:)
      real(8), allocatable :: val_expt(:,:)  !expectation values, real numbers
      integer ierr                           !error information
      logical :: all_frequencies_zero


      if (any(f == 'EL  ')) then
         ave = 0.0
      else if (any(f == 'ELGR')) then
         ave = 0.0


      else

      ! CURRENTLY ONLY SUPPORTS ELECTRIC FIELD PERTURBATION IN ADDITION TO GEOMETRIC

      if (count(f=='EL  ') > 0) then

         ! THIS CASE SHOULD NOT BE ENCOUNTERED IN OVLAVE

         num_derv = 1

         allocate(order_derv(num_derv))
         order_derv = (/count(f=='EL  ')/)

      else

         num_derv = 0
         allocate(order_derv(1))
         order_derv = (/0/)

      end if

         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
         ! gets the order of total geometric derivatives
         order_mag = count(f=='MAG ')

!          if (order_geo /= nf) then
!             call quit("interface_1el_ovlave>> only geometric derivatives implemented!")
!          end if

         ! sets the number of total geometric derivatives
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! allocates memory for expectation values
         num_expt = num_geom
         allocate(val_expt(propsize, 1), stat=ierr)
         if (ierr /= 0) call quit("interface_1el_ovlave>> failed to allocate val_expt!")
         val_expt = 0.0

! MaR: Removing a currently wrong contribution to test other parts of code
if (order_mag > 0) then
 
        ! MR: PARAMETERS FOR ADAPTED gen1int CALL BELOW ARE MAYBE SPECIFIED 
        ! INCORRECTLY (ALSO LAST ARG)
            ! calculates the expectaion values of overlap matrix
            call gen1int_host_get_expt(LONDON, INT_OVERLAP,       &
                                       0,                          &  !multipole moments
                                       0,                          &
                                       0, 0, order_mag,            &  !magnetic derivatives
                                       0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                       0, 0, 0, (/0/),             &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),             &  !geometric derivatives on ket center
                                       order_geo,                  &  !total geometric derivatives
                                       order_geo,                  &
                                       0, (/0/),                   &
                                       .false., .false., .false.,  &  !not implemented yet
                                       1, (/DFD/), propsize,       &  !expectation values
                                       val_expt, .false.,          &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
#else
                                       1, (/1, 1/),                &
#endif
                                       get_print_unit(), 0)
            val_expt = -val_expt
            ! assigns the output average

else

        ! MR: PARAMETERS FOR ADAPTED gen1int CALL BELOW ARE MAYBE SPECIFIED 
        ! INCORRECTLY (ALSO LAST ARG)
            ! calculates the expectaion values of overlap matrix
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,       &
                                       0,                          &  !multipole moments
                                       0,                          &
                                       0, 0, order_mag,            &  !magnetic derivatives
                                       0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                       0, 0, 0, (/0/),             &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),             &  !geometric derivatives on ket center
                                       order_geo,                  &  !total geometric derivatives
                                       order_geo,                  &
                                       0, (/0/),                   &
                                       .false., .false., .false.,  &  !not implemented yet
                                       1, (/DFD/), propsize,       &  !expectation values
                                       val_expt, .false.,          &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
#else
                                       1, (/1, 1/),                &
#endif
                                       get_print_unit(), 0)
            val_expt = -val_expt
            ! assigns the output average

end if

! MaR: Quick fix to get addressing in order
if (order_mag > 0) then

! write(*,*) 'val expt', val_expt(:,1)

ave(1) = val_expt(1,1)
ave(2) = val_expt(2,1)
ave(3) = val_expt(3,1)

else
            ! MR: ASSIGN DATA TO ave

            num_order = sum(order_derv)+order_geo
            ! gets the number of unique total geometric derivatives
            if (order_geo/=0) then
               call geom_total_num_derv(order_geo, min(num_atom, order_geo), num_atom, num_addr)
            else
               num_addr = 1
            end if
            ! gets the number of derivatives/powers/moments
            do ierr = 1, num_derv
               num_addr = num_addr*(order_derv(ierr)+1)*(order_derv(ierr)+2)/2
            end do
            allocate(address_list(num_order,num_addr), stat=ierr)
            if (ierr/=0) stop "failed to allocate address_list!"
            ! gets the list of addresses
            call get_address_list(num_derv, order_derv,                     &
                                  num_atom, order_geo, min(num_atom, order_geo), &
                                  num_order, num_addr, address_list)
        
            allocate(tmp_index(nf))
        
            do i = 1, num_addr
        
               tmp_index(:) = address_list(:, i)
               ave_offset = get_triang_blks_offset(nblks, nf, blk_info, blk_sizes, tmp_index)
               ave(ave_offset) = val_expt(i,1)
        
            end do
        
            deallocate(address_list)
            deallocate(tmp_index)
end if     
         deallocate(val_expt)

      end if

   end subroutine
#endif




#ifdef VAR_LSDALTON
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! MaR: This routine not finished - awaiting development in integral code
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlave_half_diff(nbra, fbra, cbra, ncbra, &
              nket, fket, cket, ncket, nblks_tuple, blks_tuple_info, blk_sizes, & 
              D, propsize, ave)


      use gen1int_host
      !> number of fields
      integer,       intent(in)    :: nbra, nket, propsize
      integer, dimension(2) :: nblks_tuple
      integer, dimension(2, nbra + nket, 3) :: blks_tuple_info
      integer, dimension(2, nbra + nket) :: blk_sizes
      !> field labels in std order
      character(4),  intent(in)    :: fbra(nf), fket(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: cbra(nf), ncbra(nf), cket(nf), ncket(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      complex(8),  intent(inout):: ave(propsize)
      real(8) :: tmp_ave(propsize)
      type(matrix) :: D
      STOP 'interface_1el_ovlave_half_diff not implemented for LSDALTON. TK'
   end subroutine
#else
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! MaR: This routine not finished - awaiting development in integral code
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlave_half_diff(nbra, fbra, cbra, ncbra, &
              nket, fket, cket, ncket, nblks_tuple, blks_tuple_info, blk_sizes, & 
              D, propsize, ave)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      !> number of fields
      integer,       intent(in)    :: nbra, nket, propsize
      integer, dimension(2) :: nblks_tuple
      integer, dimension(2, nbra + nket, 3) :: blks_tuple_info
      integer, dimension(2, nbra + nket) :: blk_sizes
      !> field labels in std order
      character(4),  intent(in)    :: fbra(nbra), fket(nket)
      !> first and number of- components in each field
      integer,       intent(in)    :: cbra(nbra), ncbra(nbra), cket(nket), ncket(nket)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      complex(8),  intent(inout):: ave(propsize)
      real(8) :: tmp_ave(propsize)
      type(matrix) :: D
      !------------------------------------------------
      integer      i


           ! MaR: Only one magnetic field supported for now
           tmp_ave = 0.0

if ((count(fbra=='MAG ') > 0) .or. (count(fket=='MAG ') > 0) ) then


           ! MaR: THIS IS THE INT CALL STRUCTURE - ADAPT TO AVE AND REMEMBER D
           ! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALL BELOW
           call gen1int_host_get_expt(LONDON, INT_OVERLAP,       &
                                     0,                          &  !multipole moments
                                     0,                          &
                                     count(fbra=='MAG '),        &
                                     count(fket=='MAG '), 0,     &  !magnetic derivatives
                                     0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                     count(fbra=='GEO '),        &  !geometric derivatives on bra center
                                     count(fbra=='GEO '),        &
                                     0, (/0/),                   &
                                     count(fket=='GEO '),        &  !geometric derivatives on ket center
                                     count(fket=='GEO '),        &
                                     0, (/0/),                   &  
                                     0,                          &  !total geometric derivatives
                                     0,                          &
                                     0, (/0/),                   &
                                     .false., .false., .false.,  &  !not implemented yet
                                     1, (/D/),                   &
                                     propsize, tmp_ave, .false.,   &  !integral matrices
#ifdef PRG_DIRAC
                                     2, (/1, 1, 2, 2/),          &
#else
                                     1, (/1, 1/),                &
#endif
                                     get_print_unit(), 0)

          ! MaR: OFFSETS AND LOCATIONS IN MEMORY POSTPONED
          ! MAKE THE CORRESPONDING OVLAVE TYPE ROUTINE

          ! Just magnetic to think about for now
          ave = tmp_ave

! if ((count(fbra=='MAG ') == 1) .or. (count(fket=='MAG ') == 1)) then
! 
! write(*,*) 'fbra', fbra
! 
! write(*,*) tmp_ave
! 
! end if

end if
   end subroutine
#endif




#ifdef VAR_LSDALTON
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average 1-electron integrals perturbed by fields f
   !>        with the (perturbed) density matrix D
   subroutine interface_1el_oneave_tr(nf, f, c, nc, D, blk_info, & 
                                      blk_sizes, propsize, ave)
      ! Gen1Int interface
      use gen1int_host
      !> number of fields
      integer,       intent(in)  :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> density matrix to average over
      type(matrix),  intent(in)  :: D
      !> output average
      complex(8),    intent(out) :: ave(propsize)
      STOP 'interface_1el_oneave_tr not implemented for LSDALTON'
   end subroutine
#else
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average 1-electron integrals perturbed by fields f
   !>        with the (perturbed) density matrix D
   subroutine interface_1el_oneave_tr(nf, f, c, nc, D, nblks, blk_info, & 
                                      blk_sizes, propsize, ave)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      !> number of fields
      integer,       intent(in)  :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      integer :: num_derv, j, k, ave_offset
      integer, allocatable, dimension(:) :: order_derv, tmp_index
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> density matrix to average over
      type(matrix),  intent(in)  :: D
      !> output average
      complex(8),    intent(out) :: ave(propsize)
      !----------------------------------------------
      integer order_mom, order_elgr          !order of Cartesian multipole moments
      integer num_mom                        !number of Cartesian multipole moments
      integer order_geo                      !order of total geometric derivatives
      integer order_mag                      !order of total magnetic derivatives
      integer num_atom                       !number of atoms
      integer num_coord                      !number of atomic coordinates
      integer num_geom                       !number of total geometric derivatives
      integer num_expt                       !number of all expectation values
      integer :: num_addr, num_order
      integer, allocatable :: address_list(:,:)
      real(8), allocatable :: val_expt(:, :) !expectation values, real numbers
      real(8), allocatable :: temp(:)
      integer ierr                           !error information
      type(matrix) :: T
      integer :: ixyz, i

      ! gets the order of Cartesian multipole moments
      order_mom = count(f=='EL  ')
      order_elgr = count(f=='ELGR')
      order_mag = count(f=='MAG ')

      ! CURRENTLY ONLY SUPPORTS ELECTRIC FIELD PERTURBATION IN ADDITION TO GEOMETRIC

      if (count(f=='EL  ') > 0) then

         num_derv = 1

         allocate(order_derv(num_derv))
         order_derv = (/count(f=='EL  ')/)

      else

         num_derv = 0
         allocate(order_derv(1))
         order_derv = (/0/)

      end if

      if (order_mom > 1) then
         ave = 0.0
      else if ((order_elgr > 0) .AND. (order_mom > 0)) then
         ave = 0.0
      else if ((order_mag > 0) .AND. (order_mom > 0)) then
         ave = 0.0
      else if ((order_elgr > 1)) then
         ave = 0.0
      else

         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
!          if (order_mom+order_geo/=nf) then
!             call quit("interface_1el_oneave>> only electric and geometric perturbations implemented!")
!          end if

         ! sets the number of operators and derivatives
         num_mom = (order_mom+1)*(order_mom+2)/2
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! allocates memory for expectation values
         num_expt = num_mom*num_geom
         allocate(val_expt(propsize,1), stat=ierr)
         if (ierr/=0) call quit("interface_1el_oneave>> failed to allocate val_expt!")
         val_expt = 0.0

         ! MR: PARAMETERS FOR ADAPTED gen1int CALLS BELOW ARE MAYBE SPECIFIED 
         ! INCORRECTLY (ALSO LAST ARG)
         ! electric perturbations
         if (order_mom/=0 .or. order_elgr/=0) then

            call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                       order_mom + 2*order_elgr,    &  !multipole moments
                                       0,                           &
                                       0, 0, order_mag,             &  !magnetic derivatives
                                       0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                       0, 0, 0, (/0/),              &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),              &  !geometric derivatives on ket center
                                       order_geo,                   &  !total geometric derivatives
                                       order_geo,                   &
                                       0, (/0/),                    &
                                       .false., .false., .false.,   &  !not implemented yet
                                       1, (/D/), propsize,          &  !expectation values
                                       val_expt, .false.,           &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),           &
#else
                                       1, (/1, 1/),                &
#endif
                                       get_print_unit(), 0)
         else if (order_mag > 0) then

            call gen1int_host_get_expt(LONDON, INT_CART_MULTIPOLE, &
                                       0,    &  !multipole moments
                                       0,                           &
                                       0, 0, order_mag,             &  !magnetic derivatives
                                       0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                       0, 0, 0, (/0/),              &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),              &  !geometric derivatives on ket center
                                       order_geo,                   &  !total geometric derivatives
                                       order_geo,                   &
                                       0, (/0/),                    &
                                       .false., .false., .false.,   &  !not implemented yet
                                       1, (/D/), propsize,          &  !expectation values
                                       val_expt, .false.,           &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),           &
#else
                                       1, (/1, 1/),                &
#endif
                                       get_print_unit(), 0)


         ! only geometric perturbations
         else

#ifdef PRG_DALTON
            call gen1int_host_get_expt(NON_LAO, INT_ONE_HAMIL,    &
                                       0,                         &  !multipole moments
                                       0,                         &
                                       0, 0, 0,                   &  !magnetic derivatives
                                       0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                       0, 0, 0, (/0/),            &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),            &  !geometric derivatives on ket center
                                       order_geo,                 &  !total geometric derivatives
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       .false., .false., .false., &  !not implemented yet
                                       1, (/D/), propsize,        &  !expectation values
                                       val_expt, .false.,         &
                                       1, (/1, 1/),               &
                                       get_print_unit(), 0)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
            allocate(temp(size(val_expt)*3))

            ! nuclear attraction
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_POT_ENERGY,   &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0, 0, (/0/),            &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),            &  !geometric derivatives on ket center
                                       order_geo,                 &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       .false., .false., .false., &
                                       1, (/D/), propsize,        &
                                       temp, .false.,             &
                                       2, (/1, 1, 2, 2/),         &
                                       get_print_unit(), 0)
            do i = 1, propsize
               val_expt(i, 1) = val_expt(i, 1) + temp(i)
            end do

            ! beta' matrix
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,      &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0, 0, (/0/),            &  !geometric derivatives on bra center
                                       0, 0, 0, (/0/),            &  !geometric derivatives on ket center
                                       order_geo,                 &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       .false., .false., .false., &
                                       1, (/D/), propsize,        &
                                       temp, .false.,             &
                                       1, (/2, 2/),               &
                                       get_print_unit(), 0)
            do i = 1, propsize
               val_expt(i, 1) = val_expt(i, 1) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*temp(i)
            end do

            ! kinetic energy
            T = 0*D
            call mat_ensure_alloc(T, only_alloc=.true.)
            do ixyz = 1, 3
               T%elms = 0.0d0
               call dcopy(D%nrow*D%ncol, D%elms(1, 1, 5-ixyz), 1, T%elms, 1)
               temp = 0.0d0
               call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                          0,                           &
                                          1,                           &
                                          0, 0, 0,                     &
                                          0, 0, 0,                     &
                                          0, 0, 0, (/0/),              &  !geometric derivatives on bra center
                                          0, 0, 0, (/0/),              &  !geometric derivatives on ket center
                                          order_geo,                   &
                                          order_geo,                   &
                                          0, (/0/),                    &
                                          .false., .false., .false.,   &
                                          1, (/T/), propsize*3,        &
                                          temp, .false.,               &
                                          2, (/1, 2, 2, 1/),           &
                                          get_print_unit(), 0)
               do i = 1, propsize
                  val_expt(i, 1) = val_expt(i, 1) + openrsp_cfg_speed_of_light*temp((i-1)*3 + ixyz)
               end do
            end do
            T = 0

            deallocate(temp)
#endif /* ifdef PRG_DIRAC */

         end if

! MaR: Quick fix to get addressing in order
if (order_mag > 0) then

ave(1) = val_expt(1,1)
ave(2) = val_expt(2,1)
ave(3) = val_expt(3,1)


elseif (order_elgr > 0) then

ave(1) = val_expt(1,1)
ave(2) = val_expt(2,1)
ave(3) = val_expt(4,1)
ave(4) = val_expt(3,1)
ave(5) = val_expt(5,1)
ave(6) = val_expt(6,1)

else

         num_order = sum(order_derv)+order_geo
         ! gets the number of unique total geometric derivatives
         if (order_geo/=0) then
            call geom_total_num_derv(order_geo, min(num_atom, order_geo), num_atom, num_addr)
         else
            num_addr = 1
         end if
         ! gets the number of derivatives/powers/moments
         do ierr = 1, num_derv
           num_addr = num_addr*(order_derv(ierr)+1)*(order_derv(ierr)+2)/2
         end do
    
         allocate(address_list(num_order,num_addr), stat=ierr)
         if (ierr/=0) stop "failed to allocate address_list!"
         ! gets the list of addresses
         call get_address_list(num_derv, order_derv,                     &
                               num_atom, order_geo, min(num_atom, order_geo), &
                               num_order, num_addr, address_list)

         allocate(tmp_index(nf))

         do i = 1, num_addr

            k = 1

            do j = order_mom + 1, order_mom + order_geo

               tmp_index(k) = address_list(j, i)
               k = k + 1

            end do

            do j = 1, order_mom

               tmp_index(k) = address_list(j, i)
               k = k + 1

            end do

            ave_offset = get_triang_blks_offset(nblks, nf, blk_info, blk_sizes, tmp_index)
            ave(ave_offset) = val_expt(i,1)

         end do

         deallocate(address_list)
         deallocate(tmp_index)

end if

         deallocate(val_expt)

      end if
   end subroutine
#endif





#ifdef VAR_LSDALTON
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ovl, w, fock)
      use gen1int_host
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> resulting overlap integral matrices (incoming content deleted)
      type(matrix),  intent(inout), optional :: ovl(propsize)
      !> frequencies of each field
      complex(8),    intent(in),    optional :: w(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout), optional :: fock(propsize)
      STOP 'interface_1el_ovlint_tr not implemented for LSDALTON. TK'
   end subroutine
#else
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ovl)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      integer :: num_derv, j, k, ave_offset, ierr
      integer, allocatable, dimension(:) :: order_derv, tmp_index
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> resulting overlap integral matrices (incoming content deleted)
      type(matrix),  intent(inout) :: ovl(propsize)
      type(matrix)   :: ovl_tmp(propsize)
      !------------------------------------------------
      integer      i
      integer order_geo  !order of total geometric derivatives
      integer num_atom   !number of atoms
      integer num_coord  !number of atomic coordinates
      integer num_geom   !number of total geometric derivatives
      integer num_ints   !number of integral matrices
      integer :: num_addr, num_order
      integer, allocatable :: address_list(:,:)
      logical :: all_frequencies_zero

      ! CURRENTLY ONLY SUPPORTS ELECTRIC FIELD PERTURBATION IN ADDITION TO GEOMETRIC

      if (count(f=='EL  ') > 0) then

         ! THIS CASE SHOULD NOT BE ENCOUNTERED IN OVLINT

         num_derv = 1

         allocate(order_derv(num_derv))
         order_derv = (/count(f=='EL  ')/)

      else

         num_derv = 0
         allocate(order_derv(1))
         order_derv = (/0/)

      end if



      if (any(f=='EL  ')) then

      else

         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
         if (order_geo/=nf) then
           call quit("interface_1el_ovlint>> only geometric derivatives implemented!")
         end if

         ! sets the number of total geometric derivatives
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! sets the number of integral matrices
         num_ints = num_geom

         ! allocates matrices
         do i = 1, propsize
            if (.not.isdef(ovl(i))) then
               call mat_init(ovl(i), nr_ao, nr_ao)
            end if
            if (.not.isdef(ovl_tmp(i))) then
               call mat_init(ovl_tmp(i), nr_ao, nr_ao)
            end if
         end do
         ! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALL BELOW
         ! calculates the overlap matrix
         call gen1int_host_get_int(NON_LAO, INT_OVERLAP,       &
                                   0,                          &  !multipole moments
                                   0,                          &
                                   0, 0, 0,                    &  !magnetic derivatives
                                   0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                   0, 0, 0, (/0/),             &  !geometric derivatives on bra center
                                   0, 0, 0, (/0/),             &  !geometric derivatives on ket center
                                   order_geo,                  &  !total geometric derivatives
                                   order_geo,                  &
                                   0, (/0/),                   &
                                   .false., .false., .false.,  &  !not implemented yet
                                   num_ints, ovl_tmp, .false.,     &  !integral matrices
#ifdef PRG_DIRAC
                                   2, (/1, 1, 2, 2/),          &
#else
                                   1, (/1, 1/),                &
#endif
                                   get_print_unit(), 0)



         num_order = sum(order_derv)+order_geo
         ! gets the number of unique total geometric derivatives
         if (order_geo/=0) then
            call geom_total_num_derv(order_geo, min(num_atom, order_geo), num_atom, num_addr)
         else
            num_addr = 1
         end if
         ! gets the number of derivatives/powers/moments
         do ierr = 1, num_derv
            num_addr = num_addr*(order_derv(ierr)+1)*(order_derv(ierr)+2)/2
         end do

         allocate(address_list(num_order,num_addr), stat=ierr) 
  
         if (ierr/=0) stop "failed to allocate address_list!"
         ! gets the list of addresses
         call get_address_list(num_derv, order_derv,                     &
                               num_atom, order_geo, min(num_atom, order_geo), &
                               num_order, num_addr, address_list)
 
         allocate(tmp_index(nf))

         do i = 1, num_addr
 
            tmp_index(:) = address_list(:, i)
            ave_offset = get_triang_blks_offset(nblks, nf, blk_info, blk_sizes, tmp_index)
            ovl(ave_offset) = ovl_tmp(i)

         end do

         deallocate(address_list)
         deallocate(tmp_index)

      end if
   end subroutine
#endif







#ifdef VAR_LSDALTON
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! MaR: This routine not finished - awaiting development in integral code
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint_half_diff(nr_ao, nbra, fbra, cbra, ncbra, &
              nket, fket, cket, ncket, nblks_tuple, blks_tuple_info, blk_sizes, & 
              propsize, fock)


      use gen1int_host
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nbra, nket, propsize
      integer, dimension(2) :: nblks_tuple
      integer, dimension(2, nbra + nket, 3) :: blks_tuple_info
      integer, dimension(2, nbra + nket) :: blk_sizes
      !> field labels in std order
      character(4),  intent(in)    :: fbra(nf), fket(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: cbra(nf), ncbra(nf), cket(nf), ncket(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout):: fock(propsize)
      STOP 'interface_1el_ovlint_half_diff not implemented for LSDALTON. TK'
   end subroutine
#else
! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! MaR: This routine not finished - awaiting development in integral code
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint_half_diff(nr_ao, nbra, fbra, cbra, ncbra, &
              nket, fket, cket, ncket, nblks_tuple, blks_tuple_info, blk_sizes, & 
              propsize, fock)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nbra, nket, propsize
      integer, dimension(2) :: nblks_tuple
      integer, dimension(2, nbra + nket, 3) :: blks_tuple_info
      integer, dimension(2, nbra + nket) :: blk_sizes
      !> field labels in std order
      character(4),  intent(in)    :: fbra(nbra), fket(nket)
      !> first and number of- components in each field
      integer,       intent(in)    :: cbra(nbra), ncbra(nbra), cket(nket), ncket(nket)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout):: fock(propsize)
      type(matrix) :: tmp_fock(propsize)
      !------------------------------------------------
      integer      i

           ! MaR: Only one magnetic field supported for now

if ((count(fbra=='MAG ') > 0) .or. (count(fket=='MAG ') > 0) ) then

           do i = 1, propsize
! MaR: It should be safe to remove this
!              if (.not.isdef(fock(i))) then
!                 call mat_init(fock(i), nr_ao, nr_ao)
!              end if
!              if (.not.isdef(tmp_fock(i))) then
!                 call mat_init(tmp_fock(i), nr_ao, nr_ao)
!              end if
                 call mat_zero_like(fock(1), tmp_fock(i))
           end do



           ! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALL BELOW
           call gen1int_host_get_int(LONDON, INT_OVERLAP,        &
                                     0,                          &  !multipole moments
                                     0,                          &
                                     count(fbra=='MAG '),        &
                                     count(fket=='MAG '), 0,     &  !magnetic derivatives
                                     0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                     count(fbra=='GEO '),        &  !geometric derivatives on bra center
                                     count(fbra=='GEO '),        &
                                     0, (/0/),                   &
                                     count(fket=='GEO '),        &  !geometric derivatives on ket center
                                     count(fket=='GEO '),        &
                                     0, (/0/),                   &
                                     0,                          &  !total geometric derivatives
                                     0,                          &
                                     0, (/0/),                   &
                                     .false., .false., .false.,  &  !not implemented yet
                                     propsize, tmp_fock, .false.,&  !integral matrices
#ifdef PRG_DIRAC
                                     2, (/1, 1, 2, 2/),          &
#else
                                     1, (/1, 1/),                &
#endif
                                     get_print_unit(), 0)

          ! MaR: OFFSETS AND LOCATIONS IN MEMORY POSTPONED
          ! MAKE THE CORRESPONDING OVLAVE TYPE ROUTINE

          ! MaR: Just magnetic field for now
          fock = tmp_fock



do i = 1, propsize

! write(*,*) tmp_fock(i)%elms
tmp_fock(i) = 0

end do

end if

   end subroutine
#endif



#ifdef VAR_LSDALTON
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   subroutine interface_1el_oneint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, oneint)
      ! Gen1Int interface
      use gen1int_host
      integer, intent(in)          :: nr_ao, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes
      integer, dimension(nblks, 3) :: blk_info
      integer :: num_derv, j, k, ave_offset
      integer, allocatable, dimension(:) :: order_derv, tmp_index
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> output perturbed integrals
      type(matrix),  intent(inout) :: oneint(propsize)
      STOP 'interface_1el_oneint_tr not implemented for LSDALTON. TK'
   end subroutine
#else
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   subroutine interface_1el_oneint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, oneint)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      integer, intent(in)          :: nr_ao, propsize
      integer :: nblks
      integer, dimension(nblks) :: blk_sizes  ! please document what these block thingies are
      integer, dimension(nblks, 3) :: blk_info
      integer :: num_derv, j, k, ave_offset, ierr
      integer, allocatable, dimension(:) :: order_derv, tmp_index
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> output perturbed integrals
      type(matrix),  intent(inout) :: oneint(propsize)
      type(matrix)   :: oneint_tmp(propsize)
      !--------------------------------------------------
      integer order_mom  !order of Cartesian multipole moments
      integer num_mom    !number of Cartesian multipole moments
      integer order_geo  !order of total geometric derivatives
      integer num_atom   !number of atoms
      integer num_coord  !number of atomic coordinates
      integer num_geom   !number of total geometric derivatives
      integer num_ints   !number of all integral matrices
      integer :: num_addr, num_order
      integer, allocatable :: address_list(:,:)
      integer imat       !incremental recorder over matrices
      integer :: i, ixyz
      type(matrix) :: A
      type(matrix), allocatable :: T(:)

      ! CURRENTLY ONLY SUPPORTS ELECTRIC FIELD PERTURBATION IN ADDITION TO GEOMETRIC

      if (count(f=='EL  ') > 0) then

         num_derv = 1

         allocate(order_derv(num_derv))
         order_derv = (/count(f=='EL  ')/)

      else

         num_derv = 0
         allocate(order_derv(1))
         order_derv = (/0/)

      end if

      if (count(f=='EL  ') > 1) then
         call mat_init(A, nr_ao, nr_ao)
         A = 0*oneint(1)
         call mat_ensure_alloc(A)

         do i = 1, propsize
            if (iszero(oneint(i))) then
               call mat_ensure_alloc(oneint(i))
               oneint(i)%elms = oneint(i)%elms + A%elms
            else
               oneint(i)%elms = oneint(i)%elms + A%elms
            end if
         end do

      else

         ! gets the order of Cartesian multipole moments
         order_mom = count(f=='EL  ')
         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
         if (order_mom+order_geo/=nf) then
           call quit("interface_1el_oneint>> only electric and geometric perturbations implemented!")
         end if

         ! sets the number of operators and derivatives
         num_mom = (order_mom+1)*(order_mom+2)/2
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo
         num_ints = num_mom*num_geom
         ! MaR: DISABLED: NOT A VALID WAY OF DISCERNING THIS ANYMORE
!          if (num_ints/=size(oneint)) then
!             call quit("interface_1el_oneint>> returning specific components is not implemented!")
!          end if

         !FIXME: it is better that we use unique components for higher order
         if (order_mom>1) then
            call quit("interface_1el_oneint>> only the first Cartesian multipole moments implemented!")
         end if

         do imat = 1, size(oneint)
            if (.not.isdef(oneint(imat))) then
               call mat_init(oneint(imat), nr_ao, nr_ao)
            end if
             if (.not.isdef(oneint_tmp(imat))) then
                call mat_init(oneint_tmp(imat), nr_ao, nr_ao)
             end if
         end do

! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALLS BELOW
! MaR: NOTE THAT NOTHING HAS BEEN DONE FOR REORDERING THE EXCLUSIVELY DIRAC CASES BELOW
! MaR: THERE MAY BE NEED FOR REORDERING OUTPUT FROM THOSE CALLS IF USING GENERAL RSP CODE
         ! electric perturbations
         if (order_mom/=0) then
            call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                      order_mom,                   &  !multipole moments
                                      0,                           &
                                      0, 0, 0,                     &  !magnetic derivatives
                                      0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                      0, 0, 0, (/0/),              &  !geometric derivatives on bra center
                                      0, 0, 0, (/0/),              &  !geometric derivatives on ket center
                                      order_geo,                   &  !total geometric derivatives
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      .false., .false., .false.,   &  !not implemented yet
                                      num_ints, oneint_tmp, .false.,   &  !integral matrices
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
#else
                                       1, (/1, 1/),                &
#endif
                                      get_print_unit(), 0)

         ! only geometric perturbations
         else
#ifdef PRG_DIRAC
            allocate(T(3*size(oneint_tmp)))
            do i = 1, 3*size(oneint_tmp)
               T(i) = 0*oneint_tmp(1)
               call mat_ensure_alloc(T(i), only_alloc=.true.)
            end do

            ! nuclear attraction
            call gen1int_host_get_int(NON_LAO, INT_POT_ENERGY,       &
                                      0,                             &
                                      0,                             &
                                      0, 0, 0,                       &
                                      0, 0, 0,                       &
                                      0, 0, 0, (/0/),                &  !geometric derivatives on bra center
                                      0, 0, 0, (/0/),                &  !geometric derivatives on ket center
                                      order_geo,                     &
                                      order_geo,                     &
                                      0, (/0/),                      &
                                      .false., .false., .false.,     &
                                      num_ints, oneint_tmp, .false., &
                                      2, (/1, 1, 2, 2/),             &
                                      get_print_unit(), 0)

            ! beta' matrix
            call gen1int_host_get_int(NON_LAO, INT_OVERLAP,      &
                                      0,                         &
                                      0,                         &
                                      0, 0, 0,                   &
                                      0, 0, 0,                   &
                                      0, 0, 0, (/0/),            &  !geometric derivatives on bra center
                                      0, 0, 0, (/0/),            &  !geometric derivatives on ket center
                                      order_geo,                 &
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      .false., .false., .false., &
                                      num_ints, T, .false.,      &
                                      1, (/2, 2/),               &
                                      get_print_unit(), 0)
            do i = 1, size(oneint_tmp)
               oneint_tmp(i) = oneint_tmp(i) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*T(i)
            end do

            ! kinetic energy
            call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                      0,                           &
                                      1,                           &
                                      0, 0, 0,                     &
                                      0, 0, 0,                     &
                                      0, 0, 0, (/0/),              &  !geometric derivatives on bra center
                                      0, 0, 0, (/0/),              &  !geometric derivatives on ket center
                                      order_geo,                   &
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      .false., .false., .false.,   &
                                      3*num_ints, T, .false.,      &
                                      2, (/1, 2, 2, 1/),           &
                                      get_print_unit(), 0)
            do i = 1, size(oneint_tmp)
               do ixyz = 1, 3
                  call daxpy(T(1)%nrow*T(1)%ncol,              &
                            -openrsp_cfg_speed_of_light,       &
                             T((i-1)*3 + ixyz)%elms,           &
                             1,                                &
                             oneint_tmp(i)%elms(1, 1, 5-ixyz), &
                             1)
               end do
            end do

            do i = 1, 3*size(oneint_tmp)
               T(i) = 0
            end do
            deallocate(T)
#else
            call gen1int_host_get_int(NON_LAO, INT_ONE_HAMIL,    &
                                      0,                         &  !multipole moments
                                      0,                         &
                                      0, 0, 0,                   &  !magnetic derivatives
                                      0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                      0, 0, 0, (/0/),            &  !geometric derivatives on bra center
                                      0, 0, 0, (/0/),            &  !geometric derivatives on ket center
                                      order_geo,                 &  !total geometric derivatives
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      .false., .false., .false., &  !not implemented yet
                                      num_ints, oneint_tmp,      &  !integral matrices
                                      .false., 1, (/1, 1/),      &
                                      get_print_unit(), 0)
#endif
         end if

         num_order = sum(order_derv)+order_geo
         ! gets the number of unique total geometric derivatives
         if (order_geo/=0) then
           call geom_total_num_derv(order_geo, min(num_atom, order_geo), num_atom, num_addr)
         else
           num_addr = 1
         end if
         ! gets the number of derivatives/powers/moments
         do ierr = 1, num_derv
           num_addr = num_addr*(order_derv(ierr)+1)*(order_derv(ierr)+2)/2
         end do
     
         allocate(address_list(num_order,num_addr), stat=ierr)

         if (ierr/=0) stop "failed to allocate address_list!"
         ! gets the list of addresses
         call get_address_list(num_derv, order_derv,                     &
                               num_atom, order_geo, min(num_atom, order_geo), &
                               num_order, num_addr, address_list)
     
         allocate(tmp_index(nf))
     
         do i = 1, num_addr
     
            k = 1
     
            do j = order_mom + 1, order_mom + order_geo
     
               tmp_index(k) = address_list(j, i)
               k = k + 1
     
            end do
     
            do j = 1, order_mom
     
               tmp_index(k) = address_list(j, i)
               k = k + 1
     
            end do
     
            ave_offset = get_triang_blks_offset(nblks, nf, blk_info, blk_sizes, tmp_index)
            oneint(ave_offset) = oneint_tmp(i)

         end do

         do i = 1, size(oneint_tmp)
            oneint_tmp = 0
         end do
        
         deallocate(address_list)
         deallocate(tmp_index)

      end if
   end subroutine
#endif



  !> \brief reorders and assigns the expectation values and/or integral matrices from
  !>        Gen1int to OpenRsp
  !> \author Bin Gao
  !> \date 2012-04-25
  !> \param num_coord is the number of atomic coordinates
  !> \param num_field is the number of fields
  !> \param first_comp constains the first component in each field ("GEO", "EL")
  !> \param num_comp contains the number of components in each field ("GEO", "EL")
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo is the order of total geometric derivatives
  !> \param val_expect contains the redudant expectation values from Gen1Int
  !> \param val_ints contains the redudant integral matrices from Gen1Int
  !> \return rsp_expect contains the redudant expectation values for OpenRsp
  !> \return rsp_ints contains the redudant integral matrices for OpenRsp
  !> \note all the expectation values and integral matrices are in the order ("EL", "GEO")
  !>       in memory; this routine is implemented temporarily since it will cost much memory
  !>       for assigning integral matrices, it would be better that the integral
  !>       code returns the required integral matrices by itself
  subroutine gen1int_reorder(num_coord, num_field, first_comp, num_comp, &
                             order_mom, order_geo, val_expect, val_ints, &
                             rsp_expect, rsp_ints)
    integer, intent(in) :: num_coord
    integer, intent(in) :: num_field
    integer, intent(in) :: first_comp(num_field)
    integer, intent(in) :: num_comp(num_field)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo
    real(8), optional, intent(in)  :: val_expect(:)
    type(matrix), optional, intent(in) :: val_ints(:)
    complex(8), optional, intent(out) :: rsp_expect(product(num_comp))
    type(matrix), optional, intent(inout) :: rsp_ints(product(num_comp))
    logical assign_expect                   !if assigning expectation values
    logical assign_ints                     !if assigning integral matrices
    integer size_rsp_comp                   !size of all components for OpenRsp
    integer, allocatable :: offset_geom(:)  !offset of total geometric derivatives
    integer, allocatable :: last_comp(:)    !last component in each field
    integer, allocatable :: idx_comp(:)     !indices of specific components
    integer addr_rsp                        !address of specific components in the results for OpenRsp
    integer addr_gen1int                    !address of specific components in the results from Gen1Int
    integer icomp                           !incremental recorder over components
    integer ymom, zmom                      !orders of yz components of Cartesian multipole moments
    logical not_found                       !indicates the next specific components are not found
    integer ierr                            !error information
    ! sets what jobs are going to do
    assign_expect = present(val_expect) .and. present(rsp_expect)
    assign_ints = present(val_ints) .and. present(rsp_ints)
    ! sets the size of all components in the fields
    if (assign_expect) then
      size_rsp_comp = size(rsp_expect)
      if (assign_ints) then
        if (size_rsp_comp/=size(rsp_ints)) &
          call quit("gen1int_reorder>> requires the same size of expectation "// &
                    "values and integral matrices!")
      end if
    else if (assign_ints) then
      size_rsp_comp = size(rsp_ints)
    else
      return
    end if
    ! sets the offset of total geometric derivatives if there are
    if (order_geo>0) then
      allocate(offset_geom(order_geo), stat=ierr)
      if (ierr/=0) call quit("gen1int_reorder>> failed to allocate offset_geom!")
      offset_geom(order_geo) = (order_mom+1)*(order_mom+2)/2  !Cartesian multipole moments come first in memory
      do icomp = order_geo-1, 1, -1
        offset_geom(icomp) = num_coord*offset_geom(icomp+1)
      end do
    end if
    ! sets the last component in each field
    allocate(last_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate last_comp!")
    do icomp = 1, num_field
      last_comp(icomp) = first_comp(icomp)+num_comp(icomp)-1
    end do
    ! indices of specific components
    allocate(idx_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate idx_comp!")
    idx_comp = first_comp
    ! gets the orders of yz components of Cartesian multipole moments
    ymom = 0
    zmom = 0
    do icomp = num_field, order_geo+1, -1
      select case(idx_comp(icomp))
      case(1)
      case(2)
        ymom = ymom+1
      case(3)
        zmom = zmom+1
      case default
        call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
      end select
    end do
    ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
    addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
    ! gets the address of total geometric derivatives in the results from Gen1Int
    do icomp = order_geo, 1, -1
      addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
    end do
    ! assigns the expectation values
    if (assign_expect) rsp_expect(1) = val_expect(addr_gen1int)
    ! assigns the integral matrices
    if (assign_ints) rsp_ints(1) = val_ints(addr_gen1int)
    ! sets other components
    do addr_rsp = 2, size_rsp_comp
      ! generates new specific components
      not_found = .true.
      do icomp = num_field, 1, -1
        ! we walk to the next component in this field, the new specific components are found
        if (idx_comp(icomp)<last_comp(icomp)) then
          idx_comp(icomp) = idx_comp(icomp)+1
          not_found = .false.
          exit
        ! we have walked to the last component of this field, back to the first component
        else
          idx_comp(icomp) = first_comp(icomp)
        end if
      end do
      ! we could not find new specific components
      if (not_found) call quit("gen1int_reorder>> can not find next components!")
      ! gets the orders of yz components of Cartesian multipole moments
      ymom = 0
      zmom = 0
    do icomp = num_field, order_geo+1, -1
        select case(idx_comp(icomp))
        case(1)
        case(2)
          ymom = ymom+1
        case(3)
          zmom = zmom+1
        case default
          call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
        end select
      end do
      ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
      addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
      ! gets the address of total geometric derivatives in the results from Gen1Int
      do icomp = order_geo, 1, -1
        addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
      end do
      ! assigns the expectation values
      if (assign_expect) rsp_expect(addr_rsp) = val_expect(addr_gen1int)
      ! assigns the integral matrices
      if (assign_ints) rsp_ints(addr_rsp) = val_ints(addr_gen1int)
    end do
    ! frees space
    if (allocated(offset_geom)) deallocate(offset_geom)
    deallocate(last_comp)
    deallocate(idx_comp)
  end subroutine



  !> \brief reorders and assigns the expectation values and/or integral matrices from
  !>        Gen1int to OpenRsp (modified by MaR for tensor symmetry nonredundancy
  !>        support. Does not work well yet but will not cause memory problems elsewhere)
  !> \author Bin Gao (small modifications by MaR)
  !> \date 2012-04-25
  !> \param num_coord is the number of atomic coordinates
  !> \param num_field is the number of fields
  !> \param first_comp constains the first component in each field ("GEO", "EL")
  !> \param num_comp contains the number of components in each field ("GEO", "EL")
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo is the order of total geometric derivatives
  !> \param val_expect contains the redudant expectation values from Gen1Int
  !> \param val_ints contains the redudant integral matrices from Gen1Int
  !> \return rsp_expect contains the redudant expectation values for OpenRsp
  !> \return rsp_ints contains the redudant integral matrices for OpenRsp
  !> \note all the expectation values and integral matrices are in the order ("EL", "GEO")
  !>       in memory; this routine is implemented temporarily since it will cost much memory
  !>       for assigning integral matrices, it would be better that the integral
  !>       code returns the required integral matrices by itself
  subroutine gen1int_reorder_tr(num_coord, num_field, first_comp, num_comp, &
                             order_mom, order_geo, propsize, val_expect, val_ints, &
                             rsp_expect, rsp_ints)
    integer, intent(in) :: num_coord
    integer, intent(in) :: num_field
    integer, intent(in) :: first_comp(num_field)
    integer, intent(in) :: num_comp(num_field)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo
    integer :: propsize
    real(8), optional, intent(in)  :: val_expect(:)
    type(matrix), optional, intent(in) :: val_ints(:)
    complex(8), optional, intent(out) :: rsp_expect(propsize)
    type(matrix), optional, intent(inout) :: rsp_ints(propsize)
    logical assign_expect                   !if assigning expectation values
    logical assign_ints                     !if assigning integral matrices
    integer size_rsp_comp                   !size of all components for OpenRsp
    integer, allocatable :: offset_geom(:)  !offset of total geometric derivatives
    integer, allocatable :: last_comp(:)    !last component in each field
    integer, allocatable :: idx_comp(:)     !indices of specific components
    integer addr_rsp                        !address of specific components in the results for OpenRsp
    integer addr_gen1int                    !address of specific components in the results from Gen1Int
    integer icomp                           !incremental recorder over components
    integer ymom, zmom                      !orders of yz components of Cartesian multipole moments
    logical not_found                       !indicates the next specific components are not found
    integer ierr                            !error information
    ! sets what jobs are going to do
    assign_expect = present(val_expect) .and. present(rsp_expect)
    assign_ints = present(val_ints) .and. present(rsp_ints)
    ! sets the size of all components in the fields
    if (assign_expect) then
      size_rsp_comp = size(rsp_expect)
      if (assign_ints) then
        if (size_rsp_comp/=size(rsp_ints)) &
          call quit("gen1int_reorder>> requires the same size of expectation "// &
                    "values and integral matrices!")
      end if
    else if (assign_ints) then
      size_rsp_comp = size(rsp_ints)
    else
      return
    end if
    ! sets the offset of total geometric derivatives if there are
    if (order_geo>0) then
      allocate(offset_geom(order_geo), stat=ierr)
      if (ierr/=0) call quit("gen1int_reorder>> failed to allocate offset_geom!")
      offset_geom(order_geo) = (order_mom+1)*(order_mom+2)/2  !Cartesian multipole moments come first in memory
      do icomp = order_geo-1, 1, -1
        offset_geom(icomp) = num_coord*offset_geom(icomp+1)
      end do
    end if
    ! sets the last component in each field
    allocate(last_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate last_comp!")
    do icomp = 1, num_field
      last_comp(icomp) = first_comp(icomp)+num_comp(icomp)-1
    end do
    ! indices of specific components
    allocate(idx_comp(num_field), stat=ierr)
    if (ierr/=0) call quit("gen1int_reorder>> failed to allocate idx_comp!")
    idx_comp = first_comp
    ! gets the orders of yz components of Cartesian multipole moments
    ymom = 0
    zmom = 0
    do icomp = num_field, order_geo+1, -1
      select case(idx_comp(icomp))
      case(1)
      case(2)
        ymom = ymom+1
      case(3)
        zmom = zmom+1
      case default
        call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
      end select
    end do
    ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
    addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
    ! gets the address of total geometric derivatives in the results from Gen1Int
    do icomp = order_geo, 1, -1
      addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
    end do
    ! assigns the expectation values
    if (assign_expect) rsp_expect(1) = val_expect(addr_gen1int)
    ! assigns the integral matrices
    if (assign_ints) rsp_ints(1) = val_ints(addr_gen1int)
    ! sets other components
    do addr_rsp = 2, size_rsp_comp
      ! generates new specific components
      not_found = .true.
      do icomp = num_field, 1, -1
        ! we walk to the next component in this field, the new specific components are found
        if (idx_comp(icomp)<last_comp(icomp)) then
          idx_comp(icomp) = idx_comp(icomp)+1
          not_found = .false.
          exit
        ! we have walked to the last component of this field, back to the first component
        else
          idx_comp(icomp) = first_comp(icomp)
        end if
      end do
      ! we could not find new specific components
      if (not_found) call quit("gen1int_reorder>> can not find next components!")
      ! gets the orders of yz components of Cartesian multipole moments
      ymom = 0
      zmom = 0
    do icomp = num_field, order_geo+1, -1
        select case(idx_comp(icomp))
        case(1)
        case(2)
          ymom = ymom+1
        case(3)
          zmom = zmom+1
        case default
          call quit("gen1int_reorder>> unknown components for Cartesian multipole moments!")
        end select
      end do
      ! the address of x^{l}y^{m}z^{n} with l+m+n=order_mom is 1+m+(2*order_mom+3-n)*n/2
      addr_gen1int = 1+ymom+(2*order_mom+3-zmom)*zmom/2
      ! gets the address of total geometric derivatives in the results from Gen1Int
      do icomp = order_geo, 1, -1
        addr_gen1int = addr_gen1int+(idx_comp(icomp)-1)*offset_geom(icomp)
      end do
      ! assigns the expectation values
      if (assign_expect) rsp_expect(addr_rsp) = val_expect(addr_gen1int)
      ! assigns the integral matrices
      if (assign_ints) rsp_ints(addr_rsp) = val_ints(addr_gen1int)
    end do
    ! frees space
    if (allocated(offset_geom)) deallocate(offset_geom)
    deallocate(last_comp)
    deallocate(idx_comp)
  end subroutine


   subroutine oneint_ave(nr_atoms, what, D, DFD, R)
      integer,           intent(in)  :: nr_atoms
      character(*),      intent(in)  :: what
      type(matrix),      intent(in)  :: D, DFD
      real(8),           intent(out) :: R(:)
#ifndef VAR_LSDALTON
#include "mxcent.h"
#include "taymol.h"
#endif
      real(8), pointer :: wrk(:)
      integer          :: lwrk
#ifdef VAR_LSDALTON
      call quit('Cannot run oneint_ave, only new integral code is compiled',-1)
#else
      call save_D_and_DFD_for_ABACUS(.false., D, DFD)
      lwrk = 50*D%nrow**2+10000*D%nrow+50000000
      call f77_memory_select(work_len=lwrk, work=wrk)
      HESMOL(:3*nr_atoms,:3*nr_atoms) = 0
      ! SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
      ! &                  NODV,DIFDIP,HFONLY,NCLONE)
      call ONEDRV(wrk,lwrk,0,.true.,len(what),.true.,.true., &
                  .true.,.false.,.true.,.false.)
      ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
      R(1:9*nr_atoms**2) = reshape(HESMOL(:3*nr_atoms,:3*nr_atoms), (/9*nr_atoms**2/))
      call f77_memory_deselect(work_len=lwrk, work=wrk)
#endif
   end subroutine

   subroutine save_D_and_DFD_for_ABACUS(anti, D, DFD)
   !ajt aug09 Same as set_dsofso in fock-eval.f90, but pertrubed
   !          D and DFD is input instead of unperturbed D and F
      logical,      intent(in) :: anti !whether the integrals D and DFD
         ! are to be contracted with are symmetric or anti-symmetric
      type(matrix), intent(in) :: D, DFD
         ! perturbed (or un-) density and
         ! 'generalized Fock matrix' (energy-weighted density matrix)
      real(8)   Dtri(D%nrow*(D%nrow+1)/2)
      real(8) DFDtri(D%nrow*(D%nrow+1)/2)
      if (iszero(D)) then
         Dtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, D%elms, Dtri)
         if (     anti) call DGETAP(D%nrow, D%elms, Dtri)
         ! scale elms by 4 if anti, 2 if symm
         Dtri = Dtri * merge(4,2,anti)
      end if
      if (iszero(DFD)) then
         DFDtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, DFD%elms, DFDtri)
         if (     anti) call DGETAP(D%nrow, DFD%elms, DFDtri)
         ! scale elms by 4 if anti, 2 if symm
         DFDtri = DFDtri * merge(4,2,anti)
      end if
      ! write fo files
#if defined(VAR_LSDALTON) || defined(PRG_DIRAC)
      call quit('Cannot call write_dsofso, only new integral code is compiled',-1)
#else
      call write_dsofso(Dtri,DFDtri)
#endif
   end subroutine

#ifdef VAR_LSDALTON
  subroutine ONEDRV_ave_ifc(fld, siz, ave, D, DFD)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix), optional, intent(in)  :: D, DFD
    STOP 'ONEDRV_ave_ifc not implemented for LSDALTON'
  end subroutine
#else
  subroutine ONEDRV_ave_ifc(fld, siz, ave, D, DFD)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix), optional, intent(in)  :: D, DFD
    !--------------------------------------------
#ifndef VAR_LSDALTON
#include "mxcent.h"
#include "taymol.h"
#endif
    real(8), allocatable :: Dtri(:)
    real(8), allocatable :: DFDtri(:)

    logical anti
    integer nc

    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    allocate(DFDtri(get_nr_ao()*(get_nr_ao()+1)/2))

    ! create triangularly packed matrices from D, DFD
    anti = (mod(count(fld=='MAG '),2) == 1)
    if (.not.present(D)) then
       Dtri = 0
    else if (anti) then
       call DGETAP(D%nrow, D%elms, Dtri)
    else !symm
       call DGEFSP(D%nrow, D%elms, Dtri)
    end if
    if (.not.present(DFD)) then
       DFDtri = 0
    else if (anti) then
       call DGETAP(DFD%nrow, DFD%elms, DFDtri)
    else !symm
       call DGEFSP(DFD%nrow, DFD%elms, DFDtri)
    end if
    ! write to files
#ifdef PRG_DIRAC
    print *, 'fix WRITE_DSOFSO'
    stop 1
#else
    call WRITE_DSOFSO(Dtri, DFDtri)
#endif
    nc = 3 * get_nr_atoms()
    HESMOL(:nc,:nc) = 0
    !  SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,
    ! &                  DIFINT,NODC,NODV,DIFDIP,DIFQDP,
    ! &                  HFONLY,NCLONE,PCM)
    call ONEDRV(f77_memory, size(f77_memory), 5, .true., size(fld), &
                .true., .true., .true., .false., .false., &
                .true., .false., .false.)
    ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
    if (size(fld)==1) &
       call quit('error in ONEDRV_ave_ifc: not implemented')
    if (size(fld)==2) &
       ave = reshape(2*HESMOL(:nc,:nc), (/nc*nc/)) !factor 2 for total dens
    deallocate(Dtri)
    deallocate(DFDtri)
  end subroutine
#endif

  !> Call GET1IN in ABACUS
#ifdef VAR_LSDALTON
  subroutine GET1IN_ave_ifc(fld, siz, ave, D)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix),           intent(in)  :: D
    STOP 'GET1IN_ave_ifc(fld, siz, ave, D) NOT implemented in LSDALTON'
  end subroutine
#else
  subroutine GET1IN_ave_ifc(fld, siz, ave, D)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix),           intent(in)  :: D
    !--------------------------------------------
    real(8)                             :: dummy(2)
#include "mxcent.h"
    real(8), allocatable :: Dtri(:)
    character(8), dimension(9*MXCENT) :: labint
    integer, dimension(9*MXCENT) :: intrep, intadr
    integer ncomp
    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    !ajt Dzero is dangerous! Should have been izero. =0 safe
    intrep = 0 !call dzero(intrep,9*MXCENT)
    intadr = 0 !call dzero(intadr,9*MXCENT)
    ! create triangularly packed matrix from D
    call DGEFSP(D%nrow, D%elms, Dtri)
    ncomp = 0
!      SUBROUTINE GET1IN(SINTMA,WORD,NCOMP,WORK,LWORK,LABINT,INTREP,
!     &                  INTADR,MPQUAD,TOFILE,KPATOM,TRIMAT,EXPVAL,
!     &                  EXP1VL,DENMAT,NPRINT)
#ifdef PRG_DIRAC
    print *, 'fix get1in call'
    stop 1
#else
    call GET1IN(dummy,'DPLGRA ',ncomp,f77_memory,size(f77_memory),labint,intrep, &
                       intadr,0,.false.,0,.false.,ave, &
                       .true.,Dtri,0)
#endif
    deallocate(Dtri)
  end subroutine
#endif

  !> \brief gets property integrals from file AOPROPER
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param prop_lab is the label of integral
  !> \param init_prop indicates if initialize the integral matrix
  !> \return prop_int contains the integral matrix
  subroutine legacy_read_integrals( prop_lab, prop_int )

    character*(8), intent(in)    :: prop_lab
    type(matrix),  intent(inout) :: prop_int

#ifdef PRG_DALTON
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include "inforb.h"
    ! IO units, use LUPROP for file AOPROPER
#include "inftap.h"
    ! external DALTON function finding the corresponding label
    logical FNDLB2
    ! information when calling subroutine FNDLB2
    character*8 RTNLBL(2)
    ! dummy stuff
    integer IDUMMY
    logical anti
#endif

#ifdef VAR_LSDALTON
    STOP 'legacy_read_integrals not implemented in LSDALTON'
#endif

#ifdef PRG_DALTON
    ! one-electron Hamiltonian
    if ( prop_lab == 'ONEHAMIL' ) then
      call QUIT( 'Not implemented!' )
      anti = .false. !radovan: was undefined, setting it to false
      call RDONEL( 'ONEHAMIL', ANTI, f77_memory(get_f77_memory_next()), NNBASX )
      call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
    else
      ! closes file AOPROPER first
      if ( LUPROP > 0 ) call GPCLOSE( LUPROP, 'KEEP' )
      call GPOPEN( LUPROP, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', IDUMMY, .false. )
      ! finds the label
      if ( FNDLB2( prop_lab, RTNLBL, LUPROP ) ) then
        ! square matrix
        if ( RTNLBL(2) == 'SQUARE' ) then
          call READT( LUPROP, N2BASX, prop_int%elms )
        ! symmetric matrix
        else if ( RTNLBL(2) == 'SYMMETRI' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        ! anti-symmetric matrix
        else if ( RTNLBL(2) == 'ANTISYMM' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DAPTGE( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        else
          call QUIT( 'Error: No symmetry label on AOPROPER!' )
        end if
        ! closes file AOPROPER
        call GPCLOSE( LUPROP, 'KEEP' )
      else
        call GPCLOSE( LUPROP, 'KEEP' )
        call QUIT( 'Integrals with label '''//prop_lab// &
                   ''' not found on file ''AOPROPER''!' )
      end if
    end if
#endif /* #ifdef PRG_DALTON */

#ifdef PRG_DIRAC
     if (.not. isdef(prop_int)) then
        call mat_init(prop_int, prop_int%nrow, prop_int%nrow)
     end if
     call read_1el_integrals(prop_lab, prop_int)
#endif /* #ifdef PRG_DIRAC */

  end subroutine

  function prefix_zeros(n, l)
    integer, intent(in) :: n, l !number, length
    character(l)        :: prefix_zeros !resulting n in ascii
    character(1), parameter :: char0to9(0:9) &
          = (/'0','1','2','3','4','5','6','7','8','9'/)
    integer :: i, k
    k = n
    do i = l, 1, -1
       prefix_zeros(i:i) = char0to9(mod(k,10))
       k = k / 10
    end do
    if (k /= 0) call quit('prefix_zeros error: Argument integer does not fit ' &
                       // 'in the specified number of ASCII caracters')
  end function

end module
