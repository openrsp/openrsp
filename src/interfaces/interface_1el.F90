module interface_1el

   use openrsp_cfg
   use matrix_defop
   use interface_molecule
   use interface_basis
   use interface_f77_memory
   use interface_pcm
   use interface_io
   use rsp_indices_and_addressing, only: get_triang_blks_offset

   implicit none

   public interface_1el_ovlave
   public interface_1el_ovlave_tr
   public interface_1el_oneave
   public interface_1el_oneave_tr
   public interface_1el_ovlint
   public interface_1el_ovlint_tr
   public interface_1el_ovlint_half_diff
   public interface_1el_oneint
   public interface_1el_oneint_tr

   public oneint_ave
   public get1in_ave_ifc
   public onedrv_ave_ifc
   public legacy_read_integrals

   private

contains

#ifdef VAR_LSDALTON
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave(nf, f, c, nc, DFD, ave, w, D)
      ! Gen1Int interface
      use gen1int_host
      !> number of fields
      integer,       intent(in)  :: nf
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> energy-weighted density matrix
      type(matrix),  intent(in)  :: DFD
      !> output average
      complex(8),    intent(out) :: ave(product(nc))
      !> field frequencies corresponding to each field
      complex(8),    intent(in), optional  :: w(nf)
      !> density matrix to contract half-differentiated overlap against
      type(matrix),  intent(in), optional  :: D
      STOP 'interface_1el_ovlave not implemted in LSDALTON. TK'      
   end subroutine
#else
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave(nf, f, c, nc, DFD, ave, w, D)
      ! Gen1Int interface
      use gen1int_api
      !> number of fields
      integer,       intent(in)  :: nf
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> energy-weighted density matrix
      type(matrix),  intent(in)  :: DFD
      !> output average
      complex(8),    intent(out) :: ave(product(nc))
      !> field frequencies corresponding to each field
      complex(8),    intent(in), optional  :: w(nf)
      !> density matrix to contract half-differentiated overlap against
      type(matrix),  intent(in), optional  :: D
      !----------------------------------------------
      type(matrix) A(2)
      real(8), parameter   :: fdistep = 2d0**(-25)
      real(8), allocatable :: tmp(:,:,:)
      real(8), allocatable :: fdi(:,:)
      integer i
      integer order_geo                      !order of total geometric derivatives
      integer num_atom                       !number of atoms
      integer num_coord                      !number of atomic coordinates
      integer num_geom                       !number of total geometric derivatives
      integer num_expt                       !number of all expectation values
      real(8), allocatable :: val_expt(:,:)  !expectation values, real numbers
      integer ierr                           !error information
      logical :: all_frequencies_zero

      if (present(w) .and. .not. present(D)) then
         call quit("error in interface_1el_ovlave: frequencies 'w' and density 'D' " &
                // 'must both be present or both absent')
      end if

      all_frequencies_zero = .true.
      if (present(w)) then
         do i = 1, nf
            if ((dabs(real(w(i))) > tiny(0.0d0))) then
               all_frequencies_zero = .false.
            end if
         end do
      end if

      if (any(f == 'EL  ')) then
         ave = 0.0
      else

         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')

         if (order_geo /= nf) then
            call quit("interface_1el_ovlave>> only geometric derivatives implemented!")
         end if

         ! sets the number of total geometric derivatives
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! allocates memory for expectation values
         num_expt = num_geom
         allocate(val_expt(num_expt, 1), stat=ierr)
         if (ierr /= 0) call quit("interface_1el_ovlave>> failed to allocate val_expt!")
         val_expt = 0.0

         !FIXME changes to call Gen1Int !FIXME \sum_{j+k=n}
         ! (-\sum_{j}w_{j}+\sum_{k}w_{k})/2 S(>)^{j}(<)^{k} ! When it comes to w, I think
         ! it's a better idea to ask the integral program for S>> and S<>, and multiply
         ! the resulting average or integral by w afterwards, rather than to send w into
         ! the integral program.  ! with field frequencies
         if (.not. all_frequencies_zero) then
            if (order_geo==1) then
               ! allocate matrices for integrals
               A(1) = mat_alloc_like(DFD)
               A(2) = mat_alloc_like(DFD)
               ! loop over nuclear coordinates
               do i = 0, nc(1)-1
                  ! (half-) perturbed overlap -i/2 Tg into A(1), Sg in A(2)
                  call legacy_read_integrals('SQHDR' // prefix_zeros(c(1)+i,3), A(1))
                  A(1) = -A(1) !SQHDR is really -dS>/dg
                  A(2) = (-w(1)/2) * (A(1) + trps(A(1)))
                  A(1) = A(1) + trps(A(1)) !=1DOVL
                  ave(1+i) = -tr(A(1),DFD)
                  ave(1+i) = ave(1+i) + tr(A(2),D)
               end do
               A(1:2) = 0 !deallocate
            else
               call quit('interface_1el_ovlave>> GEO(>1) with freqencies not implemented!')
            end if
         else
            ! calculates the expectaion values of overlap matrix
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,       &
                                       0,                          &  !multipole moments
                                       0,                          &
                                       0, 0, 0,                    &  !magnetic derivatives
                                       0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                       0, 0,                       &  !partial geometric derivatives
                                       min(2,order_geo,num_atom),  &  !total geometric derivatives
                                       order_geo,                  &
                                       0, (/0/),                   &
                                       REDUNDANT_GEO,              &
                                       .false., .false., .false.,  &  !not implemented yet
                                       1, (/DFD/), num_expt,       &  !expectation values
                                       val_expt, .false.,          &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
#else
                                       1, (/1, 1/),                &
#endif
                                       get_print_unit(), 0)
            val_expt = -val_expt
            ! assigns the output average
            if (order_geo==0) then
               ave = val_expt(:,1)
            else
               call gen1int_reorder(num_coord=num_coord, num_field=nf, &
                                    first_comp=c, num_comp=nc,         &
                                    order_mom=0, order_geo=order_geo,  &
                                    val_expect=val_expt(:,1), rsp_expect=ave)
            end if
         end if
         ! frees space
         deallocate(val_expt)

      end if
   end subroutine
#endif




#ifdef VAR_LSDALTON
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave_tr(nf, f, c, nc, DFD, nblks, blk_info, & 
                                      blk_sizes,  propsize, ave, w, D)
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
      type(matrix),  intent(in), optional  :: DFD
      !> output average
      complex(8),    intent(inout), optional :: ave(propsize)
      !> field frequencies corresponding to each field
      complex(8),    intent(in), optional  :: w(nf)
      !> density matrix to contract half-differentiated overlap against
      type(matrix),  intent(in), optional  :: D
      STOP 'interface_1el_ovlave_tr not implemented for VAR_LSDALTON'
   end subroutine
#else
   ! MaR: ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
   !> \brief host program routine to get the average f-perturbed overlap integrals
   !>        with perturbed density D and energy-weighted density DFD
   subroutine interface_1el_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ave, DFD, w, D)
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
      type(matrix),  intent(in), optional  :: DFD
      !> output average
      complex(8),    intent(inout), optional :: ave(propsize)
      !> field frequencies corresponding to each field
      complex(8),    intent(in), optional  :: w(nf)
      !> density matrix to contract half-differentiated overlap against
      type(matrix),  intent(in), optional  :: D
      !----------------------------------------------
      type(matrix) A(2)
      real(8), parameter   :: fdistep = 2d0**(-25)
      real(8), allocatable :: tmp(:,:,:)
      real(8), allocatable :: fdi(:,:)
      integer i
      integer order_geo                      !order of total geometric derivatives
      integer num_atom                       !number of atoms
      integer num_coord                      !number of atomic coordinates
      integer num_geom                       !number of total geometric derivatives
      integer num_expt                       !number of all expectation values
      integer :: num_addr, num_order
      integer, allocatable :: address_list(:,:)
      real(8), allocatable :: val_expt(:,:)  !expectation values, real numbers
      integer ierr                           !error information
      logical :: all_frequencies_zero

      if (present(w) .and. .not. present(D)) then
         call quit("error in interface_1el_ovlave: frequencies 'w' and density 'D' " &
                // 'must both be present or both absent')
      end if

      all_frequencies_zero = .true.
      if (present(w)) then
         do i = 1, nf
            if ((dabs(real(w(i))) > tiny(0.0d0))) then
               all_frequencies_zero = .false.
            end if
         end do
      end if

      if (any(f == 'EL  ')) then
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

         if (order_geo /= nf) then
            call quit("interface_1el_ovlave>> only geometric derivatives implemented!")
         end if

         ! sets the number of total geometric derivatives
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! allocates memory for expectation values
         num_expt = num_geom
         allocate(val_expt(propsize, 1), stat=ierr)
         if (ierr /= 0) call quit("interface_1el_ovlave>> failed to allocate val_expt!")
         val_expt = 0.0

         !FIXME changes to call Gen1Int !FIXME \sum_{j+k=n}
         ! (-\sum_{j}w_{j}+\sum_{k}w_{k})/2 S(>)^{j}(<)^{k} ! When it comes to w, I think
         ! it's a better idea to ask the integral program for S>> and S<>, and multiply
         ! the resulting average or integral by w afterwards, rather than to send w into
         ! the integral program.  ! with field frequencies
         if (.not. all_frequencies_zero) then
            if (order_geo==1) then
               ! allocate matrices for integrals
               A(1) = mat_alloc_like(DFD)
               A(2) = mat_alloc_like(DFD)
               ! loop over nuclear coordinates
               do i = 0, nc(1)-1
                  ! (half-) perturbed overlap -i/2 Tg into A(1), Sg in A(2)
                  call legacy_read_integrals('SQHDR' // prefix_zeros(c(1)+i,3), A(1))
                  A(1) = -A(1) !SQHDR is really -dS>/dg
                  A(2) = (-w(1)/2) * (A(1) + trps(A(1)))
                  A(1) = A(1) + trps(A(1)) !=1DOVL
                  ave(1+i) = -tr(A(1),DFD)
                  ave(1+i) = ave(1+i) + tr(A(2),D)
               end do
               A(1:2) = 0 !deallocate
            else
               call quit('interface_1el_ovlave>> GEO(>1) with freqencies not implemented!')
            end if
         else
        ! MR: PARAMETERS FOR ADAPTED gen1int CALL BELOW ARE MAYBE SPECIFIED 
        ! INCORRECTLY (ALSO LAST ARG)
            ! calculates the expectaion values of overlap matrix
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,       &
                                       0,                          &  !multipole moments
                                       0,                          &
                                       0, 0, 0,                    &  !magnetic derivatives
                                       0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                       0, 0,                       &  !partial geometric derivatives
                                       order_geo,                  &  !total geometric derivatives
                                       order_geo,                  &
                                       0, (/0/),                   &
                                       UNIQUE_GEO,              &
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
   !> \brief host program routine to get the average 1-electron integrals perturbed by fields f
   !>        with the (perturbed) density matrix D
   subroutine interface_1el_oneave(nf, f, c, nc, D, ave)
      ! Gen1Int interface
      use gen1int_host
      !> number of fields
      integer,       intent(in)  :: nf
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> density matrix to average over
      type(matrix),  intent(in)  :: D
      !> output average
      complex(8),    intent(out) :: ave(product(nc))
      STOP 'interface_1el_oneave not implemented for LSDALTON'
   end subroutine
#else
   !> \brief host program routine to get the average 1-electron integrals perturbed by fields f
   !>        with the (perturbed) density matrix D
   subroutine interface_1el_oneave(nf, f, c, nc, D, ave)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      !> number of fields
      integer,       intent(in)  :: nf
      !> field labels in std order
      character(4),  intent(in)  :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)  :: c(nf), nc(nf)
      !> density matrix to average over
      type(matrix),  intent(in)  :: D
      !> output average
      complex(8),    intent(out) :: ave(product(nc))
      !----------------------------------------------
      integer order_mom                      !order of Cartesian multipole moments
      integer num_mom                        !number of Cartesian multipole moments
      integer order_geo                      !order of total geometric derivatives
      integer num_atom                       !number of atoms
      integer num_coord                      !number of atomic coordinates
      integer num_geom                       !number of total geometric derivatives
      integer num_expt                       !number of all expectation values
      real(8), allocatable :: val_expt(:, :) !expectation values, real numbers
      real(8), allocatable :: temp(:)
      integer ierr                           !error information
      type(matrix) :: T
      type(matrix) :: P
      integer :: ixyz, i

      if (count(f == 'PNC ') > 0) then

#ifdef PRG_DIRAC
         order_geo = count(f == 'GEO ')
         if (order_geo > 1) then
            print *, 'error in oneave: pnc int geo > 1 not implemented'
            stop 1
         end if

         ave = 0.0d0

         P = mat_alloc_like(D)
         do i = 1, nc(1)
            call get_fc_integrals(P%nrow, P%elms_alpha, openrsp_cfg_pnc_center, i, 0)
            ave(i) = dot(P, D)
         end do
         P = 0
#endif /* ifdef PRG_DIRAC */

      else if (count(f == 'EL  ') > 1) then

         ave = 0.0d0

      else

         ! gets the order of Cartesian multipole moments
         order_mom = count(f=='EL  ')
         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
         if (order_mom+order_geo/=nf) then
            call quit("interface_1el_oneave>> only electric and geometric perturbations implemented!")
         end if

         ! sets the number of operators and derivatives
         num_mom = (order_mom+1)*(order_mom+2)/2
         num_atom = get_nr_atoms()
         num_coord = 3*num_atom
         num_geom = num_coord**order_geo

         ! allocates memory for expectation values
         num_expt = num_mom*num_geom
         allocate(val_expt(num_expt,1), stat=ierr)
         if (ierr/=0) call quit("interface_1el_oneave>> failed to allocate val_expt!")
         val_expt = 0.0

         ! electric perturbations
         if (order_mom/=0) then
            call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                       order_mom,                   &  !multipole moments
                                       0,                           &
                                       0, 0, 0,                     &  !magnetic derivatives
                                       0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                       0, 0,                        &  !partial geometric derivatives
                                       min(3,order_geo,num_atom),   &  !total geometric derivatives
                                       order_geo,                   &
                                       0, (/0/),                    &
                                       REDUNDANT_GEO,               &
                                       .false., .false., .false.,   &  !not implemented yet
                                       1, (/D/), num_expt,          &  !expectation values
                                       val_expt, .false.,           &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
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
                                       0, 0,                      &  !partial geometric derivatives
                                       min(3,order_geo,num_atom), &  !total geometric derivatives
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       REDUNDANT_GEO,             &
                                       .false., .false., .false., &  !not implemented yet
                                       1, (/D/), num_expt,        &  !expectation values
                                       val_expt, .false.,         &
                                       1, (/1, 1/),               &
                                       get_print_unit(), 0)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
            allocate(temp(num_expt*3))

            ! nuclear attraction
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_POT_ENERGY,   &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0,                      &
                                       min(3,order_geo,num_atom), &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       REDUNDANT_GEO,             &
                                       .false., .false., .false., &
                                       1, (/D/), num_expt,        &
                                       temp, .false.,             &
                                       2, (/1, 1, 2, 2/),         &
                                       get_print_unit(), 0)
            val_expt(:, 1) = val_expt(:, 1) + temp(1:num_expt)

            ! beta' matrix
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,      &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0,                      &
                                       min(3,order_geo,num_atom), &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       REDUNDANT_GEO,             &
                                       .false., .false., .false., &
                                       1, (/D/), num_expt,        &
                                       temp, .false.,             &
                                       1, (/2, 2/),               &
                                       get_print_unit(), 0)
            val_expt(:, 1) = val_expt(:, 1) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*temp(1:num_expt)

            ! kinetic energy
            T = mat_alloc_like(D)
            do ixyz = 1, 3
               T%elms_alpha = 0.0d0
               call dcopy(D%nrow*D%ncol, D%elms_alpha(1, 1, 5-ixyz), 1, T%elms_alpha, 1)
               temp = 0.0d0
               call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                          0,                           &
                                          1,                           &
                                          0, 0, 0,                     &
                                          0, 0, 0,                     &
                                          0, 0,                        &
                                          min(3,order_geo,num_atom),   &
                                          order_geo,                   &
                                          0, (/0/),                    &
                                          REDUNDANT_GEO,               &
                                          .false., .false., .false.,   &
                                          1, (/T/), num_expt*3,        &
                                          temp, .false.,               &
                                          2, (/1, 2, 2, 1/),           &
                                          get_print_unit(), 0)
               do i = 1, num_expt
                  val_expt(i, 1) = val_expt(i, 1) + openrsp_cfg_speed_of_light*temp((i-1)*3 + ixyz)
               end do
            end do
            T = 0

            deallocate(temp)

#endif /* ifdef PRG_DIRAC */

         end if

         ! assigns the output average
         call gen1int_reorder(num_coord=num_coord, num_field=nf,        &
                              first_comp=c, num_comp=nc,                &
                              order_mom=order_mom, order_geo=order_geo, &
                              val_expect=val_expt(:,1), rsp_expect=ave)

         ! frees space
         deallocate(val_expt)

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
      integer order_mom                      !order of Cartesian multipole moments
      integer num_mom                        !number of Cartesian multipole moments
      integer order_geo                      !order of total geometric derivatives
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
      else

         ! gets the order of total geometric derivatives
         order_geo = count(f=='GEO ')
         if (order_mom+order_geo/=nf) then
            call quit("interface_1el_oneave>> only electric and geometric perturbations implemented!")
         end if

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
         if (order_mom/=0) then
            call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                       order_mom,                   &  !multipole moments
                                       0,                           &
                                       0, 0, 0,                     &  !magnetic derivatives
                                       0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                       0, 0,                        &  !partial geometric derivatives
                                       order_geo,                   &  !total geometric derivatives
                                       order_geo,                   &
                                       0, (/0/),                    &
                                       UNIQUE_GEO,                  &
                                       .false., .false., .false.,   &  !not implemented yet
                                       1, (/D/), propsize,          &  !expectation values
                                       val_expt, .false.,           &
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
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
                                       0, 0,                      &  !partial geometric derivatives
                                       order_geo,                 &  !total geometric derivatives
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       UNIQUE_GEO,                &
                                       .false., .false., .false., &  !not implemented yet
                                       1, (/D/), propsize,        &  !expectation values
                                       val_expt, .false.,         &
                                       1, (/1, 1/),               &
                                       get_print_unit(), 0)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
            ! MaR: ALLOCATION SIZE WAS CHANGED HERE IN THE PATTERN OF EARLIER CHANGES
            ! IT WAS CHANGED TO THE TENSOR SYMMETRY NONREDUNDANT SIZE
            allocate(temp(propsize*3))

            ! nuclear attraction
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_POT_ENERGY,   &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0,                      &
                                       order_geo,                 &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       REDUNDANT_GEO,             &
                                       .false., .false., .false., &
                                       1, (/D/), propsize,        &
                                       temp, .false.,             &
                                       2, (/1, 1, 2, 2/),         &
                                       get_print_unit(), 0)
            val_expt(:, 1) = val_expt(:, 1) + temp(1:num_expt)

            ! beta' matrix
            temp = 0.0d0
            call gen1int_host_get_expt(NON_LAO, INT_OVERLAP,      &
                                       0,                         &
                                       0,                         &
                                       0, 0, 0,                   &
                                       0, 0, 0,                   &
                                       0, 0,                      &
                                       order_geo,                 &
                                       order_geo,                 &
                                       0, (/0/),                  &
                                       UNIQUE_GEO,                &
                                       .false., .false., .false., &
                                       1, (/D/), propsize,        &
                                       temp, .false.,             &
                                       1, (/2, 2/),               &
                                       get_print_unit(), 0)
            val_expt(:, 1) = val_expt(:, 1) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*temp(1:num_expt)

            ! kinetic energy
            T = mat_alloc_like(D)
            do ixyz = 1, 3
               T%elms_alpha = 0.0d0
               call dcopy(D%nrow*D%ncol, D%elms_alpha(1, 1, 5-ixyz), 1, T%elms_alpha, 1)
               temp = 0.0d0
               call gen1int_host_get_expt(NON_LAO, INT_CART_MULTIPOLE, &
                                          0,                           &
                                          1,                           &
                                          0, 0, 0,                     &
                                          0, 0, 0,                     &
                                          0, 0,                        &
                                          order_geo,                   &
                                          order_geo,                   &
                                          0, (/0/),                    &
                                          UNIQUE_GEO,                  &
                                          .false., .false., .false.,   &
                                          1, (/T/), propsize*3,        &
                                          temp, .false.,               &
                                          2, (/1, 2, 2, 1/),           &
                                          get_print_unit(), 0)
               do i = 1, num_expt
                  val_expt(i, 1) = val_expt(i, 1) + openrsp_cfg_speed_of_light*temp((i-1)*3 + ixyz)
               end do
            end do
            T = 0

            deallocate(temp)
            deallocate(num_derv)

#endif /* ifdef PRG_DIRAC */

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
            ave(ave_offset) = val_expt(i,1)

         end do

         deallocate(address_list)
         deallocate(tmp_index)
         deallocate(val_expt)

      end if
   end subroutine
#endif



#ifdef VAR_LSDALTON
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint(nr_ao, nf, f, c, nc, ovl, w, fock)
      ! Gen1Int interface
      use gen1int_host
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> resulting overlap integral matrices (incoming content deleted)
      type(matrix),  intent(inout) :: ovl(product(nc))
      !> frequencies of each field
      complex(8),    intent(in),    optional :: w(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout), optional :: fock(product(nc))
      STOP 'interface_1el_ovlint not implemented for LSDALTON'
   end subroutine
#else
   !> \brief host program routine to compute differentiated overlap matrices, and optionally
   !>        add half-differentiated overlap contribution to Fock matrices
   subroutine interface_1el_ovlint(nr_ao, nf, f, c, nc, ovl, w, fock)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> resulting overlap integral matrices (incoming content deleted)
      type(matrix),  intent(inout) :: ovl(product(nc))
      !> frequencies of each field
      complex(8),    intent(in),    optional :: w(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout), optional :: fock(product(nc))
      !------------------------------------------------
      integer      i
      integer order_geo  !order of total geometric derivatives
      integer num_atom   !number of atoms
      integer num_coord  !number of atomic coordinates
      integer num_geom   !number of total geometric derivatives
      integer num_ints   !number of integral matrices
      logical :: all_frequencies_zero

      if (present(w) .and. .not. present(fock)) then
         call quit("error in interface_1el_ovlint: frequencies 'w' and Fock matrix 'fock' "// &
                   "must both be present or both absent")
      end if

      all_frequencies_zero = .true.
      if (present(w)) then
         do i = 1, nf
            if ((dabs(real(w(i))) > tiny(0.0d0))) then
               all_frequencies_zero = .false.
            end if
         end do
      end if

      if (any(f=='EL  ')) then

         do i = 1, product(nc)
            call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         end do

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

         if (num_ints/=size(ovl)) then
            call quit("interface_1el_ovlint>> returning specific components not implemented!")
         end if

         !FIXME changes to call Gen1Int
         ! with field frequencies
         if (.not. all_frequencies_zero) then
           if (order_geo==1) then
             ! loop over nuclear coordinates
             do i = 0, nc(1)-1
                ! allocate, if needed
                if (.not.isdef(ovl(1+i))) then
                   call mat_init(ovl(1+i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
                end if
                ! overlap into ovl, half-perturbed overlap -i/2 Tg added to fock
                call legacy_read_integrals('SQHDR' // prefix_zeros(c(1)+i,3), ovl(1+i))
                ovl(1+i)  = -ovl(1+i) !SQHDR is really -dS>/dg
                fock(1+i) = fock(1+i) - w(1)/2 * ovl(1+i)
                fock(1+i) = fock(1+i) + w(1)/2 * trps(ovl(1+i))
                ovl(1+i)  = ovl(1+i)  + trps(ovl(1+i)) !=dS/dg=-1DOVL
             end do
     !FIXME to Andreas: do we need higher order geometric derivatives of overlap integrals with frequencies?
           else
             call quit('interface_1el_ovlint>> GEO(>1) with freqencies not implemented!')
           end if
         else
           ! allocates matrices
           do i = 1, num_ints
             if (.not.isdef(ovl(i))) then
               call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             end if
           end do
           ! calculates the overlap matrix
           call gen1int_host_get_int(NON_LAO, INT_OVERLAP,       &
                                     0,                          &  !multipole moments
                                     0,                          &
                                     0, 0, 0,                    &  !magnetic derivatives
                                     0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                     0, 0,                       &  !partial geometric derivatives
                                     min(2,order_geo,num_atom),  &  !total geometric derivatives
                                     order_geo,                  &
                                     0, (/0/),                   &
                                     REDUNDANT_GEO,              &
                                     .false., .false., .false.,  &  !not implemented yet
                                     num_ints, ovl, .false.,     &  !integral matrices
#ifdef PRG_DIRAC
                                     2, (/1, 1, 2, 2/),          &
#else
                                     1, (/1, 1/),                &
#endif
                                     get_print_unit(), 0)
         end if

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
                                      blk_sizes, propsize, ovl, w, fock)
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
      type(matrix),  intent(inout), optional :: ovl(propsize)
      type(matrix)   :: ovl_tmp(propsize)
      !> frequencies of each field
      complex(8),    intent(in),    optional :: w(nf)
      !> Fock matrices to which the half-differentiated overlap
      !> contribution is ADDED
      type(matrix),  intent(inout), optional :: fock(propsize)
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

      if (present(w) .and. .not. present(fock)) then
         call quit("error in interface_1el_ovlint: frequencies 'w' and Fock matrix 'fock' "// &
                   "must both be present or both absent")
      end if

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



      all_frequencies_zero = .true.
      if (present(w)) then
         do i = 1, nf
            if ((dabs(real(w(i))) > tiny(0.0d0))) then
               all_frequencies_zero = .false.
            end if
         end do
      end if

      if (any(f=='EL  ')) then

         do i = 1, propsize
            call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         end do

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

         ! MaR: DISABLED: NOT A VALID WAY OF DISCERNING THIS ANYMORE
!          if (num_ints/=size(ovl)) then
!             call quit("interface_1el_ovlint>> returning specific components not implemented!")
!          end if

         !FIXME changes to call Gen1Int
         ! with field frequencies
         if (.not. all_frequencies_zero) then
           if (order_geo==1) then
             ! loop over nuclear coordinates
             do i = 0, nc(1)-1
                ! allocate, if needed
                if (.not.isdef(ovl(1+i))) then
                   call mat_init(ovl(1+i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
                end if
                ! overlap into ovl, half-perturbed overlap -i/2 Tg added to fock
                call legacy_read_integrals('SQHDR' // prefix_zeros(c(1)+i,3), ovl(1+i))
                ovl(1+i)  = -ovl(1+i) !SQHDR is really -dS>/dg
                fock(1+i) = fock(1+i) - w(1)/2 * ovl(1+i)
                fock(1+i) = fock(1+i) + w(1)/2 * trps(ovl(1+i))
                ovl(1+i)  = ovl(1+i)  + trps(ovl(1+i)) !=dS/dg=-1DOVL
             end do
     !FIXME to Andreas: do we need higher order geometric derivatives of overlap integrals with frequencies?
           else
             call quit('interface_1el_ovlint>> GEO(>1) with freqencies not implemented!')
           end if
         else
           ! allocates matrices
           do i = 1, propsize
             if (.not.isdef(ovl(i))) then
               call mat_init(ovl(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             end if
             if (.not.isdef(ovl_tmp(i))) then
                call mat_init(ovl_tmp(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             end if
           end do
           ! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALL BELOW
           ! calculates the overlap matrix
           call gen1int_host_get_int(NON_LAO, INT_OVERLAP,       &
                                     0,                          &  !multipole moments
                                     0,                          &
                                     0, 0, 0,                    &  !magnetic derivatives
                                     0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                     0, 0,                       &  !partial geometric derivatives
                                     order_geo,                  &  !total geometric derivatives
                                     order_geo,                  &
                                     0, (/0/),                   &
                                     UNIQUE_GEO,                 &
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
      STOP 'interface_1el_ovlint_tr not implemented for LSDALTON. TK'
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

           do i = 1, propsize
             if (.not.isdef(fock(i))) then
               call mat_init(fock(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             end if
             if (.not.isdef(tmp_fock(i))) then
                call mat_init(tmp_fock(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             end if
           end do
           ! MaR: MAY NEED ANOTHER LOOK AT THE VERY LAST ARGUMENT OF THE GEN1INT CALL BELOW
           ! calculates the overlap matrix
           call gen1int_host_get_int(LONDON, INT_OVERLAP,       &
                                     0,                          &  !multipole moments
                                     0,                          &
                                     count(fbra=='MAG '),        &
                                     count(fket=='MAG '), 0,     &  !magnetic derivatives
                                     0, 0, 0,                    &  !derivatives w.r.t. total RAM
                                     count(fbra=='GEO '),        &
                                     count(fket=='GEO '),        &  !partial geometric derivatives
                                     0,                          &  !total geometric derivatives
                                     0,                          &
                                     0, (/0/),                   &
                                     UNIQUE_GEO,                 &
                                     .false., .false., .false.,  &  !not implemented yet
                                     propsize, tmp_fock, .false.,   &  !integral matrices
#ifdef PRG_DIRAC
                                     2, (/1, 1, 2, 2/),          &
#else
                                     1, (/1, 1/),                &
#endif
                                     get_print_unit(), 0)

          ! MaR: OFFSETS AND LOCATIONS IN MEMORY POSTPONED
          ! MAKE THE CORRESPONDING OVLAVE TYPE ROUTINE

   end subroutine
#endif






















#ifdef VAR_LSDALTON
   subroutine interface_1el_oneint(nr_ao, nf, f, c, nc, oneint)
      ! Gen1Int interface
      use gen1int_host
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> output perturbed integrals
      type(matrix),  intent(inout) :: oneint(product(nc))
      STOP 'interface_1el_oneint not implemented for LSDALTON. TK'
   end subroutine
#else
   subroutine interface_1el_oneint(nr_ao, nf, f, c, nc, oneint)
      ! Gen1Int interface
#ifdef VAR_LSDALTON
      use gen1int_host
#else
      use gen1int_api
#endif
      integer, intent(in)          :: nr_ao
      !> number of fields
      integer,       intent(in)    :: nf
      !> field labels in std order
      character(4),  intent(in)    :: f(nf)
      !> first and number of- components in each field
      integer,       intent(in)    :: c(nf), nc(nf)
      !> output perturbed integrals
      type(matrix),  intent(inout) :: oneint(product(nc))
      !--------------------------------------------------
      integer order_mom  !order of Cartesian multipole moments
      integer num_mom    !number of Cartesian multipole moments
      integer order_geo  !order of total geometric derivatives
      integer num_atom   !number of atoms
      integer num_coord  !number of atomic coordinates
      integer num_geom   !number of total geometric derivatives
      integer num_ints   !number of all integral matrices
      integer imat       !incremental recorder over matrices
      integer :: i, ixyz
      type(matrix) :: A
      type(matrix), allocatable :: T(:)

      if (count(f == 'PNC ') > 0) then

#ifdef PRG_DIRAC
         order_geo = count(f == 'GEO ')
         if (order_geo > 0) then
            print *, 'error in oneint: pnc int geo > 0 not implemented'
            stop 1
         end if

         if (.not. isdef(oneint(1))) then
            call mat_init(oneint(1), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         end if
         call get_fc_integrals(nr_ao, oneint(1)%elms_alpha, openrsp_cfg_pnc_center, 0, 0)
#endif /* ifdef PRG_DIRAC */

      else if (count(f == 'EL  ') > 1) then

         ! radovan: this code does not make sense to me
         !          what is it supposed to do?
         call mat_init(A, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         do i = 1, product(nc)
            if (iszero(oneint(i))) then
               call mat_ensure_alloc(oneint(i))
               oneint(i)%elms_alpha = oneint(i)%elms_alpha + A%elms_alpha
            else
               oneint(i)%elms_alpha = oneint(i)%elms_alpha + A%elms_alpha
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
         if (num_ints/=size(oneint)) then
            call quit("interface_1el_oneint>> returning specific components is not implemented!")
         end if

         !FIXME: it is better that we use unique components for higher order
         if (order_mom>1) then
            call quit("interface_1el_oneint>> only the first Cartesian multipole moments implemented!")
         end if

         do imat = 1, size(oneint)
            if (.not.isdef(oneint(imat))) then
               call mat_init(oneint(imat), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
            end if
         end do

         ! electric perturbations
         if (order_mom/=0) then
            call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                      order_mom,                   &  !multipole moments
                                      0,                           &
                                      0, 0, 0,                     &  !magnetic derivatives
                                      0, 0, 0,                     &  !derivatives w.r.t. total RAM
                                      0, 0,                        &  !partial geometric derivatives
                                      min(3,order_geo,num_atom),   &  !total geometric derivatives
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      REDUNDANT_GEO,               &
                                      .false., .false., .false.,   &  !not implemented yet
                                      num_ints, oneint, .false.,   &  !integral matrices
#ifdef PRG_DIRAC
                                       2, (/1, 1, 2, 2/),          &
#else
                                       1, (/1, 1/),                &
#endif
                                      get_print_unit(), 0)

         ! only geometric perturbations
         else
#ifdef PRG_DIRAC
            allocate(T(3*size(oneint)))
            do i = 1, 3*size(oneint)
               T(i) = mat_alloc_like(oneint(1))
            end do

            ! nuclear attraction
            call gen1int_host_get_int(NON_LAO, INT_POT_ENERGY,   &
                                      0,                         &
                                      0,                         &
                                      0, 0, 0,                   &
                                      0, 0, 0,                   &
                                      0, 0,                      &
                                      min(3,order_geo,num_atom), &
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      REDUNDANT_GEO,             &
                                      .false., .false., .false., &
                                      num_ints, oneint, .false., &
                                      2, (/1, 1, 2, 2/),         &
                                      get_print_unit(), 0)

            ! beta' matrix
            call gen1int_host_get_int(NON_LAO, INT_OVERLAP,      &
                                      0,                         &
                                      0,                         &
                                      0, 0, 0,                   &
                                      0, 0, 0,                   &
                                      0, 0,                      &
                                      min(3,order_geo,num_atom), &
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      REDUNDANT_GEO,             &
                                      .false., .false., .false., &
                                      num_ints, T, .false.,      &
                                      1, (/2, 2/),               &
                                      get_print_unit(), 0)
            do i = 1, size(oneint)
               oneint(i) = oneint(i) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*T(i)
            end do

!           ! kinetic energy
            call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                      0,                           &
                                      1,                           &
                                      0, 0, 0,                     &
                                      0, 0, 0,                     &
                                      0, 0,                        &
                                      min(3,order_geo,num_atom),   &
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      REDUNDANT_GEO,               &
                                      .false., .false., .false.,   &
                                      3*num_ints, T, .false.,      &
                                      2, (/1, 2, 2, 1/),           &
                                      get_print_unit(), 0)
            do i = 1, size(oneint)
               do ixyz = 1, 3
                  call daxpy(T(1)%nrow*T(1)%ncol,                &
                            -openrsp_cfg_speed_of_light,         &
                             T((i-1)*3 + ixyz)%elms_alpha,       &
                             1,                                  &
                             oneint(i)%elms_alpha(1, 1, 5-ixyz), &
                             1)
               end do
            end do

            do i = 1, 3*size(oneint)
               T(i) = 0
            end do
            deallocate(T)
#else
            call gen1int_host_get_int(NON_LAO, INT_ONE_HAMIL,    &
                                      0,                         &  !multipole moments
                                      0,                         &
                                      0, 0, 0,                   &  !magnetic derivatives
                                      0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                      0, 0,                      &  !partial geometric derivatives
                                      min(3,order_geo,num_atom), &  !total geometric derivatives
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      REDUNDANT_GEO,             &
                                      .false., .false., .false., &  !not implemented yet
                                      num_ints, oneint, .false., &  !integral matrices
                                      1, (/1, 1/),               &
                                      get_print_unit(), 0)
#endif
         end if

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
      integer, dimension(nblks) :: blk_sizes
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
         call mat_init(A, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
         do i = 1, propsize
            if (iszero(oneint(i))) then
               call mat_ensure_alloc(oneint(i))
               oneint(i)%elms_alpha = oneint(i)%elms_alpha + A%elms_alpha
            else
               oneint(i)%elms_alpha = oneint(i)%elms_alpha + A%elms_alpha
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
               call mat_init(oneint(imat), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
            end if
             if (.not.isdef(oneint_tmp(imat))) then
                call mat_init(oneint_tmp(imat), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
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
                                      0, 0,                        &  !partial geometric derivatives
                                      order_geo,                   &  !total geometric derivatives
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      UNIQUE_GEO,                  &
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
            allocate(T(3*size(oneint)))
            do i = 1, 3*size(oneint)
               T(i) = mat_alloc_like(oneint(1))
            end do

            ! nuclear attraction
            call gen1int_host_get_int(NON_LAO, INT_POT_ENERGY,   &
                                      0,                         &
                                      0,                         &
                                      0, 0, 0,                   &
                                      0, 0, 0,                   &
                                      0, 0,                      &
                                      order_geo,                 &
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      UNIQUE_GEO,                &
                                      .false., .false., .false., &
                                      num_ints, oneint, .false., &
                                      2, (/1, 1, 2, 2/),         &
                                      get_print_unit(), 0)

            ! beta' matrix
            call gen1int_host_get_int(NON_LAO, INT_OVERLAP,      &
                                      0,                         &
                                      0,                         &
                                      0, 0, 0,                   &
                                      0, 0, 0,                   &
                                      0, 0,                      &
                                      order_geo,                 &
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      UNIQUE_GEO,                &
                                      .false., .false., .false., &
                                      num_ints, T, .false.,      &
                                      1, (/2, 2/),               &
                                      get_print_unit(), 0)
            do i = 1, size(oneint)
               oneint(i) = oneint(i) - 2.0d0*(openrsp_cfg_speed_of_light**2.0d0)*T(i)
            end do

!           ! kinetic energy
            call gen1int_host_get_int(NON_LAO, INT_CART_MULTIPOLE, &
                                      0,                           &
                                      1,                           &
                                      0, 0, 0,                     &
                                      0, 0, 0,                     &
                                      0, 0,                        &
                                      order_geo,                   &
                                      order_geo,                   &
                                      0, (/0/),                    &
                                      UNIQUE_GEO,                  &
                                      .false., .false., .false.,   &
                                      3*num_ints, T, .false.,      &
                                      2, (/1, 2, 2, 1/),           &
                                      get_print_unit(), 0)
            do i = 1, size(oneint)
               do ixyz = 1, 3
                  call daxpy(T(1)%nrow*T(1)%ncol,                &
                            -openrsp_cfg_speed_of_light,         &
                             T((i-1)*3 + ixyz)%elms_alpha,       &
                             1,                                  &
                             oneint(i)%elms_alpha(1, 1, 5-ixyz), &
                             1)
               end do
            end do

            do i = 1, 3*size(oneint)
               T(i) = 0
            end do
            deallocate(T)
#else
            call gen1int_host_get_int(NON_LAO, INT_ONE_HAMIL,    &
                                      0,                         &  !multipole moments
                                      0,                         &
                                      0, 0, 0,                   &  !magnetic derivatives
                                      0, 0, 0,                   &  !derivatives w.r.t. total RAM
                                      0, 0,                      &  !partial geometric derivatives
                                      order_geo,                 &  !total geometric derivatives
                                      order_geo,                 &
                                      0, (/0/),                  &
                                      UNIQUE_GEO,                &
                                      .false., .false., .false., &  !not implemented yet
                                      num_ints, oneint_tmp, .false., &  !integral matrices
                                      1, (/1, 1/),               &
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
         if (.not.anti) call DGEFSP(D%nrow, D%elms_alpha, Dtri)
         if (     anti) call DGETAP(D%nrow, D%elms_alpha, Dtri)
         ! scale elms_alpha by 4 if anti, 2 if symm
         Dtri = Dtri * merge(4,2,anti)
      end if
      if (iszero(DFD)) then
         DFDtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, DFD%elms_alpha, DFDtri)
         if (     anti) call DGETAP(D%nrow, DFD%elms_alpha, DFDtri)
         ! scale elms_alpha by 4 if anti, 2 if symm
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
       call DGETAP(D%nrow, D%elms_alpha, Dtri)
    else !symm
       call DGEFSP(D%nrow, D%elms_alpha, Dtri)
    end if
    if (.not.present(DFD)) then
       DFDtri = 0
    else if (anti) then
       call DGETAP(DFD%nrow, DFD%elms_alpha, DFDtri)
    else !symm
       call DGEFSP(DFD%nrow, DFD%elms_alpha, DFDtri)
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
    call DGEFSP(D%nrow, D%elms_alpha, Dtri)
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

#ifdef VAR_LSDALTON
  subroutine legacy_read_integrals( prop_lab, prop_int )
    character*(8), intent(in) :: prop_lab
    type(matrix), intent(inout) :: prop_int
    STOP 'legacy_read_integrals not implemented in LSDALTON'
  end subroutine
#else
  !> \brief gets property integrals from file AOPROPER
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param prop_lab is the label of integral
  !> \param init_prop indicates if initialize the integral matrix
  !> \return prop_int contains the integral matrix
  subroutine legacy_read_integrals( prop_lab, prop_int )
    character*(8), intent(in) :: prop_lab
    type(matrix), intent(inout) :: prop_int
#ifndef VAR_LSDALTON
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include "inforb.h"
    ! IO units, use LUPROP for file AOPROPER
#include "inftap.h"
#endif
    ! external DALTON function finding the corresponding label
    logical FNDLB2
    ! information when calling subroutine FNDLB2
    character*8 RTNLBL(2)
    ! dummy stuff
    integer IDUMMY
    logical anti
    ! one-electron Hamiltonian
    if ( prop_lab == 'ONEHAMIL' ) then
      call QUIT( 'Not implemented!' )
#ifdef PRG_DIRAC
    print *, 'fix rdonel call'
    stop 1
#else
      anti = .false. !radovan: was undefined, setting it to false
      call RDONEL( 'ONEHAMIL', ANTI, f77_memory(get_f77_memory_next()), NNBASX )
#endif
      call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms_alpha )
    else
      ! closes file AOPROPER first
      if ( LUPROP > 0 ) call GPCLOSE( LUPROP, 'KEEP' )
      call GPOPEN( LUPROP, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', IDUMMY, .false. )
      ! finds the label
      if ( FNDLB2( prop_lab, RTNLBL, LUPROP ) ) then
        ! square matrix
        if ( RTNLBL(2) == 'SQUARE' ) then
          call READT( LUPROP, N2BASX, prop_int%elms_alpha )
        ! symmetric matrix
        else if ( RTNLBL(2) == 'SYMMETRI' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms_alpha )
        ! anti-symmetric matrix
        else if ( RTNLBL(2) == 'ANTISYMM' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DAPTGE( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms_alpha )
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
  end subroutine
#endif
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
