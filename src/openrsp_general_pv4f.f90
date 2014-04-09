subroutine openrsp_general_pv4f(num_atoms, S, D, F)

   use rsp_field_tuple, only: p_tuple, p_tuple_standardorder
   use matrix_defop, matrix => openrsp_matrix
   use rsp_sdf_caching, only: SDF, sdf_setup_datatype
   use rsp_general, only: rsp_prop
   use openrsp_cfg
   use matrix_lowlevel,  only: mat_init
   use vib_pv_contribs
   use rsp_indices_and_addressing, only: mat_init_like_and_zero
   use legacy_vibrational_properties, only: load_vib_modes

   implicit none

!  input/output
!  -----------------------------------------------------------------------------
   integer      :: num_atoms
   type(matrix) :: S, D, F

!  local variables
!  -----------------------------------------------------------------------------
    integer       :: kn(2)
    type(p_tuple) :: perturbation_tuple
    type(matrix) :: zeromat_already
    type(SDF), pointer :: F_already, D_already, S_already
    real(8), dimension(3) :: fld_dum
    real(8), dimension(3*num_atoms) :: geo_dum
    real(8), dimension(3) :: dm
    real(8), dimension(3) :: dm_orig
    character(20) :: format_line
    integer       :: n_nm, h, j, k, m, n, p, q, ierr, i
    real(8), allocatable, dimension(:) :: nm_freq, nm_freq_b
    real(8), allocatable, dimension(:,:) :: T
    real(8), parameter :: aunm = 45.56335d0
    complex(8), allocatable, dimension(:,:) :: egf_cart, egf_nm, egg_cart, egg_nm,  ff_pv
    complex(8), allocatable, dimension(:,:,:) :: egff_cart, egff_nm, fff_pv
    complex(8), allocatable, dimension(:,:,:) :: eggf_cart, eggf_nm, eggg_cart, eggg_nm
    complex(8), allocatable, dimension(:,:,:,:) :: egfff_cart, egfff_nm, ffff_pv
    complex(8), allocatable, dimension(:,:,:,:) :: eggff_cart, eggff_nm, egggf_cart, egggf_nm
    complex(8), allocatable, dimension(:,:,:,:) :: egggg_cart, egggg_nm
    complex(8), allocatable, dimension(:,:,:,:,:) :: eggfff_cart, eggfff_nm, eggggg_cart
    complex(8), allocatable, dimension(:,:,:,:,:) :: egggff_cart, egggff_nm, eggggg_nm
    complex(8), allocatable, dimension(:,:,:,:,:,:) :: egggfff_cart, egggfff_nm
    complex(8), allocatable, dimension(:,:,:,:,:,:) :: egggggg_cart, egggggg_nm
    complex(8), allocatable, dimension(:,:,:,:,:,:,:) :: egggffff_cart, egggffff_nm


       ! Calculate dipole moment

       kn = (/0,0/)

       perturbation_tuple%n_perturbations = 1
       allocate(perturbation_tuple%pdim(1))
       allocate(perturbation_tuple%plab(1))
       allocate(perturbation_tuple%pid(1))
       allocate(perturbation_tuple%freq(1))

       perturbation_tuple%plab = (/'EL  '/)
       perturbation_tuple%pdim = (/3/)
       perturbation_tuple%pid = (/1/)
       perturbation_tuple%freq = (/0.0d0/)

       call rsp_prop(perturbation_tuple, kn, F, D, S, file_id='Ef')

       ! Read dipole moment from file

       open(unit = 258, file='rsp_tensor_Ef', status='old', action='read', iostat=ierr)
       read(258,*) fld_dum
          dm = fld_dum
       close(258)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


       fld_dum = 0.0

       ! Get normal mode transformation matrix
       ! Get normal mode frequencies

       allocate(T(3*num_atoms, 3*num_atoms))

       allocate(nm_freq_b(3*num_atoms))

       nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

       nm_freq = 0.0

       nm_freq = nm_freq_b(1:n_nm)
       deallocate(nm_freq_b)

       allocate(ff_pv(3, 3))
       allocate(fff_pv(3, 3, 3))
       allocate(ffff_pv(3, 3, 3, 3))
       allocate(egf_cart(3*num_atoms, 3))
       allocate(egf_nm(n_nm, 3))
       allocate(egff_cart(3*num_atoms, 3, 3))
       allocate(egff_nm(n_nm, 3, 3))
       allocate(egfff_cart(3*num_atoms, 3, 3, 3))
       allocate(egfff_nm(n_nm, 3, 3, 3))

       ff_pv = 0.0
       fff_pv = 0.0
       ffff_pv = 0.0
       egf_cart = 0.0
       egf_nm = 0.0
       egff_cart = 0.0
       egff_nm = 0.0
       egfff_cart = 0.0
       egfff_nm = 0.0

       ! Calculate gradient of dipole moment

       kn = (/0,1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(2))
       allocate(perturbation_tuple%plab(2))
       allocate(perturbation_tuple%pid(2))
       allocate(perturbation_tuple%freq(2))

       perturbation_tuple%plab = (/'GEO ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3/)
       perturbation_tuple%pid = (/1, 2/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2/)

       ! ASSUME CLOSED SHELL
       call mat_init(zeromat_already, S%nrow, S%ncol, is_zero=.true.)
       call mat_init_like_and_zero(S, zeromat_already)


       call sdf_setup_datatype(S_already, S)
       call sdf_setup_datatype(D_already, D)
       call sdf_setup_datatype(F_already, F)


       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egf')

       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Egf', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          read(258,*) fld_dum
          egf_cart(i, :) = fld_dum
       end do

       close(258)

       egf_nm = trans_cartnc_1w1d(3*num_atoms, n_nm, egf_cart, T(:,1:n_nm))

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

      ! Calculate other properties for given orders of anharmonicity


       allocate(eggf_cart(3*num_atoms, 3*num_atoms, 3))
       allocate(eggf_nm(n_nm, n_nm, 3))
       allocate(eggff_cart(3*num_atoms, 3*num_atoms, 3, 3))
       allocate(eggff_nm(n_nm, n_nm, 3, 3))
       allocate(eggfff_cart(3*num_atoms, 3*num_atoms, 3, 3, 3))
       allocate(eggfff_nm(n_nm, n_nm, 3, 3, 3))
       allocate(egggf_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3))
       allocate(egggf_nm(n_nm, n_nm, n_nm, 3))
       allocate(egggff_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3))
       allocate(egggff_nm(n_nm, n_nm, n_nm, 3, 3))
       allocate(egggfff_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3, 3))
       allocate(egggfff_nm(n_nm, n_nm, n_nm, 3, 3, 3))
       allocate(eggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms))
       allocate(eggg_nm(n_nm, n_nm, n_nm))
       allocate(egggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3*num_atoms))
       allocate(egggg_nm(n_nm, n_nm, n_nm, n_nm))



       if (openrsp_cfg_general_pv_el_anh > 0) then

          ! Calculate Hessian of dipole moment

          kn = (/1,1/)

          perturbation_tuple%n_perturbations = 3
          allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

          perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  '/)
          perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3/)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
          perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, file_id='Eggf')

          open(unit = 258, file='rsp_tensor_Eggf', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                read(258,*) fld_dum
                eggf_cart(i, j, :) = fld_dum
             end do
          end do

          close(258)

          eggf_nm = trans_cartnc_1w2d(3*num_atoms, n_nm, eggf_cart, T(:,1:n_nm))

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)

         ! Calculate Hessian of polarizability

          kn = (/1,2/)

          perturbation_tuple%n_perturbations = 4
          allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

          perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  '/)
          perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3, 3/)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
          perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, file_id='Eggff')

          open(unit = 258, file='rsp_tensor_Eggff', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3
                   read(258,*) fld_dum
                   eggff_cart(i, j, k, :) = (-1.0) * fld_dum
                end do
             end do
          end do

          close(258)

          eggff_nm = trans_cartnc_2w2d(3*num_atoms, n_nm, eggff_cart, T(:,1:n_nm))

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)


          ! Calculate Hessian of first hyperpolarizability

          kn = (/1,3/)

          perturbation_tuple%n_perturbations = 5
          allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

          perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
          perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3, 3, 3/)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
          perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, file_id='Eggfff')

          open(unit = 258, file='rsp_tensor_Eggfff', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3
                   do m = 1, 3
                      read(258,*) fld_dum
                      eggfff_cart(i, j, k, m, :) = fld_dum
                   end do
                end do
             end do
          end do

          close(258)

          eggfff_nm = trans_cartnc_3w2d(3*num_atoms, n_nm, eggfff_cart, T(:,1:n_nm))

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)

          if (openrsp_cfg_general_pv_el_anh > 1) then

             ! Calculate cubic force field of dipole moment

             kn = (/1,2/)

             perturbation_tuple%n_perturbations = 4
             allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

             perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  '/)
             perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3/)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
             perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

             perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

             call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egggf')

             open(unit = 258, file='rsp_tensor_Egggf', status='old', action='read', iostat=ierr)

             do i = 1, 3*num_atoms
                do j = 1, 3*num_atoms
                   do k = 1, 3*num_atoms
                      read(258,*) fld_dum
                      egggf_cart(i, j, k, :) = fld_dum
                   end do
                end do
             end do

             close(258)

             egggf_nm = trans_cartnc_1w3d(3*num_atoms, n_nm, egggf_cart, T(:,1:n_nm))

             deallocate(perturbation_tuple%pdim)
             deallocate(perturbation_tuple%plab)
             deallocate(perturbation_tuple%pid)
             deallocate(perturbation_tuple%freq)


             ! Calculate cubic force field of polarizability

             kn = (/2,2/)

             perturbation_tuple%n_perturbations = 5
             allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

             perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  '/)
             perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3/)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
             perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

             perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

             call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egggff')

             open(unit = 258, file='rsp_tensor_Egggff', status='old', action='read', iostat=ierr)

             do i = 1, 3*num_atoms
                do j = 1, 3*num_atoms
                   do k = 1, 3*num_atoms
                      do m = 1, 3
                         read(258,*) fld_dum
                         egggff_cart(i, j, k, m, :) = (-1.0) * fld_dum
                      end do
                   end do
                end do
             end do

             close(258)

             egggff_nm = trans_cartnc_2w3d(3*num_atoms, n_nm, egggff_cart, T(:,1:n_nm))

             deallocate(perturbation_tuple%pdim)
             deallocate(perturbation_tuple%plab)
             deallocate(perturbation_tuple%pid)
             deallocate(perturbation_tuple%freq)

             ! Calculate cubic force field of first hyperpolarizability

             kn = (/2,3/)

             perturbation_tuple%n_perturbations = 6
             allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

             perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
             perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3, 3/)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
             perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

             perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

             call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egggfff')

             open(unit = 258, file='rsp_tensor_Egggfff', status='old', action='read', iostat=ierr)

             do i = 1, 3*num_atoms
                do j = 1, 3*num_atoms
                   do k = 1, 3*num_atoms
                      do m = 1, 3
                         do n = 1, 3
                            read(258,*) fld_dum
                            egggfff_cart(i, j, k, m, n, :) = fld_dum
                         end do
                      end do
                   end do
                end do
             end do

             close(258)

             egggfff_nm = trans_cartnc_3w3d(3*num_atoms, n_nm, egggfff_cart, T(:,1:n_nm))

             deallocate(perturbation_tuple%pdim)
             deallocate(perturbation_tuple%plab)
             deallocate(perturbation_tuple%pid)
             deallocate(perturbation_tuple%freq)


          end if

       end if

       if (openrsp_cfg_general_pv_mech_anh > 0) then

          ! Calculate cubic force field

          kn = (/1,1/)

          perturbation_tuple%n_perturbations = 3
          allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
          allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

          perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO '/)
          perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms/)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
          perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, file_id='Eggg')

          open(unit = 258, file='rsp_tensor_Eggg', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                read(258,*) geo_dum
                eggg_cart(i, j, :) = geo_dum
             end do
          end do

          close(258)

          eggg_nm = trans_cartnc_0w3d(3*num_atoms, n_nm, eggg_cart, T(:,1:n_nm))

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)


          if (openrsp_cfg_general_pv_mech_anh > 1) then

             ! Calculate quartic force field

             kn = (/1,2/)

             perturbation_tuple%n_perturbations = 4
             allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
             allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

             perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO '/)
             perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3*num_atoms/)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)
             perturbation_tuple%freq = (/(0.0d0*i, i = 1, perturbation_tuple%n_perturbations)/)

             perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
             perturbation_tuple%pid = (/(i, i = 1, perturbation_tuple%n_perturbations)/)

             call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egggg')

             open(unit = 258, file='rsp_tensor_Egggg', status='old', action='read', iostat=ierr)

             do i = 1, 3*num_atoms
                do j = 1, 3*num_atoms
                   do k = 1, 3*num_atoms
                      read(258,*) geo_dum
                      egggg_cart(i, j, k, :)  = geo_dum
                   end do
                end do
             end do

             close(258)

             egggg_nm = trans_cartnc_0w4d(3*num_atoms, n_nm, egggg_cart, T(:,1:n_nm))

             deallocate(perturbation_tuple%pdim)
             deallocate(perturbation_tuple%plab)
             deallocate(perturbation_tuple%pid)
             deallocate(perturbation_tuple%freq)

          end if

       end if

! End calculation of anharmonicity-related response properties


       ! Normalize dipole moment - follows procedure by AJT

       dm_orig = dm

       if ( ((sum(dm * dm))**0.5) < 0.00001 ) then

          dm = (/0d0, 0d0, 1d0/)

       else

          dm = dm/((sum(dm * dm))**0.5)

       end if

       do k = 1, openrsp_cfg_nr_freq_tuples

          do h = 1, 3

             ! Calculate PV contribution to polarizability
             if (openrsp_cfg_general_pv_total_anh > 0) then

                ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0) * &
                        openrsp_cfg_real_freqs((k - 1) * 3 + h), &
                        openrsp_cfg_real_freqs((k - 1) * 3 + h) /),  &
                        dm_1d = egf_nm, dm_2d = eggf_nm, &
                        dm_3d = egggf_nm, fcc = eggg_nm, fcq = egggg_nm, &
                        rst_order = openrsp_cfg_general_pv_total_anh)

             else

                ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0) * &
                        openrsp_cfg_real_freqs((k - 1) * 3 + h), &
                        openrsp_cfg_real_freqs((k - 1) * 3 + h) /),  &
                        dm_1d = egf_nm, dm_2d = eggf_nm, &
                        dm_3d = egggf_nm, fcc = eggg_nm, fcq = egggg_nm, &
                        rst_elec=openrsp_cfg_general_pv_el_anh, &
                        rst_mech=openrsp_cfg_general_pv_mech_anh)

             end if

             ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0) * &
                              openrsp_cfg_real_freqs((k - 1) * 3 + h), &
                              openrsp_cfg_real_freqs((k - 1) * 3 + h) /), dm_1d = egf_nm) 


             if (k * h == 1) then

                open(unit=259, file='alpha_pv', status='replace', action='write') 

             else

                open(unit=259, file='alpha_pv', status='old', action='write', position='append') 

             end if

             write(259,*) 'Pure vibrational output'
             write(259,*) '======================='
             write(259,*) ' '
             if (openrsp_cfg_general_pv_total_anh > 0) then
                write(259,*) 'Total order of anharmonicity (mechanical and electrical):', &
                openrsp_cfg_general_pv_total_anh
             else
                write(259,*) 'Order of electrical anharmonicity:', &
                openrsp_cfg_general_pv_el_anh
                write(259,*) 'Order of mechanical anharmonicity:', &
                openrsp_cfg_general_pv_mech_anh
             end if
             write(259,*) ' '
             write(259, *) 'Dipole moment (Debye): ', dm_orig*0.393456
             write(259, *) 'Conversion factor: 0.393456 (http://www.theochem.ru.nl/units.html)'
             write(259, *) 'Dipole moment length (normalization factor) (Debye):', &
                            ((sum(dm_orig * dm_orig))**0.5)*0.393456
             write(259, *) ' '
             write(259,*) 'Frequency combination', (k - 1) * 3 + h
             write(259,*) ' '
             write(259,*) 'Frequencies w0, w1 are (a.u.)', (-1.0)* &
                           openrsp_cfg_real_freqs((k - 1) * 3 + h), ' ,', &
                           openrsp_cfg_real_freqs((k - 1) * 3 + h)
             write(259,*) 'Wavelengths for w0, w1 are (nm)', (-1.0)* &
                           aunm / openrsp_cfg_real_freqs((k - 1) * 3 + h), ' ,', &
                           aunm / openrsp_cfg_real_freqs((k - 1) * 3 + h)
             write(259,*) ' '
             write(259,*) 'PV contribution to polarizability'
             write(259,*) '================================='
             write(259,*) ' '

             format_line = '(      f20.8)'
             write(format_line(2:7), '(i6)') size(ff_pv,1)

             do i = 1, 3
                write(259, format_line) real(ff_pv(i,:))
             end do
             write(259,*) ' '
             write(259,*) 'Isotropic (a.u.):', real(((1.0)/(3.0)) * &
                                        (ff_pv(1,1) + ff_pv(2,2) + ff_pv(3,3)))
             write(259,*) 'Isotropic (10^-25 esu):', real(((1.0)/(3.0)) * &
                                        (ff_pv(1,1) + ff_pv(2,2) + ff_pv(3,3))) * 1.481847
             ! Follows method by AJT
             write(259,*) 'Dipole^2 (a.u.):', real(sum( (/ ((ff_pv(i,j) * dm(i) * dm(j), &
                                      i = 1, 3), j = 1, 3) /) ))
             write(259,*) 'Dipole^2 (10^-25 esu):', real(sum( (/ ((ff_pv(i,j) * dm(i) * dm(j), &
                                      i = 1, 3), j = 1, 3) /) )) * 1.481847
             write(259,*) ' '

             close(259)

          end do

       end do


       ! Calculate gradient of polarizability

       kn = (/0,2/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(3))
       allocate(perturbation_tuple%plab(3))
       allocate(perturbation_tuple%pid(3))
       allocate(perturbation_tuple%freq(3))

       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3/)

       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egff')

       ! Read polarizability gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Egff', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          do j = 1, 3
             read(258,*) fld_dum
             egff_cart(i, j, :) = (-1.0) * fld_dum
          end do
       end do

       close(258)

       egff_nm = trans_cartnc_2w1d(3*num_atoms, n_nm, egff_cart, T(:,1:n_nm))

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


       do k = 1, openrsp_cfg_nr_freq_tuples

          do h = 1, 3

             ! Calculate PV contribution to 1st hyperpolarizability

             if (h == 1) then

                if (openrsp_cfg_general_pv_total_anh > 0) then

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 1), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_order = openrsp_cfg_general_pv_total_anh)

                else

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 1), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_elec=openrsp_cfg_general_pv_el_anh, &
                            rst_mech=openrsp_cfg_general_pv_mech_anh)

                end if

             elseif (h == 2) then

                if (openrsp_cfg_general_pv_total_anh > 0) then

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 1), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_order = openrsp_cfg_general_pv_total_anh)

                else

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 1), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_elec=openrsp_cfg_general_pv_el_anh, &
                            rst_mech=openrsp_cfg_general_pv_mech_anh)

                end if

             elseif (h == 3) then

                if (openrsp_cfg_general_pv_total_anh > 0) then

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 2) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_order = openrsp_cfg_general_pv_total_anh)

                else

                   fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0) * &
                            (openrsp_cfg_real_freqs((k - 1) * 3 + 2) + &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3)), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 2), &
                            openrsp_cfg_real_freqs((k - 1) * 3 + 3) /), &
                            dm_1d = egf_nm, dm_2d = eggf_nm, &
                            dm_3d = egggf_nm, po_1d = egff_nm, &
                            po_2d = eggff_nm, po_3d = egggff_nm, &
                            fcc = eggg_nm, fcq = egggg_nm, &
                            rst_elec=openrsp_cfg_general_pv_el_anh, &
                            rst_mech=openrsp_cfg_general_pv_mech_anh)

                end if

             end if

             if (k * h == 1) then

                open(unit=259, file='beta_pv', status='replace', action='write') 

             else

                open(unit=259, file='beta_pv', status='old', action='write', position='append') 

             end if

             write(259,*) 'Pure vibrational output'
             write(259,*) '======================='
             write(259,*) ' '
             if (openrsp_cfg_general_pv_total_anh > 0) then
                write(259,*) 'Total order of anharmonicity (mechanical and electrical):', &
                openrsp_cfg_general_pv_total_anh
             else
                write(259,*) 'Order of electrical anharmonicity:', &
                openrsp_cfg_general_pv_el_anh
                write(259,*) 'Order of mechanical anharmonicity:', &
                openrsp_cfg_general_pv_mech_anh
             end if
             write(259,*) ' '
             write(259, *) 'Dipole moment (Debye): ', dm_orig*0.393456
             write(259, *) 'Conversion factor: 0.393456 (http://www.theochem.ru.nl/units.html)'
             write(259, *) 'Dipole moment length (normalization factor) (Debye):', &
                            ((sum(dm_orig * dm_orig))**0.5)*0.393456
             write(259, *) ' '
             write(259,*) 'Frequency combination', (k - 1) * 3 + h
             write(259,*) ' '

             if (h == 1) then

                write(259,*) 'Frequencies w0, w1, w2 are (a.u.)', (-1.0) * &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 2)), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 1), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 2)
                write(259,*) 'Wavelengths for w0, w1, w2 are (nm)', (-1.0) * aunm / &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 2)), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 1), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 2)

             elseif (h == 2) then

                write(259,*) 'Frequencies w0, w1, w2 are (a.u.)', (-1.0) * &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 1), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)
                write(259,*) 'Wavelengths for w0, w1, w2 are (nm)', (-1.0) * aunm / &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 1) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 1), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 3)


             elseif (h == 3) then

                write(259,*) 'Frequencies w0, w1, w2 are (a.u.)', (-1.0) * &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 2) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 2), ' ,', &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)
                write(259,*) 'Wavelengths for w0, w1, w2 are (nm)', (-1.0) * aunm / &
                             (openrsp_cfg_real_freqs((k - 1) * 3 + 2) + &
                             openrsp_cfg_real_freqs((k - 1) * 3 + 3)), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 2), ' ,', &
                             aunm / openrsp_cfg_real_freqs((k - 1) * 3 + 3)

             end if


             write(259,*) ' '
             write(259,*) 'PV contribution to 1st hyperpolarizability'
             write(259,*) '=========================================='
             write(259,*) ' '

             format_line = '(      f20.8)'
             write(format_line(2:7), '(i6)') size(ff_pv,1)

             do i = 1, 3
                do j = 1, 3
                   write(259, format_line) real(fff_pv(i,j,:))
                end do
                write(259,*) ' '
             end do
             ! Follows method by AJT
             write(259,*) 'Isotropic (a.u.):', real(((1.0)/(5.0)) * sum ( (/ (((fff_pv(i,i,j) + &
                                        fff_pv(i,j,i) + fff_pv(j,i,i)) * dm(j), &
                                        i = 1, 3), j = 1, 3) /)))
             write(259,*) 'Isotropic (10^-30 esu):', real(((1.0)/(5.0)) * sum ( (/ (((fff_pv(i,i,j) + &
                                        fff_pv(i,j,i) + fff_pv(j,i,i)) * dm(j), &
                                        i = 1, 3), j = 1, 3) /)))*0.00863922
             ! Follows method by AJT
             write(259,*) 'Dipole^3 (a.u.):', real(sum( (/ (((fff_pv(i,j,k) * dm(i) * dm(j) * dm(k), &
                                      i = 1, 3), j = 1, 3), k = 1, 3) /) ))
             write(259,*) 'Dipole^3 (10^-30 esu):', real(sum( (/ (((fff_pv(i,j,k) * dm(i) * dm(j) * dm(k), &
                                      i = 1, 3), j = 1, 3), k = 1, 3) /) ))*0.00863922
             write(259,*) ' '

             close(259)

          end do

       end do


       ! Calculate gradient of first hyperpolarizability

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Egfff')

       ! Read 1st hyperpolarizability gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Egfff', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          do j = 1, 3
             do k = 1, 3
                read(258,*) fld_dum
                egfff_cart(i, j, k, :) = fld_dum
             end do
          end do
       end do

       close(258)

       egfff_nm = trans_cartnc_3w1d(3*num_atoms, n_nm, egfff_cart, T(:,1:n_nm))

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


       ! Calculate PV contribution to 2nd hyperpolarizability


       do k = 1, openrsp_cfg_nr_freq_tuples

          if (openrsp_cfg_general_pv_total_anh > 0) then

             ffff_pv = gamma_pv(n_nm, nm_freq, (/ (-1.0) * &
                       sum(openrsp_cfg_real_freqs((k-1)*3 + 1: (k-1)*3 + 3)), &
                       openrsp_cfg_real_freqs((k-1)*3 + 1), &
                       openrsp_cfg_real_freqs((k-1)*3 + 2), &
                       openrsp_cfg_real_freqs((k-1)*3 + 3)/), &
                       dm_1d = egf_nm, dm_2d = eggf_nm, &
                       dm_3d = egggf_nm, po_1d = egff_nm, &
                       po_2d = eggff_nm, po_3d = egggff_nm, &
                       hp_1d = egfff_nm, hp_2d = eggfff_nm, &
                       hp_3d = egggfff_nm, fcc = eggg_nm, fcq = egggg_nm, &
                       rst_order = openrsp_cfg_general_pv_total_anh)

          else

             ffff_pv = gamma_pv(n_nm, nm_freq, (/ (-1.0) * &
                       sum(openrsp_cfg_real_freqs((k-1)*3 + 1: (k-1)*3 + 3)), &
                       openrsp_cfg_real_freqs((k-1)*3 + 1), &
                       openrsp_cfg_real_freqs((k-1)*3 + 2), &
                       openrsp_cfg_real_freqs((k-1)*3 + 3)/), &
                       dm_1d = egf_nm, dm_2d = eggf_nm, &
                       dm_3d = egggf_nm, po_1d = egff_nm, &
                       po_2d = eggff_nm, po_3d = egggff_nm, &
                       hp_1d = egfff_nm, hp_2d = eggfff_nm, &
                       hp_3d = egggfff_nm, fcc = eggg_nm, fcq = egggg_nm, &
                       rst_elec=openrsp_cfg_general_pv_el_anh, &
                       rst_mech=openrsp_cfg_general_pv_mech_anh)

          end if

          if (k == 1) then

             open(unit=259, file='gamma_pv', status='replace', action='write') 

          else

             open(unit=259, file='gamma_pv', status='old', action='write', position='append') 

          end if

          write(259,*) 'Pure vibrational output'
          write(259,*) '======================='
          write(259,*) ' '
          if (openrsp_cfg_general_pv_total_anh > 0) then
             write(259,*) 'Total order of anharmonicity (mechanical and electrical):', &
             openrsp_cfg_general_pv_total_anh
          else
             write(259,*) 'Order of electrical anharmonicity:', &
             openrsp_cfg_general_pv_el_anh
             write(259,*) 'Order of mechanical anharmonicity:', &
             openrsp_cfg_general_pv_mech_anh
          end if
          write(259,*) ' '

          write(259, *) 'Dipole moment (Debye): ', dm_orig*0.393456
          write(259, *) 'Conversion factor: 0.393456 (http://www.theochem.ru.nl/units.html)'
          write(259, *) 'Dipole moment length (normalization factor) (Debye):', &
                         ((sum(dm_orig * dm_orig))**0.5)*0.393456
          write(259, *) ' '
          write(259,*) 'Frequency combination', k
          write(259,*) ' '
          write(259,*) 'Frequencies w0, w1, w2, w3 are (a.u.)', (-1.0) * &
                       sum(openrsp_cfg_real_freqs((k-1)*3 + 1: (k-1)*3 + 3)), ' ,', &
                       openrsp_cfg_real_freqs((k-1)*3 + 1), ' ,', &
                       openrsp_cfg_real_freqs((k-1)*3 + 2), ' ,', &
                       openrsp_cfg_real_freqs((k-1)*3 + 3)
          write(259,*) 'Wavelengths for w0, w1, w2, w3 are (nm)', (-1.0) * &
                       aunm / sum(openrsp_cfg_real_freqs((k-1)*3 + 1: (k-1)*3 + 3)), ' ,', &
                       aunm / openrsp_cfg_real_freqs((k-1)*3 + 1), ' ,', &
                       aunm / openrsp_cfg_real_freqs((k-1)*3 + 2), ' ,', &
                       aunm / openrsp_cfg_real_freqs((k-1)*3 + 3)
          write(259,*) ' '
          write(259,*) 'PV contribution to 2nd hyperpolarizability'
          write(259,*) '=========================================='
          write(259,*) ' '
          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') size(ff_pv,1)

          do i = 1, 3
             do j = 1, 3
                do m = 1, 3
                   write(259, format_line) real(ffff_pv(i,j,m,:))
                end do 
                write(259,*) ' '
             end do
             write(259,*) ' '
          end do
          ! Follows method by AJT
          write(259,*) 'Isotropic (a.u.):', real(((1.0)/(15.0)) * sum ((/ ((ffff_pv(i,i,j,j) + &
                       ffff_pv(i,j,i,j) + ffff_pv(i,j,j,i), i = 1, 3), j = 1, 3) /)))
          ! Follows method by AJT
          write(259,*) 'Dipole^4 (a.u.):', real(sum( (/ ((((ffff_pv(i,j,k,m) * dm(i) * dm(j) * dm(k) * &
                                   dm(m), i = 1, 3), j = 1, 3), k = 1, 3), m = 1, 3) /) ))
          write(259,*) ' '

          close(259)

       end do

       deallocate(T)
       deallocate(nm_freq)
       deallocate(ff_pv)
       deallocate(fff_pv)
       deallocate(ffff_pv)
       deallocate(egf_cart)
       deallocate(egf_nm)
       deallocate(egff_cart)
       deallocate(egff_nm)
       deallocate(egfff_cart)
       deallocate(egfff_nm)

end subroutine
