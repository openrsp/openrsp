!!  openrsp: interface for solving response equations using atomic basis
!!  Copyright 2009 Bin Gao
!!
!!  This file is part of openrsp.
!!
!!  openrsp is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  openrsp is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with openrsp. If not, see <http://www.gnu.org/licenses/>.
!!
!!  openrsp is an interface for solving response equations using
!!  Andreas J. Thorvaldsen's codes, and the solver in DALTON
!!
!!  NOTICE: DALTON is distributed under its own License Agreement
!!          You need to first obtain such agreement and DALTON
!!          from http://www.kjemi.uio.no/software/dalton/dalton.html
!!
!!  NOTICE: Andreas J. Thorvaldsen's codes are also distributed under
!!          GNU Lesser General Public License, you may also first
!!          obtain his code
!!
!!  This file is the main module of openrsp, provides external user interfaces,
!!  also contains the interface of calling Andreas J. Thorvaldsen's codes
!!
!!  2009-12-08, Bin Gao:
!!  * first version

!> \brief main module of openrsp
!> \details contains the external user interfaces,
!>          and the interface of calling Andreas J. Thorvaldsen's codes
!> \author Bin Gao
!> \date 2009-12-08
module openrsp

  use interface_molecule, only: get_nr_atoms, interface_molecule_init
  use interface_io, only: get_print_unit, interface_io_init
  use interface_xc, only: is_ks_calculation, interface_xc_init, xcint_wakeup_workers
  use interface_scf, only: interface_scf_init,  &
                           interface_scf_get_s, &
                           interface_scf_get_d, &
                           interface_scf_get_g, &
                           interface_scf_get_h1
  use interface_f77_memory, only: interface_f77_memory_init, interface_f77_memory_finalize
  use interface_rsp_solver, only: rsp_mosolver_init, rsp_mosolver_finalize
  use interface_1el, only: interface_1el_init
  use interface_basis, only: get_nr_ao, interface_basis_init
  use interface_pelib, only: pe_add_full_operator
  use pe_variables, only: peqm
  use rsp_field_tuple, only: p_tuple, p_tuple_standardorder
  use rsp_general, only: rsp_prop, openrsp_get_property_2014
  use openrsp_cfg
  use matrix_defop, matrix => openrsp_matrix
  use matrix_lowlevel,  only: mat_init
  use eri_contractions, only: ctr_arg
  use eri_basis_loops,  only: unopt_geodiff_loop
  use legacy_properties, only: magnetizability, vcd_aat
  use legacy_vibrational_properties, only: load_vib_modes
  use vib_pv_contribs
  use rsp_sdf_caching, only: SDF, sdf_setup_datatype
  use rsp_indices_and_addressing, only: mat_init_like_and_zero
  use iso_c_binding
  use xcint_fortran_interface

  implicit none

  public openrsp_setup
  public openrsp_calc
  public openrsp_finalize

  ! NEW 2014
  
  public openrsp_calculation_setup
  public write_rsp_tensor
  
  ! END NEW 2014
  
  private

  ! ------------- SCF state and settings ------------
  ! overlap matrix
  type(matrix) S
  ! AO density matrix
  type(matrix), target :: D
  ! Fock matrix
  type(matrix) F

  ! --------------- geometry and basis settings -----------
  ! number of atoms
  integer num_atoms
  ! number of cgto shell blocks
  integer num_cgto_blocks
  ! number of exponents and contraction coefficients
  integer num_exp_and_ctr
  ! array to hold exponents and contraction coefs
  real(8), allocatable :: exp_and_ctr(:)

  ! print level to be used here and in solver
  integer :: print_level = 0

contains

! NEW 2014

  subroutine openrsp_calculation_setup(num_pert, pert_dims, pert_first_comp, pert_labels, &
             pert_freqs, scf_routine, rsp_solver_routine, nucpot_routine, oneel_routine, &
             twoel_routine, xc_routine, prop_size, rsp_tensor, file_id)

    implicit none
    
    integer :: num_pert, id_outp, prop_size
    integer, allocatable, dimension(:) :: pert_labels, pert_dims, pert_first_comp, pert_names
    integer , dimension(2) :: kn
    character, dimension(20) :: file_id
    complex(8), allocatable, dimension(:) :: pert_freqs, rsp_tensor
    external :: scf_routine, rsp_solver_routine
    external :: nucpot_routine, oneel_routine, twoel_routine, xc_routine
    external :: mat_dim, max_mem_mat, unused_iounit, use_disk_for_mat, save_sdf_end
    
    call openrsp_get_property_2014(num_pert, pert_names, pert_labels, pert_dims, pert_freqs, &
                              kn, scf_routine, rsp_solver_routine, nucpot_routine, oneel_routine, &
                              twoel_routine, xc_routine, id_outp, prop_size, rsp_tensor, file_id)
                              
    call write_rsp_tensor(prop_size, rsp_tensor)
    
    
  end subroutine
  
  
  subroutine write_rsp_tensor(prop_size, rsp_tensor)
  
    implicit none
   
    integer :: prop_size
    complex(8), dimension(prop_size) :: rsp_tensor

    
  end subroutine



! END NEW 2014

  subroutine openrsp_setup(LWORK, WORK)
    integer, intent(in)    :: LWORK
    integer                :: nbast, lupri
    real(8), intent(inout) :: WORK(LWORK)
    type(matrix)           :: H1 !one electron Hamiltonian
    type(matrix)           :: G  !two electron Hamiltonian
    real(c_double), allocatable   :: xc_dmat(:)
    real(c_double), allocatable   :: xc_mat(:)
    integer                :: mat_dim
    real(c_double)         :: xc_energy(1)
    real(c_double)         :: num_electrons
    real(8), target        :: temp(1)
    real(8)                :: ave(100)
    real(8), allocatable   :: pe_dmat(:,:)
    real(8), allocatable   :: pe_fmat(:,:)
    real(8)                :: pe_energy(1)
    integer :: ierr

    type(ctr_arg) :: arg(1)

    call interface_molecule_init()
    call interface_1el_init()
    call interface_io_init()
    call interface_xc_init()
    call interface_scf_init()
    call interface_basis_init()

    nbast = get_nr_ao()
    lupri = get_print_unit()

    ! prints the header and license information
    call TITLER('OpenRSP: Response functions computed using AO basis', '*', -1)
    write (LUPRI,*) '>> -- solving equations with the MO response solver in DALTON.'
    write (LUPRI,*)
    write (LUPRI,*) '>> Original reference:'
    write (LUPRI,*) '>>  Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,'// &
                        ' Poul Jorgensen, and Sonia Coriani, '
    write (LUPRI,*) '>>  J. Chem. Phys. 129, 214108 (2008).'
    write (LUPRI,*)

    call interface_f77_memory_init(work_len=lwork, work=work)

    ! initialize the MO response solver
    call rsp_mosolver_init(openrsp_cfg_solver_maxitr, &
                           openrsp_cfg_solver_maxphp, &
                           openrsp_cfg_solver_maxred, &
                           openrsp_cfg_solver_thresh, &
                           openrsp_cfg_solver_optorb, &
                           openrsp_cfg_cpp_damping,   &
                           openrsp_cfg_cpp_used)

    ! initialize and allocate overlap matrix
    call mat_init(S, NBAST, NBAST)

    ! get the overlap
    call interface_scf_get_s(S)

    ! get the AO density matrix and divide by 2
    D = 0*S
    call mat_ensure_alloc(D, only_alloc=.true.)
    call interface_scf_get_d(D)
    D = 0.5d0*D

    ! get the one electron Hamiltonian
    H1 = 0*S
    call mat_ensure_alloc(H1, only_alloc=.true.)
    call interface_scf_get_h1(H1)

    ! get the two electron contribution (G) to Fock matrix
    G = 0*S
    call mat_ensure_alloc(G, only_alloc=.true.)
    call interface_scf_get_g(D, G)

    ! Fock matrix F = H1 + G
    F = H1 + G

    H1 = 0
    G  = 0

    if (is_ks_calculation()) then

       ! add xc contribution to the fock matrix
       mat_dim = D%nrow
       allocate(xc_dmat(mat_dim*mat_dim))
       xc_dmat = 0.0d0
       call daxpy(mat_dim*mat_dim, 1.0d0, D%elms, 1, xc_dmat, 1)
       allocate(xc_mat(mat_dim*mat_dim))
       call xcint_wakeup_workers()
       call xcint_integrate(XCINT_MODE_RKS, &
                            0,              &
                            (/0/),          &
                            (/0/),          &
                            1,              &
                            (/0/),          &
                            (/1/),          &
                            xc_dmat,        &
                            0,              &
                            xc_energy,      &
                            1,              &
                            xc_mat,         &
                            num_electrons)
       call daxpy(mat_dim*mat_dim, 1.0d0, xc_mat, 1, F%elms, 1)
       deallocate(xc_dmat)
       deallocate(xc_mat)
    end if

    if (peqm) then
        mat_dim = D%nrow
        allocate(pe_fmat(mat_dim,mat_dim))
        allocate(pe_dmat(mat_dim,mat_dim))
        pe_fmat = 0.0d0
        pe_dmat = 0.0d0
        call daxpy(mat_dim*mat_dim, 2.0d0, D%elms, 1, pe_dmat, 1)
        call pe_add_full_operator(pe_dmat, pe_fmat, pe_energy)
        call daxpy(mat_dim*mat_dim, 1.0d0, pe_fmat, 1, F%elms, 1)
        deallocate(pe_fmat, pe_dmat)
    end if

    num_atoms = get_nr_atoms()

end subroutine

  subroutine openrsp_calc()

    logical       :: reduced_nc
    integer       :: kn(2)
    integer       :: i, setup_i
    type(p_tuple) :: perturbation_tuple
    type(matrix) :: zeromat_already
    type(SDF), pointer :: F_already, D_already, S_already
    real(8), dimension(3) :: fld_dum
    real(8), dimension(num_atoms) :: masses
    real(8), dimension(3*num_atoms) :: geo_dum
    real(8), dimension(3) :: dm
    real(8), dimension(3) :: dm_orig
    character(len=4) :: file_id
    character(20) :: format_line
    integer       :: n_nm, h, j, k, m, n, p, q, ierr
    real(8), allocatable, dimension(:) :: nm_freq, nm_freq_b
    real(8), allocatable, dimension(:) :: hr_intensities_hv, hr_intensities_vv
    real(8), allocatable, dimension(:,:) :: T
    real(8), dimension(3,3) :: ff
    real(8), dimension(3,3,3) :: fff
    real(8) :: c0, reccm, aunm, xfam
    complex(8), dimension(3) :: CID
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
    complex(8), dimension(3,3,3,3) :: Effff, Effmfww, Effmfw2w
    complex(8), dimension(3,3,6,3) :: Effqfww, Effqfw2w


    ! Source: http://www.webqc.org/unitconverters-js.php
    aunm = 45.56335
    reccm = 219474.6313705

    if (openrsp_cfg_magnetizability) then
       call magnetizability(S, D, F)
    end if

    if (openrsp_cfg_vcd) then
       call vcd_aat(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_general_specify) then

       kn = openrsp_cfg_specify_kn

       perturbation_tuple%n_perturbations = openrsp_cfg_specify_order
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%plab = openrsp_cfg_specify_plab

       do i = 1, openrsp_cfg_specify_order

          if (openrsp_cfg_specify_plab(i) == 'GEO ') then

             perturbation_tuple%pdim(i) = 3*num_atoms

          elseif (openrsp_cfg_specify_plab(i) == 'EL  ') then

             perturbation_tuple%pdim(i) = 3

          elseif (openrsp_cfg_specify_plab(i) == 'MAG ') then

             perturbation_tuple%pdim(i) = 3


          elseif (openrsp_cfg_specify_plab(i) == 'MAG0') then

             perturbation_tuple%pdim(i) = 3

          elseif (openrsp_cfg_specify_plab(i) == 'ELGR') then

             perturbation_tuple%pdim(i) = 6


          else

             write(*,*) 'ERROR: Unrecognized field', openrsp_cfg_specify_plab(i)
             write(*,*) 'encountered when setting up customized response property calculation'

          end if

       end do

       perturbation_tuple%pid = (/(i, i = 1, openrsp_cfg_specify_order)/)

       ! Requires that openrsp_cfg_nr_freq_tuples is set to 1 if no frequencies specified
       ! That is currently the default
       ! Requires that openrsp_cfg_nr_real_freqs is used as freqs. per tuple
       ! That is also currently the default

       do k = 1, openrsp_cfg_nr_freq_tuples

          if(allocated(openrsp_cfg_real_freqs)) then

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order - &
             openrsp_cfg_nr_real_freqs - 1), (-1.0) * &
             sum(openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs)), &
                 openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs) /)

          else

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order)/)

          end if

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/(i, i = 1, openrsp_cfg_specify_order)/)

file_id = '    '


          write(file_id, '(I4)') k


          if (openrsp_cfg_nr_freq_tuples == 1) then

             call rsp_prop(perturbation_tuple, kn, F_unpert=F, D_unpert=D, S_unpert=S)

          else

             if (k == 1) then

                ! ASSUME CLOSED SHELL
                call mat_init(zeromat_already, S%nrow, S%ncol, is_zero=.true.)
                call mat_init_like_and_zero(S, zeromat_already)


                call sdf_setup_datatype(S_already, S)
                call sdf_setup_datatype(D_already, D)
                call sdf_setup_datatype(F_already, F)

             end if

             call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, &
                           file_id = 'freq_tuple_' // trim(adjustl(file_id)))

          end if

       end do


    end if

if (openrsp_cfg_general_trans_cartnc) then


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

       ! Normalize dipole moment - follows procedure by AJT

       dm_orig = dm

       if ( ((sum(dm * dm))**0.5) < 0.00001 ) then

          dm = (/0d0, 0d0, 1d0/)

       else

          dm = dm/((sum(dm * dm))**0.5)

       end if

       ! Get vibrational information

       allocate(T(3*num_atoms, 3*num_atoms))
       allocate(nm_freq_b(3*num_atoms))

       nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

       nm_freq = 0.0

       nm_freq = nm_freq_b(1:n_nm)
       deallocate(nm_freq_b)

       reduced_nc = .TRUE.

       if (reduced_nc) then


          open(unit = 259, file='nm_freq', status='replace', action='write', iostat=ierr)

          do i = 1, n_nm
                      write(259,*) nm_freq(i)* reccm
          end do

          close(259)

          do i = 1, n_nm

             write(*,*) 'Mode:', i, ', frequency in cm^-1', nm_freq(i)* reccm
             T(:, i) = T(:, i) / ((nm_freq(i)**0.5))

          end do

          T = T !* (reccm)**0.5

       end if

       if ((openrsp_cfg_general_cartnc_order_field == 4) .AND. &
           (openrsp_cfg_general_cartnc_order_geo == 3)) then

          allocate(egggffff_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3, 3, 3))
          allocate(egggffff_nm(n_nm, n_nm, n_nm, 3, 3, 3, 3))


          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3*num_atoms
                   do m = 1, 3
                      do n = 1, 3
                         do p = 1, 3
                            read(258,*) fld_dum
                            egggffff_cart(i, j, k, m, n, p, :) = fld_dum
                         end do
                      end do
                   end do
                end do
             end do
          end do

          close(258)

          egggffff_nm = trans_cartnc_4w3d(3*num_atoms, n_nm, egggffff_cart, T)

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3


          do i = 1, n_nm
             do j = i, n_nm
                do k = j, n_nm
                   do m = 1, 3
                      do n = 1, 3
                         do p = 1, 3
                            write(*,format_line) real(egggffff_nm(i, j, k, m, n, p, :))
                         end do
                         write(*,*) ' '
                      end do
                      write(*,*) ' '
                   end do
                   write(*,*) ' '
                end do
                write(*,*) ' '
             end do
             write(*,*) ' '
          end do

          write(*,*) 'Tensor output:'
          write(*,*) ' '

           do n = 1, n_nm
             do p = n, n_nm
                do q = p, n_nm

                   ! Follows method by AJT
                   write(*,*) n, p, q, real(((1.0)/(15.0)) * &
                              sum   ((/ ((egggffff_nm(n,p,q,i,i,j,j) + &
                                          egggffff_nm(n,p,q,i,j,i,j) + &
                                          egggffff_nm(n,p,q,i,j,j,i), i = 1, 3), j = 1, 3) /)))


!                    write(*,*) 'Normal modes: ', n, p, q
!                    ! Follows method by AJT
!                    write(*,*) 'Isotropic (a.u.):', real(((1.0)/(15.0)) * &
!                               sum   ((/ ((egggffff_nm(n,p,q,i,i,j,j) + &
!                                           egggffff_nm(n,p,q,i,j,i,j) + &
!                                           egggffff_nm(n,p,q,i,j,j,i), i = 1, 3), j = 1, 3) /)))
!                    ! Follows method by AJT
!                    write(*,*) 'Dipole^4 (a.u.):', real(sum( (/ &
!                               ((((egggffff_nm(n,p,q,i,j,k,m) * &
!                               dm(i) * dm(j) * dm(k) * &
!                               dm(m), i = 1, 3), j = 1, 3), k = 1, 3), m = 1, 3) /) ))
!                    write(*,*) ' '
                end do
             end do
          end do

          write(*,*) ' '
          write(*,*) 'End tensor output:'


          deallocate(egggffff_cart)
          deallocate(egggffff_nm)


      else if (openrsp_cfg_general_cartnc_order_geo == 2) then

          allocate(egg_cart(3*num_atoms, 3*num_atoms))
          allocate(egg_nm(n_nm, n_nm))
!           allocate(geo_dum(3*num_atoms))

          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
                   read(258,*) geo_dum
                   egg_cart(i, :) = geo_dum
          end do

          close(258)

          egg_nm = trans_cartnc_0w2d(3*num_atoms, n_nm, egg_cart, T(:,1:n_nm))

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3

do i = 1, n_nm

write(*,*) 'Diagonal', i, ' of transformed Hessian:', real(egg_nm(i,i))*reccm

end do


          open(unit = 259, file='rsp_tensor_nm', status='replace', action='write', iostat=ierr)

          do i = 1, n_nm
                      write(259,*) real(egg_nm(i, :))
          end do

          close(259)

          deallocate(egg_cart)
          deallocate(egg_nm)




       else if (openrsp_cfg_general_cartnc_order_geo == 3) then


          allocate(eggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms))
          allocate(eggg_nm(n_nm, n_nm, n_nm))
!           allocate(geo_dum(3*num_atoms))

          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                   read(258,*) geo_dum
                   eggg_cart(i, j, :) = geo_dum
             end do
          end do

          close(258)


          eggg_nm = trans_cartnc_0w3d(3*num_atoms, n_nm, eggg_cart, T(:,1:n_nm))

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3

          open(unit = 259, file='rsp_tensor_nm', status='replace', action='write', iostat=ierr)

          do i = 1, n_nm
             do j = 1, n_nm
                      write(259,*) real(eggg_nm(i, j, :))*reccm
             end do
             write(259,*) ' '
          end do

          close(259)

          deallocate(eggg_cart)
          deallocate(eggg_nm)




       else if (openrsp_cfg_general_cartnc_order_geo == 4) then

          allocate(egggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, 3*num_atoms))
          allocate(egggg_nm(n_nm, n_nm, n_nm, n_nm))
!           allocate(geo_dum(3*num_atoms))

          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3*num_atoms
                   read(258,*) geo_dum
                   egggg_cart(i, j, k, :) = geo_dum
                end do
             end do
          end do

          close(258)


          egggg_nm = trans_cartnc_0w4d(3*num_atoms, n_nm, egggg_cart, T(:,1:n_nm))

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3

          open(unit = 259, file='rsp_tensor_nm', status='replace', action='write', iostat=ierr)

          do i = 1, n_nm
             do j = 1, n_nm
                do k = 1, n_nm
                      write(259,*) real(egggg_nm(i, j, k, :))*reccm
                end do
                write(259,*) ' '
             end do
             write(259,*) ' '
          end do

          close(259)

          deallocate(egggg_cart)
          deallocate(egggg_nm)



      else if ((openrsp_cfg_general_cartnc_order_field == 0) .AND. &
        (openrsp_cfg_general_cartnc_order_geo == 5)) then

          allocate(eggggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                 3*num_atoms, 3*num_atoms))
          allocate(eggggg_nm(n_nm, n_nm, n_nm, n_nm, n_nm))
!           allocate(geo_dum(3*num_atoms))

          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3*num_atoms
                   do m = 1, 3*num_atoms
                      read(258,*) geo_dum
                      eggggg_cart(i, j, k, m, :) = geo_dum
                   end do
                end do
             end do
          end do

          close(258)


          eggggg_nm = trans_cartnc_0w5d(3*num_atoms, n_nm, eggggg_cart, T(:,1:n_nm))

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3


          do i = 1, n_nm
             do j = i, n_nm
                do k = j, n_nm
                   do m = k, n_nm
                      do n = m, n_nm
                         write(*,*) i, j, k, m, n, real(eggggg_nm(i, j, k, m, n))*reccm
                      end do
                   end do
                end do
             end do
          end do


          deallocate(eggggg_cart)
          deallocate(eggggg_nm)



       else if ((openrsp_cfg_general_cartnc_order_field == 0) .AND. &
        (openrsp_cfg_general_cartnc_order_geo == 6)) then

          allocate(egggggg_cart(3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                 3*num_atoms, 3*num_atoms, 3*num_atoms))
          allocate(egggggg_nm(n_nm, n_nm, n_nm, n_nm, n_nm, n_nm))
!           allocate(geo_dum(3*num_atoms))

          open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

          do i = 1, 3*num_atoms
             do j = 1, 3*num_atoms
                do k = 1, 3*num_atoms
                   do m = 1, 3*num_atoms
                      do n = 1, 3*num_atoms
                         read(258,*) geo_dum
                         egggggg_cart(i, j, k, m, n, :) = geo_dum
                      end do
                   end do
                end do
             end do
          end do

          close(258)


          egggggg_nm = trans_cartnc_0w6d(3*num_atoms, n_nm, egggggg_cart, T(:,1:n_nm))

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') 3


          do i = 1, n_nm
             do j = i, n_nm
                do k = j, n_nm
                   do m = k, n_nm
                      do n = m, n_nm
                         do p = n, n_nm
                            write(*,*) i, j, k, m, n, p, real(egggggg_nm(i, j, k, m, n, p))*reccm
                         end do
                      end do
!                       write(*,*) ' '
                   end do
!                    write(*,*) ' '
                end do
!                 write(*,*) ' '
             end do
!              write(*,*) ' '
          end do


          deallocate(egggggg_cart)
          deallocate(egggggg_nm)
!           deallocate(geo_dum)



       else

          write(*,*) 'Unsupported field/geo combination - nothing was done'
          write(*,*) openrsp_cfg_general_cartnc_order_geo
          write(*,*) openrsp_cfg_general_cartnc_order_field

       end if



    end if

    if (openrsp_cfg_general_shg) then


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

       call rsp_prop(perturbation_tuple, kn, F_unpert=F, D_unpert=D, S_unpert=S, file_id='Ef')

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

       ! Read dipole moment from file

       open(unit = 258, file='rsp_tensor_Ef', status='old', action='read', iostat=ierr)
       read(258,*) fld_dum
          dm = fld_dum
       close(258)

       ! Normalize dipole moment - follows procedure by AJT

       open(unit = 320, file='shg_output', status='replace', action='write', iostat=ierr)
       write(320, *) 'Second harmonic generation output'
       write(320, *) '================================='
       write(320, *) ' '
       write(320, *) 'Dipole moment (Debye): ', dm*0.393456
       write(320, *) 'Conversion factor: 0.393456 (http://www.theochem.ru.nl/units.html)'
       write(320, *) 'Dipole moment length (normalization factor) (Debye):',  &
                     ((sum(dm * dm))**0.5)*0.393456
       write(320, *) ' '


       if ( ((sum(dm * dm))**0.5) < 0.00001 ) then

          dm = (/0d0, 0d0, 1d0/)

       else

          dm = dm/((sum(dm * dm))**0.5)

       end if



       ! Requires that openrsp_cfg_nr_freq_tuples is set to 1 if no frequencies specified
       ! That is currently the default
       ! Requires that openrsp_cfg_nr_real_freqs is used as freqs. per tuple
       ! That is also currently the default




       do k = 1, openrsp_cfg_nr_freq_tuples

          write(320,*) 'Frequency combination', k
          write(320,*) ' '

          if (k == 1) then

             ! ASSUME CLOSED SHELL
             call mat_init(zeromat_already, S%nrow, S%ncol, is_zero=.true.)
             call mat_init_like_and_zero(S, zeromat_already)

             call sdf_setup_datatype(S_already, S)
             call sdf_setup_datatype(D_already, D)
             call sdf_setup_datatype(F_already, F)

          end if


! alpha(-2w; 2w)


          kn = (/0,1/)

          perturbation_tuple%n_perturbations = 2
          allocate(perturbation_tuple%pdim(2))
          allocate(perturbation_tuple%plab(2))
          allocate(perturbation_tuple%pid(2))
          allocate(perturbation_tuple%freq(2))

          perturbation_tuple%plab = (/'EL  ',  'EL  '/)
          perturbation_tuple%pdim = (/3, 3/)
          perturbation_tuple%pid = (/1, 2/)

          if(allocated(openrsp_cfg_real_freqs)) then

             perturbation_tuple%freq = (/ (-1.0) * &
             sum(openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs)), &
             sum(openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs))/)

          else

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order)/)

          end if

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/1, 2/)

          write(file_id, '(I4)'), k


          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, &
                        file_id = 'ff2w_freq_tuple_' // trim(adjustl(file_id)))


          open(unit = 258, file='rsp_tensor_' // 'ff2w_freq_tuple_' // trim(adjustl(file_id)), &
               status='old', action='read', iostat=ierr)


          write(320,*) ' '
          write(320,*) 'Polarizability (-2w;2w)'
          write(320,*) '======================='
          write(320,*) ' '
          write(320,*) 'Frequencies w0, w1 are (a.u.)', perturbation_tuple%freq(1), ' ,', &
                        perturbation_tuple%freq(2)
          write(320,*) 'Wavelengths for w0, w1 are (nm)', aunm/perturbation_tuple%freq(1), ' ,', &
                        aunm/perturbation_tuple%freq(2)
          write(320,*) ' '
          write(320,*) 'Polarizability:'
          write(320,*) ' '


          ff = 0.0

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') size(ff,1)

          do i = 1, 3
                read(258,*) fld_dum
                ff(i, :) = fld_dum
                write(320, format_line) ff(i,:)
          end do


          write(320,*) ' '
          write(320,*) 'Isotropic (a.u.):', real(((1.0)/(3.0)) * &
                                     (ff(1,1) + ff(2,2) + ff(3,3)))
          write(320,*) 'Isotropic (10^-25 esu):', real(((1.0)/(3.0)) * &
                                     (ff(1,1) + ff(2,2) + ff(3,3))) * 1.481847
          ! Follows method by AJT
          write(320,*) 'Dipole^2 (a.u.):', real(sum( (/ ((ff(i,j) * dm(i) * dm(j), &
                                   i = 1, 3), j = 1, 3) /) ))
          write(320,*) 'Dipole^2 (10^-25 esu):', real(sum( (/ ((ff(i,j) * dm(i) * dm(j), &
                                   i = 1, 3), j = 1, 3) /) )) * 1.481847
          write(320,*) ' '



          close(258)

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)

! alpha(-w; w)


          kn = (/0,1/)

          perturbation_tuple%n_perturbations = 2
          allocate(perturbation_tuple%pdim(2))
          allocate(perturbation_tuple%plab(2))
          allocate(perturbation_tuple%pid(2))
          allocate(perturbation_tuple%freq(2))

          perturbation_tuple%plab = (/'EL  ',  'EL  '/)
          perturbation_tuple%pdim = (/3, 3/)
          perturbation_tuple%pid = (/1, 2/)

          if(allocated(openrsp_cfg_real_freqs)) then

             perturbation_tuple%freq = (/ (-1.0) * &
             openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1), &
             openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1)/)

          else

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order)/)

          end if

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/1, 2/)

          write(file_id, '(I4)'), k


          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, &
                        file_id = 'ffw_freq_tuple_' // trim(adjustl(file_id)))


          open(unit = 258, file='rsp_tensor_' // 'ffw_freq_tuple_' // trim(adjustl(file_id)), &
               status='old', action='read', iostat=ierr)


          write(320,*) ' '
          write(320,*) 'Polarizability (-w;w)'
          write(320,*) '======================='
          write(320,*) ' '

          write(320,*) 'Frequencies w0, w1 are (a.u.)', perturbation_tuple%freq(1), ' ,', &
                        perturbation_tuple%freq(2)
          write(320,*) 'Wavelengths for w0, w1 are (nm)', aunm/perturbation_tuple%freq(1), ' ,', &
                        aunm/perturbation_tuple%freq(2)
          write(320,*) ' '
          write(320,*) 'Polarizability:'
          write(320,*) ' '


          ff = 0.0

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') size(ff,1)

          do i = 1, 3
                read(258,*) fld_dum
                ff(i, :) = fld_dum
                write(320, format_line) ff(i,:)
          end do


          write(320,*) ' '
          write(320,*) 'Isotropic (a.u.):', real(((1.0)/(3.0)) * &
                                     (ff(1,1) + ff(2,2) + ff(3,3)))
          write(320,*) 'Isotropic (10^-25 esu):', real(((1.0)/(3.0)) * &
                                     (ff(1,1) + ff(2,2) + ff(3,3))) * 1.481847
          ! Follows method by AJT
          write(320,*) 'Dipole^2 (a.u.):', real(sum( (/ ((ff(i,j) * dm(i) * dm(j), &
                                   i = 1, 3), j = 1, 3) /) ))
          write(320,*) 'Dipole^2 (10^-25 esu):', real(sum( (/ ((ff(i,j) * dm(i) * dm(j), &
                                   i = 1, 3), j = 1, 3) /) )) * 1.481847
          write(320,*) ' '



          close(258)

          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)


! beta(-w; w)


          kn = (/1,1/)

          perturbation_tuple%n_perturbations = 3
          allocate(perturbation_tuple%pdim(3))
          allocate(perturbation_tuple%plab(3))
          allocate(perturbation_tuple%pid(3))
          allocate(perturbation_tuple%freq(3))

          perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  '/)
          perturbation_tuple%pdim = (/3, 3, 3/)
          perturbation_tuple%pid = (/1, 2, 3/)



          if(allocated(openrsp_cfg_real_freqs)) then

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order - &
             openrsp_cfg_nr_real_freqs - 1), (-1.0) * &
             sum(openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs)), &
                 openrsp_cfg_real_freqs((k - 1) * openrsp_cfg_nr_real_freqs + 1: &
                                              k * openrsp_cfg_nr_real_freqs) /)

          else

             perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order)/)

          end if

          perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
          perturbation_tuple%pid = (/1, 2, 3/)

          write(file_id, '(I4)'), k


          call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                        S_already=S_already, zeromat_already=zeromat_already, &
                        file_id = 'fff_freq_tuple_' // trim(adjustl(file_id)))



          write(320,*) ' '
          write(320,*) 'Hyperpolarizability (-2w;w,w)'
          write(320,*) '============================='
          write(320,*) ' '
          write(320,*) 'Frequencies w0, w1, w2 are (a.u.)', perturbation_tuple%freq(1), ' ,', &
                        perturbation_tuple%freq(2), perturbation_tuple%freq(3)
          write(320,*) 'Wavelengths for w0, w1, w2 are (nm)', aunm/perturbation_tuple%freq(1), ' ,', &
                        aunm/perturbation_tuple%freq(2), aunm/perturbation_tuple%freq(3)
          write(320,*) ' '
          write(320,*) 'Hyperpolarizability:'
          write(320,*) ' '

          open(unit = 258, file='rsp_tensor_' // 'fff_freq_tuple_' // trim(adjustl(file_id)), &
               status='old', action='read', iostat=ierr)

          fff = 0.0

          format_line = '(      f20.8)'
          write(format_line(2:7), '(i6)') size(ff,1)

          do i = 1, 3
             do j = 1, 3
                read(258,*) fld_dum
                fff(i, j, :) = fld_dum
                write(320, format_line) fff(i,j,:)
             end do
          end do

          close(258)

          write(320,*) ' '
          ! Follows method by AJT
          write(320,*) 'Isotropic (a.u.):', real((1.0/5.0) * sum ( (/ (( (fff(i,i,j) + &
                       fff(i,j,i) + fff(j,i,i)) * dm(j), i = 1, 3), j = 1, 3) /)))
          write(320,*) 'Isotropic (10^-30 esu):', real((1.0/5.0) * sum ( (/ (( (fff(i,i,j) + &
                       fff(i,j,i) + fff(j,i,i)) * dm(j), i = 1, 3), j = 1, 3) /))) * 0.00863922
          ! Follows method by AJT
          write(320,*) 'Dipole^3 (a.u.):', real(sum( (/ (((fff(i,j,k) * dm(i) * dm(j) * dm(k), &
                                   i = 1, 3), j = 1, 3), k = 1, 3) /) ))
          write(320,*) 'Dipole^3 (10^-30 esu):', real(sum( (/ (((fff(i,j,k) * dm(i) * dm(j) * dm(k), &
                                   i = 1, 3), j = 1, 3), k = 1, 3) /) )) * 0.00863922
          write(320,*) ' '




          deallocate(perturbation_tuple%pdim)
          deallocate(perturbation_tuple%plab)
          deallocate(perturbation_tuple%pid)
          deallocate(perturbation_tuple%freq)


       end do

       close(320)

    end if

    if (openrsp_cfg_general_pv2f) then
       call openrsp_general_pv2f(num_atoms, S, D, F)
    end if

    if (openrsp_cfg_general_pv3f) then
       call openrsp_general_pv3f(num_atoms, S, D, F)
    end if

    if (openrsp_cfg_general_pv4f) then
       call openrsp_general_pv4f(num_atoms, S, D, F)
    end if

    if (openrsp_cfg_general_hyper_raman) then

       ! Verify the value of the conversion factor
       reccm = 219474.6313705d0
       fld_dum = 0.0

       ! Get normal mode transformation matrix
       ! Get normal mode frequencies

       allocate(T(3*num_atoms, 3*num_atoms))

       allocate(nm_freq_b(3*num_atoms))

       nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

       nm_freq = nm_freq_b(1:n_nm) * openrsp_cfg_general_hypram_freqscale

       deallocate(nm_freq_b)

       allocate(egfff_cart(3*num_atoms, 3, 3, 3))
       allocate(egfff_nm(n_nm, 3, 3, 3))
       allocate(hr_intensities_hv(n_nm))
       allocate(hr_intensities_vv(n_nm))

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
       perturbation_tuple%freq = (/0.0d0, -2.0d0* openrsp_cfg_real_freqs(1),&
                       openrsp_cfg_real_freqs(1), openrsp_cfg_real_freqs(1)/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       ! ASSUME CLOSED SHELL
       call mat_init(zeromat_already, S%nrow, S%ncol, is_zero=.true.)
       call mat_init_like_and_zero(S, zeromat_already)


       call sdf_setup_datatype(S_already, S)
       call sdf_setup_datatype(D_already, D)
       call sdf_setup_datatype(F_already, F)



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

       open(unit = 12, file='hr_output', status='replace', action='write', iostat=ierr)

       write(12,*) 'HYPER-RAMAN INTENSITY OUTPUT'
       write(12,*) '============================'
       write(12,*) ' '
       ! Write something about time and date
       ! Adapt frequency to multiple frequency tuple setup
       write(12,100) openrsp_cfg_real_freqs(1)
       100 format('The Hyper-Raman intensities are calculated using an electronic frequency of ', &
       F10.7, ' a.u.')
       write(12,*) ' '

       write(12,101) openrsp_cfg_temperature
       101 format('The temperature is set at', F10.5 ' K.')
       write(12,*) ' '

       write(12,102) n_nm
       102 format('The number of normal modes considered in this system is', I5, ' .')
       write(12,*) 'These modes have the following frequencies:'
       write(12,*) ' '

       write(12,*) 'Mode      Frequency (a.u.)'
       write(12,*) '=========================='
       write(12,*) ' '

       do i = 1, n_nm
          write(12,103) i, nm_freq(i)
          103 format(I6, 4X, F10.7)
       end do


       write(12,*) ' '

       write(12,120) openrsp_cfg_general_hypram_freqscale
       120 format('The frequencies are all scaled by a factor of', F10.5, ' .')
       write(12,*) ' '


       ! Calculate the intensities mode by mode

       ! Adapt frequency to multiple frequency tuple setup
       do i = 1, n_nm
          hr_intensities_hv(i) = hyp_raman_hv(n_nm, openrsp_cfg_temperature, &
                                 openrsp_cfg_real_freqs(1), nm_freq(i), dreal(egfff_nm(i, :, :, :)))
       end do

       write(12,*) 'HR intensities using HV polarization:'
       write(12,*) ' '
       write(12,*) 'Mode        Freq. (a.u.)    Freq (cm^-1)      Intensity '
       write(12,*) '====================================================================='
       write(12,*) ' '

       do i = 1, n_nm
          write(12,114) i, nm_freq(i), nm_freq(i)*reccm, hr_intensities_hv(i)
          114 format(I5, 4X, F13.7, 6X, F13.7, 4X, F25.10)
       end do
       write(12,*) ' '

       ! Adapt frequency to multiple frequency tuple setup
       do i = 1, n_nm
          hr_intensities_vv(i) = hyp_raman_vv(n_nm, openrsp_cfg_temperature, &
                                 openrsp_cfg_real_freqs(1), nm_freq(i), dreal(egfff_nm(i, :, :, :)))
       end do

       write(12,*) 'HR intensities using VV polarization:'
       write(12,*) ' '
       write(12,*) 'Mode        Freq. (a.u.)    Freq (cm^-1)      Intensity '
       write(12,*) '====================================================================='
       write(12,*) ' '

       do i = 1, n_nm
          write(12,115) i, nm_freq(i), nm_freq(i)*reccm, hr_intensities_vv(i)
          115 format(I5, 4X, F13.7, 6X, F13.7, 4X, F25.10)
       end do
       write(12,*) ' '
       close(12)

       deallocate(hr_intensities_hv)
       deallocate(hr_intensities_vv)
       deallocate(nm_freq)
       deallocate(T)

    end if

    if (openrsp_cfg_general_efishgcid) then
       call openrsp_efishgcid(num_atoms, S, D, F)
    end if

  end subroutine

  subroutine openrsp_finalize()
    integer lupri
    lupri = get_print_unit()
    ! free

    S = 0
    D = 0
    F = 0

    call rsp_mosolver_finalize()
    call interface_f77_memory_finalize()
    ! stamp with date, time and hostname
    call TSTAMP(' ', lupri)
    write (lupri,*)
    write (lupri,*) '>> ** End of OpenRSP Section'
    write (lupri,*)
  end subroutine

end module

!> \brief driver used for DALTON
!> \author Bin Gao
!> \date 2009-12-08
!> \param WORK contains the work memory
!> \param LWORK is the size of the work memory
subroutine openrsp_driver(work, lwork)

   use openrsp

   implicit none

   integer, intent(in)    :: lwork
   real(8), intent(inout) :: work(lwork)

   call openrsp_setup(lwork, work)
   call openrsp_calc()
   call openrsp_finalize()

end subroutine
