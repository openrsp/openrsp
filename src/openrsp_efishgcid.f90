subroutine openrsp_efishgcid(num_atoms, S, D, F)

   use rsp_field_tuple, only: p_tuple, p_tuple_standardorder
   use matrix_defop, matrix => openrsp_matrix
   use rsp_sdf_caching, only: SDF, sdf_setup_datatype
   use rsp_general, only: rsp_prop
   use openrsp_cfg
   use matrix_lowlevel,  only: mat_init
   use rsp_mag_prop, only: M_efishg, D_efishg, L_efishg
   use rsp_indices_and_addressing, only: mat_init_like_and_zero

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
   integer       :: i, j, k, m, ierr
   real(8) :: c0
   complex(8), dimension(3) :: CID
   complex(8), dimension(3,3,3,3) :: Effff, Effmfww, Effmfw2w
   complex(8), dimension(3,3,6,3) :: Effqfww, Effqfw2w




       ! FIND CORRECT VALUE: AROUND 137?
       c0 = 137.035999074

       fld_dum = 0.0


       Effff = 0.0
       Effmfww = 0.0
       Effmfw2w = 0.0
       Effqfww = 0.0
       Effqfw2w = 0.0

       ! ASSUME CLOSED SHELL
       call mat_init(zeromat_already, S%nrow, S%ncol, is_zero=.true.)
       call mat_init_like_and_zero(S, zeromat_already)

       call sdf_setup_datatype(S_already, S)
       call sdf_setup_datatype(D_already, D)
       call sdf_setup_datatype(F_already, F)

do k = 1, openrsp_cfg_nr_freq_tuples

       ! Calculate Effff(w, w, 0)

       kn = (/1,2/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/-2.0d0 * openrsp_cfg_real_freqs(k), &
       openrsp_cfg_real_freqs(k), openrsp_cfg_real_freqs(k), 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)


       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, file_id='Effff')


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Effff', status='old', action='read', iostat=ierr)

       do i = 1, 3
          do j = 1, 3
             do m = 1, 3
                read(258,*) fld_dum
                Effff(i, :, m, j) = -(1.0) * fld_dum
             end do
          end do
       end do

       close(258)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)



       ! Calculate Effqf(w, w, 0)

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'ELGR', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/6, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/openrsp_cfg_real_freqs(k), &
       openrsp_cfg_real_freqs(k), -2.0d0 * openrsp_cfg_real_freqs(k), 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)


       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, &
                           file_id='Effqfww')


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Effqfww', status='old', action='read', iostat=ierr)

! IMPORTANT: FIND OUT HOW TUPLES ARE SORTED AND READ IN ACCORDINGLY (IT MAY ALREADY BE CORRECT)

! Remaining: (1.0) * fld_dum/2.0, (1.0) * fld_dum!/2.0
       do i = 1, 6
          do j = 1, 3
             do m = 1, 3
                read(258,*) fld_dum
                Effqfww(j, :, i, m) = (1.0) * fld_dum/2.0
             end do
          end do
       end do


       close(258)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

       ! Calculate Effqf(w, -2w, 0)

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'ELGR', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/6, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/-2.0d0 * openrsp_cfg_real_freqs(k), &
       openrsp_cfg_real_freqs(k), openrsp_cfg_real_freqs(k), 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)


       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, &
                           file_id='Effqfw2w')


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Effqfw2w', status='old', action='read', iostat=ierr)

! IMPORTANT: FIND OUT HOW TUPLES ARE SORTED AND READ IN ACCORDINGLY (IT MAY ALREADY BE CORRECT)
       do i = 1, 6
          do j = 1, 3
             do m = 1, 3
                read(258,*) fld_dum
                Effqfw2w(:, m, i, j) = (1.0) * fld_dum/2.0
             end do
          end do
       end do


       close(258)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


! THIS AND Effmf((w,-2w, 0): Find out about freq. ordering in calculation and ordering when reading
! tensors in and factors -1 and 1/2 

       ! Calculate Effmf(w, w, 0)

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'MAG ', 'EL  '/)
       perturbation_tuple%pdim = (/3, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/-2.0d0 * openrsp_cfg_real_freqs(k), &
       openrsp_cfg_real_freqs(k),  openrsp_cfg_real_freqs(k), 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)


       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, &
                           file_id='Emfffww')


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Emfffww', status='old', action='read', iostat=ierr)

! IMPORTANT: FIND OUT HOW TUPLES ARE SORTED AND READ IN ACCORDINGLY
       do i = 1, 3
          do j = 1, 3
             do m = 1, 3
                read(258,*) fld_dum
                Effmfww(j, :, i, m) = (1.0) * fld_dum!/2.0
             end do
          end do
       end do


       close(258)

       
                 write(*,*) 'Isotropic (a.u.):', real(((1.0)/(15.0)) * sum ((/ ((Effmfww(i,i,j,j) + &
                       Effmfww(i,j,i,j) + Effmfww(i,j,j,i), i = 1, 3), j = 1, 3) /)))
       
       
       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

       ! Calculate Effmf(w, -2w, 0)

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'MAG ', 'EL  '/)
       perturbation_tuple%pdim = (/3, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/openrsp_cfg_real_freqs(k), &
       openrsp_cfg_real_freqs(k), -2.0d0 *  openrsp_cfg_real_freqs(k), 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       call rsp_prop(perturbation_tuple, kn, F_already=F_already, D_already=D_already, &
                           S_already=S_already, zeromat_already=zeromat_already, &
                           file_id='Effmfw2w')


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor_Effmfw2w', status='old', action='read', iostat=ierr)

! IMPORTANT: FIND OUT HOW TUPLES ARE SORTED AND READ IN ACCORDINGLY
       do i = 1, 3
          do j = 1, 3
             do m = 1, 3
                read(258,*) fld_dum
                Effmfw2w(:, m, i, j) = (1.0) * fld_dum!/2.0
             end do
          end do
       end do

       close(258)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


do i = 1, 3

write(*,*) 'M', i, ' is', M_efishg(i, Effff, Effmfww, Effmfw2w)
write(*,*) 'L', i, ' is', L_efishg(i, Effff, Effqfww, Effqfw2w)
write(*,*) 'D', i, ' is', D_efishg(i, Effff)

 CID(i) = (1.0/c0) * ( M_efishg(i, Effff, Effmfww, Effmfw2w) + &
          openrsp_cfg_real_freqs(k) *  L_efishg(i, Effff, Effqfww, Effqfw2w) ) / &
          ( D_efishg(i, Effff) )

end do


! Add output printing

write(*,*) ' '
write(*,*) 'CID values'
write(*,*) '=========='
write(*,*) 'CID(1) =', CID(1)
write(*,*) 'CID(2) =', CID(2)
write(*,*) 'CID(3) =', CID(3)
write(*,*) ' '



end do

end subroutine
