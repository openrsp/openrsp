module fermi_contact

#ifdef PRG_DIRAC
!  radovan: routine to calculate FC "integrals"
!           and derivatives with respect to nuclear centers
!           factor dfac is missing since these are directly used
!           for E_PNC calculations and not for FC

   use interface_ao
   use ao_eval
   use matrix_defop

   implicit none

   public get_fc_integrals

   private

contains

   subroutine get_fc_integrals(            &
                               M,          &
                               center_pnc, &
                               component1, &
                               component2  &
                              )

!     --------------------------------------------------------------------------
      type(matrix)                  :: M
      integer, intent(in)           :: center_pnc
      integer, intent(in), optional :: component1
      integer, intent(in), optional :: component2
!     --------------------------------------------------------------------------
      integer                       :: ixyz
      integer                       :: jxyz
      real(8), allocatable          :: ao(:, :)
      real(8), allocatable          :: buffer(:, :)
      integer                       :: iblock
      integer                       :: i, j, k, l
      integer                       :: nr1, nr2, st1, st2
      real(8)                       :: d
      integer                       :: order
      integer                       :: center_i
      integer                       :: center_j
      integer                       :: center_d1
      integer                       :: center_d2
      integer                       :: center_k
      integer                       :: center_l
      logical                       :: ij_different_center

      integer, parameter :: pos(9) = (/1, 2, 3, 2, 4, 5, 3, 5, 6/)
!     --------------------------------------------------------------------------

      M%elms_alpha = 0.0d0

!     if you have 4 atoms, then there are 12 components
!     code below figures out the center and direction based on component
      order = 0
      if (present(component1)) then
         center_d1 = (component1 + 2)/3
         if (mod(component1, 3) == 1) ixyz = 1
         if (mod(component1, 3) == 2) ixyz = 2
         if (mod(component1, 3) == 0) ixyz = 3
         order = order + 1
      end if
      if (present(component2)) then
         center_d2 = (component2 + 2)/3
         if (mod(component2, 3) == 1) jxyz = 1
         if (mod(component2, 3) == 2) jxyz = 2
         if (mod(component2, 3) == 0) jxyz = 3
         order = order + 1
      end if

      call interface_ao_read(.false.)
      call ao_eval_init(order, 0, .false.)

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0

      call get_ao(1,                         &
                  center_xyz(center_pnc, 1), &
                  center_xyz(center_pnc, 2), &
                  center_xyz(center_pnc, 3), &
                  ao,                        &
                  buffer)

      select case (order)

         case (0)
            do iblock = 1, nr_ao_blocks
               nr1 = ao_block_nr(iblock)
               st1 = ao_block_start(iblock)
               nr2 = ao_block_nr(lssl_block_partner(iblock, 0, 0))
               st2 = ao_block_start(lssl_block_partner(iblock, 0, 0))
               if (st1 < st2) then
                  do i = st1, st1 + nr1 - 1
                     do j = st2, st2 + nr2 - 1
                        d = ao(1, i)*ao(1, j)
                        M%elms_alpha(j, i, 1) = d
                        M%elms_alpha(i, j, 1) = d
                     end do
                  end do
               end if
            end do

         case (1)
            do iblock = 1, nr_ao_blocks
               nr1 = ao_block_nr(iblock)
               st1 = ao_block_start(iblock)
               nr2 = ao_block_nr(lssl_block_partner(iblock, 0, 0))
               st2 = ao_block_start(lssl_block_partner(iblock, 0, 0))
               if (st1 < st2) then
                  do i = st1, st1 + nr1 - 1
                     center_i = ao_center(i)
                     do j = st2, st2 + nr2 - 1
                        center_j = ao_center(j)
                        ij_different_center = (center_i /= center_j)

                        d = 0.0d0

                        ! i_1 j
                        k = i
                        l = j
                        center_k = center_i
                        if (center_d1 == center_k) then
                           if (center_d1 == center_pnc) then
                              if (ij_different_center) then
                                 d = d + ao(1, ao_off_g1_m0(ixyz, 0) + l)*ao(1, k)
                              end if
                           else
                              d = d - ao(1, ao_off_g1_m0(ixyz, 0) + k)*ao(1, l)
                           end if
                        end if

                        ! j_1 i
                        k = j
                        l = i
                        center_k = center_j
                        if (center_d1 == center_k) then
                           if (center_d1 == center_pnc) then
                              if (ij_different_center) then
                                 d = d + ao(1, ao_off_g1_m0(ixyz, 0) + l)*ao(1, k)
                              end if
                           else
                              d = d - ao(1, ao_off_g1_m0(ixyz, 0) + k)*ao(1, l)
                           end if
                        end if

                        M%elms_alpha(j, i, 1) = d
                        M%elms_alpha(i, j, 1) = d
                     end do
                  end do
               end if
            end do

         case (2)
            do iblock = 1, nr_ao_blocks
               nr1 = ao_block_nr(iblock)
               st1 = ao_block_start(iblock)
               nr2 = ao_block_nr(lssl_block_partner(iblock, 0, 0))
               st2 = ao_block_start(lssl_block_partner(iblock, 0, 0))
               if (st1 < st2) then
                  do i = st1, st1 + nr1 - 1
                     center_i = ao_center(i)
                     do j = st2, st2 + nr2 - 1
                        center_j = ao_center(j)
                        ij_different_center = (center_i /= center_j)

                        d = 0.0d0

                        k = i
                        l = j
                        ! i_12 j
                        center_k = center_i
                        center_l = center_i
                        if (center_d1 == center_k) then
                           if (center_d2 == center_l) then
                              if (center_d1 == center_pnc) then
                                 if (ij_different_center) then
                                    d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + l)*ao(1, k)
                                 end if
                              else
                                 d = d + ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + k)*ao(1, l)
                              end if
                           end if
                        end if
                        ! i_1  j_2
                        center_k = center_i
                        center_l = center_j
                        if (center_d1 == center_k) then
                           if (center_d2 == center_l) then
                            ! if ((center_d1 == center_pnc) .and. (center_d2 == center_pnc)) then
                            ! end if
                              if ((center_d1 == center_pnc) .and. (center_d2 /= center_pnc)) then
                                 d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + l)*ao(1, k)
                              end if
                              if ((center_d1 /= center_pnc) .and. (center_d2 == center_pnc)) then
                                 d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + k)*ao(1, l)
                              end if
                              if ((center_d1 /= center_pnc) .and. (center_d2 /= center_pnc)) then
                                 d = d + ao(1, ao_off_g1_m0(ixyz, 0) + k)*ao(1, ao_off_g1_m0(jxyz, 0) + l)
                              end if
                           end if
                        end if

                        k = j
                        l = i
                        ! j_12 i
                        center_k = center_j
                        center_l = center_j
                        if (center_d1 == center_k) then
                           if (center_d2 == center_l) then
                              if (center_d1 == center_pnc) then
                                 if (ij_different_center) then
                                    d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + l)*ao(1, k)
                                 end if
                              else
                                 d = d + ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + k)*ao(1, l)
                              end if
                           end if
                        end if
                        ! j_1  i_2
                        center_k = center_j
                        center_l = center_i
                        if (center_d1 == center_k) then
                           if (center_d2 == center_l) then
                            ! if ((center_d1 == center_pnc) .and. (center_d2 == center_pnc)) then
                            ! end if
                              if ((center_d1 == center_pnc) .and. (center_d2 /= center_pnc)) then
                                 d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + l)*ao(1, k)
                              end if
                              if ((center_d1 /= center_pnc) .and. (center_d2 == center_pnc)) then
                                 d = d - ao(1, ao_off_g2_m0(pos((ixyz-1)*3 + jxyz), 0) + k)*ao(1, l)
                              end if
                              if ((center_d1 /= center_pnc) .and. (center_d2 /= center_pnc)) then
                                 d = d + ao(1, ao_off_g1_m0(ixyz, 0) + k)*ao(1, ao_off_g1_m0(jxyz, 0) + l)
                              end if
                           end if
                        end if

                        M%elms_alpha(j, i, 1) = d
                        M%elms_alpha(i, j, 1) = d
                     end do
                  end do
               end if
            end do

         case default
            print *, 'error: order too hight in fermi_contact'
            stop 1

      end select

      deallocate(ao)
      deallocate(buffer)

   end subroutine
#endif /* ifdef PRG_DIRAC */

end module
