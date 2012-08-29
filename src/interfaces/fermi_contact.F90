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

   subroutine get_fc_integrals(           &
                               M,         &
                               center,    &
                               component  &
                              )

!     --------------------------------------------------------------------------
      type(matrix)                  :: M
      integer, intent(in)           :: center
      integer, intent(in), optional :: component
!     --------------------------------------------------------------------------
      logical                       :: take_derv
      integer                       :: center_d
      integer                       :: ixyz_d
      real(8), allocatable          :: ao(:, :)
      real(8), allocatable          :: buffer(:, :)
      integer                       :: i, j, iblock
      integer                       :: nr1, nr2, st1, st2
      real(8)                       :: d
!     --------------------------------------------------------------------------

      M%elms_alpha = 0.0d0

!     if you have 4 atoms, then there are 12 components
!     code below figures out the center and direction based on component
      if (present(component)) then
         center_d = (component + 2)/3
         if (mod(component, 3) == 1) ixyz_d = 1
         if (mod(component, 3) == 2) ixyz_d = 2
         if (mod(component, 3) == 0) ixyz_d = 3
         take_derv = .true.
      else
         take_derv = .false.
      end if

      call interface_ao_read(.false.)
      call ao_eval_init(1, 0, .false.)

      allocate(ao(1,     nr_ao_slices*nr_ao_cartesian))
      ao = 0.0d0
      allocate(buffer(1, nr_ao_slices*nr_ao_cartesian))
      buffer = 0.0d0

      call get_ao(1, center_xyz(center, 1), center_xyz(center, 2), center_xyz(center, 3), ao, buffer)

      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(lssl_block_partner(iblock, 0, 0))
         st2 = ao_block_start(lssl_block_partner(iblock, 0, 0))
         do i = st1, st1 + nr1 - 1
            do j = st2, st2 + nr2 - 1
               if (take_derv) then
                  d = 0.0d0
                  if (ao_center(i) == center_d) then
                     if (ao_center(i) /= center) then
                        d = d - ao(1, ao_off_g1_m0(ixyz_d, 0) + i)*ao(1, j)
                     else
                        if (ao_center(i) /= ao_center(j)) then
                           d = d + ao(1, ao_off_g1_m0(ixyz_d, 0) + j)*ao(1, i)
                        end if
                     end if
                  end if
                  if (ao_center(j) == center_d) then
                     if (ao_center(j) /= center) then
                        d = d - ao(1, ao_off_g1_m0(ixyz_d, 0) + j)*ao(1, i)
                     else
                        if (ao_center(i) /= ao_center(j)) then
                           d = d + ao(1, ao_off_g1_m0(ixyz_d, 0) + i)*ao(1, j)
                        end if
                     end if
                  end if
               else
                  d = ao(1, i)*ao(1, j)
               end if
               M%elms_alpha(j, i, 1) = d
            end do
         end do
      end do

      deallocate(ao)
      deallocate(buffer)

   end subroutine
#endif /* ifdef PRG_DIRAC */

end module
