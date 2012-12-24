module interface_dirac_gen1int

#ifndef VAR_LSDALTON
   use gen1int_api

   implicit none

   public get_1el_integrals

   private

contains
   subroutine get_1el_integrals(                 &
                                M,               &
                                prop_name,       &
                                num_ints,        &
                                order_mom,       &
                                order_elec,      &
                                order_geo_total, &
                                max_num_cent,    &
                                blocks,          &
                                print_unit       &
                               )

      type(matrix)                :: M(*)
      character(*), intent(in)    :: prop_name
      integer,      intent(in)    :: num_ints
      integer,      intent(in)    :: order_mom
      integer,      intent(in)    :: order_elec
      integer,      intent(in)    :: order_geo_total
      integer,      intent(in)    :: max_num_cent
      integer,      intent(in)    :: blocks(:)
      integer,      intent(in)    :: print_unit

      logical,      parameter     :: write_integrals_to_file = .false.
      integer                     :: i

      do i = 1, num_ints
         M(i)%elms = 0.0d0
      end do

      call gen1int_host_get_int(non_lao, trim(prop_name),      &
                                order_mom,                     &  !multipole moments
                                order_elec,                    &
                                0, 0, 0,                       &  !magnetic derivatives
                                0, 0, 0,                       &  !derivatives w.r.t. total ram
                                0, 0,                          &  !partial geometric derivatives
                                max_num_cent, order_geo_total, &  !total geometric derivatives
                                0, (/0/), redundant_geo,       &  !total geometric derivatives
                                .false., .false., .false.,     &  !not implemented
                                num_ints, M,                   &
                                write_integrals_to_file,       &
                                size(blocks)/2, (/blocks/),    &
                                print_unit, 0)

   end subroutine
#endif
end module
