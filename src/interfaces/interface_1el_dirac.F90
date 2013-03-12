module interface_1el_dirac

   use matrix_defop
   use matrix_lowlevel, only: mat_print
   use interface_molecule, only: get_nr_atoms

   implicit none

#ifdef PRG_DIRAC
   public order_1el_integrals
   public read_1el_integrals
   public print_prp_common_blocks

   private

   integer, parameter :: max_aoproper_indices = 900
   character(16)      :: op_name_to_aoproper_index(max_aoproper_indices) = '----------------'

#include "dcbbas.h"
#include "dcbcls.h"
#include "dcbgen.h"
#include "dcbprl.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"

contains

   subroutine order_1el_integrals()

!     --------------------------------------------------------------------------
      character(1) :: ss
      character(4) :: blocks
      character(8) :: s, s1, s2, s3
      real(8)      :: f
      integer      :: i, ij, iatom
      integer      :: irep(18)
      character(2) ::  i_iter(3) = (/ 'X',  'Y',  'Z'/)
      character(2) :: ij_iter(6) = (/'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'/)
!     --------------------------------------------------------------------------
!     todo: 1. adapt for symmetry
!              only diplen is right at the moment
!           2. generalize for other hamiltonians
!              not every combination hamiltonian--property is tested
!     --------------------------------------------------------------------------

      logical :: need_1el_f = .false.
      logical :: need_1el_g = .false.
      logical :: need_1el_v = .false.
      logical :: need_1el_q = .false.
      logical :: need_1el_b = .false.
      logical :: need_1el_o = .true.

      logical :: diamagnetic_via_pn = .false.
      logical :: nonzero_small_metric = .false.
      logical :: openrsp_cfg_mllsus = .false.


!     decide what to do with the small-small metric
!     =============================================

      ss = '0'
      if (nonzero_small_metric) ss = '+'
      ss = '+' !fixme hardcoded


!     add to dcbcls
!     =============

      if (need_1el_f) then
         call add_class('DIPLEN ', '+00'//ss)
         if (need_1el_g) then
            call add_class('DPLGRA ', '+00'//ss)
         end if
      end if

      if (need_1el_v) then
         call add_class('DIPVEL ', '+00'//ss)
      end if

      if (need_1el_q) then
         call add_class('THETA  ', '+00'//ss)
         if (need_1el_g) then
            call add_class('QUAGRA ', '+00'//ss)
         end if
      end if

      if (need_1el_o .or. need_1el_b) then
         call add_class('DIPLEN ', '0++0')
         if (need_1el_g) then
            call add_class('DPLGRA ', '0++0')
         end if
         call add_class('DSUSCGO', '+00'//ss)
      end if

      if (need_1el_b) then
         call add_class('S1MAG  ', '000+')
         call add_class('S2MAG  ', '000+')
         call add_class('RM1H2  ', '+00'//ss)
         call add_class('RM2H2  ', '+00'//ss)
         call add_class('RM1H3  ', '0++0')
         call add_class('RM2H3  ', '0++0')
         call add_class('RM1RN  ', '0++0')
         call add_class('S1MAG  ', '+00'//ss)
         call add_class('S2MAG  ', '+00'//ss)
         call add_class('S1MAGR ', '+00'//ss)
         call add_class('S2MKET ', '+00'//ss)
         call add_class('S2MMIX ', '+00'//ss)
         if (openrsp_cfg_mllsus) then
            ! reproduces Dalton
            call add_class('RDSUSLL', '00+0')
         else
            ! reproduces DC RKB large c
            call add_class('RDSUSLL', '0++0')
         end if
         call add_class('DSUSNOL', '+00'//ss)
         if (need_1el_f) then
            call add_class('CM1    ', '+00'//ss)
            call add_class('CM2    ', '+00'//ss)
         end if
         if (need_1el_q) then
            call add_class('QDBINT ', '+00'//ss)
         end if
      end if

      if (need_1el_g) then
         call add_class('G1O    ', '+00'//ss)
         call add_class('G1N    ', '+00'//ss)
         call add_class('G1B    ', '000+')
         call add_class('G1KX   ', '0++0')
         call add_class('G1KY   ', '0++0')
         call add_class('G1KZ   ', '0++0')
      end if


!     add to dcbxpr and dcbprl
!     ========================

      f = 1.0d0
      blocks = '+00'//ss

      if (need_1el_f) then

         irep(1) = isymax(1, 1)
         irep(2) = isymax(2, 1)
         irep(3) = isymax(3, 1)

         do i = 1, size(i_iter)
            s      = ' DIPLEN '
            s(1:1) = i_iter(i)
            call add_operator(s // '        ', &
                              blocks,          &
                              1,               &
                              1,               &
                              1 + irep(i),     &
                              (/1 + irep(i)/), &
                              (/f/),           &
                              (/s/))
         end do

         if (need_1el_g) then
            do iatom = 1, get_nr_atoms()*3
               do i = 1, size(i_iter)
                  s = prefix_zeros(iatom, 2) // ' DPG ' // i_iter(i)
                  call add_operator(s // '        ', &
                                    blocks,          &
                                    1,               &
                                    1,               &
                                    1,               &
                                    (/1/),           &
                                    (/f/),           &
                                    (/s/))
               end do
            end do
         end if
      end if

      if (need_1el_v) then

         irep(1) = isymax(1, 1)
         irep(2) = isymax(2, 1)
         irep(3) = isymax(3, 1)

         do i = 1, size(i_iter)
            s      = ' DIPVEL '
            s(1:1) = i_iter(i)
            call add_operator(s // '        ', &
                              blocks,          &
                              1,               &
                             -1,               &
                              1 + irep(i),     &
                              (/1 + irep(i)/), &
                              (/f/),           &
                              (/s/))
         end do
      end if

      if (need_1el_q) then

         irep(1) = 0
         irep(2) = ieor(isymax(1, 1), isymax(2, 1))
         irep(3) = ieor(isymax(1, 1), isymax(3, 1))
         irep(4) = 0
         irep(5) = ieor(isymax(2, 1), isymax(3, 1))
         irep(6) = 0

         do ij = 1, size(ij_iter)
            s = ij_iter(ij) // 'THETA '
            call add_operator(s // '        ',  &
                              blocks,           &
                              1,                &
                              1,                &
                              1 + irep(ij),     &
                              (/1 + irep(ij)/), &
                              (/f/),            &
                              (/s/))
         end do

         if (need_1el_g) then
            do iatom = 1, get_nr_atoms()*3
               do ij = 1, size(ij_iter)
                  s = prefix_zeros(iatom, 2) // 'QDG ' // ij_iter(ij)
                  call add_operator(s // '        ', &
                                    blocks,          &
                                    1,               &
                                    1,               &
                                    1,               &
                                    (/1/),           &
                                    (/f/),           &
                                    (/s/))
               end do
            end do
         end if
      end if

      if (need_1el_o .or. need_1el_b) then

         f = 1.0d0*cval
         blocks = '0++0'

         call add_operator('XANGMOM         ',                     &
                           blocks,                                 &
                           5,                                      &
                          -1,                                      &
                           (1 + ieor(isymax(3, 1), isymax(2, 1))), &
                           (/1 + isymax(3, 1), 1 + isymax(2, 1)/), &
                           (/f, f/),                               &
                           (/'ZDIPLEN ', 'YDIPLEN '/))
         call add_operator('YANGMOM         ',                     &
                           blocks,                                 &
                           6,                                      &
                          -1,                                      &
                           (1 + ieor(isymax(1, 1), isymax(3, 1))), &
                           (/1 + isymax(1, 1), 1 + isymax(3, 1)/), &
                           (/f, f/),                               &
                           (/'XDIPLEN ', 'ZDIPLEN '/))
         call add_operator('ZANGMOM         ',                     &
                           blocks,                                 &
                           7,                                      &
                          -1,                                      &
                           (1 + ieor(isymax(2, 1), isymax(1, 1))), &
                           (/1 + isymax(2, 1), 1 + isymax(1, 1)/), &
                           (/f, f/),                               &
                           (/'YDIPLEN ', 'XDIPLEN '/))

         if (need_1el_g) then

            f = -cval
            blocks = '0++0'

            do iatom = 1, get_nr_atoms()*3
               s1 = prefix_zeros(iatom, 2) // ' DPG Z'
               s2 = prefix_zeros(iatom, 2) // ' DPG Y'
               call add_operator(prefix_zeros(iatom, 3) // 'AMDRX        ', &
                                 blocks,                                    &
                                 5,                                         &
                                -1,                                         &
                                 1,                                         &
                                 (/1,  1/),                                 &
                                 (/f,  f/),                                 &
                                 (/s1, s2/))
               s1 = prefix_zeros(iatom, 2) // ' DPG X'
               s2 = prefix_zeros(iatom, 2) // ' DPG Z'
               call add_operator(prefix_zeros(iatom, 3) // 'AMDRY        ', &
                                 blocks,                                    &
                                 6,                                         &
                                -1,                                         &
                                 1,                                         &
                                 (/1,  1/),                                 &
                                 (/f,  f/),                                 &
                                 (/s1, s2/))
               s1 = prefix_zeros(iatom, 2) // ' DPG Y'
               s2 = prefix_zeros(iatom, 2) // ' DPG X'
               call add_operator(prefix_zeros(iatom, 3) // 'AMDRZ        ', &
                                 blocks,                                    &
                                 7,                                         &
                                -1,                                         &
                                 1,                                         &
                                 (/1,  1/),                                 &
                                 (/f,  f/),                                 &
                                 (/s1, s2/))
            end do
         end if

         blocks = '+00'//ss
         if (diamagnetic_via_pn) then
            f = 0.0d0
         else
            f = 1.0d0
         end if

         irep(1) = 0
         irep(2) = ieor(isymax(1, 1), isymax(2, 1))
         irep(3) = ieor(isymax(1, 1), isymax(3, 1))
         irep(4) = 0
         irep(5) = ieor(isymax(2, 1), isymax(3, 1))
         irep(6) = 0

         do ij = 1, size(ij_iter)
            s = ij_iter(ij) // 'SUSCGO'
            call add_operator(s // '        ',  &
                              blocks,           &
                              1,                &
                              1,                &
                              1 + irep(ij),     &
                              (/1 + irep(ij)/), &
                              (/f/),            &
                              (/s/))
         end do
      end if

      if (need_1el_b) then

         f = -2.0d0*cval*cval
         blocks = '000+'

         do i = 1, size(i_iter)
            s = 'dS/dB   '
            s(6:6) = i_iter(i)
            call add_operator('dbet/dB' // i_iter(i) // '        ', &
                              blocks,                               &
                              1,                                    &
                             -1,                                    &
                              1,                                    &
                              (/1/),                                &
                              (/f/),                                &
                              (/s/))
         end do

         do ij = 1, size(ij_iter)
            s = 'dS/dB2' // ij_iter(ij)
            call add_operator('dbet/dB2' // ij_iter(ij) // '      ', &
                              blocks,                                &
                              1,                                     &
                              1,                                     &
                              1,                                     &
                              (/1/),                                 &
                              (/f/),                                 &
                              (/s/))
         end do

         f = 1.0d0
         blocks = '+00'//ss

         do i = 1, size(i_iter)
            s = ' RM1H2  '
            s(1:1) = i_iter(i)
            call add_operator('dnuc/dB' // i_iter(i) // '        ', &
                              blocks,                               &
                              1,                                    &
                             -1,                                    &
                              1,                                    &
                              (/1/),                                &
                              (/f/),                                &
                              (/s/))
         end do

         do ij = 1, size(ij_iter)
            s = '  RM2H2 '
            s(1:2) = ij_iter(ij)
            call add_operator('dnuc/dB2' // ij_iter(ij) // '      ', &
                              blocks,                                &
                              1,                                     &
                              1,                                     &
                              1,                                     &
                              (/1/),                                 &
                              (/f/),                                 &
                              (/s/))
         end do

         f = cval
         blocks = '0++0'

         irep(1) = ieor(isymax(2, 1), isymax(3, 1))
         irep(2) = ieor(isymax(3, 1), isymax(1, 1))
         irep(3) = ieor(isymax(1, 1), isymax(2, 1))

         do i = 1, size(i_iter)
            s1 = ' RM1H3X '
            s2 = ' RM1H3Y '
            s3 = ' RM1H3Z '
            s1(1:1) = i_iter(i)
            s2(1:1) = i_iter(i)
            s3(1:1) = i_iter(i)
            call add_operator('dkin/dB' // i_iter(i) // '        ',     &
                              blocks,                                   &
                              8,                                        &
                              0,                                        &
                              1 + irep(i),                              &
                            (/1 + ieor(isymax(i, 2), isymax(1, 1)),     &
                              1 + ieor(isymax(i, 2), isymax(2, 1)),     &
                              1 + ieor(isymax(i, 2), isymax(3, 1))/),   &
                              (/f,  f,  f/),                            &
                              (/s1, s2, s3/))
         end do

         do ij = 1, size(ij_iter)
            s1 = '  RM2H3X'
            s2 = '  RM2H3Y'
            s3 = '  RM2H3Z'
            s1(1:2) = ij_iter(ij)
            s2(1:2) = ij_iter(ij)
            s3(1:2) = ij_iter(ij)
            call add_operator('dkin/dB2' // ij_iter(ij) // '      ', &
                              blocks,                                &
                              8,                                     &
                              0,                                     &
                              1,                                     &
                              (/1,  1,  1/),                         &
                              (/f,  f,  f/),                         &
                              (/s1, s2, s3/))
         end do

         f = 0.5d0*cval
         blocks = '0++0'

         call add_operator('dvec/dBX        ', &
                           blocks,             &
                           5,                  &
                           0,                  &
                           1,                  &
                           (/1, 1/),           &
                           (/f, f/),           &
                           (/'ZRM1RN  ', 'YRM1RN  '/))
         call add_operator('dvec/dBY        ', &
                           blocks,             &
                           6,                  &
                           0,                  &
                           1,                  &
                           (/1, 1/),           &
                           (/f, f/),           &
                           (/'XRM1RN  ', 'ZRM1RN  '/))
         call add_operator('dvec/dBZ        ', &
                           blocks,             &
                           7,                  &
                           0,                  &
                           1,                  &
                           (/1, 1/),           &
                           (/f, f/),           &
                           (/'YRM1RN  ', 'XRM1RN  '/))

         f = 1.0d0
         blocks = '+00'//ss

         irep(1) = ieor(isymax(2, 1), isymax(3, 1))
         irep(2) = ieor(isymax(3, 1), isymax(1, 1))
         irep(3) = ieor(isymax(1, 1), isymax(2, 1))

         do i = 1, size(i_iter)
            s = 'dS/dB   '
            s(6:6) = i_iter(i)
            call add_operator('dS/dB' // i_iter(i) // '          ', &
                              blocks,                               &
                              1,                                    &
                             -1,                                    &
                              1 + irep(i),                          &
                              (/1 + irep(i)/),                      &
                              (/f/),                                &
                              (/s/))
         end do

         do ij = 1, size(ij_iter)
            s = 'dS/dB2' // ij_iter(ij)
            call add_operator('dS/dB2' // ij_iter(ij) // '        ', &
                              blocks,                                &
                              1,                                     &
                              1,                                     &
                              1,                                     &
                              (/1/),                                 &
                              (/f/),                                 &
                              (/s/))
         end do

         do ij = 1, size(ij_iter)
            s = '>>S/B2' // ij_iter(ij)
            call add_operator(s // '        ', blocks, 1, 0, 1, (/1/), (/f/), (/s/))
         end do

         do ij = 1, size(ij_iter)
            s = '<>S/B2' // ij_iter(ij)
            call add_operator(s // '        ', blocks, 1, 0, 1, (/1/), (/f/), (/s/))
         end do

         irep(1) = ieor(isymax(2, 1), isymax(3, 1))
         irep(2) = ieor(isymax(3, 1), isymax(1, 1))
         irep(3) = ieor(isymax(1, 1), isymax(2, 1))

         do i = 1, size(i_iter)
            s = 'd|S>/dB '
            s(8:8) = i_iter(i)
            call add_operator('d|S>/dB' // i_iter(i) // '        ', &
                              blocks,                               &
                              1,                                    &
                              0,                                    &
                              1 + irep(i),                          &
                              (/1 + irep(i)/),                      &
                              (/f/),                                &
                              (/s/))
         end do

         if (openrsp_cfg_mllsus) then
            ! reproduces Dalton
            blocks = '00+0'
            f = -2.0d0*cval
         else
            ! reproduces DC RKB large c
            blocks = '0++0'
            f = -1.0d0*cval
         end if

         call add_operator('XXrdsusll       ',                      &
                           blocks,                                  &
                           5,                                       &
                           0,                                       &
                           (1 + ieor(isymax(1, 1), isymax(1, 1))),  &
                          (/1 + ieor(isymax(1, 2), isymax(3, 1)),   &
                            1 + ieor(isymax(1, 2), isymax(2, 1))/), &
                           (/f, f/),                                &
                           (/'XXRDSULZ', 'XXRDSULY'/))
         call add_operator('YXrdsusll       ',                      &
                           blocks,                                  &
                           6,                                       &
                           0,                                       &
                           (1 + ieor(isymax(2, 1), isymax(1, 1))),  &
                          (/1 + ieor(isymax(1, 2), isymax(1, 1)),   &
                            1 + ieor(isymax(1, 2), isymax(3, 1))/), &
                           (/f, f/),                                &
                           (/'XYRDSULX', 'XYRDSULZ'/))
         call add_operator('ZXrdsusll       ',                      &
                           blocks,                                  &
                           7,                                       &
                           0,                                       &
                           (1 + ieor(isymax(3, 1), isymax(1, 1))),  &
                          (/1 + ieor(isymax(1, 2), isymax(2, 1)),   &
                            1 + ieor(isymax(1, 2), isymax(1, 1))/), &
                           (/f, f/),                                &
                           (/'XZRDSULY', 'XZRDSULX'/))

         call add_operator('XYrdsusll       ',                      &
                           blocks,                                  &
                           5,                                       &
                           0,                                       &
                           (1 + ieor(isymax(1, 1), isymax(2, 1))),  &
                          (/1 + ieor(isymax(2, 2), isymax(3, 1)),   &
                            1 + ieor(isymax(2, 2), isymax(2, 1))/), &
                           (/f, f/),                                &
                           (/'YXRDSULZ', 'YXRDSULY'/))
         call add_operator('YYrdsusll       ',                      &
                           blocks,                                  &
                           6,                                       &
                           0,                                       &
                           (1 + ieor(isymax(2, 1), isymax(2, 1))),  &
                          (/1 + ieor(isymax(2, 2), isymax(1, 1)),   &
                            1 + ieor(isymax(2, 2), isymax(3, 1))/), &
                           (/f, f/),                                &
                           (/'YYRDSULX', 'YYRDSULZ'/))
         call add_operator('ZYrdsusll       ',                      &
                           blocks,                                  &
                           7,                                       &
                           0,                                       &
                           (1 + ieor(isymax(3, 1), isymax(2, 1))),  &
                          (/1 + ieor(isymax(2, 2), isymax(2, 1)),   &
                            1 + ieor(isymax(2, 2), isymax(1, 1))/), &
                           (/f, f/),                                &
                           (/'YZRDSULY', 'YZRDSULX'/))

         call add_operator('XZrdsusll       ',                      &
                           blocks,                                  &
                           5,                                       &
                           0,                                       &
                           (1 + ieor(isymax(1, 1), isymax(3, 1))),  &
                          (/1 + ieor(isymax(3, 2), isymax(3, 1)),   &
                            1 + ieor(isymax(3, 2), isymax(2, 1))/), &
                           (/f, f/),                                &
                           (/'ZXRDSULZ', 'ZXRDSULY'/))
         call add_operator('YZrdsusll       ',                      &
                           blocks,                                  &
                           6,                                       &
                           0,                                       &
                           (1 + ieor(isymax(2, 1), isymax(3, 1))),  &
                          (/1 + ieor(isymax(3, 2), isymax(1, 1)),   &
                            1 + ieor(isymax(3, 2), isymax(3, 1))/), &
                           (/f, f/),                                &
                           (/'ZYRDSULX', 'ZYRDSULZ'/))
         call add_operator('ZZrdsusll       ',                      &
                           blocks,                                  &
                           7,                                       &
                           0,                                       &
                           (1 + ieor(isymax(3, 1), isymax(3, 1))),  &
                          (/1 + ieor(isymax(3, 2), isymax(2, 1)),   &
                            1 + ieor(isymax(3, 2), isymax(1, 1))/), &
                           (/f, f/),                                &
                           (/'ZZRDSULY', 'ZZRDSULX'/))

         blocks = '+00'//ss
         if (diamagnetic_via_pn) then
            f = 0.0d0
         else
            f = 1.0d0
         end if

         irep(1) = 0
         irep(2) = ieor(isymax(1, 1), isymax(2, 1))
         irep(3) = ieor(isymax(1, 1), isymax(3, 1))
         irep(4) = 0
         irep(5) = ieor(isymax(2, 1), isymax(3, 1))
         irep(6) = 0

         do ij = 1, size(ij_iter)
            s      = '  DSUSNL'
            s(1:2) = ij_iter(ij)
            call add_operator(ij_iter(ij) // 'DSUSNL        ', &
                              blocks,                          &
                              1,                               &
                              0,                               &
                              1 + irep(ij),                    &
                              (/1 + irep(ij)/),                &
                              (/f/),                           &
                              (/s/))
         end do

         if (need_1el_f) then

            f = 1.0d0
            blocks = '+00'//ss

            call add_operator('X-CM1 X         ', blocks, 1, -1, 1, (/1/), (/f/), (/'X-CM1 X '/))
            call add_operator('X-CM1 Y         ', blocks, 1, -1, 1, (/1/), (/f/), (/'X-CM1 Y '/))
            call add_operator('X-CM1 Z         ', blocks, 1, -1, 1, (/1/), (/f/), (/'X-CM1 Z '/))
            call add_operator('Y-CM1 X         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Y-CM1 X '/))
            call add_operator('Y-CM1 Y         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Y-CM1 Y '/))
            call add_operator('Y-CM1 Z         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Y-CM1 Z '/))
            call add_operator('Z-CM1 X         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Z-CM1 X '/))
            call add_operator('Z-CM1 Y         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Z-CM1 Y '/))
            call add_operator('Z-CM1 Z         ', blocks, 1, -1, 1, (/1/), (/f/), (/'Z-CM1 Z '/))

            call add_operator('X-CM2XX         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2XX '/))
            call add_operator('X-CM2XY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2XY '/))
            call add_operator('X-CM2XZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2XZ '/))
            call add_operator('X-CM2YY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2YY '/))
            call add_operator('X-CM2YZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2YZ '/))
            call add_operator('X-CM2ZZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'X-CM2ZZ '/))
            call add_operator('Y-CM2XX         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2XX '/))
            call add_operator('Y-CM2XY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2XY '/))
            call add_operator('Y-CM2XZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2XZ '/))
            call add_operator('Y-CM2YY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2YY '/))
            call add_operator('Y-CM2YZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2YZ '/))
            call add_operator('Y-CM2ZZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Y-CM2ZZ '/))
            call add_operator('Z-CM2XX         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2XX '/))
            call add_operator('Z-CM2XY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2XY '/))
            call add_operator('Z-CM2XZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2XZ '/))
            call add_operator('Z-CM2YY         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2YY '/))
            call add_operator('Z-CM2YZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2YZ '/))
            call add_operator('Z-CM2ZZ         ', blocks, 1, 1, 1, (/1/), (/f/), (/'Z-CM2ZZ '/))

         end if
         if (need_1el_q) then
            f = 1.0d0
            blocks = '+00'//ss

            irep( 1) = isymax(1, 2)
            irep( 2) = isymax(2, 2)
            irep( 3) = isymax(3, 2)
            irep( 4) = isymax(1, 2)
            irep( 5) = isymax(2, 2)
            irep( 6) = isymax(3, 2)
            irep( 7) = isymax(1, 2)
            irep( 8) = isymax(2, 2)
            irep( 9) = isymax(3, 2)
            irep(10) = ieor(isymax(1,2),ieor(isymax(1,1),isymax(2,1)))
            irep(11) = ieor(isymax(2,2),ieor(isymax(1,1),isymax(2,1)))
            irep(12) = ieor(isymax(3,2),ieor(isymax(1,1),isymax(2,1)))
            irep(13) = ieor(isymax(1,2),ieor(isymax(1,1),isymax(3,1)))
            irep(14) = ieor(isymax(2,2),ieor(isymax(1,1),isymax(3,1)))
            irep(15) = ieor(isymax(3,2),ieor(isymax(1,1),isymax(3,1)))
            irep(16) = ieor(isymax(1,2),ieor(isymax(2,1),isymax(3,1)))
            irep(17) = ieor(isymax(2,2),ieor(isymax(2,1),isymax(3,1)))
            irep(18) = ieor(isymax(3,2),ieor(isymax(2,1),isymax(3,1)))

            call add_operator('XX-QDB X        ', blocks, 1, -1, 1+irep( 1), (/1+irep( 1)/), (/f/), (/'XX-QDB X'/))
            call add_operator('XX-QDB Y        ', blocks, 1, -1, 1+irep( 2), (/1+irep( 2)/), (/f/), (/'XX-QDB Y'/))
            call add_operator('XX-QDB Z        ', blocks, 1, -1, 1+irep( 3), (/1+irep( 3)/), (/f/), (/'XX-QDB Z'/))
            call add_operator('YY-QDB X        ', blocks, 1, -1, 1+irep( 4), (/1+irep( 4)/), (/f/), (/'YY-QDB X'/))
            call add_operator('YY-QDB Y        ', blocks, 1, -1, 1+irep( 5), (/1+irep( 5)/), (/f/), (/'YY-QDB Y'/))
            call add_operator('YY-QDB Z        ', blocks, 1, -1, 1+irep( 6), (/1+irep( 6)/), (/f/), (/'YY-QDB Z'/))
            call add_operator('ZZ-QDB X        ', blocks, 1, -1, 1+irep( 7), (/1+irep( 7)/), (/f/), (/'ZZ-QDB X'/))
            call add_operator('ZZ-QDB Y        ', blocks, 1, -1, 1+irep( 8), (/1+irep( 8)/), (/f/), (/'ZZ-QDB Y'/))
            call add_operator('ZZ-QDB Z        ', blocks, 1, -1, 1+irep( 9), (/1+irep( 9)/), (/f/), (/'ZZ-QDB Z'/))
            call add_operator('XY-QDB X        ', blocks, 1, -1, 1+irep(10), (/1+irep(10)/), (/f/), (/'XY-QDB X'/))
            call add_operator('XY-QDB Y        ', blocks, 1, -1, 1+irep(11), (/1+irep(11)/), (/f/), (/'XY-QDB Y'/))
            call add_operator('XY-QDB Z        ', blocks, 1, -1, 1+irep(12), (/1+irep(12)/), (/f/), (/'XY-QDB Z'/))
            call add_operator('XZ-QDB X        ', blocks, 1, -1, 1+irep(13), (/1+irep(13)/), (/f/), (/'XZ-QDB X'/))
            call add_operator('XZ-QDB Y        ', blocks, 1, -1, 1+irep(14), (/1+irep(14)/), (/f/), (/'XZ-QDB Y'/))
            call add_operator('XZ-QDB Z        ', blocks, 1, -1, 1+irep(15), (/1+irep(15)/), (/f/), (/'XZ-QDB Z'/))
            call add_operator('YZ-QDB X        ', blocks, 1, -1, 1+irep(16), (/1+irep(16)/), (/f/), (/'YZ-QDB X'/))
            call add_operator('YZ-QDB Y        ', blocks, 1, -1, 1+irep(17), (/1+irep(17)/), (/f/), (/'YZ-QDB Y'/))
            call add_operator('YZ-QDB Z        ', blocks, 1, -1, 1+irep(18), (/1+irep(18)/), (/f/), (/'YZ-QDB Z'/))
         end if
      end if

      if (need_1el_g) then

         f = 1.0d0
         blocks = '+00'//ss

         do iatom = 1, get_nr_atoms()*3
            s = 'G1O' // prefix_zeros(iatom, 3) // '  '
            call add_operator(s // '        ', &
                              blocks,          &
                              1,               &
                              1,               &
                              1,               &
                              (/1/),           &
                              (/f/),           &
                              (/s/))
         end do

         do iatom = 1, get_nr_atoms()*3
            s = 'G1N' // prefix_zeros(iatom, 3) // '  '
            call add_operator(s // '        ', &
                              blocks,          &
                              1,               &
                              1,               &
                              1,               &
                              (/1/),           &
                              (/f/),           &
                              (/s/))
         end do

         f = 2.0d0*cval*cval
         blocks = '000+'

         do iatom = 1, get_nr_atoms()*3
            s = 'G1B' // prefix_zeros(iatom, 3) // '  '
            call add_operator(s // '        ', &
                              blocks,          &
                              1,               &
                              1,               &
                              1,               &
                              (/1/),           &
                              (/f/),           &
                              (/s/))
         end do

         f = cval
         blocks = '0++0'

         do iatom = 1, get_nr_atoms()*3
            do i = 1, size(i_iter)
               s      = 'G1K     '
               s(4:4) = i_iter(i)
               s(5:7) = prefix_zeros(iatom, 3)
               call add_operator(s // '        ',    &
                                 blocks,             &
                                 1,                  &
                                -1,                  &
                                 1,                  &
                                 (/1/),              &
                                 (/f/),              &
                                 (/s/))
            end do
         end do
      end if

   end subroutine

   subroutine read_1el_integrals(op_name, P)

!     --------------------------------------------------------------------------
      character(*),  intent(in)    :: op_name
      type(matrix),  intent(inout) :: P
!     --------------------------------------------------------------------------
      real(8), allocatable         :: ptri(:)
      real(8), allocatable         :: op1int(:)
      logical, allocatable         :: first(:)
      logical                      :: lexpst(2) = .false.
      integer                      :: op_index, file_unit, idimension
      real(8)                      :: dummy_dp
      logical                      :: file_exists
      logical                      :: debug_me = .false.
      real(8), allocatable         :: work(:)
      integer                      :: lwork
!     --------------------------------------------------------------------------

      file_unit = lu1int
      call opnfil(file_unit, 'AOPROPER', 'UNKNOWN', 'read_1el_integrals')
      rewind file_unit

      op_index = index_on_aoproper(op_name)

      idimension = nnbbasx
      if (iprptim(op_index) == 0) then
         idimension = n2bbasx
      end if

      allocate(ptri(nz*idimension))
      allocate(op1int(idimension))
      allocate(first(nz))

      P%elms = 0.0d0

      P%pg_sym = iprpsym(op_index)
      P%ih_sym = iprptim(op_index)

      lwork = n2bbasxq*2
      allocate(work(lwork))

      call prpex2(op_index,   &
                  dummy_dp,   &
                  .false.,    &
                  lexpst,     &
                  idimension, &
                  ptri,       &
                  P%elms,     &
                  op1int,     &
                  first,      &
                  work, 1, lwork, 0)

      deallocate(work)
      close(file_unit, status = 'keep')

      if (debug_me) then
         call mat_print(P, label = 'debug P in AO basis')
      end if

      deallocate(ptri)
      deallocate(op1int)
      deallocate(first)

   end subroutine

   subroutine add_operator(op_name,           &
                           blocks,            &
                           typ,               &
                           tr_sym,            &
                           pg_sym,            &
                           pg_sym_components, &
                           prefactors,        &
                           labels)

!     --------------------------------------------------------------------------
      character(16), intent(in) :: op_name
      character(4),  intent(in) :: blocks
      integer,       intent(in) :: typ
      integer,       intent(in) :: tr_sym
      integer,       intent(in) :: pg_sym
      integer,       intent(in) :: pg_sym_components(:)
      real(8),       intent(in) :: prefactors(:)
      character(8),  intent(in) :: labels(:)
!     --------------------------------------------------------------------------
      integer                   :: i, ixpr, iprl
!     --------------------------------------------------------------------------

      ixpr = nprps + 1
      iprl = nprplbl

      if (ixpr > maxprps) then
         call quit('add_operator: ixpr > maxprps')
      end if

      prpnam(ixpr)  = op_name
      iprptyp(ixpr) = typ
      iprpsym(ixpr) = pg_sym
      iprptim(ixpr) = tr_sym

      do i = 1, size(prefactors)

         iprl = iprl + 1

         if (iprl > maxprplbl) then
            call quit('add_operator: iprl > maxprplbl')
         end if

         prplbl(iprl)  = labels(i)
         pdoint(iprl)  = blocks
         iprlrep(iprl) = pg_sym_components(i) - 1
         iprltyp(iprl) = tr_sym

         facprp(i, ixpr)  = prefactors(i)
         iprplbl(i, ixpr) = iprl

      end do

      nprps   = ixpr
      nprplbl = iprl

!     private
      op_name_to_aoproper_index(ixpr) = op_name

   end subroutine

   subroutine print_prp_common_blocks()

!     --------------------------------------------------------------------------
      integer :: i
!     --------------------------------------------------------------------------

      do i = 1, nprps
         write(*, *) 'dcbxpr:  --------------'
         write(*, *) 'dcbxpr:  i             ', i
         write(*, *) 'dcbxpr:  prpnam(i)     ', prpnam(i)
         write(*, *) 'dcbxpr:  iprptyp(i)    ', iprptyp(i)
         write(*, *) 'dcbxpr:  iprpsym(i)    ', iprpsym(i)
         write(*, *) 'dcbxpr:  iprptim(i)    ', iprptim(i)
         write(*, *) 'dcbxpr:  facprp(1, i)  ', facprp(1, i)
         write(*, *) 'dcbxpr:  facprp(2, i)  ', facprp(2, i)
         write(*, *) 'dcbxpr:  facprp(3, i)  ', facprp(3, i)
         write(*, *) 'dcbxpr:  iprplbl(1, i) ', iprplbl(1, i)
         write(*, *) 'dcbxpr:  iprplbl(2, i) ', iprplbl(2, i)
         write(*, *) 'dcbxpr:  iprplbl(3, i) ', iprplbl(3, i)
         write(*, *) 'dcbxpr:  --------------'
      end do

      do i = 1, nprplbl
         write(*, *) 'dcbprl:  -----------'
         write(*, *) 'dcbprl:  i          ', i
         write(*, *) 'dcbprl:  prplbl(i)  ', prplbl(i)
         write(*, *) 'dcbprl:  pdoint(i)  ', pdoint(i)
         write(*, *) 'dcbprl:  iprlrep(i) ', iprlrep(i)
         write(*, *) 'dcbprl:  iprltyp(i) ', iprltyp(i)
         write(*, *) 'dcbprl:  -----------'
      end do

      do i = 1, nprpcls
         write(*, *) 'dcbcls:  ----------'
         write(*, *) 'dcbcls:  i         ', i
         write(*, *) 'dcbcls:  clscal(i) ', clscal(i)
         write(*, *) 'dcbcls:  clsint(i) ', clsint(i)
         write(*, *) 'dcbcls:  clscmb(i) ', clscmb(i)
         write(*, *) 'dcbcls:  ----------'
      end do

   end subroutine

   subroutine add_class(label, blocks)

!     --------------------------------------------------------------------------
      character(7), intent(in) :: label
      character(4), intent(in) :: blocks
!     --------------------------------------------------------------------------
      integer                  :: i
!     --------------------------------------------------------------------------

      i = nprpcls + 1

      if (i > maxcls) then
         call quit('add_class: i > maxcls')
      end if

      clscal(i) = .true.
      clsint(i) = label
      clscmb(i) = blocks

!     important only for cartesian or spherical moment integrals
!     here just set it to something
      iordcl(i) = 0

      nprpcls = i

   end subroutine

   function index_on_aoproper(op_name)

!     --------------------------------------------------------------------------
      integer                   :: index_on_aoproper
      character(*), intent(in)  :: op_name
!     --------------------------------------------------------------------------
      integer                   :: i
!     --------------------------------------------------------------------------

      index_on_aoproper = 0
      do i = 1, max_aoproper_indices
         if (op_name_to_aoproper_index(i) == op_name(1:len_trim(op_name))) then
            index_on_aoproper = i
         end if
      end do

      if (index_on_aoproper == 0) then
         call quit('index_on_aoproper: index for ' &
                   // op_name(1:len_trim(op_name)) &
                   // ' not found')
      end if

   end function

   function prefix_zeros(i, n)

!     prefix_zeros(137, 6) returns '000137'

!     --------------------------------------------------------------------------
      integer,      intent(in) :: i
      integer,      intent(in) :: n
      character(n)             :: prefix_zeros
!     --------------------------------------------------------------------------
      integer                  :: k
      character(1)             :: c09(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
!     --------------------------------------------------------------------------

      do k = 1, n
         prefix_zeros(n-k+1:n-k+1) = c09(mod(i, 10**k)/10**(k-1))
      end do

   end function

#endif /* ifdef PRG_DIRAC */
end module
