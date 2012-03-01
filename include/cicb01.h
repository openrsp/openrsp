      integer :: NCIRED, NJCR, NKCR, JCROOT, KCROOT, JCONV_C
      real(8) :: THRCIT, residual_croot, energy_root, emy_ci

      integer, parameter :: MXCIRT = 150 ! max # of roots just like maxrts!!!

      COMMON /CICB01/ THRCIT,residual_croot(MXCIRT),energy_root(MXCIRT),&
     &                emy_ci, NCIRED, NJCR, NKCR, JCROOT(MXCIRT),       &
     &                KCROOT(MXCIRT), JCONV_C
