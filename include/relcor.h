!
!     DARWN : Dawin correction
!     RMASSV: Mass-velocity correction
!     EKIN  : electronic kinetic energy,
!             may be used for one of the diagonal adiabatic corrections
!
#if defined (SYS_CRAY) 
      REAL  DARWN, RMASSV, EKIN
#else
      DOUBLE PRECISION DARWN, RMASSV, EKIN
#endif
      COMMON /RELCOR/ DARWN, RMASSV, EKIN
