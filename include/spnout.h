! -- FILE: spnout.h --
      LOGICAL DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,
     &        FCFIN, SPNISO, NCSPNI,
     &        SOS, SOSSPN, SOSOCC, SOSOCS
      COMMON /SPNOUT/ ABUND,                                        ! real*8
     &        ISPPRI, ISOTPS(MXCENT),                               ! integer
     &        NSTATS, NSTATT, NSTATI, NSTATF, NITRST, NUCSPI,
     &        DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,     ! logical
     &        FCFIN, SPNISO, NCSPNI(MXCENT),
     &        SOS, SOSSPN, SOSOCC, SOSOCS
! -- end of spnout.h --
