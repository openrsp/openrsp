#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using VAR_INT64 /Jan-2007 hjaaj
#endif
#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#if defined(__CVERSION__)
#ifdef COMMENT
! NOTE : MXSHEL should have value of MXSHEL in maxorb.h
#endif
#define MXSHEL 750
#define MAXNOD 999
#define NPARI  ((MAXNOD + 1) + 6)
#define MAXTSK ( MXSHEL * (MXSHEL + 1) / 2 )
extern struct common_infpar {
#if defined (VAR_INT64)
    long iprpar, ntask, ncode, ndegdi, master, mynum, mytid;
    long nodtot, nodeid[MAXNOD+1], nfmat, mtottk, parher, debug, pario;
    long timing, slave;
#else
    int  iprpar, ntask, ncode, ndegdi, master, mynum, mytid;
    int  nodtot, nodeid[MAXNOD+1], nfmat, mtottk, parher, debug, pario;
    int  timing, slave;
#endif
    char nodnam[MAXNOD][20], myname[20];
} daltoninfpar_;
#else
! File: infpar.h for Dalton; special information for parallel calculations
!
!     Parameters NPARI must be updated after changes (for parallelization)
!
!     NOTE: Integers  (IPRPAR,...,MASTER,...,MYTID)
!           Logicals  (TIMING,SLAVE)
!           Character (NODNAM,MYNAME) should NOT be sent to slaves
!     THUS: NPARI is length from NODTOT,...,PARIO
!
      INTEGER   MAXNOD, NPARI, MAXTSK
      PARAMETER ( MAXNOD = 999, NPARI = (MAXNOD + 1) + 6 )
      PARAMETER ( MAXTSK = (MXSHEL * (MXSHEL + 1))/2 )
      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID
      INTEGER NODTOT, NODEID(0:MAXNOD), NFMAT, MTOTTK
      LOGICAL PARHER, PARIO, DEBUG,     TIMING, SLAVE
      CHARACTER*20   NODNAM(0:MAXNOD), MYNAME
      COMMON /DALTONINFPAR/                                              &
     &        IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID         &
     &       ,NODTOT, NODEID, NFMAT, MTOTTK, PARHER, DEBUG, PARIO        &
     &       ,TIMING, SLAVE , NODNAM, MYNAME
#ifdef PRG_DIRAC
      EQUIVALENCE (NODTOT, NUMNOD)
!     ... for some obscure reason NODTOT in Dalton is called NUMNOD in Dirac.
!         This equivalence is important such that both e.g. abacus routines using NODTOT
!         and dirac routines using NUMNOD are working correctly together /hjaaj Dec 07
#endif

#if defined (VAR_INT64)
!     integer array my_STATUS contains MPI_SOURCE information.
!     Proper use of my_STATUS on 64-bit machines in 
!     combination with VAR_INT64 requires explicit declaration 
!     as INTEGER*4 /March-2007 sk - This is only true for 
!     32-bit MPI libraries.
!     64-bit MPI libraries require INTEGER my_STATUS (which will then be INTEGER*8)
!
#endif
      INTEGER my_STATUS

! -- end of infpar.h --
#endif
