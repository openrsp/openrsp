c  -------------------------------------------------------------------
c          PVM version 3.3:  Parallel Virtual Machine System
c  -------------------------------------------------------------------
c
c     Definitions to be included with
c     User's Fortran application
c     ----------------------------------

      integer PVMTASKDEFAULT, PVMTASKHOST, PVMTASKARCH, PVMTASKDEBUG
      integer PVMTASKTRACE, PVMMPPFRONT, PVMHOSTCOMPL
      integer PVMHOST, PVMARCH, PVMDEBUG, PVMTRACE
      integer PVMDATADEFAULT, PVMDATARAW, PVMDATAINPLACE
      integer PVMDEFAULT, PVMRAW, PVMINPLACE
      integer PVMTASKEXIT, PVMHOSTDELETE, PVMHOSTADD
      integer PVMROUTE, PVMDEBUGMASK, PVMAUTOERR
      integer PVMOUTPUTTID, PVMOUTPUTCODE, PVMRESVTIDS
      integer PVMTRACETID, PVMTRACECODE, PVMFRAGSIZE
      integer PVMDONTROUTE, PVMALLOWDIRECT, PVMROUTEDIRECT
      integer STRING, BYTE1, INTEGER2, INTEGER4
      integer REAL4, COMPLEX8, REAL8, COMPLEX16

      integer PvmOk, PvmSysErr, PvmBadParam, PvmMismatch
      integer PvmNoData, PvmNoHost, PvmNoFile, PvmNoMem
      integer PvmBadMsg, PvmNoBuf, PvmNoSuchBuf
      integer PvmNullGroup, PvmDupGroup, PvmNoGroup
      integer PvmNotInGroup, PvmNoinst, PvmHostFail, PvmNoParent
      integer PvmNotImpl, PvmDSysErr, PvmBadVersion, PvmOutOfRes
      integer PvmDupHost, PvmCantStart, PvmAlready, PvmNoTask
      integer PvmNoEntry, PvmDupEntry

c     --------------------
c     spawn 'flag' options
c     --------------------
      parameter( PVMTASKDEFAULT  =  0)
      parameter( PVMTASKHOST     =  1)
      parameter( PVMTASKARCH     =  2)
      parameter( PVMTASKDEBUG    =  4)
      parameter( PVMTASKTRACE    =  8)
      parameter( PVMMPPFRONT     = 16)
      parameter( PVMHOSTCOMPL    = 32)
c     --------------------------------
c     old option names still supported
c     --------------------------------
      parameter( PVMHOST  =  1)
      parameter( PVMARCH  =  2)
      parameter( PVMDEBUG =  4)
      parameter( PVMTRACE =  8)

c     -------------------------
c     buffer 'encoding' options
c     -------------------------
      parameter( PVMDATADEFAULT = 0)
      parameter( PVMDATARAW     = 1)
      parameter( PVMDATAINPLACE = 2)
c     --------------------------------
c     old option names still supported
c     --------------------------------
      parameter( PVMDEFAULT = 0)
      parameter( PVMRAW     = 1)
      parameter( PVMINPLACE = 2)

c     ----------------------
c     notify 'about' options
c     ----------------------
      parameter( PVMTASKEXIT   = 1 )
      parameter( PVMHOSTDELETE = 2 )
      parameter( PVMHOSTADD    = 3 )

c     --------------------------------
c     packing/unpacking 'what' options
c     --------------------------------
      parameter( STRING   = 0)
      parameter( BYTE1    = 1)
      parameter( INTEGER2 = 2)
      parameter( INTEGER4 = 3)
      parameter( REAL4    = 4)
      parameter( COMPLEX8 = 5)
      parameter( REAL8    = 6)
      parameter( COMPLEX16= 7)

c     --------------------------------
c     setopt/getopt options for 'what'
c     --------------------------------
      parameter( PVMROUTE      = 1)
      parameter( PVMDEBUGMASK  = 2)
      parameter( PVMAUTOERR    = 3)
      parameter( PVMOUTPUTTID  = 4)
      parameter( PVMOUTPUTCODE = 5)
      parameter( PVMTRACETID   = 6)
      parameter( PVMTRACECODE  = 7)
      parameter( PVMFRAGSIZE   = 8)
      parameter( PVMRESVTIDS   = 9)

c     --------------------------------------------
c     routing options for 'how' in setopt function
c     --------------------------------------------
      parameter( PVMDONTROUTE  = 1)
      parameter( PVMALLOWDIRECT= 2)
      parameter( PVMROUTEDIRECT= 3)

c     --------------------------
c     error 'info' return values
c     --------------------------
      parameter( PvmOk         =   0)
      parameter( PvmBadParam   =  -2)
      parameter( PvmMismatch   =  -3)
      parameter( PvmNoData     =  -5)
      parameter( PvmNoHost     =  -6)
      parameter( PvmNoFile     =  -7)
      parameter( PvmNoMem      = -10)
      parameter( PvmBadMsg     = -12)
      parameter( PvmSysErr     = -14)
      parameter( PvmNoBuf      = -15)
      parameter( PvmNoSuchBuf  = -16)
      parameter( PvmNullGroup  = -17)
      parameter( PvmDupGroup   = -18)
      parameter( PvmNoGroup    = -19)
      parameter( PvmNotInGroup = -20)
      parameter( PvmNoInst     = -21)
      parameter( PvmHostFail   = -22)
      parameter( PvmNoParent   = -23)
      parameter( PvmNotImpl    = -24)
      parameter( PvmDSysErr    = -25)
      parameter( PvmBadVersion = -26)
      parameter( PvmOutOfRes   = -27)
      parameter( PvmDupHost    = -28)
      parameter( PvmCantStart  = -29)
      parameter( PvmAlready    = -30)
      parameter( PvmNoTask     = -31)
      parameter( PvmNoEntry    = -32)
      parameter( PvmDupEntry   = -33)
