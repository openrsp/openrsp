!
! file: infvar.h (originally from SIRIUS)
!       contains INFormation on number and symmetry of mcscf VARiables
!
      INTEGER MAXWOP, JWOP,                                             &
     &        NCONF, NWOPT, NVAR, JWOPSY, NWOPH, NVARH, NCDETS
!
!     MAXWOP = maximum number of orbital rotations (dimension of JWOP)
!     JWOP(1:2,*) = /P, Q/, rotation from orbital P to orbital Q
!     NCONF  = number of configurations (CSFs or determinants)
!     NCDETS = number of determinants
!     NWOPT  = number of orbital rotations = number of non-redundant elements in JWOP
!     NVAR   = NCONF + NWOPT = total number of variables
!     JWOPSY = symmetry of orbital rotations
!     NWOPH  = NWOPT + no. of active-active rotations needed for Hessian
!     NVARH  = NCONF + NWOPH
!
      PARAMETER ( MAXWOP = 230000 )
      COMMON /INFVAR/ NCONF,NWOPT,NVAR,JWOPSY,NWOPH,NVARH,NCDETS,       &
     &                JWOP(2,MAXWOP)
