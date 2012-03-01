#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C CRAY: do not vectorize this subroutine
CDIR$ NOVECTOR
#endif
#if defined (SYS_NEC)
C We would like to avoid vectorizing the entire subroutine, but a
C valid NEC directive is badly needed.
*VDIR NOVECTOR
#endif
