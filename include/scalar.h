#if defined (SYS_ALLIANT)
CVD$ NOVECTOR
#endif
#if defined (SYS_HAL)
C$DIR SCALAR
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
CDIR$ NEXTSCALAR
#endif
#if defined (SYS_NEC)
*VDIR NOVECTOR
#endif
