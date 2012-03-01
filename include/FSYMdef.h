/* file FSYMdef.h:
 * Match Fortran name mangling. If the Fortran compiler does not
 * mangle names, define FUNDERSCORE=0 in CFLAGS.  g77 and compaq fort
 * (cryptically referred to with HAVE_GCPP below) for linux-alpha both
 * insert a second underscore if routine name contains at least one
 * underscore /hjaaj Oct04 */
#if defined(NO_UNDERSCORE) || (defined(FUNDERSCORE) &&FUNDERSCORE == 0)
#define FSYM(a) a
#define FSYM2(a) a
#else
#define FSYM(a) a ## _
#if (defined(FUNDERSCORE) && FUNDERSCORE == 2)
#define FSYM2(a) a ## __
#else
#define FSYM2(a) a ## _
#endif
#endif
