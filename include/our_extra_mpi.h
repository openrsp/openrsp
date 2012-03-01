/* FILE our_extra_mpi.h :
 * to handle if the Fortran integers are integer*4 or integer*8
 * */

#ifdef VAR_INT64
#define fortran_MPI_INT MPI_LONG_LONG_INT
#else
#define fortran_MPI_INT MPI_INT
#endif

/* end of our_extra_mpi.h */
