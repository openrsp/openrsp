module interface_1el

   use matrix_defop
   use dalton_ifc !fixme ugly

   implicit none

   public oneint_ave

   private

contains

   subroutine oneint_ave(nr_atoms, what, D, DFD, R)
      !> structure containing the integral program settings
      integer,           intent(in)  :: nr_atoms
      character(*),      intent(in)  :: what
      type(matrix),      intent(in)  :: D, DFD
      real(8),           intent(out) :: R(:)
#ifndef LSDALTON_ONLY
#include <mxcent.h>
#include <taymol.h>
#endif
      real(8), pointer :: wrk(:)
      integer          :: lwrk
#ifdef LSDALTON_ONLY
      call quit('Cannot run oneint_ave, only new integral code is compiled',-1)
#else
      call save_D_and_DFD_for_ABACUS(.false., D, DFD)
      lwrk = 50*D%nrow**2+10000*D%nrow+50000000
      call di_select_wrk(wrk, lwrk)
      HESMOL(:3*nr_atoms,:3*nr_atoms) = 0
      ! SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
      ! &                  NODV,DIFDIP,HFONLY,NCLONE)
      call ONEDRV(wrk,lwrk,0,.true.,len(what),.true.,.true., &
                  .true.,.false.,.true.,.false.)
      ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
      R(1:9*nr_atoms**2) = reshape(HESMOL(:3*nr_atoms,:3*nr_atoms), (/9*nr_atoms**2/))
      call di_deselect_wrk(wrk, lwrk)
#endif
   end subroutine

   subroutine save_D_and_DFD_for_ABACUS(anti, D, DFD)
   !ajt aug09 Same as set_dsofso in fock-eval.f90, but pertrubed
   !          D and DFD is input instead of unperturbed D and F
      logical,      intent(in) :: anti !whether the integrals D and DFD
         ! are to be contracted with are symmetric or anti-symmetric
      type(matrix), intent(in) :: D, DFD
         ! perturbed (or un-) density and
         ! 'generalized Fock matrix' (energy-weighted density matrix)
      real(8)   Dtri(D%nrow*(D%nrow+1)/2)
      real(8) DFDtri(D%nrow*(D%nrow+1)/2)
      if (iszero(D)) then
         Dtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, D%elms, Dtri)
         if (     anti) call DGETAP(D%nrow, D%elms, Dtri)
         ! scale elms by 4 if anti, 2 if symm
         Dtri = Dtri * merge(4,2,anti)
      end if
      if (iszero(DFD)) then
         DFDtri = 0
      else
         if (.not.anti) call DGEFSP(D%nrow, DFD%elms, DFDtri)
         if (     anti) call DGETAP(D%nrow, DFD%elms, DFDtri)
         ! scale elms by 4 if anti, 2 if symm
         DFDtri = DFDtri * merge(4,2,anti)
      end if
      ! write fo files
#ifdef LSDALTON_ONLY
      call quit('Cannot call write_dsofso, only new integral code is compiled',-1)
#else
      call write_dsofso(Dtri,DFDtri)
#endif
   end subroutine

end module
