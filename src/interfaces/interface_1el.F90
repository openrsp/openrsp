module interface_1el

   use matrix_defop
   use interface_molecule
   use interface_f77_memory
   use interface_pcm

   implicit none

   public oneint_ave
   public get1in_ave_ifc
   public onedrv_ave_ifc
   public di_get_overlap_and_h1
   public di_read_operator_int

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
      call f77_memory_select(work_len=lwrk, work=wrk)
      HESMOL(:3*nr_atoms,:3*nr_atoms) = 0
      ! SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
      ! &                  NODV,DIFDIP,HFONLY,NCLONE)
      call ONEDRV(wrk,lwrk,0,.true.,len(what),.true.,.true., &
                  .true.,.false.,.true.,.false.)
      ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
      R(1:9*nr_atoms**2) = reshape(HESMOL(:3*nr_atoms,:3*nr_atoms), (/9*nr_atoms**2/))
      call f77_memory_deselect(work_len=lwrk, work=wrk)
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

  subroutine ONEDRV_ave_ifc(fld, siz, ave, D, DFD)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix), optional, intent(in)  :: D, DFD
    !--------------------------------------------
#include <mxcent.h>
#include <taymol.h>

    real(8), allocatable :: Dtri(:)
    real(8), allocatable :: DFDtri(:)

    logical anti
    integer nc

    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    allocate(DFDtri(get_nr_ao()*(get_nr_ao()+1)/2))

    ! create triangularly packed matrices from D, DFD
    anti = (mod(count(fld=='MAG '),2) == 1)
    if (.not.present(D)) then
       Dtri = 0
    else if (anti) then
       call DGETAP(D%nrow, D%elms, Dtri)
    else !symm
       call DGEFSP(D%nrow, D%elms, Dtri)
    end if
    if (.not.present(DFD)) then
       DFDtri = 0
    else if (anti) then
       call DGETAP(DFD%nrow, DFD%elms, DFDtri)
    else !symm
       call DGEFSP(DFD%nrow, DFD%elms, DFDtri)
    end if
    ! write to files
    call WRITE_DSOFSO(Dtri, DFDtri)
    nc = 3 * get_nr_atoms()
    HESMOL(:nc,:nc) = 0
    !  SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,
    ! &                  DIFINT,NODC,NODV,DIFDIP,DIFQDP,
    ! &                  HFONLY,NCLONE,PCM)
    call ONEDRV(f77_memory, size(f77_memory), 5, .true., size(fld), &
                .true., .true., .true., .false., .false., &
                .true., .false., .false.)
    ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
    if (size(fld)==1) &
       call quit('error in ONEDRV_ave_ifc: not implemented')
    if (size(fld)==2) &
       ave = reshape(2*HESMOL(:nc,:nc), (/nc*nc/)) !factor 2 for total dens
    deallocate(Dtri)
    deallocate(DFDtri)
  end subroutine

  !> Call GET1IN in ABACUS
  subroutine GET1IN_ave_ifc(fld, siz, ave, D)
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix),           intent(in)  :: D
    !--------------------------------------------
    real(8)                             :: dummy(2)
#include <mxcent.h>
    real(8), allocatable :: Dtri(:)
    character(8), dimension(9*MXCENT) :: labint
    integer, dimension(9*MXCENT) :: intrep, intadr
    integer ncomp
    allocate(Dtri(  get_nr_ao()*(get_nr_ao()+1)/2))
    !ajt Dzero is dangerous! Should have been izero. =0 safe
    intrep = 0 !call dzero(intrep,9*MXCENT)
    intadr = 0 !call dzero(intadr,9*MXCENT)
    ! create triangularly packed matrix from D
    call DGEFSP(D%nrow, D%elms, Dtri)
    ncomp = 0
!      SUBROUTINE GET1IN(SINTMA,WORD,NCOMP,WORK,LWORK,LABINT,INTREP,
!     &                  INTADR,MPQUAD,TOFILE,KPATOM,TRIMAT,EXPVAL,
!     &                  EXP1VL,DENMAT,NPRINT)
    call GET1IN(dummy,'DPLGRA ',ncomp,f77_memory,size(f77_memory),labint,intrep, &
                       intadr,0,.false.,0,.false.,ave, &
                       .true.,Dtri,0)
    deallocate(Dtri)
  end subroutine

  !> \brief gets the overlap and one-electron Hamiltonian matrices
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return S contains the overlap matrix
  !> \return H1 contains the one-electron Hamiltonian matrix
  subroutine di_get_overlap_and_H1( S, H1 )
!   radovan: fixme it's not good for other codes to do S and H1 at the same time
!            split this in two
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(inout) :: S
    type(matrix), intent(inout) :: H1
    ! uses NBAST, NNBASX, N2BASX
#include <inforb.h>
    ! IO units, use LUPROP for file AOPROPER
#include<inftap.h>
    ! starts of the overlap and one electron Hamiltonian matrices in work memory
    integer work_ovlp, work_ham1
    ! PCM one-electron contributions
    integer work_pcm
    ! integer constants
    real(8), parameter :: one = 1.0D+00
    ! dummy stuff
    integer idummy
    ! external DALTON function finding the corresponding label
    logical FNDLAB
    ! reads overlap and one electron Hamiltonian matrices by calling RDONEL
    work_ovlp = get_f77_memory_next()
    work_ham1 = work_ovlp + NNBASX
    call set_f77_memory_next(work_ham1 + NNBASX)

    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', get_f77_memory_next()-1, get_f77_memory_total() )
    call RDONEL( 'OVERLAP', .true., f77_memory(work_ovlp), NNBASX )
    call RDONEL( 'ONEHAMIL', .true., f77_memory(work_ham1), NNBASX )
    ! PCM one-electron contributions
    if ( get_is_pcm_calculation() ) then
      work_pcm = get_f77_memory_next()
      call set_f77_memory_next(work_pcm + NNBASX)
      if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', get_f77_memory_next()-1, get_f77_memory_total() )
#ifdef USE_WAVPCM
      call pcm_ao_rsp_1elfock( f77_memory(work_pcm) )
      ! adds to one-electron Hamiltonian
      f77_memory( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    = f77_memory( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    + f77_memory( work_pcm  : work_pcm  + NNBASX - 1 )
#endif
      ! cleans
      call set_f77_memory_next(work_pcm)
    end if
    ! fills the data into matrices S and H1
    !N N2BASX = NBAST * NBAST
    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'DSPTSI', get_f77_memory_next()+N2BASX-1, get_f77_memory_total() )
    ! gets S
    call DSPTSI( NBAST, f77_memory(work_ovlp), S%elms )
    ! gets H1
    call DSPTSI( NBAST, f77_memory(work_ham1), H1%elms )
    ! clean
    call set_f77_memory_next(work_ovlp)
  end subroutine



  !> \brief gets property integrals from file AOPROPER
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param prop_lab is the label of integral
  !> \param init_prop indicates if initialize the integral matrix
  !> \return prop_int contains the integral matrix
  subroutine di_read_operator_int( prop_lab, prop_int )
    implicit integer (i,m-n)
#include <implicit.h>
    character*(8), intent(in) :: prop_lab
    type(matrix), intent(inout) :: prop_int
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
    ! IO units, use LUPROP for file AOPROPER
#include<inftap.h>
    ! external DALTON function finding the corresponding label
    logical FNDLB2
    ! information when calling subroutine FNDLB2
    character*8 RTNLBL(2)
    ! dummy stuff
    integer IDUMMY
    ! one-electron Hamiltonian
    if ( prop_lab == 'ONEHAMIL' ) then
      call QUIT( 'Not implemented!' )
      call RDONEL( 'ONEHAMIL', ANTI, f77_memory(get_f77_memory_next()), NNBASX )
      call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
    else
      ! closes file AOPROPER first
      if ( LUPROP > 0 ) call GPCLOSE( LUPROP, 'KEEP' )
      call GPOPEN( LUPROP, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', IDUMMY, .false. )
      ! finds the label
      if ( FNDLB2( prop_lab, RTNLBL, LUPROP ) ) then
        ! square matrix
        if ( RTNLBL(2) == 'SQUARE' ) then
          call READT( LUPROP, N2BASX, prop_int%elms )
        ! symmetric matrix
        else if ( RTNLBL(2) == 'SYMMETRI' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DSPTSI( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        ! anti-symmetric matrix
        else if ( RTNLBL(2) == 'ANTISYMM' ) then
          call READT( LUPROP, NNBASX, f77_memory(get_f77_memory_next()) )
          call DAPTGE( NBAST, f77_memory(get_f77_memory_next()), prop_int%elms )
        else
          call QUIT( 'Error: No symmetry label on AOPROPER!' )
        end if
        ! closes file AOPROPER
        call GPCLOSE( LUPROP, 'KEEP' )
      else
        call GPCLOSE( LUPROP, 'KEEP' )
        call QUIT( 'Integrals with label '''//prop_lab// &
                   ''' not found on file ''AOPROPER''!' )
      end if
    end if
  end subroutine


end module
