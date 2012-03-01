!
!VERSION : $Revision: 10420 $
! DATE    : $Date: 2010-03-22 00:46:02 +0100 (Mo, 22 Mär 2010) $
! FILE    : parluci.h
!
!                     PARALLEL INFORMATION
!
!     common blocks which are general to LUCITA and LUCIAREL
!     ======================================================
!
!
!     definitions, communication handles and group information
!
      INTEGER LUCI_MASTER,LUCI_MYPROC,LUCI_NMPROC,                      &
     &        MYNEW_ID, MYNEW_COMM, NEWCOMM_PROC, JCOMCOM, ICOMM,       &
     &        ICOMM_ID, ICOMM_SIZE, MY_GROUPN, N_MASTER,                &
     &        MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,              &
     &        N_MASTER_SM, MYNEW_ID_SM_C, MYNEW_COMM_SM_C,              &
     &        NEWCOMM_PROC_SM_C, N_MASTER_SM_C, IXCOMM, IXCOMM_SZ,      &
     &        IXCOMM_RK, IBLOCKCOMM, NUM_BLOCKS, NUM_BLOCKS2,           &
     &        NPTEST_VAR, IRC_SAVE,         JVEC_SF, LBLOCK_LUCI,       &
     &        IAM_NOT_INV,                    NCBLOCKS_USED, IIOMOD,    &
     &        NFLGRPS, NROOT_INFO, IDISTROUTE, IPRINTSTAT,I_NZERO_LEN,  &
     &        MSLVOUT_REL, LUWRT, IIOMOD_REL, NFLGRPS_REL, NROOT_SAVE,  &
     &        NVEC_SAVE, L_COMBI, ISHARE_DENSM, IPARALLELIO,            &
     &        LF2_ZERO, ISMEMFAC, IIJKL_ROD_LL, LEVEL_SM,IT_SHL,IC_SHL, &
     &        ICOUNTABLK_C, ISI_CALC_BL,IBI_MULT_BL,L_COMBI_S,          &
     &        L_COMBI_MAX, I_NZERO_LEN_S,I_NZERO_LEN_C, IFERM_SYM_T,    &
     &        LU_INFO
!
      PARAMETER (MSLVOUT_REL = 7)
!
      COMMON/LUCIPAR/ LUCI_MASTER,LUCI_MYPROC,LUCI_NMPROC,              &
     &                MY_GROUPN,MYNEW_COMM,MYNEW_ID,NEWCOMM_PROC,       &
     &                JCOMCOM,ICOMM,ICOMM_ID,ICOMM_SIZE,                &
     &                IBLOCKCOMM,N_MASTER,                              &
     &                MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,      &
     &                N_MASTER_SM, MYNEW_ID_SM_C, MYNEW_COMM_SM_C,      &
     &                NEWCOMM_PROC_SM_C, N_MASTER_SM_C,IXCOMM,          &
     &                IXCOMM_SZ,IXCOMM_RK, NUM_BLOCKS, NUM_BLOCKS2,     &
     &                NPTEST_VAR, IRC_SAVE,         JVEC_SF,            &
     &                LBLOCK_LUCI, IAM_NOT_INV,                         &
     &                NCBLOCKS_USED, IIOMOD, NFLGRPS, NROOT_INFO,       &
     &                IDISTROUTE, IPRINTSTAT, LUWRT, IIOMOD_REL,        &
     &                NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE, L_COMBI,      &
     &                I_NZERO_LEN, ISHARE_DENSM, IPARALLELIO,           &
     &                LF2_ZERO, ISMEMFAC, IIJKL_ROD_LL, LEVEL_SM,       &
     &                IT_SHL, IC_SHL, ICOUNTABLK_C, ISI_CALC_BL,        &
     &                IBI_MULT_BL,L_COMBI_S,L_COMBI_MAX,I_NZERO_LEN_S,  &
     &                I_NZERO_LEN_C, IFERM_SYM_T,                       &
     &                LU_INFO

!     currently we work with 9 scratch files
!     --------------------------------------
      integer, parameter :: nr_files = 9

!     MPI file lists and MPI file handles
!     --------------------------------------
      INTEGER IDIA,ILU1,ILU2,ILU3,ILU4,ILU5,ILU6,ILU7,ILUC
      INTEGER IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4, IALL_LU5,         &
     &        IALL_LU6, IALL_LU7, IALL_LUC, IIJKL_ROD, ILPRP_X,         &
     &        IJKL_PRE
      COMMON/LUCIAPFILE/ IDIA,ILU1,ILU2,ILU3,ILU4,ILU5,ILU6,ILU7,ILUC,  &
     &                   IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4,        &
     &                   IALL_LU5, IALL_LU6, IALL_LU7, IALL_LUC,        &
     &                   IIJKL_ROD, ILPRP_X, IJKL_PRE
#if defined (VAR_MPI)
      INTEGER MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL, FILE_INFO_OBJ
      INTEGER*8 MY_VEC1_IOFF, MY_VEC2_IOFF, MY_DIA_OFF, MY_LU1_OFF,     &
     &          MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF, MY_LU5_OFF,         &
     &          MY_LU6_OFF, MY_LU7_OFF, MY_LUC_OFF
      COMMON/LUCIFILEOFF/ MY_VEC1_IOFF, MY_VEC2_IOFF,                   &
     &                    MY_DIA_OFF, MY_LU1_OFF,                       &
     &                    MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF,           &
     &                    MY_LU5_OFF, MY_LU6_OFF, MY_LU7_OFF,           &
     &                    MY_LUC_OFF,                                   &
     &                    MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL,     &
     &                    FILE_INFO_OBJ
#endif
!     integer*8 block
      INTEGER*8           LEN_ALL_INT, IS_LENGTH_TT, NINT_TOTAL,        &
     &                    LEN_T_BUFF_NZ, MY_T_LEN, MY_T_LEN_D,          &
     &                    LMEMFREE_PTR
      COMMON/LUCIPAR_I8/  LEN_ALL_INT, IS_LENGTH_TT, NINT_TOTAL,        &
     &                    LEN_T_BUFF_NZ, MY_T_LEN, MY_T_LEN_D,          &
     &                    LMEMFREE_PTR
!     double precision block
      DOUBLE PRECISION  RNORM_FAC, TRUNC_FAC, mem_cp_time, xreadtimebi, &
     &                  xcomputesi
      COMMON/LUCI_REAL/ RNORM_FAC, TRUNC_FAC, mem_cp_time, xreadtimebi, &
     &                  xcomputesi
!     logical block
      LOGICAL CSCREEN, SHARED_M, REORD_IJKL, CHECK_TC, USE_EX_IDIA,     &
     &        SPLIT_IJKL, NO_INDIF, NO_CDIAF, CXPROPRUN, COMPDISTL,     &
     &        COMPRCCTOS, CHECKPOINT_LUCIX, CPBLCK_FILE_LUCIX,          &
     &        CKRCIONLY, LOWSRT_IJKL, TIMING_par
!
      COMMON/LUCI_LOG/                                                  &
     &        CSCREEN, SHARED_M, REORD_IJKL, CHECK_TC, USE_EX_IDIA,     &
     &        SPLIT_IJKL, NO_INDIF, NO_CDIAF, CXPROPRUN, COMPDISTL,     &
     &        COMPRCCTOS, CHECKPOINT_LUCIX, CPBLCK_FILE_LUCIX,          &
     &        CKRCIONLY, LOWSRT_IJKL, TIMING_par
!     character block
      CHARACTER*9 NIIJKL_ROD
        CHARACTER*1 SYMFLABEL, symflab_data(8)
      data symflab_data/'a','b','c','d','e','f','g','h'/
      COMMON/LUCI_CHAR/ NIIJKL_ROD, SYMFLABEL
! --- end of parluci.h ---
