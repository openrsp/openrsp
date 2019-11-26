/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud

  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

*/

#if defined(QCMATRIX_SINGLE_PRECISION)
#define QCREAL_FMT "f"
#else
#define QCREAL_FMT "lf"
#endif

#include "OpenRSPTestData.h"

/* <function name='OpenRSPTestReadMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-18'>
     Reads a dat file and puts the data into matrices
     <param name='file_mat' direction='in'>
       Name of the dat file
     </param>
     <param name='num_skip_mat' direction='in'>
       Number of the first few matrices to skip
     </param>
     <param name='num_read_mat' direction='in'>
       Number of matrices to read
     </param>
     <param name='mat_read' direction='inout'>
       Matrices read
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPTestReadMat(const char *file_mat,
                              const QInt num_skip_mat,
                              const QInt num_read_mat,
                              QcMat *mat_read[])
{
    FILE *fp_dat;          /* file pointer */
    QInt dat_num_read;     /* number of arguments read from the dat file */
    QInt dat_num_mat;      /* number of matrices in the dat file */
    QInt dat_num_row;      /* number of rows of matrices */
    QInt dat_num_col;      /* number of columns of matrices */
    QInt dat_num_val;      /* number of values of matrices */
    QReal *dat_real_elms;  /* elements of real part */
    QReal *dat_imag_elms;  /* elements of imaginary part */
    QInt imat;             /* incremental recorder over matrices */
    QInt ival;             /* incremental recorder over values */
    QErrorCode ierr;       /* error information */
    /* parameters for setting matrices */
    const QInt NUM_BLOCKS=1;
#if defined(OPENRSP_ZERO_BASED)
    const QInt IDX_BLOCK_ROW[]={0};
    const QInt IDX_BLOCK_COL[]={0};
    const QInt IDX_FIRST_ROW=0;
    const QInt IDX_FIRST_COL=0;
#else
    const QInt IDX_BLOCK_ROW[]={1};
    const QInt IDX_BLOCK_COL[]={1};
    const QInt IDX_FIRST_ROW=1;
    const QInt IDX_FIRST_COL=1;
#endif
    QcSymType sym_type;
    QcDataType data_type[1];
    fp_dat = fopen(file_mat, "r");
    if (fp_dat==NULL) {
        printf("OpenRSPTestReadMat>> file: %s\n", file_mat);
        QErrorExit(FILE_AND_LINE, "open the dat file");
    }
    /* reads the number of matrices */
    if (fscanf(fp_dat, "%"QINT_FMT"", &dat_num_mat)!=1) {
        QErrorExit(FILE_AND_LINE, "read number of matrices");
    }
    /* checks the number of matrices to read */
    if (num_skip_mat+num_read_mat>dat_num_mat) {
        printf("OpenRSPTestReadMat>> number of skipping matrices %"QINT_FMT"\n",
               num_skip_mat);
        printf("OpenRSPTestReadMat>> number of reading matrices %"QINT_FMT"\n",
               num_read_mat);
        printf("OpenRSPTestReadMat>> number of matrices on the file %"QINT_FMT"\n",
               dat_num_mat);
        QErrorExit(FILE_AND_LINE, "not enough matrices on the file");
    }
    /* skips the first few matrices, including the '\n' of the first line */
    dat_num_read = fgetc(fp_dat);
    for (imat=0; imat<num_skip_mat; imat++) {
        /*fscanf(fp_dat, "%*[^\n]\n", NULL);*/
        /* reads the information of matrix */
        dat_num_read = fscanf(fp_dat,
                              "%d %d %"QINT_FMT" %"QINT_FMT"",
                              &sym_type, &data_type[0], &dat_num_row, &dat_num_col);
        if (dat_num_read!=4) {
            printf("OpenRSPTestReadMat>> number of arguments read %"QINT_FMT"\n",
                   dat_num_read);
            QErrorExit(FILE_AND_LINE, "read information of matrix");
        }
        dat_num_read = fgetc(fp_dat);
        do {
            dat_num_read = fgetc(fp_dat);
        } while (dat_num_read!='\n');
        /* skips imaginary part */
        if (data_type[0]!=QREALMAT) {
            do {
                dat_num_read = fgetc(fp_dat);
            } while (dat_num_read!='\n');
        }
    }
    /* initializes the elments of imaginary part */
    dat_imag_elms = NULL;
    for (imat=0; imat<num_read_mat; imat++) {
        /* reads the information of matrix */
        dat_num_read = fscanf(fp_dat,
                              "%d %d %"QINT_FMT" %"QINT_FMT"",
                              &sym_type, &data_type[0], &dat_num_row, &dat_num_col);
        if (dat_num_read!=4) {
            printf("OpenRSPTestReadMat>> number of arguments read %"QINT_FMT"\n",
                   dat_num_read);
            QErrorExit(FILE_AND_LINE, "read information of matrix");
        }
        /* allocates memory for elements of real part */
        dat_num_val = dat_num_row*dat_num_col;
        dat_real_elms = (QReal *)malloc(dat_num_val*sizeof(QReal));
        if (dat_real_elms==NULL) {
            printf("OpenRSPTestReadMat>> matrix %"QINT_FMT"\n",
                   imat+num_skip_mat);
            printf("OpenRSPTestReadMat>> number of elements %"QINT_FMT"\n",
                   dat_num_val);
            QErrorExit(FILE_AND_LINE, "allocates memory for elements of real part");
        }
        /* reads elements of real part */
        for (ival=0; ival<dat_num_val; ival++) {
            if (fscanf(fp_dat, "%"QCREAL_FMT"", &dat_real_elms[ival])!=1) {
                printf("OpenRSPTestReadMat>> %"QINT_FMT"-th element\n", ival);
                QErrorExit(FILE_AND_LINE, "read elements of real part");
            }
        }
        if (data_type[0]!=QREALMAT) {
            /* allocates memory for elements of imaginary part */
            dat_imag_elms = (QReal *)malloc(dat_num_val*sizeof(QReal));
            if (dat_imag_elms==NULL) {
                printf("OpenRSPTestReadMat>> matrix %"QINT_FMT"\n",
                       imat+num_skip_mat);
                printf("OpenRSPTestReadMat>> number of elements %"QINT_FMT"\n",
                       dat_num_val);
                QErrorExit(FILE_AND_LINE, "allocates memory for elements of imaginary part");
            }
            /* reads elements of imaginary part */
            for (ival=0; ival<dat_num_val; ival++) {
                if (fscanf(fp_dat, "%"QCREAL_FMT"", &dat_imag_elms[ival])!=1) {
                    printf("OpenRSPTestReadMat>> %"QINT_FMT"-th element\n", ival);
                    QErrorExit(FILE_AND_LINE, "read elements of imaginary part");
                }
            }
        }
        /* creates the matrix */
        ierr = QcMatCreate(mat_read[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatCreate()");
        ierr = QcMatBlockCreate(mat_read[imat], NUM_BLOCKS);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate()");
        ierr = QcMatSetSymType(mat_read[imat], sym_type);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetSymType()");
        ierr = QcMatSetDataType(mat_read[imat],
                                NUM_BLOCKS,
                                IDX_BLOCK_ROW,
                                IDX_BLOCK_COL,
                                data_type);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType()");
        ierr = QcMatSetDimMat(mat_read[imat], dat_num_row, dat_num_col);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat()");
        ierr = QcMatAssemble(mat_read[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble()");
        /* saves the values into the matrix */
        ierr = QcMatSetValues(mat_read[imat],
                              IDX_BLOCK_ROW[0],
                              IDX_BLOCK_COL[0],
                              IDX_FIRST_ROW,
                              dat_num_row,
                              IDX_FIRST_COL,
                              dat_num_col,
                              dat_real_elms,
                              dat_imag_elms);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetValues()");
        free(dat_real_elms);
        dat_real_elms = NULL;
        if (dat_imag_elms!=NULL) {
            free(dat_imag_elms);
            dat_imag_elms = NULL;
        }
    }
    fclose(fp_dat);
    return QSUCCESS;
}
/* <function name='OpenRSPTestZeroMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-18'>
     Zero out a few matrices
     <param name='num_mat' direction='in'>
       Number of matrices to be zeroed out
     </param>
     <param name='mat_zero' direction='inout'>
       Matrices to be zeroed out
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPTestZeroMat(const QInt num_mat, QcMat *mat_zero[])
{
    QInt imat;        /* incremental recorder over matrices */
    QErrorCode ierr;  /* error information */
    /* parameters for setting matrices */
    const QInt NUM_BLOCKS=1;
#if defined(OPENRSP_ZERO_BASED)
    const QInt IDX_BLOCK_ROW[]={0};
    const QInt IDX_BLOCK_COL[]={0};
#else
    const QInt IDX_BLOCK_ROW[]={1};
    const QInt IDX_BLOCK_COL[]={1};
#endif
    QcDataType DATA_TYPE[1]={QREALMAT};
    const QInt DAT_NUM_ROW=32;
    const QInt DAT_NUM_COL=32;
    for (imat=0; imat<num_mat; imat++) {
        ierr = QcMatCreate(mat_zero[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatCreate()");
        ierr = QcMatBlockCreate(mat_zero[imat], NUM_BLOCKS);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate()");
        ierr = QcMatSetSymType(mat_zero[imat], QSYMMAT);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetSymType()");
        ierr = QcMatSetDataType(mat_zero[imat],
                                NUM_BLOCKS,
                                IDX_BLOCK_ROW,
                                IDX_BLOCK_COL,
                                DATA_TYPE);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType()");
        ierr = QcMatSetDimMat(mat_zero[imat], DAT_NUM_ROW, DAT_NUM_COL);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat()");
        ierr = QcMatAssemble(mat_zero[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble()");
        ierr = QcMatZeroEntries(mat_zero[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatZeroEntries()");
    }
    return QSUCCESS;
}

