/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

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
    FILE *fp_dat;       /* file pointer */
    QInt dat_num_read;  /* number of arguments read from the dat file */
    QInt dat_num_mat;   /* number of matrices in the dat file */
    QInt dat_num_row;   /* number of rows of matrices */
    QInt dat_num_col;   /* number of columns of matrices */
    QInt dat_num_val;   /* number of values of matrices */
    QReal *dat_values;  /* values of matrices */
    QInt imat;          /* incremental recorder over matrices */
    QInt ival;          /* incremental recorder over values */
    QBool assembled;    /* matrix assembled or not */
    QErrorCode ierr;    /* error information */
    /* parameters for setting matrices */
    const QInt NUM_BLOCKS=1;
#if defined(ZERO_BASED_NUMBERING)
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
    QcDataType DATA_TYPE[1]={QREALMAT};
    fp_dat = fopen(file_mat, "r");
    if (fp_dat==NULL) {
        printf("OpenRSPTestReadMat>> file: %s\n", file_mat);
        QErrorExit(FILE_AND_LINE, "open the dat file");
    }
    dat_num_read = fscanf(fp_dat,
                          "%"QINT_FMT" %"QINT_FMT" %"QINT_FMT"",
                          &dat_num_mat, &dat_num_row, &dat_num_col);
    if (dat_num_read!=3) {
        printf("OpenRSPTestReadMat>> number of arguments read %"QINT_FMT"\n",
               dat_num_read);
        QErrorExit(FILE_AND_LINE, "read information of data");
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
    dat_num_val = dat_num_row*dat_num_col;
    dat_values = (QReal *)malloc(dat_num_val*sizeof(QReal));
    if (dat_values==NULL) {
        printf("OpenRSPTestReadMat>> number of values %"QINT_FMT"\n",
               dat_num_val);
        QErrorExit(FILE_AND_LINE, "allocates memory for values");
    }
    /* skips the first few lines, including the '\n' of the first line */
    for (imat=0; imat<=num_skip_mat; imat++) {
        /*fscanf(fp_dat, "%*[^\n]\n", NULL);*/
        do {
            dat_num_read = fgetc(fp_dat);
        } while (dat_num_read!='\n');
    }
    for (imat=0; imat<num_read_mat; imat++) {
        for (ival=0; ival<dat_num_val; ival++) {
            if (fscanf(fp_dat, "%"QCREAL_FMT"", &dat_values[ival])!=1) {
                printf("OpenRSPTestReadMat>> %"QINT_FMT"-th element\n", ival);
                QErrorExit(FILE_AND_LINE, "read values");
            }
        }
        /* checks if the matrix is assembled or not */
        ierr = QcMatIsAssembled(mat_read[imat], &assembled);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatIsAssembled()");
        if (assembled==QFALSE) {
            ierr = QcMatBlockCreate(mat_read[imat], NUM_BLOCKS);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatBlockCreate()");
            ierr = QcMatSetDataType(mat_read[imat],
                                    NUM_BLOCKS,
                                    IDX_BLOCK_ROW,
                                    IDX_BLOCK_COL,
                                    DATA_TYPE);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDataType()");
            ierr = QcMatSetDimMat(mat_read[imat], dat_num_row, dat_num_col);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetDimMat()");
            ierr = QcMatAssemble(mat_read[imat]);
            QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatAssemble()");
        }
        /* saves the values into the matrix */
        ierr = QcMatSetValues(mat_read[imat],
                              IDX_BLOCK_ROW[0],
                              IDX_BLOCK_COL[0],
                              IDX_FIRST_ROW,
                              dat_num_row,
                              IDX_FIRST_COL,
                              dat_num_col,
                              dat_values,
                              NULL);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatSetValues()");
    }
    free(dat_values);
    dat_values = NULL;
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
    QInt imat;          /* incremental recorder over matrices */
    QBool assembled;    /* matrix assembled or not */
    QErrorCode ierr;    /* error information */
    /* parameters for setting matrices */
    const QInt NUM_BLOCKS=1;
#if defined(ZERO_BASED_NUMBERING)
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
        /* checks if the matrix is assembled or not */
        ierr = QcMatIsAssembled(mat_zero[imat], &assembled);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatIsAssembled()");
        if (assembled==QFALSE) {
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
        }
        /* zeroes out the matrix */
        ierr = QcMatZeroEntries(mat_zero[imat]);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling QcMatZeroEntries()");
    }
    return QSUCCESS;
}

