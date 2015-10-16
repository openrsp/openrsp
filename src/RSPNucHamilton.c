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

#include "RSPNucHamilton.h"

/* <function name='RSPNucHamiltonCreate'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Create the context of nuclear Hamiltonian, should be called at first
     <param name='nuc_hamilton' direction='inout'>
       The context of nuclear Hamiltonian
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the nuclear Hamiltonian
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback-function context
     </param>
     <param name='get_nuc_contrib' direction='in'>
       User-specified function for calculating contribution of the
       nuclear Hamiltonian and its derivatives
     </param>
     <param name='num_atoms' direction='in'>
       Number of atoms
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPNucHamiltonCreate(RSPNucHamilton *nuc_hamilton,
                                const QInt num_pert_lab,
                                const QcPertInt *pert_labels,
                                const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                void *user_ctx,
#endif
                                const GetNucContrib get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                const QInt num_atoms)
{
    QInt ilab;  /* incremental recorders over perturbation labels */
    QInt jlab;
    if (num_pert_lab<0) {
        printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPNucHamiltonCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    nuc_hamilton->num_pert_lab = num_pert_lab;
    if (nuc_hamilton->num_pert_lab>0) {
        nuc_hamilton->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (nuc_hamilton->pert_max_orders==NULL) {
            printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        nuc_hamilton->nuc_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (nuc_hamilton->nuc_pert_orders==NULL) {
            printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on nuclear Hamiltonian");
        }
        nuc_hamilton->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (nuc_hamilton->pert_labels==NULL) {
            printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        nuc_hamilton->nuc_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (nuc_hamilton->nuc_pert_labels==NULL) {
            printf("RSPNucHamiltonCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on nuclear Hamiltonian");
        }
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPNucHamiltonCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPNucHamiltonCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPNucHamiltonCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPNucHamiltonCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            nuc_hamilton->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPNucHamiltonCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPNucHamiltonCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ilab]);
                QErrorExit(FILE_AND_LINE, "only positive order allowed");
            }
            nuc_hamilton->pert_max_orders[ilab] = pert_max_orders[ilab];
        }
    }
    else {
        nuc_hamilton->pert_max_orders = NULL;
        nuc_hamilton->nuc_pert_orders = NULL;
        nuc_hamilton->pert_labels = NULL;
        nuc_hamilton->nuc_pert_labels = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    nuc_hamilton->user_ctx = user_ctx;
#endif
    nuc_hamilton->get_nuc_contrib = get_nuc_contrib;
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    nuc_hamilton->num_atoms = num_atoms;
    return QSUCCESS;
}

/* <function name='RSPNucHamiltonAssemble'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Assembles the context of nuclear Hamiltonian
     <param name='nuc_hamilton' direction='inout'>
       The context of nuclear Hamiltonian
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPNucHamiltonAssemble(RSPNucHamilton *nuc_hamilton,
                                  const RSPPert *rsp_pert)
{
    QErrorCode ierr;  /* error information */
    if (nuc_hamilton->num_pert_lab>0 &&
        (nuc_hamilton->pert_labels==NULL || nuc_hamilton->pert_max_orders==NULL)) {
        QErrorExit(FILE_AND_LINE, "perturbations of nuclear Hamiltonian not set");
    }
    if (nuc_hamilton->get_nuc_contrib==NULL) {
        QErrorExit(FILE_AND_LINE, "callback function of nuclear Hamiltonian not set");
    }
    /* checks perturbation labels and allowed maximal orders against
       all known perturbations */
    ierr = RSPPertValidateLabelOrder(rsp_pert,
                                     nuc_hamilton->num_pert_lab,
                                     nuc_hamilton->pert_labels,
                                     nuc_hamilton->pert_max_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertValidateLabelOrder()");
    return QSUCCESS;
}

/* <function name='RSPNucHamiltonWrite'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Writes the context of nuclear Hamiltonian
     <param name='nuc_hamilton' direction='in'>
       The context of nuclear Hamiltonian
     </param>
     <param name='fp_nuc' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton *nuc_hamilton,
                               FILE *fp_nuc)
{
    QInt ilab;  /* incremental recorder over perturbation labels */
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> number of pert. labels that nuclear Hamiltonian depends on %"QINT_FMT"\n",
            nuc_hamilton->num_pert_lab);
    fprintf(fp_nuc, "RSPNucHamiltonWrite>> label           maximum-order\n");
    for (ilab=0; ilab<nuc_hamilton->num_pert_lab; ilab++) {
        fprintf(fp_nuc,
                "RSPNucHamiltonWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                nuc_hamilton->pert_labels[ilab],
                nuc_hamilton->pert_max_orders[ilab]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (nuc_hamilton->user_ctx!=NULL) {
        fprintf(fp_nuc, "RSPNucHamiltonWrite>> user-defined function context given\n");
    }
#endif
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> number of atoms %"QINT_FMT"\n",
            nuc_hamilton->num_atoms);
    return QSUCCESS;
}

/* <function name='RSPNucHamiltonGetContributions'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates contribution of the nuclear Hamiltonian
     <param name='nuc_hamilton' direction='inout'>
       The context of nuclear Hamiltonian
     </param>
     <param name='nuc_len_tuple' direction='in'>
       Length of the perturbation tuple on the nuclear Hamiltonian
     </param>
     <param name='nuc_pert_tuple' direction='in'>
       Perturbation tuple on the nuclear Hamiltonian
     </param>
     <param name='size_pert' direction='in'>
       Size of the perturbations on the nuclear Hamiltonian
     </param>
     <param name='val_nuc' direction='inout'>
       The contribution of the nuclear Hamiltonian
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPNucHamiltonGetContributions(RSPNucHamilton *nuc_hamilton,
                                          const QInt nuc_len_tuple,
                                          const QcPertInt *nuc_pert_tuple,
                                          const QInt size_pert,
                                          QReal *val_nuc)
{
    QErrorCode ierr;  /* error information */
    /* gets perturbation labels and corresponding orders out of the internal
       perturbation tuple on the nuclear Hamiltonian */
    ierr = RSPPertInternTupleToHostLabelOrder(nuc_len_tuple,
                                              nuc_pert_tuple,
                                              nuc_hamilton->num_pert_lab,
                                              nuc_hamilton->pert_labels,
                                              nuc_hamilton->pert_max_orders,
                                              &nuc_hamilton->nuc_num_pert,
                                              nuc_hamilton->nuc_pert_labels,
                                              nuc_hamilton->nuc_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    /* checks if the perturbations on the nuclear Hamiltonian
       result in zero values */
    if (nuc_hamilton->nuc_num_pert<0) return QSUCCESS;
    /* calculates contribution of nuclear Hamiltonian using the
       callback function */
    nuc_hamilton->get_nuc_contrib(nuc_hamilton->nuc_num_pert,
                                  nuc_hamilton->nuc_pert_labels,
                                  nuc_hamilton->nuc_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  nuc_hamilton->user_ctx,
#endif
                                  size_pert,
                                  val_nuc);
    return QSUCCESS;
}

/*% \brief gets the number of atoms
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{in} nuc_hamilton the context of nuclear Hamiltonian
    \param[QInt:int]{out} num_atoms number of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonGetNumAtoms(const RSPNucHamilton *nuc_hamilton,
                                     QInt *num_atoms)
{
    *num_atoms = nuc_hamilton->num_atoms;
    return QSUCCESS;
}

/* <function name='RSPNucHamiltonDestroy'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Destroys the context of nuclear Hamiltonian, should be called at the end
     <param name='nuc_hamilton' direction='inout'>
       The context of nuclear Hamiltonian
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPNucHamiltonDestroy(RSPNucHamilton *nuc_hamilton)
{
    if (nuc_hamilton->pert_max_orders!=NULL) {
        free(nuc_hamilton->pert_max_orders);
        nuc_hamilton->pert_max_orders = NULL;
    }
    if (nuc_hamilton->nuc_pert_orders!=NULL) {
        free(nuc_hamilton->nuc_pert_orders);
        nuc_hamilton->nuc_pert_orders = NULL;
    }
    if (nuc_hamilton->pert_labels!=NULL) {
        free(nuc_hamilton->pert_labels);
        nuc_hamilton->pert_labels = NULL;
    }
    if (nuc_hamilton->nuc_pert_labels!=NULL) {
        free(nuc_hamilton->nuc_pert_labels);
        nuc_hamilton->nuc_pert_labels = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    nuc_hamilton->user_ctx = NULL;
#endif
    nuc_hamilton->get_nuc_contrib = NULL;
    return QSUCCESS;
}

