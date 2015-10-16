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

#include "OpenRSP.h"

/* <function name='OpenRSPCreate' author='Bin Gao' date='2014-01-28'>
     Creates the OpenRSP context
     <param name='open_rsp' direction='inout'>The OpenRSP context</param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPCreate(OpenRSP *open_rsp)
{
    open_rsp->assembled = QFALSE;
    open_rsp->rsp_pert = NULL;
    /*open_rsp->elec_wav = NULL;*/
    /*open_rsp->elec_wav_type = ELEC_AO_D_MATRIX;*/
    open_rsp->overlap = NULL;
    open_rsp->one_oper = NULL;
    open_rsp->two_oper = NULL;
    open_rsp->xc_fun = NULL;
    open_rsp->nuc_hamilton = NULL;
    open_rsp->rsp_solver = NULL;
    return QSUCCESS;
}

/* <function name='OpenRSPAssemble' author='Bin Gao' date='2014-07-30'>
     Assembles the OpenRSP context
     <param name='open_rsp' direction='inout'>The OpenRSP context</param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPAssemble(OpenRSP *open_rsp)
{
    QErrorCode ierr;  /* error information */
    open_rsp->assembled = QFALSE;
    /* assembles host program perturbations */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertAssemble(open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertAssemble()");
    }
    else {
        QErrorExit(FILE_AND_LINE, "perturbations not set by OpenRSPSetPerturbations()");
    }
/*FIXME: to implement ierr = xxAssemble(open_rsp->elec_eom); */
    /* assembles overlap integrals */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapAssemble(open_rsp->overlap,
                                  open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapAssemble()");
    }
    /* assembles one-electron operators */
    if (open_rsp->one_oper!=NULL) {
        ierr = RSPOneOperAssemble(open_rsp->one_oper,
                                  open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperAssemble()");
    }
    /* assembles two-electron operators */
    if (open_rsp->two_oper!=NULL) {
        ierr = RSPTwoOperAssemble(open_rsp->two_oper,
                                  open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperAssemble()");
    }
    /* assembles XC functionals */
    if (open_rsp->xc_fun!=NULL) {
        ierr = RSPXCFunAssemble(open_rsp->xc_fun,
                                open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunAssemble()");
    }
    /* assembles nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton!=NULL) {
        ierr = RSPNucHamiltonAssemble(open_rsp->nuc_hamilton,
                                      open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonAssemble()");
    }
    /* assembles linear response equation solver */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverAssemble(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverAssemble()");
    }
    else {
        QErrorExit(FILE_AND_LINE, "solver not set by OpenRSPSetSolver()");
    }
    open_rsp->assembled = QTRUE;
    return QSUCCESS;
}

/* <function name='OpenRSPWrite' author='Bin Gao' date='2014-07-30'>
     Writes the OpenRSP context
     <param name='open_rsp' direction='in'>The OpenRSP context</param>
     <param name='file_name' direction='in'>File to write the context</param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPWrite(const OpenRSP *open_rsp, const QChar *file_name)
{
    FILE *fp_rsp;     /* file pointer */
    QErrorCode ierr;  /* error information */
    /* opens the file */
    fp_rsp = fopen(file_name, "a");
    if (fp_rsp==NULL) {
        printf("OpenRSPWrite>> file: %s\n", file_name);
        QErrorExit(FILE_AND_LINE, "failed to open the file in appending mode");
    }
    fprintf(fp_rsp, "\nOpenRSP library compiled at %s, %s\n", __TIME__, __DATE__);
    /* context of the (electronic) wave function */
    /*FIXME: ierr = xxWrite(open_rsp->elec_eom); */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertWrite(open_rsp->rsp_pert, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertWrite()");
    }
    if (open_rsp->overlap!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> overlap integrals\n");
        ierr = RSPOverlapWrite(open_rsp->overlap, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapWrite()");
    }
    if (open_rsp->one_oper!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> linked list of one-electron operators\n");
        ierr = RSPOneOperWrite(open_rsp->one_oper, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperWrite()");
    }
    if (open_rsp->two_oper!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> linked list of two-electron operators\n");
        ierr = RSPTwoOperWrite(open_rsp->two_oper, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperWrite()");
    }
    if (open_rsp->xc_fun!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> linked list of XC functionals\n");
        ierr = RSPXCFunWrite(open_rsp->xc_fun, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunWrite()");
    }
    if (open_rsp->nuc_hamilton!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> nuclear Hamiltonian\n");
        ierr = RSPNucHamiltonWrite(open_rsp->nuc_hamilton, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonWrite()");
    }
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverWrite(open_rsp->rsp_solver, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverWrite()");
    }
    /* closes the file */
    fclose(fp_rsp);
    return QSUCCESS;
}

/* <function name='OpenRSPDestroy' author='Bin Gao' date='2014-01-28'>
     Destroys the OpenRSP context
     <param name='open_rsp' direction='inout'>The OpenRSP context</param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPDestroy(OpenRSP *open_rsp)
{
    QErrorCode ierr;  /* error information */
    open_rsp->assembled = QFALSE;
//    if (open_rsp->elec_eom!=NULL) {
///*FIXME: to implement ierr = xxDestroy(open_rsp->elec_eom); */
//        free(open_rsp->elec_eom);
//        open_rsp->elec_eom = NULL;
//    }
    /* destroys the context of all perturbations involved in calculations */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertDestroy(open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertDestroy()");
        free(open_rsp->rsp_pert);
        open_rsp->rsp_pert = NULL;
    }
    /* destroys the context of overlap integrals */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapDestroy(open_rsp->overlap);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapDestroy()");
        free(open_rsp->overlap);
        open_rsp->overlap = NULL;
    }
    /* destroys the linked list of one-electron operators */
    if (open_rsp->one_oper!=NULL) {
        ierr = RSPOneOperDestroy(&open_rsp->one_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperDestroy()");
    }
    /* destroys the linked list of two-electron operators */
    if (open_rsp->two_oper!=NULL) {
        ierr = RSPTwoOperDestroy(&open_rsp->two_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperDestroy()");
    }
    /* destroys the linked list of exchange-correlation functionals */
    if (open_rsp->xc_fun!=NULL) {
        ierr = RSPXCFunDestroy(&open_rsp->xc_fun);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunDestroy()");
    }
    /* destroys the context of nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton!=NULL) {
        ierr = RSPNucHamiltonDestroy(open_rsp->nuc_hamilton);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonDestroy()");
        free(open_rsp->nuc_hamilton);
        open_rsp->nuc_hamilton = NULL;
    }
    /* destroys the context of linear response equation sovler */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverDestroy(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverDestroy()");
        free(open_rsp->rsp_solver);
        open_rsp->rsp_solver = NULL;
    }
    return QSUCCESS;
}

/* <function name='OpenRSPSetPerturbations' author='Bin Gao' date='2015-06-29'>
     Sets all perturbations involved in response theory calculations
     <param name='open_rsp' direction='inout'>The OpenRSP context</param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels involved in calculations
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='pert_num_comps' direction='in'>
       Number of components of a perturbation described by exactly one of
       the above different labels, up to the allowed maximal order, size
       is therefore the sum of <pert_max_orders>
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_pert_concatenation' direction='in'>
       User specified function for getting the ranks of components of
       sub-perturbation tuples (with the same perturbation label) for given
       components of the corresponding concatenated perturbation tuple
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPSetPerturbations(OpenRSP *open_rsp,
                                   const QInt num_pert_lab,
                                   const QcPertInt *pert_labels,
                                   const QInt *pert_max_orders,
                                   const QInt *pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid *user_ctx,
#endif
                                   const GetPertCat get_pert_concatenation)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of all perturbations involved in calculations */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertDestroy(open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertDestroy()");
    }
    else {
        open_rsp->rsp_pert = (RSPPert *)malloc(sizeof(RSPPert));
        if (open_rsp->rsp_pert==NULL) {
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbations");
        }
    }
    ierr = RSPPertCreate(open_rsp->rsp_pert,
                         num_pert_lab,
                         pert_labels,
                         pert_max_orders,
                         pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                         user_ctx,
#endif
                         get_pert_concatenation);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertCreate()");
    return QSUCCESS;
}

/* <function name='OpenRSPSetOverlap' author='Bin Gao' date='2014-07-30'>
     Set the overlap operator
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       overlap operator
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_overlap_mat' direction='in'>
       User-specified callback function to calculate integral matrices of
       overlap operator as well as its derivatives with respect to different
       perturbations
     </param>
     <param name='get_overlap_exp' direction='in'>
       User-specified callback function to calculate expectation values of
       overlap operator as well as its derivatives with respect to different
       perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPSetOverlap(OpenRSP *open_rsp,
                             const QInt num_pert_lab,
                             const QcPertInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             const GetOverlapMat get_overlap_mat,
                             const GetOverlapExp get_overlap_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of overlap operator */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapDestroy(open_rsp->overlap);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapDestroy()");
    }
    else {
        open_rsp->overlap = (RSPOverlap *)malloc(sizeof(RSPOverlap));
        if (open_rsp->overlap==NULL) {
            QErrorExit(FILE_AND_LINE, "allocates memory for overlap");
        }
    }
    ierr = RSPOverlapCreate(open_rsp->overlap,
                            num_pert_lab,
                            pert_labels,
                            pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            user_ctx,
#endif
                            get_overlap_mat,
                            get_overlap_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapCreate()");
    return QSUCCESS;
}

/* <function name='OpenRSPAddOneOper' author='Bin Gao' date='2014-07-30'>
     Add a one-electron operator to the Hamiltonian
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       one-electron operator
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_one_oper_mat' direction='in'>
       User-specified callback function to calculate integral matrices of
       one-electron operator as well as its derivatives with respect to
       different perturbations
     </param>
     <param name='get_one_oper_exp' direction='in'>
       User-specified callback function to calculate expectation values of
       one-electron operator as well as its derivatives with respect to
       different perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPAddOneOper(OpenRSP *open_rsp,
                             const QInt num_pert_lab,
                             const QcPertInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             const GetOneOperMat get_one_oper_mat,
                             const GetOneOperExp get_one_oper_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of one-electron operators */
    if (open_rsp->one_oper==NULL) {
        ierr = RSPOneOperCreate(&open_rsp->one_oper,
                                num_pert_lab,
                                pert_labels,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                get_one_oper_mat,
                                get_one_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperCreate()");
    }
    /* adds the one-electron operator to the linked list */
    else {
        ierr = RSPOneOperAdd(open_rsp->one_oper,
                             num_pert_lab,
                             pert_labels,
                             pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             user_ctx,
#endif
                             get_one_oper_mat,
                             get_one_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperAdd()");
    }
    return QSUCCESS;
}

/* <function name='OpenRSPAddTwoOper' author='Bin Gao' date='2014-08-05'>
     Add a two-electron operator to the Hamiltonian
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       two-electron operator
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_two_oper_mat' direction='in'>
       User-specified callback function to calculate integral matrices of
       two-electron operator as well as its derivatives with respect to
       different perturbations
     </param>
     <param name='get_two_oper_exp' direction='in'>
       User-specified callback function to calculate expectation values of
       two-electron operator as well as its derivatives with respect to
       different perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPAddTwoOper(OpenRSP *open_rsp,
                             const QInt num_pert_lab,
                             const QcPertInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             const GetTwoOperMat get_two_oper_mat,
                             const GetTwoOperExp get_two_oper_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of two-electron operators */
    if (open_rsp->two_oper==NULL) {
        ierr = RSPTwoOperCreate(&open_rsp->two_oper,
                                num_pert_lab,
                                pert_labels,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                get_two_oper_mat,
                                get_two_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperCreate()");
    }
    /* adds the two-electron operator to the linked list */
    else {
        ierr = RSPTwoOperAdd(open_rsp->two_oper,
                             num_pert_lab,
                             pert_labels,
                             pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             user_ctx,
#endif
                             get_two_oper_mat,
                             get_two_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperAdd()");
    }
    return QSUCCESS;
}

/* <function name='OpenRSPAddXCFun' author='Bin Gao' date='2015-06-23'>
     Add an XC functional to the Hamiltonian
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       XC functional
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_xc_fun_mat' direction='in'>
       User-specified callback function to calculate integral matrices of
       XC functional as well as its derivatives with respect to
       different perturbations
     </param>
     <param name='get_xc_fun_exp' direction='in'>
       User-specified callback function to calculate expectation values of
       XC functional as well as its derivatives with respect to
       different perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPAddXCFun(OpenRSP *open_rsp,
                           const QInt num_pert_lab,
                           const QcPertInt *pert_labels,
                           const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                           QVoid *user_ctx,
#endif
                           const GetXCFunMat get_xc_fun_mat,
                           const GetXCFunExp get_xc_fun_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of XC functionals */
    if (open_rsp->xc_fun==NULL) {
        ierr = RSPXCFunCreate(&open_rsp->xc_fun,
                              num_pert_lab,
                              pert_labels,
                              pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                              user_ctx,
#endif
                              get_xc_fun_mat,
                              get_xc_fun_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunCreate()");
    }
    /* adds the XC functional to the linked list */
    else {
        ierr = RSPXCFunAdd(open_rsp->xc_fun,
                           num_pert_lab,
                           pert_labels,
                           pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                           user_ctx,
#endif
                           get_xc_fun_mat,
                           get_xc_fun_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunAdd()");
    }
    return QSUCCESS;
}

/* <function name='OpenRSPSetNucHamilton' author='Bin Gao' date='2015-02-12'>
     Set the context of nuclear Hamiltonian
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       nuclear Hamiltonian
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_nuc_contrib' direction='in'>
       User-specified callback function to calculate nuclear contributions
     </param>
     <param name='num_atoms' direction='in'>
       Number of atoms
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPSetNucHamilton(OpenRSP *open_rsp,
                                 const QInt num_pert_lab,
                                 const QcPertInt *pert_labels,
                                 const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                 QVoid *user_ctx,
#endif

                                 const GetNucContrib get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                 const QInt num_atoms)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton!=NULL) {
        ierr = RSPNucHamiltonDestroy(open_rsp->nuc_hamilton);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonDestroy()");
    }
    else {
        open_rsp->nuc_hamilton = (RSPNucHamilton *)malloc(sizeof(RSPNucHamilton));
        if (open_rsp->nuc_hamilton==NULL) {
            QErrorExit(FILE_AND_LINE, "allocates memory for nuclear Hamiltonian");
        }
    }
    ierr = RSPNucHamiltonCreate(open_rsp->nuc_hamilton,
                                num_pert_lab,
                                pert_labels,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                num_atoms);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonCreate()");
    return QSUCCESS;
}

/* <function name='OpenRSPSetLinearRSPSolver' author='Bin Gao' date='2014-08-06'>
     Set the context of linear response equation solver
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_linear_rsp_solution' direction='in'>
       User-specified callback function of linear response equation solver
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPSetLinearRSPSolver(OpenRSP *open_rsp,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     QVoid *user_ctx,
#endif
                                     const GetLinearRSPSolution get_linear_rsp_solution)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of response equation solver */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverDestroy(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverDestroy()");
    }
    else {
        open_rsp->rsp_solver = (RSPSolver *)malloc(sizeof(RSPSolver));
        if (open_rsp->rsp_solver==NULL) {
            QErrorExit(FILE_AND_LINE, "allocates memory for solver");
        }
    }
    ierr = RSPSolverCreate(open_rsp->rsp_solver,
#if defined(OPENRSP_C_USER_CONTEXT)
                           user_ctx,
#endif
                           get_linear_rsp_solution);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverCreate()");
    return QSUCCESS;
}

