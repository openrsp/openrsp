/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

*/

#include "OpenRSP.h"

/* <function name='OpenRSPCreate' author='Bin Gao' date='2014-01-28'>
     Creates the OpenRSP context
     <param name='open_rsp' direction='inout'>The OpenRSP context</param>
     <param name='num_atoms' direction='in'>
       Number of atoms
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPCreate(OpenRSP *open_rsp, const QInt num_atoms)
{
    open_rsp->assembled = QFALSE;
    open_rsp->rsp_pert = NULL;
    /*open_rsp->elec_wav = NULL;*/
    /*open_rsp->elec_wav_type = ELEC_AO_D_MATRIX;*/
    open_rsp->overlap = NULL;
    open_rsp->one_oper = NULL;
    open_rsp->two_oper = NULL;
    open_rsp->xc_fun = NULL;
    open_rsp->zero_oper = NULL;
    open_rsp->rsp_solver = NULL;
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    open_rsp->num_atoms = num_atoms;
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
    if (open_rsp->zero_oper!=NULL) {
        ierr = RSPZeroOperAssemble(open_rsp->zero_oper,
                                      open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperAssemble()");
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
     <param name='fp_rsp' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPWrite(const OpenRSP *open_rsp, FILE *fp_rsp)
{
    QErrorCode ierr;  /* error information */
    fprintf(fp_rsp,
            "\nOpenRSP library compiled at %s, %s\n",
            __TIME__,
            __DATE__);
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
    if (open_rsp->zero_oper!=NULL) {
        fprintf(fp_rsp, "OpenRSPWrite>> nuclear Hamiltonian\n");
        ierr = RSPZeroOperWrite(open_rsp->zero_oper, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperWrite()");
    }
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverWrite(open_rsp->rsp_solver, fp_rsp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverWrite()");
    }
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    fprintf(fp_rsp,
            "OpenRSPWrite>> number of atoms %"QINT_FMT"\n",
            open_rsp->num_atoms);
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
    if (open_rsp->zero_oper!=NULL) {
        ierr = RSPZeroOperDestroy(&open_rsp->zero_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperDestroy()");
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
                                   void *user_ctx,
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
                             void *user_ctx,
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
                             void *user_ctx,
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
                             void *user_ctx,
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
                           void *user_ctx,
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

/* <function name='OpenRSPAddZeroOper' author='Bin Gao' date='2015-02-12'>
     Add a zero-electron operator to the Hamiltonian
     <param name='open_rsp' direction='inout'>
       The context of response theory calculations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act on the
       zero-electron operator
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
     <param name='get_zero_oper_contrib' direction='in'>
       User-specified callback function to calculate contributions from the
       zero-electron operator
     </param>
     <return>Error information</return>
   </function> */
QErrorCode OpenRSPAddZeroOper(OpenRSP *open_rsp,
                              const QInt num_pert_lab,
                              const QcPertInt *pert_labels,
                              const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                              void *user_ctx,
#endif

                              const GetZeroOperContrib get_zero_oper_contrib)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of zero-electron operators */
    if (open_rsp->zero_oper==NULL) {
        ierr = RSPZeroOperCreate(&open_rsp->zero_oper,
                                 num_pert_lab,
                                 pert_labels,
                                 pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                 user_ctx,
#endif
                                 get_zero_oper_contrib);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperCreate()");
    }
    /* adds the zero-electron operator to the linked list */
    else {
        ierr = RSPZeroOperAdd(open_rsp->zero_oper,
                              num_pert_lab,
                              pert_labels,
                              pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                              user_ctx,
#endif
                              get_zero_oper_contrib);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperAdd()");
    }
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
                                     void *user_ctx,
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

/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
void OpenRSPGetRSPFun_f(const QInt num_atoms,
                        const QInt num_props,
                        const QInt *len_tuple,
                        const QcPertInt *pert_tuple,
                        const QInt *num_freq_configs,
                        const QReal *pert_freqs,
                        const QInt *kn_rules,
                        const QcMat *ref_ham,
                        const QcMat *ref_overlap,
                        const QcMat *ref_state,
                        RSPSolver *rsp_solver,
                        RSPZeroOper *zero_oper,
                        RSPOverlap *overlap,
                        RSPOneOper *one_oper,
                        RSPTwoOper *two_oper,
                        RSPXCFun *xc_fun,
                        const QInt r_flag,
                        const QReal write_threshold,
                        const QInt size_rsp_funs,
                        QReal *rsp_funs);

/*@% \brief gets the response functions for given perturbations
     \author Bin Gao
     \date 2014-07-31
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QcMat:struct]{in} ref_ham Hamiltonian of referenced state
     \param[QcMat:struct]{in} ref_state electronic state of referenced state
     \param[QcMat:struct]{in} ref_overlap overlap integral matrix of referenced state
     \param[QInt:int]{in} num_props number of properties to calculate
     \param[QInt:int]{in} len_tuple length of perturbation tuple for each property
     \param[QInt:int]{in} pert_tuple ordered list of perturbation labels
         for each property
     \param[QInt:int]{in} num_freq_configs number of different frequency
         configurations for each property
     \param[QReal:real]{in} pert_freqs complex frequencies of each perturbation label
         (except for the perturbation a) over all frequency configurations
     \param[QInt:int]{in} kn_rules number k for the kn rule for each property
     \param[QInt:int]{in} size_rsp_funs size of the response functions
     \param[QReal:real]{out} rsp_funs the response functions
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPGetRSPFun(OpenRSP *open_rsp,
                            const QcMat *ref_ham,
                            const QcMat *ref_state,
                            const QcMat *ref_overlap,
                            const QInt num_props,
                            const QInt *len_tuple,
                            const QcPertInt *pert_tuple,
                            const QInt *num_freq_configs,
                            const QReal *pert_freqs,
                            const QInt *kn_rules,
                            const QInt r_flag,
                            const QReal write_threshold,
                            const QInt size_rsp_funs,
                            QReal *rsp_funs)
{
    //QErrorCode ierr;  /* error information */
    if (open_rsp->assembled==QFALSE) {
        QErrorExit(FILE_AND_LINE, "OpenRSPAssemble() should be called before calculations");
    }
    //switch (open_rsp->elec_wav_type) {
    ///* density matrix-based response theory */
    //case ELEC_AO_D_MATRIX:
        OpenRSPGetRSPFun_f(open_rsp->num_atoms,
                           num_props,
                           len_tuple,
                           pert_tuple,
                           num_freq_configs,
                           pert_freqs,
                           kn_rules,
                           ref_ham,
                           ref_overlap,
                           ref_state,
                           open_rsp->rsp_solver,
                           open_rsp->zero_oper,
                           open_rsp->overlap,
                           open_rsp->one_oper,
                           open_rsp->two_oper,
                           open_rsp->xc_fun,
                           //id_outp,
                           r_flag,
                           write_threshold,
                           size_rsp_funs,
                           rsp_funs);
    //    break;
    ///* molecular orbital (MO) coefficient matrix-based response theory */
    //case ELEC_MO_C_MATRIX:
    //    break;
    ///* couple cluster-based response theory */
    //case ELEC_COUPLED_CLUSTER:
    //    break;
    //default:
    //    printf("OpenRSPGetRSPFun>> type of (electronic) wave function %d\n",
    //           open_rsp->elec_wav_type);
    //    QErrorExit(FILE_AND_LINE, "invalid type of (electronic) wave function");
    //}
    return QSUCCESS;
}

/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
void OpenRSPGetResidue_f(const QInt num_atoms,
                         const QInt num_props,
                         const QInt *len_tuple,
                         const QcPertInt *pert_tuple,
                         const QInt *residue_num_pert,
                         const QInt *residue_idx_pert,
                         const QInt *num_freq_configs,
                         const QReal *pert_freqs,
                         const QInt *kn_rules,
                         const QcMat *ref_ham,
                         const QcMat *ref_overlap,
                         const QcMat *ref_state,
                         const QInt order_residue,
                         const QInt num_excit,
                         const QReal *excit_energy,
                         QcMat *eigen_vector[],
                         RSPSolver *rsp_solver,
                         //RSPZeroOper *zero_oper,
                         RSPOverlap *overlap,
                         RSPOneOper *one_oper,
                         RSPTwoOper *two_oper,
                         RSPXCFun *xc_fun,
                         const QInt r_flag,
                         const QReal write_threshold,
                         const QInt size_residues,
                         QReal *residues);

/*@% \brief gets the residues for given perturbations
     \author Bin Gao
     \date 2014-07-31
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QcMat:struct]{in} ref_ham Hamiltonian of referenced state
     \param[QcMat:struct]{in} ref_state electronic state of referenced state
     \param[QcMat:struct]{in} ref_overlap overlap integral matrix of referenced state
     \param[QInt:int]{in} order_residue order of residues, that is also the length of
         each excitation tuple
     \param[QInt:int]{in} num_excit number of excitation tuples that will be used for
         residue calculations
     \param[QReal:real]{in} excit_energy excitation energies of all tuples, size is
         ``order_residue`` :math:`\times` ``num_excit``, and arranged
         as ``[num_excit][order_residue]``; that is, there will be
         ``order_residue`` frequencies of perturbation labels (or sums
         of frequencies of perturbation labels) respectively equal to
         the ``order_residue`` excitation energies per tuple
         ``excit_energy[i][:]`` (``i`` runs from ``0`` to ``num_excit-1``)
     \param[QcMat:struct]{in} eigen_vector eigenvectors (obtained from the generalized
         eigenvalue problem) of all excitation tuples, size is ``order_residue``
         :math:`\times` ``num_excit``, and also arranged in memory
         as ``[num_excit][order_residue]`` so that each eigenvector has
         its corresponding excitation energy in ``excit_energy``
     \param[QInt:int]{in} num_props number of properties to calculate
     \param[QInt:int]{in} len_tuple length of perturbation tuple for each property
     \param[QInt:int]{in} pert_tuple ordered list of perturbation labels
         for each property
     \param[QInt:int]{in} residue_num_pert for each property and each excitation energy
         in the tuple, the number of perturbation labels whose sum of
         frequencies equals to that excitation energy, size is ``order_residue``
         :math:`\times` ``num_props``, and arragned as ``[num_props][order_residue]``;
         a negative ``residue_num_pert[i][j]`` (``i`` runs from ``0`` to
         ``num_props-1``) means that the sum of frequencies of perturbation
         labels equals to ``-excit_energy[:][j]``
     \param[QInt:int]{in} residue_idx_pert for each property and each excitation energy
         in the tuple, the indices of perturbation labels whose sum of
         frequencies equals to that excitation energy, size is
         ``sum(residue_num_pert)``, and arranged as ``[residue_num_pert]``
     \param[QInt:int]{in} num_freq_configs number of different frequency
         configurations for each property
     \param[QReal:real]{in} pert_freqs complex frequencies of each perturbation
         label (except for the perturbation a) over all frequency configurations
         and excitation tuples
     \param[QInt:int]{in} kn_rules number k for the kn rule for each property
     \param[QInt:int]{in} size_residues size of the residues
     \param[QReal:real]{out} residues the residues
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPGetResidue(OpenRSP *open_rsp,
                             const QcMat *ref_ham,
                             const QcMat *ref_state,
                             const QcMat *ref_overlap,
                             const QInt order_residue,
                             const QInt num_excit,
                             const QReal *excit_energy,
                             QcMat *eigen_vector[],
                             const QInt num_props,
                             const QInt *len_tuple,
                             const QcPertInt *pert_tuple,
                             const QInt *residue_num_pert,
                             const QInt *residue_idx_pert,
                             const QInt *num_freq_configs,
                             const QReal *pert_freqs,
                             const QInt *kn_rules,
                             const QInt r_flag,
                             const QReal write_threshold,
                             const QInt size_residues,
                             QReal *residues)
{
    //QErrorCode ierr;  /* error information */
    if (open_rsp->assembled==QFALSE) {
        QErrorExit(FILE_AND_LINE, "OpenRSPAssemble() should be invoked before any calculation");
    }
    //switch (open_rsp->elec_wav_type) {
    ///* density matrix-based response theory */
    //case ELEC_AO_D_MATRIX:
        OpenRSPGetResidue_f(open_rsp->num_atoms,
                            num_props,
                            len_tuple,
                            pert_tuple,
                            residue_num_pert,
                            residue_idx_pert,
                            num_freq_configs,
                            pert_freqs,
                            kn_rules,
                            ref_ham,
                            ref_overlap,
                            ref_state,
                            order_residue,
                            num_excit,
                            excit_energy,
                            eigen_vector,
                            open_rsp->rsp_solver,
                            //open_rsp->zero_oper,
                            open_rsp->overlap,
                            open_rsp->one_oper,
                            open_rsp->two_oper,
                            open_rsp->xc_fun,
                            //id_outp,
                            r_flag,
                            write_threshold,
                            size_residues,
                            residues);
    //    break;
    ///* molecular orbital (MO) coefficient matrix-based response theory */
    //case ELEC_MO_C_MATRIX:
    //    break;
    ///* couple cluster-based response theory */
    //case ELEC_COUPLED_CLUSTER:
    //    break;
    //default:
    //    printf("OpenRSPGetResidue>> type of (electronic) wave function %d\n",
    //           open_rsp->elec_wav_type);
    //    QErrorExit(FILE_AND_LINE, "invalid type of (electronic) wave function");
    //}
    return QSUCCESS;
}

