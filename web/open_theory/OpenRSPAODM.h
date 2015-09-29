/* atomic-orbital density matrix */

/*#include "OpenRSPClass.h"*/
#include "OpenRSP.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*AODMOneOperMat)(const QInt,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const QInt,
                                QcMat*[]);
typedef QVoid (*AODMOneOperExp)(const QInt,
                                const QInt*,
                                const QInt*,
                                const QInt,
                                QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const QInt,
                                QReal*);

/* linked list of one-electron operators */
typedef struct AODMOneOper AODMOneOper;
struct AODMOneOper {
    QInt num_pert;                   /* number of different perturbation labels that
                                        can act as perturbations on the one-electron operator */
    QInt *pert_labels;               /* all the different perturbation labels */
    QInt *pert_max_orders;           /*  maximum allowed order of each perturbation (label) */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                 /* user-defined callback function context */
#endif
    AODMOneOperMat get_one_oper_mat;  /* user specified function for getting integral matrices */
    AODMOneOperExp get_one_oper_exp;  /* user specified function for getting expectation values */
    AODMOneOper *next_oper;           /* pointer to the next one-electron operator */
};

typedef struct {
    AODMOneOper *one_oper;          /* linked list of one-electron operators */
} RSPAODM;

extern QErrorCode RSPAODMCreate(OpenRSP*);
extern QErrorCode RSPAODMAddOneOper(OpenRSP*,
                                    const QInt,
                                    const QInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    QVoid*,
#endif
                                    va_list);
extern QErrorCode RSPAODMDestroy(OpenRSP*);
