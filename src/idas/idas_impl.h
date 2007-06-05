/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2007-06-05 21:03:55 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file (private version) for the main IDAS solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDAS_IMPL_H
#define _IDAS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdarg.h>

#include <idas/idas.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

/* 
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */


/* Basic IDA constants */
  
#define HMAX_INV_DEFAULT RCONST(0.0) /* hmax_inv default value          */
#define MAXORD_DEFAULT   5           /* maxord default value            */
#define MXORDP1          6           /* max. number of N_Vectors in phi */
#define MXSTEP_DEFAULT   500         /* mxstep default value            */

/* itol */
#define IDA_NN               0
#define IDA_SS               1
#define IDA_SV               2
#define IDA_WF               3
#define IDA_EE               4 

/*
 * -----------------------------------------------------------------
 * Types: struct IDAMemRec, IDAMem
 * -----------------------------------------------------------------
 * The type IDAMem is type pointer to struct IDAMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */

typedef struct IDAMemRec {

  realtype ida_uround;    /* machine unit roundoff */

  /*--------------------------
    Problem Specification Data 
    --------------------------*/

  IDAResFn       ida_res;            /* F(t,y(t),y'(t))=0; the function F  */
  void          *ida_user_data;      /* user pointer passed to res         */

  int            ida_itol;           /* itol = IDA_SS, IDA_SV or IDA_WF    */
  realtype       ida_rtol;           /* relative tolerance                 */
  realtype       ida_Satol;          /* scalar absolute tolerance          */  
  N_Vector       ida_Vatol;          /* vector absolute tolerance          */  
  IDAEwtFn       ida_efun;           /* function to set ewt                */
  void          *ida_edata;          /* user pointer passed to efun        */
  
  /*-----------------------
    Quadrature Related Data 
    -----------------------*/

  booleantype    ida_quadr;

  IDAQuadRhsFn   ida_rhsQ;
  void          *ida_user_dataQ;

  booleantype    ida_errconQ;

  int            ida_itolQ;
  realtype       ida_rtolQ;
  realtype       ida_SatolQ;    /* scalar absolute tolerance for quadratures    */
  N_Vector       ida_VatolQ;    /* vector absolute tolerance for quadratures    */

  /*------------------------
    Sensitivity Related Data
    ------------------------*/

  booleantype    ida_sensi;
  int            ida_Ns;
  int            ida_ism;

  IDASensResFn   ida_resS;
  void          *ida_user_dataS;
  booleantype    ida_resSDQ;

  realtype      *ida_p;
  realtype      *ida_pbar;
  int           *ida_plist;
  int            ida_DQtype;
  realtype       ida_DQrhomax;

  booleantype    ida_errconS;        /* TRUE if sensitivities in err. control  */

  int            ida_itolS;
  realtype       ida_rtolS;        /* relative tolerance for sensitivities   */
  realtype       *ida_SatolS;      /* scalar absolute tolerances for sensi.  */
  N_Vector       *ida_VatolS;      /* vector absolute tolerances for sensi.  */

  /*-----------------------------------------------
    Divided differences array and associated arrays
    -----------------------------------------------*/

  N_Vector ida_phi[MXORDP1];   /* phi = (maxord+1) arrays of divided differences */

  realtype ida_psi[MXORDP1];   /* differences in t (sums of recent step sizes)   */
  realtype ida_alpha[MXORDP1]; /* ratios of current stepsize to psi values       */
  realtype ida_beta[MXORDP1];  /* ratios of current to previous product of psi's */
  realtype ida_sigma[MXORDP1]; /* product successive alpha values and factorial  */
  realtype ida_gamma[MXORDP1]; /* sum of reciprocals of psi values               */

  /*-------------------------
    N_Vectors for integration
    -------------------------*/

  N_Vector ida_ewt;         /* error weight vector                           */
  N_Vector ida_yy;          /* work space for y vector (= user's yret)       */
  N_Vector ida_yp;          /* work space for y' vector (= user's ypret)     */
  N_Vector ida_delta;       /* residual vector                               */
  N_Vector ida_id;          /* bit vector for diff./algebraic components     */
  N_Vector ida_constraints; /* vector of inequality constraint options       */
  N_Vector ida_savres;      /* saved residual vector (= tempv1)              */
  N_Vector ida_ee;          /* accumulated corrections to y vector, but
			       set equal to estimated local errors upon
			       successful return                             */
  N_Vector ida_mm;          /* mask vector in constraints tests (= tempv2)   */
  N_Vector ida_tempv1;      /* work space vector                             */
  N_Vector ida_tempv2;      /* work space vector                             */
  N_Vector ida_ynew;        /* work vector for y in IDACalcIC (= tempv2)     */
  N_Vector ida_ypnew;       /* work vector for yp in IDACalcIC (= ee)        */
  N_Vector ida_delnew;      /* work vector for delta in IDACalcIC (= phi[2]) */
  N_Vector ida_dtemp;       /* work vector in IDACalcIC (= phi[3])           */

  /*----------------------------
    Quadrature Related N_Vectors 
    ----------------------------*/

  N_Vector ida_phiQ[MXORDP1];
  N_Vector ida_yyQ;
  N_Vector ida_ypQ;
  N_Vector ida_ewtQ;
  N_Vector ida_eeQ;

  /*---------------------------
    Sensitivity Related Vectors 
    ---------------------------*/

  N_Vector *ida_phiS[MXORDP1];
  N_Vector *ida_ewtS;

  N_Vector *ida_yS0;        /* initial yS vector (user-supplied)            */
  N_Vector *ida_ypS0;       /* initial yS' vector (user-supplied)           */

  N_Vector *ida_eeS;        /* cumulative sensitivity corrections           */

  N_Vector *ida_yyS;        /* allocated and used for:                      */
  N_Vector *ida_ypS;        /*                 ism = SIMULTANEOUS           */
  N_Vector *ida_deltaS;     /*                 ism = STAGGERED              */

  N_Vector ida_tmpS1;       /* work space vectors  | tmpS1 = tempv1         */
  N_Vector ida_tmpS2;       /* for resS            | tmpS2 = tempv2         */
  N_Vector ida_tmpS3;       /*                     | tmpS3 = allocated      */    

  /*------------------------------ 
    Variables for use by IDACalcIC
    ------------------------------*/

  realtype ida_t0;          /* initial t                                     */
  N_Vector ida_yy0;         /* initial y vector (user-supplied).             */
  N_Vector ida_yp0;         /* initial y' vector (user-supplied).            */

  int ida_icopt;            /* IC calculation user option                    */
  booleantype ida_lsoff;    /* IC calculation linesearch turnoff option      */
  int ida_maxnh;            /* max. number of h tries in IC calculation      */
  int ida_maxnj;            /* max. number of J tries in IC calculation      */
  int ida_maxnit;           /* max. number of Netwon iterations in IC calc.  */
  int ida_nbacktr;          /* number of IC linesearch backtrack operations  */
  int ida_sysindex;         /* computed system index (0 or 1)                */
  realtype ida_epiccon;     /* IC nonlinear convergence test constant        */
  realtype ida_steptol;     /* minimum Newton step size in IC calculation    */
  realtype ida_tscale;      /* time scale factor = abs(tout1 - t0)           */

  /* Tstop information */

  booleantype ida_tstopset;
  realtype ida_tstop;

  /* Step Data */

  int ida_kk;        /* current BDF method order                              */
  int ida_knew;      /* order for next step from order decrease decision      */
  int ida_phase;     /* flag to trigger step doubling in first few steps      */
  int ida_ns;        /* counts steps at fixed stepsize and order              */

  realtype ida_hin;      /* initial step                                      */
  realtype ida_hh;       /* current step size h                               */
  realtype ida_rr;       /* rr = hnext / hused                                */
  realtype ida_tn;       /* current internal value of t                       */
  realtype ida_tretlast; /* value of tret previously returned by IDASolve     */
  realtype ida_cj;       /* current value of scalar (-alphas/hh) in Jacobian  */
  realtype ida_cjlast;   /* cj value saved from last successful step          */
  realtype ida_cjold;    /* cj value saved from last call to lsetup           */
  realtype ida_cjratio;  /* ratio of cj values: cj/cjold                      */
  realtype ida_ss;       /* scalar used in Newton iteration convergence test  */
  realtype ida_epsNewt;  /* test constant in Newton convergence test          */
  realtype ida_epcon;    /* coeficient of the Newton covergence test          */
  realtype ida_toldel;   /* tolerance in direct test on Newton corrections    */

  realtype ida_ssS;      /* scalar ss for staggered sensitivities             */

  /*------
    Limits
    ------*/

  int ida_maxncf;        /* max numer of convergence failures                 */
  int ida_maxcor;        /* max number of Newton corrections                  */
  int ida_maxnef;        /* max number of error test failures                 */

  int ida_maxord;        /* max value of method order k:                      */
  int ida_maxord_alloc;  /* value of maxord used when allocating memory       */
  long int ida_mxstep;   /* max number of internal steps for one user call    */
  realtype ida_hmax_inv; /* inverse of max. step size hmax (default = 0.0)    */

  int ida_maxcorS;       /* max number of Newton corrections for sensitivity
			    systems (staggered method)                        */

  /*--------
    Counters
    --------*/

  long int ida_nst;      /* number of internal steps taken                    */

  long int ida_nre;      /* number of function (res) calls                    */
  long int ida_nrQe;
  long int ida_nrSe;
  long int ida_nreS;

  long int ida_ncfn;     /* number of corrector convergence failures          */
  long int ida_ncfnQ;
  long int ida_ncfnS;

  long int ida_netf;     /* number of error test failures                     */
  long int ida_netfQ;
  long int ida_netfS;

  long int ida_nni;      /* number of Newton iterations performed             */
  long int ida_nniS;

  long int ida_nsetups;  /* number of lsetup calls                            */
  long int ida_nsetupsS;
  
  /*---------------------------
    Space requirements for IDAS
    ---------------------------*/

  long int ida_lrw1;     /* no. of realtype words in 1 N_Vector               */
  long int ida_liw1;     /* no. of integer words in 1 N_Vector                */
  long int ida_lrw1Q;
  long int ida_liw1Q;
  long int ida_lrw;      /* number of realtype words in IDA work vectors      */
  long int ida_liw;      /* no. of integer words in IDA work vectors          */


  /*-------------------------------------------
    Error handler function and error ouput file 
    -------------------------------------------*/

  IDAErrHandlerFn ida_ehfun; /* Error messages are handled by ehfun           */
  void *ida_eh_data;         /* dats pointer passed to ehfun                  */
  FILE *ida_errfp;           /* IDA error messages are sent to errfp          */

  /* Flags to verify correct calling sequence */
    
  booleantype ida_SetupDone;     /* set to FALSE by IDAInit and IDAReInit
				    set to TRUE by IDACalcIC or IDASolve      */

  booleantype ida_VatolMallocDone;
  booleantype ida_constraintsMallocDone;
  booleantype ida_idMallocDone;

  booleantype ida_MallocDone;    /* set to FALSE by IDACreate
				    set to TRUE by IDAInit
				    tested by IDAReInit and IDASolve          */

  booleantype ida_VatolQMallocDone;
  booleantype ida_quadMallocDone;

  booleantype ida_VatolSMallocDone;
  booleantype ida_SatolSMallocDone;
  booleantype ida_sensMallocDone;


  /*------------------
    Linear Solver Data
    ------------------*/

  /* Linear Solver functions to be called */

  int (*ida_linit)(struct IDAMemRec *idamem);

  int (*ida_lsetup)(struct IDAMemRec *idamem, N_Vector yyp, 
		    N_Vector ypp, N_Vector resp, 
		    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3); 

  int (*ida_lsolve)(struct IDAMemRec *idamem, N_Vector b, N_Vector weight,
		    N_Vector ycur, N_Vector ypcur, N_Vector rescur);

  int (*ida_lperf)(struct IDAMemRec *idamem, int perftask);

  int (*ida_lfree)(struct IDAMemRec *idamem);

  /* Linear Solver specific memory */

  void *ida_lmem;           

  /* Flag to request a call to the setup routine */
  
  booleantype ida_forceSetup;

  /* Flag to indicate successful ida_linit call */

  booleantype ida_linitOK;

  /*------------
    Saved Values
    ------------*/

  booleantype    ida_setupNonNull;   /* Does setup do something?              */
  booleantype    ida_constraintsSet; /* constraints vector present            */
  booleantype    ida_suppressalg;    /* true means suppress algebraic vars
					in local error tests                  */
  int ida_kused;         /* method order used on last successful step         */
  realtype ida_h0u;      /* actual initial stepsize                           */
  realtype ida_hused;    /* step size used on last successful step            */
  realtype ida_tolsf;    /* tolerance scale factor (saved value)              */

  /*----------------
    Rootfinding Data
    ----------------*/

  IDARootFn ida_gfun;    /* Function g for roots sought                     */
  int ida_nrtfn;         /* number of components of g                       */
  int *ida_iroots;       /* array for root information                      */
  int *ida_rootdir;      /* array specifying direction of zero-crossing     */
  realtype ida_tlo;      /* nearest endpoint of interval in root search     */
  realtype ida_thi;      /* farthest endpoint of interval in root search    */
  realtype ida_trout;    /* t return value from rootfinder routine          */
  realtype *ida_glo;     /* saved array of g values at t = tlo              */
  realtype *ida_ghi;     /* saved array of g values at t = thi              */
  realtype *ida_grout;   /* array of g values at t = trout                  */
  realtype ida_toutc;    /* copy of tout (if NORMAL mode)                   */
  realtype ida_ttol;     /* tolerance on root location                      */
  int ida_taskc;         /* copy of parameter itask                         */
  int ida_irfnd;         /* flag showing whether last step had a root       */
  long int ida_nge;      /* counter for g evaluations                       */


} *IDAMem;

/* 
 * =================================================================
 *   A D J O I N T   M O D U L E    M E M O R Y    B L O C K
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Forward references for pointers to various structures 
 * -----------------------------------------------------------------
 */

typedef struct IDAadjMemRec *IDAadjMem;
typedef struct CkpntMemRec *CkpntMem;
typedef struct DtpntMemRec *DtpntMem;

/*
 * -----------------------------------------------------------------
 * Types for functions provided by an interpolation module
 * -----------------------------------------------------------------
 * IDAAGetYFn:   Function type for a function that returns the 
 *               interpolated forward solution.
 * IDAAStorePnt: Function type for a function that stores a new
 *               point in the structure d
 * -----------------------------------------------------------------
 */

typedef int (*IDAAGetYFn)(IDAadjMem IDAADJ_mem, realtype t, N_Vector y, N_Vector yp);
typedef int (*IDAAStorePntFn)(IDAMem IDA_mem, DtpntMem d);

/*
 * -----------------------------------------------------------------
 * Types : struct CkpntMemRec, CkpntMem
 * -----------------------------------------------------------------
 * The type CkpntMem is type pointer to struct CkpntMemRec.
 * This structure contains fields to store all information at a
 * check point that is needed to 'hot' start IDAS.
 * -----------------------------------------------------------------
 */

struct CkpntMemRec {

  /* Integration limits */
  realtype ck_t0;
  realtype ck_t1;
    
  /* Nordsieck History Array */
  N_Vector ck_phi[MXORDP1];

  /* Nordsieck History Array for quadratures */
  N_Vector ck_phiQ[MXORDP1];
    
  /* Do we need to carry quadratures? */
  booleantype ck_quadr;

  realtype ck_psi[MXORDP1];
  realtype ck_alpha[MXORDP1];
  realtype ck_beta[MXORDP1];
  realtype ck_sigma[MXORDP1];
  realtype ck_gamma[MXORDP1];
    
  /* Step data */
  long int     ck_nst;
  long int     ck_ns;
  int          ck_kk;
  int          ck_kused;
  int          ck_knew;
  int          ck_phase;
    
  realtype     ck_hh;
  realtype     ck_hused;
  realtype     ck_rr;
  realtype     ck_cj;
  realtype     ck_cjlast;
  realtype     ck_cjold;
  realtype     ck_cjratio;
        
  /* Pointer to next structure in list */
  struct CkpntMemRec *ck_next;
    
};

/*
 * -----------------------------------------------------------------
 * Type : struct DtpntMemRec
 * -----------------------------------------------------------------
 * This structure contains fields to store all information at a
 * data point that is needed to interpolate solution of forward
 * simulations. Its content field is interpType-dependent.
 * -----------------------------------------------------------------
 */
  
struct DtpntMemRec {
  realtype t;    /* time */
  void *content; /* interpType-dependent content */
};

/* Data for cubic Hermite interpolation */
typedef struct HermiteDataMemRec {
  N_Vector y;
  N_Vector yd;
} *HermiteDataMem;

/* Data for polynomial interpolation */
typedef struct PolynomialDataMemRec {
  N_Vector y;
  int order;
} *PolynomialDataMem;


/*
 * -----------------------------------------------------------------
 * Type : struct IDAadjMemRec
 * -----------------------------------------------------------------
 * The type IDAadjMem is type pointer to struct IDAadjMemRec.
 * This structure contins fields to store all information
 * necessary for adjoint sensitivity analysis.
 * -----------------------------------------------------------------
 */

struct IDAadjMemRec {
    
  /* IDAS memory for forward runs */
  struct IDAMemRec *IDA_mem;
    
  /* IDAS memory for backward run */
  struct IDAMemRec *IDAB_mem;
    
  /* Storage for check point information */
  struct CkpntMemRec *ck_mem;

  /* Interpolation type */
  int ca_interpType;

  /* Storage for data from forward runs */
  struct DtpntMemRec **dt_mem;
    
  /* Functions set by the interpolation module */
  IDAAStorePntFn ia_storePnt; /* store a new interpolation point */
  IDAAGetYFn     ia_getY;     /* interpolate forward solution    */


  /* Residual function (resB) for backward run */
  IDAResFnB ia_resB;
    
  /* Right hand side quadrature function (fQB) for backward run */
  IDAQuadRhsFnB ia_rhsQB;

  /* User user_dataB */
  void *ia_user_dataB;
    
  /* User user_dataQB */
  void *ia_user_dataQB;

  /* Memory block for a linear solver's interface to IDAA */
  void *ia_lmemB;

  /* Function to free any memory allocated by the linear solver */
  void (*ia_lfreeB)(IDAadjMem IDAADJ_mem);

  /* Memory block for a preconditioner's module interface to IDAA */ 
  void *ia_pmemB;
    
  /* Unit roundoff */
  realtype ia_uround;
    
  /* Integration interval */
  realtype ia_tinitial, ia_tfinal;
    
  /* Time at which to extract quadratures */
  realtype ia_t_for_quad;
    
  /* Number of check points */
  int ia_nckpnts;
    
  /* Number of steps between 2 check points */
  long int ia_nsteps;
    
  /* Flag to indicate that data in dt_mem is new */
  booleantype ia_newData;
    
  /* address of the check point structure for which data is available */
  struct CkpntMemRec *ia_ckpntData;
    
  /* Actual number of data points saved in current dt_mem */
  /* Commonly, np = nsteps+1                              */
  long int ia_np;
    
  /* Workspace used by the Hermite interpolation */
  N_Vector ia_Y0, ia_Y1;       /* pointers to phi[0] and phi[1] */

  /* Workspace for polynomial interpolation */
  N_Vector ia_Y[MXORDP1];      /* pointers to phi[i] */
  realtype ia_T[MXORDP1];

  /* Workspace for wrapper functions */
  N_Vector ia_ytmp, ia_yptmp;
    
};


/*
 * =================================================================
 *     I N T E R F A C E   T O    L I N E A R   S O L V E R S     
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * int (*ida_linit)(IDAMem IDA_mem);                               
 * -----------------------------------------------------------------
 * The purpose of ida_linit is to allocate memory for the          
 * solver-specific fields in the structure *(idamem->ida_lmem) and 
 * perform any needed initializations of solver-specific memory,   
 * such as counters/statistics. An (*ida_linit) should return      
 * 0 if it has successfully initialized the IDA linear solver and 
 * a non-zero value otherwise. If an error does occur, an 
 * appropriate message should be issued.
 * ----------------------------------------------------------------
 */                                                                 

/*
 * -----------------------------------------------------------------
 * int (*ida_lsetup)(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,   
 *                  N_Vector resp,                                 
 *            N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);  
 * -----------------------------------------------------------------
 * The job of ida_lsetup is to prepare the linear solver for       
 * subsequent calls to ida_lsolve. Its parameters are as follows:  
 *                                                                 
 * idamem - problem memory pointer of type IDAMem. See the big     
 *          typedef earlier in this file.                          
 *                                                                 
 *                                                                 
 * yyp   - the predicted y vector for the current IDA internal     
 *         step.                                                   
 *                                                                 
 * ypp   - the predicted y' vector for the current IDA internal    
 *         step.                                                   
 *                                                                 
 * resp  - F(tn, yyp, ypp).                                        
 *                                                                 
 * tempv1, tempv2, tempv3 - temporary N_Vectors provided for use   
 *         by ida_lsetup.                                          
 *                                                                 
 * The ida_lsetup routine should return 0 if successful,
 * a positive value for a recoverable error, and a negative value 
 * for an unrecoverable error.
 * -----------------------------------------------------------------
 */                                                                 

/*
 * -----------------------------------------------------------------
 * int (*ida_lsolve)(IDAMem IDA_mem, N_Vector b, N_Vector weight,  
 *               N_Vector ycur, N_Vector ypcur, N_Vector rescur);  
 * -----------------------------------------------------------------
 * ida_lsolve must solve the linear equation P x = b, where        
 * P is some approximation to the system Jacobian                  
 *                  J = (dF/dy) + cj (dF/dy')                      
 * evaluated at (tn,ycur,ypcur) and the RHS vector b is input.     
 * The N-vector ycur contains the solver's current approximation   
 * to y(tn), ypcur contains that for y'(tn), and the vector rescur 
 * contains the N-vector residual F(tn,ycur,ypcur).                
 * The solution is to be returned in the vector b. 
 *                                                                 
 * The ida_lsolve routine should return 0 if successful,
 * a positive value for a recoverable error, and a negative value 
 * for an unrecoverable error.
 * -----------------------------------------------------------------
 */                                                                 

/*
 * -----------------------------------------------------------------
 * int (*ida_lperf)(IDAMem IDA_mem, int perftask);                 
 * -----------------------------------------------------------------
 * ida_lperf is called two places in IDAS where linear solver       
 * performance data is required by IDAS. For perftask = 0, an       
 * initialization of performance variables is performed, while for 
 * perftask = 1, the performance is evaluated.                     
 * -----------------------------------------------------------------
 */                                                                 

/*
 * -----------------------------------------------------------------
 * int (*ida_lfree)(IDAMem IDA_mem);                               
 * -----------------------------------------------------------------
 * ida_lfree should free up any memory allocated by the linear     
 * solver. This routine is called once a problem has been          
 * completed and the linear solver is no longer needed.            
 * -----------------------------------------------------------------
 */                                                                 

/*
 * =================================================================
 *   I D A S    I N T E R N A L   F U N C T I O N S
 * =================================================================
 */

/* Prototype of internal ewtSet function */

int IDAEwtSet(N_Vector ycur, N_Vector weight, void *data);

/* High level error handler */

void IDAProcessError(IDAMem IDA_mem, 
		     int error_code, const char *module, const char *fname, 
		     const char *msgfmt, ...);

/* Prototype of internal errHandler function */

void IDAErrHandler(int error_code, const char *module, const char *function, 
		   char *msg, void *data);

/* Prototype for internal sensitivity residual DQ function */

int IDASensResDQ(int Ns, realtype t, 
		 N_Vector yy, N_Vector yp, N_Vector resval,
		 N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
		 void *user_dataS,
		 N_Vector ytemp, N_Vector yptemp, N_Vector restemp);

/*
 * =================================================================
 *    I D A S    E R R O R    M E S S A G E S
 * =================================================================
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME "t = %Lg, "
#define MSG_TIME_H "t = %Lg and h = %Lg, "
#define MSG_TIME_INT "t = %Lg is not between tcur - hu = %Lg and tcur = %Lg."
#define MSG_TIME_TOUT "tout = %Lg"
#define MSG_TIME_TSTOP "tstop = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME "t = %lg, "
#define MSG_TIME_H "t = %lg and h = %lg, "
#define MSG_TIME_INT "t = %lg is not between tcur - hu = %lg and tcur = %lg."
#define MSG_TIME_TOUT "tout = %lg"
#define MSG_TIME_TSTOP "tstop = %lg"

#else

#define MSG_TIME "t = %g, "
#define MSG_TIME_H "t = %g and h = %g, "
#define MSG_TIME_INT "t = %g is not between tcur - hu = %g and tcur = %g."
#define MSG_TIME_TOUT "tout = %g"
#define MSG_TIME_TSTOP "tstop = %g"

#endif

/* General errors */

#define MSG_MEM_FAIL       "A memory request failed."
#define MSG_NO_MEM         "ida_mem = NULL illegal."
#define MSG_NO_MALLOC      "Attempt to call before IDAMalloc."
#define MSG_BAD_NVECTOR    "A required vector operation is not implemented."

/* Initialization errors */

#define MSG_Y0_NULL        "y0 = NULL illegal."
#define MSG_YP0_NULL       "yp0 = NULL illegal."
#define MSG_BAD_ITOL       "Illegal value for itol. The legal values are IDA_SS, IDA_SV, and IDA_WF."
#define MSG_RES_NULL       "res = NULL illegal."
#define MSG_BAD_RTOL       "rtol < 0 illegal."
#define MSG_ATOL_NULL      "atol = NULL illegal."
#define MSG_BAD_ATOL       "Some atol component < 0.0 illegal."
#define MSG_ROOT_FUNC_NULL "g = NULL illegal."

#define MSG_MISSING_ID     "id = NULL but suppressalg option on."
#define MSG_NO_TOLS        "No integration tolerances have been specified."
#define MSG_FAIL_EWT       "The user-provide EwtSet function failed."
#define MSG_BAD_EWT        "Some initial ewt component = 0.0 illegal."
#define MSG_Y0_FAIL_CONSTR "y0 fails to satisfy constraints."
#define MSG_LSOLVE_NULL    "The linear solver's solve routine is NULL."
#define MSG_LINIT_FAIL     "The linear solver's init routine failed."

#define MSG_NO_QUAD        "Illegal attempt to call before calling IDAQuadInit."
#define MSG_BAD_EWTQ       "Initial ewtQ has component(s) equal to zero (illegal)."
#define MSG_BAD_ITOLQ      "Illegal value for itolQ. The legal values are IDA_SS and IDA_SV."
#define MSG_NO_TOLQ        "No integration tolerances for quadrature variables have been specified."
#define MSG_NULL_ATOLQ     "atolQ = NULL illegal."
#define MSG_BAD_RTOLQ      "rtolQ < 0 illegal."
#define MSG_BAD_ATOLQ      "atolQ has negative component(s) (illegal)."  

#define MSG_NO_SENSI       "Illegal attempt to call before calling IDASensInit."
#define MSG_BAD_EWTS       "Initial ewtS has component(s) equal to zero (illegal)."
#define MSG_BAD_ITOLS      "Illegal value for itolS. The legal values are IDA_SS, IDA_SV, and IDA_EE."
#define MSG_NULL_ATOLS     "atolS = NULL illegal."
#define MSG_BAD_RTOLS      "rtolS < 0 illegal."
#define MSG_BAD_ATOLS      "atolS has negative component(s) (illegal)."  
#define MSG_BAD_PBAR       "pbar has zero component(s) (illegal)."
#define MSG_BAD_PLIST      "plist has negative component(s) (illegal)."
#define MSG_BAD_NS         "NS <= 0 illegal."
#define MSG_NULL_YYS0      "yyS0 = NULL illegal."
#define MSG_NULL_YPS0      "ypS0 = NULL illegal."
#define MSG_BAD_ISM        "Illegal value for ism. Legal values are: IDA_SIMULTANEOUS and IDA_STAGGERED."
#define MSG_BAD_IS         "Illegal value for is."
#define MSG_NULL_DKYA      "dkyA = NULL illegal."
#define MSG_BAD_DQTYPE     "Illegal value for DQtype. Legal values are: IDA_CENTERED and IDA_FORWARD."
#define MSG_BAD_DQRHO      "DQrhomax < 0 illegal."

/* IDACalcIC error messages */

#define MSG_IC_BAD_ICOPT   "icopt has an illegal value."
#define MSG_IC_MISSING_ID  "id = NULL conflicts with icopt."
#define MSG_IC_TOO_CLOSE   "tout1 too close to t0 to attempt initial condition calculation."
#define MSG_IC_BAD_ID      "id has illegal values."
#define MSG_IC_BAD_EWT     "Some initial ewt component = 0.0 illegal."
#define MSG_IC_RES_NONREC  "The residual function failed unrecoverably. "
#define MSG_IC_RES_FAIL    "The residual function failed at the first call. "
#define MSG_IC_SETUP_FAIL  "The linear solver setup failed unrecoverably."
#define MSG_IC_SOLVE_FAIL  "The linear solver solve failed unrecoverably."
#define MSG_IC_NO_RECOVERY "The residual routine or the linear setup or solve routine had a recoverable error, but IDACalcIC was unable to recover."
#define MSG_IC_FAIL_CONSTR "Unable to satisfy the inequality constraints."
#define MSG_IC_FAILED_LINS "The linesearch algorithm failed with too small a step."
#define MSG_IC_CONV_FAILED "Newton/Linesearch algorithm failed to converge."

/* IDASolve error messages */

#define MSG_YRET_NULL      "yret = NULL illegal."
#define MSG_YPRET_NULL     "ypret = NULL illegal."
#define MSG_TRET_NULL      "tret = NULL illegal."
#define MSG_BAD_ITASK      "itask has an illegal value."
#define MSG_TOO_CLOSE      "tout too close to t0 to start integration."
#define MSG_BAD_HINIT      "Initial step is not towards tout."
#define MSG_BAD_TSTOP      "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME "in the direction of integration."
#define MSG_BAD_INIT_ROOT  "Root found at and very near initial t."
#define MSG_CLOSE_ROOTS    "Root found at and very near " MSG_TIME "."
#define MSG_MAX_STEPS      "At " MSG_TIME ", mxstep steps taken before reaching tout." 
#define MSG_EWT_NOW_FAIL   "At " MSG_TIME "the user-provide EwtSet function failed."
#define MSG_EWT_NOW_BAD    "At " MSG_TIME "some ewt component has become <= 0.0."
#define MSG_TOO_MUCH_ACC   "At " MSG_TIME "too much accuracy requested."

#define MSG_BAD_T          "Illegal value for t." MSG_TIME_INT
#define MSG_BAD_TOUT       "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration."

#define MSG_BAD_K        "Illegal value for k."
#define MSG_NULL_DKY     "dky = NULL illegal."
#define MSG_NULL_DKYP    "dkyp = NULL illegal."

#define MSG_ERR_FAILS      "At " MSG_TIME_H "the error test failed repeatedly or with |h| = hmin."
#define MSG_CONV_FAILS     "At " MSG_TIME_H "the corrector convergence failed repeatedly or with |h| = hmin."
#define MSG_SETUP_FAILED   "At " MSG_TIME "the linear solver setup failed unrecoverably."
#define MSG_SOLVE_FAILED   "At " MSG_TIME "the linear solver solve failed unrecoverably."
#define MSG_REP_RES_ERR    "At " MSG_TIME "repeated recoverable residual errors."
#define MSG_RES_NONRECOV   "At " MSG_TIME "the residual function failed unrecoverably."
#define MSG_FAILED_CONSTR  "At " MSG_TIME "unable to satisfy inequality constraints."
#define MSG_RTFUNC_FAILED  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable manner."
#define MSG_NO_ROOT        "Rootfinding was not initialized."

#define MSG_EWTQ_NOW_BAD "At " MSG_TIME ", a component of ewtQ has become <= 0."
#define MSG_QRHSFUNC_FAILED "At " MSG_TIME ", the quadrature right-hand side routine failed in an unrecoverable manner."
#define MSG_QRHSFUNC_UNREC "At " MSG_TIME ", the quadrature right-hand side failed in a recoverable manner, but no recovery is possible."
#define MSG_QRHSFUNC_REPTD "At " MSG_TIME "repeated recoverable quadrature right-hand side function errors."
#define MSG_QRHSFUNC_FIRST "The quadrature right-hand side routine failed at the first call."

#define MSG_NULL_P "p = NULL when using internal DQ for sensitivity residual illegal."
#define MSG_EWTS_NOW_BAD "At " MSG_TIME ", a component of ewtS has become <= 0."
#define MSG_SRHSFUNC_FAILED "At " MSG_TIME ", the sensitivity residual routine failed in an unrecoverable manner."
#define MSG_SRHSFUNC_UNREC "At " MSG_TIME ", the sensitivity residual failed in a recoverable manner, but no recovery is possible."
#define MSG_SRHSFUNC_REPTD "At " MSG_TIME "repeated recoverable sensitivity residual function errors."

/* IDASet* / IDAGet* error messages */

#define MSG_NEG_MAXORD     "maxord<=0 illegal."
#define MSG_BAD_MAXORD     "Illegal attempt to increase maximum order."
#define MSG_NEG_MXSTEPS    "mxsteps < 0 illegal."
#define MSG_NEG_HMAX       "hmax < 0 illegal."
#define MSG_NEG_EPCON      "epcon < 0.0 illegal."
#define MSG_BAD_CONSTR     "Illegal values in constraints vector."
#define MSG_BAD_EPICCON    "epiccon < 0.0 illegal."
#define MSG_BAD_MAXNH      "maxnh < 0 illegal."
#define MSG_BAD_MAXNJ      "maxnj < 0 illegal."
#define MSG_BAD_MAXNIT     "maxnit < 0 illegal."
#define MSG_BAD_STEPTOL    "steptol < 0.0 illegal."

#define MSG_TOO_LATE       "IDAGetConsistentIC can only be called before IDASolve."

/*
 * =================================================================
 *    I D A A    E R R O R    M E S S A G E S
 * =================================================================
 */

#define MSGAM_NULL_IDAMEM  "ida_mem = NULL illegal."
#define MSGAM_NULL_CAMEM   "idaadj_mem = NULL illegal."
#define MSGAM_BAD_STEPS    "Steps nonpositive illegal."
#define MSGAM_MEM_FAIL     "A memory request failed."
#define MSGAM_BAD_INTERP   "Illegal value for interp."
#define MSGAM_BAD_ITASKB   "Illegal value for itaskB. Legal values are IDA_NORMAL and IDA_ONE_STEP."
#define MSGAM_BAD_TB0      "The initial time tB0 is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_TBOUT    "The final time tBout is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_T        "Bad t for interpolation."
#define MSGAM_WRONG_INTERP "This function cannot be called for the specified interp type."


#ifdef __cplusplus
}
#endif

#endif
