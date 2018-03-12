/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the interface file for the main ARKODE integrator.
 *---------------------------------------------------------------
 * ARKODE is used to solve numerically the ordinary initial value
 * problem:
 *              M y'(t) = fe(t,y) + fi(t,y),
 *             y(t0) = y0,
 * where t0, y0 in R^N, y: R -> R^N, and fe,fi: R x R^N -> R^N 
 * are given.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_H
#define _ARKODE_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
                         ARKODE CONSTANTS
===============================================================*/

/*---------------------------------------------------------------
 Enumerations for inputs to ARKodeCreate and ARKode.
-----------------------------------------------------------------
 Symbolic constants for the iter parameter to 
 ARKodeCreate and the input parameter itask to ARKode, are given 
 below.
 
 itask: The itask input parameter to ARKode indicates the job
        of the solver for the next user step. The ARK_NORMAL
        itask is to have the solver take internal steps until
        it has reached or just passed the user specified tout
        parameter. The solver then interpolates in order to
        return an approximate value of y(tout). The ARK_ONE_STEP
        option tells the solver to just take one internal step
        and return the solution at the point reached by that 
        step.  If fixed time steps are desired, then the user 
        can run in ARK_ONE_STEP mode, and set hmin and hmax to
        be the desired time step -- this will effectively 
        disable all temporal error control, and the user may
        need to increase the iteration counts for the nonlinear 
        and linear solvers to guarantee convergence for larger 
        step sizes.
---------------------------------------------------------------*/

/* itask */
#define ARK_NORMAL         1
#define ARK_ONE_STEP       2


/* ARKODE return flags */
#define ARK_SUCCESS               0
#define ARK_TSTOP_RETURN          1
#define ARK_ROOT_RETURN           2

#define ARK_WARNING              99

#define ARK_TOO_MUCH_WORK        -1
#define ARK_TOO_MUCH_ACC         -2
#define ARK_ERR_FAILURE          -3
#define ARK_CONV_FAILURE         -4

#define ARK_LINIT_FAIL           -5
#define ARK_LSETUP_FAIL          -6
#define ARK_LSOLVE_FAIL          -7
#define ARK_RHSFUNC_FAIL         -8
#define ARK_FIRST_RHSFUNC_ERR    -9
#define ARK_REPTD_RHSFUNC_ERR    -10
#define ARK_UNREC_RHSFUNC_ERR    -11
#define ARK_RTFUNC_FAIL          -12
#define ARK_LFREE_FAIL           -13
#define ARK_MASSINIT_FAIL        -14
#define ARK_MASSSETUP_FAIL       -15
#define ARK_MASSSOLVE_FAIL       -16
#define ARK_MASSFREE_FAIL        -17
#define ARK_MASSMULT_FAIL        -18

#define ARK_MEM_FAIL             -20
#define ARK_MEM_NULL             -21
#define ARK_ILL_INPUT            -22
#define ARK_NO_MALLOC            -23
#define ARK_BAD_K                -24
#define ARK_BAD_T                -25
#define ARK_BAD_DKY              -26
#define ARK_TOO_CLOSE            -27

#define ARK_POSTPROCESS_FAIL     -28


/*===============================================================
                          FUNCTION TYPES
===============================================================*/

/*---------------------------------------------------------------
 Type : ARKRhsFn
-----------------------------------------------------------------
 The f functions which define the right hand side of the ODE
 system M*y' = fe(t,y) + fi(t,y) must have type ARKRhsFn.
 f takes as input the independent variable value t, and the
 dependent variable vector y.  It stores the result of fe(t,y) 
 or fi(t,y) in the vector ydot.  The y and ydot arguments are of 
 type N_Vector.
 (Allocation of memory for ydot is handled within ARKODE)
 The user_data parameter is the same as the user_data
 parameter set by the user through the ARKodeSetUserData routine.
 This user-supplied pointer is passed to the user's fe or fi
 function every time it is called.

 A ARKRhsFn should return 0 if successful, a negative value if
 an unrecoverable error occured, and a positive value if a 
 recoverable error (e.g. invalid y values) occured. 
 If an unrecoverable occured, the integration is halted. 
 If a recoverable error occured, then (in most cases) ARKODE
 will try to correct and retry.
---------------------------------------------------------------*/
typedef int (*ARKRhsFn)(realtype t, N_Vector y,
                        N_Vector ydot, void *user_data);

/*---------------------------------------------------------------
 Type : ARKRootFn
-----------------------------------------------------------------
 A function g, which defines a set of functions g_i(t,y) whose
 roots are sought during the integration, must have type 
 ARKRootFn. The function g takes as input the independent 
 variable value t, and the dependent variable vector y.  It 
 stores the nrtfn values g_i(t,y) in the realtype array gout.
 (Allocation of memory for gout is handled within ARKODE.)
 The user_data parameter is the same as that passed by the user
 to the ARKodeSetUserData routine.  This user-supplied pointer 
 is passed to the user's g function every time it is called.

 A ARKRootFn should return 0 if successful or a non-zero value
 if an error occured (in which case the integration will be 
 halted).
---------------------------------------------------------------*/
typedef int (*ARKRootFn)(realtype t, N_Vector y, 
                         realtype *gout, void *user_data);

/*---------------------------------------------------------------
 Type : ARKEwtFn
-----------------------------------------------------------------
 A function e, which sets the error weight vector ewt, must have
 type ARKEwtFn.  The function e takes as input the current 
 dependent variable y. It must set the vector of error weights 
 used in the WRMS norm:
 
   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]

 Typically, the vector ewt has components:
 
   ewt_i = 1 / (reltol * |y_i| + abstol_i)

 The user_data parameter is the same as that passed by the user
 to the ARKodeSetUserData routine.  This user-supplied pointer 
 is passed to the user's e function every time it is called.
 A ARKEwtFn e must return 0 if the error weight vector has been
 successfuly set and a non-zero value otherwise.
---------------------------------------------------------------*/
typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

/*---------------------------------------------------------------
 Type : ARKRwtFn
-----------------------------------------------------------------
 A function r, which sets the residual weight vector rwt, must 
 have type ARKRwtFn.  The function r takes as input the current 
 dependent variable y.  It must set the vector of residual 
 weights used in the WRMS norm:
 
   ||v||_WRMS = sqrt [ 1/N * sum ( rwt_i * v_i)^2 ]

 Typically, the vector rwt has components:
 
   rwt_i = 1 / (reltol * |(M*y)_i| + rabstol_i)

 The user_data parameter is the same as that passed by the user
 to the ARKodeSetUserData routine.  This user-supplied pointer 
 is passed to the user's r function every time it is called.
 A ARKRwtFn e must return 0 if the residual weight vector has 
 been successfuly set and a non-zero value otherwise.
---------------------------------------------------------------*/
typedef int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void *user_data);

/*---------------------------------------------------------------
 Type : ARKErrHandlerFn
-----------------------------------------------------------------
 A function eh, which handles error messages, must have type
 ARKErrHandlerFn.  The function eh takes as input the error code, 
 the name of the module reporting the error, the error message, 
 and a pointer to user data, the same as that passed to 
 ARKodeSetUserData.
 
 All error codes are negative, except ARK_WARNING which indicates 
 a warning (the solver continues).

 A ARKErrHandlerFn has no return value.
---------------------------------------------------------------*/
typedef void (*ARKErrHandlerFn)(int error_code, const char *module, 
                                const char *function, char *msg, 
                                void *user_data); 

/*---------------------------------------------------------------
 Type : ARKAdaptFn
-----------------------------------------------------------------
 A function which sets the new time step h, must have type 
 ARKAdaptFn.  The function takes as input the current dependent 
 variable y, the current time t, the last 3 step sizes h, 
 the last 3 error estimates, the method order q, and a 
 pointer to user data. The function must set the scalar step size 
 for the upcoming time step.  This value will subsequently be 
 bounded by the user-supplied values for the minimum and maximum 
 allowed time step, and the time step satisfying the explicit 
 stability restriction.  The user_data parameter is the same as 
 that passed by the user to the ARKodeSetUserData routine.  This 
 user-supplied pointer is passed to the function every time it 
 is called.

 A ARKAdaptFn must return 0 if the new time step has been
 successfuly set and a non-zero value otherwise.
---------------------------------------------------------------*/
typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1, 
                          realtype h2, realtype h3, 
                          realtype e1, realtype e2, 
                          realtype e3, int q, realtype *hnew,
                          void *user_data);

/*---------------------------------------------------------------
 Type : ARKExpStabFn
-----------------------------------------------------------------
 A function which returns the time step satisfying the stability 
 restriction for the explicit portion of the ODE.  The function 
 takes as input the current dependent variable y, the current 
 time t, and a pointer to user data. The function must set the 
 scalar step size satisfying the stability restriction for the 
 upcoming time step.  This value will subsequently be bounded by 
 the user-supplied values for the minimum and maximum allowed 
 time step, and the accuracy-based time step.  The user_data 
 parameter is the same as that passed by the user to the 
 ARKodeSetUserData routine.  This user-supplied pointer is passed
 to the function every time it is called.

 If this function is not supplied (NULL), or if it returns a
 negative time step size, then ARKode will assume that there is
 no explicit stability restriction on the time step size.

 A ARKExpStabFn must return 0 if the step size limit has been
 successfuly set and a non-zero value otherwise.
---------------------------------------------------------------*/
typedef int (*ARKExpStabFn)(N_Vector y, realtype t, 
                            realtype *hstab, void *user_data);

/*---------------------------------------------------------------
 Type : ARKVecResizeFn
-----------------------------------------------------------------
 When calling ARKodeResize, the user may specify a vector resize
 function to be used to convert any existing N_Vectors in the 
 ARKode memory structure to the new problem size.  This would 
 typically be used if there is a user-supplied N_Vector module 
 that allows dynamic resizing of the vector data structures 
 without the need to delete/allocate memory on each call.  

 The default behavior will be to delete the vector memory and 
 re-clone from the new vector; if this is the desired behavior 
 then specification of the ARKVecResizeFn is not recommended.

 The first argument, 'y', is the vector to be resized.

 The second argument, 'ytemplate', is the user-provided vector 
 with the "new" size, that may be used as a template.

 The third argument, 'user_data', is a user-provided data 
 structure to ARKodeResize, in case additional data is 
 necessary for the resize operation.

 A ARKVecResizeFn should return 0 if successful, and a nonzero 
 value if an error occurred.
---------------------------------------------------------------*/
typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate, 
                              void *user_data);

/*---------------------------------------------------------------
 Type : ARKPostProcessStepFn
-----------------------------------------------------------------
 A function that is used to process the results of each timestep
 solution, in preparation for subsequent steps.  A routine of 
 this type is designed for tasks such as inter-processor 
 communication, computation of derived quantities, etc..

 IF THIS IS USED TO MODIFY ANY OF THE ACTIVE STATE DATA, THEN ALL
 THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND STABILITY ARE 
 LOST.

 Inputs: 
   t          current time of ARKode solution
   y          current ARKode solution N_Vector for processing
   user_data  the structure passed by the user to the 
              ARKodeSetUserData routine.
---------------------------------------------------------------*/
typedef int (*ARKPostProcessStepFn)(realtype t, N_Vector y, 
                                    void *user_data);


/*===============================================================
                       USER-CALLABLE ROUTINES
===============================================================*/

/*---------------------------------------------------------------
 Function : ARKodeCreate
-----------------------------------------------------------------
 ARKodeCreate creates an internal memory block for a problem to
 be solved by ARKODE.

 If successful, ARKodeCreate returns a pointer to initialized
 problem memory. This pointer should be passed to ARKodeInit.
 If an initialization error occurs, ARKodeCreate prints an error
 message to standard err and returns NULL.
---------------------------------------------------------------*/
SUNDIALS_EXPORT void *ARKodeCreate();

/*---------------------------------------------------------------
 Integrator optional input specification functions
-----------------------------------------------------------------
 The following functions can be called to set optional inputs
 to values other than the defaults given below.  


 Function                 |  Optional input / [ default value ]
-----------------------------------------------------------------
 ARKodeSetDefaults        | resets all optional inputs to ARKode 
                          | default values.  Does not change 
                          | problem-defining function pointers 
                          | fe and fi or user_data pointer.  Also
                          | leaves alone any data structures or 
                          | options related to root-finding 
                          | (those can be reset using 
                          ! ARKodeRootInit).
                          | [internal]
                          |
 ARKodeSetDenseOrder      | polynomial order to be used for dense 
                          | output.  Allowed values are between 0
                          | and min(q,5) (where q is the order of
                          | the integrator)
                          | [3]
                          |
 ARKodeSetErrHandlerFn    | user-provided ErrHandler function.
                          | [internal]
                          |
 ARKodeSetErrFile         | the file pointer for an error file
                          | where all ARKODE warning and error
                          | messages will be written if the 
                          | default internal error handling 
                          | function is used. This parameter can 
                          | be stdout (standard output), stderr 
                          | (standard error), or a file pointer 
                          | (corresponding to a user error file 
                          | opened for writing) returned by fopen.
                          | If not called, then all messages will
                          | be written to stderr.
                          | [stderr]
                          |
 ARKodeSetUserData        | a pointer to user data that will be
                          | passed to the user's f function every
                          | time f is called.
                          | [NULL]
                          |
 ARKodeSetDiagnostics     | the file pointer for a diagnostics file 
                          | where all ARKODE step adaptivity and solver 
                          | information is written.  This parameter can 
                          | be stdout or stderr, though the preferred 
                          | approach is to specify a file pointer 
                          | (corresponding to a user diagnostics file 
                          | opened for writing) returned by fopen.  If 
                          | not called, or if called with a NULL file
                          | pointer, all diagnostics output is disabled.
                          | NOTE: when run in parallel, only one process 
                          | should set a non-NULL value for this pointer, 
                          | since statistics from all processes would be 
                          | identical.
                          | [NULL]
                          |
 ARKodeSetMaxNumSteps     | maximum number of internal steps to be
                          | taken by the solver in its attempt to
                          | reach tout.
                          | [500]
                          |
 ARKodeSetMaxHnilWarns    | maximum number of warning messages
                          | issued by the solver that t+h==t on 
                          | the next internal step. A value of -1 
                          | means no such messages are issued.
                          | [10]
                          |
 ARKodeSetInitStep        | initial step size.
                          | [estimated by ARKODE]
                          |
 ARKodeSetMinStep         | minimum absolute value of step size
                          | allowed.
                          | [0.0]
                          |
 ARKodeSetMaxStep         | maximum absolute value of step size
                          | allowed.
                          | [infinity]
                          |
 ARKodeSetStopTime        | the independent variable value past
                          | which the solution is not to proceed.
                          | [infinity]
                          |
 ARKodeSetFixedStep       | specifies to use a fixed step size 
                          | throughout integration
                          | [off]
-----------------------------------------------------------------
 ARKodeSetRootDirection      | Specifies the direction of zero
                             | crossings to be monitored
                             | [both directions]
                             |
 ARKodeSetNoInactiveRootWarn | disable warning about possible
                             | g==0 at beginning of integration
-----------------------------------------------------------------
 Return flag:
   ARK_SUCCESS   if successful
   ARK_MEM_NULL  if the arkode memory is NULL
   ARK_ILL_INPUT if an argument has an illegal value
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeSetDefaults(void *arkode_mem);
SUNDIALS_EXPORT int ARKodeSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int ARKodeSetErrHandlerFn(void *arkode_mem, 
                                          ARKErrHandlerFn ehfun, 
                                          void *eh_data);
SUNDIALS_EXPORT int ARKodeSetErrFile(void *arkode_mem, 
                                     FILE *errfp);
SUNDIALS_EXPORT int ARKodeSetUserData(void *arkode_mem, 
                                      void *user_data);
SUNDIALS_EXPORT int ARKodeSetDiagnostics(void *arkode_mem, 
                                         FILE *diagfp);
SUNDIALS_EXPORT int ARKodeSetMaxNumSteps(void *arkode_mem, 
                                         long int mxsteps);
SUNDIALS_EXPORT int ARKodeSetMaxHnilWarns(void *arkode_mem, 
                                          int mxhnil);
SUNDIALS_EXPORT int ARKodeSetInitStep(void *arkode_mem, 
                                      realtype hin);
SUNDIALS_EXPORT int ARKodeSetMinStep(void *arkode_mem, 
                                     realtype hmin);
SUNDIALS_EXPORT int ARKodeSetMaxStep(void *arkode_mem, 
                                     realtype hmax);
SUNDIALS_EXPORT int ARKodeSetStopTime(void *arkode_mem, 
                                      realtype tstop);
SUNDIALS_EXPORT int ARKodeSetFixedStep(void *arkode_mem, 
                                       realtype hfixed);

SUNDIALS_EXPORT int ARKodeSetRootDirection(void *arkode_mem, 
                                           int *rootdir);
SUNDIALS_EXPORT int ARKodeSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int ARKodeSetPostprocessStepFn(void *arkode_mem,
                                               ARKPostProcessStepFn ProcessStep);

/*---------------------------------------------------------------
 Function : ARKodeResize
-----------------------------------------------------------------
 ARKodeResize re-initializes ARKODE's memory for a problem with a
 changing vector size.  It is assumed that the problem dynamics 
 before and after the vector resize will be comparable, so that 
 all time-stepping heuristics prior to calling ARKodeResize 
 remain valid after the call.  If instead the dynamics should be 
 re-calibrated, the ARKode memory structure should be deleted 
 with a call to ARKodeFree, and re-created with calls to 
 ARKodeCreate and ARKodeInit.

 To aid in the vector-resize operation, the user can supply a 
 vector resize function, that will take as input an N_Vector with
 the previous size, and return as output a corresponding vector 
 of the new size.  If this function (of type ARKVecResizeFn) is 
 not supplied (i.e. is set to NULL), then all existing N_Vectors 
 will be destroyed and re-cloned from the input vector.

 In the case that the dynamical time scale should be modified 
 slightly from the previous time scale, an input "hscale" is 
 allowed, that will re-scale the upcoming time step by the 
 specified factor.  If a value <= 0 is specified, the default of 
 1.0 will be used.

 Other arguments:
   arkode_mem       Existing ARKode memory data structure.
   ynew             The newly-sized solution vector, holding 
                    the current dependent variable values.
   t0               The current value of the independent 
                    variable.
   resize_data      User-supplied data structure that will be 
                    passed to the supplied resize function.

 The return value of ARKodeResize is equal to ARK_SUCCESS = 0 if
 there were no errors; otherwise it is a negative int equal to:
   ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                    ARKodeCreate has not been called).
   ARK_NO_MALLOC    indicating that arkode_mem has not been
                    allocated (i.e., ARKodeInit has not been
                    called).
   ARK_ILL_INPUT    indicating an input argument was illegal
                    (including an error from the supplied 
                    resize function).
 In case of an error return, an error message is also printed.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeResize(void *arkode_mem, N_Vector ynew, 
                                 realtype hscale, realtype t0,
                                 ARKVecResizeFn resize, 
                                 void *resize_data);

/*---------------------------------------------------------------
 Functions : ARKodeSStolerances
             ARKodeSVtolerances
             ARKodeWFtolerances
-----------------------------------------------------------------
 These functions specify the integration tolerances. One of them
 SHOULD be called before the first call to ARKode; otherwise 
 default values of reltol=1e-4 and abstol=1e-9 will be used, 
 which may be entirely incorrect for a specific problem.

 ARKodeSStolerances specifies scalar relative and absolute 
   tolerances.
 ARKodeSVtolerances specifies scalar relative tolerance and a 
   vector absolute tolerance (a potentially different absolute 
   tolerance for each vector component).
 ARKodeWFtolerances specifies a user-provided function (of type 
   ARKEwtFn) which will be called to set the error weight vector.

 The tolerances reltol and abstol define a vector of error 
 weights, ewt, with components
   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in SS case), or
   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in SV case).
 This vector is used in all error and convergence tests, which
 use a weighted RMS norm on all error-like vectors v:
    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
 where N is the problem dimension.

 The return value of these functions is equal to ARK_SUCCESS=0 
 if there were no errors; otherwise it is a negative int equal 
 to:
   ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                    ARKodeCreate has not been called).
   ARK_NO_MALLOC    indicating that arkode_mem has not been
                    allocated (i.e., ARKodeInit has not been
                    called).
   ARK_ILL_INPUT    indicating an input argument was illegal
                    (e.g. a negative tolerance)
 In case of an error return, an error message is also printed.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeSStolerances(void *arkode_mem, 
                                       realtype reltol, 
                                       realtype abstol);
SUNDIALS_EXPORT int ARKodeSVtolerances(void *arkode_mem, 
                                       realtype reltol, 
                                       N_Vector abstol);
SUNDIALS_EXPORT int ARKodeWFtolerances(void *arkode_mem, 
                                       ARKEwtFn efun);

/*---------------------------------------------------------------
 Functions : ARKodeResStolerance
             ARKodeResVtolerance
             ARKodeResFtolerance
-----------------------------------------------------------------
 These functions specify the absolute residual tolerance. 
 Specification of the absolute residual tolerance is only 
 necessary for problems with non-identity mass matrices in which
 the units of the solution vector y dramatically differ from the 
 units of My, where M is the user-supplied mass matrix.  If this
 occurs, one of these routines SHOULD be called before the first 
 call to ARKode; otherwise the default value of rabstol=1e-9 will 
 be used, which may be entirely incorrect for a specific problem.

 ARKodeResStolerances specifies a scalar residual tolerance.

 ARKodeResVtolerances specifies a vector residual tolerance 
  (a potentially different absolute residual tolerance for 
   each vector component).

 ARKodeResFtolerances specifies a user-provides function (of 
   type ARKRwtFn) which will be called to set the residual 
   weight vector.

 The tolerances reltol (defined for both the solution and 
 residual) and rabstol define a vector of residual weights, 
 rwt, with components
   rwt[i] = 1/(reltol*abs(My[i]) + abstol)     (in S case), or
   rwt[i] = 1/(reltol*abs(My[i]) + abstol[i])  (in V case).
 This vector is used in all solver convergence tests, which
 use a weighted RMS norm on all residual-like vectors v:
    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*rwt[i])^2 ),
 where N is the problem dimension.

 The return value of these functions is equal to ARK_SUCCESS=0 
 if there were no errors; otherwise it is a negative int equal 
 to:
   ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                    ARKodeCreate has not been called).
   ARK_NO_MALLOC    indicating that arkode_mem has not been
                    allocated (i.e., ARKodeInit has not been
                    called).
   ARK_ILL_INPUT    indicating an input argument was illegal
                    (e.g. a negative tolerance)
 In case of an error return, an error message is also printed.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeResStolerance(void *arkode_mem, 
                                        realtype rabstol);
SUNDIALS_EXPORT int ARKodeResVtolerance(void *arkode_mem, 
                                        N_Vector rabstol);
SUNDIALS_EXPORT int ARKodeResFtolerance(void *arkode_mem, 
                                        ARKRwtFn rfun);

/*---------------------------------------------------------------
 Function : ARKodeRootInit
-----------------------------------------------------------------
 ARKodeRootInit initializes a rootfinding problem to be solved
 during the integration of the ODE system.  It must be called
 after ARKodeCreate, and before ARKode.  The arguments are:

 arkode_mem = pointer to ARKODE memory returned by ARKodeCreate.

 nrtfn      = number of functions g_i, an integer >= 0.

 g          = name of user-supplied function, of type ARKRootFn,
              defining the functions g_i whose roots are sought.

 If a new problem is to be solved with a call to either
 ARKStepReInit or ERKStepReInit, where the new problem has no
 root functions but the prior one did, then call ARKodeRootInit
 again with nrtfn = 0.

 The return value of ARKodeRootInit is ARK_SUCCESS = 0 if there 
 were no errors; otherwise it is a negative int equal to:
   ARK_MEM_NULL    indicating arkode_mem was NULL, or
   ARK_MEM_FAIL    indicating a memory allocation failed.
                   (including an attempt to increase maxord).
   ARK_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
 In case of an error return, an error message is also printed.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeRootInit(void *arkode_mem, 
                                   int nrtfn, 
                                   ARKRootFn g);

/*---------------------------------------------------------------
 Function : ARKode
-----------------------------------------------------------------
 ARKode integrates the ODE over an interval in t.

 ARKode may be run in one of two modes (ARK_NORMAL or 
 ARK_ONE_STEP), as determined by the itask argument:

 If itask is ARK_NORMAL, then the solver integrates from its
 current internal t value to a point at or beyond tout, then
 interpolates to t = tout and returns y(tout) in the user-
 allocated vector yout.  This interpolation is typically less 
 accurate than the full time step solutions produced by the 
 solver, since the interpolating polynomial relies on the 
 internal stage solutions, that may have reduced accuracy in 
 comparison with the full time step solutions.  If the user 
 wishes that this returned value have full method accuracy, they 
 may issue a call to ARKodeSetStopTime before the call to ARKode 
 to specify a fixed stop time to end the time step and return to 
 the user.  Once the integrator returns at a tstop time, any 
 future testing for tstop is disabled (and can be reenabled only 
 though a new call to ARKodeSetStopTime).

 If itask is ARK_ONE_STEP, then the solver takes one internal 
 time step and returns in yout the value of y at the new internal 
 time. In this case, tout is used only during the first call to 
 ARKode to determine the direction of integration and the rough 
 scale of the t variable.  As with the ARK_NORMAL mode, a user 
 may specify a specific stop time for output of this step, 
 assuming that the requested step is smaller than the step taken 
 by the method. 

 The time reached by the solver is placed in (*tret). The
 user is responsible for allocating the memory for this value.

 arkode_mem is the pointer to ARKODE memory returned by
            ARKodeCreate.

 tout  is the next time at which a computed solution is desired.

 yout  is the computed solution vector. In ARK_NORMAL mode with no
       errors and no roots found, yout=y(tout).

 tret  is a pointer to a real location. ARKode sets (*tret) to
       the time reached by the solver and returns
       yout=y(*tret).

 itask is ARK_NORMAL or ARK_ONE_STEP, as described above.

 Here is a brief description of each return value:

 ARK_SUCCESS:      ARKode succeeded and no roots were found.

 ARK_ROOT_RETURN:  ARKode succeeded, and found one or more roots.
                   If nrtfn > 1, call ARKodeGetRootInfo to see
                   which g_i were found to have a root at (*tret).

 ARK_TSTOP_RETURN: ARKode succeeded and returned at tstop.

 ARK_MEM_NULL:     The arkode_mem argument was NULL.

 ARK_NO_MALLOC:    arkode_mem was not allocated.

 ARK_ILL_INPUT:    One of the inputs to ARKode is illegal. This
                   includes the situation when a component of the
                   error weight vectors becomes < 0 during
                   internal time-stepping.  It also includes the
                   situation where a root of one of the root
                   functions was found both at t0 and very near t0.
                   The ILL_INPUT flag will also be returned if the
                   linear solver routine ARK--- (called by the user
                   after calling ARKodeCreate) failed to set one of
                   the linear solver-related fields in arkode_mem or
                   if the linear solver's init routine failed. In
                   any case, the user should see the printed
                   error message for more details.

 ARK_TOO_MUCH_WORK: The solver took mxstep internal steps but
                   could not reach tout prior to the maximum number
                   of steps (ark_mxstep).

 ARK_TOO_MUCH_ACC: The solver could not satisfy the accuracy
                   demanded by the user for some internal step.

 ARK_ERR_FAILURE:  Error test failures occurred too many times
                   (= ark_maxnef) during one internal time step 
                   or occurred with |h| = hmin.

 ARK_CONV_FAILURE: Convergence test failures occurred too many
                   times (= ark_maxncf) during one internal time
                   step or occurred with |h| = hmin.

 ARK_LINIT_FAIL:   The linear solver's initialization function 
                   failed.

 ARK_LSETUP_FAIL:  The linear solver's setup routine failed in an
                   unrecoverable manner.

 ARK_LSOLVE_FAIL:  The linear solver's solve routine failed in an
                   unrecoverable manner.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKode(void *arkode_mem, realtype tout, 
                           N_Vector yout, realtype *tret, int itask);

/*---------------------------------------------------------------
 Function : ARKodeGetDky
-----------------------------------------------------------------
 ARKodeGetDky computes the kth derivative of the y function at
 time t, where tn-hu <= t <= tn, tn denotes the current
 internal time reached, and hu is the last internal step size
 successfully used by the solver. The user may request
 k=0, 1, ..., d, where d = min(5,q), with q the order of accuracy 
 for the time integration method. The derivative vector is 
 returned in dky. This vector must be allocated by the caller. It 
 is only legal to call this function after a successful return 
 from ARKode.

 arkode_mem is the pointer to ARKODE memory returned by
            ARKodeCreate.

 t   is the time at which the kth derivative of y is evaluated.
     The legal range for t is [tn-hu,tn] as described above.

 k   is the order of the derivative of y to be computed. The
     legal range for k is [0,min(q,3)] as described above.

 dky is the output derivative vector [((d/dy)^k)y](t).

 The return value for ARKodeGetDky is one of:

   ARK_SUCCESS:  ARKodeGetDky succeeded.

   ARK_BAD_K:    k is not in the range 0, 1, ..., s-1.

   ARK_BAD_T:    t is not in the interval [tn-hu,tn].

   ARK_BAD_DKY:  The dky argument was NULL.

   ARK_MEM_NULL: The arkode_mem argument was NULL.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeGetDky(void *arkode_mem, realtype t, 
                                 int k, N_Vector dky);

/*---------------------------------------------------------------
 Integrator optional output extraction functions
-----------------------------------------------------------------
 The following functions can be called to get optional outputs
 and statistics related to the main integrator.
-----------------------------------------------------------------
 ARKodeGetNumSteps returns the cumulative number of internal
                   steps taken by the solver

 ARKodeGetActualInitStep returns the actual initial step size
                         used by ARKODE

 ARKodeGetLastStep returns the step size for the last internal
                   step

 ARKodeGetCurrentStep returns the step size to be attempted on
                      the next internal step

 ARKodeGetCurrentTime returns the current internal time reached
                      by the solver

 ARKodeGetTolScaleFactor returns a suggested factor by which the
                         user's tolerances should be scaled when
                         too much accuracy has been requested for
                         some internal step

 ARKodeGetErrWeights returns the current error weight vector.
                     The user must allocate space for eweight.

 ARKodeGetResWeights returns the current residual weight vector.
                     The user must allocate space for rweight.

 ARKodeGetWorkSpace returns the ARKODE real and integer workspaces

 ARKodeGetNumGEvals returns the number of calls to the user's
                    g function (for rootfinding)

 ARKodeGetRootInfo returns the indices for which g_i was found to 
                   have a root. The user must allocate space for 
                   rootsfound. For i = 0 ... nrtfn-1, 
                   rootsfound[i] = 1 if g_i has a root, and = 0 
                   if not.

 ARKodeGet* return values:
   ARK_SUCCESS   if succesful
   ARK_MEM_NULL  if the arkode memory was NULL
---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKodeGetWorkSpace(void *arkode_mem, 
                                       long int *lenrw, 
                                       long int *leniw);
SUNDIALS_EXPORT int ARKodeGetNumSteps(void *arkode_mem, 
                                      long int *nsteps);
SUNDIALS_EXPORT int ARKodeGetActualInitStep(void *arkode_mem, 
                                            realtype *hinused);
SUNDIALS_EXPORT int ARKodeGetLastStep(void *arkode_mem, 
                                      realtype *hlast);
SUNDIALS_EXPORT int ARKodeGetCurrentStep(void *arkode_mem, 
                                         realtype *hcur);
SUNDIALS_EXPORT int ARKodeGetCurrentTime(void *arkode_mem, 
                                         realtype *tcur);
SUNDIALS_EXPORT int ARKodeGetTolScaleFactor(void *arkode_mem, 
                                            realtype *tolsfac);
SUNDIALS_EXPORT int ARKodeGetErrWeights(void *arkode_mem, 
                                        N_Vector eweight);
SUNDIALS_EXPORT int ARKodeGetResWeights(void *arkode_mem, 
                                        N_Vector rweight);
SUNDIALS_EXPORT int ARKodeGetNumGEvals(void *arkode_mem, 
                                       long int *ngevals);
SUNDIALS_EXPORT int ARKodeGetRootInfo(void *arkode_mem, 
                                      int *rootsfound);

/*---------------------------------------------------------------
 As a convenience, the following function provides many of the
 optional outputs as a group.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeGetStepStats(void *arkode_mem, 
                                       long int *nsteps,
                                       realtype *hinused, 
                                       realtype *hlast, 
                                       realtype *hcur, 
                                       realtype *tcur);

/*---------------------------------------------------------------
 The following function returns the name of the constant 
 associated with a ARKODE return flag
---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ARKodeGetReturnFlagName(long int flag);

/*---------------------------------------------------------------
 Function : ARKodeWriteParameters
-----------------------------------------------------------------
 ARKodeWriteParameters outputs all solver parameters to the 
 provided file pointer.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKodeWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
 Function : ARKodeFree
-----------------------------------------------------------------
 ARKodeFree frees the problem memory arkode_mem allocated by
 ARKodeCreate and ARKodeInit. Its only argument is the pointer
 arkode_mem returned by ARKodeCreate.
---------------------------------------------------------------*/
SUNDIALS_EXPORT void ARKodeFree(void **arkode_mem);

#ifdef __cplusplus
}
#endif

#endif
