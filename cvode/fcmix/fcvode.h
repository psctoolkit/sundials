/******************************************************************
 * File          : fcvode.h                                       *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 30 MArch 2003                                  *
 *----------------------------------------------------------------*
 * This is the header file for FCVODE, the Fortran interface to   *
 * the CVODE package.                                             *
 ******************************************************************/

#ifndef _fcvode_h
#define _fcvode_h


/******************************************************************************

                  FCVODE Interface Package

The FCVODE Interface Package is a package of C functions which support
the use of the CVODE solver, for the solution of ODE systems 
dy/dt = f(t,y), in a mixed Fortran/C setting.  While CVODE is written
in C, it is assumed here that the user's calling program and
user-supplied problem-defining routines are written in Fortran.  
This package provides the necessary interface to CVODE for both the
serial and the parallel NVECTOR implementations.

The user-callable functions, with the corresponding CVODE functions,
are as follows:

  FMENVINITS and FMENVINITP interface to M_EnvInit_Serial and
             M_EnvInitParallel (defined by nvector_serial
             and nvector_parallel, respectively)

  FCVMALLOC  interfaces to CVodeMalloc

  FCVREINIT  interfaces to CVReInit

  FCVDIAG    interfaces to CVDiag

  FCVDENSE0, FCVDENSE1, FCVREINDENSE0, FCVREINDENSE1
             interface to CVDense or CVReInitDense for the various options

  FCVBAND0, FCVBAND1, FCVREINBAND0, FCVREINBAND1
             interface to CVBand or CVReInitBand for the various options

  FCVSPGMR00, FCVSPGMR01, FCVSPGMR10, FCVSPGMR11, FCVSPGMR20, FCVSPGMR21,
  FCVREINSPGMR00, FCVREINSPGMR01, FCVREINSPGMR10, FCVREINSPGMR11,
  FCVREINSPGMR20, FCVREINSPGMR21,
             interface to CVSpgmr or CVReInitSpgmr for the various options

  FCVODE     interfaces to CVode

  FCVDKY     interfaces to CVodeDky

  FCVFREE    interfaces to CVodeFree

  FMENVFREES and FMENVFREEP interface to M_EnvFree_Serial and
             M_EnvFree_parallel(defined by nvector_serial
             and nvector_parallel, respectively)

The user-supplied functions, each listed with the corresponding interface
function which calls it (and its type within CVODE), are as follows:
  CVFUN    is called by the interface function CVf of type RhsFn
  CVDJAC   is called by the interface function CVDenseJac of type CVDenseJacFn
  CVBJAC   is called by the interface function CVBandJac of type CVBandJacFn
  CVPSOL   is called by the interface function CVPSol of type CVSpgmrSolveFn
  CVPRECO  is called by the interface function CVPreco of type CVSpgmrPrecondFn
  CVJTIMES is called by the interface function CVJtimes of type CVSpgmrJtimesFn
In contrast to the case of direct use of CVODE, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program.

Important note on portability.
In this package, the names of the interface functions, and the names of
the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header file fcvode.h.  Those mapping definitions depend in turn on a
pair of parameters, CRAY and UNDERSCORE, defined in the header file
fcmixpar.h, which is machine-dependent.  The names into which the dummy
names are mapped are either upper or lower case, and may or may not have
an underscore appended, depending on these parameters.  Check, and if 
necessary modify, the file fcmixpar.h for a given machine environment.

===============================================================================

                 Usage of the FPVODE Interface Package

The usage of FCVODE requires calls to six or seven interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These
function calls and user routines are summarized separately below.

Some details are omitted, and the user is referred to the user documents
on CVODE and PVODE for more complete documentation.  Information on the
arguments of any given user-callable interface routine, or of a given
user-supplied function called by an interface function, can be found in
the documentation on the corresponding function in the CVODE package.

The number labels on the instructions below end with s for instructions
that apply to the serial version of CVODE only, and end with p for
those that apply to the parallel version only.


(1) User-supplied right-hand side routine: CVFUN
The user must in all cases supply the following Fortran routine
      SUBROUTINE CVFUN (T, Y, YDOT)
      DIMENSION Y(*), YDOT(*)
It must set the YDOT array to f(t,y), the right-hand side of the ODE
system, as function of T = t and the array Y = y.  Here Y and YDOT
are distributed vectors.

(2) Optional user-supplied dense Jacobian approximation routine: CVDJAC
As an option when using the DENSE linear solver, the user may supply a
routine that computes a dense approximation of the system Jacobian 
J = df/dy. If supplied, it must have the following form:
      SUBROUTINE CVDJAC (NEQ, T, Y, FY, DJAC, WK1, WK2, WK3)
      DIMENSION Y(*), FY(*), EWT(*), DJAC(NEQ,*), WK1(*), WK2(*), WK3(*)
Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
the Jacobian and store it columnwise in DJAC.

(3) Optional user-supplied band Jacobian approximation routine: CVBJAC
As an option when using the BAND linear solver, the user may supply a
routine that computes a band approximation of the system Jacobian 
J = df/dy. If supplied, it must have the following form:
      SUBROUTINE CVBJAC (NEQ, MU, ML, MDIM, T, Y, FY,
     1                   BJAC, WK1, WK2, WK3)
      DIMENSION Y(*), FY(*), EWT(*), BJAC(MDIM,*), WK1(*), WK2(*), WK3(*)
Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. 
It must load the MDIM by N array BJAC with the Jacobian matrix at the
current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.

(4) Optional user-supplied Jacobian-vector product routine: CVJTIMES
As an option when using the SPGMR linear solver, the user may supply a 
routine that computes the product of the system Jacobian J = df/dy and 
a given vector v.  If supplied, it must have the following form:
      SUBROUTINE CVJTIMES (V, FJV, T, Y, FY, WORK, IER)
      DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), WORK(*)
Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
compute the product vector Jv, where the vector v is stored in V, and store
the product in FJV.  On return, set IER = 0 if CVJTIMES was successful,
and nonzero otherwise.

(5) Initialization:  FMENVINITS / FMENVINITP , FCVMALLOC, FCVREINIT

(5.1s) To initialize the serial machine environment, the user must make
the following call:
       CALL FMENVINITS (NEQ, IER)
The arguments are:
NEQ     = size of vectors
IER     = return completion flag. Values are 0 = success, -1 = failure.

(5.1p) To initialize the parallel machine environment, the user must make 
the following call:
       CALL FMENVINITP (NLOCAL, NGLOBAL, IER)
The arguments are:
NLOCAL  = local size of vectors on this processor
NGLOBAL = the system size, and the global size of vectors (the sum 
          of all values of NLOCAL)
IER     = return completion flag. Values are 0 = success, -1 = failure.
Note: If MPI was initialized by the user, the communicator must be
set to MPI_COMM_WORLD.  If not, this routine initializes MPI and sets
the communicator equal to MPI_COMM_WORLD.

(5.2) To set various problem and solution parameters and allocate
internal memory, make the following call:
      CALL FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1               IOPT, ROPT, IER)
The arguments are:
T0     = initial value of t
Y0     = array of initial conditions
METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
ITMETH = nonlinear iteration method: 1 = functional iteration, 2 = Newton iter.
IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array
RTOL   = relative tolerance (scalar)
ATOL   = absolute tolerance (scalar or array)
INOPT  = optional input flag: 0 = none, 1 = inputs used
IOPT   = array of length 40 for integer optional inputs and outputs
         (declare as INTEGER*4 or INTEGER*8 according to C type long int)
ROPT   = array of length 40 for real optional inputs and outputs
         The optional inputs are MAXORD, MXSTEP, MXHNIL, SLDET, H0, HMAX,
         HMIN, stored in IOPT(1), IOPT(2), IOPT(3), IOPT(14), ROPT(1),
         ROPT(2), ROPT(3), respectively.  If any of these optional inputs
         are used, set the others to zero to indicate default values.
         The optional outputs are NST, NFE, NSETUPS, NNI, NCFN, NETF, QU, QCUR,
         LENRW, LENIW, NOR, HU,HCUR, TCUR, TOLSF, stored in IOPT(4) .. IOPT(13),
         IOPT(15), ROPT(4) .. ROPT(7), resp.  See the CVODE manual for details. 
IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
         See printed message for details in case of failure.

(5.3) To re-initialize the CVODE solver for the solution of a new problem
of the same size as one already solved, make the following call:
      CALL FCVREINIT(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
     1               IOPT, ROPT, IER)
The arguments have the same names and meanings as those of FCVMALLOC,
except that NEQ has been omitted from the argument list (being unchanged
for the new problem).  FCVREINIT performs the same initializations as
FCVMALLOC, but does no memory allocation, using instead the existing
internal memory created by the previous FCVMALLOC call.  The call to
specify the linear system solution method may or may not be needed;
see paragraph (6) below.

(6) Specification of linear system solution method.
In the case of a stiff system, the implicit BDF method involves the solution
of linear systems related to the Jacobian J = df/dy of the ODE system.
CVODE presently includes four choices for the treatment of these systems,
and the user of FCVODE must call a routine with a specific name to make the
desired choice.

(6.1) Diagonal approximate Jacobian.
This choice is appropriate when the Jacobian can be well approximated by
a diagonal matrix.  The user must make the call:
      CALL FCVDIAG (IER)
IER is an error return flag: 0 = success, -1 = memory failure.
There is no additional user-supplied routine.  Optional outputs specific
to the approximate diagonal Jacobian case are LRW and LIW, stored in
IOPT(16) and IOPT(17), respectively.  (See the CVODE manual for descriptions.)

(6.2) DENSE treatment of the linear system.
The user must make one of the following two calls:
      CALL FCVDENSE0(NEQ, IER)
          if CVDJAC is not supplied 
or
      CALL FCVDENSE1(NEQ, IER)
          if CVDJAC is supplied 

In both cases, the argument is:
IER = error return flag: 0 = success , -1 = memory allocation failure,
                        -2 = illegal input. 

  In the case FCVDENSE1, the user program must include the CVDJAC routine 
for the evaluation of the dense approximation to the Jacobian.

     Optional outputs specific to the DENSE case are NJE, LRW, and LIW,
stored in IOPT(16), IOPT(17), and IOPT(18), respectively.  (See the CVODE
manual for descriptions.)

     If a sequence of problems of the same size is being solved using the
DENSE linear solver, then following the call to FCVREINIT, a call to specify
the linear solver may or may not be needed.  If the choice of user vs internal
Jacobian routine is unchanged, no such call is needed, and the user program can
proceed to the FCVODE calls.  However, if there is a change in this choice, then
a call to FCVREINDENSE0(IER) or FCVREINDENSE1(IER) must be made, for the case
of a CVDJAC routine not supplied or supplied, respectively.  This call
reinitializes the DENSE linear solver without reallocating its memory.
These routines have a return argument IER, with the same meaning as for the
FCVDENSE* routines.

(6.3) BAND treatment of the linear system
The user must make one of the following two calls:
      CALL FCVBAND0(NEQ, MU, ML, IER)
          if CVBJAC is not supplied
or
      CALL FCVBAND1(NEQ, MU, ML, IER)
          if CVBJAC is supplied

In both cases, the arguments are:
MU  = upper bandwidth
ML  = lower bandwidth
IER = error return flag: 0 = success , -1 = memory allocation failure,
                        -2 = illegal input.     

  In the case FCVBAND1, the user program must include the CVBJAC routine 
for the evaluation of the banded approximation to the Jacobian.

     Optional outputs specific to the BAND case are NJE, LRW, and LIW,
stored in IOPT(16), IOPT(17), and IOPT(18), respectively.  (See the CVODE
manual for descriptions.)

     If a sequence of problems of the same size is being solved using the
BAND linear solver, then following the call to FCVREINIT, a call to specify
the linear solver may or may not be needed.  If the values of MU and ML and
the choice of user vs internal Jacobian routine are unchanged, no such call is
needed, and the user program can proceed to the FCVODE calls.  If MU and ML
are unchanged, but the choice of Jacobian routine is changed, then a call to
FCVREINBAND0 or FCVREINBAND1 must be made, for the case of a CVBJAC routine
not supplied or supplied, respectively.  This call reinitializes the BAND
linear solver without reallocating its memory.  These routines have the same
arguments (with the same meanings) as for the FCVBAND* routines.  If there
is a change to MU or ML, then one of the FCVBAND* routines must be called.

(6.4) SPGMR treatment of the linear systems.
For the Scaled Preconditioned GMRES solution of the linear systems,
the user must make one of the following six calls:
      CALL FCVSPGMR00 (IGSTYPE, MAXL, DELT, IER)              
          if no preconditioning is to be done and CVJTIMES is not supplied;

      CALL FCVSPGMR10 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
          if the preconditioner involves no data setup and CVJTIMES is not
          supplied;

      CALL FCVSPGMR20 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
          if the preconditioner involves data setup and CVJTIMES is not
          supplied;

      CALL FCVSPGMR01 (IGSTYPE, MAXL, DELT, IER)
          if the preconditioning is to be done but CVJTIMES is supplied;

      CALL FCVSPGMR11 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
          if the preconditioner involves no data setup but CVJTIMES is supplied;

      CALL FCVSPGMR21 (IPRETYPE, IGSTYPE, MAXL, DELT, IER)
          if the preconditioner involves data setup and CVJTIMES is supplied.

(In the two-digit suffix on the name above, the first digit is the number of
preconditioner routines supplied, and the second digit is 1 or 0 according as
CVJTIMES is supplied or not.)

In all cases, the arguments are:
IPRETYPE = preconditioner type: 1 = left only, 2 = right only, 3 = both sides.
IGSTYPE  = Gram-schmidt process type: 0 = modified G-S, 1 = classical G-S.
MAXL     = maximum Krylov subspace dimension; 0 indicates default.
DELT     = linear convergence tolerance factor; 0.0 indicates default.
IER      = error return flag: 0 = success, -1 = memory allocation failure,
           -2 = illegal input.

     In the cases FCVSPGMR10, FCVSPGMR11, FCVSPGMR20, and FCVSPGMR21, the user
program must include the following routine for solution of the preconditioner
linear system:
      SUBROUTINE CVPSOL (T, Y,FY, VT, GAMMA, EWT,DELTA, NFE, R, LR, Z, IER)
      DIMENSION Y(*), FY(*), VT(*), EWT(*), R(*), Z(*),
Typically this routine will use only NEQ, T, Y, GAMMA, R, LR, and Z.  It
must solve the preconditioner linear system Pz = r, where r = R is input, 
and store the solution z in Z.  Here P is the left preconditioner if LR = 1
and the right preconditioner if LR = 2.  The preconditioner (or the product
of the left and right preconditioners if both are nontrivial) should be an 
approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).

     In the cases FCVSPGMR20 and FCVSPGMR21, the user program must also include
the following routine for the evaluation and preprocessing of the preconditioner:
      SUBROUTINE CVPRECO (T, Y, FY, JOK, JCUR, GAMMA, EWT, H, UROUND, 
     1                   NFE, V1, V2, V3, IER)
      DIMENSION Y(*), FY(*), EWT(*), V1(*), V2(*), V3(*) 
Typically this routine will use only NEQ, T, Y, JOK, and GAMMA. It must
perform any evaluation of Jacobian-related data and preprocessing needed
for the solution of the preconditioner linear systems by CVPSOL.
The JOK argument allows for Jacobian data to be saved and reused:  If 
JOK = 0, this data should be recomputed from scratch.  If JOK = 1, a saved
copy of it may be reused, and the preconditioner constructed from it.
On return, set JCUR = 1 if Jacobian data was computed, and 0 otherwise.
Also on return, set IER = 0 if CVPRECO was successful, set IER positive if a 
recoverable error occurred, and set IER negative if a non-recoverable error
occurred.

     Optional outputs specific to the SPGMR case are NPE, NLI, NPS, NCFL,
LRW, and LIW, stored in IOPT(16) ... IOPT(21), respectively.  (See the CVODE
manual for descriptions.)

     If a sequence of problems of the same size is being solved using the SPGMR
linear solver, then following the call to FCVREINIT, a call to the FCVSPGMR**
routine may or may not be needed.  First, if the choice among the six SPGMR
options is the same and the input arguments are the same, no FCVSPGMR** call
is needed.  If a different choice of options is desired, or there is a change
in input arguments other than MAXL, then the user program should call one of
the routines FCVREINSPGMR00, FCVREINSPGMR01, FCVREINSPGMR10, FCVREINSPGMR11,
FCVREINSPGMR20, or FCVREINSPGMR21.  In this case, the FCVREINSPGMR** routine
reinitializes the SPGMR linear solver, but without reallocating its memory.
The arguments of each FCVREINSPGMR** routine have the same names and meanings
as the corresponding FCVSPGMR** routine.  Finally, if the value of MAXL is
being changed, then a call to one of the six FCVSPGMR** routines must be made,
where again a different choice of that routine is allowed.  

(7) The integrator: FCVODE
Carrying out the integration is accomplished by making calls as follows:
      CALL FCVODE (TOUT, T, Y, ITASK, IER)
The arguments are:
TOUT  = next value of t at which a solution is desired (input)
T     = value of t reached by the solver on output
Y     = array containing the computed solution on output
ITASK = task indicator: 0 = normal mode (overshoot TOUT and interpolate)
        1 = one-step mode (return after each internal step taken)
IER   = completion flag: 0 = success, values -1 ... -8 are various
        failure modes (see CVODE manual).
The current values of the optional outputs are available in IOPT and ROPT.

(8) Computing solution derivatives: FCVDKY
To obtain a derivative of the solution, of order up to the current method
order, make the following call:
      CALL FCVDKY (T, K, DKY, IER
The arguments are:
T   = value of t at which solution derivative is desired
K   = derivative order (0 .le. K .le. QU)
DKY = array containing computed K-th derivative of y on return
IER = return flag: = 0 for success, < 0 for illegal argument.

(9) Memory freeing: FCVFREE and FMENVFREES / FMENVFREEP
To the free the internal memory created by the calls to FCVMALLOC and
FMENVINITS or FMENVINITP, depending on the version (serial/parallel), make
the following calls, in this order:
      CALL FCVFREE
      CALL FMENVFREES or CALL FMENVFREEP  


******************************************************************************/

#include "fcmixpar.h"   /* parameters for function name definitions */

/* Definitions of interface function names */

#if (CRAY)

#define FCV_MALLOC  FCVMALLOC
#define FCV_REINIT  FCVREINIT
#define FCV_DIAG    FCVDIAG
#define FCV_DENSE0  FCVDENSE0
#define FCV_DENSE1  FCVDENSE1
#define FCV_REINDENSE0 FCVREINDENSE0
#define FCV_REINDENSE1 FCVREINDENSE1
#define FCV_BAND0   FCVBAND0
#define FCV_BAND1   FCVBAND1
#define FCV_REINBAND0 FCVREINBAND0
#define FCV_REINBAND1 FCVREINBAND1
#define FCV_SPGMR00 FCVSPGMR00
#define FCV_SPGMR10 FCVSPGMR10
#define FCV_SPGMR20 FCVSPGMR20
#define FCV_SPGMR01 FCVSPGMR01
#define FCV_SPGMR11 FCVSPGMR11
#define FCV_SPGMR21 FCVSPGMR21
#define FCV_REINSPGMR00 FCVREINSPGMR00
#define FCV_REINSPGMR10 FCVREINSPGMR10
#define FCV_REINSPGMR20 FCVREINSPGMR20
#define FCV_REINSPGMR01 FCVREINSPGMR01
#define FCV_REINSPGMR11 FCVREINSPGMR11
#define FCV_REINSPGMR21 FCVREINSPGMR21
#define FCV_CVODE   FCVODE
#define FCV_DKY     FCVDKY
#define FCV_FREE    FCVFREE
#define FCV_FUN     CVFUN
#define FCV_DJAC    CVDJAC
#define FCV_BJAC    CVBJAC
#define FCV_PSOL    CVPSOL
#define FCV_PRECO   CVPRECO
#define FCV_JTIMES  CVJTIMES

#elif (UNDERSCORE)

#define FCV_MALLOC  fcvmalloc_
#define FCV_REINIT  fcvreinit_
#define FCV_DIAG    fcvdiag_
#define FCV_DENSE0  fcvdense0_
#define FCV_DENSE1  fcvdense1_
#define FCV_REINDENSE0 fcvreindense0_
#define FCV_REINDENSE1 fcvreindense1_
#define FCV_BAND0   fcvband0_
#define FCV_BAND1   fcvband1_
#define FCV_REINBAND0 fcvreinband0_
#define FCV_REINBAND1 fcvreinband1_
#define FCV_SPGMR00 fcvspgmr00_
#define FCV_SPGMR10 fcvspgmr10_
#define FCV_SPGMR20 fcvspgmr20_
#define FCV_SPGMR01 fcvspgmr01_
#define FCV_SPGMR11 fcvspgmr11_
#define FCV_SPGMR21 fcvspgmr21_
#define FCV_REINSPGMR00 fcvreinspgmr00_
#define FCV_REINSPGMR10 fcvreinspgmr10_
#define FCV_REINSPGMR20 fcvreinspgmr20_
#define FCV_REINSPGMR01 fcvreinspgmr01_
#define FCV_REINSPGMR11 fcvreinspgmr11_
#define FCV_REINSPGMR21 fcvreinspgmr21_
#define FCV_CVODE   fcvode_
#define FCV_DKY     fcvdky_
#define FCV_FREE    fcvfree_
#define FCV_FUN     cvfun_
#define FCV_DJAC    cvdjac_
#define FCV_BJAC    cvbjac_
#define FCV_PSOL    cvpsol_
#define FCV_PRECO   cvpreco_
#define FCV_JTIMES  cvjtimes_

#else

#define FCV_MALLOC  fcvmalloc
#define FCV_REINIT  fcvreinit
#define FCV_DIAG    fcvdiag
#define FCV_DENSE0  fcvdense0
#define FCV_DENSE1  fcvdense1
#define FCV_REINDENSE0 fcvreindense0
#define FCV_REINDENSE1 fcvreindense1
#define FCV_BAND0   fcvband0
#define FCV_BAND1   fcvband1
#define FCV_REINBAND0 fcvreinband0
#define FCV_REINBAND1 fcvreinband1
#define FCV_SPGMR00 fcvspgmr00
#define FCV_SPGMR10 fcvspgmr10
#define FCV_SPGMR20 fcvspgmr20
#define FCV_SPGMR01 fcvspgmr01
#define FCV_SPGMR11 fcvspgmr11
#define FCV_SPGMR21 fcvspgmr21
#define FCV_REINSPGMR00 fcvreinspgmr00
#define FCV_REINSPGMR10 fcvreinspgmr10 
#define FCV_REINSPGMR20 fcvreinspgmr20 
#define FCV_REINSPGMR01 fcvreinspgmr01 
#define FCV_REINSPGMR11 fcvreinspgmr11 
#define FCV_REINSPGMR21 fcvreinspgmr21 
#define FCV_CVODE   fcvode
#define FCV_DKY     fcvdky
#define FCV_FREE    fcvfree
#define FCV_FUN     cvfun
#define FCV_DJAC    cvdjac
#define FCV_BJAC    cvbjac
#define FCV_PSOL    cvpsol
#define FCV_PRECO   cvpreco
#define FCV_JTIMES  cvjtimes

#endif


/* CVODE header files  */

#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "cvode.h"         /* definition of type RHSFn                      */
#include "nvector.h"       /* definition of type N_Vector, machEnvType      */
#include "dense.h"         /* definition of DenseMat                        */
#include "band.h"          /* definition of BandMat                         */

/* Prototypes: Functions Called by the CVODE Solver */

void CVf(realtype t, N_Vector y, N_Vector ydot, void *f_data);

void CVDenseJac(integertype N, DenseMat J, realtype t, 
                N_Vector y, N_Vector fy, void *jac_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

void CVBandJac(integertype N, integertype mupper, integertype mlower,
               BandMat J, realtype t, N_Vector y, N_Vector fy,
               void *jac_data,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int CVPreco(realtype tn, N_Vector y,N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma, N_Vector ewt, realtype h,
            realtype uround, long int *nfePtr, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int CVPSol(realtype tn, N_Vector y, N_Vector fy, N_Vector vtemp,
           realtype gamma, N_Vector ewt, realtype delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z);

int CVJtimes(N_Vector v, N_Vector Jv, realtype t, 
             N_Vector y, N_Vector fy,
             void *jac_data, N_Vector work);


/* Declarations for global variables, shared among various routines */

void *CV_cvodemem;
N_Vector CV_yvec;

#endif
