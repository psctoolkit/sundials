/* clang-format off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * We consider the Kepler problem. We choose one body to be the center of our
 * coordinate system and then we use the coordinates q = (q1, q2) to represent
 * the position of the second body relative to the first (center). This yields
 * the ODE:
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) ]
 * with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
 *
 * By default We solve the problem by letting y = [ q, p ]^T then using a 4th
 * order symplectic integrator via the SPRKStep time-stepper of ARKODE with a
 * fixed time-step size.
 *
 * The program also accepts command line arguments to change the method
 * used and time-stepping strategy. The program has the following CLI arguments:
 * 
 *   --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)
 *   --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)
 *   --method <string>           which method to use (default ARKODE_SYMPLECTIC_MCLACHLAN_4_4)
 *   --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)
 *   --tf <Real>                 the final time for the simulation (default 100)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where applicable
 * 
 * References:
 *    Ernst Hairer, Christain Lubich, Gerhard Wanner
 *    Geometric Numerical Integration: Structure-Preserving
 *    Algorithms for Ordinary Differential Equations
 *    Springer, 2006,
 *    ISSN 0179-3632
 * --------------------------------------------------------------------------*/
/* clang-format on */

#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for MRIStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "arkode/arkode.h"

typedef struct
{
  sunrealtype ecc;
} * UserData;

typedef struct
{
  int step_mode;
  int stepper;
  int use_compsums;
  sunrealtype dt;
  sunrealtype tf;
  const char* method_name;
} ProgramArgs;

static int check_retval(void* returnvalue, const char* funcname, int opt);
static int ParseArgs(int argc, char* argv[], ProgramArgs* args);
static void PrintHelp();

static void InitialConditions(N_Vector y0, sunrealtype ecc);
static sunrealtype Hamiltonian(N_Vector yvec);
static sunrealtype AngularMomentum(N_Vector y);

static int dydt(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

static sunrealtype Q(N_Vector yvec, sunrealtype alpha);
static sunrealtype G(N_Vector yvec, sunrealtype alpha);
static int Adapt(N_Vector y, sunrealtype t, sunrealtype h1, sunrealtype h2,
                 sunrealtype h3, sunrealtype e1, sunrealtype e2, sunrealtype e3,
                 int q, int p, sunrealtype* hnew, void* user_data);

int main(int argc, char* argv[])
{
  ProgramArgs args;
  void* arkode_mem       = NULL;
  SUNContext sunctx      = NULL;
  N_Vector y             = NULL;
  SUNNonlinearSolver NLS = NULL;
  UserData udata         = NULL;
  sunrealtype tout       = NAN;
  sunrealtype tret       = NAN;
  sunrealtype H0         = NAN;
  sunrealtype L0         = NAN;
  FILE* conserved_fp     = NULL;
  FILE* solution_fp      = NULL;
  FILE* times_fp         = NULL;
  int rootsfound         = 0;
  int argi               = 0;
  int iout               = 0;
  int retval             = 0;

  /* CLI args */
  if (ParseArgs(argc, argv, &args)) { return 1; };
  const int step_mode     = args.step_mode;
  const int stepper       = args.stepper;
  const char* method_name = args.method_name;
  const int use_compsums  = args.use_compsums;
  const sunrealtype dt    = args.dt;
  sunrealtype Tf          = args.tf;

  /* Default problem parameters */
  const sunrealtype T0       = SUN_RCONST(0.0);
  const sunrealtype ecc      = SUN_RCONST(0.6);
  const sunrealtype dTout    = 10 * dt; // SUN_RCONST(1.);
  const int num_output_times = (int)ceil(Tf / dTout);

  printf("\n   Begin Kepler Problem\n\n");

  /* Allocate and fill udata structure */
  udata      = (UserData)malloc(sizeof(*udata));
  udata->ecc = ecc;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* Allocate our state vector */
  y = N_VNew_Serial(4, sunctx);

  /* Fill the initial conditions */
  InitialConditions(y, ecc);

  /* Create SPRKStep integrator */
  if (stepper == 0)
  {
    arkode_mem = SPRKStepCreate(force, velocity, T0, y, sunctx);

    /* Optional: enable temporal root-finding */
    SPRKStepRootInit(arkode_mem, 1, rootfn);
    if (check_retval(&retval, "SPRKStepRootInit", 1)) return 1;

    retval = SPRKStepSetMethodName(arkode_mem, method_name);
    if (check_retval(&retval, "SPRKStepSetMethodName", 1)) return 1;

    retval = SPRKStepSetUseCompensatedSums(arkode_mem, use_compsums);
    if (check_retval(&retval, "SPRKStepSetUseCompensatedSums", 1)) return 1;

    if (step_mode == 0)
    {
      retval = SPRKStepSetFixedStep(arkode_mem, dt);
      if (check_retval(&retval, "SPRKStepSetFixedStep", 1)) return 1;

      retval = SPRKStepSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 1);
      if (check_retval(&retval, "SPRKStepSetMaxNumSteps", 1)) return 1;
    }
    else
    {
      fprintf(stderr,
              "ERROR: adaptive time-steps are not supported with SPRKStep\n");
      return 1;
    }

    retval = SPRKStepSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "SPRKStepSetUserData", 1)) return 1;
  }
  else if (stepper == 1)
  {
    arkode_mem = ARKStepCreate(dydt, NULL, T0, y, sunctx);

    retval = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method_name);
    if (check_retval(&retval, "ARKStepSetTableName", 1)) return 1;

    ARKStepRootInit(arkode_mem, 1, rootfn);
    if (check_retval(&retval, "ARKStepRootInit", 1)) return 1;

    retval = ARKStepSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

    retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
    if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return 1;

    if (step_mode == 0) { retval = ARKStepSetFixedStep(arkode_mem, dt); }
    else
    {
      retval = ARKStepSStolerances(arkode_mem, dt, dt);
      if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
    }
  }

  /* Open output files */
  if (stepper == 0)
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }
  else
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }

  /* Print out starting energy, momentum before integrating */
  tret = T0;
  tout = T0 + dTout;
  H0   = Hamiltonian(y);
  L0   = AngularMomentum(y);
  fprintf(stdout, "t = %.4Lf, H(p,q) = %.16Lf, L(p,q) = %.16Lf\n", tret, H0, L0);
  fprintf(times_fp, "%.16Lf\n", tret);
  fprintf(conserved_fp, "%.16Lf, %.16Lf\n", H0, L0);
  N_VPrintFile(y, solution_fp);

  /* Do integration */
  if (stepper == 0)
  {
    for (iout = 0; iout < num_output_times; iout++)
    {
      sunrealtype hlast = SUN_RCONST(0.0);

      /* Optional: if the stop time is not set, then its possible that the the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      SPRKStepSetStopTime(arkode_mem, tout);
      retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        fprintf(stdout, "ROOT RETURN:\t");
        SPRKStepGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "t = %.4Lf g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg\n", tret,
                rootsfound, N_VGetArrayPointer(y)[0], N_VGetArrayPointer(y)[1]);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                tret, Hamiltonian(y) - H0, AngularMomentum(y) - L0);

        /* Continue to tout */
        retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
      }

      /* Output current integration status */
      fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
              tret, Hamiltonian(y) - H0, AngularMomentum(y) - L0);
      fprintf(times_fp, "%.16Lf\n", tret);
      fprintf(conserved_fp, "%.16Lf, %.16Lf\n", Hamiltonian(y),
              AngularMomentum(y));

      N_VPrintFile(y, solution_fp);

      /* Check if the solve was successful, if so, update the time and continue
       */
      if (retval >= 0)
      {
        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
      }
      else
      {
        fprintf(stderr, "Solver failure, stopping integration\n");
        break;
      }
    }
  }
  else
  {
    for (iout = 0; iout < num_output_times; iout++)
    {
      /* Optional: if the stop time is not set, then its possible that the the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      ARKStepSetStopTime(arkode_mem, tout);
      retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        fprintf(stdout, "ROOT RETURN:\t");
        ARKStepGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "t = %.4Lf g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg\n", tret,
                rootsfound, N_VGetArrayPointer(y)[0], N_VGetArrayPointer(y)[1]);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
                tret, Hamiltonian(y) - H0, AngularMomentum(y) - L0);

        /* Continue to tout */
        retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
      }

      /* Output current integration status */
      fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
              tret, Hamiltonian(y) - H0, AngularMomentum(y) - L0);
      fprintf(times_fp, "%.16Lf\n", tret);
      fprintf(conserved_fp, "%.16Lf, %.16Lf\n", Hamiltonian(y),
              AngularMomentum(y));
      N_VPrintFile(y, solution_fp);

      /* Check if the solve was successful, if so, update the time and continue
       */
      if (retval >= 0)
      {
        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
      }
      else
      {
        fprintf(stderr, "ERROR: Solver failure, stopping integration\n");
        break;
      }
    }
  }

  free(udata);
  fclose(times_fp);
  fclose(conserved_fp);
  fclose(solution_fp);
  if (NLS) { SUNNonlinSolFree(NLS); }
  N_VDestroy(y);
  if (stepper == 0)
  {
    SPRKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    SPRKStepFree(&arkode_mem);
  }
  else
  {
    ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    ARKStepFree(&arkode_mem);
  }

  SUNContext_Free(&sunctx);

  return 0;
}

void InitialConditions(N_Vector y0vec, sunrealtype ecc)
{
  const sunrealtype zero = SUN_RCONST(0.0);
  const sunrealtype one  = SUN_RCONST(1.0);
  sunrealtype* y0        = N_VGetArrayPointer(y0vec);

  y0[0] = one - ecc;
  y0[1] = zero;
  y0[2] = zero;
  y0[3] = SUNRsqrt((one + ecc) / (one - ecc));
}

void Solution(N_Vector yvec, UserData user_data) {}

sunrealtype Hamiltonian(N_Vector yvec)
{
  sunrealtype H              = 0.0;
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  const sunrealtype sqrt_qTq = SUNRsqrt(y[0] * y[0] + y[1] * y[1]);
  const sunrealtype pTp      = y[2] * y[2] + y[3] * y[3];

  H = SUN_RCONST(0.5) * pTp - SUN_RCONST(1.0) / sqrt_qTq;

  return H;
}

sunrealtype AngularMomentum(N_Vector yvec)
{
  sunrealtype L        = 0.0;
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  L = q1 * p2 - q2 * p1;

  return L;
}

int dydt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  int retval = 0;

  retval += force(t, yvec, ydotvec, user_data);
  retval += velocity(t, yvec, ydotvec, user_data);

  return retval;
}

int velocity(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  sunrealtype* ydot    = N_VGetArrayPointer(ydotvec);
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  ydot[0] = p1;
  ydot[1] = p2;
  // ydot[2] = ydot[3] = SUN_RCONST(0.0);

  return 0;
}

int force(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  UserData udata             = (UserData)user_data;
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  sunrealtype* ydot          = N_VGetArrayPointer(ydotvec);
  const sunrealtype q1       = y[0];
  const sunrealtype q2       = y[1];
  const sunrealtype sqrt_qTq = SUNRsqrt(q1 * q1 + q2 * q2);

  // ydot[0] = ydot[1] = SUN_RCONST(0.0);
  ydot[2] = -q1 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));
  ydot[3] = -q2 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));

  return 0;
}

int rootfn(sunrealtype t, N_Vector yvec, sunrealtype* gout, void* user_data)
{
  UserData udata       = (UserData)user_data;
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];

  /* We want to know when the body crosses the position (0.36, -0.22) */
  gout[0] = (q1 - SUN_RCONST(0.36)) + (q2 + SUN_RCONST(0.22));

  return 0;
}

/* Functions needed to implement the step density integrating controller
   proposed by Hairer and Soderlind in https://doi.org/10.1137/04060699. */

sunrealtype G(N_Vector yvec, sunrealtype alpha)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  const sunrealtype pTq = p1 * q1 + p2 * q2;
  const sunrealtype qTq = q1 * q1 + q2 * q2;

  return (-alpha * pTq / qTq);
}

sunrealtype Q(N_Vector yvec, sunrealtype alpha)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];

  const sunrealtype qTq = q1 * q1 + q2 * q2;

  return SUNRpowerR(qTq, alpha / SUN_RCONST(2.0));
}

int ParseArgs(int argc, char* argv[], ProgramArgs* args)
{
  args->step_mode    = 0;
  args->stepper      = 0;
  args->method_name  = NULL;
  args->use_compsums = 0;
  args->dt           = SUN_RCONST(1e-2);
  args->tf           = SUN_RCONST(100.);

  for (int argi = 1; argi < argc; argi++)
  {
    if (!strcmp(argv[argi], "--step-mode"))
    {
      argi++;
      if (!strcmp(argv[argi], "fixed")) { args->step_mode = 0; }
      else if (!strcmp(argv[argi], "adapt")) { args->step_mode = 1; }
      else
      {
        fprintf(stderr, "ERROR: --step-mode must be 'fixed' or 'adapt'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--stepper"))
    {
      argi++;
      if (!strcmp(argv[argi], "SPRK")) { args->stepper = 0; }
      else if (!strcmp(argv[argi], "ERK")) { args->stepper = 1; }
      else
      {
        fprintf(stderr, "ERROR: --stepper must be 'SPRK' or 'ERK'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--method"))
    {
      argi++;
      args->method_name = argv[argi];
    }
    else if (!strcmp(argv[argi], "--dt"))
    {
      argi++;
      args->dt = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--tf"))
    {
      argi++;
      args->tf = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--use-compensated-sums"))
    {
      args->use_compsums = 1;
    }
    else if (!strcmp(argv[argi], "--help"))
    {
      PrintHelp();
      return 1;
    }
    else
    {
      fprintf(stderr, "ERROR: unrecognized argument %s\n", argv[argi]);
      PrintHelp();
      return 1;
    }
  }

  if (!args->method_name)
  {
    if (args->stepper == 0)
    {
      args->method_name = "ARKODE_SYMPLECTIC_MCLACHLAN_4_4";
    }
    else if (args->stepper == 1)
    {
      args->method_name = "ARKODE_ZONNEVELD_5_3_4";
    }
  }

  return 0;
}

void PrintHelp()
{
  fprintf(stderr, "ark_kepler: an ARKODE example demonstrating the SPRKStep "
                  "time-stepping module solving the Kepler problem\n");
  fprintf(stderr, "  --step-mode <fixed, adapt>  should we use a fixed "
                  "time-step or adaptive time-step (default fixed)\n");
  fprintf(stderr,
          "  --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep "
          "with an ERK method (default SPRK)\n");
  fprintf(stderr, "  --method <string>           which method to use (default "
                  "ARKODE_SYMPLECTIC_MCLACHLAN_4_4)\n");
  fprintf(stderr, "  --dt <Real>                 the fixed-time step size to "
                  "use if fixed "
                  "time stepping is turned on (default 0.01)\n");
  fprintf(stderr, "  --tf <Real>                 the final time for the "
                  "simulation (default 100)\n");
  fprintf(stderr,
          "  --use-compensated-sums      turns on compensated summation in "
          "ARKODE where applicable\n");
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval = NULL;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
