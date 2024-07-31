/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test for GetUserData functions
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"

static void err_fn(int line, const char* func, const char* file, const char* msg,
                   SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  fprintf(stderr, "Error at line %i of %s in %s: %s\n", line, func, file, msg);
  exit(err_code);
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(SUN_RCONST(0.0), ydot);
  return 0;
}

static int check_step(void* arkode_mem, N_Vector y, sunrealtype h_expected)
{
  static int step = 0;
  step++;

  sunrealtype tret;
  ARKodeEvolve(arkode_mem, SUN_RCONST(1.0), y, &tret, ARK_ONE_STEP);

  long int local_err_fails;
  ARKodeGetNumErrTestFails(arkode_mem, &local_err_fails);
  if (local_err_fails != 0)
  {
    fprintf(stderr, "Expected 0 local error failures at step %i but is %li\n",
            step, local_err_fails);
  }

  const N_Vector err = N_VClone(y);
  ARKodeGetEstLocalErrors(arkode_mem, err);
  sunrealtype err_norm = N_VMaxNorm(err);
  N_VDestroy(err);

  if (err_norm != 0)
  {
    fprintf(stderr, "Expected local error at step %i to be 0 but is %g\n", step,
            err_norm);
    return 1;
  }

  sunrealtype h_actual;
  ARKodeGetCurrentStep(arkode_mem, &h_actual);
  if (h_expected != h_actual)
  {
    fprintf(stderr, "Expected h at step %i to be %g but is %g\n", step,
            h_expected, h_actual);
    return 1;
  }

  return 0;
}

int main()
{
  SUNContext sunctx;
  int retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  retval = SUNContext_PushErrHandler(sunctx, err_fn, NULL);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_PushErrHandler returned %i\n", retval);
    return 1;
  }

  const N_Vector y = N_VNew_Serial(1, sunctx);
  N_VConst(SUN_RCONST(1.0), y);

  void* arkode_mem = ERKStepCreate(f, SUN_RCONST(0.0), y, sunctx);

  const sunrealtype h0           = SUN_RCONST(1.0e-4);
  const sunrealtype first_growth = SUN_RCONST(1234.0);
  const sunrealtype growth       = SUN_RCONST(3.0);

  ARKodeSetInitStep(arkode_mem, h0);
  ARKodeSetMaxFirstGrowth(arkode_mem, first_growth);
  ARKodeSetMaxGrowth(arkode_mem, growth);

  sunrealtype h_expect = first_growth * h0;

  for (int i = 0; i < 4; i++)
  {
    retval = check_step(arkode_mem, y, h_expect);
    if (retval != 0) { return retval; }
    h_expect *= growth;
  }

  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}
