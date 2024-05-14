/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNAdjointSolver class definition.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_SOLVER_H
#define _SUNADJOINT_SOLVER_H

#include <sunadjoint/sunadjoint_checkpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#include "sundials/sundials_types.h"

struct SUNAdjointSolver_
{
  SUNStepper stepper;
  SUNJacFn JacFn, JacPFn;
  SUNMatrix Jac, JacP;
  SUNJacTimesFn Jvp, JPvp, vJp, vJPp;
  SUNAdjointCheckpointScheme checkpoint_scheme;
  sunrealtype tf;
  SUNContext sunctx;
};

typedef struct SUNAdjointSolver_* SUNAdjointSolver;

#ifdef __cplusplus
extern "C" {
#endif

// TODO(CJB): I think this should be a private function that is only used
// within the package CreateAdjointSolver routines.
SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Create(SUNStepper stepper,
                                   sunindextype num_cost_fns, N_Vector sf, sunrealtype tf,
                                   SUNAdjointCheckpointScheme checkpoint_scheme,
                                   SUNContext sunctx,
                                   SUNAdjointSolver* adj_solver);

/*
  Solves the adjoint system.

  :param adj_solver: The adjoint solver object.
  :param tout: The time at which the adjoint solution is desired.
  :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
  :param tret: On return, the time reached by the adjoint solver.
  :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

  :returns: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Solve(SUNAdjointSolver adj_solver,
                                  sunrealtype tout, N_Vector sens,
                                  sunrealtype* tret, int* stop_reason);


/*
  Evolves the adjoint system backwards one step.

  :param adj_solver: The adjoint solver object.
  :param tout: The time at which the adjoint solution is desired.
  :param sens: The vector of sensitivity solutions dg/dy0 and dg/dp.
  :param tret: On return, the time reached by the adjoint solver.
  :param stop_reason: On return, an integer code that indicates why the adjoint solver stopped.

  :returns: A SUNErrCode indicating failure or success.
 */
SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Step(SUNAdjointSolver adj_solver,
                                 sunrealtype tout, N_Vector sens,
                                 sunrealtype* tret, int* stop_reason);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetJacFn(SUNAdjointSolver, SUNJacFn JacFn, SUNMatrix Jac);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetJacPFn(SUNAdjointSolver, SUNJacFn JacPFn, SUNMatrix JacP);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetJacTimesVecFn(SUNAdjointSolver, SUNJacTimesFn Jvp,
                                             SUNJacTimesFn JPvp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_SetVecTimesJacFn(SUNAdjointSolver, SUNJacTimesFn vJp,
                                             SUNJacTimesFn vJPp);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointSolver_Destroy(SUNAdjointSolver*);

#ifdef __cplusplus
}
#endif
#endif /* _SUNADJOINT_SOLVER_H */
