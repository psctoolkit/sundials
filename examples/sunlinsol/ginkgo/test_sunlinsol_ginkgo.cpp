/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_ginkgo.hpp>

#include "test_sunlinsol.h"

#if defined(USE_HIP)
#define OMP_OR_HIP_OR_CUDA(a, b, c) b
#elif defined(USE_CUDA)
#define OMP_OR_HIP_OR_CUDA(a, b, c) c
#else
#define OMP_OR_HIP_OR_CUDA(a, b, c) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#else
#include <nvector/nvector_serial.h>
#endif

using GkoMatrixType = gko::matrix::Csr<sunrealtype>;
using SUNMatrixType = sundials::ginkgo::Matrix<GkoMatrixType>;
using GkoSolverType = gko::solver::Gmres<sunrealtype>;
using SUNLinearSolverType =
    sundials::ginkgo::LinearSolver<GkoSolverType, SUNMatrixType>;

/* ----------------------------------------------------------------------
 * SUNLinSol_Ginkgo Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;           /* counter for test failures  */
  sunindextype cols, rows; /* matrix columns, rows       */
  sunindextype nblocks;    /* number of matrix blocks    */
  N_Vector x, y, b;        /* test vectors               */
  int print_timing;
  realtype* xdata;
  sundials::Context sunctx;

  auto gko_exec =
      OMP_OR_HIP_OR_CUDA(gko::OmpExecutor::create(),
                         gko::HipExecutor::create(0, gko::OmpExecutor::create(),
                                                  true),
                         gko::CudaExecutor::create(0, gko::OmpExecutor::create(),
                                                   true));

  /* check input and set matrix dimensions */
  if (argc < 4)
  {
    printf("ERROR: THREE (3) Inputs required: matrix cols, number of blocks, "
           "print timing \n");
    return (-1);
  }

  cols = (sunindextype)atol(argv[1]);
  if (cols <= 0)
  {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return (-1);
  }
  rows = cols;

  nblocks = (sunindextype)atol(argv[2]);
  if (nblocks <= 0)
  {
    printf("ERROR: number of blocks must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\nGinkgo linear solver test: size %ld, blocks %ld\n\n",
         (long int)cols, (long int)nblocks);

  /* Create vectors and matrices */
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, rows);

  x = OMP_OR_HIP_OR_CUDA(N_VNew_Serial(cols, sunctx), N_VNew_Hip(cols, sunctx),
                         N_VNew_Cuda(cols, sunctx));
  b = N_VClone(x);
  y = N_VClone(x);

  auto matrix_dim = gko::dim<2>(rows, cols);
  auto gko_matdata =
      gko::matrix_data<sunrealtype>(matrix_dim, distribution, generator);
  auto gko_matrixA = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_matrixB = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
  auto gko_ident   = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));

  /* Fill matrices */
  gko_matrixA->read(gko_matdata);
  gko_matrixB->read(gko_matdata);
  gko_ident->read(gko::matrix_data<sunrealtype>::diag(matrix_dim, 1.0));

  /* Fill x with random data */
  xdata = N_VGetArrayPointer(x);
  for (sunindextype i = 0; i < cols; i++)
  {
    xdata[i] = distribution(generator);
  }
  OMP_OR_HIP_OR_CUDA(, N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  /* Wrap ginkgo matrices for SUNDIALS.
     sundials::ginkgo::Matrix is overloaded to a SUNMatrix. */
  SUNMatrixType A{gko_matrixA, sunctx};
  SUNMatrixType B{gko_matrixB, sunctx};

  /* You can also create the SUNMatrix via cloning */
  SUNMatrix I = SUNMatClone(A);

  /* Copy A and x into B and y to print in case of solver failure */
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails)
  {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    SUNMatDestroy(I);
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);

    return (1);
  }

  /*
   * Create linear solver.
   */

  /* Use default stopping critieria */
  auto crit = sundials::ginkgo::DefaultStop::build() //
                  .with_max_iters(cols)              //
                  .on(gko_exec);

  auto precon = gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build() //
                                .on(gko_exec);

  auto gko_solver_factory = GkoSolverType::build()               //
                                .with_criteria(gko::share(crit)) //
                                .with_preconditioner(gko::share(precon)) //
                                .on(gko_exec);

  SUNLinearSolverType LS{gko::share(gko_solver_factory), sunctx};
  LS.logResNorm(true); // turn on residual norm logging

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS.get(), SUNLINEARSOLVER_GINKGO, 0);
  fails += Test_SUNLinSolGetType(LS.get(), SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS.get(), 0);
  fails += Test_SUNLinSolSetup(LS.get(), A, 0);
  fails += Test_SUNLinSolSolve(LS.get(), A, x, b, SUN_RCONST(1e-10), SUNTRUE, 0);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nx (original) =\n");
    N_VPrint(y);
    printf("\nx (computed) =\n");
    N_VPrint(x);
    printf("\nb =\n");
    N_VPrint(b);
  }
  else
  {
    printf("\nSUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Print solve information */
  printf("Number of linear solver iterations: %ld\n",
         static_cast<long int>(SUNLinSolNumIters(LS)));
  printf("Final residual norm: %g\n", SUNLinSolResNorm(LS));

  /* Free solver, matrix and vectors */
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector expected, N_Vector actual, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(actual),
                     N_VCopyFromDevice_Cuda(actual));
  OMP_OR_HIP_OR_CUDA(, N_VCopyFromDevice_Hip(expected),
                     N_VCopyFromDevice_Cuda(expected));

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);

  if (xldata != yldata)
  {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
  {
    printf("Check_vector failures:\n");
    for (i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int)i, xdata[i],
               ydata[i], SUNRabs(xdata[i] - ydata[i]));
  }

  return failure > 0;
}

void sync_device()
{
  OMP_OR_HIP_OR_CUDA(, hipDeviceSynchronize(), cudaDeviceSynchronize());
}
