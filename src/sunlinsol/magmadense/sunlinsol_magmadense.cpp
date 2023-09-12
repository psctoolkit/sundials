/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#include <sunmatrix/sunmatrix_magmadense.h>

/* Interfaces to match 'realtype' with the correct MAGMA functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgesv          magma_dgesv_gpu
#define xgetrf         magma_dgetrf_gpu
#define xgetrf_batched magma_dgetrf_batched
#define xgetrs         magma_dgetrs_gpu
#define xgetrs_batched magma_dgetrs_batched
#define xset_pointer   magma_dset_pointer
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgesv          magma_sgesv_gpu
#define xgetrf         magma_sgetrf_gpu
#define xgetrf_batched magma_sgetrf_batched
#define xgetrs         magma_sgetrs_gpu
#define xgetrs_batched magma_sgetrs_batched
#define xset_pointer   magma_sset_pointer
#else
#error Incompatible realtype for MAGMA
#endif

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * MAGMADENSE solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define MAGMADENSE_CONTENT(S) ((SUNLinearSolverContent_MagmaDense)(S->content))
#define MHELP(S)              (MAGMADENSE_CONTENT(S)->memhelp)
#define QUEUE(S)              (MAGMADENSE_CONTENT(S)->q)
#define STRATEGY(S)           (MAGMADENSE_CONTENT(S)->strategy)
#define INITIALIZED(S)        (MAGMADENSE_CONTENT(S)->initialized)
#define ORDER(S)              (MAGMADENSE_CONTENT(S)->N)
#define NBLOCKS(S)            (MAGMADENSE_CONTENT(S)->nblocks)
#define PIVOTS(S)             ((sunindextype*)MAGMADENSE_CONTENT(S)->pivots->ptr)
#define DPIVOTS(S)            ((sunindextype*)MAGMADENSE_CONTENT(S)->pivots->ptr)
#define PIVOTSARRAY(S)        ((sunindextype**)MAGMADENSE_CONTENT(S)->pivotsarr->ptr)
#define RHSARRAY(S)           ((realtype**)MAGMADENSE_CONTENT(S)->rhsarr->ptr)
#define WORKDOUBLE(S)         ((double*)MAGMADENSE_CONTENT(S)->workd->ptr)
#define WORKSINGLE(S)         ((float*)MAGMADENSE_CONTENT(S)->works->ptr)
#define INFOARRAY(S)          ((sunindextype*)MAGMADENSE_CONTENT(S)->infoarr->ptr)
#define LASTFLAG(S)           (MAGMADENSE_CONTENT(S)->last_flag)
#define ASYNCHRONOUS(S)       (MAGMADENSE_CONTENT(S)->async)

/*
 * ----------------------------------------------------------------------------
 * Implementation specific routines
 * ----------------------------------------------------------------------------
 */

/*
 * Constructor functions
 */

SUNLinearSolver SUNLinSol_MagmaDense(N_Vector y, SUNMatrix Amat, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_MagmaDense content;
  SUNMatrixContent_MagmaDense A;
  sunindextype M, nblocks;

  /* Check inputs */
  if (y == NULL || Amat == NULL) return (NULL);

  if (y->ops == NULL || Amat->ops == NULL) return (NULL);

  if (y->ops->nvgetlength == NULL || y->ops->nvgetdevicearraypointer == NULL ||
      Amat->ops->getid == NULL)
    return (NULL);

  /* Check compatibility with supplied SUNMatrix */
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return (NULL);

  if (Amat->content == NULL) return (NULL);

  A = (SUNMatrixContent_MagmaDense)Amat->content;

  /* Check that the matrix is square */
  if (A->M != A->N) return (NULL);

  M       = A->M;
  nblocks = A->nblocks;

  /* Check that the matirx and vector dimensions agree */
  if (M * nblocks != N_VGetLength(y)) return (NULL);

  /* Create the linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) return (NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_MagmaDense;
  S->ops->getid      = SUNLinSolGetID_MagmaDense;
  S->ops->initialize = SUNLinSolInitialize_MagmaDense;
  S->ops->setup      = SUNLinSolSetup_MagmaDense;
  S->ops->solve      = SUNLinSolSolve_MagmaDense;
  S->ops->lastflag   = SUNLinSolLastFlag_MagmaDense;
  S->ops->space      = SUNLinSolSpace_MagmaDense;
  S->ops->free       = SUNLinSolFree_MagmaDense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_MagmaDense)malloc(sizeof(*content));
  if (content == NULL)
  {
    SUNLinSolFree(S);
    return (NULL);
  }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag   = 0;
  content->initialized = 0;
  content->strategy    = SUN_MAGMA_STRATEGY_DEFAULT;
  content->async       = SUNTRUE;
  content->N           = M;
  content->nblocks     = nblocks;
  content->pivots      = NULL;
  content->dpivots     = NULL;
  content->pivotsarr   = NULL;
  content->infoarr     = NULL;
  content->rhsarr      = NULL;
  content->works       = NULL;
  content->workd       = NULL;
  content->memhelp     = A->memhelp;
  content->q           = A->q;

  return (S);
}

/*
 * Set functions
 */

int SUNLinSol_MagmaDense_SetAsync(SUNLinearSolver S, booleantype onoff)
{
  if (S == NULL) return SUNLS_MEM_NULL;
  ASYNCHRONOUS(S) = onoff;
  return SUNLS_SUCCESS;
}

int SUNLinSol_MagmaDense_SetStrategy(SUNLinearSolver S, int strategy)
{
  if (S == NULL) return SUNLS_MEM_NULL;
  if (INITIALIZED(S)) return SUNLS_ILL_INPUT;
#if !defined(SUNDIALS_SINGLE_PRECISION)
  if (strategy == SUN_MAGMA_STRATEGY_HTGETRF) { return SUNLS_ILL_INPUT; }
#endif
  if (NBLOCKS(S) > 1 && STRATEGY(S) != SUN_MAGMA_STRATEGY_GETRF)
  {
    return SUNLS_ILL_INPUT;
  }
  STRATEGY(S) = strategy;
  return SUNLS_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Implementation of generic SUNLinearSolver operations.
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_MagmaDense(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_MagmaDense(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_MAGMADENSE);
}

int SUNLinSolInitialize_MagmaDense(SUNLinearSolver S)
{
  int retval = 0;

  /* All solver-specific memory has already been allocated */
  if (S == NULL) return SUNLS_MEM_NULL;

  if (!INITIALIZED(S))
  {
    /* Allocate data */

    /* The pivots need to be in host memory when calling the non-batched
       methods, but in device memory for the batched methods. */
    retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->pivots,
                                   ORDER(S) * NBLOCKS(S) * sizeof(sunindextype),
                                   SUNMEMTYPE_HOST, nullptr);
    if (retval)
    {
      SUNLinSolFree(S);
      return (SUNLS_MEM_FAIL);
    }

    /* If we have multiple blocks, then we need to allocate some extra
      pointer arrays needed when calling MAGMA batched methods. */
    if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF && NBLOCKS(S) > 1)
    {
      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->pivots,
                                     ORDER(S) * NBLOCKS(S) * sizeof(sunindextype),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }

      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->pivotsarr,
                                     NBLOCKS(S) * sizeof(sunindextype*),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }

      /* Set the pivots array on the device */
      magma_iset_pointer(PIVOTSARRAY(S), /* 2D output array */
                         PIVOTS(S),      /* 1D input array */
                         1,              /* leading dimension of input */
                         0,              /* row */
                         0,              /* column */
                         ORDER(S),       /* rows in a block */
                         NBLOCKS(S),     /* number of blocks */
                         QUEUE(S));

      /* We use pinned memory for the info array because we are going to
        check its values on the host and we need it to have fast transfers. */
      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->infoarr,
                                     NBLOCKS(S) * sizeof(sunindextype),
                                     SUNMEMTYPE_PINNED, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }

      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->rhsarr,
                                     NBLOCKS(S) * sizeof(realtype*),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }
    }

    if (STRATEGY(S) == SUN_MAGMA_STRATEGY_HGESV_REFINE ||
        STRATEGY(S) == SUN_MAGMA_STRATEGY_SGESV_REFINE)
    {
      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->dpivots,
                                     ORDER(S) * NBLOCKS(S) * sizeof(sunindextype),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }

      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->workd,
                                     ORDER(S) * sizeof(double),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }

      retval = SUNMemoryHelper_Alloc(MHELP(S), &MAGMADENSE_CONTENT(S)->works,
                                     ORDER(S) * ORDER(S) * sizeof(float),
                                     SUNMEMTYPE_DEVICE, nullptr);
      if (retval)
      {
        SUNLinSolFree(S);
        return (SUNLS_MEM_FAIL);
      }
    }
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (SUNLS_SUCCESS);
}

int SUNLinSolSetup_MagmaDense(SUNLinearSolver S, SUNMatrix A)
{
  /* Check for valid inputs */
  if (S == NULL) return SUNLS_MEM_NULL;

  if (A == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_NULL;
    return (SUNLS_MEM_NULL);
  }

  /* Ensure that A is a magma dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_MAGMADENSE)
  {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return (SUNLS_ILL_INPUT);
  }

  sunindextype ier = 0;
  sunindextype M   = SUNMatrix_MagmaDense_BlockRows(A);

  if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF && NBLOCKS(S) > 1)
  {
    xgetrf_batched(M, /* number of rows per block */
                   M, /* number of columns per block */
                   SUNMatrix_MagmaDense_BlockData(A), M, /* leading dimension of
                                                            each block */
                   PIVOTSARRAY(S), INFOARRAY(S),
                   SUNMatrix_MagmaDense_NumBlocks(A), QUEUE(S));

    if (!ASYNCHRONOUS(S))
    {
      magma_queue_sync(QUEUE(S));
      /* Check if there were any failures when factoring */
      for (sunindextype k = 0; k < SUNMatrix_MagmaDense_NumBlocks(A); k++)
      {
        if (INFOARRAY(S)[k] < 0) ier = INFOARRAY(S)[k];
        if (INFOARRAY(S)[k] > 0)
        {
          ier = INFOARRAY(S)[k];
          break;
        }
      }
    }
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF)
  {
    if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF)
    {
      xgetrf(M,                               /* number of rows */
             M,                               /* number of columns */
             SUNMatrix_MagmaDense_Data(A), M, /* leading dimension of A */
             PIVOTS(S), &ier);
    }
    if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_HTGETRF)
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    magma_htgetrf_gpu(M,                               /* number of rows */
                      M,                               /* number of columns */
                      SUNMatrix_MagmaDense_Data(A), M, /* leading dimension of A
                                                        */
                      PIVOTS(S), &ier);

    if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));
#else
    ier = SUNLS_LUFACT_FAIL;
#endif
  }

  LASTFLAG(S) = ier;
  if (ier > 0) return (SUNLS_LUFACT_FAIL);
  if (ier < 0) return (SUNLS_PACKAGE_FAIL_UNREC);
  return (SUNLS_SUCCESS);
}

int SUNLinSolSolve_MagmaDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                              N_Vector b, realtype tol)
{
  /* Check for valid inputs */
  if (S == NULL) return (SUNLS_MEM_NULL);

  if ((A == NULL) || (x == NULL) || (b == NULL))
  {
    LASTFLAG(S) = SUNLS_MEM_NULL;
    return (SUNLS_MEM_NULL);
  }

  /* Ensure that A is a magma dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_MAGMADENSE)
  {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return (SUNLS_ILL_INPUT);
  }

  int ier              = 0;
  sunindextype M       = SUNMatrix_MagmaDense_BlockRows(A);
  sunindextype nblocks = SUNMatrix_MagmaDense_NumBlocks(A);
  realtype* bdata      = N_VGetDeviceArrayPointer(b);
  if (bdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return (SUNLS_MEM_FAIL);
  }
  realtype* xdata = N_VGetDeviceArrayPointer(x);
  if (xdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return (SUNLS_MEM_FAIL);
  }

  /* Copy b into x for getrs which uses x as input and output */
  if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GESV ||
      STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF ||
      STRATEGY(S) == SUN_MAGMA_STRATEGY_HTGETRF)
  {
    N_VScale(ONE, b, x);
  }

  /* Call MAGMA to solve the linear system */
  if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF && NBLOCKS(S) > 1)
  {
    /* First, set pointers to RHS blocks */
    xset_pointer(RHSARRAY(S), /* 2D output array */
                 xdata,       /* 1D input array */
                 1,           /* leading dimension of input */
                 0,           /* rows */
                 0,           /* cols */
                 M,           /* number of rows in block */
                 nblocks, QUEUE(S));

    /* Now, solve the batch system */
    xgetrs_batched(MagmaNoTrans, M, /* order of the matrix */
                   1,               /* number of right hand sides */
                   SUNMatrix_MagmaDense_BlockData(A), M, /* leading dimension of
                                                            A */
                   PIVOTSARRAY(S), RHSARRAY(S), /* right hand side (input),
                                                   solution (output) */
                   M,                           /* leading dimension of b */
                   nblocks, QUEUE(S));
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GETRF)
  {
    xgetrs(MagmaNoTrans, M,                 /* order of the matrix */
           1,                               /* number of right hand sides */
           SUNMatrix_MagmaDense_Data(A), M, /* leading dimension of A */
           PIVOTS(S), xdata, /* right hand side (input), solution (output) */
           M,                /* leading dimension of x */
           &ier);
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_GESV)
  {
    xgesv(M, 1, SUNMatrix_MagmaDense_Data(A), M, PIVOTS(S), xdata, M, &ier);
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_SGESV_REFINE)
  {
    int iter = 0;

    magma_dsgesv_iteref_gpu(MagmaNoTrans, M, /* order of the matrix */
                            1,               /* number of right hand sides */
                            SUNMatrix_MagmaDense_Data(A), M, /* lddA */
                            PIVOTS(S), DPIVOTS(S), bdata,    /* right hand side
                                                                (input) */
                            M,                               /* lddB */
                            xdata, /* solution (output) */
                            M,     /* lddX */
                            WORKDOUBLE(S), WORKSINGLE(S), &iter, &ier);

    fprintf(stderr, "Num. refinement iters = %d\n", iter);
  }
  else if (STRATEGY(S) == SUN_MAGMA_STRATEGY_HGESV_REFINE)
  {
    int iter = 0;

    magma_dhgesv_iteref_gpu(MagmaNoTrans, M, /* order of the matrix */
                            1,               /* number of right hand sides */
                            SUNMatrix_MagmaDense_Data(A), M, /* lddA */
                            PIVOTS(S), DPIVOTS(S), bdata,    /* right hand side
                                                                (input) */
                            M,                               /* lddB */
                            xdata, /* solution (output) */
                            M,     /* lddX */
                            WORKDOUBLE(S), WORKSINGLE(S), &iter, &ier);

    fprintf(stderr, "Num. refinement iters = %d\n", iter);
  }
  else { return -1; }

  if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));

  LASTFLAG(S) = ier;
  return ((ier < 0) ? SUNLS_PACKAGE_FAIL_UNREC : SUNLS_SUCCESS);
}

sunindextype SUNLinSolLastFlag_MagmaDense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return (-1);
  return (LASTFLAG(S));
}

int SUNLinSolSpace_MagmaDense(SUNLinearSolver S, long int* lenrwLS,
                              long int* leniwLS)
{
  *lenrwLS = 0;
  /* TODO: the strategy updates make this inaccurate */
  *leniwLS = 2 + MAGMADENSE_CONTENT(S)->N;
  return (SUNLS_SUCCESS);
}

int SUNLinSolFree_MagmaDense(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL) return (SUNLS_SUCCESS);

  /* delete items from contents, then delete generic structure */
  if (S->content)
  {
    if (MAGMADENSE_CONTENT(S)->pivots)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->pivots, nullptr);
    if (MAGMADENSE_CONTENT(S)->dpivots)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->dpivots, nullptr);
    if (MAGMADENSE_CONTENT(S)->pivotsarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->pivotsarr,
                              nullptr);
    if (MAGMADENSE_CONTENT(S)->infoarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->infoarr, nullptr);
    if (MAGMADENSE_CONTENT(S)->rhsarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->rhsarr, nullptr);
    if (MAGMADENSE_CONTENT(S)->works)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->works, nullptr);
    if (MAGMADENSE_CONTENT(S)->workd)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->workd, nullptr);
    free(S->content);
    S->content = NULL;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return (SUNLS_SUCCESS);
}
