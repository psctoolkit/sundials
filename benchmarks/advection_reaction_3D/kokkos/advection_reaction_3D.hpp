/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                David J. Gardner, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------*/

#ifndef ADVECTION_REACTION_3D_HPP
#define ADVECTION_REACTION_3D_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mpi.h>

#include <sundials/sundials_context.h>
#include <nvector/nvector_mpiplusx.h>
#include "nvector/nvector_kokkos.hpp"
#include "check_retval.h"
#include "ParallelGrid.hpp"

/* Set SUNDIALS Kokkos vector shortcut */
using SUNVector = sundials::kokkos::Vector<ExecSpace>;

using sundials_tools::ParallelGrid;
using sundials_tools::BoundaryType;
using sundials_tools::StencilType;
using std::string;

/* Number of dimensions */
constexpr int NDIMS = 3;

/* Maximum size of output directory string */
constexpr int MXSTR = 2048;

/*
 * Data structure for problem options
 */

struct UserOptions
{
  int      npxyz[3]; /* number of processors in x,y,z */
  sunindextype npts; /* number of spatial mesh points */
  realtype t0;       /* initial time                  */
  realtype tf;       /* final time                    */
  realtype rtol;     /* relative tolerance            */
  realtype atol;     /* absolute tolerance            */
  int      order;    /* method order                  */
  string   method;   /* method string                 */
  string   nls;      /* nonlinear solver to use       */
  int      fpaccel;  /* number of fixedpoint vectors  */
  int      precond;  /* to precondition or not        */
  int      fused;    /* use fused vector ops          */
  int      nout;     /* number of outputs             */
  int      save;     /* save solution to disk         */
  char*    outputdir;
};


/*
 * Data structure for problem specific data
 */

struct UserData
{
  SUNContext ctx;
  SUNProfiler prof;

  /* MPI data */
  MPI_Comm    comm;
  int         myid;
  int         nprocs;
  MPI_Request req[2];

  /* Should reactions be added to the advection or not */
  bool add_reactions;

  /* File handles for output */
  FILE*  TFID;     /* time output file pointer     */
  FILE*  UFID;     /* solution output file pointer */
  FILE*  VFID;
  FILE*  WFID;

  /* Solution masks */
  N_Vector umask;
  N_Vector vmask;
  N_Vector wmask;

  /* Problem parameters */
  realtype  xmax; /* maximum x value              */
  realtype  A;    /* concentration of species A   */
  realtype  B;    /* w source rate                */
  realtype  k1;   /* reaction rates               */
  realtype  k2;
  realtype  k3;
  realtype  k4;
  realtype  k5;
  realtype  k6;
  realtype  c;    /* advection coefficient        */

  /* Parallel mesh */
  ParallelGrid<sunindextype,NDIMS>* grid;

  /* Count of implicit function evals by the task local nonlinear solver */
  long int nnlfi;

  /* Integrator options */
  UserOptions* uopt;

  /* Constructor that takes the context */
  UserData(SUNContext ctx) : ctx(ctx) {
    SUNContext_GetProfiler(ctx, &prof);
    umask = nullptr;
    vmask = nullptr;
    wmask = nullptr;
    TFID = nullptr;
    UFID = nullptr;
    VFID = nullptr;
    WFID = nullptr;
    uopt = nullptr;
  }

  /* destructor frees the problem data */
  ~UserData();
};


/*
 * Functions to evolve the solution (defined by the drivers)
 */

/* function that does ARKStep setup and evolves the solution with a DIRK method */
extern int EvolveProblemDIRK(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does ARKStep setup and evolves the solution with an IMEX method */
extern int EvolveProblemIMEX(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does ERKStep setup and evolves the solution */
extern int EvolveProblemExplicit(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does CVODE BDF setup and evolves the solution */
extern int EvolveProblemBDF(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does CVODE Adams setup and evolves the solution */
extern int EvolveProblemAdams(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does IDA BDF setup and evolves the solution */
extern int EvolveDAEProblem(N_Vector y, UserData* udata, UserOptions* uopt);


/*
 * Helper functions
 */

/* function to set initial condition */
int SetIC(N_Vector y, UserData* udata);

/* function to fill neighbor data */
int FillSendBuffers(N_Vector y, UserData* udata);

/* functions for processing command line args */
int SetupProblem(int argc, char *argv[], UserData* udata, UserOptions* uopt,
                 SUNContext ctx);
void InputError(char *name);
int ComponentMask(N_Vector mask, const int component, const UserData* udata);

/* function to write solution to disk */
int WriteOutput(realtype t, N_Vector y, UserData* udata, UserOptions* uopt);

#endif
