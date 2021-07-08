/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *                David Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a steady-state anisotropic 2D heat equation
 * with an additional nonlinear term,
 *
 *   b = kx u_xx + ky u_yy + c(u),
 *
 * for (x,y) in [0, 1]^2, with boundary conditions
 *
 *   u = 0,
 *
 * and the forcing term
 *
 *   b(x,y) = kx 2 pi^2 * (cos^2(pi x) - sin^2(pi x)) sin^2(pi y)
 *            ky 2 pi^2 * (cos^2(pi y) - sin^2(pi y)) sin^2(pi x)
 *            + c(u_exact).
 *
 * Under this setup, the problem has the analytical solution
 *
 *   u(x,y) = u_exact = sin^2(pi x ) sin^2(pi y).
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid.
 * The problem is solved via Fixed Point Iteration by adding u to both
 * sides of the equation to setup the fixed point function. 
 * This matrix-free setup has the following form
 * 
 *   u = c3 * u_{ij} + c1*(u_{i-1,j} + u_{i+1,j}) + c2*(u_{i,j-1} + u_{i,j+1})
                 + c(u)_{ij} - b_{ij} + u_{ij}  
 *
 * where the constants c1, c2, and c3 are defined as follows
 *
 *   c1 = kx/hx^2 
 *   c2 = ky/hy^2 
 *   c3  = -2 * (c1 + c2)
 *
 * Several command line options are available to change the problem parameters
 * and KINSOL settings. Use the flag --help for more information. 
 * ---------------------------------------------------------------------------*/

#include "kin_heat2D_nonlin_p.hpp"     // header file containing UserData
                                       //   and function declartions

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int flag;                   // reusable error-checking flag
  UserData *udata    = NULL;  // user data structure
  N_Vector u         = NULL;  // vector for storing solution
  N_Vector scale     = NULL;  // vector for scaling initial value
  void *kin_mem      = NULL;  // KINSOL memory structure

  // Timing variables
  double t1 = 0.0;
  double t2 = 0.0;

  // MPI variables
  MPI_Comm comm_w = MPI_COMM_WORLD; // MPI communicator
  int myid;                         // MPI process ID

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  flag = MPI_Comm_rank(comm_w, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

  // Set output process flag
  bool outproc = (myid == 0);

  // ------------------------------------------
  // Setup UserData and parallel decomposition
  // ------------------------------------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) return 1;

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata, outproc);
  if (flag != 0) return 1;

  // Setup parallel decomposition
  flag = SetupDecomp(comm_w, udata);
  if (check_flag(&flag, "SetupDecomp", 1)) return 1;

  // Output problem setup/options
  if (outproc)
  {
    flag = PrintUserData(udata);
    if (check_flag(&flag, "PrintUserData", 1)) return 1;
  }

  // ------------------------
  // Create parallel vectors
  // ------------------------

  // Create vector for solution
  u = N_VNew_Parallel(udata->comm_c, udata->nodes_loc, udata->nodes);
  if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;

  // Create vector for scaling initial value
  scale = N_VNew_Parallel(udata->comm_c, udata->nodes_loc, udata->nodes);
  if (check_flag((void *) scale, "N_VNew_Parallel", 0)) return 1;
  N_VConst(ONE, scale);

  // Set initial condition
  N_VConst(ONE, u);
 
  // Create vector for error
  udata->e = N_VClone(u);
  if (check_flag((void *) (udata->e), "N_VClone", 0)) return 1;

  // Create vector for b
  udata->b = N_VClone(u);
  if (check_flag((void *) (udata->b), "N_VClone", 0)) return 1;

  // Create temp vector for FPFunction evaluation
  udata->vtemp = N_VClone(u);
  if (check_flag((void *) (udata->vtemp), "N_VClone", 0)) return 1;
  
  // -----------------------------------------------
  // Set b and c(u) for RHS evaluation in FPFunction
  // -----------------------------------------------
  flag = SetupRHS(udata);
  if (check_flag(&flag, "SetupRHS", 1)) return 1;

  // --------------
  // Setup KINSOL 
  // --------------

  // Initialize KINSOL memory
  kin_mem = KINCreate();
  if (check_flag((void *) kin_mem, "KINCreate", 0)) return 1;

  // Set number of prior residuals used in Anderson Accleration
  flag = KINSetMAA(kin_mem, udata->maa);
  if (check_flag(&flag, "KINSetMAA", 0)) return 1;
  
  // Set orthogonlization routine used in Anderson Accleration
  flag = KINSetOrthAA(kin_mem, udata->orthaa);
  if (check_flag(&flag, "KINSetOrthAA", 0)) return 1;

  // Set Fixed Point Function
  flag = KINInit(kin_mem, FPFunction, u);
  if (check_flag(&flag, "KINInit", 1)) return 1;

  // Specify tolerances
  flag = KINSetFuncNormTol(kin_mem, udata->rtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return 1;

  // Set maximum number of iterations
  flag = KINSetNumMaxIters(kin_mem, udata->maxits);
  if (check_flag(&flag, "KINSetMaxNumIters", 1)) return 1;

  // Set Anderson Acceleration damping parameter
  flag = KINSetDampingAA(kin_mem, udata->damping);
  if (check_flag(&flag, "KINSetDampingAA", 1)) return 1;

  // Attach user data
  flag = KINSetUserData(kin_mem, (void *) udata);
  if (check_flag(&flag, "KINSetUserData", 1)) return 1;

  // ----------------------------
  // Call KINSol to solve problem
  // ----------------------------

  // No scaling used
  N_VConst(ONE, scale);
  
  // Start timer
  t1 = MPI_Wtime();

  // Call main solver
  flag = KINSol(kin_mem,        // KINSol memory block
                u,              // inital guess on input; solution vector
                KIN_FP,         // global strategy choice
                scale,          // scaling vector, for the variable u 
                scale);         // scaling vector for function values fval
  if (check_flag(&flag, "KINSol", 1)) return(1);

  // Stop timer
  t2 = MPI_Wtime();

  // Update timer
  udata->totaltime = t2 - t1;

  // -----------------------
  // Get solver statistics
  // -----------------------
  if (udata->output > 0 && outproc)
  {
    cout << "Final statistics:" << endl;
    flag = OutputStats(kin_mem, udata);
    if (check_flag(&flag, "OutputStats", 1)) return 1;
  }

  // ------------------------
  // Calculate solution error
  // ------------------------

  // Output final error
  flag = SolutionError(u, udata->e, udata);
  if (check_flag(&flag, "SolutionError", 1)) return 1;

  realtype maxerr = N_VMaxNorm(udata->e);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << "  Max error = " << maxerr << endl;
    cout << endl;
  }

  // ------------------------------
  // Print timing
  // ------------------------------
  if (udata->timing)
  {
    flag = OutputTiming(udata);
    if (check_flag(&flag, "OutputTiming", 1)) return 1;
  }
  
  // ------------------------------
  // Free memory 
  // ------------------------------

  KINFree(&kin_mem);         // Free solver memory
  N_VDestroy(u);             // Free vectors
  N_VDestroy(scale);
  FreeUserData(udata);       // Free user data
  delete udata;
  flag = MPI_Finalize();     // Finalize MPI
  return 0;

}

// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData *udata)
{
  int flag;

  // Check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL ||
      udata->Srecv != NULL || udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  flag = MPI_Comm_size(comm_w, &(udata->nprocs_w));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Check the processor grid
  if ((udata->npx * udata->npy) != udata->nprocs_w)
  {
    cerr << "Error: npx * npy != nproc" << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = udata->npx;
  dims[1] = udata->npy;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  flag = MPI_Cart_create(comm_w, 2, dims, periods, 0, &(udata->comm_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(udata->comm_c, &(udata->myid_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(udata->comm_c, 2, dims, periods, coords);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << flag << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = udata->nx / dims[0];
  sunindextype rx = udata->nx % dims[0];

  udata->is = qx * idx + (idx < rx ? idx : rx);
  udata->ie = udata->is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (udata->ie > (udata->nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = udata->ny / dims[1];
  sunindextype ry = udata->ny % dims[1];

  udata->js = qy * idy + (idy < ry ? idy : ry);
  udata->je = udata->js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (udata->je > (udata->ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  udata->nx_loc = (udata->ie) - (udata->is) + 1;
  udata->ny_loc = (udata->je) - (udata->js) + 1;

  // Initialize global and local vector lengths
  udata->nodes     = udata->nx * udata->ny;
  udata->nodes_loc = udata->nx_loc * udata->ny_loc;

  // Determine if this proc has neighbors
  udata->HaveNbrW = (udata->is != 0);
  udata->HaveNbrE = (udata->ie != udata->nx-1);
  udata->HaveNbrS = (udata->js != 0);
  udata->HaveNbrN = (udata->je != udata->ny-1);

  // Allocate exchange buffers if necessary
  if (udata->HaveNbrW)
  {
    udata->Wrecv = new realtype[udata->ny_loc];
    udata->Wsend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrE)
  {
    udata->Erecv = new realtype[udata->ny_loc];
    udata->Esend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrS)
  {
    udata->Srecv = new realtype[udata->nx_loc];
    udata->Ssend = new realtype[udata->nx_loc];
  }
  if (udata->HaveNbrN)
  {
    udata->Nrecv = new realtype[udata->nx_loc];
    udata->Nsend = new realtype[udata->nx_loc];
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (udata->HaveNbrW)
  {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // East neighbor
  if (udata->HaveNbrE)
  {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // South neighbor
  if (udata->HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // North neighbor
  if (udata->HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    flag = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the solver 
// -----------------------------------------------------------------------------

// Fixed point function to compute G(u) 
static int FPFunction(N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) return -1;

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Constants for computing diffusion term
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain and compute rhs forcing term
  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;
  realtype cos_sqr_x, cos_sqr_y;

  // Iterate over subdomain interior and add rhs diffusion term
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] +
        cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) return -1;

  // Iterate over subdomain boundaries and add rhs diffusion term
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face (updates south-west and north-west corners if necessary)
  if (udata->HaveNbrW)
  {
    i = 0;
    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // East face (updates south-east and north-east corners if necessary)
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;
    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // South face (excludes corners)
  if (udata->HaveNbrS)
  {
    j = 0;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // North face (excludes corners)
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)] + 
        cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }
 
  // Add c(u)
  flag = udata->c(u, udata->vtemp, user_data);
  if (check_flag(&flag, "c(u)", 1)) return 1;
  N_VLinearSum(ONE, udata->vtemp, ONE, f, f);

  // Add u
  N_VLinearSum(ONE, u, ONE, f, f);

  // Subtract b 
  N_VLinearSum(-ONE, udata->b, ONE, f, f);

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->fevaltime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// RHS helper functions
// -----------------------------------------------------------------------------

// Set nonlinear function c(u)
static int SetC(UserData *udata)
{
  if (0 < udata->c_int && udata->c_int < 18)
  {
    if (udata->c_int == 1) udata->c = c1;
    else if (udata->c_int == 2) udata->c = c2;
    else if (udata->c_int == 3) udata->c = c3;
    else if (udata->c_int == 4) udata->c = c4;
    else if (udata->c_int == 5) udata->c = c5;
    else if (udata->c_int == 6) udata->c = c6;
    else if (udata->c_int == 7) udata->c = c7;
    else if (udata->c_int == 8) udata->c = c8;
    else if (udata->c_int == 9) udata->c = c9;
    else if (udata->c_int == 10) udata->c = c10;
    else if (udata->c_int == 11) udata->c = c11;
    else if (udata->c_int == 12) udata->c = c12;
    else if (udata->c_int == 13) udata->c = c13;
    else if (udata->c_int == 14) udata->c = c14;
    else if (udata->c_int == 15) udata->c = c15;
    else if (udata->c_int == 16) udata->c = c16;
    else if (udata->c_int == 17) udata->c = c17;
  }
  else return 1;

  // Return success
  return 0;
}

// Set nonlinear function c(u)
static int SetupRHS(void *user_data)
{
  int           flag;
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      :          1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      :          1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // ------------------------------------------------------
  // Setup function c(u) for FPFunction based on user input
  // ------------------------------------------------------
  flag = SetC(udata);
  if (check_flag(&flag, "SetC", 1)) return 1;
  
  // ------------------------------------------------------
  // Setup b for FPFunction 
  // ------------------------------------------------------

  // Access data array
  realtype *barray = N_VGetArrayPointer(udata->b);
  if (check_flag((void *) barray, "N_VGetArrayPointer", 0)) return 1;

  // Initialize rhs vector to zero (handles boundary conditions)
  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;
  realtype cos_sqr_x, cos_sqr_y;

  realtype bx = (udata->kx) * TWO * PI * PI;
  realtype by = (udata->ky) * TWO * PI * PI;

  // Iterate over subdomain and compute forcing term (b)
  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      x = (udata->is + i) * udata->dx;
      y = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      cos_sqr_x = cos(PI * x) * cos(PI * x);
      cos_sqr_y = cos(PI * y) * cos(PI * y);

      barray[IDX(i,j,nx_loc)] = 
        bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y
        + by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x;
    }
  }

  // Calculate c(u_exact) and add to forcing term (b)
  flag = Solution(udata->e, udata);
  if (check_flag(&flag, "rhs", 1)) return 1;

  flag = udata->c(udata->e, udata->vtemp, user_data);
  if (check_flag(&flag, "c(u)", 1)) return 1;

  // b = kx u_xx (u_exact) + ky u_yy (u_exact) + c(u_exact)
  N_VLinearSum(ONE, udata->vtemp, ONE, udata->b, udata->b);

  // Return success
  return 0;
}

// Post exchange receives
static int PostRecv(UserData *udata)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Open Irecv buffers
  if (udata->HaveNbrW)
  {
    flag = MPI_Irecv(udata->Wrecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, MPI_ANY_TAG, udata->comm_c, &(udata->reqRW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Irecv(udata->Erecv, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, MPI_ANY_TAG, udata->comm_c, &(udata->reqRE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Irecv(udata->Srecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, MPI_ANY_TAG, udata->comm_c, &(udata->reqRS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Irecv(udata->Nrecv, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, MPI_ANY_TAG, udata->comm_c, &(udata->reqRN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Send exchange data
static int SendData(N_Vector y, UserData *udata)
{
  int flag, i;
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Start timer
  double t1 = MPI_Wtime();

  // Access data array
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  // Send data
  if (udata->HaveNbrW)
  {
    for (i = 0; i < ny_loc; i++) udata->Wsend[i] = Y[IDX(0,i,nx_loc)];
    flag = MPI_Isend(udata->Wsend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, 0, udata->comm_c, &(udata->reqSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    for (i = 0; i < ny_loc; i++) udata->Esend[i] = Y[IDX(nx_loc-1,i,nx_loc)];
    flag = MPI_Isend(udata->Esend, (int) udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, 1, udata->comm_c, &(udata->reqSE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    for (i = 0; i < nx_loc; i++) udata->Ssend[i] = Y[IDX(i,0,nx_loc)];
    flag = MPI_Isend(udata->Ssend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, 2, udata->comm_c, &(udata->reqSS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    for (i = 0; i < nx_loc; i++) udata->Nsend[i] = Y[IDX(i,ny_loc-1,nx_loc)];
    flag = MPI_Isend(udata->Nsend, (int) udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, 3, udata->comm_c, &(udata->reqSN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Wait for exchange data
static int WaitRecv(UserData *udata)
{
  // Local variables
  int flag;
  MPI_Status stat;

  // Start timer
  double t1 = MPI_Wtime();

  // Wait for messages to finish
  if (udata->HaveNbrW)
  {
    flag = MPI_Wait(&(udata->reqRW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Wait(&(udata->reqRE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Wait(&(udata->reqRS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Wait(&(udata->reqRN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  // Diffusion coefficient
  udata->kx = ONE;
  udata->ky = ONE;

    // Upper bounds in x and y directions
  udata->xu = ONE;
  udata->yu = ONE;

  // Global number of nodes in the x and y directions
  udata->nx    = 35;
  udata->ny    = 35;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Locals number of nodes in the x and y directions (set in SetupDecomp)
  udata->nx_loc    = 0;
  udata->ny_loc    = 0;
  udata->nodes_loc = 0;

  // Global indices of this subdomain (set in SetupDecomp)
  udata->is = 0;
  udata->ie = 0;
  udata->js = 0;
  udata->je = 0;

  // MPI variables (set in SetupDecomp)
  udata->comm_c = MPI_COMM_NULL;

  udata->nprocs_w = 1;
  udata->npx      = 1;
  udata->npy      = 1;

  udata->myid_c = 0;

  // Flags denoting neighbors (set in SetupDecomp)
  udata->HaveNbrW = true;
  udata->HaveNbrE = true;
  udata->HaveNbrS = true;
  udata->HaveNbrN = true;

  // Exchange receive buffers (allocated in SetupDecomp)
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;

  // Exchange send buffers (allocated in SetupDecomp)
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  // Neighbors IDs (set in SetupDecomp)
  udata->ipW = -1;
  udata->ipE = -1;
  udata->ipS = -1;
  udata->ipN = -1;

  // Integrator settings
  udata->rtol        = RCONST(1.e-8);   // relative tolerance
  udata->maa         = 60;              // 60 vectors in Anderson Acceleration space 
  udata->damping     = ONE;             // no damping for Anderson Acceleration
  udata->orthaa      = 0;               // use MGS for Anderson Acceleration 
  udata->maxits      = 200;             // max number of fixed point iterations

  // c function
  udata->c     = NULL;
  udata->c_int = 1; 

  // Vectors
  udata->b     = NULL;
  udata->vtemp = NULL;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->e      = NULL;

  // Timing variables
  udata->timing       = true;
  udata->totaltime    = 0.0;
  udata->fevaltime    = 0.0;
  udata->exchangetime = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // Free exchange buffers
  if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
  if (udata->Wsend != NULL)  delete[] udata->Wsend;
  if (udata->Erecv != NULL)  delete[] udata->Erecv;
  if (udata->Esend != NULL)  delete[] udata->Esend;
  if (udata->Srecv != NULL)  delete[] udata->Srecv;
  if (udata->Ssend != NULL)  delete[] udata->Ssend;
  if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
  if (udata->Nsend != NULL)  delete[] udata->Nsend;

  // Free MPI Cartesian communicator
  if (udata->comm_c != MPI_COMM_NULL)
    MPI_Comm_free(&(udata->comm_c));

  // Free b vector
  if (udata->b)
  {
    N_VDestroy(udata->b);
    udata->b = NULL;
  }
  
  // Free temporary vector
  if (udata->vtemp)
  {
    N_VDestroy(udata->vtemp);
    udata->vtemp = NULL;
  }

  // Free error vector
  if (udata->e)
  {
    N_VDestroy(udata->e);
    udata->e = NULL;
  }

  // Return success
  return 0;
}

// Read command line inputs
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc)
{
  // Check for input args
  int arg_idx = 1;

  while (arg_idx < (*argc))
  {
    string arg = (*argv)[arg_idx++];

    // Mesh points
    if (arg == "--mesh")
    {
      udata->nx = stoi((*argv)[arg_idx++]);
      udata->ny = stoi((*argv)[arg_idx++]);
    }
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
    }
    // Domain upper bounds
    else if (arg == "--domain")
    {
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--k")
    {
      udata->kx = stod((*argv)[arg_idx++]);
      udata->ky = stod((*argv)[arg_idx++]);
    }
    // Fixed Point settings
    else if (arg == "--rtol")
    {
      udata->rtol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--maa")
    {
      udata->maa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--damping")
    {
      udata->damping = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--orthaa")
    {
      udata->orthaa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--maxits")
    {
      udata->maxits = stoi((*argv)[arg_idx++]);
    }
    // RHS settings
    else if (arg == "--c")
    {
      udata->c_int = stoi((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output")
    {
      udata->output = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--timing")
    {
      udata->timing = true;
    }
    // Help
    else if (arg == "--help")
    {
      if (outproc) InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      if (outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
      return -1;
    }
  }

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(N_Vector u, UserData *udata)
{
  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;

  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y;
    }
  }

  return 0;
}

// Compute the solution error
static int SolutionError(N_Vector u, N_Vector e, UserData *udata)
{
  // Check absolute error between output u and G(u)
  int flag = Solution(e, udata);
  if (flag != 0) return -1;
  
  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --mesh <nx> <ny>        : mesh points in the x and y directions" << endl;
  cout << "  --np <npx> <npy>        : number of MPI processes in the x and y directions" << endl;
  cout << "  --domain <xu> <yu>      : domain upper bound in the x and y direction" << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --maa                   : size of Anderson Acceleration subspace" << endl;
  cout << "  --damping               : damping for Anderson Acceleration" << endl;
  cout << "  --orthaa                : orthogonalization routined used in Anderson Acceleration" << endl;
  cout << "  --maxits <iterations>   : max fixed point iterations" << endl;
  cout << "  --c <c_int>             : nonlinear function choice" << endl;
  cout << "  --output                : output nonlinear solver statistics" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "Stationary 2D Heat PDE with Nonlinear Terms:"  << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  nprocs         = " << udata->nprocs_w        << endl;
  cout << "  npx            = " << udata->npx             << endl;
  cout << "  npy            = " << udata->npy             << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  kx             = " << udata->kx              << endl;
  cout << "  ky             = " << udata->ky              << endl;
  cout << "  xu             = " << udata->xu              << endl;
  cout << "  yu             = " << udata->yu              << endl;
  cout << "  nx             = " << udata->nx              << endl;
  cout << "  ny             = " << udata->ny              << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc          << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc          << endl;
  cout << "  dx             = " << udata->dx              << endl;
  cout << "  dy             = " << udata->dy              << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  rtol           = " << udata->rtol            << endl;
  cout << "  maa            = " << udata->maa             << endl;
  cout << "  damping        = " << udata->damping         << endl;
  cout << "  orthaa         = " << udata->orthaa          << endl;
  cout << "  maxits         = " << udata->maxits          << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  c              = " << udata->c_int           << endl;
  cout << " --------------------------------- "           << endl;
  cout << "  output         = " << udata->output          << endl;
  cout << " --------------------------------- "           << endl;
  cout << endl;

  return 0;
}

// Print nonlinear solver statistics
static int OutputStats(void *kinsol_mem, UserData* udata)
{
  int flag;

  // Get solver stats
  long int nfe, nni;
  flag = KINGetNumNonlinSolvIters(kinsol_mem, &nni);
  if (check_flag(&flag, "KINGetNumNonlinSolvIters", 1)) return(1);
  flag = KINGetNumFuncEvals(kinsol_mem, &nfe);
  if (check_flag(&flag, "KINGetNumFuncEvals", 1)) return(1);

  cout << setprecision(6);

  cout << "  Func evals       = " << nfe     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << endl;

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->totaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Total time                = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->fevaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Function evaluation time  = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->exchangetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Boundary exchange time    = " << maxtime << " sec" << endl;
    cout << endl;
  }

  return 0;
}

// Check function return value
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return flag value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int *) flagvalue);
    if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl << "ERROR: " << funcname << " returned with flag = "
           << errflag << endl << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl << "ERROR: check_flag called with an invalid option value"
         << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
