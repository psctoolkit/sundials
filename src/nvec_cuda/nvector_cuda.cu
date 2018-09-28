/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a MPI+CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <nvector/cuda/Vector.hpp>
#include <nvector/cuda/VectorKernels.cuh>
#include <nvector/cuda/VectorArrayKernels.cuh>
#include <sundials/sundials_mpi.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

using namespace suncudavec;

/*
 * Type definitions
 */

typedef suncudavec::Vector<realtype, sunindextype> vector_type;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Cuda(N_Vector v)
{
  return SUNDIALS_NVEC_CUDA;
}

N_Vector N_VNewEmpty_Cuda(sunindextype length)
{
  N_Vector v;
  N_Vector_Ops ops;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Cuda;
  ops->nvclone           = N_VClone_Cuda;
  ops->nvcloneempty      = N_VCloneEmpty_Cuda;
  ops->nvdestroy         = N_VDestroy_Cuda;
  ops->nvspace           = N_VSpace_Cuda;
  ops->nvgetarraypointer = NULL;
  ops->nvsetarraypointer = NULL;

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_Cuda;
  ops->nvconst        = N_VConst_Cuda;
  ops->nvprod         = N_VProd_Cuda;
  ops->nvdiv          = N_VDiv_Cuda;
  ops->nvscale        = N_VScale_Cuda;
  ops->nvabs          = N_VAbs_Cuda;
  ops->nvinv          = N_VInv_Cuda;
  ops->nvaddconst     = N_VAddConst_Cuda;
  ops->nvdotprod      = N_VDotProd_Cuda;
  ops->nvmaxnorm      = N_VMaxNorm_Cuda;
  ops->nvwrmsnormmask = N_VWrmsNormMask_Cuda;
  ops->nvwrmsnorm     = N_VWrmsNorm_Cuda;
  ops->nvmin          = N_VMin_Cuda;
  ops->nvwl2norm      = N_VWL2Norm_Cuda;
  ops->nvl1norm       = N_VL1Norm_Cuda;
  ops->nvcompare      = N_VCompare_Cuda;
  ops->nvinvtest      = N_VInvTest_Cuda;
  ops->nvconstrmask   = N_VConstrMask_Cuda;
  ops->nvminquotient  = N_VMinQuotient_Cuda;

  /* fused vector operations */
  ops->nvlinearcombination = N_VLinearCombination_Cuda;
  ops->nvscaleaddmulti     = N_VScaleAddMulti_Cuda;
  ops->nvdotprodmulti      = NULL; // N_VDotProdMulti_Cuda;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Cuda;
  ops->nvscalevectorarray             = N_VScaleVectorArray_Cuda;
  ops->nvconstvectorarray             = N_VConstVectorArray_Cuda;
  ops->nvwrmsnormvectorarray          = NULL; // N_VWrmsNormVectorArray_Cuda;
  ops->nvwrmsnormmaskvectorarray      = NULL; // N_VWrmsNormMaskVectorArray_Cuda;
  ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Cuda;
  ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;

  /* Attach ops and set content to NULL */
  v->content = NULL;
  v->ops     = ops;

  return(v);
}

#if SUNDIALS_MPI_ENABLED
N_Vector N_VNew_Cuda(MPI_Comm comm,
                     sunindextype local_length,
                     sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(local_length);
  if (v == NULL)
    return(NULL);

  v->content = new vector_type(comm, local_length, global_length);

  return(v);
}
#else
N_Vector N_VNew_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(length);
  if (v == NULL)
    return(NULL);

  v->content = new vector_type(SUNMPI_COMM_WORLD, length, length);

  return(v);
}
#endif


N_Vector N_VMake_Cuda(N_VectorContent_Cuda c)
{
  N_Vector v;
  vector_type* x = static_cast<vector_type*>(c);
  sunindextype length = x->size();

  v = NULL;
  v = N_VNewEmpty_Cuda(length);
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}

/* -----------------------------------------------------------------
 * Function to return the length of the vector.
 */
sunindextype N_VGetLength_Cuda(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return xd->size();
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyFromDev();
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to stdout
 */

void N_VPrint_Cuda(N_Vector x)
{
  N_VPrintFile_Cuda(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to outfile
 */

void N_VPrintFile_Cuda(N_Vector x, FILE *outfile)
{
  sunindextype i;
  vector_type* xd = static_cast<vector_type*>(x->content);

  for (i = 0; i < xd->size(); i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd->host()[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd->host()[i]);
#else
    fprintf(outfile, "%11.8g\n", xd->host()[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}


/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Cuda(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;

  /* standard vector operations */
  ops->nvlinearsum    = w->ops->nvlinearsum;
  ops->nvconst        = w->ops->nvconst;
  ops->nvprod         = w->ops->nvprod;
  ops->nvdiv          = w->ops->nvdiv;
  ops->nvscale        = w->ops->nvscale;
  ops->nvabs          = w->ops->nvabs;
  ops->nvinv          = w->ops->nvinv;
  ops->nvaddconst     = w->ops->nvaddconst;
  ops->nvdotprod      = w->ops->nvdotprod;
  ops->nvmaxnorm      = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
  ops->nvmin          = w->ops->nvmin;
  ops->nvwl2norm      = w->ops->nvwl2norm;
  ops->nvl1norm       = w->ops->nvl1norm;
  ops->nvcompare      = w->ops->nvcompare;
  ops->nvinvtest      = w->ops->nvinvtest;
  ops->nvconstrmask   = w->ops->nvconstrmask;
  ops->nvminquotient  = w->ops->nvminquotient;

  /* fused vector operations */
  ops->nvlinearcombination = w->ops->nvlinearcombination;
  ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /* Create content */
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Cuda(N_Vector w)
{
  N_Vector v;
  vector_type* wdat = static_cast<vector_type*>(w->content);
  vector_type* vdat = new vector_type(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Cuda(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Cuda(N_Vector v)
{
  vector_type* x = static_cast<vector_type*>(v->content);
  if (x != NULL) {
    delete x;
    v->content = NULL;
  }

  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Cuda(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  int npes;

  SUNMPI_Comm_size(comm, &npes);

  *lrw = getGlobalSize<realtype, sunindextype>(X);
  *liw = 2*npes;
}

void N_VConst_Cuda(realtype a, N_Vector X)
{
  vector_type *xvec = extract<realtype, sunindextype>(X);
  setConst(a, *xvec);
}

void N_VLinearSum_Cuda(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *yvec = extract<realtype, sunindextype>(Y);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  linearSum(a, *xvec, b, *yvec, *zvec);
}

void N_VProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *yvec = extract<realtype, sunindextype>(Y);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  prod(*xvec, *yvec, *zvec);
}

void N_VDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *yvec = extract<realtype, sunindextype>(Y);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  div(*xvec, *yvec, *zvec);
}

void N_VScale_Cuda(realtype a, N_Vector X, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  scale(a, *xvec, *zvec);
}

void N_VAbs_Cuda(N_Vector X, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  absVal(*xvec, *zvec);
}

void N_VInv_Cuda(N_Vector X, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  inv(*xvec, *zvec);
}

void N_VAddConst_Cuda(N_Vector X, realtype b, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  addConst(b, *xvec, *zvec);
}

realtype N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *yvec = extract<realtype, sunindextype>(Y);

  realtype sum = dotProd(*xvec, *yvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);

  realtype locmax = maxNorm(*xvec);

  realtype globmax = SUNMPI_Allreduce_scalar(locmax, 2, comm);
  return globmax;
}

realtype N_VWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const sunindextype Nglob = getGlobalSize<realtype,sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *wvec = extract<realtype, sunindextype>(W);

  realtype sum = wL2NormSquare(*xvec, *wvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum/Nglob);
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const sunindextype Nglob = getGlobalSize<realtype,sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *wvec = extract<realtype, sunindextype>(W);
  const vector_type *ivec = extract<realtype, sunindextype>(Id);

  realtype sum = wL2NormSquareMask(*xvec, *wvec, *ivec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum/Nglob);
}

realtype N_VMin_Cuda(N_Vector X)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);

  realtype locmin = findMin(*xvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return globmin;
}

realtype N_VWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  const vector_type *wvec = extract<realtype, sunindextype>(W);

  realtype sum = wL2NormSquare(*xvec, *wvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum);
}

realtype N_VL1Norm_Cuda(N_Vector X)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);

  realtype sum = L1Norm(*xvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  compare(c, *xvec, *zvec);
}

booleantype N_VInvTest_Cuda(N_Vector X, N_Vector Z)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *zvec = extract<realtype, sunindextype>(Z);
  realtype locmin = invTest(*xvec, *zvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return (globmin < HALF);
}

/*
 * Creates mask for variables violating constraints
 */
booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X);
  const vector_type *cvec = extract<realtype, sunindextype>(C);
  const vector_type *xvec = extract<realtype, sunindextype>(X);
  vector_type *mvec = extract<realtype, sunindextype>(M);

  realtype locsum = constrMask(*cvec, *xvec, *mvec);

  realtype globsum = SUNMPI_Allreduce_scalar(locsum, 1, comm);
  return (globsum < HALF);
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(num);
  const vector_type *numvec = extract<realtype, sunindextype>(num);
  const vector_type *denvec = extract<realtype, sunindextype>(denom);

  realtype locmin = minQuotient(*numvec, *denvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return globmin;
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type*  Zv;

  Zv = extract<realtype, sunindextype>(Z);

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = extract<realtype, sunindextype>(X[i]);

  err = linearCombination(nvec, c, Xv, Zv);

  delete[] Xv;

  return err == cudaSuccess ? 0 : -1;
}

int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X, N_Vector* Y,
                           N_Vector* Z)
{
  cudaError_t err;
  vector_type*  Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = extract<realtype, sunindextype>(X);

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = extract<realtype, sunindextype>(Y[i]);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = extract<realtype, sunindextype>(Z[i]);

  err = scaleAddMulti(nvec, c, Xv, Yv, Zv);

  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  cudaError_t err;
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(x);
  sunindextype N = getGlobalSize<realtype, sunindextype>(x);
  vector_type*  Xv;
  vector_type** Yv;

  Xv = extract<realtype, sunindextype>(x);

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = extract<realtype, sunindextype>(Y[i]);

  err = dotProdMulti(nvec, Xv, Yv, dotprods);

  delete[] Yv;

  SUNMPI_Allreduce(dotprods, nvec, 1, comm);

  return err == cudaSuccess ? 0 : -1;
}



/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Cuda(int nvec, realtype a, N_Vector* X, realtype b,
                                 N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = extract<realtype, sunindextype>(X[i]);

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = extract<realtype, sunindextype>(Y[i]);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = extract<realtype, sunindextype>(Z[i]);

  err = linearSumVectorArray(nvec, a, Xv, b, Yv, Zv);

  delete[] Xv;
  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = extract<realtype, sunindextype>(X[i]);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = extract<realtype, sunindextype>(Z[i]);

  err = scaleVectorArray(nvec, c, Xv, Zv);

  delete[] Xv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Zv;

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = extract<realtype, sunindextype>(Z[i]);

  err = constVectorArray(nvec, c, Zv);

  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  cudaError_t err;
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X[0]);
  sunindextype N = getGlobalSize<realtype, sunindextype>(X[0]);
  vector_type** Xv;
  vector_type** Wv;

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = extract<realtype, sunindextype>(X[k]);

  Wv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Wv[k] = extract<realtype, sunindextype>(W[k]);

  err = wL2NormSquareVectorArray(nvec, Xv, Wv, norms);

  delete[] Xv;
  delete[] Wv;

  SUNMPI_Allreduce(norms, nvec, 1, comm);

  for (int k=0; k<nvec; ++k) {
    norms[k] = std::sqrt(norms[k]/N);
  }

  return err == cudaSuccess ? 0 : -1;
}


int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  cudaError_t err;
  SUNMPI_Comm comm = getMPIComm<realtype, sunindextype>(X[0]);
  sunindextype N = getGlobalSize<realtype, sunindextype>(X[0]);
  vector_type** Xv;
  vector_type** Wv;
  vector_type*  IDv;

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = extract<realtype, sunindextype>(X[k]);

  Wv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Wv[k] = extract<realtype, sunindextype>(W[k]);

  IDv = extract<realtype, sunindextype>(id);

  err = wL2NormSquareMaskVectorArray(nvec, Xv, Wv, IDv, norms);

  delete[] Xv;
  delete[] Wv;

  SUNMPI_Allreduce(norms, nvec, 1, comm);

  for (int k=0; k<nvec; ++k) {
    norms[k] = std::sqrt(norms[k]/N);
  }

  return err == cudaSuccess ? 0 : -1;
}


int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                      N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = extract<realtype, sunindextype>(X[k]);

  Yv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Yv[k*nsum+j] = extract<realtype, sunindextype>(Y[j][k]);

  Zv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Zv[k*nsum+j] = extract<realtype, sunindextype>(Z[j][k]);

  err = scaleAddMultiVectorArray(nvec, nsum, c, Xv, Yv, Zv);

  delete[] Xv;
  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Zv;

  Xv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Xv[k*nsum+j] = extract<realtype, sunindextype>(X[j][k]);

  Zv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Zv[k] = extract<realtype, sunindextype>(Z[k]);

  err = linearCombinationVectorArray(nvec, nsum, c, Xv, Zv);

  delete[] Xv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}

} // extern "C"
