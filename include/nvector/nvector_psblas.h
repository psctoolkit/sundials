/* -----------------------------------------------------------------
 * Programmer(s): Fabio Durastante @ IAC-CNR
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the main header file for the PSBLAS-enabled implementation
 * of the NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type sunbooleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_PSBLAS(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PSBLAS
#define _NVECTOR_PSBLAS

#include <stdio.h>
#include <mpi.h>
#include "psb_base_cbind.h"
#include "psb_c_dbase.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_mpi_types.h>

#undef I

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PSBLAS implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_PSBLAS {
  sunbooleantype own_data;        /* ownership of data                */
  psb_c_descriptor *cdh;       /* descriptor for data distribution */
  psb_c_dvector *pvec;	       /* PSBLAS vector                    */
  psb_c_ctxt *cctxt;           /* PSBLAS communicator              */
};

typedef struct _N_VectorContent_PSBLAS *N_VectorContent_PSBLAS;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_P, NV_DESCRIPTOR_P, NV_OWN_DATA_P,
 *    NV_PVEC_P, NV_CCTXT_P
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_P(v)    ( (N_VectorContent_PSBLAS)(v->content) )

#define NV_DESCRIPTOR_P(v)  ( NV_CONTENT_P(v)->cdh )

#define NV_OWN_DATA_P(v)   ( NV_CONTENT_P(v)->own_data )

#define NV_PVEC_P(v)       ( NV_CONTENT_P(v)->pvec )

#define NV_CCTXT_P(v)      ( NV_CONTENT_P(v)->cctxt )


/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_PSBLAS
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_PSBLAS(psb_c_ctxt *cctxt, psb_c_descriptor *cdh);

SUNDIALS_EXPORT N_Vector N_VNewEmpty_PSBLAS(psb_c_ctxt *cctxt, psb_c_descriptor *cdh);

SUNDIALS_EXPORT N_Vector N_VMake_PSBLAS(psb_c_ctxt *cctxt, psb_c_descriptor *cdh,
  psb_i_t m, psb_l_t *irow,double *val);

SUNDIALS_EXPORT void N_VAsb_PSBLAS(N_Vector v);

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_PSBLAS(int count, N_Vector w);

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_PSBLAS(int count, N_Vector w);

SUNDIALS_EXPORT void N_VDestroyVectorArray_PSBLAS(N_Vector *vs, int count);

SUNDIALS_EXPORT sunindextype N_VGetLength_PSBLAS(N_Vector v);

SUNDIALS_EXPORT sunindextype N_VGetLocalLength_PSBLAS(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_PSBLAS(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_PSBLAS(N_Vector v, FILE *outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_PSBLAS(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_PSBLAS(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_PSBLAS(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_PSBLAS(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_PSBLAS(N_Vector v, sunindextype *lrw,
                                       sunindextype *liw);
SUNDIALS_EXPORT sunrealtype *N_VGetArrayPointer_PSBLAS(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_PSBLAS(sunrealtype *v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_PSBLAS(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_PSBLAS(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_PSBLAS(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_PSBLAS(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_PSBLAS(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_PSBLAS(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_PSBLAS(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_PSBLAS(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_PSBLAS(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_PSBLAS(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_PSBLAS(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_PSBLAS(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_PSBLAS(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_PSBLAS(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_PSBLAS(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_PSBLAS(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_PSBLAS(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_PSBLAS(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_PSBLAS(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_PSBLAS(int nvec, sunrealtype* c, N_Vector* V,
                                                  N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_PSBLAS(int nvec, sunrealtype* a, N_Vector x,
                                              N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_PSBLAS(int nvec, N_Vector x,
                                             N_Vector *Y, sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_PSBLAS(int nvec,
                                                     sunrealtype a, N_Vector* X,
                                                     sunrealtype b, N_Vector* Y,
                                                     N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_PSBLAS(int nvec, sunrealtype* c,
                                                 N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_PSBLAS(int nvecs, sunrealtype c,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_PSBLAS(int nvecs, N_Vector* X,
                                                    N_Vector* W, sunrealtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_PSBLAS(int nvec, N_Vector* X,
                                                        N_Vector* W, N_Vector id,
                                                        sunrealtype* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_PSBLAS(int nvec, int nsum,
                                                         sunrealtype* a,
                                                         N_Vector* X,
                                                         N_Vector** Y,
                                                         N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_PSBLAS(int nvec, int nsum,
                                                             sunrealtype* c,
                                                             N_Vector** X,
                                                             N_Vector* Z);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_PSBLAS(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_PSBLAS(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_PSBLAS(N_Vector v, sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
