/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.0
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

/* ---------------------------------------------------------------
 * Programmer(s): Auto-generated by swig.
 * ---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -------------------------------------------------------------*/

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* qualifier for exported *const* global data variables*/
#ifndef SWIGEXTERN
# ifdef __cplusplus
#   define SWIGEXTERN extern
# else
#   define SWIGEXTERN
# endif
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13




#include <assert.h>
#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
 { printf("In " DECL ": " MSG); assert(0); RETURNNULL; }


#include <stdio.h>
#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(_WATCOM)
# ifndef snprintf
#  define snprintf _snprintf
# endif
#endif


/* Support for the `contract` feature.
 *
 * Note that RETURNNULL is first because it's inserted via a 'Replaceall' in
 * the fortran.cxx file.
 */
#define SWIG_contract_assert(RETURNNULL, EXPR, MSG) \
 if (!(EXPR)) { SWIG_exception_impl("$decl", SWIG_ValueError, MSG, RETURNNULL); } 


#define SWIGVERSION 0x040000 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) (void *)((const void *)(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),(void**)(a)) 


#include "arkode/arkode_lsrkstep.h"

SWIGEXPORT void * _wrap_FLSRKStepCreate(ARKRhsFn farg1, ARKRhsFn farg2, double const *farg3, N_Vector farg4, void *farg5) {
  void * fresult ;
  ARKRhsFn arg1 = (ARKRhsFn) 0 ;
  ARKRhsFn arg2 = (ARKRhsFn) 0 ;
  sunrealtype arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  SUNContext arg5 = (SUNContext) 0 ;
  void *result = 0 ;
  
  arg1 = (ARKRhsFn)(farg1);
  arg2 = (ARKRhsFn)(farg2);
  arg3 = (sunrealtype)(*farg3);
  arg4 = (N_Vector)(farg4);
  arg5 = (SUNContext)(farg5);
  result = (void *)LSRKStepCreate(arg1,arg2,arg3,arg4,arg5);
  fresult = result;
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepReInit(void *farg1, ARKRhsFn farg2, ARKRhsFn farg3, double const *farg4, N_Vector farg5) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKRhsFn arg2 = (ARKRhsFn) 0 ;
  ARKRhsFn arg3 = (ARKRhsFn) 0 ;
  sunrealtype arg4 ;
  N_Vector arg5 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKRhsFn)(farg2);
  arg3 = (ARKRhsFn)(farg3);
  arg4 = (sunrealtype)(*farg4);
  arg5 = (N_Vector)(farg5);
  result = (int)LSRKStepReInit(arg1,arg2,arg3,arg4,arg5);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetMethod(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKODE_LSRKMethodType arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKODE_LSRKMethodType)(*farg2);
  result = (int)LSRKStepSetMethod(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetDomEigFn(void *farg1, ARKDomEigFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKDomEigFn arg2 = (ARKDomEigFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKDomEigFn)(farg2);
  result = (int)LSRKStepSetDomEigFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetDomEigFrequency(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)LSRKStepSetDomEigFrequency(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetMaxStageNum(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)LSRKStepSetMaxStageNum(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetDomEigSafetyFactor(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  sunrealtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (sunrealtype)(*farg2);
  result = (int)LSRKStepSetDomEigSafetyFactor(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepSetSSPStageNum(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)LSRKStepSetSSPStageNum(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepGetNumRhsEvals(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)LSRKStepGetNumRhsEvals(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepGetNumDomEigUpdates(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)LSRKStepGetNumDomEigUpdates(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepGetMaxStageNum(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)LSRKStepGetMaxStageNum(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepGetAverStageNum(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)LSRKStepGetAverStageNum(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FLSRKStepGetTimestepperStats(void *farg1, long *farg2, long *farg3, long *farg4, long *farg5, long *farg6, long *farg7, long *farg8, long *farg9, double *farg10, double *farg11) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  long *arg4 = (long *) 0 ;
  long *arg5 = (long *) 0 ;
  long *arg6 = (long *) 0 ;
  long *arg7 = (long *) 0 ;
  long *arg8 = (long *) 0 ;
  long *arg9 = (long *) 0 ;
  sunrealtype *arg10 = (sunrealtype *) 0 ;
  sunrealtype *arg11 = (sunrealtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  arg4 = (long *)(farg4);
  arg5 = (long *)(farg5);
  arg6 = (long *)(farg6);
  arg7 = (long *)(farg7);
  arg8 = (long *)(farg8);
  arg9 = (long *)(farg9);
  arg10 = (sunrealtype *)(farg10);
  arg11 = (sunrealtype *)(farg11);
  result = (int)LSRKStepGetTimestepperStats(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
  fresult = (int)(result);
  return fresult;
}


