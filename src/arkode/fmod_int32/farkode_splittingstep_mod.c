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


enum {
    SWIG_MEM_OWN = 0x01,
    SWIG_MEM_RVALUE = 0x02,
    SWIG_MEM_CONST = 0x04
};


#define SWIG_check_mutable(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
    if ((SWIG_CLASS_WRAPPER).cmemflags & SWIG_MEM_CONST) { \
        SWIG_exception_impl(FUNCNAME, SWIG_TypeError, \
            "Cannot pass const " TYPENAME " (class " FNAME ") " \
            "as a mutable reference", \
            RETURNNULL); \
    }


#define SWIG_check_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
  if (!(SWIG_CLASS_WRAPPER).cptr) { \
    SWIG_exception_impl(FUNCNAME, SWIG_TypeError, \
                        "Cannot pass null " TYPENAME " (class " FNAME ") " \
                        "as a reference", RETURNNULL); \
  }


#define SWIG_check_mutable_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
    SWIG_check_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL); \
    SWIG_check_mutable(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL);


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


#include "arkode/arkode_splittingstep.h"


typedef struct {
    void* cptr;
    int cmemflags;
} SwigClassWrapper;


SWIGINTERN SwigClassWrapper SwigClassWrapper_uninitialized() {
    SwigClassWrapper result;
    result.cptr = NULL;
    result.cmemflags = 0;
    return result;
}


#include <stdlib.h>
#ifdef _MSC_VER
# ifndef strtoull
#  define strtoull _strtoui64
# endif
# ifndef strtoll
#  define strtoll _strtoi64
# endif
#endif


#include <string.h>


SWIGINTERN void SWIG_assign(SwigClassWrapper* self, SwigClassWrapper other) {
  if (self->cptr == NULL) {
    /* LHS is unassigned */
    if (other.cmemflags & SWIG_MEM_RVALUE) {
      /* Capture pointer from RHS, clear 'moving' flag */
      self->cptr = other.cptr;
      self->cmemflags = other.cmemflags & (~SWIG_MEM_RVALUE);
    } else {
      /* Become a reference to the other object */
      self->cptr = other.cptr;
      self->cmemflags = other.cmemflags & (~SWIG_MEM_OWN);
    }
  } else if (other.cptr == NULL) {
    /* Replace LHS with a null pointer */
    free(self->cptr);
    *self = SwigClassWrapper_uninitialized();
  } else {
    if (self->cmemflags & SWIG_MEM_OWN) {
      free(self->cptr);
    }
    self->cptr = other.cptr;
    if (other.cmemflags & SWIG_MEM_RVALUE) {
      /* Capture RHS */
      self->cmemflags = other.cmemflags & ~SWIG_MEM_RVALUE;
    } else {
      /* Point to RHS */
      self->cmemflags = other.cmemflags & ~SWIG_MEM_OWN;
    }
  }
}


typedef struct {
    void* data;
    size_t size;
} SwigArrayWrapper;


SWIGINTERN SwigArrayWrapper SwigArrayWrapper_uninitialized() {
  SwigArrayWrapper result;
  result.data = NULL;
  result.size = 0;
  return result;
}

SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_alpha_set(SwigClassWrapper const *farg1, double *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  sunrealtype *arg2 = (sunrealtype *) 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::alpha", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (sunrealtype *)(farg2);
  if (arg1) (arg1)->alpha = arg2;
}


SWIGEXPORT double * _wrap_SplittingStepCoefficientsMem_alpha_get(SwigClassWrapper const *farg1) {
  double * fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  sunrealtype *result = 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::alpha", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (sunrealtype *) ((arg1)->alpha);
  fresult = result;
  return fresult;
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_beta_set(SwigClassWrapper const *farg1, void *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  sunrealtype ***arg2 = (sunrealtype ***) 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::beta", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (sunrealtype ***)(farg2);
  if (arg1) (arg1)->beta = arg2;
}


SWIGEXPORT void * _wrap_SplittingStepCoefficientsMem_beta_get(SwigClassWrapper const *farg1) {
  void * fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  sunrealtype ***result = 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::beta", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (sunrealtype ***) ((arg1)->beta);
  fresult = result;
  return fresult;
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_sequential_methods_set(SwigClassWrapper const *farg1, int const *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::sequential_methods", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (int)(*farg2);
  if (arg1) (arg1)->sequential_methods = arg2;
}


SWIGEXPORT int _wrap_SplittingStepCoefficientsMem_sequential_methods_get(SwigClassWrapper const *farg1) {
  int fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::sequential_methods", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (int) ((arg1)->sequential_methods);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_stages_set(SwigClassWrapper const *farg1, int const *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::stages", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (int)(*farg2);
  if (arg1) (arg1)->stages = arg2;
}


SWIGEXPORT int _wrap_SplittingStepCoefficientsMem_stages_get(SwigClassWrapper const *farg1) {
  int fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::stages", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (int) ((arg1)->stages);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_partitions_set(SwigClassWrapper const *farg1, int const *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::partitions", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (int)(*farg2);
  if (arg1) (arg1)->partitions = arg2;
}


SWIGEXPORT int _wrap_SplittingStepCoefficientsMem_partitions_get(SwigClassWrapper const *farg1) {
  int fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::partitions", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (int) ((arg1)->partitions);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_order_set(SwigClassWrapper const *farg1, int const *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::order", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  arg2 = (int)(*farg2);
  if (arg1) (arg1)->order = arg2;
}


SWIGEXPORT int _wrap_SplittingStepCoefficientsMem_order_get(SwigClassWrapper const *farg1) {
  int fresult ;
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  int result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::order", return 0);
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  result = (int) ((arg1)->order);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_new_SplittingStepCoefficientsMem() {
  SwigClassWrapper fresult ;
  struct SplittingStepCoefficientsMem *result = 0 ;
  
  result = (struct SplittingStepCoefficientsMem *)calloc(1, sizeof(struct SplittingStepCoefficientsMem));
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (1 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void _wrap_delete_SplittingStepCoefficientsMem(SwigClassWrapper *farg1) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  
  SWIG_check_mutable(*farg1, "struct SplittingStepCoefficientsMem *", "SplittingStepCoefficientsMem", "SplittingStepCoefficientsMem::~SplittingStepCoefficientsMem()", return );
  arg1 = (struct SplittingStepCoefficientsMem *)(farg1->cptr);
  free((char *) arg1);
}


SWIGEXPORT void _wrap_SplittingStepCoefficientsMem_op_assign__(SwigClassWrapper *farg1, SwigClassWrapper const *farg2) {
  struct SplittingStepCoefficientsMem *arg1 = (struct SplittingStepCoefficientsMem *) 0 ;
  struct SplittingStepCoefficientsMem *arg2 = 0 ;
  
  (void)sizeof(arg1);
  (void)sizeof(arg2);
  SWIG_assign(farg1, *farg2);
  
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_Alloc(int const *farg1, int const *farg2, int const *farg3) {
  SwigClassWrapper fresult ;
  int arg1 ;
  int arg2 ;
  int arg3 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  arg2 = (int)(*farg2);
  arg3 = (int)(*farg3);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_Alloc(arg1,arg2,arg3);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_Create(int const *farg1, int const *farg2, int const *farg3, int const *farg4, double *farg5, double *farg6) {
  SwigClassWrapper fresult ;
  int arg1 ;
  int arg2 ;
  int arg3 ;
  int arg4 ;
  sunrealtype *arg5 = (sunrealtype *) 0 ;
  sunrealtype *arg6 = (sunrealtype *) 0 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  arg2 = (int)(*farg2);
  arg3 = (int)(*farg3);
  arg4 = (int)(*farg4);
  arg5 = (sunrealtype *)(farg5);
  arg6 = (sunrealtype *)(farg6);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_Create(arg1,arg2,arg3,arg4,arg5,arg6);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void _wrap_FSplittingStepCoefficients_Free(SwigClassWrapper const *farg1) {
  SplittingStepCoefficients arg1 = (SplittingStepCoefficients) 0 ;
  
  SWIG_check_mutable(*farg1, "SplittingStepCoefficients", "SplittingStepCoefficientsMem", "SplittingStepCoefficients_Free(SplittingStepCoefficients)", return );
  arg1 = (SplittingStepCoefficients)(farg1->cptr);
  SplittingStepCoefficients_Free(arg1);
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_Copy(SwigClassWrapper const *farg1) {
  SwigClassWrapper fresult ;
  SplittingStepCoefficients arg1 = (SplittingStepCoefficients) 0 ;
  SplittingStepCoefficients result;
  
  SWIG_check_mutable(*farg1, "SplittingStepCoefficients", "SplittingStepCoefficientsMem", "SplittingStepCoefficients_Copy(SplittingStepCoefficients)", return SwigClassWrapper_uninitialized());
  arg1 = (SplittingStepCoefficients)(farg1->cptr);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_Copy(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void _wrap_FSplittingStepCoefficients_Write(SwigClassWrapper const *farg1, void *farg2) {
  SplittingStepCoefficients arg1 = (SplittingStepCoefficients) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  
  SWIG_check_mutable(*farg1, "SplittingStepCoefficients", "SplittingStepCoefficientsMem", "SplittingStepCoefficients_Write(SplittingStepCoefficients,FILE *)", return );
  arg1 = (SplittingStepCoefficients)(farg1->cptr);
  arg2 = (FILE *)(farg2);
  SplittingStepCoefficients_Write(arg1,arg2);
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_LoadCoefficients(int const *farg1) {
  SwigClassWrapper fresult ;
  ARKODE_SplittingCoefficientsID arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (ARKODE_SplittingCoefficientsID)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_LoadCoefficients(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_LoadCoefficientsByName(SwigArrayWrapper *farg1) {
  SwigClassWrapper fresult ;
  char *arg1 = (char *) 0 ;
  SplittingStepCoefficients result;
  
  arg1 = (char *)(farg1->data);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_LoadCoefficientsByName((char const *)arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigArrayWrapper _wrap_FSplittingStepCoefficients_IDToName(int const *farg1) {
  SwigArrayWrapper fresult ;
  ARKODE_SplittingCoefficientsID arg1 ;
  char *result = 0 ;
  
  arg1 = (ARKODE_SplittingCoefficientsID)(*farg1);
  result = (char *)SplittingStepCoefficients_IDToName(arg1);
  fresult.size = strlen((const char*)(result));
  fresult.data = (char *)(result);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_LieTrotter(int const *farg1) {
  SwigClassWrapper fresult ;
  int arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_LieTrotter(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_Strang(int const *farg1) {
  SwigClassWrapper fresult ;
  int arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_Strang(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_Parallel(int const *farg1) {
  SwigClassWrapper fresult ;
  int arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_Parallel(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_SymmetricParallel(int const *farg1) {
  SwigClassWrapper fresult ;
  int arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_SymmetricParallel(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_ThirdOrderSuzuki(int const *farg1) {
  SwigClassWrapper fresult ;
  int arg1 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_ThirdOrderSuzuki(arg1);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_TripleJump(int const *farg1, int const *farg2) {
  SwigClassWrapper fresult ;
  int arg1 ;
  int arg2 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  arg2 = (int)(*farg2);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_TripleJump(arg1,arg2);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_FSplittingStepCoefficients_SuzukiFractal(int const *farg1, int const *farg2) {
  SwigClassWrapper fresult ;
  int arg1 ;
  int arg2 ;
  SplittingStepCoefficients result;
  
  arg1 = (int)(*farg1);
  arg2 = (int)(*farg2);
  result = (SplittingStepCoefficients)SplittingStepCoefficients_SuzukiFractal(arg1,arg2);
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (0 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void * _wrap_FSplittingStepCreate(SwigClassWrapper const *farg1, int const *farg2, double const *farg3, N_Vector farg4, void *farg5) {
  void * fresult ;
  SUNStepper *arg1 = (SUNStepper *) 0 ;
  int arg2 ;
  sunrealtype arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  SUNContext arg5 = (SUNContext) 0 ;
  void *result = 0 ;
  
  SWIG_check_mutable(*farg1, "SUNStepper *", "SWIGTYPE_p_SUNStepper", "SplittingStepCreate(SUNStepper *,int,sunrealtype,N_Vector,SUNContext)", return 0);
  arg1 = (SUNStepper *)(farg1->cptr);
  arg2 = (int)(*farg2);
  arg3 = (sunrealtype)(*farg3);
  arg4 = (N_Vector)(farg4);
  arg5 = (SUNContext)(farg5);
  result = (void *)SplittingStepCreate(arg1,arg2,arg3,arg4,arg5);
  fresult = result;
  return fresult;
}


SWIGEXPORT int _wrap_FSplittingStep_SetCoefficients(void *farg1, SwigClassWrapper const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  SplittingStepCoefficients arg2 = (SplittingStepCoefficients) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  SWIG_check_mutable(*farg2, "SplittingStepCoefficients", "SplittingStepCoefficientsMem", "SplittingStep_SetCoefficients(void *,SplittingStepCoefficients)", return 0);
  arg2 = (SplittingStepCoefficients)(farg2->cptr);
  result = (int)SplittingStep_SetCoefficients(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSplittingStep_GetNumEvolves(void *farg1, int const *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  arg3 = (long *)(farg3);
  result = (int)SplittingStep_GetNumEvolves(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}



