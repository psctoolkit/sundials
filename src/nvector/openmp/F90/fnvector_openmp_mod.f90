! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.
module fnvector_openmp_mod
 use, intrinsic :: ISO_C_BINDING
 use fsundials_types
 use fnvector_mod
 implicit none
 private

 ! PUBLIC METHODS AND TYPES

  public :: FN_VGetData_OpenMP
  
 public :: FN_VNew_OpenMP
 public :: FN_VNewEmpty_OpenMP
 public :: FN_VMake_OpenMP
 public :: FN_VCloneVectorArray_OpenMP
 public :: FN_VCloneVectorArrayEmpty_OpenMP
 public :: FN_VDestroyVectorArray_OpenMP
 public :: FN_VGetLength_OpenMP
 public :: FN_VPrint_OpenMP
 public :: FN_VPrintFile_OpenMP
 public :: FN_VGetVectorID_OpenMP
 public :: FN_VCloneEmpty_OpenMP
 public :: FN_VClone_OpenMP
 public :: FN_VDestroy_OpenMP
 public :: FN_VSpace_OpenMP
 public :: FN_VGetArrayPointer_OpenMP
 public :: FN_VSetArrayPointer_OpenMP
 public :: FN_VLinearSum_OpenMP
 public :: FN_VConst_OpenMP
 public :: FN_VProd_OpenMP
 public :: FN_VDiv_OpenMP
 public :: FN_VScale_OpenMP
 public :: FN_VAbs_OpenMP
 public :: FN_VInv_OpenMP
 public :: FN_VAddConst_OpenMP
 public :: FN_VDotProd_OpenMP
 public :: FN_VMaxNorm_OpenMP
 public :: FN_VWrmsNorm_OpenMP
 public :: FN_VWrmsNormMask_OpenMP
 public :: FN_VMin_OpenMP
 public :: FN_VWL2Norm_OpenMP
 public :: FN_VL1Norm_OpenMP
 public :: FN_VCompare_OpenMP
 public :: FN_VInvTest_OpenMP
 public :: FN_VConstrMask_OpenMP
 public :: FN_VMinQuotient_OpenMP
 public :: FN_VLinearCombination_OpenMP
 public :: FN_VScaleAddMulti_OpenMP
 public :: FN_VDotProdMulti_OpenMP
 public :: FN_VLinearSumVectorArray_OpenMP
 public :: FN_VScaleVectorArray_OpenMP
 public :: FN_VConstVectorArray_OpenMP
 public :: FN_VWrmsNormVectorArray_OpenMP
 public :: FN_VWrmsNormMaskVectorArray_OpenMP
 public :: FN_VScaleAddMultiVectorArray_OpenMP
 public :: FN_VLinearCombinationVectorArray_OpenMP
 public :: FN_VEnableFusedOps_OpenMP
 public :: FN_VEnableLinearCombination_OpenMP
 public :: FN_VEnableScaleAddMulti_OpenMP
 public :: FN_VEnableDotProdMulti_OpenMP
 public :: FN_VEnableLinearSumVectorArray_OpenMP
 public :: FN_VEnableScaleVectorArray_OpenMP
 public :: FN_VEnableConstVectorArray_OpenMP
 public :: FN_VEnableWrmsNormVectorArray_OpenMP
 public :: FN_VEnableWrmsNormMaskVectorArray_OpenMP
 public :: FN_VEnableScaleAddMultiVectorArray_OpenMP
 public :: FN_VEnableLinearCombinationVectorArray_OpenMP

 ! WRAPPER DECLARATIONS
 interface
function FN_VNew_OpenMP(vec_length, num_threads) &
bind(C, name="N_VNew_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT64_T), value :: vec_length
integer(C_INT), value :: num_threads
type(C_PTR) :: fresult
end function

function FN_VNewEmpty_OpenMP(vec_length, num_threads) &
bind(C, name="N_VNewEmpty_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT64_T), value :: vec_length
integer(C_INT), value :: num_threads
type(C_PTR) :: fresult
end function

function FN_VMake_OpenMP(vec_length, v_data, num_threads) &
bind(C, name="N_VMake_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT64_T), value :: vec_length
real(C_DOUBLE), dimension(*) :: v_data
integer(C_INT), value :: num_threads
type(C_PTR) :: fresult
end function

function FN_VCloneVectorArray_OpenMP(count, w) &
bind(C, name="N_VCloneVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: count
type(C_PTR), value :: w
type(C_PTR) :: fresult
end function

function FN_VCloneVectorArrayEmpty_OpenMP(count, w) &
bind(C, name="N_VCloneVectorArrayEmpty_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: count
type(C_PTR), value :: w
type(C_PTR) :: fresult
end function

subroutine FN_VDestroyVectorArray_OpenMP(vs, count) &
bind(C, name="N_VDestroyVectorArray_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: vs
integer(C_INT), value :: count
end subroutine

function FN_VGetLength_OpenMP(v) &
bind(C, name="N_VGetLength_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
integer(C_INT64_T) :: fresult
end function

subroutine FN_VPrint_OpenMP(v) &
bind(C, name="N_VPrint_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
end subroutine

subroutine FN_VPrintFile_OpenMP(v, outfile) &
bind(C, name="N_VPrintFile_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
type(C_PTR), value :: outfile
end subroutine

function FN_VGetVectorID_OpenMP(v) &
bind(C, name="N_VGetVectorID_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
integer(C_INT) :: fresult
end function

function FN_VCloneEmpty_OpenMP(w) &
bind(C, name="N_VCloneEmpty_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: w
type(C_PTR) :: fresult
end function

function FN_VClone_OpenMP(w) &
bind(C, name="N_VClone_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: w
type(C_PTR) :: fresult
end function

subroutine FN_VDestroy_OpenMP(v) &
bind(C, name="N_VDestroy_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
end subroutine

subroutine FN_VSpace_OpenMP(v, lrw, liw) &
bind(C, name="N_VSpace_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
integer(C_INT64_T) :: lrw
integer(C_INT64_T) :: liw
end subroutine

function FN_VGetArrayPointer_OpenMP(v) &
bind(C, name="N_VGetArrayPointer_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
type(C_PTR) :: fresult
end function

subroutine FN_VSetArrayPointer_OpenMP(v_data, v) &
bind(C, name="N_VSetArrayPointer_OpenMP")
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), dimension(*) :: v_data
type(C_PTR), value :: v
end subroutine

subroutine FN_VLinearSum_OpenMP(a, x, b, y, z) &
bind(C, name="N_VLinearSum_OpenMP")
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), value :: a
type(C_PTR), value :: x
real(C_DOUBLE), value :: b
type(C_PTR), value :: y
type(C_PTR), value :: z
end subroutine

subroutine FN_VConst_OpenMP(c, z) &
bind(C, name="N_VConst_OpenMP")
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), value :: c
type(C_PTR), value :: z
end subroutine

subroutine FN_VProd_OpenMP(x, y, z) &
bind(C, name="N_VProd_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: y
type(C_PTR), value :: z
end subroutine

subroutine FN_VDiv_OpenMP(x, y, z) &
bind(C, name="N_VDiv_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: y
type(C_PTR), value :: z
end subroutine

subroutine FN_VScale_OpenMP(c, x, z) &
bind(C, name="N_VScale_OpenMP")
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), value :: c
type(C_PTR), value :: x
type(C_PTR), value :: z
end subroutine

subroutine FN_VAbs_OpenMP(x, z) &
bind(C, name="N_VAbs_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: z
end subroutine

subroutine FN_VInv_OpenMP(x, z) &
bind(C, name="N_VInv_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: z
end subroutine

subroutine FN_VAddConst_OpenMP(x, b, z) &
bind(C, name="N_VAddConst_OpenMP")
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
real(C_DOUBLE), value :: b
type(C_PTR), value :: z
end subroutine

function FN_VDotProd_OpenMP(x, y) &
bind(C, name="N_VDotProd_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: y
real(C_DOUBLE) :: fresult
end function

function FN_VMaxNorm_OpenMP(x) &
bind(C, name="N_VMaxNorm_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
real(C_DOUBLE) :: fresult
end function

function FN_VWrmsNorm_OpenMP(x, w) &
bind(C, name="N_VWrmsNorm_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: w
real(C_DOUBLE) :: fresult
end function

function FN_VWrmsNormMask_OpenMP(x, w, id) &
bind(C, name="N_VWrmsNormMask_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: w
type(C_PTR), value :: id
real(C_DOUBLE) :: fresult
end function

function FN_VMin_OpenMP(x) &
bind(C, name="N_VMin_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
real(C_DOUBLE) :: fresult
end function

function FN_VWL2Norm_OpenMP(x, w) &
bind(C, name="N_VWL2Norm_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: w
real(C_DOUBLE) :: fresult
end function

function FN_VL1Norm_OpenMP(x) &
bind(C, name="N_VL1Norm_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
real(C_DOUBLE) :: fresult
end function

subroutine FN_VCompare_OpenMP(c, x, z) &
bind(C, name="N_VCompare_OpenMP")
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE), value :: c
type(C_PTR), value :: x
type(C_PTR), value :: z
end subroutine

function FN_VInvTest_OpenMP(x, z) &
bind(C, name="N_VInvTest_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: x
type(C_PTR), value :: z
logical(C_BOOL) :: fresult
end function

function FN_VConstrMask_OpenMP(c, x, m) &
bind(C, name="N_VConstrMask_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: c
type(C_PTR), value :: x
type(C_PTR), value :: m
logical(C_BOOL) :: fresult
end function

function FN_VMinQuotient_OpenMP(num, denom) &
bind(C, name="N_VMinQuotient_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: num
type(C_PTR), value :: denom
real(C_DOUBLE) :: fresult
end function

function FN_VLinearCombination_OpenMP(nvec, c, v, z) &
bind(C, name="N_VLinearCombination_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
real(C_DOUBLE) :: c
type(C_PTR), value :: v
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VScaleAddMulti_OpenMP(nvec, a, x, y, z) &
bind(C, name="N_VScaleAddMulti_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
real(C_DOUBLE) :: a
type(C_PTR), value :: x
type(C_PTR), value :: y
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VDotProdMulti_OpenMP(nvec, x, y, dotprods) &
bind(C, name="N_VDotProdMulti_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
type(C_PTR), value :: x
type(C_PTR), value :: y
real(C_DOUBLE) :: dotprods
integer(C_INT) :: fresult
end function

function FN_VLinearSumVectorArray_OpenMP(nvec, a, x, b, y, z) &
bind(C, name="N_VLinearSumVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
real(C_DOUBLE), value :: a
type(C_PTR), value :: x
real(C_DOUBLE), value :: b
type(C_PTR), value :: y
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VScaleVectorArray_OpenMP(nvec, c, x, z) &
bind(C, name="N_VScaleVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
real(C_DOUBLE) :: c
type(C_PTR), value :: x
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VConstVectorArray_OpenMP(nvecs, c, z) &
bind(C, name="N_VConstVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvecs
real(C_DOUBLE), value :: c
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VWrmsNormVectorArray_OpenMP(nvecs, x, w, nrm) &
bind(C, name="N_VWrmsNormVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvecs
type(C_PTR), value :: x
type(C_PTR), value :: w
real(C_DOUBLE) :: nrm
integer(C_INT) :: fresult
end function

function FN_VWrmsNormMaskVectorArray_OpenMP(nvecs, x, w, id, nrm) &
bind(C, name="N_VWrmsNormMaskVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvecs
type(C_PTR), value :: x
type(C_PTR), value :: w
type(C_PTR), value :: id
real(C_DOUBLE) :: nrm
integer(C_INT) :: fresult
end function

function FN_VScaleAddMultiVectorArray_OpenMP(nvec, nsum, a, x, y, z) &
bind(C, name="N_VScaleAddMultiVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
integer(C_INT), value :: nsum
real(C_DOUBLE) :: a
type(C_PTR), value :: x
type(C_PTR), value :: y
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VLinearCombinationVectorArray_OpenMP(nvec, nsum, c, x, z) &
bind(C, name="N_VLinearCombinationVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
integer(C_INT), value :: nvec
integer(C_INT), value :: nsum
real(C_DOUBLE) :: c
type(C_PTR), value :: x
type(C_PTR), value :: z
integer(C_INT) :: fresult
end function

function FN_VEnableFusedOps_OpenMP(v, tf) &
bind(C, name="N_VEnableFusedOps_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableLinearCombination_OpenMP(v, tf) &
bind(C, name="N_VEnableLinearCombination_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableScaleAddMulti_OpenMP(v, tf) &
bind(C, name="N_VEnableScaleAddMulti_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableDotProdMulti_OpenMP(v, tf) &
bind(C, name="N_VEnableDotProdMulti_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableLinearSumVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableLinearSumVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableScaleVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableScaleVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableConstVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableConstVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableWrmsNormVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableWrmsNormVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableWrmsNormMaskVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableWrmsNormMaskVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableScaleAddMultiVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableScaleAddMultiVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

function FN_VEnableLinearCombinationVectorArray_OpenMP(v, tf) &
bind(C, name="N_VEnableLinearCombinationVectorArray_OpenMP") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: v
logical(C_BOOL), value :: tf
integer(C_INT) :: fresult
end function

 end interface


contains
 ! FORTRAN PROXY CODE

  subroutine FN_VGetData_OpenMP(vec, vdata)

      use, intrinsic :: iso_c_binding
      implicit none

      type(C_PTR)        :: vec
      integer(C_INT64_T) :: len
      type(C_PTR)        :: cptr
      real(C_DOUBLE), dimension(:), pointer :: vdata

      len = FN_VGetLength_OpenMP(vec)
      cptr = FN_VGetArrayPointer_OpenMP(vec)

      call c_f_pointer(cptr, vdata, (/len/))

  end subroutine FN_VGetData_OpenMP
  

end module
