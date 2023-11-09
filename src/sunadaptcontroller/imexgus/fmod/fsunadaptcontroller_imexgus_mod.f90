! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.

! ---------------------------------------------------------------
! Programmer(s): Auto-generated by swig.
! ---------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2023, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ---------------------------------------------------------------

module fsunadaptcontroller_imexgus_mod
 use, intrinsic :: ISO_C_BINDING
 use fsundials_adaptcontroller_mod
 use fsundials_types_mod
 use fsundials_context_mod
 implicit none
 private

 ! DECLARATION CONSTRUCTS
 public :: FSUNAdaptController_ImExGus
 public :: FSUNAdaptController_SetParams_ImExGus
 public :: FSUNAdaptController_GetType_ImExGus
 public :: FSUNAdaptController_EstimateStep_ImExGus
 public :: FSUNAdaptController_Reset_ImExGus
 public :: FSUNAdaptController_SetDefaults_ImExGus
 public :: FSUNAdaptController_Write_ImExGus
 public :: FSUNAdaptController_SetErrorBias_ImExGus
 public :: FSUNAdaptController_UpdateH_ImExGus
 public :: FSUNAdaptController_Space_ImExGus

! WRAPPER DECLARATIONS
interface
function swigc_FSUNAdaptController_ImExGus(farg1) &
bind(C, name="_wrap_FSUNAdaptController_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNAdaptController_SetParams_ImExGus(farg1, farg2, farg3, farg4, farg5) &
bind(C, name="_wrap_FSUNAdaptController_SetParams_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE), intent(in) :: farg2
real(C_DOUBLE), intent(in) :: farg3
real(C_DOUBLE), intent(in) :: farg4
real(C_DOUBLE), intent(in) :: farg5
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_GetType_ImExGus(farg1) &
bind(C, name="_wrap_FSUNAdaptController_GetType_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_EstimateStep_ImExGus(farg1, farg2, farg3, farg4, farg5) &
bind(C, name="_wrap_FSUNAdaptController_EstimateStep_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE), intent(in) :: farg2
integer(C_INT), intent(in) :: farg3
real(C_DOUBLE), intent(in) :: farg4
type(C_PTR), value :: farg5
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_Reset_ImExGus(farg1) &
bind(C, name="_wrap_FSUNAdaptController_Reset_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_SetDefaults_ImExGus(farg1) &
bind(C, name="_wrap_FSUNAdaptController_SetDefaults_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_Write_ImExGus(farg1, farg2) &
bind(C, name="_wrap_FSUNAdaptController_Write_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_SetErrorBias_ImExGus(farg1, farg2) &
bind(C, name="_wrap_FSUNAdaptController_SetErrorBias_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_UpdateH_ImExGus(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNAdaptController_UpdateH_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE), intent(in) :: farg2
real(C_DOUBLE), intent(in) :: farg3
integer(C_INT) :: fresult
end function

function swigc_FSUNAdaptController_Space_ImExGus(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNAdaptController_Space_ImExGus") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_PTR), value :: farg3
integer(C_INT) :: fresult
end function

end interface


contains
 ! MODULE SUBPROGRAMS
function FSUNAdaptController_ImExGus(sunctx) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
type(SUNAdaptController), pointer :: swig_result
type(C_PTR) :: sunctx
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 

farg1 = sunctx
fresult = swigc_FSUNAdaptController_ImExGus(farg1)
call c_f_pointer(fresult, swig_result)
end function

function FSUNAdaptController_SetParams_ImExGus(c, k1e, k2e, k1i, k2i) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
real(C_DOUBLE), intent(in) :: k1e
real(C_DOUBLE), intent(in) :: k2e
real(C_DOUBLE), intent(in) :: k1i
real(C_DOUBLE), intent(in) :: k2i
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
real(C_DOUBLE) :: farg2 
real(C_DOUBLE) :: farg3 
real(C_DOUBLE) :: farg4 
real(C_DOUBLE) :: farg5 

farg1 = c_loc(c)
farg2 = k1e
farg3 = k2e
farg4 = k1i
farg5 = k2i
fresult = swigc_FSUNAdaptController_SetParams_ImExGus(farg1, farg2, farg3, farg4, farg5)
swig_result = fresult
end function

function FSUNAdaptController_GetType_ImExGus(c) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(SUNAdaptController_Type) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(c)
fresult = swigc_FSUNAdaptController_GetType_ImExGus(farg1)
swig_result = fresult
end function

function FSUNAdaptController_EstimateStep_ImExGus(c, h, p, dsm, hnew) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
real(C_DOUBLE), intent(in) :: h
integer(C_INT), intent(in) :: p
real(C_DOUBLE), intent(in) :: dsm
real(C_DOUBLE), target, intent(inout) :: hnew
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
real(C_DOUBLE) :: farg2 
integer(C_INT) :: farg3 
real(C_DOUBLE) :: farg4 
type(C_PTR) :: farg5 

farg1 = c_loc(c)
farg2 = h
farg3 = p
farg4 = dsm
farg5 = c_loc(hnew)
fresult = swigc_FSUNAdaptController_EstimateStep_ImExGus(farg1, farg2, farg3, farg4, farg5)
swig_result = fresult
end function

function FSUNAdaptController_Reset_ImExGus(c) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(c)
fresult = swigc_FSUNAdaptController_Reset_ImExGus(farg1)
swig_result = fresult
end function

function FSUNAdaptController_SetDefaults_ImExGus(c) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(c)
fresult = swigc_FSUNAdaptController_SetDefaults_ImExGus(farg1)
swig_result = fresult
end function

function FSUNAdaptController_Write_ImExGus(c, fptr) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
type(C_PTR) :: fptr
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 

farg1 = c_loc(c)
farg2 = fptr
fresult = swigc_FSUNAdaptController_Write_ImExGus(farg1, farg2)
swig_result = fresult
end function

function FSUNAdaptController_SetErrorBias_ImExGus(c, bias) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
real(C_DOUBLE), intent(in) :: bias
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
real(C_DOUBLE) :: farg2 

farg1 = c_loc(c)
farg2 = bias
fresult = swigc_FSUNAdaptController_SetErrorBias_ImExGus(farg1, farg2)
swig_result = fresult
end function

function FSUNAdaptController_UpdateH_ImExGus(c, h, dsm) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
real(C_DOUBLE), intent(in) :: h
real(C_DOUBLE), intent(in) :: dsm
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
real(C_DOUBLE) :: farg2 
real(C_DOUBLE) :: farg3 

farg1 = c_loc(c)
farg2 = h
farg3 = dsm
fresult = swigc_FSUNAdaptController_UpdateH_ImExGus(farg1, farg2, farg3)
swig_result = fresult
end function

function FSUNAdaptController_Space_ImExGus(c, lenrw, leniw) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNAdaptController), target, intent(inout) :: c
integer(C_LONG), dimension(*), target, intent(inout) :: lenrw
integer(C_LONG), dimension(*), target, intent(inout) :: leniw
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_PTR) :: farg3 

farg1 = c_loc(c)
farg2 = c_loc(lenrw(1))
farg3 = c_loc(leniw(1))
fresult = swigc_FSUNAdaptController_Space_ImExGus(farg1, farg2, farg3)
swig_result = fresult
end function


end module
