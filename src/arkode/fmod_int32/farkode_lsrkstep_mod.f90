! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.

! ---------------------------------------------------------------
! Programmer(s): Auto-generated by swig.
! ---------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2024, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ---------------------------------------------------------------

module farkode_lsrkstep_mod
 use, intrinsic :: ISO_C_BINDING
 use farkode_mod
 use fsundials_core_mod
 implicit none
 private

 ! DECLARATION CONSTRUCTS
 public :: FLSRKStepCreate
 public :: FLSRKodeSetMethod
 public :: FLSRKodeSetSprRadFn
 public :: FLSRKodeSetConstJac
 public :: FLSRKodeSetSprRadFrequency
 public :: FLSRKodeSetMaxStageNum
 public :: FLSRKodeSetSprRadSafetyFactor
 public :: FLSRKStepReInit
 public :: FLSRKStepGetNumRhsEvals
 public :: FLSRKStepGetTimestepperStats

! WRAPPER DECLARATIONS
interface
function swigc_FLSRKStepCreate(farg1, farg2, farg3, farg4, farg5) &
bind(C, name="_wrap_FLSRKStepCreate") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_FUNPTR), value :: farg1
type(C_FUNPTR), value :: farg2
real(C_DOUBLE), intent(in) :: farg3
type(C_PTR), value :: farg4
type(C_PTR), value :: farg5
type(C_PTR) :: fresult
end function

function swigc_FLSRKodeSetMethod(farg1, farg2) &
bind(C, name="_wrap_FLSRKodeSetMethod") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKodeSetSprRadFn(farg1, farg2) &
bind(C, name="_wrap_FLSRKodeSetSprRadFn") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_FUNPTR), value :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKodeSetConstJac(farg1) &
bind(C, name="_wrap_FLSRKodeSetConstJac") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FLSRKodeSetSprRadFrequency(farg1, farg2) &
bind(C, name="_wrap_FLSRKodeSetSprRadFrequency") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKodeSetMaxStageNum(farg1, farg2) &
bind(C, name="_wrap_FLSRKodeSetMaxStageNum") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKodeSetSprRadSafetyFactor(farg1, farg2) &
bind(C, name="_wrap_FLSRKodeSetSprRadSafetyFactor") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKStepReInit(farg1, farg2, farg3, farg4, farg5) &
bind(C, name="_wrap_FLSRKStepReInit") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_FUNPTR), value :: farg2
type(C_FUNPTR), value :: farg3
real(C_DOUBLE), intent(in) :: farg4
type(C_PTR), value :: farg5
integer(C_INT) :: fresult
end function

function swigc_FLSRKStepGetNumRhsEvals(farg1, farg2) &
bind(C, name="_wrap_FLSRKStepGetNumRhsEvals") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
integer(C_INT) :: fresult
end function

function swigc_FLSRKStepGetTimestepperStats(farg1, farg2, farg3, farg4, farg5, farg6, farg7, farg8, farg9, farg10, farg11) &
bind(C, name="_wrap_FLSRKStepGetTimestepperStats") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_PTR), value :: farg3
type(C_PTR), value :: farg4
type(C_PTR), value :: farg5
type(C_PTR), value :: farg6
type(C_PTR), value :: farg7
type(C_PTR), value :: farg8
type(C_PTR), value :: farg9
type(C_PTR), value :: farg10
type(C_PTR), value :: farg11
integer(C_INT) :: fresult
end function

end interface


contains
 ! MODULE SUBPROGRAMS
function FLSRKStepCreate(fe, fi, t0, y0, sunctx) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
type(C_PTR) :: swig_result
type(C_FUNPTR), intent(in), value :: fe
type(C_FUNPTR), intent(in), value :: fi
real(C_DOUBLE), intent(in) :: t0
type(N_Vector), target, intent(inout) :: y0
type(C_PTR) :: sunctx
type(C_PTR) :: fresult 
type(C_FUNPTR) :: farg1 
type(C_FUNPTR) :: farg2 
real(C_DOUBLE) :: farg3 
type(C_PTR) :: farg4 
type(C_PTR) :: farg5 

farg1 = fe
farg2 = fi
farg3 = t0
farg4 = c_loc(y0)
farg5 = sunctx
fresult = swigc_FLSRKStepCreate(farg1, farg2, farg3, farg4, farg5)
swig_result = fresult
end function

function FLSRKodeSetMethod(arkode_mem, method) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_INT), intent(in) :: method
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = arkode_mem
farg2 = method
fresult = swigc_FLSRKodeSetMethod(farg1, farg2)
swig_result = fresult
end function

function FLSRKodeSetSprRadFn(arkode_mem, spr) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
type(C_FUNPTR), intent(in), value :: spr
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_FUNPTR) :: farg2 

farg1 = arkode_mem
farg2 = spr
fresult = swigc_FLSRKodeSetSprRadFn(farg1, farg2)
swig_result = fresult
end function

function FLSRKodeSetConstJac(arkode_mem) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = arkode_mem
fresult = swigc_FLSRKodeSetConstJac(farg1)
swig_result = fresult
end function

function FLSRKodeSetSprRadFrequency(arkode_mem, nsteps) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_INT), intent(in) :: nsteps
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = arkode_mem
farg2 = nsteps
fresult = swigc_FLSRKodeSetSprRadFrequency(farg1, farg2)
swig_result = fresult
end function

function FLSRKodeSetMaxStageNum(arkode_mem, stagemaxlimit) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_INT), intent(in) :: stagemaxlimit
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = arkode_mem
farg2 = stagemaxlimit
fresult = swigc_FLSRKodeSetMaxStageNum(farg1, farg2)
swig_result = fresult
end function

function FLSRKodeSetSprRadSafetyFactor(arkode_mem, sprsfty) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
real(C_DOUBLE), intent(in) :: sprsfty
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
real(C_DOUBLE) :: farg2 

farg1 = arkode_mem
farg2 = sprsfty
fresult = swigc_FLSRKodeSetSprRadSafetyFactor(farg1, farg2)
swig_result = fresult
end function

function FLSRKStepReInit(arkode_mem, fe, fi, t0, y0) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
type(C_FUNPTR), intent(in), value :: fe
type(C_FUNPTR), intent(in), value :: fi
real(C_DOUBLE), intent(in) :: t0
type(N_Vector), target, intent(inout) :: y0
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_FUNPTR) :: farg2 
type(C_FUNPTR) :: farg3 
real(C_DOUBLE) :: farg4 
type(C_PTR) :: farg5 

farg1 = arkode_mem
farg2 = fe
farg3 = fi
farg4 = t0
farg5 = c_loc(y0)
fresult = swigc_FLSRKStepReInit(farg1, farg2, farg3, farg4, farg5)
swig_result = fresult
end function

function FLSRKStepGetNumRhsEvals(arkode_mem, nfevals) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_LONG), dimension(*), target, intent(inout) :: nfevals
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 

farg1 = arkode_mem
farg2 = c_loc(nfevals(1))
fresult = swigc_FLSRKStepGetNumRhsEvals(farg1, farg2)
swig_result = fresult
end function

function FLSRKStepGetTimestepperStats(arkode_mem, expsteps, accsteps, attempts, fevals, sprfevals, netfails, stagemax, &
  nsprupdates, sprmax, sprmin) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(C_PTR) :: arkode_mem
integer(C_LONG), dimension(*), target, intent(inout) :: expsteps
integer(C_LONG), dimension(*), target, intent(inout) :: accsteps
integer(C_LONG), dimension(*), target, intent(inout) :: attempts
integer(C_LONG), dimension(*), target, intent(inout) :: fevals
integer(C_LONG), dimension(*), target, intent(inout) :: sprfevals
integer(C_LONG), dimension(*), target, intent(inout) :: netfails
integer(C_LONG), dimension(*), target, intent(inout) :: stagemax
integer(C_LONG), dimension(*), target, intent(inout) :: nsprupdates
real(C_DOUBLE), dimension(*), target, intent(inout) :: sprmax
real(C_DOUBLE), dimension(*), target, intent(inout) :: sprmin
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_PTR) :: farg3 
type(C_PTR) :: farg4 
type(C_PTR) :: farg5 
type(C_PTR) :: farg6 
type(C_PTR) :: farg7 
type(C_PTR) :: farg8 
type(C_PTR) :: farg9 
type(C_PTR) :: farg10 
type(C_PTR) :: farg11 

farg1 = arkode_mem
farg2 = c_loc(expsteps(1))
farg3 = c_loc(accsteps(1))
farg4 = c_loc(attempts(1))
farg5 = c_loc(fevals(1))
farg6 = c_loc(sprfevals(1))
farg7 = c_loc(netfails(1))
farg8 = c_loc(stagemax(1))
farg9 = c_loc(nsprupdates(1))
farg10 = c_loc(sprmax(1))
farg11 = c_loc(sprmin(1))
fresult = swigc_FLSRKStepGetTimestepperStats(farg1, farg2, farg3, farg4, farg5, farg6, farg7, farg8, farg9, farg10, farg11)
swig_result = fresult
end function


end module
