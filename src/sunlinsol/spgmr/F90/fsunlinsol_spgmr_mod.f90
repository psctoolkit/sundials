! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.

! ---------------------------------------------------------------
! Programmer(s): Auto-generated by swig.
! ---------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ---------------------------------------------------------------

module fsunlinsol_spgmr_mod
 use, intrinsic :: ISO_C_BINDING
 use fsunlinsol_mod
 use fsundials_types
 use fnvector_mod
 use fsundials_types
 implicit none
 private

 ! PUBLIC METHODS AND TYPES
 public :: FSUNLinSol_SPGMR
 public :: FSUNLinSol_SPGMRSetPrecType
 public :: FSUNLinSol_SPGMRSetGSType
 public :: FSUNLinSol_SPGMRSetMaxRestarts
 public :: FSUNSPGMR
 public :: FSUNSPGMRSetPrecType
 public :: FSUNSPGMRSetGSType
 public :: FSUNSPGMRSetMaxRestarts
 public :: FSUNLinSolGetType_SPGMR
 public :: FSUNLinSolInitialize_SPGMR
 public :: FSUNLinSolSetATimes_SPGMR
 public :: FSUNLinSolSetPreconditioner_SPGMR
 public :: FSUNLinSolSetScalingVectors_SPGMR
 public :: FSUNLinSolSetup_SPGMR
 public :: FSUNLinSolSolve_SPGMR
 public :: FSUNLinSolNumIters_SPGMR
 public :: FSUNLinSolResNorm_SPGMR
 public :: FSUNLinSolResid_SPGMR
 public :: FSUNLinSolLastFlag_SPGMR
 public :: FSUNLinSolSpace_SPGMR
 public :: FSUNLinSolFree_SPGMR

 ! PARAMETERS
 integer(C_INT), parameter, public :: SUNSPGMR_MAXL_DEFAULT = 5_C_INT
 integer(C_INT), parameter, public :: SUNSPGMR_MAXRS_DEFAULT = 0_C_INT

 ! WRAPPER DECLARATIONS
 interface
function FSUNLinSol_SPGMR(y, pretype, maxl) &
bind(C, name="SUNLinSol_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
integer(C_INT), value :: pretype
integer(C_INT), value :: maxl
type(C_PTR) :: fresult
end function

function FSUNLinSol_SPGMRSetPrecType(s, pretype) &
bind(C, name="SUNLinSol_SPGMRSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: pretype
integer(C_INT) :: fresult
end function

function FSUNLinSol_SPGMRSetGSType(s, gstype) &
bind(C, name="SUNLinSol_SPGMRSetGSType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: gstype
integer(C_INT) :: fresult
end function

function FSUNLinSol_SPGMRSetMaxRestarts(s, maxrs) &
bind(C, name="SUNLinSol_SPGMRSetMaxRestarts") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: maxrs
integer(C_INT) :: fresult
end function

function FSUNSPGMR(y, pretype, maxl) &
bind(C, name="SUNSPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
integer(C_INT), value :: pretype
integer(C_INT), value :: maxl
type(C_PTR) :: fresult
end function

function FSUNSPGMRSetPrecType(s, pretype) &
bind(C, name="SUNSPGMRSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: pretype
integer(C_INT) :: fresult
end function

function FSUNSPGMRSetGSType(s, gstype) &
bind(C, name="SUNSPGMRSetGSType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: gstype
integer(C_INT) :: fresult
end function

function FSUNSPGMRSetMaxRestarts(s, maxrs) &
bind(C, name="SUNSPGMRSetMaxRestarts") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: maxrs
integer(C_INT) :: fresult
end function

function FSUNLinSolGetType_SPGMR(s) &
bind(C, name="SUNLinSolGetType_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolInitialize_SPGMR(s) &
bind(C, name="SUNLinSolInitialize_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolSetATimes_SPGMR(s, a_data, atimes) &
bind(C, name="SUNLinSolSetATimes_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a_data
type(C_FUNPTR), value :: atimes
integer(C_INT) :: fresult
end function

function FSUNLinSolSetPreconditioner_SPGMR(s, p_data, pset, psol) &
bind(C, name="SUNLinSolSetPreconditioner_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: p_data
type(C_FUNPTR), value :: pset
type(C_FUNPTR), value :: psol
integer(C_INT) :: fresult
end function

function FSUNLinSolSetScalingVectors_SPGMR(s, s1, s2) &
bind(C, name="SUNLinSolSetScalingVectors_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: s1
type(C_PTR), value :: s2
integer(C_INT) :: fresult
end function

function FSUNLinSolSetup_SPGMR(s, a) &
bind(C, name="SUNLinSolSetup_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a
integer(C_INT) :: fresult
end function

function FSUNLinSolSolve_SPGMR(s, a, x, b, tol) &
bind(C, name="SUNLinSolSolve_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a
type(C_PTR), value :: x
type(C_PTR), value :: b
real(C_DOUBLE), value :: tol
integer(C_INT) :: fresult
end function

function FSUNLinSolNumIters_SPGMR(s) &
bind(C, name="SUNLinSolNumIters_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolResNorm_SPGMR(s) &
bind(C, name="SUNLinSolResNorm_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
real(C_DOUBLE) :: fresult
end function

function FSUNLinSolResid_SPGMR(s) &
bind(C, name="SUNLinSolResid_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR) :: fresult
end function

function FSUNLinSolLastFlag_SPGMR(s) &
bind(C, name="SUNLinSolLastFlag_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: fresult
end function

function FSUNLinSolSpace_SPGMR(s, lenrwls, leniwls) &
bind(C, name="SUNLinSolSpace_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: lenrwls
integer(C_LONG) :: leniwls
integer(C_INT) :: fresult
end function

function FSUNLinSolFree_SPGMR(s) &
bind(C, name="SUNLinSolFree_SPGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

 end interface


end module