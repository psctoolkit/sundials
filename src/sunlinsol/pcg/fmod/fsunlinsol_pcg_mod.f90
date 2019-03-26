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

module fsunlinsol_pcg_mod
 use, intrinsic :: ISO_C_BINDING
 use fsunlinsol_mod
 use fsundials_types
 use fnvector_mod
 use fsundials_types
 implicit none
 private

 ! DECLARATION CONSTRUCTS
 integer(C_INT), parameter, public :: SUNPCG_MAXL_DEFAULT = 5_C_INT
 public :: FSUNLinSol_PCG
 public :: FSUNLinSol_PCGSetPrecType
 public :: FSUNLinSol_PCGSetMaxl
 public :: FSUNPCG
 public :: FSUNPCGSetPrecType
 public :: FSUNPCGSetMaxl
 public :: FSUNLinSolGetType_PCG
 public :: FSUNLinSolInitialize_PCG
 public :: FSUNLinSolSetATimes_PCG
 public :: FSUNLinSolSetPreconditioner_PCG
 public :: FSUNLinSolSetScalingVectors_PCG
 public :: FSUNLinSolSetup_PCG
 public :: FSUNLinSolSolve_PCG
 public :: FSUNLinSolNumIters_PCG
 public :: FSUNLinSolResNorm_PCG
 public :: FSUNLinSolResid_PCG
 public :: FSUNLinSolLastFlag_PCG
 public :: FSUNLinSolSpace_PCG
 public :: FSUNLinSolFree_PCG

! WRAPPER DECLARATIONS
interface
function FSUNLinSol_PCG(y, pretype, maxl) &
bind(C, name="SUNLinSol_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
integer(C_INT), value :: pretype
integer(C_INT), value :: maxl
type(C_PTR) :: fresult
end function

function FSUNLinSol_PCGSetPrecType(s, pretype) &
bind(C, name="SUNLinSol_PCGSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: pretype
integer(C_INT) :: fresult
end function

function FSUNLinSol_PCGSetMaxl(s, maxl) &
bind(C, name="SUNLinSol_PCGSetMaxl") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: maxl
integer(C_INT) :: fresult
end function

function FSUNPCG(y, pretype, maxl) &
bind(C, name="SUNPCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: y
integer(C_INT), value :: pretype
integer(C_INT), value :: maxl
type(C_PTR) :: fresult
end function

function FSUNPCGSetPrecType(s, pretype) &
bind(C, name="SUNPCGSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: pretype
integer(C_INT) :: fresult
end function

function FSUNPCGSetMaxl(s, maxl) &
bind(C, name="SUNPCGSetMaxl") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT), value :: maxl
integer(C_INT) :: fresult
end function

function FSUNLinSolGetType_PCG(s) &
bind(C, name="SUNLinSolGetType_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolInitialize_PCG(s) &
bind(C, name="SUNLinSolInitialize_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolSetATimes_PCG(s, a_data, atimes) &
bind(C, name="SUNLinSolSetATimes_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: a_data
type(C_FUNPTR), value :: atimes
integer(C_INT) :: fresult
end function

function FSUNLinSolSetPreconditioner_PCG(s, p_data, pset, psol) &
bind(C, name="SUNLinSolSetPreconditioner_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: p_data
type(C_FUNPTR), value :: pset
type(C_FUNPTR), value :: psol
integer(C_INT) :: fresult
end function

function FSUNLinSolSetScalingVectors_PCG(s, s1, nul) &
bind(C, name="SUNLinSolSetScalingVectors_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: s1
type(C_PTR), value :: nul
integer(C_INT) :: fresult
end function

function FSUNLinSolSetup_PCG(s, nul) &
bind(C, name="SUNLinSolSetup_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: nul
integer(C_INT) :: fresult
end function

function FSUNLinSolSolve_PCG(s, nul, x, b, tol) &
bind(C, name="SUNLinSolSolve_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR), value :: nul
type(C_PTR), value :: x
type(C_PTR), value :: b
real(C_DOUBLE), value :: tol
integer(C_INT) :: fresult
end function

function FSUNLinSolNumIters_PCG(s) &
bind(C, name="SUNLinSolNumIters_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

function FSUNLinSolResNorm_PCG(s) &
bind(C, name="SUNLinSolResNorm_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
real(C_DOUBLE) :: fresult
end function

function FSUNLinSolResid_PCG(s) &
bind(C, name="SUNLinSolResid_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
type(C_PTR) :: fresult
end function

function FSUNLinSolLastFlag_PCG(s) &
bind(C, name="SUNLinSolLastFlag_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: fresult
end function

function FSUNLinSolSpace_PCG(s, lenrwls, leniwls) &
bind(C, name="SUNLinSolSpace_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_LONG) :: lenrwls
integer(C_LONG) :: leniwls
integer(C_INT) :: fresult
end function

function FSUNLinSolFree_PCG(s) &
bind(C, name="SUNLinSolFree_PCG") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: s
integer(C_INT) :: fresult
end function

end interface


end module
