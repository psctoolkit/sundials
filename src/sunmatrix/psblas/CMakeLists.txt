# ---------------------------------------------------------------
# Programmer:  F. Durastante @ IAC-CNR
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the psblas SUNMatrix library

INSTALL(CODE "MESSAGE(\"\nInstall SUNMATRIX_PSBLAS\n\")")






INCLUDE(${PROJECT_SOURCE_DIR}/src/sunmatrix/psblas/makeincinputcmake)

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable sunmatrixpsblas_SOURCES with the sources for the SUNMATRIXPSBLAS lib
SET(sunmatrixpsblas_SOURCES sunmatrix_psblas.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXPSBLAS library
SET(shared_SOURCES
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable sunmatrixpsblas_HEADERS with the exported SUNMATRIXPSBLAS header files
SET(sunmatrixpsblas_HEADERS
  ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_psblas.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)





# PSBLAS LINK DIRECTORIES
#SET(PSBLAS_LIBS  -L${AMG_DIR}/lib -lamg_cbind -lamg_prec -L/usr/local/psblas/lib -lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base -lpsb_cbind)
LINK_DIRECTORIES(${PSBLAS_LIBRARY_DIR})
LINK_DIRECTORIES(${PSCTOOLKIT_DIR} ${AMG_DIR})

LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)
#LINK_DIRECTORIES(${PSBLAS_LIBRARY_DIR})
#LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib)
LINK_LIBRARIES(${PSBLDLIBS} ${PSBLAS_LIBS} -lgfortran)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXPSBLAS library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXPSBLAS library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunmatrixpsblas_static STATIC ${sunmatrixpsblas_SOURCES} ${shared_SOURCES})
  TARGET_COMPILE_DEFINITIONS(sundials_sunmatrixpsblas_static PUBLIC ${PSBCDEFINES} ${PSBFDEFINES})
  TARGET_LINK_LIBRARIES(sundials_sunmatrixpsblas_static ${LINKED_LIBRARIES})
  #TARGET_INCLUDE_DIRECTORIES(sundials_sunmatrixpsblas_static PUBLIC ${AMG_INCLUDE} ${AMG_MODULES})
  
  TARGET_INCLUDE_DIRECTORIES(sundials_sunmatrixpsblas_static PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include ${PSBLAS_INCLUDE} ${AMG_INCLUDE} ${PSBLAS_MODULES} ${AMG_MODULES} ${PSCTOOLKIT_DIR} ${AMG_DIR} ${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)

  SET_TARGET_PROPERTIES(sundials_sunmatrixpsblas_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixpsblas CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunmatrixpsblas_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
  MESSAGE(STATUS "Added to SUNMATRIX_PSBLAS static libs --- ${PSBLAS_INCLUDE} ${PSBLAS_MODULES} ${AMG_INCLUDE} ${AMG_MODULES};")
  MESSAGE(STATUS "Added to SUNMATRIX_PSBLAS static libs ------  ${PSBCDEFINES} ${PSBFDEFINES}) -- ${PSBLAS_LIBRARY_DIR}; --- ${PSBLDLIBS} ${PSBLAS_LIBS}; __ ${LINKED_LIBRARIES};")
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXPSBLAS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXPSBLAS library
IF(BUILD_SHARED_LIBS AND FALSE)
  ADD_LIBRARY(sundials_sunmatrixpsblas_shared SHARED ${sunmatrixpsblas_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunmatrixpsblas_shared m)
  ENDIF()

  TARGET_COMPILE_DEFINITIONS(sundials_sunmatrixpsblas_shared PUBLIC ${PSBCDEFINES} ${PSBFDEFINES})
  TARGET_LINK_LIBRARIES(sundials_sunmatrixpsblas_shared ${LINKED_LIBRARIES})
  TARGET_INCLUDE_DIRECTORIES(sundials_sunmatrixpsblas_shared PUBLIC ${PSBLAS_INCLUDE} ${PSBLAS_MODULES})

  SET_TARGET_PROPERTIES(sundials_sunmatrixpsblas_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixpsblas CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunmatrixpsblas_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  INSTALL(TARGETS sundials_sunmatrixpsblas_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS AND FALSE)

# Install the SUNMATRIXPSBLAS header files
INSTALL(FILES ${sunmatrixpsblas_HEADERS} DESTINATION include/sunmatrix)

#
MESSAGE(STATUS "Added SUNMATRIX_PSBLAS module")
