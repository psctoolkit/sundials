# ---------------------------------------------------------------
# Programmer:  F. Durastante @ IAC CNR
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
# CMakeLists.txt file for the PSBLASFCG SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_PSBLAS\n\")")


# This part will be removed when PSBLAS/MLD2P4 and the other libraries will be
# installed by using CMake, this works as a patch for including the information
# contained in the Make.inc files (Compile flags etc.), I know it is very
# horrible indeed
#EXECUTE_PROCESS(COMMAND ./make2cmakeset.sh ${PSBLAS_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src/sunlinsol/psblas WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/sunlinsol/psblas OUTPUT_QUIET ERROR_QUIET)
INCLUDE(${PROJECT_SOURCE_DIR}/src/sunlinsol/psblas/makeincinputcmake)


IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable sunlinsolpsblas_SOURCES with the sources for the SUNLINSOLPSBLASFCG lib
SET(sunlinsolpsblas_SOURCES sunlinsol_psblas.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOL_PSBLAS library
SET(shared_SOURCES
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_math.c
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${PROJECT_SOURCE_DIR}/src/sundials/sundials_iterative.c
	)

# Add variable sunlinsolpsblas_HEADERS with the exported SUNLINSOL_PSBLAS header files
SET(sunlinsolpsblas_HEADERS
  ${PROJECT_SOURCE_DIR}/include/sunlinsol/sunlinsol_psblas.h
  ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_psblas.h
  ${PROJECT_SOURCE_DIR}/include/nvector/nvector_psblas.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY ${PSBCDEFINES})

# PSBLAS LINK DIRECTORIES
#SET(PSBLAS_LIBS -lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base -L/usr/local/psblas/lib -lpsb_cbind)
#SET(PSBLAS_LIBRARY_DIR /usr/local/psblas/include)
#SET(AMG_LIBRARY_DIR /usr/local/amg4psblas/include)
# Find all static libraries in the lib directory
#file(GLOB PSBLAS_LIBS1 "/usr/local/psblas/lib/*.a")  # Adjust the pattern if needed

#LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)
#LINK_DIRECTORIES(${PSCTOOLKIT_DIR} ${AMG_DIR})
#LINK_DIRECTORIES(${PSBLAS_LIBRARY_DIR})
#LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib)
#LINK_LIBRARIES(${PSBLDLIBS} ${PSBLAS_LIBS})
#  MESSAGE(STATUS "Added to SUNLINSOL_PSBLAS static libs ------  ${PSBLDLIBS} __ ${PSBLAS_LIBS}) -- ${PSBLAS_LIBS1};")

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOL_PSBLAS library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOL_PSBLAS library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolpsblas_static STATIC ${sunlinsolpsblas_SOURCES} ${shared_SOURCES})
  TARGET_COMPILE_DEFINITIONS(sundials_sunlinsolpsblas_static PUBLIC ${PSBCDEFINES} ${PSBDEFINES})

 
  #TARGET_INCLUDE_DIRECTORIES(sundials_sunlinsolpsblas_static PUBLIC ${AMG_INCLUDE} ${AMG_MODULES})
  TARGET_INCLUDE_DIRECTORIES(sundials_sunlinsolpsblas_static PUBLIC ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_BINARY_DIR}/include ${PSBLAS_INCLUDE} ${PSBLAS_MODULES} ${AMG_INCLUDE} ${AMG_MODULES})
  TARGET_INCLUDE_DIRECTORIES(sundials_sunlinsolpsblas_static PUBLIC ${PSCTOOLKIT_DIR} ${AMG_DIR} ${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)
  #TARGET_LINK_LIBRARIES(sundials_sunlinsolpsblas_static ${LINKED_LIBRARIES})
    TARGET_LINK_LIBRARIES(sundials_sunlinsolpsblas_static ${LINKED_LIBRARIES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolpsblas_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpsblas CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolpsblas_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
  MESSAGE(STATUS "Added to SUNLINSOL_PSBLAS static libs --- ${PSBLAS_INCLUDE} ${PSBLAS_MODULES} ${AMG_INCLUDE} ${AMG_MODULES};")
  MESSAGE(STATUS "Added to SUNLINSOL_PSBLAS static libs ------  ${PSBCDEFINES} ${PSBFDEFINES}) -- ${PSBLAS_LIBRARY_DIR}; --- ${PSBLDLIBS} *** ${PSBLAS_LIBS}; Last dirs ${PSCTOOLKIT_DIR} ${AMG_DIR};")
  
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOL_PSBLAS library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOL_PSBLAS library
IF(BUILD_SHARED_LIBS AND FALSE)
  ADD_LIBRARY(sundials_sunlinsolpsblas_shared SHARED ${sunlinsolpsblas_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_sunlinsolpsblas_shared m)
  ENDIF()

  TARGET_COMPILE_DEFINITIONS(sundials_sunlinsolpsblas_shared PUBLIC ${PSBCDEFINES} ${PSBFDEFINES})
  TARGET_LINK_LIBRARIES(sundials_sunlinsolpsblas_shared ${LINKED_LIBRARIES})
  TARGET_INCLUDE_DIRECTORIES(sundials_sunlinsolpsblas_shared PUBLIC ${PSBLAS_INCLUDE} ${PSBLAS_MODULES})

  SET_TARGET_PROPERTIES(sundials_sunlinsolpsblas_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpsblas CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolpsblas_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolpsblas_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS AND FALSE)

# Install the SUNLINSOL_PSBLAS header files
INSTALL(FILES ${sunlinsolpsblas_HEADERS} DESTINATION include/sunlinsol)

#
MESSAGE(STATUS "Added SUNLINSOL_PSBLAS module")
