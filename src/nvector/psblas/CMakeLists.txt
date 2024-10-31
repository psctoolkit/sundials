# ---------------------------------------------------------------
# Programmer:  Fabio Durastante @ IAC-CNR
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
# CMakeLists.txt file for the PSBLAS NVECTOR library
# ---------------------------------------------------------------

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_PSBLAS\n\")")


set(sundials_SOURCE_DIR ${PROJECT_SOURCE_DIR})

#EXECUTE_PROCESS(COMMAND ./make2cmakeset.sh ${PSBLAS_INCLUDE_DIR} ${sundials_SOURCE_DIR}/src/nvector/psblas WORKING_DIRECTORY ${sundials_SOURCE_DIR}/src/nvector/psblas OUTPUT_QUIET ERROR_QUIET)
INCLUDE(${sundials_SOURCE_DIR}/src/nvector/psblas/makeincinputcmake)

  message(STATUS "Checking for PSCTOOLKIT ${sundials_SOURCE_DIR};")

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable nvecpsblas_SOURCES with the sources for the nvecpsblas lib
SET(nvecpsblas_SOURCES nvector_psblas.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the nvecpsblas library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecpsblas_HEADERS with the exported nvecpsblas header files
SET(nvecpsblas_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_psblas.h
  )

  message(STATUS "Checking for shared incl${nvecpsblas_HEADERS};")

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY and PSBLAS flags
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)



file(GLOB PSBLAS_LIBS1 "/usr/local/psblas/lib/*.a")  # Adjust the pattern if needed

# PSBLAS LINK DIRECTORIES
#SET(PSBLAS_LIBS -lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base -L/usr/local/psblas/lib -lpsb_cbind)
LINK_DIRECTORIES(${PSBLAS_LIBRARY_DIR})
LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)
LINK_DIRECTORIES(${PSCTOOLKIT_DIR} ${AMG_DIR})
LINK_LIBRARIES(${PSBLDLIBS} ${PSBLAS_LIBS})

#LINK_DIRECTORIES(${PSBLAS_LIBRARY_DIR})
#LINK_DIRECTORIES(${PSCTOOLKIT_DIR}/lib)
LINK_LIBRARIES(${PSBLDLIBS} ${PSBLAS_LIBS} -lgfortran)
# Rules for building and installing the static library:
#  - Add the build target for the nvecpsblas library
#  - Set the library name and make sure it is not deleted
#  - Install the nvecpsblas library
IF(BUILD_STATIC_LIBS)
  MESSAGE(STATUS "Found STATIC --")
  ADD_LIBRARY(sundials_nvecpsblas_static
              STATIC ${nvecpsblas_SOURCES} ${shared_SOURCES})
  TARGET_COMPILE_DEFINITIONS(sundials_nvecpsblas_static PUBLIC ${PSBCDEFINES} ${PSBFDEFINES})
  TARGET_LINK_LIBRARIES(sundials_nvecpsblas_static ${LINKED_LIBRARIES} ${PSBLAS_LIBS1})
  TARGET_INCLUDE_DIRECTORIES(sundials_nvecpsblas_static PUBLIC ${sundials_SOURCE_DIR}/include/ ${PROJECT_BINARY_DIR}/include ${PSBLAS_INCLUDE} ${PSBLAS_MODULES} ${AMG_INCLUDE} ${AMG_MODULES} ${PSCTOOLKIT_DIR} ${AMG_DIR} ${PSCTOOLKIT_DIR}/lib ${AMG_DIR}/lib)
  #TARGET_INCLUDE_DIRECTORIES(sundials_nvecpsblas_static PUBLIC ${AMG_INCLUDE} ${AMG_MODULES})

  SET_TARGET_PROPERTIES(sundials_nvecpsblas_static
    PROPERTIES OUTPUT_NAME sundials_nvecpsblas CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecpsblas_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  MESSAGE(STATUS "Found STATIC ${CMAKE_INSTALL_LIBDIR};")
  
  MESSAGE(STATUS "Added to NVECTOR_PSBLAS static libs ------  ${PSBCDEFINES} ${PSBFDEFINES}) -- ${PSBLAS_LIBRARY_DIR}; --- ${PSBLDLIBS} ${PSBLAS_LIBS}; __ ${LINKED_LIBRARIES};")
  
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the nvecpsblas library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the nvecpsblas library
IF(BUILD_SHARED_LIBS AND FALSE)
  MESSAGE(STATUS "Found SHARED --")
  ADD_LIBRARY(sundials_nvecpsblas_shared
              SHARED ${nvecpsblas_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecpsblas_shared m)
  ENDIF()

  TARGET_COMPILE_DEFINITIONS(sundials_nvecpsblas_shared PUBLIC ${PSBCDEFINES} ${PSBFDEFINES})
  TARGET_LINK_LIBRARIES(sundials_nvecpsblas_shared ${LINKED_LIBRARIES})
  TARGET_INCLUDE_DIRECTORIES(sundials_nvecpsblas_shared PUBLIC ${PSBLAS_INCLUDE} ${PSBLAS_MODULES})

  SET_TARGET_PROPERTIES(sundials_nvecpsblas_shared
    PROPERTIES OUTPUT_NAME sundials_nvecpsblas CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecpsblas_shared
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecpsblas_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS AND FALSE)

# Install the nvecpsblas header files
INSTALL(FILES ${nvecpsblas_HEADERS} DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_PSBLAS module")
