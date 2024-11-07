#!/bin/bash
# This script reads as input the Make.inc file for the PSBLAS/AMG2P4 library
# and produces a file that can be included from the CMake installer of KINSOL

liblocation=$1
wheretoputfile=$2

echo "The Make.inc files are in "$liblocation

cat <<EOF > $liblocation/Makefile
include $liblocation/psblas/include/Make.inc.psblas
include $liblocation/amg4psblas/include/Make.inc.amg4psblas
include $liblocation/psblas-ext/include/Make.inc.ext
all:
	@echo "" > $wheretoputfile/makeincinputcmake
	@echo "SET(PSBDEFINES "\${AMGFDEFINES}")\n" >> $wheretoputfile/makeincinputcmake
	@echo "SET(PSBCDEFINES "\${AMGCDEFINES}")\n" >> $wheretoputfile/makeincinputcmake
	@echo "SET(PSBLDLIBS "\${AMGLDLIBS}")\n" >> $wheretoputfile/makeincinputcmake
	@echo "SET(PSBLAS_LIBS "\${PSBLAS_LIBS} "-lpsb_cbind)\n" >> $wheretoputfile/makeincinputcmake
	@echo "SET(PSBRSBLDLIBS "\${LIBRSB_LIBS}")\n" >> $wheretoputfile/makeincinputcmake
	@echo "SET(PSBGPULDLIBS "\${SPGPU_LIBS} \${CUDA_LIBS}")\n" >> $wheretoputfile/makeincinputcmake
EOF

echo "Creating input file for CMake : makeinputcmake"

(cd $liblocation ; make all)

(cd $liblocation ; rm Makefile)
