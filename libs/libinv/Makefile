include $(SEPINC)/SEP.top
include ../../makefile

SRCDIR = .
OBJDIR = .

LIBNAME=${LIBDIR}/libAinv.a
SOURCES= $(notdir $(wildcard ${SRCDIR}/*.r90 ${SRCDIR}/*.f90 ${SRCDIR}/*.f ${SRCDIR}/*.c))

lbfgs.o: lbfgs.f
	${F90C} -O${FCOLEVEL} -Bstatic -c  -o lbfgs.o ./lbfgs.f

blas.o: blas.f
	${F90C} -O${FCOLEVEL} -Bstatic -c  -o blas.o ./blas.f

default: all_Lib

all_Lib: ${LIBNAME}
	-${TOUCH} junk.${F90MODSUFFIX} 
	-${CP} *.${F90MODSUFFIX} ${INCDIR}
	@echo "#########################################"
	@echo "Done making ${LIBNAME}"
	@echo "#########################################"

aclean: jclean
	-${TOUCH} ${INCDIR}/junk.${F90MODSUFFIX}
	-${TOUCH} ${INCDIR}/junk.${F90MODSUFFIX}
	-${TOUCH} ${LIBNAME} ${LIBNAME}.bak 
#	-${RM} ${INCDIR}/${SOURCES}.${F90MODSUFFIX}
	-${RM} ${INCDIR}/*.${F90MODSUFFIX}
	-${RM} ${LIBNAME} ${LIBNAME}.bak

include $(SEPINC)/SEP.bottom
