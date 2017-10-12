include $(SEPINC)/SEP.top
include ../../makefile

SRCDIR = .
OBJDIR = .

LIBNAME=${LIBDIR}/libAnlinv.a
SOURCES= $(notdir $(wildcard ${SRCDIR}/*.r90 ${SRCDIR}/*.f90 ${SRCDIR}/*.f ${SRCDIR}/*.c))

nlinv_lbfgs.o: nlinv_lbfgs.f
	${F90C} -O${FCOLEVEL} -Bstatic -c  -o nlinv_lbfgs.o ./nlinv_lbfgs.f

nlinv_blas.o: nlinv_blas.f
	${F90C} -O${FCOLEVEL} -Bstatic -c  -o nlinv_blas.o ./nlinv_blas.f

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
