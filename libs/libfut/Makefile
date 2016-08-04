include $(SEPINC)/SEP.top
include ../../makefile

SRCDIR = .
OBJDIR = .

LIBNAME=${LIBDIR}/libAfut.a
SOURCES= $(notdir $(wildcard ${SRCDIR}/*.f90 ${SRCDIR}/*.f ${SRCDIR}/*.c))

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
	-${RM} ${INCDIR}/*.${F90MODSUFFIX}
	-${RM} ${LIBNAME} ${LIBNAME}.bak

include $(SEPINC)/SEP.bottom
