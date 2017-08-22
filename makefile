ADIR=../../

UF90INCLUDES= -I/${MKLROOT}/include/intel64/lp64/ -qopenmp
UF90LIBS= -L/u/st/by/aguitton/bins/fftw/lib/ -qopenmp  -lbei -lsepfft -lsep2df90 -lsep3df90 -lsep3d -lsepf90 -lsep -lsupersetf90 -lsuperset  -lfftw3f  -L/${MKLROOT}/lib/intel64 -qopenmp -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -lsepgeef90 -lAdec -lAfut -lAwave -lAflt -lAinv -lsepgeef90
UF77LIBS:= ${UF77LIBS} -qopenmp -lbei

###########################################
# Add -g 2 lines below if debugger needed
OF90FLAGS:= -axCOMMON-AVX512 -qopenmp -fPIC -FR  -D_GNU_SOURCE -I${SEPINCDIR} -c
OF77FLAGS:= -axCOMMON-AVX512 -qopenmp -fPIC -FR                               -c 
###########################################

LIBDIR  = ${ADIR}/lib/$(MTYPE)
INCDIR  = ${ADIR}/inc/$(MTYPE)
BINDIR  = ${ADIR}/bin/${MTYPE}

UF90INCLUDES := ${UF90INCLUDES} ${F90INCFLAG}${INCDIR}
UF90LIBDIRS  := ${UF90LIBDIRS} ${LIBDIR}
F90MODSUFFIX = mod

dirstruct:
	\rm -rf bin lib inc
	mkdir bin
	mkdir bin/LINUX
	mkdir lib
	mkdir lib/LINUX
	mkdir inc
	mkdir inc/LINUX

clean: 
	(cd libs/libfut; make aclean)
	(cd libs/libdec; make aclean)
	(cd libs/libinv; make aclean)
	(cd libs/libflt; make aclean)
	(cd libs/libwave; make aclean)
	(cd progs/logdecon; make aclean)
	(cd progs/futterman; make aclean)
	(cd progs/modeling; make aclean)
	(cd progs/rtm; make aclean)
	(cd progs/afwi; make aclean)
	(cd progs/match; make aclean)
	(cd progs/procs; make aclean)
	(cd progs/inversion; make aclean)


all: 
#	make dirstruct
	(cd libs/libfut; make)
	(cd libs/libdec; make)
	(cd libs/libinv; make)
	(cd libs/libflt; make)
	(cd libs/libwave; make)
#	(cd progs/logdecon; make)
#	(cd progs/match; make)
#	(cd progs/futterman; make)
	(cd progs/modeling; make )
	(cd progs/rtm; make)
	(cd progs/afwi; make)
	(cd progs/procs; make)
	(cd progs/inversion; make)

listfiles='ls  */*/*90'

listdo:
	for files in ${listfiles}; do\
		echo $$files ;\
		sed -i '1i ! ' $$files ;\
		sed -i '1i ! -----------------------------------------------' $$files ;\
		sed -i '1i ! Copyright (c) 2016-2017 Bellevue Geophysics LLC' $$files ;\
		sed -i '1i ! -----------------------------------------------' $$files ;\
		sed -i '1i ! ' $$files ;\
	done
