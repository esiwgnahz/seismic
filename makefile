ADIR=../../

UF90INCLUDES= -I/${MKLROOT}/include/intel64/lp64/ -qopenmp
UF90LIBS= -L/u/st/by/aguitton/bins/fftw/lib/ -qopenmp  -lbei -lsepfft -lsep2df90 -lsep3df90 -lsep3d -lsepf90 -lsep -lsupersetf90 -lsuperset  -lfftw3f  -L/${MKLROOT}/lib/intel64 -qopenmp -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -lAdec -lAinv -lAfut -lAwave
UF77LIBS:= ${UF77LIBS} -qopenmp -lbei

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
#	(cd libs/libfut; make aclean)
#	(cd libs/libdec; make aclean)
#	(cd libs/libinv; make aclean)
	(cd libs/libwave; make aclean)
#	(cd progs/logdecon; make aclean)
#	(cd progs/futterman; make aclean)
	(cd progs/modeling; make aclean)
	(cd progs/rtm; make aclean)


all: 
#	make dirstruct
#	(cd libs/libfut; make)
#	(cd libs/libdec; make)
#	(cd libs/libinv; make)
	(cd libs/libwave; make)
#	(cd progs/logdecon; make)
#	(cd progs/futterman; make)
	(cd progs/modeling; make)
	(cd progs/rtm; make)

