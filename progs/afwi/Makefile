include ${SEPINC}/SEP.top
include ../../makefile

default: ${BINDIR}/2DAFWI.x ${BINDIR}/3DAFWI.x ${BINDIR}/TestRegularization.x

${BINDIR}/2DAFWI.x : ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a 

aclean: jclean       
	@-${RM} ${BINDIR}/*.x 

include ${SEPINC}/SEP.bottom
