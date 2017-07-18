include ${SEPINC}/SEP.top
include ../../makefile

#default: ${BINDIR}/Acoustic_RTM_sep.x ${BINDIR}/Acoustic_RTM.x

default: ${BINDIR}/3DRTM.x ${BINDIR}/Acoustic_RTM_opt_memory.x ${BINDIR}/LaplaceFilter.x #${BINDIR}/Acoustic_RTM_opt.x ${BINDIR}/Acoustic_RTM_opt_memory.x ${BINDIR}/LaplaceFilter.x

${BINDIR}/3DRTM.x: ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a
${BINDIR}/Acoustic_RTM_opt_memory.x: ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a
aclean: jclean       
	@-${RM} ${BINDIR}/*.x

include ${SEPINC}/SEP.bottom
