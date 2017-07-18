include ${SEPINC}/SEP.top
include ../../makefile

default: ${BINDIR}/3DBORNMOD.x ${BINDIR}/Acoustic_modeling.x ${BINDIR}/Acoustic_modeling_sep.x #${BINDIR}/Acoustic_Born_modeling.x ${BINDIR}/Acoustic_modeling_sep.x ${BINDIR}/Acoustic_Born_modeling_opt.x ${BINDIR}/Acoustic_Born_modeling_opt_memory.x

${BINDIR}/3DBORNMOD.x:  ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a
${BINDIR}/Acoustic_modeling.x: ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a
${BINDIR}/Acoustic_modeling_sep.x: ../../lib/LINUX/libAinv.a ../../lib/LINUX/libAwave.a
aclean: jclean       
	@-${RM} ${BINDIR}/*.x 

include ${SEPINC}/SEP.bottom
