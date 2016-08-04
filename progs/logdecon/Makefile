include ${SEPINC}/SEP.top
include ../../makefile

default: ${BINDIR}/Apply_filter.x \
	${BINDIR}/Apply_sparse_decon_filter.x \
	${BINDIR}/Burg_pef.x \
	${BINDIR}/Center.x \
	${BINDIR}/Cliptomax.x \
	${BINDIR}/Compute_wavelet.x \
	${BINDIR}/Convo.x \
	${BINDIR}/Phase_rotation.x \
	${BINDIR}/Rickerdecon.x \
	${BINDIR}/Sparse_decon_bfgs.x \
	${BINDIR}/Sparse_decon_bfgs_guessinit.x \
	${BINDIR}/Sparse_decon_bfgs_noise_model.x \
	${BINDIR}/Taper.x 

aclean: jclean       
	@-${RM} ${BINDIR}/*.x 

include ${SEPINC}/SEP.bottom
