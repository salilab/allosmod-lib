include ../Makefile.include

.PHONY: install

FILES=${DATADIR}/BMA.rsr ${DATADIR}/FUC.rsr ${DATADIR}/getavgpdb2.sh \
      ${DATADIR}/get_glyc_restraint.sh ${DATADIR}/get_loopadjres.sh \
      ${DATADIR}/getmax ${DATADIR}/getmin ${DATADIR}/get_rest.sh \
      ${DATADIR}/getsegpdb ${DATADIR}/GLB.rsr ${DATADIR}/MAN.rsr \
      ${DATADIR}/NAG.rsr ${DATADIR}/NAN.rsr ${DATADIR}/NGA.rsr \
      ${DATADIR}/par_all_glyco.lib ${DATADIR}/pm_fast_cm.py.in \
      ${DATADIR}/pm_moderate_am.py.in ${DATADIR}/pm_moderate_cm.py.in \
      ${DATADIR}/pm_script.py.in ${DATADIR}/qsub.sh.in \
      ${DATADIR}/randomize_list_2.pl ${DATADIR}/README_user \
      ${DATADIR}/restyp.dat ${DATADIR}/salign_sub.sh \
      ${DATADIR}/top_all_glyco.lib

install: ${FILES}

${DATADIR}/%: %
	@if [ ! -d ${DATADIR} ]; then mkdir -p ${DATADIR}; fi
	cp $< $@
