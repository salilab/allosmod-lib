include ../Makefile.include

.PHONY: install

FILES=${BINDIR}/allosmod

install: ${FILES}

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@
