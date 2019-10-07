include ../Makefile.include

.PHONY: install clean

allosmod: allosmod.in
	FULL_PYTHON=`which ${PYTHON}`; sed -e "s^#!.*^#!$$FULL_PYTHON^" $< > $@ && chmod a+x $@

FILES=${BINDIR}/allosmod

install: ${FILES}

clean:
	rm -f allosmod

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@
