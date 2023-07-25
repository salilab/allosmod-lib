include ../../../Makefile.include

.PHONY: install
PY=${PYTHONDIR}/allosmod/util

FILES=${PY}/__init__.py ${PY}/align.py

install: ${FILES}

${PY}/%: %
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install -m 0644 $< $@
	${PYTHON} -m py_compile $@
