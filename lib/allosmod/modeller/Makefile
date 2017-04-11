include ../../../Makefile.include

.PHONY: clean install

PY=${PYTHONDIR}/allosmod/modeller

FILES=${PY}/__init__.py ${PY}/forms.py ${PY}/_truncated_gaussian.so

_truncated_gaussian.so: truncated_gaussian_wrap.o truncated_gaussian.o
	${CC} -o $@ ${LDFLAGS} ${LIBS} $^

%.o: %.c
	${CC} -o $@ ${SHFLAGS} ${CFLAGS} -I../../../src -c $<

truncated_gaussian.o: ../../../src/truncated_gaussian.c
	${CC} -o $@ ${SHFLAGS} ${CFLAGS} -I../../../src -c $<

truncated_gaussian_wrap.c: ../../../src/truncated_gaussian.i
	${SWIG} -o $@ -python -noproxy $<

clean:
	rm -f *.o *_wrap.c *.so

install: ${FILES}

${PY}/_truncated_gaussian.so: _truncated_gaussian.so
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install $< $@

${PY}/%.py: %.py
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install -m 0644 $< $@
	python -m py_compile $@
