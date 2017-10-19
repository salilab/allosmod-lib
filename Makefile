include Makefile.include

.PHONY: install test pyext clean

pyext:
	${MAKE} -C lib/allosmod/modeller

clean:
	${MAKE} -C lib/allosmod/modeller clean

install: pyext
	${MAKE} -C bin install
	${MAKE} -C lib/allosmod install
	${MAKE} -C data install

test: pyext
	nosetests --processes=8 test
