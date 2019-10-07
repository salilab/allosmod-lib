include Makefile.include

.PHONY: install test pyext clean all bin

all: pyext bin

bin:
	${MAKE} -C bin

pyext:
	${MAKE} -C lib/allosmod/modeller

clean:
	${MAKE} -C lib/allosmod/modeller clean
	${MAKE} -C bin clean

install: pyext bin
	${MAKE} -C bin install
	${MAKE} -C lib/allosmod install
	${MAKE} -C data install

test: pyext bin
	nosetests --processes=8 test
