include ../../Makefile.include

SUBDIRS=util modeller
.PHONY: install ${SUBDIRS}

PY=${PYTHONDIR}/allosmod

FILES=${PY}/__init__.py ${PY}/bin_data.py ${PY}/contpres.py \
      ${PY}/count_alignments.py ${PY}/edit_restraints.py \
      ${PY}/get_add_restraint.py ${PY}/get_allosteric_site.py \
      ${PY}/get_auto_align.py ${PY}/getavgpdb.py ${PY}/getcofm.py \
      ${PY}/get_contacts.py ${PY}/get_glyc_restraint.py \
      ${PY}/get_inter_contacts.py ${PY}/get_loopadjres.py \
      ${PY}/get_pm_glyc.py ${PY}/get_pm_initialstruct.py ${PY}/get_q_ca.py \
      ${PY}/get_qiavg_ca.py ${PY}/get_qmatrix.py ${PY}/get_rest.py \
      ${PY}/getrofg.py ${PY}/get_ss.py ${PY}/__init__.py \
      ${PY}/make_mod_inputs.py ${PY}/make_pm_script.py ${PY}/min_rmsd.py \
      ${PY}/pdb2ali.py ${PY}/pdb_fix_res.py ${PY}/rotatepdb.py \
      ${PY}/salign0.py ${PY}/setchain.py ${PY}/setup.py ${PY}/spline.py \
      ${PY}/translatepdb.py ${PY}/config/__init__.py

install: ${FILES} ${SUBDIRS}

${SUBDIRS}:
	${MAKE} -C $@ install

${PY}/config/__init__.py::
	@if [ ! -d ${PY}/config ]; then mkdir -p ${PY}/config; fi
	echo "datadir = r'${DATADIR}'" > $@
	echo "local_scratch = r'${LOCAL_SCRATCH}'" >> $@
	echo "global_scratch = r'${GLOBAL_SCRATCH}'" >> $@
	${PYTHON} -m py_compile $@

${PY}/%: %
	@if [ ! -d ${PY} ]; then mkdir -p ${PY}; fi
	install -m 0644 $< $@
	${PYTHON} -m py_compile $@
