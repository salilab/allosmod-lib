Installation
************

In the Sali Lab
===============

If you are working in the Sali lab, you don't need to build and install
AllosMod - it is already set up for you as a module. Just run
``module load allosmod`` to load it.

Dependencies
============

* `Python <https://www.python.org>`_ 3.6 or later.

* `MODELLER <https://salilab.org/modeller/>`_. This library expects to be able
  to directly import the Modeller Python package (via ``import modeller``)
  and to run the Modeller binary (e.g. :command:`mod9.15`). If you installed
  Modeller from the ``.tar.gz`` package, you will need to set the
  ``PYTHONPATH``, ``PATH`` and ``LD_LIBRARY_PATH`` environment variables
  to facilitate this.

* `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_. It is expected that the
  :command:`mkdssp` binary is in the ``PATH``.

* `nose <https://nose.readthedocs.io/en/latest/>`_ is also needed to run the
  test suite (recommended but not essential).

If you will also be using SAXS profiles (e.g. for the AllosMod-FoXS web service)
then you will need:

* `ProFit <http://www.bioinf.org.uk/programs/profit/>`_. The
  :command:`profit` binary needs to be in the ``PATH``.

In the Sali lab, running ``module load modeller dssp profit`` will get all
of these dependencies.


Building
========

Use ``make test`` to test the library, and ``make install`` to install it.
In most cases you will need to tell ``make`` where to install (if running on
a Linux cluster, AllosMod will need to be installed on a network-accessible
filesystem) and where any temporary/scratch disks usable by AllosMod are.
(AllosMod needs a 'local' scratch disk, accessible by individual jobs on each
node, as well as a 'global' scratch disk, on a network-accessible filesystem.)
Do this with something like
``make PREFIX=/shared/allosmod GLOBAL_SCRATCH=/scratch install``. See
``Makefile.include`` for all make variables that can be configured.
