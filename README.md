[![python coverage](https://salilab.org/coverage/stat/?s=allosmod-lib&t=python)](http://salilab.org/coverage/allosmod-lib/python/)

[![docs](https://readthedocs.org/projects/allosmod/?badge=latest)](http://allosmod.readthedocs.org/)

This is a library of utility functions used by the
[AllosMod](http://salilab.org/allosmod/)
and [AllosMod-FoXS](http://salilab.org/allosmod-foxs/) web services.

Dependencies
============

- [MODELLER](http://salilab.org/modeller/). This library expects to be able
  to directly import the Modeller Python package (via `import modeller`) and to
  run the Modeller binary (e.g. `mod9.15`). If you installed Modeller from the
  `.tar.gz` package, you will need to set the `PYTHONPATH`, `PATH` and
  `LD_LIBRARY_PATH` environment variables to facilitate this.

- [DSSP](http://swift.cmbi.ru.nl/gv/dssp/). It is expected that the `dssp`
  binary is in the `PATH`.

If you will also be using SAXS profiles (e.g. for the AllosMod-FoXS web service)
then you will need:

- [ProFit](http://www.bioinf.org.uk/programs/profit/). The `profit` binary
  needs to be in the `PATH`.

In the Sali lab, running `module load modeller dssp profit` will get all
of these dependencies.


Building
========

Use `scons test` to test the library, and `scons install` to install it.
See `scons -h` for options to control the build.
