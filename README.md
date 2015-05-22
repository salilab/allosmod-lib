[![python coverage](https://salilab.org/coverage/stat/?s=allosmod-lib&t=python)](http://salilab.org/coverage/allosmod-lib/python/)

This is a library of utility functions used by the
[AllosMod](http://salilab.org/allosmod/)
and [AllosMod-FoXS](http://salilab.org/allosmod-foxs/) web services.

Dependencies
============

- [MODELLER](http://salilab.org/modeller/). This library expects to be able
  to directly import the Modeller Python package (via `import modeller`). If you
  installed Modeller from the `.tar.gz` package, you will need to set the
  `PYTHONPATH` and `LD_LIBRARY_PATH` environment variables to facilitate this.

Building
========

Use `scons test` to test the library, and `scons install` to install it.
See `scons -h` for options to control the build.
