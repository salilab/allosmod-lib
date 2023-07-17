import tempfile
import os
import sys
import shutil
import subprocess
import contextlib


def set_search_paths(topdir):
    """Set search paths so that we can run binaries and import
       Python modules"""
    os.environ['PATH'] = os.path.join(topdir, 'bin') + ':' + os.environ['PATH']
    os.environ['PYTHONPATH'] = \
        os.path.join(topdir, 'lib') + ':' + os.environ.get('PYTHONPATH', '')
    sys.path.append(os.path.join(topdir, 'lib'))
    return os.path.join(topdir, 'test')


@contextlib.contextmanager
def temporary_directory():
    _tmpdir = tempfile.mkdtemp()
    yield _tmpdir
    shutil.rmtree(_tmpdir, ignore_errors=True)


@contextlib.contextmanager
def mock_method(cls, method_name, replacement=None):
    """Temporarily replace the given method in the given class."""
    def mock(*args, **keys):
        mock.args = args
        mock.keys = keys
    if replacement is None:
        replacement = mock
    old_method = getattr(cls, method_name)
    setattr(cls, method_name, replacement)
    yield mock
    setattr(cls, method_name, old_method)


def check_output(args, stderr=None, retcode=0, input=None, *other, **keys):
    """Run a subprocess and return its output.
       If the return code from the subprocess does not match `retcode`, an
       `OSError` exception is raised.

       Note: this is similar to `subprocess.check_output` but allows for the
       checked-for return code to be non-zero.
    """
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr,
                         stdin=subprocess.PIPE if input else None,
                         *other, **keys)
    stdout, stderr = p.communicate(input)
    if p.returncode != retcode:
        raise OSError("Process %s exited with code %d, output %s"
                      % (" ".join(args), p.returncode, stdout))
    return stdout


if 'coverage' in sys.modules:
    import atexit
    # Collect coverage information from subprocesses
    __site_tmpdir = tempfile.mkdtemp()
    with open(os.path.join(__site_tmpdir, 'sitecustomize.py'), 'w') as fh:
        fh.write("""
import coverage
import atexit

_cov = coverage.coverage(branch=True, data_suffix=True, auto_data=True,
                         data_file='%s/.coverage')
_cov.start()

def _coverage_cleanup(c):
    c.stop()
atexit.register(_coverage_cleanup, _cov)
""" % os.getcwd())

    os.environ['PYTHONPATH'] = \
        __site_tmpdir + ':' + os.environ.get('PYTHONPATH', '')

    def __cleanup(d):
        shutil.rmtree(d, ignore_errors=True)
    atexit.register(__cleanup, __site_tmpdir)
