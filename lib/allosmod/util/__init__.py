import contextlib
import tempfile
import shutil

def read_templates(template_file):
    """Read list of templates from the given file. The AllosMod list format
       is one template per file; each line lists the PDB file, chain,
       starting residue, and ending residue, separated by spaces."""
    with open(template_file) as fh:
        return [line.rstrip('\r\n').split()[0] for line in fh]

@contextlib.contextmanager
def temporary_directory():
    """Make a temporary directory"""
    tempd = tempfile.mkdtemp()
    yield tempd
    shutil.rmtree(tempd)

def get_modeller_environ(opts, **keys):
    """Get a new Modeller environment, and set up logging."""
    import modeller
    e = modeller.environ(**keys)
    if opts.verbose:
        modeller.log.verbose()
    else:
        modeller.log.none()
    return e
