import contextlib
import optparse
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

class ModellerOptionParser(optparse.OptionParser):
    """Add options to control the amount of Modeller logging to optparse"""

    def __init__(self, *args, **keys):
        optparse.OptionParser.__init__(self, *args, **keys)
        self.add_option("-v", "--verbose", action="store_true", dest="verbose",
                        help="verbose output")

    def parse_args(self):
        import modeller
        opts, args = optparse.OptionParser.parse_args(self)
        if opts.verbose:
            modeller.log.verbose()
        else:
            modeller.log.none()
        return opts, args
