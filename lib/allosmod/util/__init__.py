from __future__ import print_function
import contextlib
import optparse
import tempfile
import subprocess
import shutil
import os
import re


def check_output(args, stderr=None, retcode=0, input=None, *other, **keys):
    """Run a subprocess and return its output.
       If the return code from the subprocess does not match `retcode`, an
       `OSError` exception is raised.

       Note: this is similar to `subprocess.check_output` but that requires
       Python 2.7.
    """
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr,
                         stdin=subprocess.PIPE if input else None,
                         *other, **keys)
    stdout, stderr = p.communicate(input)
    if p.returncode != retcode:
        raise OSError("Process %s exited with code %d, output %s"
                      % (" ".join(args), p.returncode, stdout))
    return stdout


def subst_file(fh_in, fh_out, subs):
    """Make file fh_out by substituting variables in fh_in.
       Substitutions are of the form @FOO@, and will be replaced by
       the value of subs['FOO']. @@ is replaced by @.
    """
    r = re.compile(r'@(?P<key>\w*?)@')

    def repl(match):
        key = match.group('key')
        if not key:
            return "@"    # Replace @@ with @
        elif key in subs:
            return subs[key]
        else:
            raise ValueError("Unknown substitution %s" % key)
    for line in fh_in:
        fh_out.write(r.sub(repl, line))


def get_data_file(fname):
    """Return the full path to a file in the data directory"""
    import allosmod.config
    return os.path.join(allosmod.config.datadir, fname)


def fix_newlines(fname):
    """Remove any \r characters from the file, in place"""
    with open(fname) as fh:
        contents = fh.read()
    if '\r' in contents:
        with open(fname, 'w') as fh:
            fh.write(contents.replace('\r', ''))


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


class FileFormatError(Exception):
    pass


class Sequence(object):
    """Representation of a single amino acid sequence"""

    def __init__(self):
        self.prottyp = 'sequence'
        self.range = [['', ''], ['', '']]
        self.atom_file = self.name = self.source = ''
        self.resolution = self.rfactor = ''

    def get_residues(self):
        """Get residues in this sequence (same as 'primary' but without
           gaps or chain breaks)"""
        return self.primary.replace('-', '').replace('/', '')


class PIRFile(object):
    """Representation of a PIR-format file"""

    def _parse_pir_header(self, num, line, seq):
        seq.primary = ''
        spl = line.rstrip().split(':')
        if len(spl) != 10:
            raise FileFormatError(
                      "Invalid PIR header at line %d (expecting 10 fields "
                      "split by colons): %s" % (num + 1, line))
        (seq.prottyp, seq.atom_file, seq.range[0][0], seq.range[0][1],
         seq.range[1][0], seq.range[1][1], seq.name, seq.source,
         seq.resolution, seq.rfactor) = spl
        if seq.prottyp == '':
            seq.prottyp = 'sequence'

    def read(self, fh):
        """Read sequences from the given stream in PIR format. A list of
           the sequences is returned, as :class:`Sequence` objects."""
        seq = None
        terminator = re.compile(r'\*\s*$')
        for (num, line) in enumerate(fh):
            if line.startswith('C;') or line.startswith('R;'):
                # Skip comment lines
                continue
            elif line.startswith('>P1;'):
                if seq:
                    raise FileFormatError(
                        "PIR sequence without terminating * at line %d: %s"
                        % (num + 1, line))
                seq = Sequence()
                seq.primary = None
                seq.code = line[4:].strip()
            elif seq and seq.primary is None:
                self._parse_pir_header(num, line, seq)
            else:
                line = line.rstrip()
                if line:
                    if seq is None:
                        raise FileFormatError(
                             "PIR sequence found without a preceding header "
                             "at line %d: %s" % (num + 1, line))
                    (line, count) = terminator.subn("", line)
                    seq.primary += line
                    # See if this was the last line in the sequence
                    if count == 1:
                        yield seq
                        seq = None
        if seq:
            raise FileFormatError(
                     "PIR sequence without terminating * at end of file")

    def write(self, fh, seq, width=70):
        """Write a single :class:`Sequence` object to the given stream in
           PIR format."""
        print(">P1;" + seq.code, file=fh)
        start, end = seq.range
        print(":".join(str(x) for x in [seq.prottyp, seq.atom_file,
                                        start[0], start[1], end[0],
                                        end[1], seq.name, seq.source,
                                        seq.resolution, seq.rfactor]), file=fh)
        for pos in range(0, len(seq.primary), width):
            print(seq.primary[pos:pos+width], file=fh)
        print('*', file=fh)


class PDBParser(object):
    def __init__(self, filter):
        self.filter = filter

    def parse(self, fh):
        for line in fh:
            if self.filter(line):
                yield line


def atom_filter(line):
    return line.startswith('ATOM')


def atom_hetatm_filter(line):
    return line.startswith('ATOM') or line.startswith('HETATM')
