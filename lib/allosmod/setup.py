"""Check inputs and do initial setup."""

import argparse
from configparser import ConfigParser
import allosmod.util
import allosmod.getcofm
import allosmod.config
import collections
import os
import sys


AlignInfo = collections.namedtuple('AlignInfo', ['codes', 'maxres'])


class _FakeSectionHead(object):
    def __init__(self, fp):
        self.fp = fp
        self.sechead = '[main]\n'

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if line == '':
            raise StopIteration
        else:
            return line

    def readline(self):
        if self.sechead:
            try:
                return self.sechead
            finally:
                self.sechead = None
        else:
            line = self.fp.readline()
            # Hide custom restraints from ConfigParser
            if len(line) > 5 and line[:5] in ('HARM ', 'LOBD ', 'UPBD '):
                return '# ' + line
            else:
                return line


class ConfigFile(object):
    """AllosMod config file (input.dat).
       After calling parse(), config can be looked up like a dict, e.g.

       c = ConfigFile(templates)
       c.parse('input.dat')
       print c['delEmax']

       Field names are case insensitive."""

    def __init__(self, templates, glyco=False, maxres=0):
        self._d = {}
        if templates:
            self._first_template = templates[0]
            self.num_templates = len(templates)
        else:
            self._first_template = None
            self.num_templates = 0
        self.glyco = glyco
        self.maxres = maxres

    def _set_defaults(self):
        defaults = {'DELEMAX': 'CALC', 'DEVIATION': 1.0, 'RAS': 1000,
                    'SAMPLING': 'simulation', 'MDTEMP': 300.0,
                    'REPEAT_OPTIMIZATION': 1, 'ATTACH_GAPS': True,
                    'SCRAPP': False, 'BREAK': False, 'SCLBREAK': 0.1,
                    'ZCUTOFF': 3.5, 'CHEMFR': 'cdensity', 'COARSE': False,
                    'LOCALRIGID': False, 'PW': False,
                    'LIGPDB': self._first_template,
                    'ASPDB': self._first_template}
        # set default if key not present or was explicitly set to empty string
        for k, v in defaults.items():
            if self._d.get(k, '') == '':
                self._d[k] = v

    def __getitem__(self, key):
        return self._d[key.upper()]

    def __setitem__(self, key, val):
        self._d[key.upper()] = val

    def get(self, key, default):
        if key in self:
            return self[key]
        else:
            return default

    def __contains__(self, key):
        return key.upper() in self._d

    def keys(self):
        return self._d.keys()

    def write(self, fname):
        """Write configuration to the given file."""
        with open(fname, 'w') as fh:
            for key, val in self._d.items():
                fh.write('%s=%s\n' % (key, str(val)))

    def parse(self, fname):
        """Read in configuration from the given file.
           Yields any parsing errors."""
        cp = ConfigParser()
        with open(fname) as fh:
            cp.read_file(_FakeSectionHead(fh))
        for name, value in cp.items('main'):
            self._d[name.upper()] = value

        def parse_delemax(val):
            upval = val.upper()
            if upval == 'CALC':
                return upval
            else:
                return float(val)

        def parse_mdtemp(val):
            if isinstance(val, float):
                return val
            lowval = val.lower()
            if lowval == 'scan':
                return lowval
            else:
                return float(val)

        class ParseChoice(object):
            def __init__(self, choices):
                self.choices = choices

            def __call__(self, val):
                val = val.lower()
                if val in self.choices:
                    return val
                else:
                    raise ValueError("not one of %s" % ", ".join(self.choices))

        def parse_boolean(val):
            if isinstance(val, bool):
                return val
            val = val.lower()
            true_types = ('1', 'yes', 'true', 'on')
            false_types = ('0', 'no', 'false', 'off')
            if val in true_types:
                return True
            elif val in false_types:
                return False
            else:
                raise ValueError("not one of %s"
                                 % ", ".join(true_types + false_types))

        def float_or_int(val):
            if isinstance(val, (int, float)):
                return val
            if '.' in val:
                return float(val)
            else:
                return int(val)
        converters = {'DELEMAX': parse_delemax, 'NRUNS': int,
                      'DEVIATION': float, 'RAS': float_or_int,
                      'REPEAT_OPTIMIZATION': int, 'MDTEMP': parse_mdtemp,
                      'SAMPLING': ParseChoice(('simulation', 'moderate_cm',
                                               'moderate_am', 'fast_cm')),
                      'COARSE': parse_boolean, 'ATTACH_GAPS': parse_boolean,
                      'SCRAPP': parse_boolean, 'BREAK': parse_boolean,
                      'PW': parse_boolean, 'SCLBREAK': float,
                      'ZCUTOFF': float,
                      'CHEMFR': ParseChoice(('cdensity', 'charge')),
                      'LOCALRIGID': parse_boolean}
        required_fields = ['NRUNS']
        for k in required_fields:
            if k not in self:
                yield "Missing variable in %s: %s" % (fname, k)
        self._set_defaults()
        for k, converter in converters.items():
            if k in self._d:
                try:
                    self._d[k] = converter(self._d[k])
                except ValueError as err:
                    yield "Invalid variable in %s: %s: %s" \
                          % (fname, k, str(err))
        # If doing glycosylation, force sampling and no deviation
        if self.glyco:
            self['SAMPLING'] = 'moderate_cm'
            self['DEVIATION'] = 0.
        # If only one structure, force rAS=1000
        if self.num_templates == 1 and self['SAMPLING'] != 'moderate_cm':
            self['RAS'] = 1000
        # If number of residues > 1500, use coarse landscape
        if self.maxres > 1500:
            self['COARSE'] = True


class ErrorAccumulator(object):
    error = False

    def report(self, msg):
        print(msg)
        self.error = True


class Setup(object):
    align_file = 'align.ali'
    list_file = 'list'
    dat_file = 'input.dat'
    target = 'pm.pdb'
    ligand = 'lig.pdb'

    def __init__(self):
        self.err = ErrorAccumulator()

    def do_setup(self):
        self.check_input_files_exist()
        align_info = self.check_alignment()
        templates = self.check_list_file(align_info)
        self.config = self.check_input_dat(align_info, templates)
        self.make_ligand(templates)
        self.templates = templates

    def with_glyc2(self):
        """Return True iff doing glycosylation option 2"""
        if os.path.exists('allosmod.py'):
            with open('allosmod.py') as fh:
                for line in fh:
                    if 'self.patch(residue_type=' in line:
                        return True

    def get_other_pdb(self):
        """Get the first template that is not ASPDB (if any)"""
        aspdb = self.config['ASPDB']
        for template in self.templates:
            if template != aspdb:
                return template
        return aspdb

    def make_script_file(self, script_file):
        with open(script_file, 'w') as fh_out:
            self.make_script_file_header(fh_out)
            self.substitute_script_file(fh_out)

    def make_script_file_header(self, fh):
        fh.write('#!/bin/bash\n')
        fh.write(
            'TASK=( null \\\n'
            + ''.join('%d \\\n' % i for i in range(self.config['NRUNS']))
            + ')\n')

    def substitute_script_file(self, fh_out):
        subs = {}
        for k in self.config.keys():
            val = self.config[k]
            subs[k] = str(val).lower() if isinstance(val, bool) else str(val)
        subs.update({'LOCAL_SCRATCH': allosmod.config.local_scratch,
                     'GLOBAL_SCRATCH': allosmod.config.global_scratch,
                     'DATADIR': allosmod.config.datadir,
                     'GLYC1': '1' if self.config.glyco else '0',
                     'GLYC2': '1' if self.with_glyc2() else '0',
                     'COARSE': '--coarse' if self.config['COARSE'] else '',
                     'LOCALRIGID': '--locrigid'
                                   if self.config['LOCALRIGID'] else '',
                     'OTHPDB': self.get_other_pdb()})
        template = allosmod.util.get_data_file('qsub.sh.in')
        allosmod.util.subst_file(open(template), fh_out, subs)

    def make_ligand(self, templates):
        """If no ligand provided, use center of mass of first template"""
        if not os.path.exists(self.ligand) and templates \
           and os.path.exists(templates[0]):
            c = allosmod.getcofm.CenterOfMassPDBParser(
                                           filter=allosmod.util.atom_filter)
            cofm = c.get_cofm(open(templates[0]))
            with open(self.ligand, 'w') as fh:
                fh.write("ATOM      1  XX  ALA A   1    %8.3f%8.3f%8.3f  "
                         "1.00 99.99           C\n" % cofm)

    def error(self):
        return self.err.error

    def check_input_dat(self, align_info, templates):
        if not os.path.exists(self.dat_file):
            return
        d = ConfigFile(templates, os.path.exists('glyc.dat'),
                       align_info.maxres if align_info else 0)
        for e in d.parse(self.dat_file):
            self.err.report(e)
        if align_info:
            for varname in 'LIGPDB', 'ASPDB':
                code = d.get(varname, None)
                if code and code not in align_info.codes:
                    self.err.report("Missing %s in %s for file: %s"
                                    % (varname, self.align_file, code))
        return d

    def check_list_file(self, align_info):
        if not os.path.exists(self.list_file):
            return
        files = []
        with open(self.list_file) as fh:
            for line in fh:
                fil = line.rstrip('\r\n')
                files.append(fil)
                if not os.path.exists(fil):
                    self.err.report("Missing file: %s" % fil)
                if align_info is not None and fil not in align_info.codes:
                    self.err.report("Missing sequence in %s "
                                    "for file: %s" % (self.align_file, fil))
        return files

    def check_alignment(self):
        if not os.path.exists(self.align_file):
            return
        p = allosmod.util.PIRFile()
        try:
            seqs = list(p.read(open(self.align_file)))
        except allosmod.util.FileFormatError as exc:
            self.err.report(str(exc))
            return
        for i in range(1, len(seqs)):
            if len(seqs[i].primary) != len(seqs[0].primary):
                self.err.report("Sequences in %s are not properly aligned"
                                % self.align_file)
        maxres = max(len(s.get_residues()) for s in seqs)
        codes = [s.code for s in seqs]
        if self.target not in codes:
            self.err.report("Missing sequence in %s for file: %s"
                            % (self.align_file, self.target))
        return AlignInfo(codes, maxres)

    def check_input_files_exist(self):
        """Make sure all files exist and don't contain Windows-style
           line endings, which can confuse later parts of the pipeline"""
        files = [self.dat_file, self.align_file, self.list_file]
        for f in files:
            if not os.path.exists(f):
                self.err.report("Missing file: %s" % f)
            else:
                allosmod.util.fix_newlines(f)


def parse_args():
    parser = argparse.ArgumentParser(description="""
Check inputs and do initial setup.
If all the inputs are OK, a script file qsub.sh is generated. This is designed
to be run on an SGE cluster, and will generate MODELLER input files for the
AllosMod protocol. If the SCRAPP option is turned on in input.dat, the files
will automatically be run and the outputs deposited on the global scratch disk.
""")
    _ = parser.parse_args()


def main():
    parse_args()
    s = Setup()
    s.do_setup()
    if s.error():
        sys.exit(1)
    else:
        s.make_script_file('qsub.sh')


if __name__ == '__main__':
    main()
