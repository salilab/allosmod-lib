"""Generate Modeller inputs to model with glycosylation"""

from __future__ import print_function, absolute_import
import optparse
import allosmod.util
import allosmod.util.align
import allosmod.config
import string
import os


class BondTypeError(Exception):
    """Error raised for an invalid bond type."""
    pass


class InvalidResidueError(Exception):
    """Error raised for an invalid residue number."""
    pass


class NoSugarsError(Exception):
    """Error raised for an empty glyc.dat."""
    pass


def read_template_file(template_file):
    """Read the list of templates from a file and return it."""
    with open(template_file) as fh:
        return [x.rstrip('\r\n') for x in fh]


def make_ini_model(rand, target, templates, align, fh):
    """Make a Modeller script to generate an initial model"""
    svars = {'knowns': ", ".join(repr(x) for x in templates),
             'data': allosmod.config.datadir,
             'target': target, 'align': align, 'rand': rand}
    fh.write("""
from modeller import *
from modeller.scripts import complete_pdb
import allosmod

env =environ(rand_seed=%(rand)d, restyp_lib_file='%(data)s/restyp.dat',
             copy=None)

# Read in HETATM records from template PDBs
env.io.atom_files_directory = ['.', '../atom_files']
env.libs.topology.read(file='%(data)s/top_all_glyco.lib')
env.libs.parameters.read(file='%(data)s/par_all_glyco.lib')
env.io.hetatm = True
env.io.hydrogen = True

a = allosmod.AllosModel(env, deviation=None, alnfile='%(align)s',
                        knowns=[%(knowns)s],
                        sequence='%(target)s')
# Very thorough VTFM optimization:
a.library_schedule = allosmod.MDnone
a.max_var_iterations = 500

# MD optimization:
a.md_level = allosmod.none

# Repeat the whole cycle 1 time and do not stop unless obj.func. > 1E9
a.repeat_optimization = 1
a.max_molpdf = 1e9

a.starting_model = 1
a.ending_model = 1
a.make()
""" % svars)


def make_ini_model_with_sugar(rand, target, templates):
    """Make a Modeller script to generate an initial model with sugar"""
    with open('model_ini.py', 'w') as fh:
        make_ini_model(rand, target, templates, 'align2.ali', fh)


def make_ini_model_protein_only(rand, target, templates):
    """Make a Modeller script to generate an initial model with protein only"""
    with open('model_ini0.py', 'w') as fh:
        make_ini_model(rand, target, templates, 'align.ali', fh)


def make_glyc_models(rand, target, templates, rep_opt):
    """Make a Modeller script to generate glycosylated models"""
    svars = {'knowns': ", ".join(repr(x) for x in templates),
             'data': allosmod.config.datadir,
             'target': target, 'rep_opt': rep_opt, 'rand': rand}
    with open('model_glyc.py', 'w') as fh:
        fh.write("""
from modeller import *
from modeller.scripts import complete_pdb
import allosmod

env =environ(rand_seed=%(rand)d, restyp_lib_file='%(data)s/restyp.dat',
             copy=None)

# Read in HETATM records from template PDBs
env.io.atom_files_directory = ['.', '../atom_files']
env.libs.topology.read(file='%(data)s/top_all_glyco.lib')
env.libs.parameters.read(file='%(data)s/par_all_glyco.lib')
env.io.hetatm = True
env.io.hydrogen = True

a = allosmod.AllosModel(env, csrfile='converted.rsr', deviation=None,
                        alnfile='align2.ali',
                        knowns=[%(knowns)s],
                        sequence='%(target)s')
# Very thorough VTFM optimization:
a.library_schedule = allosmod.MDopt
a.max_var_iterations = 500

# MD optimization:
a.md_level = allosmod.moderate

# Repeat the whole cycle 1 time and do not stop unless obj.func. > 1E9
a.repeat_optimization = %(rep_opt)d
a.max_molpdf = 1e9

a.starting_model = 1
a.ending_model = 1
a.make()
""" % svars)


class Sugar(object):
    _one_letter_map = {'NAG': '1', 'MAN': '2', 'BMA': '3', 'GLB': '4',
                       'FUC': '5', 'NAN': '8', 'NGA': '9'}
    _connect_atom_map = {'NGLA': 'ND2', 'NGLB': 'ND2',
                         'SGPA': 'OG', 'SGPB': 'OG',
                         'TGPA': 'OG1', 'TGPB': 'OG1',
                         '16ab': 'O6', '16fu': 'O6', '14bb': 'O4',
                         '13ab': 'O3', '13bb': 'O3', '12aa': 'O2',
                         '12ba': 'O2', 'sa23': 'O3', 'sa26': 'O6'}

    def __init__(self, monomer, bond_type, attach_res):
        if bond_type not in self._connect_atom_map.keys():
            raise BondTypeError(
                "Invalid O1 bond type %s. Valid types are: %s"
                % (bond_type,
                   ", ".join(sorted(self._connect_atom_map.keys()))))
        if not attach_res.isdigit():
            raise InvalidResidueError(
                "Invalid residue index for sugar "
                "attachment point: %s. Note that residues are numbered "
                "sequentially starting from 1, with no chain ID or "
                "insertion code." % attach_res)
        self.monomer, self.bond_type = monomer, bond_type
        self.attach_res = int(attach_res)
        self.one_letter_code = self._one_letter_map[monomer]

    def get_connect_atom(self):
        """Return the type of the atom this sugar is bonded to"""
        return self._connect_atom_map[self.bond_type]


class SugarChain(list):
    """Store a chain of Sugars"""
    # todo: assert that sugar-sugar bonds are sane
    pass


def read_glyc_file(fname):
    """Parse the glyc.dat file and return it."""
    s = []
    chain = None
    with open(fname) as fh:
        for line in fh:
            line = line.rstrip('\r\n ')
            # Skip empty lines
            if len(line) == 0:
                continue
            monomer, bond_type, attach_res = line.split()
            if bond_type in ('NGLA', 'NGLB', 'SGPA', 'SGPB', 'TGPA', 'TGPB'):
                chain = SugarChain()
                s.append(chain)
            chain.append(Sugar(monomer, bond_type, attach_res))
    if len(s) == 0:
        raise NoSugarsError("You have provided a glyc.dat that contains "
                            "no valid sugars.")
    return s


def count_residues(alnfile, code):
    """Return the number of residues in the given sequence"""
    p = allosmod.util.PIRFile()
    with open(alnfile) as fh:
        for seq in p.read(fh):
            if seq.code == code:
                return len(seq.get_residues())


def get_residue_chains(pdbfile):
    """Get a mapping from residue number to chain ID"""
    res = {}
    p = allosmod.util.PDBParser(allosmod.util.atom_hetatm_filter)
    with open(pdbfile) as fh:
        for line in p.parse(fh):
            res[line[22:27].strip()] = line[21]
    return res


def get_first_unused_chain(chain_for_res):
    """Return a chain ID that's not already used"""
    all_chains = set(chain_for_res.values())
    for chain_id in string.ascii_uppercase:
        if chain_id not in all_chains:
            return chain_id
    raise ValueError("Cannot find an unused chain ID")


class _Connection(object):
    """Represent a connection between residues"""
    def __init__(self, patch_type, residues, chains, atom_type):
        self.patch_type, self.residues = patch_type, residues
        self.chains, self.atom_type = chains, atom_type

    def write_patch(self, fh):
        """Write a Modeller patch for this connection"""
        fh.write("        self.patch(residue_type='%s', "
                 "residues=(self.residues['%s:%s'],\n"
                 % (self.patch_type, self.residues[0], self.chains[0]))
        fh.write("                                                  "
                 "self.residues['%s:%s']))\n"
                 % (self.residues[1], self.chains[1]))

    def write_rest_in(self, fh):
        """Write an entry to get_rest.in for this connection"""
        fh.write("%s %s %s\n" % (self.residues[0], self.atom_type,
                                 self.patch_type))


def _check_attachments(sugar_chains, chain_for_res):
    """Make sure the attachment point for each sugar chain exists."""
    attachments = [str(sc[0].attach_res) for sc in sugar_chains]
    bad_attachments = [r for r in attachments if r not in chain_for_res]
    if bad_attachments:
        raise InvalidResidueError(
            "Sugar attachment point(s) not found in "
            "PDB file: %s. Note that residues are numbered "
            "sequentially starting from 1, with no chain ID or "
            "insertion code." % (", ".join(bad_attachments)))


def add_glycosidic_bonds(target, glycpm, sugar_chains):
    """Add bonds between protein and sugar chains"""
    last_res = count_residues('align.ali', target)
    chain_for_res = get_residue_chains(glycpm)
    _check_attachments(sugar_chains, chain_for_res)
    sugar_chain = get_first_unused_chain(chain_for_res)
    start_res = last_res + 1
    for sc in sugar_chains:
        # Add bond between first monomer and protein
        protein_chain = chain_for_res[str(sc[0].attach_res)]
        yield _Connection(patch_type=sc[0].bond_type,
                          residues=(sc[0].attach_res, start_res),
                          chains=(protein_chain, sugar_chain),
                          atom_type=sc[0].get_connect_atom())
        # Add bonds between monomers in the chain
        for n_s, s in enumerate(sc[1:]):
            yield _Connection(patch_type=s.bond_type,
                              residues=(start_res - 1 + s.attach_res,
                                        start_res + n_s + 1),
                              chains=(sugar_chain, sugar_chain),
                              atom_type=s.get_connect_atom())
        start_res += len(sc)


def make_alignment_with_sugar(target, templates, sugar_chains):
    """Make a Modeller alignment to include sugars in the model."""
    one_letters = []
    for sc in sugar_chains:
        one_letters.extend(s.one_letter_code for s in sc)
    p = allosmod.util.PIRFile()
    with open('align.ali') as fh_in:
        with open('align2.ali', 'w') as fh_out:
            for seq in p.read(fh_in):
                if seq.code == target:
                    seq.range[1][0] = int(seq.range[1][0]) + len(one_letters)
                    seq.primary += '/' + ''.join(one_letters)
                else:
                    seq.primary += '/' + '-' * len(one_letters)
                p.write(fh_out, seq)


def get_gaps(sugar_chains):
    """Get a list of places to insert alignment gaps"""
    gaps = []
    for chain in sugar_chains:
        r = chain[0].attach_res
        gaps.append((max(1, r - 2), r + 2))
    return sorted(gaps)


def remove_overlaps(gaps):
    """Combine any overlapping gaps into a single gap"""
    result = []
    current_lb, current_ub = -1, -1
    for lb, ub in gaps:
        if lb > current_ub:
            result.append((lb, ub))
            current_lb, current_ub = lb, ub
        else:
            # overlap, so update the last gap
            # list is sorted, so current_lb <= lb
            current_ub = max(current_ub, ub)
            result[-1] = (current_lb, current_ub)
    return result


def add_alignment_gaps(sugar_chains):
    """Add gaps to the alignment near insertion sites"""
    gaps = remove_overlaps(get_gaps(sugar_chains))
    for gap in gaps:
        with open('align2.ali.tmp', 'w') as fh:
            allosmod.util.align.insert_gap("align2.ali", 0, gap[0]-1,
                                           gap[1]-1, fh)
        os.rename('align2.ali.tmp', 'align2.ali')


def get_pm_glyc(target, template_file, rand, rep_opt, att_gap, glycpm):
    templates = read_template_file(template_file)
    make_ini_model_with_sugar(rand, target, templates)
    make_ini_model_protein_only(rand, target, templates)
    make_glyc_models(rand, target, templates, rep_opt)
    if glycpm == "script":
        return
    sugar_chains = read_glyc_file('glyc.dat')
    with open('allosmod.py', 'a') as py:
        py.write('    # define bonds between sugar monomers\n')
        py.write('    def special_patches(self, aln):\n')
        with open('get_rest.in', 'w') as rest:
            for conn in add_glycosidic_bonds(target, glycpm, sugar_chains):
                conn.write_patch(py)
                conn.write_rest_in(rest)
    make_alignment_with_sugar(target, templates, sugar_chains)
    if att_gap.upper() == 'TRUE':
        add_alignment_gaps(sugar_chains)


def parse_args():
    usage = """%prog [opts] <target> <templates> <rand> <rep_opt>
                 <att_gap> <glycpm>

Generate Modeller scripts and alignments to model with glycosylation.

<target> alignment code of the target (usually pm.pdb)
<template_file> file containing a list of all templates
<rand> random seed
<rep_opt> option to allow repeat optimization
<att_gap> if set to "TRUE", inserts gaps to allow flexibility at
          glycosylation sites
<glycpm> either "script" or the name of a PDB file
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 6:
        parser.error("incorrect number of arguments")
    return args


def main():
    target, template_file, rand, rep_opt, att_gap, glycpm = parse_args()
    get_pm_glyc(target, template_file, int(rand), int(rep_opt),
                att_gap, glycpm)


if __name__ == '__main__':
    main()
