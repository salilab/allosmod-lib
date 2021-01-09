"""Generate Modeller distance restraints."""

from __future__ import print_function, absolute_import, division
import optparse
import re
import collections


Restraint = collections.namedtuple('Restraint', ['distance', 'stddev',
                                                 'resind1', 'resind2'])


def get_restraints(dat_file, restr_type):
    """Yield Restraint objects of the given type from the .dat file"""
    with open(dat_file) as fh:
        for line in fh:
            s = line.rstrip('\n\r').split()
            if len(s) == 4 and s[0] == restr_type:
                # Remove any trailing comma, then split on commas
                resind = re.sub(r',\s*$', '', s[3]).split(',')
                for i in range(0, len(resind) - 1, 2):
                    yield Restraint(float(s[1]), float(s[2]),
                                    resind[i], resind[i + 1])


def get_atom_indexes(pdb_file):
    """Get a mapping from residue to atom indexes"""
    atind = 0
    atinds = {}
    attyps = {}
    pref_attyp = {'CA': 1, 'N': 2, 'P': 3, 'C': 4, 'O': 5}
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atind += 1
                attyp = line[12:16].strip()
                resind = line[22:27].strip()
                new_pref = pref_attyp.get(attyp, 9999)
                old_pref = pref_attyp.get(attyps.get(resind, None), 9999)
                if resind not in attyps or new_pref < old_pref:
                    attyps[resind] = attyp
                    atinds[resind] = atind
    return atinds


def get_add_restraint(dat_file, pdb_file, restr_type):
    modeller_form = {'HARM': 3, 'UPBD': 2, 'LOBD': 1}[restr_type]
    atinds = get_atom_indexes(pdb_file)
    for r in get_restraints(dat_file, restr_type):
        if r.resind1 in atinds and r.resind2 in atinds:
            print("R  %3d   1   1  27   2   2   1 %5d %5d     %8.4f  %8.4f"
                  % (modeller_form, atinds[r.resind1], atinds[r.resind2],
                     r.distance, r.stddev))


def parse_args():
    usage = """%prog [opts] <input.dat>
                                 <PDB file> <restraint type>

Generate Modeller distance restraints.

Any residue-residue distance restraints specified in <input.dat> of
<restraint type> are converted to Modeller restraint format and printed out.

Restraints in <input.dat> have a simple format; each line is of the form:
RESTRAINT_TYPE DIST STDDEV INDEX1,INDEX2

where RESTRAINT_TYPE can be HARM (harmonic), UPBD (upper bound),
      or LOBD (lower bound);
DIST is the desired distance (in angstroms);
STDDEV is the standard deviation (in angstroms);
INDEX1,INDEX2 are the residue indexes.

Since Modeller restraints are atom-atom, each restraint is mapped onto atoms
with the following preference: CA, N, P, C, O. If the residue has none of those
atom types, the first atom in the residue is used for the restraint.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return (args[0], args[1], args[2])


def main():
    dat_file, pdb_file, restr_type = parse_args()
    get_add_restraint(dat_file, pdb_file, restr_type)


if __name__ == '__main__':
    main()
