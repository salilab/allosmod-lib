"""Get all residues involved in charge contacts."""

import modeller
import math
import allosmod.util


charged_residues = dict.fromkeys(('ARG', 'HIS', 'LYS', 'HSD', 'HSE',
                                  'HSP', 'HID', 'HIE', 'HIP', 'ASP', 'GLU'))


def get_restrained_atoms(mdl, rsr_file):
    for line in rsr_file:
        spl = line.split()
        if len(spl) > 10 and spl[5] == '2' and spl[4] != '1':
            yield mdl.atoms[int(spl[8])-1], mdl.atoms[int(spl[9])-1]


def charged_ca_pair(a1, a2):
    return a1.name == 'CA' and a2.name == 'CA' \
           and a1.residue.pdb_name in charged_residues \
           and a2.residue.pdb_name in charged_residues


class ChargedContactFinder:
    def __init__(self, env, rsr_file, pdb_file):
        self.env, self.rsr_file, self.pdb_file = env, rsr_file, pdb_file

    def find(self):
        m = self._m = modeller.Model(self.env, file=self.pdb_file)
        charge = [0] * len(m.residues)
        total = [0] * len(m.residues)
        for a1, a2 in get_restrained_atoms(m, open(self.rsr_file)):
            r1 = a1.residue.index - 1
            r2 = a2.residue.index - 1
            total[r1] += 1
            total[r2] += 1
            if charged_ca_pair(a1, a2):
                charge[r1] += 1
                charge[r2] += 1
        self._total = total
        self._contacts = []
        for n, z in enumerate(zip(charge, total)):
            c, t = z
            self._contacts.append((n+1, c if t > 143 else 0))

    def print_contacts(self, fh):
        for n, c in self._contacts:
            print("%d %d" % (n, c), file=fh)

    def print_all_buried(self, sclbreak, fh):
        for n, c in enumerate(self._total):
            if c > 143 and self._m.residues[n].pdb_name in charged_residues:
                print("%d %s" % (n+1, sclbreak), file=fh)

    def print_cdensity(self, cutoff, sclbreak, fh):
        xavg = sum(x[1] for x in self._contacts) / len(self._contacts)
        xsd = math.sqrt(sum((x[1]-xavg)**2
                        for x in self._contacts) / len(self._contacts))
        if xsd > 0.:
            for n, c in self._contacts:
                v = (c - xavg) / xsd
                if abs(v) > cutoff:
                    print("%d %s" % (n, sclbreak), file=fh)
        elif cutoff < 0:
            for n, c in self._contacts:
                print("%d %s" % (n, sclbreak), file=fh)


def parse_args():
    usage = """%prog [opts] <restraint file> <PDB file> <sclbreak>

Output all residues from <PDB file> that are involved in charge contacts
(as specified by <restraint file>). Charged residues will also be output
to break.dat (with the given <sclbreak> scale factor).

Note that residue indices (starting from 1) are used, not PDB residue numbers.
"""
    parser = allosmod.util.ModellerOptionParser(usage)
    parser.add_option("--cdensity_cutoff", type=float, default=None,
                      dest="cdensity_cutoff", metavar='FLOAT',
                      help="If given, only residues that are buried and with "
                           "a high charge density (> cutoff standard "
                           "deviations above the mean) will be output; "
                           "otherwise, all charged residues will be output.")

    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return args[0], args[1], args[2], opts


def main():
    rsr_file, pdb_file, sclbreak, opts = parse_args()
    e = modeller.Environ()
    e.io.hetatm = True
    a = ChargedContactFinder(e, rsr_file, pdb_file)
    a.find()
    a.print_contacts(open('contpres.dat', 'w'))
    if opts.cdensity_cutoff is not None:
        a.print_cdensity(opts.cdensity_cutoff, sclbreak,
                         open('break.dat', 'a'))
    else:
        a.print_all_buried(sclbreak, open('break.dat', 'a'))


if __name__ == '__main__':
    main()
