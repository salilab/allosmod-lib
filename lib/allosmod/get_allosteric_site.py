"""Get the allosteric site (all residues near the ligand)"""

from __future__ import print_function, absolute_import
import modeller
import os
import allosmod.util
from allosmod.salign0 import salign0
from allosmod.get_inter_contacts import get_inter_contacts
import sys


class AllostericSiteError(Exception):
    pass


def get_fit_filename(pdb):
    if pdb.endswith('.pdb'):
        return pdb[:-4] + '_fit.pdb'
    else:
        return pdb + '_fit.pdb'


class AllostericSiteFinder(object):
    def __init__(self, env, pdb1, ligand, pdb2, rcut):
        self.__allosteric_site = None
        self.env, self.pdb1, self.ligand, self.pdb2, self.rcut \
            = env, pdb1, ligand, pdb2, rcut

    def find(self):
        """Return a Modeller selection corresponding to the allosteric site.
           @raise AllostericSiteError on error."""
        if self.__allosteric_site is None:
            # align PDB2 to PDB1 and superimpose antigen
            try:
                salign0(self.env, self.pdb1, self.pdb2)
            except modeller.ModellerError as err:
                raise AllostericSiteError(
                    "Could not align %s with %s: %s. "
                    "This is usually due to a poor alignment."
                    % (self.pdb2, self.pdb1, str(err)))
            pmfit = get_fit_filename(self.pdb2)

            # determine residues in PDB2 that contact LIG1
            self.__pmfit = modeller.Model(self.env, file=pmfit)
            lig1 = modeller.Model(self.env, file=self.ligand)
            self.__allosteric_site = \
                modeller.Selection(
                    [ri for ri, rj, dist
                     in get_inter_contacts(self.env, self.__pmfit, lig1,
                                           self.rcut)])
            os.unlink(pmfit)
            os.unlink(get_fit_filename(self.pdb1))
            if len(self.__allosteric_site) == 0:
                raise AllostericSiteError("No allosteric site found")
        return self.__allosteric_site

    def write_atom_list(self, atomlist):
        allosteric_site = self.find()
        for n, atom in enumerate(self.__pmfit.atoms):
            atomlist.write('%d ' % (n+1))
            atomlist.write('AS\n' if atom in allosteric_site else 'RS\n')


def parse_args():
    usage = """%prog [opts] <PDB 1> <ligand> <PDB 2> <cutoff>

Get the allosteric site; this is all residues in <PDB 2> within <cutoff>
of <ligand> (after superposition of <PDB 2> onto <PDB 1>).
"""
    parser = allosmod.util.ModellerOptionParser(usage)
    parser.add_option("--output_pdb", default=None, dest="output_pdb",
                      metavar='FILE',
                      help="Output the allosteric site as a PDB file, "
                           "with the given name.")
    parser.add_option("--atom_list", default=None, dest="atom_list",
                      metavar='FILE',
                      help="Make an atom list file, identifying each atom "
                           "in <PDB 2> as in the allosteric site (AS) or the "
                           "regulated site (RS). This is used by "
                           "editrestraints.")

    opts, args = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments")
    return args[0], args[1], args[2], float(args[3]), opts


def main():
    pdb1, ligand, pdb2, rcut, opts = parse_args()
    e = modeller.Environ()
    e.io.hetatm = True
    a = AllostericSiteFinder(e, pdb1, ligand, pdb2, rcut)
    try:
        if opts.output_pdb:
            a.find().write(file=opts.output_pdb)
        if opts.atom_list:
            a.write_atom_list(open(opts.atom_list, 'w'))
    except AllostericSiteError as err:
        print(str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
