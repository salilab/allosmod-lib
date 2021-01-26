"""Make a list of all contacts between two structures."""

from __future__ import print_function, absolute_import
import allosmod.util
from allosmod.get_contacts import _get_average_aa, get_contact_type
from allosmod.get_contacts import get_contact_dist


def get_inter_contacts(env, mdl1, mdl2, rcut):
    av1 = [_get_average_aa(r) for r in mdl1.residues]
    av2 = [_get_average_aa(r) for r in mdl2.residues]

    rcut2 = rcut * rcut
    for i in av1:
        for j in av2:
            dist = get_contact_dist(i, j, rcut2)
            if dist is not None:
                yield i.r, j.r, dist


def parse_args():
    usage = """%prog <PDB file 1> <PDB file 2> <cutoff>

Make a list of all residues in contact, i.e. all residues in <PDB file 1>
which are less than <cutoff> angstroms from a residue in <PDB file 2>.

Distance is considered to be between representative centers of each residue.
This center is the mass center of all heavy atoms in the amino acid sidechain,
plus CB (for all residues except GLY) or CA (for GLY). If these atoms don't
exist, the CA, O or N atom (in that order) is used as the center.
"""
    parser = allosmod.util.ModellerOptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return args[0], args[1], float(args[2])


def main():
    import modeller
    file1, file2, rcut = parse_args()
    e = modeller.Environ()
    e.io.hetatm = True
    mdl1 = modeller.Model(e, file=file1)
    mdl2 = modeller.Model(e, file=file2)
    for ri, rj, dist in get_inter_contacts(e, mdl1, mdl2, rcut):
        print("  %6s  %2s  %6s  %2s  %3s  %3s%3d%11.3f  %1d  %1d"
              % (ri.num, ri.chain.name, rj.num, rj.chain.name,
                 ri.pdb_name, rj.pdb_name, get_contact_type(ri, rj),
                 dist, 6, 6))


if __name__ == '__main__':
    main()
