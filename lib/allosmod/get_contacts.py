"""Make a list of all residues in contact."""

from __future__ import print_function, absolute_import
import optparse
import math
import collections


Residue = collections.namedtuple('Residue', ['r', 'average'])


def _get_average(to_average):
    x = sum(a.x for a in to_average)
    y = sum(a.y for a in to_average)
    z = sum(a.z for a in to_average)
    return (x / len(to_average), y / len(to_average),
            z / len(to_average))


def _get_average_hem(r):
    def get_av(names):
        return _get_average([a for a in r.atoms if a.name in names])
    return Residue(r, (get_av(['NA', 'C1A', 'C2A', 'C3A', 'C4A']),
                       get_av(['NB', 'C1B', 'C2B', 'C3B', 'C4B']),
                       get_av(['NC', 'C1C', 'C2C', 'C3C', 'C4C']),
                       get_av(['ND', 'C1D', 'C2D', 'C3D', 'C4D'])))


def _get_average_aa(r):
    main_chain_atom_names = dict.fromkeys(['O', 'N', 'C', 'OT', 'CA', 'CB'])
    to_average = []
    for a in r.atoms:
        if a.name not in main_chain_atom_names \
          or (a.name == 'CA' and r.pdb_name == 'GLY') \
          or (a.name == 'CB' and r.pdb_name != 'GLY'):
            to_average.append(a)
    if to_average:
        return Residue(r, (_get_average(to_average),))
    else:
        for typ in ('CA', 'O', 'N'):
            try:
                a = r.atoms[typ]
                return Residue(r, ((a.x, a.y, a.z),))
            except KeyError:
                pass
    raise ValueError("no average")


def get_average_coordinate(r):
    if r.pdb_name == 'HEM':
        return _get_average_hem(r)
    else:
        return _get_average_aa(r)


def get_contact_type(r1, r2):
    hphob = 1
    cgpol = 2
    restyp = {'ALA': cgpol, 'GLY': cgpol, 'PRO': cgpol, 'SER': cgpol,
              'THR': cgpol, 'ASN': cgpol, 'ASP': cgpol, 'GLN': cgpol,
              'GLU': cgpol, 'ARG': cgpol, 'HIS': cgpol, 'LYS': cgpol,
              'HSD': cgpol, 'CYS': hphob, 'ILE': hphob, 'LEU': hphob,
              'MET': hphob, 'PHE': hphob, 'TRP': hphob, 'TYR': hphob,
              'VAL': hphob}
    c1 = restyp.get(r1.pdb_name, None)
    c2 = restyp.get(r2.pdb_name, None)
    if c1 == hphob and c2 == hphob:
        return 1  # hydrophobic-hydrophobic
    elif c1 == cgpol and c2 == cgpol:
        return 2  # charge/polar-charge/polar
    elif c1 is None or c2 is None:
        return 0  # other
    else:
        return 3  # hydrophobic-charge/polar


def get_contact_dist(r1, r2, rcut2):
    for av1 in r1.average:
        for av2 in r2.average:
            dist = sum((a-b) ** 2 for a, b in zip(av1, av2))
            if dist < rcut2:
                return math.sqrt(dist)


def get_contacts(pdb_file, rcut):
    import modeller

    modeller.log.none()
    e = modeller.Environ()
    e.io.hetatm = True

    m = modeller.Model(e, file=pdb_file)

    rcut2 = rcut * rcut
    av = [get_average_coordinate(r) for r in m.residues]
    for i in range(len(av) - 3):
        for j in range(i + 3, len(av)):
            ri = av[i].r
            rj = av[j].r
            if ri.hetatm and rj.hetatm:
                continue  # do not print het-het contacts
            dist = get_contact_dist(av[i], av[j], rcut2)
            if dist is not None:
                yield ri, rj, dist


def parse_args():
    usage = """%prog <PDB file> <cutoff>

Make a list of all residues in contact, i.e. all residues in <PDB file>
which are less than <cutoff> angstroms from another residue.

Distance is considered to be between representative centers of each residue.
For most residues, this center is the mass center of all heavy atoms in
the amino acid sidechain, plus CB (for all residues except GLY) or CA
(for GLY). If these atoms don't exist, the CA, O or N atom (in that order)
is used as the center.

HEM residues are considered to have four centers, corresponding to the centers
of each of the four pyrrole rings.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return args[0], float(args[1])


def main():
    pdb_file, rcut = parse_args()
    for ri, rj, dist in get_contacts(pdb_file, rcut):
        print("  %6s  %2s  %6s  %2s  %3s  %3s%3d%11.3f  %1d  %1d"
              % (ri.num, ri.chain.name, rj.num, rj.chain.name,
                 ri.pdb_name, rj.pdb_name, get_contact_type(ri, rj),
                 dist, 6, 6))


if __name__ == '__main__':
    main()
