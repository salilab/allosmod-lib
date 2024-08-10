"""Calculate qiavg for CA atoms.

Qi is a distance similarity metric in which a score of 1 implies identity
to a reference structure and a score of 0 implies all contacts are different.
"""

import optparse
import math


def get_coordinates(m, fname):
    m.read(file=fname)

    def get_coord(res):
        if 'CA' in res.atoms:
            return res.atoms['CA']
        elif 'P' in res.atoms:
            return res.atoms['P']
    return [get_coord(r) for r in m.residues]


def get_distance(ci, cj):
    if ci is not None and cj is not None:
        dx = ci.x - cj.x
        dy = ci.y - cj.y
        dz = ci.z - cj.z
        return dx * dx + dy * dy + dz * dz


def get_distances(coord, rcut):
    rcut2 = rcut * rcut
    dist = {}
    for i in range(len(coord) - 1):
        for j in range(i + 1, len(coord)):
            d = get_distance(coord[i], coord[j])
            if d is not None and d < rcut2:
                dist[(i, j)] = dist[(j, i)] = math.sqrt(d)
    return dist


def get_qi_ca(m, len_coord, dist, template, avg_qi_cut, navg_qi_cut, fh):
    coord = get_coordinates(m, template)
    if len(coord) != len_coord:
        raise ValueError("different numbers of residues")
    qi_cut = [0.] * len(coord)
    nqi_cut = [0] * len(coord)
    for i in range(len(coord)):
        for j in range(len(coord)):
            if abs(i - j) < 2 or (i, j) not in dist:
                continue
            d = get_distance(coord[i], coord[j])
            if d is not None:
                delta = (dist[(i, j)] - math.sqrt(d)) / (abs(j - i) ** 0.15)
                qi_cut[i] += math.exp(-delta * delta * 0.5)
                nqi_cut[i] += 1
    for i in range(len(coord)):
        if nqi_cut[i] > 0:
            qi_cut[i] /= nqi_cut[i]
            avg_qi_cut[i] += qi_cut[i]
            navg_qi_cut[i] += 1
        else:
            qi_cut[i] = 1.1
        print("%6d %9.4f" % (i+1, qi_cut[i]), file=fh)


def get_qiavg_ca(target, templates, rcut):
    import modeller

    modeller.log.none()
    e = modeller.Environ()
    e.io.hetatm = False

    m = modeller.Model(e)
    coord = get_coordinates(m, target)
    dist = get_distances(coord, rcut)
    avg_qi_cut = [0.] * len(coord)
    navg_qi_cut = [0] * len(coord)
    for n, template in enumerate(templates):
        get_qi_ca(m, len(coord), dist, template, avg_qi_cut, navg_qi_cut,
                  open('qi_%d.dat' % (n+1), 'w'))
    with open('qi_avg.dat', 'w') as fh:
        for i in range(len(coord)):
            if navg_qi_cut[i] > 0:
                avg_qi_cut[i] /= navg_qi_cut[i]
            else:
                avg_qi_cut[i] = 1.1
            print("%6d %9.4f" % (i+1, avg_qi_cut[i]), file=fh)


def parse_args():
    usage = """%prog <target PDB> <cutoff> <template PDB> [...]

Calculate qiavg for CA atoms.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) < 3:
        parser.error("incorrect number of arguments")
    return args[0], float(args[1]), args[2:]


def main():
    target, cutoff, templates = parse_args()
    get_qiavg_ca(target, templates, cutoff)


if __name__ == '__main__':
    main()
