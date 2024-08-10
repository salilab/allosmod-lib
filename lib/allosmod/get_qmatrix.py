"""Calculate q matrix for given proteins."""

import optparse
import math
import random
from allosmod.get_qiavg_ca import get_coordinates, get_distance
from allosmod.get_q_ca import QScore


def get_squared_distances(coord):
    dist = {}
    for i in range(len(coord) - 2):
        for j in range(i + 2, len(coord)):
            d = get_distance(coord[i], coord[j])
            if d is not None:
                dist[(i, j)] = d
    return dist


def get_template_distances(m, numres, template):
    coord = get_coordinates(m, template)
    return get_squared_distances(coord)


def write_q_output(q_cuts, templates, mat_fh, avg_fh):
    def get_q(q_cuts, temp1, temp2):
        return q_cuts.get((temp1, temp2), None) \
               or q_cuts.get((temp2, temp1), None) or 1.0
    qavg = QScore()
    for itg, template in enumerate(templates):
        qs = [get_q(q_cuts, itg, jtg) for jtg in range(len(templates))]
        for jtg, q in enumerate(qs):
            if jtg != itg:
                qavg.add(q)
        print(template + "".join(" %.4f" % q for q in qs), file=mat_fh)
    print("Qa,b: %.2f" % qavg.average(), file=avg_fh)


def make_matrix_from_dists(templates, dists, numres, rcut2):
    q_cuts = {}
    for itg in range(len(templates) - 1):
        for jtg in range(itg + 1, len(templates)):
            q_cuts[(itg, jtg)] = q_cut = QScore()
            for i in range(numres - 2):
                for j in range(i + 2, numres):
                    ditg = dists[itg].get((i, j), None)
                    djtg = dists[jtg].get((i, j), None)
                    if ditg is None or djtg is None:
                        continue
                    if ditg < rcut2 or djtg < rcut2:
                        delta = (math.sqrt(ditg) - math.sqrt(djtg)) \
                                / (abs(j-i) ** 0.15)
                        q_cut.add(math.exp(-delta * delta * 0.5))
    for k in q_cuts.keys():
        q_cuts[k] = q_cuts[k].average()
    return q_cuts


def get_qmatrix(target, templates, rcut):
    import modeller

    rcut2 = rcut * rcut
    modeller.log.none()
    e = modeller.Environ()
    e.io.hetatm = False

    m = modeller.Model(e, file=target)
    numres = len(m.residues)
    dists = []
    for template in templates:
        dists.append(get_template_distances(m, numres, template))
    q_cuts = make_matrix_from_dists(templates, dists, numres, rcut2)
    write_q_output(q_cuts, templates, open('qmatrix.dat', 'w'),
                   open('cq_aq_qavg_qsd.dat', 'w'))


def shuffle_templates(templates):
    random.shuffle(templates)
    del templates[500:]


def parse_args():
    usage = """%prog <target PDB> <cutoff> <template PDB> [...]

Calculate q matrix for templates.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) < 3:
        parser.error("incorrect number of arguments")
    return args[0], float(args[1]), args[2:]


def main():
    target, cutoff, templates = parse_args()
    shuffle_templates(templates)
    get_qmatrix(target, templates, cutoff)


if __name__ == '__main__':
    main()
