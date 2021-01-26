"""Calculate q for CA atoms."""

from __future__ import print_function, absolute_import
import optparse
import math
from allosmod.get_qiavg_ca import get_distance, get_coordinates


class QScore(object):
    def __init__(self):
        self.q = 0.
        self.count = 0

    def add(self, val):
        self.q += val
        self.count += 1

    def average(self):
        if self.count == 0:
            return 0.
        else:
            return self.q / self.count


class QScores(object):
    def __init__(self):
        self.q = QScore()
        self.qs = QScore()
        self.qm = QScore()
        self.ql = QScore()

    def average(self):
        return [x.average() for x in (self.q, self.qs, self.qm, self.ql)]

    def add(self, i, j, dist):
        dij = abs(j-i)
        delta = dist / (dij ** 0.15)
        qcont = math.exp(-delta * delta * 0.5)
        self.q.add(qcont)
        if dij < 5:
            self.qs.add(qcont)
        elif dij < 13:
            self.qm.add(qcont)
        else:
            self.ql.add(qcont)


def get_distances(coord, resdelta=1):
    dist = {}
    for i in range(len(coord) - resdelta):
        for j in range(i + resdelta, len(coord)):
            d = get_distance(coord[i], coord[j])
            if d is not None:
                dist[(i, j)] = math.sqrt(d)
    return dist


def get_qi_ca(m, len_coord, dist, template, rcut):
    q_tot = QScores()
    q_cut = QScores()
    coord = get_coordinates(m, template)
    if len(coord) != len_coord:
        raise ValueError("different numbers of residues")
    for i in range(len(coord) - 2):
        for j in range(i + 2, len(coord)):
            if (i, j) not in dist:
                continue
            d = get_distance(coord[i], coord[j])
            if d is not None:
                delta = dist[(i, j)] - math.sqrt(d)
                q_tot.add(i, j, delta)
                if dist[(i, j)] < rcut:
                    q_cut.add(i, j, delta)
    return q_tot, q_cut


def write_q_scores(qs, fh):
    for n, q in enumerate(qs):
        print("%6d %9.4f  %9.4f  %9.4f  %9.4f" % ((n+1,) + tuple(q.average())),
              file=fh)


def get_q_ca(target, templates, rcut):
    import modeller

    modeller.log.none()
    e = modeller.Environ()
    e.io.hetatm = False

    m = modeller.Model(e)
    coord = get_coordinates(m, target)
    dist = get_distances(coord)
    q_tot = []
    q_cut = []
    for template in templates:
        t, c = get_qi_ca(m, len(coord), dist, template, rcut)
        q_tot.append(t)
        q_cut.append(c)
    write_q_scores(q_tot, open('qscore1to%d.dat' % len(coord), 'w'))
    write_q_scores(q_cut, open('qs_cut1to%d.dat' % len(coord), 'w'))


def parse_args():
    usage = """%prog <target PDB> <cutoff> <template PDB> [...]

Calculate q for CA atoms.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) < 3:
        parser.error("incorrect number of arguments")
    return args[0], float(args[1]), args[2:]


def main():
    target, cutoff, templates = parse_args()
    get_q_ca(target, templates, cutoff)


if __name__ == '__main__':
    main()
