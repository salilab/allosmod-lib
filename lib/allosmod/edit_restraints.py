"""Create AllosMod-specific restraints"""

from __future__ import print_function, absolute_import, division
import optparse
import modeller
import math
import sys
import allosmod.util
import allosmod.get_contacts
from allosmod.get_ss import get_ss

class Sigmas(object):
    def __init__(self, ntotal, sig_AS, sig_RS, sig_inter):
        self.ntotal = ntotal
        self.sig_AS = sig_AS
        self.sig_RS = sig_RS
        self.sig_inter = sig_inter

    def get(self, atoms):
        """Get sigma most appropriate for the given atoms"""
        isSC1, isSC2 = atoms[0].isSC, atoms[1].isSC
        isAS1, isAS2 = atoms[0].isAS, atoms[1].isAS

        if not isSC1 and isSC2: # BB-BB
            sig_scale = 1.0
        elif isSC1 != isSC2: # SC-BB
            sig_scale=1.5
        else: # SC-SC
            sig_scale=1.5*1.5

        if isAS1 and isAS2:
            if self.ntotal == 1: # if one template, avoid treating like AS
                return self.sig_AS*sig_scale
            else: # if allosteric site, do not scale side chain interactions
                return self.sig_AS
        elif not isAS1 and not isAS2:
            return self.sig_RS*sig_scale
        else: # interface
            return self.sig_inter*sig_scale

class TruncatedGaussianParameters(object):
    def __init__(self, delEmax, slope, scl_delx, breaks):
        self.delEmax, self.slope, self.scl_delx = delEmax, slope, scl_delx
        self.breaks = breaks

    def get_dele(self, atoms):
        bscale = 1.0
        is_break = False
        for a in atoms:
            ri = a.a.residue.index
            if ri in self.breaks:
                is_break = True
                bscale *= self.breaks[ri]
        if is_break:
            return bscale * self.delEmax, bscale * self.delEmax
        else:
            return self.delEmax, self.delEmax * 10.0

class Restraint(object):
    def __init__(self, line, atoms):
        self.line = line
        spl = line.split()
        self.form = int(spl[1])
        self.modal = int(spl[2])
        self.feat = int(spl[3])
        self.group = int(spl[4])
        natom = int(spl[5])
        self.nparam = int(spl[6])
        self.nfeat = int(spl[7])
        ind = 7 + min(self.nfeat, 1)
        self.atoms = [atoms[int(i)-1] for i in spl[ind:ind+natom]]
        self.handle_parameters(spl[ind+natom:])

    def handle_parameters(self, params):
        raise ValueError("Could not handle %d" % self.form)

    def write(self, fh=sys.stdout):
        fh.write('R %4d %4d %4d %4d %4d %4d %4d '
                 % (self.form, self.modal, self.feat, self.group,
                    len(self.atoms), self.nparam, 0))
        fh.write(' '.join('%6d' % x.a.index for x in self.atoms))
        fh.write(' ')
        self.write_parameters(fh)
        fh.write('\n')

    def is_intrahet(self):
        """Return True iff all atoms in this restraint are HETATM"""
        for a in self.atoms:
            if not a.a.residue.hetatm:
                return False
        return True

    def is_ca_cb_interaction(self):
        for a in self.atoms:
            if not a.isNUC and not a.isCA and not a.isCB:
                return False
        return True

    def is_sidechain_sidechain_interaction(self):
        for a in self.atoms:
            if not a.isSC and not a.isCB:
                return False
        return True

    def is_beta_beta_interaction(self, beta_structure):
        for a in self.atoms:
            if a.a.residue.index not in beta_structure:
                return False
        return True

    def is_intra_protein_interaction(self):
        for a in self.atoms:
            if a.isNUC:
                return False
        return True

    def is_intra_dna_interaction(self):
        for a in self.atoms:
            if not a.isNUC or not a.torestr:
                return False
        return True

    def is_protein_dna_interaction(self):
        dna = [a for a in self.atoms if a.isNUC]
        if len(dna) == 0 or len(dna) == len(self.atoms):
            return False
        for a in dna:
            if a.torestr:
                return True
        return False

    def is_allosteric_interaction(self):
        for a in self.atoms:
            if not a.isAS:
                return False
        return True


class GaussianRestraint(Restraint):
    def handle_parameters(self, params):
        self.mean = self.firstmean = float(params[0])
        self.stdev = float(params[1])
    def write_parameters(self, fh):
        fh.write(' '.join('%9.4f' % x for x in (self.mean, self.stdev)))

    def any_mean_below(self, threshold):
        return self.mean < threshold

    def transform(self, tgparams, modal, stdev, truncated=True,
                  fh=sys.stdout):
        """Convert this restraint into a multigaussian, and write out"""
        if truncated:
            delE, delEmaxLOC = tgparams.get_dele(self.atoms)
            parameters = [delEmaxLOC, tgparams.slope, tgparams.scl_delx]
        else:
            parameters = []
        parameters.extend([1.0/modal]*modal + [self.mean] * modal
                          + [stdev]*modal)
        fh.write('R %4d %4d %4d %4d %4d %4d %4d '
                 % (50 if truncated else 4, modal, self.feat, self.group,
                    len(self.atoms), len(parameters), 0))
        fh.write(' '.join('%6d' % x.a.index for x in self.atoms))
        fh.write(' ')
        fh.write(' '.join('%9.4f' % x for x in parameters))
        fh.write('\n')

class MultiGaussianRestraint(Restraint):
    def handle_parameters(self, params):
        self.weights = float(params[:self.modal])
        self.means = float(params[self.modal:self.modal*2])
        self.firstmean = self.means[0]
        self.stdevs = float(params[self.model*2:self.modal*3])
    def write_parameters(self, fh):
        fh.write(' '.join('%9.4f' % x for x in
                                     self.weights + self.means + self.stdevs))

    def any_mean_below(self, threshold):
        for m in self.means:
            if m < threshold:
                return True
        return False

    def transform(self, tgparams, modal, stdev, truncated=True,
                  fh=sys.stdout):
        """Convert this restraint into a multigaussian, and write out"""
        # Note that modal is ignored
        modal = self.modal
        if truncated:
            delE, delEmaxLOC = tgparams.get_dele(self.atoms)
            parameters = [delEmaxLOC, tgparams.slope, tgparams.scl_delx]
        else:
            parameters = []
        parameters.extend([1.0/modal]*modal + self.means + [stdev]*modal)
        fh.write('R %4d %4d %4d %4d %4d %4d %4d '
                 % (50 if truncated else 4, modal, self.feat, self.group,
                    len(self.atoms), len(parameters), 0))
        fh.write(' '.join('%6d' % x.a.index for x in self.atoms))
        fh.write(' ')
        fh.write(' '.join('%9.4f' % x for x in parameters))
        fh.write('\n')

class CosineRestraint(Restraint):
    def handle_parameters(self, params):
        self.phase, self.force = [float(x) for x in params]
    def write_parameters(self, fh):
        fh.write(' '.join('%9.4f' % x for x in (self.phase, self.force)))

class BinormalRestraint(Restraint):
    def handle_parameters(self, params):
        # Keep as is
        self._params = params
    def write_parameters(self, fh):
        fh.write(' '.join('%9s' % x for x in self._params))

class SplineRestraint(Restraint):
    def handle_parameters(self, params):
        # Keep as is
        self._params = params
    def write_parameters(self, fh):
        fh.write(' '.join('%9s' % x for x in self._params))

def filter_rs_rs(atoms):
    """Allow only RS-RS contacts"""
    for atom in atoms:
        if atom.isAS:
            return False
    return True

def filter_not_rs_rs(atoms):
    """Allow AS-AS and AS-RS contacts"""
    natomRS = sum(0 if atom.isAS else 1 for atom in atoms)
    return natomRS != len(atoms)

def add_ca_boundary_restraints(atoms, fh=sys.stdout):
    """Add restraints to CA's enforce boundary conditions: cube soft boundary"""
    for a in atoms:
        if a.isCA:
            for feat in range(9,12):
                for (form, mean) in ((2, 100.0), (1, -100.0)):
                    fh.write('R %4d %4d %4d %4d %4d %4d %4d %6d    '
                             '%9.4f %9.4f\n'
                             % (form, 0, feat, 27, 1, 2, 1, a.a.index, mean,
                                10.0))

def parse_restraints_file(fh, atoms, filter=None):
    restraint_from_form = {3: GaussianRestraint,
                           4: MultiGaussianRestraint,
                           7: CosineRestraint,
                           9: BinormalRestraint,
                           10: SplineRestraint}
    for line in fh:
        if line.startswith('R'):
            typ, form, rest = line.split(None, 2)
            r = restraint_from_form.get(int(form), Restraint)(line, atoms)
            if filter is None or filter(r.atoms):
                yield r

class Atom(object):
    isAS = isNUC = isSC = isCA = isCB = torestr = False
    def __init__(self, a):
        self.a = a

class ContactMap(object):
    def __init__(self):
        self.__d = {}

    def __getitem__(self, key):
        i,j = key
        if isinstance(i, Atom):
            i = i.a.residue.index
        if isinstance(j, Atom):
            j = j.a.residue.index
        i, j = min(i,j), max(i,j)
        return (i,j) in self.__d
    def __setitem__(self, key, val):
        i,j = key
        i, j = min(i,j), max(i,j)
        return self.__d.__setitem__((i, j), True)

def parse_atomlist_asrs(atomlist_asrs):
    retval = {'AS': True, 'RS': False}
    for line in atomlist_asrs:
        num, allos_type = line.rstrip('\n\r').split()
        yield retval[allos_type]

def get_contacts(contacts_pdbs, rcut):
    contacts = ContactMap()
    for f in contacts_pdbs:
        for ri, rj, dist in allosmod.get_contacts.get_contacts('pm_' + f, rcut):
            contacts[(ri.index, rj.index)] = True
    return contacts

def get_breaks(fh):
    breaks = {}
    for line in fh:
        resind, scale = line.rstrip('\r\n').split()
        breaks[int(resind)] = float(scale)
    return breaks

def get_beta(pdb_file):
    beta = {}
    dssp = list(get_ss(pdb_file))
    beta_fraction = dssp.count('E') / len(dssp)
    helix_fraction = dssp.count('H') / len(dssp)
    beta_structure = beta_fraction > 0.20 and helix_fraction < 0.05
    if beta_structure:
        for n,x in dssp:
            if x == 'E':
                # Residue indices start at 1
                beta[n+1] = True
    return beta

def get_nuc_restrained(atm_name, res_name):
    if res_name in ('ADE', 'A', 'DA'):
        return atm_name in ('N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7',
                            'C8', 'N9')
    elif res_name in ('THY', 'T', 'DT'):
        return atm_name in ('N1', 'C2', 'O2', 'N1', 'N3', 'C4', 'O4', 'C5',
                            'C6', 'C7')
    elif res_name in ('URA', 'U', 'DU'):
        return atm_name in ('N1', 'C2', 'O2', 'N1', 'N3', 'C4', 'O4', 'C5',
                            'C6')
    elif res_name in ('GUA', 'G', 'DG'):
        return atm_name in ('N1', 'N2', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6',
                            'O6', 'N7', 'C8', 'N9')
    else:
        return atm_name in ('N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6')


def test(listoth_rsr, listas_rsr, pdb_file, contacts_pdbs, atomlist_asrs, sig_AS, sig_RS, sig_inter, rcut, ntotal, delEmax, break_file, coarse, locrigid):
    # if empty_AS, then only use: cov bonds, angles, dihedrals for AS
    # (ignore nonbonded contacts within AS site, RS and interface OK)
    empty_AS=False
    # if tgauss_AS, then AS nonbonded contacts are converted to truncated
    # Gaussian, instead of being constrained
    tgauss_AS=True
    #scale intraHET contacts to prevent exploding 
    HETscale=1.0
    #set nucleotide energy params
    delEmaxNUC=0.12
    rcutNUC=8.0
    distco_scsc=5.0
    sigmas = Sigmas(ntotal, sig_AS, sig_RS, sig_inter)

    modeller.log.none()
    e = modeller.environ()
    e.io.hetatm = True
    m = modeller.model(e, file=pdb_file)
    atoms = [Atom(a) for a in m.atoms]
    contacts = get_contacts(contacts_pdbs, rcut)
    breaks = get_breaks(open(break_file)) if break_file else {}
    tgparams = TruncatedGaussianParameters(delEmax=delEmax, slope=4.0,
                                           scl_delx=0.7, breaks=breaks)
    beta_structure = get_beta(pdb_file)
    NUCLEIC_ACIDS = dict.fromkeys(['ADE', 'A', 'DA', 'THY', 'T', 'DT', 'URA', 'U', 'DU', 'GUA', 'G', 'DG', 'CYT', 'C', 'DC'])
    BACKBONE_ATOMS = dict.fromkeys(['CA', 'CB', 'O', 'N', 'C', 'OT', 'NA', 'NB', 'NC', 'ND', 'C1A', 'C2A', 'C3A', 'C4A', 'C1B', 'C2B', 'C3B', 'C4B', 'C1C', 'C2C', 'C3C', 'C4C', 'C1D', 'C2D', 'C3D', 'C4D'])
    for a in atoms:
        r = a.a.residue
        if r.pdb_name in NUCLEIC_ACIDS:
            a.isNUC = True
            a.torestr = get_nuc_restrained(a.a.name, r.pdb_name)
            for rj in range(1, len(m.residues) + 1):
                contacts[(r.index,rj)] = True
        if a.a.name in BACKBONE_ATOMS or r.pdb_name in NUCLEIC_ACIDS:
            a.isSC = False
            a.isCA = a.a.name == 'CA'
            a.isCB = a.a.name == 'CB'
        else:
            a.isSC = a.a.name != 'H'
    for a, asrs in zip(atoms, parse_atomlist_asrs(open(atomlist_asrs))):
        a.isAS = asrs

    if coarse:
        ndist = ndistCACB = 0
        for r in parse_restraints_file(open(listas_rsr), atoms):
            if (isinstance(r, GaussianRestraint)
                and len(r.atoms) == 2 and r.group != 1) \
               or (isinstance(r, MultiGaussianRestraint)
                   and len(r.atoms) == 2):
                if contacts[(r.atoms[0], r.atoms[1])]:
                    if (not r.atoms[0].isSC and not r.atoms[0].isCB) \
                       or (not r.atoms[1].isSC and not r.atoms[1].isCB) \
                       or r.firstmean < distco_scsc:
                        ndist += 1
                    if (r.atoms[0].isCA or r.atoms[0].isCB) \
                       and (r.atoms[1].isCA or r.atoms[1].isCB):
                        ndistCACB += 1
        delEmax = (6.5 / 7.8) * (ndist / ndistCACB) * delEmax
        delEmaxNUC = (6.5 / 7.8) * (ndist / ndistCACB) * delEmaxNUC

    print("MODELLER5 VERSION: MODELLER FORMAT")
    # get restraints
    for (fname, filter) in ((listoth_rsr, filter_rs_rs),
                            (listas_rsr, filter_not_rs_rs)):
        for r in parse_restraints_file(open(fname), atoms, filter):
            # gaussian; bond, angle or torsion
            # keep as is for prot, scale for HET
            if isinstance(r, GaussianRestraint) and \
               ((len(r.atoms) == 2 and r.group == 1) \
                or len(r.atoms) == 3 or len(r.atoms) == 4):
                if r.is_intrahet():
                    r.stdev *= HETscale
                r.write()
            # multigaussian; angle or torsion
            # keep as is for prot, scale for HET
            elif isinstance(r, MultiGaussianRestraint) and \
                len(r.atoms) >= 3:
                if r.is_intrahet():
                    r.stdevs = [x * HETscale for x in r.stdevs]
                r.write()
            # cosine dihedrals for backbone dihedrals and 
            # side chain dihedrals with few alignments
            elif isinstance(r, CosineRestraint):
                if r.is_intrahet():
                    r.force *= HETscale
                r.write()
            # Keep as is
            elif isinstance(r, SplineRestraint):
                r.write()
            # Gaussian or MultiGaussian distance restraint
            elif len(r.atoms) == 2 \
                 and ((isinstance(r, GaussianRestraint) and r.group != 1) \
                      or (isinstance(r, MultiGaussianRestraint))):
                if coarse and not r.is_ca_cb_interaction():
                    continue
                # omit side chain interactions > 5 Ang
                if r.is_sidechain_sidechain_interaction() \
                   and not r.any_mean_below(distco_scsc):
                    continue
                seqdst = abs(r.atoms[0].a.residue.index
                             - r.atoms[1].a.residue.index)
                if locrigid and 2 <= seqdst <= 5 and r.is_ca_cb_interaction():
                    r.transform(tgparams, modal=2, stdev=2.0)
                elif locrigid and 6 <= seqdst <= 12 \
                   and r.is_ca_cb_interaction() and r.any_mean_below(6.0):
                    r.transform(tgparams, modal=2, stdev=2.0)
                elif 2 <= seqdst and r.any_mean_below(6.0) \
                   and r.is_ca_cb_interaction() \
                   and r.is_beta_beta_interaction(beta_structure):
                    r.transform(tgparams, modal=2, stdev=2.0)
                elif contacts[(r.atoms[0], r.atoms[1])]:
                    if isinstance(r, MultiGaussianRestraint):
                        sig = sigmas.get(r.atoms) * ntotal * ntotal
                        r.transform(tgparams, modal=r.modal,
                                    stdev=sig, truncated=delEmax != 0.)
                    elif r.is_intra_protein_interaction():
                        sig = sigmas.get(r.atoms) * ntotal * ntotal
                        if r.is_allosteric_interaction():
                            if not empty_AS:
                                if tgauss_AS:
                                    r.transform(tgparams, modal=2,
                                                stdev=sig)
                                else:
                                    r.stdev = sig
                                    r.write()
                        else: # RS or interface
                            if delEmax == 0.:
                                r.transform(tgparams, modal=1,
                                            stdev=sig, truncated=False)
                            else:
                                r.transform(tgparams, modal=2,
                                            stdev=sig)
                    elif r.is_protein_dna_interaction() \
                         and r.any_mean_below(rcutNUC):
                        sig = sigmas.get(r.atoms)
                        r.transform(tgparams, modal=2, stdev=sig)
                    elif r.is_intra_dna_interaction() \
                         and r.any_mean_below(rcutNUC):
                        r.stdev = 1.0
                        r.write()
            # add intra heme contacts to maintain geometry (not included
            # above due to contmap), keep all as is
            elif r.is_intrahet() and len(r.atoms) == 2:
                if isinstance(r, GaussianRestraint) and r.group != 1:
                    r.stdev *= HETscale
                    r.write()
                elif isinstance(r, MultiGaussianRestraint):
                    r.stdevs = [x * HETscale for x in r.stdevs]
                    r.write()
    add_ca_boundary_restraints(atoms)

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
    return test('listOTH.rsr', 'listAS.rsr', 'pm_3UWP.pdb', ['3UWP.pdb'], 'atomlistASRS', 2.0, 2.0, 2.0, 11.0, 1, 0.01, None, False, False)

if __name__ == '__main__':
    main()
