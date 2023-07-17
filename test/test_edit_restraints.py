import unittest
import subprocess
import os
import sys
if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from io import BytesIO as StringIO
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)

import allosmod.get_contacts  # noqa: E402
import allosmod.get_ss  # noqa: E402
import allosmod.edit_restraints  # noqa: E402


class TruncatedGaussianRestraint(allosmod.edit_restraints.Restraint):
    def handle_parameters(self, params):
        self.delE, self.slope, self.scl_delx = [float(x) for x in params[:3]]
        del params[:3]
        self.weights = [float(x) for x in params[:self.modal]]
        self.means = [float(x) for x in params[self.modal:self.modal*2]]
        self.stdevs = [float(x) for x in params[self.modal*2:self.modal*3]]


class MockRestraintEditor(allosmod.edit_restraints.RestraintEditor):
    def __init__(self):
        sigmas = allosmod.edit_restraints.Sigmas(2, 2.0, 3.0, 4.0)
        allosmod.edit_restraints.RestraintEditor.__init__(
            self, "dummyoth.rsr", "dummyas.rsr", "dummypdbfile", ["contact"],
            "dummyatomlist", sigmas, 10.0, 0.2, None, False, False)
        self.contacts = allosmod.edit_restraints.ContactMap()
        self.beta_structure = {}
        self.HETscale *= 4.0  # fail if HETscale not in parent

    def check_parse_restraint(self, r, delEmax=10.0):
        from allosmod.edit_restraints import TruncatedGaussianParameters
        r_from_form = {3: allosmod.edit_restraints.GaussianRestraint,
                       4: allosmod.edit_restraints.MultiGaussianRestraint,
                       7: allosmod.edit_restraints.CosineRestraint,
                       9: allosmod.edit_restraints.BinormalRestraint,
                       10: allosmod.edit_restraints.SplineRestraint,
                       50: TruncatedGaussianRestraint}

        tgparams = TruncatedGaussianParameters(delEmax, delEmaxNUC=20.0,
                                               slope=3.0, scl_delx=4.0,
                                               breaks={})
        fh = StringIO()
        self.parse_restraint(tgparams, r, fh)
        fh.seek(0)
        for line in fh:
            if line.startswith('R'):
                typ, form, rest = line.split(None, 2)
                yield r_from_form[int(form)](line, r.atoms)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to edit_restraints"""
        for args in ([], [''] * 6):
            check_output(['allosmod', 'edit_restraints'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.edit_restraints'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_setup_atoms(self):
        """Test setup_atoms()"""
        import modeller
        from allosmod.edit_restraints import RestraintEditor
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        with open('atomlistASRS', 'w') as fh:
            for i in range(64):
                fh.write("%d %s\n" % (i+1, "AS" if i < 10 else "RS"))
        with open('break.dat', 'w') as fh:
            fh.write('1 20')
        env = modeller.environ()
        env.io.hetatm = True

        e = RestraintEditor(
            'test.rsr', 'test.rsr',
            os.path.join(test_dir, 'input', 'test_editrsr.pdb'),
            ['test.pdb'], 'atomlistASRS', sigmas, 10.0, 0.1,
            'break.dat', False, False)

        class Residue(object):
            pass

        def mock_get_cont(fname, rcut):
            ri = Residue()
            ri.index = 1
            rj = Residue()
            rj.index = 4
            yield ri, rj, 10.0

        def mock_get_ss(pdb_file):
            return []
        with utils.mock_method(allosmod.get_contacts, 'get_contacts',
                               mock_get_cont):
            with utils.mock_method(allosmod.get_ss, 'get_ss', mock_get_ss):
                e.setup_atoms(env)
        contacts = sorted(e.contacts.keys())
        # Should have the 1-4 interaction from get_contacts, plus the two
        # nucleic acids (#6 and #7) should interact with everything
        self.assertEqual(contacts, [(1, 4), (1, 6), (1, 7), (2, 6), (2, 7),
                                    (3, 6), (3, 7), (4, 6), (4, 7), (5, 6),
                                    (5, 7), (6, 6), (6, 7), (7, 7)])
        self.assertEqual(e.breaks, {1: 20.0})
        # First 10 atoms should be allosteric site
        self.assertEqual([a.isAS for a in e.atoms], [True]*10 + [False]*54)
        # One CA and CB for each residue
        self.assertEqual(len([a for a in e.atoms if a.isCA]), 4)
        self.assertEqual(len([a for a in e.atoms if a.isCB]), 4)
        self.assertTrue(e.atoms[1].isCA)
        self.assertTrue(e.atoms[2].isCB)
        self.assertEqual(len([a for a in e.atoms if a.isSC]), 5)
        self.assertTrue(e.atoms[3].isSC)
        self.assertEqual(len([a for a in e.atoms if a.isNUC]), 39)
        self.assertEqual(len([a for a in e.atoms if a.torestr]), 18)
        os.unlink('break.dat')
        os.unlink('atomlistASRS')

    def test_sigmas(self):
        """Test Sigmas class"""
        class Atom(object):
            def __init__(self, isAS, isSC):
                self.isAS, self.isSC = isAS, isSC
        # Test one template, allosteric site, BB-BB
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=True, isSC=False),
                        Atom(isAS=True, isSC=False)))
        self.assertAlmostEqual(g, 1.0 * 1.0, places=1)
        # Test one template, allosteric site, SC-BB
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=True, isSC=True),
                        Atom(isAS=True, isSC=False)))
        self.assertAlmostEqual(g, 1.5 * 1.0, places=1)
        g = sigmas.get((Atom(isAS=True, isSC=False),
                        Atom(isAS=True, isSC=True)))
        self.assertAlmostEqual(g, 1.5 * 1.0, places=1)
        # Test one template, allosteric site, SC-SC
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=True, isSC=True),
                        Atom(isAS=True, isSC=True)))
        self.assertAlmostEqual(g, 1.5 * 1.5 * 1.0, places=1)
        # Test two templates, allosteric site, SC-SC
        sigmas = allosmod.edit_restraints.Sigmas(2, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=True, isSC=True),
                        Atom(isAS=True, isSC=True)))
        self.assertAlmostEqual(g, 1.0, places=1)
        g = sigmas.get_scaled((Atom(isAS=True, isSC=True),
                               Atom(isAS=True, isSC=True)))
        self.assertAlmostEqual(g, 4.0, places=1)
        # Test one template, regulated site, SC-SC
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=False, isSC=True),
                        Atom(isAS=False, isSC=True)))
        self.assertAlmostEqual(g, 1.5 * 1.5 * 2.0, places=1)
        # Test one template, interface, SC-SC
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        g = sigmas.get((Atom(isAS=True, isSC=True),
                        Atom(isAS=False, isSC=True)))
        self.assertAlmostEqual(g, 1.5 * 1.5 * 3.0, places=1)
        g = sigmas.get((Atom(isAS=False, isSC=True),
                        Atom(isAS=True, isSC=True)))
        self.assertAlmostEqual(g, 1.5 * 1.5 * 3.0, places=1)

    def test_truncated_gaussian_parameters(self):
        """Test TruncatedGaussianParameters class"""
        from allosmod.edit_restraints import TruncatedGaussianParameters

        class Residue(object):
            pass

        class Atom(object):
            def __init__(self, ri):
                self.a = self
                self.a.residue = Residue()
                self.a.residue.index = ri

        tgparams = TruncatedGaussianParameters(delEmax=1.0, delEmaxNUC=2.0,
                                               slope=3.0, scl_delx=4.0,
                                               breaks={})
        # nuc
        e = tgparams.get_dele(None, local=False, nuc=True)
        self.assertAlmostEqual(e, 2.0, places=1)
        # non-local, no breaks
        e = tgparams.get_dele((Atom(1), Atom(2)), local=False, nuc=False)
        self.assertAlmostEqual(e, 1.0, places=1)
        # local, no breaks
        e = tgparams.get_dele((Atom(1), Atom(2)), local=True, nuc=False)
        self.assertAlmostEqual(e, 10.0, places=1)

        # breaks
        tgparams = TruncatedGaussianParameters(delEmax=2.0, delEmaxNUC=6.0,
                                               slope=3.0, scl_delx=4.0,
                                               breaks={1: 20.0})
        for local in True, False:
            e = tgparams.get_dele((Atom(1), Atom(2)), local=local, nuc=False)
            self.assertAlmostEqual(e, 40.0, places=1)

    def test_restraint(self):
        """Test Restraint base class"""
        from allosmod.edit_restraints import Restraint

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        self.assertRaises(ValueError, Restraint,
                          "R 3 1 9 12 2 2 1 3 2 10.00 20.00",
                          [Atom(i) for i in range(1, 10)])

    def make_restraint(self, modify_atom_func, args, natom, cls, fmt):
        from allosmod.edit_restraints import Atom

        class ModellerResidue(object):
            hetatm = False

        class ModellerAtom(object):
            def __init__(self, ind):
                self.index = ind
                self.residue = ModellerResidue()
        atoms = [Atom(ModellerAtom(i+1)) for i in range(natom)]
        if modify_atom_func:
            modify_atom_func(atoms, args)
        r = cls(fmt % (len(atoms),
                       ' '.join('%d' % (x+1) for x in range(natom))), atoms)
        return r

    def make_gaussian_restraint(self, modify_atom_func, args=None, natom=2):
        from allosmod.edit_restraints import GaussianRestraint
        return self.make_restraint(modify_atom_func, args, natom,
                                   GaussianRestraint,
                                   "R 3 1 9 12 %d 2 1 %s 10.00 20.00")

    def make_multi_gaussian_restraint(self, modify_atom_func, args=None,
                                      natom=2):
        from allosmod.edit_restraints import MultiGaussianRestraint
        return self.make_restraint(modify_atom_func, args, natom,
                                   MultiGaussianRestraint,
                                   "R 4 2 9 12 %d 6 1 %s 0.8 0.2 10.00 "
                                   "20.00 5.0 8.0")

    def make_cosine_restraint(self, modify_atom_func, args=None, natom=2):
        from allosmod.edit_restraints import CosineRestraint
        return self.make_restraint(modify_atom_func, args, natom,
                                   CosineRestraint,
                                   "R 7 2 9 12 %d 2 1 %s 20.0 30.0")

    def make_spline_restraint(self, modify_atom_func=None, args=None, natom=2):
        from allosmod.edit_restraints import SplineRestraint
        return self.make_restraint(modify_atom_func, args, natom,
                                   SplineRestraint,
                                   "R 10 22 3 13 %d 3 1 %s x y z")

    def test_is_intrahet(self):
        """Test Restraint.is_intrahet()"""
        def set_hetatm(atoms, hetatm):
            for a, h in zip(atoms, hetatm):
                a.a.residue.hetatm = h
        r = self.make_gaussian_restraint(set_hetatm, [True, True])
        self.assertTrue(r.is_intrahet())
        r = self.make_gaussian_restraint(set_hetatm, [True, False])
        self.assertFalse(r.is_intrahet())
        r = self.make_gaussian_restraint(set_hetatm, [False, False])
        self.assertFalse(r.is_intrahet())

    def test_is_ca_cb_interaction(self):
        """Test Restraint.is_ca_cb_interaction()"""
        def modify_atoms(atoms, nuc_ca_cb):
            for a, n in zip(atoms, nuc_ca_cb):
                a.isNUC, a.isCA, a.isCB = n
        # No atom is CA, CB, or a nucleotide
        r = self.make_gaussian_restraint(modify_atoms, [[False, False, False],
                                                        [True, True, True]])
        self.assertFalse(r.is_ca_cb_interaction())
        # nuc-nuc interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, False, False],
                                                        [True, False, False]])
        self.assertTrue(r.is_ca_cb_interaction())
        # CA-nuc interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, False, False],
                                                        [False, True, False]])
        self.assertTrue(r.is_ca_cb_interaction())
        # CA-CB interaction
        r = self.make_gaussian_restraint(modify_atoms, [[False, False, True],
                                                        [False, True, False]])
        self.assertTrue(r.is_ca_cb_interaction())

    def test_is_sidechain_sidechain_interaction(self):
        """Test Restraint.is_sidechain_sidechain_interaction()"""
        def modify_atoms(atoms, sc_cb):
            for a, s in zip(atoms, sc_cb):
                a.isSC, a.isCB = s
        # SC-SC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, False],
                                                        [True, False]])
        self.assertTrue(r.is_sidechain_sidechain_interaction())
        # SC-CB interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, False],
                                                        [False, True]])
        self.assertTrue(r.is_sidechain_sidechain_interaction())
        # BB-BB interaction
        r = self.make_gaussian_restraint(modify_atoms, [[False, False],
                                                        [False, False]])
        self.assertFalse(r.is_sidechain_sidechain_interaction())

    def test_is_beta_beta_interaction(self):
        """Test Restraint.is_beta_beta_interaction()"""
        def modify_atoms(atoms, rind):
            for a, r in zip(atoms, rind):
                a.a.residue.index = r
        beta_structure = {4: True}
        # beta-beta interaction
        r = self.make_gaussian_restraint(modify_atoms, [4, 4])
        self.assertTrue(r.is_beta_beta_interaction(beta_structure))
        # beta-nonbeta interaction
        r = self.make_gaussian_restraint(modify_atoms, [4, 2])
        self.assertFalse(r.is_beta_beta_interaction(beta_structure))

    def test_is_intra_protein_interaction(self):
        """Test Restraint.is_intra_protein_interaction()"""
        def modify_atoms(atoms, nuc):
            for a, n in zip(atoms, nuc):
                a.isNUC = n
        # protein-protein interaction
        r = self.make_gaussian_restraint(modify_atoms, [False, False])
        self.assertTrue(r.is_intra_protein_interaction())
        # protein-nucleotide interaction
        r = self.make_gaussian_restraint(modify_atoms, [False, True])
        self.assertFalse(r.is_intra_protein_interaction())

    def test_is_intra_dna_interaction(self):
        """Test Restraint.is_intra_dna_interaction()"""
        def modify_atoms(atoms, nuc_r):
            for a, n in zip(atoms, nuc_r):
                a.isNUC, a.torestr = n
        # restrNUC-restrNUC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, True],
                                                        [True, True]])
        self.assertTrue(r.is_intra_dna_interaction())
        # restrNUC-NUC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, True],
                                                        [True, False]])
        self.assertFalse(r.is_intra_dna_interaction())
        # restrNUC-protein interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, True],
                                                        [False, False]])
        self.assertFalse(r.is_intra_dna_interaction())

    def test_is_protein_dna_interaction(self):
        """Test Restraint.is_protein_dna_interaction()"""
        def modify_atoms(atoms, nuc_r):
            for a, n in zip(atoms, nuc_r):
                a.isNUC, a.torestr = n
        # protein-restrNUC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[False, False],
                                                        [True, True]])
        self.assertTrue(r.is_protein_dna_interaction())
        # protein-unrestrNUC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[False, False],
                                                        [True, False]])
        self.assertFalse(r.is_protein_dna_interaction())
        # protein-protein interaction
        r = self.make_gaussian_restraint(modify_atoms, [[False, False],
                                                        [False, False]])
        self.assertFalse(r.is_protein_dna_interaction())
        # restrNUC-restrNUC interaction
        r = self.make_gaussian_restraint(modify_atoms, [[True, True],
                                                        [True, True]])
        self.assertFalse(r.is_protein_dna_interaction())

    def test_is_allosteric_interaction(self):
        """Test Restraint.is_allosteric_interaction()"""
        def modify_atoms(atoms, allos):
            for a, n in zip(atoms, allos):
                a.isAS = n
        # AS-AS interaction
        r = self.make_gaussian_restraint(modify_atoms, [True, True])
        self.assertTrue(r.is_allosteric_interaction())
        # AS-RS interaction
        r = self.make_gaussian_restraint(modify_atoms, [True, False])
        self.assertFalse(r.is_allosteric_interaction())
        # RS-RS interaction
        r = self.make_gaussian_restraint(modify_atoms, [False, False])
        self.assertFalse(r.is_allosteric_interaction())

    def test_gaussian_restraint(self):
        """Test GaussianRestraint class"""
        from allosmod.edit_restraints import GaussianRestraint
        from allosmod.edit_restraints import TruncatedGaussianParameters

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = GaussianRestraint("R 3 1 9 12 2 2 1 3 2 10.00 20.00",
                              [Atom(i) for i in range(1, 10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        self.assertEqual(r.modal, 1)
        self.assertEqual(r.feat, 9)
        self.assertEqual(r.group, 12)
        self.assertAlmostEqual(r.mean, 10.0, places=1)
        self.assertAlmostEqual(r.firstmean, 10.0, places=1)
        self.assertAlmostEqual(r.stdev, 20.0, places=1)
        self.assertTrue(r.any_mean_below(15.0))
        self.assertFalse(r.any_mean_below(5.0))
        self.assertFalse(r.any_mean_below(10.0))
        s = StringIO()
        r.write(s)
        self.assertEqual(s.getvalue(), 'R    3   1   9  12   2   2   1     '
                                       '3     2      10.0000   20.0000\n')
        # transform to multigaussian
        s = StringIO()
        r.transform(None, 2, 30.0, truncated=False, fh=s)
        self.assertEqual(s.getvalue(),
                         'R    4   2   9  12   2   6   1     3     2       '
                         '0.5000    0.5000   10.0000   10.0000   '
                         '30.0000   30.0000\n')
        # transform to truncated gaussian
        s = StringIO()
        tgparams = TruncatedGaussianParameters(delEmax=2.0, delEmaxNUC=3.0,
                                               slope=4.0, scl_delx=5.0,
                                               breaks={})
        r.transform(tgparams, 2, 30.0, truncated=True, nuc=True, fh=s)
        self.assertEqual(
            s.getvalue(),
            'R   50   2   9  12   2   9   1     3     2       '
            '3.0000    4.0000    5.0000    0.5000    0.5000   '
            '10.0000   10.0000   30.0000   30.0000\n')
        r.rescale(4.0)
        self.assertAlmostEqual(r.stdev, 5.0, places=1)

    def test_multi_gaussian_restraint(self):
        """Test MultiGaussianRestraint class"""
        from allosmod.edit_restraints import MultiGaussianRestraint
        from allosmod.edit_restraints import TruncatedGaussianParameters

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = MultiGaussianRestraint("R 4 2 9 12 2 6 1 3 2 0.8 0.2 "
                                   "10.00 20.00 30.00 40.00",
                                   [Atom(i) for i in range(1, 10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        self.assertEqual(r.modal, 2)
        self.assertEqual(r.feat, 9)
        self.assertEqual(r.group, 12)
        self.assertEqual(len(r.means), 2)
        self.assertEqual(len(r.stdevs), 2)
        self.assertAlmostEqual(r.means[0], 10.0, places=1)
        self.assertAlmostEqual(r.means[1], 20.0, places=1)
        self.assertAlmostEqual(r.firstmean, 10.0, places=1)
        self.assertAlmostEqual(r.stdevs[0], 30.0, places=1)
        self.assertAlmostEqual(r.stdevs[1], 40.0, places=1)
        self.assertTrue(r.any_mean_below(15.0))
        self.assertFalse(r.any_mean_below(5.0))
        self.assertFalse(r.any_mean_below(10.0))
        s = StringIO()
        r.write(s)
        self.assertEqual(
            s.getvalue(),
            'R    4   2   9  12   2   6   1     3     2       '
            '0.8000    0.2000   10.0000   20.0000   30.0000   40.0000\n')
        # transform to multigaussian
        s = StringIO()
        r.transform(None, 8, 70.0, truncated=False, fh=s)
        self.assertEqual(s.getvalue(),
                         'R    4   2   9  12   2   6   1     3     2       '
                         '0.5000    0.5000   10.0000   20.0000   '
                         '70.0000   70.0000\n')
        # transform to truncated gaussian
        s = StringIO()
        tgparams = TruncatedGaussianParameters(delEmax=2.0, delEmaxNUC=3.0,
                                               slope=4.0, scl_delx=5.0,
                                               breaks={})
        r.transform(tgparams, 8, 70.0, truncated=True, nuc=True, fh=s)
        self.assertEqual(
            s.getvalue(),
            'R   50   2   9  12   2   9   1     3     2       '
            '3.0000    4.0000    5.0000    0.5000    0.5000   '
            '10.0000   20.0000   70.0000   70.0000\n')
        r.rescale(5.0)
        self.assertAlmostEqual(r.stdevs[0], 6.0, places=1)
        self.assertAlmostEqual(r.stdevs[1], 8.0, places=1)

    def test_cosine_restraint(self):
        """Test CosineRestraint class"""
        from allosmod.edit_restraints import CosineRestraint

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = CosineRestraint("R 7 2 9 12 2 2 1 3 2 2.0 4.0",
                            [Atom(i) for i in range(1, 10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        self.assertEqual(r.modal, 2)
        self.assertEqual(r.feat, 9)
        self.assertEqual(r.group, 12)
        self.assertAlmostEqual(r.phase, 2.0, places=1)
        self.assertAlmostEqual(r.force, 4.0, places=1)
        s = StringIO()
        r.write(s)
        self.assertEqual(
            s.getvalue(),
            'R    7   2   9  12   2   2   1     3     '
            '2       2.0000    4.0000\n')
        r.rescale(2.5)
        self.assertAlmostEqual(r.force, 10.0, places=1)

    def test_binormal_restraint(self):
        """Test BinormalRestraint class"""
        from allosmod.edit_restraints import BinormalRestraint

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = BinormalRestraint("R 9 2 9 12 2 2 1 3 2 x y z",
                              [Atom(i) for i in range(1, 10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        s = StringIO()
        r.write(s)
        self.assertEqual(s.getvalue(),
                         'R    9   2   9  12   2   2   1     '
                         '3     2            x         y         z\n')

    def test_spline_restraint(self):
        """Test SplineRestraint class"""
        from allosmod.edit_restraints import SplineRestraint

        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = SplineRestraint("R 10 2 9 12 2 2 1 3 2 x y z",
                            [Atom(i) for i in range(1, 10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        s = StringIO()
        r.write(s)
        self.assertEqual(s.getvalue(),
                         'R   10   2   9  12   2   2   1     '
                         '3     2            x         y         z\n')

    def test_restraint_filters(self):
        """Test restraint filters"""
        from allosmod.edit_restraints import filter_rs_rs, filter_not_rs_rs

        class Atom(object):
            def __init__(self, isAS):
                self.isAS = isAS
        as_as = (Atom(True), Atom(True))
        as_rs = (Atom(True), Atom(False))
        rs_as = (Atom(False), Atom(True))
        rs_rs = (Atom(False), Atom(False))

        self.assertTrue(filter_rs_rs(rs_rs))
        self.assertFalse(filter_rs_rs(rs_as))
        self.assertFalse(filter_rs_rs(as_rs))
        self.assertFalse(filter_rs_rs(as_as))

        self.assertFalse(filter_not_rs_rs(rs_rs))
        self.assertTrue(filter_not_rs_rs(rs_as))
        self.assertTrue(filter_not_rs_rs(as_rs))
        self.assertTrue(filter_not_rs_rs(as_as))

    def test_add_ca_boundary_restraints(self):
        """Test add_ca_boundary_restraints()"""
        from allosmod.edit_restraints import add_ca_boundary_restraints

        class Atom(object):
            def __init__(self, ind, isCA):
                self.isCA = isCA
                self.index = ind
                self.a = self
        atoms = [Atom(5, False), Atom(7, True)]
        s = StringIO()
        add_ca_boundary_restraints(atoms, s)
        restraints = s.getvalue().rstrip('\n').split('\n')
        self.assertEqual(len(restraints), 6)
        # Only atom #7 should be restrained
        for r in restraints:
            self.assertEqual(r[34:37], " 7 ")

    def test_parse_restraints_file(self):
        """Test parse_restraints_file()"""
        from allosmod.edit_restraints import parse_restraints_file
        from allosmod.edit_restraints import filter_rs_rs

        class Atom(object):
            def __init__(self, ind, isAS):
                self.index = ind
                self.a = self
                self.isAS = isAS
        atoms = [Atom(i, True) for i in range(1, 10)]
        # no filter
        s = StringIO("R 3 1 9 12 2 2 1 3 2 10.00 20.00\n\n")
        rs = list(parse_restraints_file(s, atoms))
        self.assertEqual(len(rs), 1)
        self.assertEqual([a.a.index for a in rs[0].atoms], [3, 2])

        # rs-rs filter
        s = StringIO("R 3 1 9 12 2 2 1 3 2 10.00 20.00\n\n")
        rs = list(parse_restraints_file(s, atoms, filter_rs_rs))
        self.assertEqual(len(rs), 0)

    def test_atom(self):
        """Test Atom class"""
        from allosmod.edit_restraints import Atom
        a = Atom('foo')
        self.assertEqual(a.isAS, False)
        self.assertEqual(a.isNUC, False)
        self.assertEqual(a.isSC, False)
        self.assertEqual(a.isCA, False)
        self.assertEqual(a.isCB, False)
        self.assertEqual(a.torestr, False)
        self.assertEqual(a.a, 'foo')

    def test_contact_map(self):
        """Test ContactMap class"""
        from allosmod.edit_restraints import ContactMap, Atom

        class ModellerResidue(object):
            pass

        class ModellerAtom(object):
            def __init__(self, ind):
                self.residue = ModellerResidue()
                self.residue.index = ind
        c = ContactMap()
        self.assertEqual(list(c.keys()), [])
        self.assertFalse(c[(1, 4)])
        c[(1, 4)] = True
        c[(5, 2)] = True
        self.assertTrue(c[(1, 4)])
        self.assertTrue(c[(4, 1)])
        self.assertTrue(c[(2, 5)])
        self.assertTrue(c[(5, 2)])
        self.assertTrue(c[(Atom(ModellerAtom(5)), 2)])
        self.assertTrue(c[(5, Atom(ModellerAtom(2)))])
        self.assertEqual(len(c.keys()), 2)

    def test_get_beta(self):
        """Test get_beta()"""
        from allosmod.edit_restraints import get_beta

        def mock_get_ss(pdb_file):
            if pdb_file == 'empty':
                return []
            elif pdb_file == 'all_beta':
                return ['E']*3
            elif pdb_file == 'all_helix':
                return ['H']*10
            elif pdb_file == 'some_beta':
                return ['', 'E', '']
        with utils.mock_method(allosmod.get_ss, 'get_ss', mock_get_ss):
            self.assertEqual(get_beta('empty'), {})
            self.assertEqual(get_beta('all_beta'), {1: True, 2: True, 3: True})
            self.assertEqual(get_beta('all_helix'), {})
            self.assertEqual(get_beta('some_beta'), {2: True})

    def test_get_nuc_restrained(self):
        """Test get_nuc_restrained()"""
        from allosmod.edit_restraints import get_nuc_restrained
        self.assertFalse(get_nuc_restrained('OP1', 'any residue'))
        self.assertFalse(get_nuc_restrained('OP2', 'any residue'))
        self.assertTrue(get_nuc_restrained("O3'", 'any residue'))
        self.assertTrue(get_nuc_restrained('N1', 'ADE'))
        self.assertTrue(get_nuc_restrained('C2', 'DT'))
        self.assertTrue(get_nuc_restrained('O2', 'U'))
        self.assertTrue(get_nuc_restrained('O6', 'G'))
        self.assertTrue(get_nuc_restrained('N1', 'CYT'))
        self.assertFalse(get_nuc_restrained('N2', 'CYT'))
        self.assertFalse(get_nuc_restrained('N2', 'URA'))

    def test_parse_other(self):
        """Test parse of restraints unknown to AllosMod"""
        class MockRestraint(object):
            atoms = []

            def is_intrahet(self):
                return False
        e = MockRestraintEditor()
        r2 = list(e.check_parse_restraint(MockRestraint()))
        # Unknown restraints are ignored
        self.assertEqual(len(r2), 0)

    def test_parse_coarse_ca_ca_intra_protein(self):
        """Test parse of coarse AS-AS CA-CA intra-protein restraint"""
        e = MockRestraintEditor()
        e.coarse = True
        e.contacts[(1, 2)] = True  # non-local interaction

        def modify_atoms(atoms, arg):
            atoms[0].isAS = atoms[1].isAS = True  # AS-AS
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
        self.assertAlmostEqual(r2[0].delE, 10.0, places=1)
        self.assertAlmostEqual(r2[0].slope, 3.0, places=1)
        self.assertAlmostEqual(r2[0].scl_delx, 4.0, places=1)
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 8.0, places=1)
        # with tgauss_AS off
        e.tgauss_AS = False
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.GaussianRestraint)
        self.assertAlmostEqual(r2[0].stdev, 8.0, places=1)
        # with empty_AS on
        e.tgauss_AS = True
        e.empty_AS = True
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 0)

    def test_parse_ca_ca_multi_intra_protein(self):
        """Test parse of CA-CA multigauss intra-protein restraint"""
        e = MockRestraintEditor()
        e.contacts[(1, 2)] = True  # non-local interaction

        def modify_atoms(atoms, arg):
            atoms[0].isAS = atoms[1].isAS = True  # AS-AS
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_multi_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
        self.assertAlmostEqual(r2[0].delE, 10.0, places=1)
        self.assertAlmostEqual(r2[0].slope, 3.0, places=1)
        self.assertAlmostEqual(r2[0].scl_delx, 4.0, places=1)
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 2.0, places=1)
        # with delEmax = 0
        r2 = list(e.check_parse_restraint(r, delEmax=0.))
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.MultiGaussianRestraint)
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 2.0, places=1)

    def test_parse_rs_ca_ca_intra_protein(self):
        """Test parse of RS-RS CA-CA intra-protein restraint"""
        e = MockRestraintEditor()
        e.coarse = True
        e.contacts[(1, 2)] = True  # non-local interaction

        def modify_atoms(atoms, arg):
            atoms[0].isAS = atoms[1].isAS = False  # RS-RS
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
        self.assertAlmostEqual(r2[0].delE, 10.0, places=1)
        self.assertAlmostEqual(r2[0].slope, 3.0, places=1)
        self.assertAlmostEqual(r2[0].scl_delx, 4.0, places=1)
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 12.0, places=1)
        # with delEmax = 0
        r2 = list(e.check_parse_restraint(r, delEmax=0.))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.MultiGaussianRestraint)
        self.assertEqual(len(r2[0].stdevs), 1)
        self.assertAlmostEqual(r2[0].stdevs[0], 12.0, places=1)

    def test_parse_bond_restraint(self):
        """Test parse of bond restraint"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            for a, h in zip(atoms, arg):
                a.a.residue.hetatm = h
        for het, scale in (False, 1.0), (True, 4.0):
            r = self.make_gaussian_restraint(modify_atoms, [het, het])
            r.group = 1
            r2 = list(e.check_parse_restraint(r))
            self.assertEqual(len(r2), 1)
            self.assertEqual(type(r2[0]),
                             allosmod.edit_restraints.GaussianRestraint)
            self.assertAlmostEqual(r2[0].mean, 10.0, places=1)
            self.assertAlmostEqual(r2[0].stdev, 20.0 / scale, places=1)

    def test_parse_angle_restraint(self):
        """Test parse of angle/dihedral restraint"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            for a, h in zip(atoms, arg):
                a.a.residue.hetatm = h
        for natom in (3, 4):
            for het, scale in (False, 1.0), (True, 4.0):
                r = self.make_gaussian_restraint(modify_atoms, [het]*natom,
                                                 natom=natom)
                r.group = 1
                r2 = list(e.check_parse_restraint(r))
                self.assertEqual(len(r2), 1)
                self.assertEqual(type(r2[0]),
                                 allosmod.edit_restraints.GaussianRestraint)
                self.assertAlmostEqual(r2[0].mean, 10.0, places=1)
                self.assertAlmostEqual(r2[0].stdev, 20.0 / scale, places=1)

    def test_parse_multi_angle_restraint(self):
        """Test parse of multigauss angle/dihedral restraint"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            for a, h in zip(atoms, arg):
                a.a.residue.hetatm = h
        for natom in (3, 4):
            for het, scale in (False, 1.0), (True, 4.0):
                r = self.make_multi_gaussian_restraint(
                    modify_atoms, [het]*natom, natom=natom)
                r2 = list(e.check_parse_restraint(r))
                self.assertEqual(len(r2), 1)
                self.assertEqual(
                    type(r2[0]),
                    allosmod.edit_restraints.MultiGaussianRestraint)
                self.assertEqual(len(r2[0].means), 2)
                self.assertAlmostEqual(r2[0].weights[0], 0.8, places=1)
                self.assertAlmostEqual(r2[0].weights[1], 0.2, places=1)
                self.assertAlmostEqual(r2[0].means[0], 10.0, places=1)
                self.assertAlmostEqual(r2[0].means[1], 20.0, places=1)
                self.assertAlmostEqual(r2[0].stdevs[0], 5.0 / scale, places=1)
                self.assertAlmostEqual(r2[0].stdevs[1], 8.0 / scale, places=1)

    def test_parse_cosine_restraint(self):
        """Test parse of cosine restraint"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            for a, h in zip(atoms, arg):
                a.a.residue.hetatm = h
        for het, scale in (False, 1.0), (True, 4.0):
            r = self.make_cosine_restraint(modify_atoms, [het, het])
            r2 = list(e.check_parse_restraint(r))
            self.assertEqual(len(r2), 1)
            self.assertEqual(type(r2[0]),
                             allosmod.edit_restraints.CosineRestraint)
            self.assertAlmostEqual(r2[0].phase, 20.0, places=1)
            self.assertAlmostEqual(r2[0].force, 30.0 * scale, places=1)

    def test_parse_spline_restraint(self):
        """Test parse of spline restraint"""
        # should pass through as-is
        e = MockRestraintEditor()
        r = self.make_spline_restraint()
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.SplineRestraint)
        self.assertEqual(r2[0]._params, ['x', 'y', 'z'])

    def test_parse_coarse_not_ca_cb(self):
        """Test parse of coarse restraint, not CA-CB"""
        e = MockRestraintEditor()
        e.coarse = True

        def modify_atoms(atoms, arg):
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        # Restraint should be omitted
        self.assertEqual(len(r2), 0)

    def test_parse_sidechain_too_long(self):
        """Test parse of sidechain-sidechain restraint, mean too big"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            atoms[0].isSC = atoms[1].isSC = True
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r.mean = 5.1
        r2 = list(e.check_parse_restraint(r))
        # Restraint should be omitted
        self.assertEqual(len(r2), 0)

    def test_local_2to5_ca_cb(self):
        """Test locrigid CA-CB restraint with res range 2-5"""
        e = MockRestraintEditor()
        e.locrigid = True

        def modify_atoms(atoms, seqdst):
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 30
            atoms[1].a.residue.index = 30 + seqdst
        for seqdst in 2, 5, -2, -5:
            r = self.make_gaussian_restraint(modify_atoms, seqdst)
            r2 = list(e.check_parse_restraint(r))
            # sigma should be set to 2.0
            self.assertEqual(len(r2), 1)
            self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
            self.assertEqual(len(r2[0].stdevs), 2)
            self.assertAlmostEqual(r2[0].stdevs[0], 2.0, places=1)
            self.assertAlmostEqual(r2[0].stdevs[1], 2.0, places=1)

    def test_local_6to12_ca_cb(self):
        """Test locrigid CA-CB restraint with res range 6-12"""
        e = MockRestraintEditor()
        e.locrigid = True

        def modify_atoms(atoms, seqdst):
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 30
            atoms[1].a.residue.index = 30 + seqdst
        for seqdst in 6, 12, -6, -12:
            r = self.make_gaussian_restraint(modify_atoms, seqdst)
            r.mean = 5.9
            r2 = list(e.check_parse_restraint(r))
            # sigma should be set to 2.0
            self.assertEqual(len(r2), 1)
            self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
            self.assertEqual(len(r2[0].stdevs), 2)
            self.assertAlmostEqual(r2[0].stdevs[0], 2.0, places=1)
            self.assertAlmostEqual(r2[0].stdevs[1], 2.0, places=1)
            # If mean >= 6.0, restraint is omitted
            r = self.make_gaussian_restraint(modify_atoms, seqdst)
            r.mean = 6.1
            r2 = list(e.check_parse_restraint(r))
            self.assertEqual(len(r2), 0)

    def test_under3_ca_cb(self):
        """Test CA-CB restraint with res range <= 2"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            atoms[0].isCA = atoms[1].isCA = True  # CA-CA
            atoms[0].a.residue.index = 30
            atoms[1].a.residue.index = 32
        r = self.make_gaussian_restraint(modify_atoms)
        e.beta_structure = {30: True, 32: True}
        r.mean = 5.9
        r2 = list(e.check_parse_restraint(r))
        # sigma should be set to 2.0
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 2.0, places=1)
        self.assertAlmostEqual(r2[0].stdevs[1], 2.0, places=1)
        # If mean >= 6.0, restraint is omitted
        r = self.make_gaussian_restraint(modify_atoms)
        r.mean = 6.1
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 0)
        # If not beta structure, restraint is omitted
        r = self.make_gaussian_restraint(modify_atoms)
        r.mean = 5.9
        e.beta_structure = {}
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 0)

    def test_parse_het_het(self):
        """Test parse of het-het restraint"""
        e = MockRestraintEditor()

        def modify_atoms(atoms, arg):
            atoms[0].a.residue.hetatm = True
            atoms[1].a.residue.hetatm = True
        r = self.make_gaussian_restraint(modify_atoms)
        r2 = list(e.check_parse_restraint(r))
        # keep as is (but scaled by HETscale, 4.0)
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.GaussianRestraint)
        self.assertAlmostEqual(r2[0].mean, 10.0, places=1)
        self.assertAlmostEqual(r2[0].stdev, 20.0 / 4.0, places=1)

    def test_protein_dna_restraint(self):
        """Test parse of protein-dna restraint"""
        e = MockRestraintEditor()
        e.contacts[(1, 2)] = True  # non-local interaction

        def modify_atoms(atoms, arg):
            atoms[1].isNUC = atoms[1].torestr = True  # protein-DNA, RS-RS
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r.mean = 7.9
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]), TruncatedGaussianRestraint)
        self.assertAlmostEqual(r2[0].delE, 20.0, places=1)  # delEmaxNUC
        self.assertEqual(len(r2[0].stdevs), 2)
        self.assertAlmostEqual(r2[0].stdevs[0], 3.0, places=1)
        self.assertAlmostEqual(r2[0].stdevs[1], 3.0, places=1)
        # If mean > rcutNUC (8.0), restraint should be omitted
        r.mean = 8.1
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 0)

    def test_intra_dna_restraint(self):
        """Test parse of intra-dna restraint"""
        e = MockRestraintEditor()
        e.contacts[(1, 2)] = True  # non-local interaction

        def modify_atoms(atoms, arg):
            atoms[0].isNUC = atoms[0].torestr = True  # DNA-DNA, RS-RS
            atoms[1].isNUC = atoms[1].torestr = True
            atoms[0].a.residue.index = 1
            atoms[1].a.residue.index = 2
        r = self.make_gaussian_restraint(modify_atoms)
        r.mean = 7.9
        # stdev should be forced to 1.0, mean not changed
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 1)
        self.assertEqual(type(r2[0]),
                         allosmod.edit_restraints.GaussianRestraint)
        self.assertAlmostEqual(r2[0].mean, 7.9, places=1)
        self.assertAlmostEqual(r2[0].stdev, 1.0, places=1)
        # If mean > rcutNUC (8.0), restraint should be omitted
        r.mean = 8.1
        r2 = list(e.check_parse_restraint(r))
        self.assertEqual(len(r2), 0)

    def test_setup_delEmax_no_coarse(self):
        """Test setup_delEmax(), coarse=False"""
        e = MockRestraintEditor()
        # coarse=False; delE* should be unchanged
        e.setup_delEmax()
        self.assertAlmostEqual(e.delEmax, 0.2, places=1)
        self.assertAlmostEqual(e.delEmaxNUC, 0.12, places=2)

    def test_setup_delEmax_no_rsr(self):
        """Test setup_delEmax(), no restraints"""
        e = MockRestraintEditor()
        # empty file; delE* should be unchanged
        e.coarse = True
        e.atoms = []
        open('dummyas.rsr', 'w').close()

        def mock_parse(fh, atoms):
            return []
        with utils.mock_method(allosmod.edit_restraints,
                               'parse_restraints_file', mock_parse):
            e.setup_delEmax()
        self.assertAlmostEqual(e.delEmax, 0.2, places=1)
        self.assertAlmostEqual(e.delEmaxNUC, 0.12, places=2)
        os.unlink('dummyas.rsr')

    def test_setup_delEmax_rsr(self):
        """Test setup_delEmax(), with some restraints"""
        e = MockRestraintEditor()
        e.contacts[(1, 2)] = True
        e.coarse = True
        e.atoms = []
        open('dummyas.rsr', 'w').close()

        def mock_parse(fh, atoms):
            def set_resind(atoms):
                atoms[0].a.residue.index = 1
                atoms[1].a.residue.index = 2

            def make_ca_ca(atoms, arg):
                set_resind(atoms)
                atoms[0].isCA = atoms[1].isCA = True  # CA-CA

            def make_cb_cb(atoms, arg):
                set_resind(atoms)
                atoms[0].isCB = atoms[1].isCB = True  # CB-CB

            def make_sc_sc(atoms, arg):
                set_resind(atoms)
                atoms[0].isSC = atoms[1].isSC = True  # SC-SC
            # CA-CA restraint
            r = self.make_gaussian_restraint(make_ca_ca)
            yield r
            # CB-CB restraint
            r = self.make_gaussian_restraint(make_cb_cb)
            yield r
            # SC-SC restraint, but short (<distco_scsc)
            r = self.make_gaussian_restraint(make_sc_sc)
            r.mean = 4.9
            yield r
            # Restraint with natoms!=2
            r = self.make_gaussian_restraint(make_ca_ca, natom=3)
            yield r
            # Restraint not in contacts
            r = self.make_gaussian_restraint(make_ca_ca)
            r.atoms[0].a.residue.index = 3
            yield r
            # Restraint of wrong type
            r = self.make_cosine_restraint(make_ca_ca)
            yield r
        with utils.mock_method(allosmod.edit_restraints,
                               'parse_restraints_file', mock_parse):
            e.setup_delEmax()
        self.assertAlmostEqual(e.delEmax, 0.08, places=2)
        self.assertAlmostEqual(e.delEmaxNUC, 0.05, places=2)
        os.unlink('dummyas.rsr')

    def test_simple(self):
        """Simple complete run of edit_restraints"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'pm_test.pdb'), 'w') as fh:
                fh.write("""
ATOM      1  N   ARG A   1     -18.387  -9.167  -1.701  1.00  0.54           N
ATOM      2  CA  ARG A   1     -17.434  -9.856  -0.787  1.00  0.45           C
ATOM      3  C   ARG A   1     -15.998  -9.610  -1.251  1.00  0.45           C
ATOM      4  O   ARG A   1     -15.130 -10.444  -1.087  1.00  0.52           O
ATOM      5  CB  ARG A   1     -17.793 -11.337  -0.899  1.00  0.46           C
""")
            with open(os.path.join(tmpdir, 'test.rsr'), 'w') as fh:
                fh.write("R 3 1 9 12 2 2 1 1 2 10.00 20.00\n")
            with open(os.path.join(tmpdir, 'list4contacts'), 'w') as fh:
                fh.write("test.pdb\n")
            with open(os.path.join(tmpdir, 'atomlistASRS'), 'w') as fh:
                fh.write("1 AS\n2 AS\n")
            check_output(['allosmod', 'edit_restraints', 'test.rsr',
                          'test.rsr', 'pm_test.pdb', 'list4contacts',
                          'atomlistASRS'], cwd=tmpdir)


if __name__ == '__main__':
    unittest.main()
