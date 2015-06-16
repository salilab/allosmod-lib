import unittest
import subprocess
import os
from io import BytesIO
import sys
from test_pdb2ali import check_output
from test_modeller import mock_method
import allosmod.get_contacts
import allosmod.get_ss

test_dir = os.path.dirname(sys.argv[0])

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to edit_restraints"""
        for args in ([], ['' * 6]):
            out = check_output(['allosmod', 'edit_restraints'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.edit_restraints'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_setup_atoms(self):
        """Test setup_atoms()"""
        import modeller
        from allosmod.edit_restraints import RestraintEditor, Sigmas
        sigmas = allosmod.edit_restraints.Sigmas(1, 1.0, 2.0, 3.0)
        with open('atomlistASRS', 'w') as fh:
            for i in range(64):
               fh.write("%d %s\n" % (i+1, "AS" if i < 10 else "RS"))
        with open('break.dat', 'w') as fh:
            fh.write('1 20')
        env = modeller.environ()
        env.io.hetatm = True

        e = RestraintEditor('test.rsr', 'test.rsr',
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
        with mock_method(allosmod.get_contacts, 'get_contacts', mock_get_cont):
            with mock_method(allosmod.get_ss, 'get_ss', mock_get_ss):
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
        self.assertEqual(len([a for a in e.atoms if a.torestr]), 22)
        os.unlink('break.dat')
        os.unlink('atomlistASRS')

    def test_sigmas(self):
        """Test Sigmas class"""
        from allosmod.edit_restraints import Sigmas
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
                          [Atom(i) for i in range(1,10)])

    def make_restraint(self, modify_atom_func, *args):
        from allosmod.edit_restraints import GaussianRestraint, Atom
        class ModellerResidue(object):
            hetatm = False
        class ModellerAtom(object):
            def __init__(self, ind):
                self.index = ind
                self.residue = ModellerResidue()
        atoms = [Atom(ModellerAtom(1)), Atom(ModellerAtom(2))]
        modify_atom_func(atoms, *args)
        r = GaussianRestraint("R 3 1 9 12 2 2 1 1 2 10.00 20.00", atoms)
        return r
        
    def test_is_intrahet(self):
        """Test Restraint.is_intrahet()"""
        def set_hetatm(atoms, hetatm):
            for a, h in zip(atoms, hetatm):
                a.a.residue.hetatm = h
        r = self.make_restraint(set_hetatm, [True, True])
        self.assertTrue(r.is_intrahet())
        r = self.make_restraint(set_hetatm, [True, False])
        self.assertFalse(r.is_intrahet())
        r = self.make_restraint(set_hetatm, [False, False])
        self.assertFalse(r.is_intrahet())

    def test_is_ca_cb_interaction(self):
        """Test Restraint.is_ca_cb_interaction()"""
        def modify_atoms(atoms, nuc_ca_cb):
            for a, n in zip(atoms, nuc_ca_cb):
                a.isNUC, a.isCA, a.isCB = n
        # No atom is CA, CB, or a nucleotide
        r = self.make_restraint(modify_atoms, [[False, False, False],
                                               [True, True, True]])
        self.assertFalse(r.is_ca_cb_interaction())
        # nuc-nuc interaction
        r = self.make_restraint(modify_atoms, [[True, False, False],
                                               [True, False, False]])
        self.assertTrue(r.is_ca_cb_interaction())
        # CA-nuc interaction
        r = self.make_restraint(modify_atoms, [[True, False, False],
                                               [False, True, False]])
        self.assertTrue(r.is_ca_cb_interaction())
        # CA-CB interaction
        r = self.make_restraint(modify_atoms, [[False, False, True],
                                               [False, True, False]])
        self.assertTrue(r.is_ca_cb_interaction())

    def test_is_sidechain_sidechain_interaction(self):
        """Test Restraint.is_sidechain_sidechain_interaction()"""
        def modify_atoms(atoms, sc_cb):
            for a, s in zip(atoms, sc_cb):
                a.isSC, a.isCB = s
        # SC-SC interaction
        r = self.make_restraint(modify_atoms, [[True, False], [True, False]])
        self.assertTrue(r.is_sidechain_sidechain_interaction())
        # SC-CB interaction
        r = self.make_restraint(modify_atoms, [[True, False], [False, True]])
        self.assertTrue(r.is_sidechain_sidechain_interaction())
        # BB-BB interaction
        r = self.make_restraint(modify_atoms, [[False, False], [False, False]])
        self.assertFalse(r.is_sidechain_sidechain_interaction())

    def test_is_beta_beta_interaction(self):
        """Test Restraint.is_beta_beta_interaction()"""
        def modify_atoms(atoms, rind):
            for a, r in zip(atoms, rind):
                a.a.residue.index = r
        beta_structure = {4: True}
        # beta-beta interaction
        r = self.make_restraint(modify_atoms, [4, 4])
        self.assertTrue(r.is_beta_beta_interaction(beta_structure))
        # beta-nonbeta interaction
        r = self.make_restraint(modify_atoms, [4, 2])
        self.assertFalse(r.is_beta_beta_interaction(beta_structure))

    def test_is_intra_protein_interaction(self):
        """Test Restraint.is_intra_protein_interaction()"""
        def modify_atoms(atoms, nuc):
            for a, n in zip(atoms, nuc):
                a.isNUC = n
        # protein-protein interaction
        r = self.make_restraint(modify_atoms, [False, False])
        self.assertTrue(r.is_intra_protein_interaction())
        # protein-nucleotide interaction
        r = self.make_restraint(modify_atoms, [False, True])
        self.assertFalse(r.is_intra_protein_interaction())

    def test_is_intra_dna_interaction(self):
        """Test Restraint.is_intra_dna_interaction()"""
        def modify_atoms(atoms, nuc_r):
            for a, n in zip(atoms, nuc_r):
                a.isNUC, a.torestr = n
        # restrNUC-restrNUC interaction
        r = self.make_restraint(modify_atoms, [[True, True], [True, True]])
        self.assertTrue(r.is_intra_dna_interaction())
        # restrNUC-NUC interaction
        r = self.make_restraint(modify_atoms, [[True, True], [True, False]])
        self.assertFalse(r.is_intra_dna_interaction())
        # restrNUC-protein interaction
        r = self.make_restraint(modify_atoms, [[True, True], [False, False]])
        self.assertFalse(r.is_intra_dna_interaction())

    def test_is_protein_dna_interaction(self):
        """Test Restraint.is_protein_dna_interaction()"""
        def modify_atoms(atoms, nuc_r):
            for a, n in zip(atoms, nuc_r):
                a.isNUC, a.torestr = n
        # protein-restrNUC interaction
        r = self.make_restraint(modify_atoms, [[False, False], [True, True]])
        self.assertTrue(r.is_protein_dna_interaction())
        # protein-unrestrNUC interaction
        r = self.make_restraint(modify_atoms, [[False, False], [True, False]])
        self.assertFalse(r.is_protein_dna_interaction())
        # protein-protein interaction
        r = self.make_restraint(modify_atoms, [[False, False], [False, False]])
        self.assertFalse(r.is_protein_dna_interaction())
        # restrNUC-restrNUC interaction
        r = self.make_restraint(modify_atoms, [[True, True], [True, True]])
        self.assertFalse(r.is_protein_dna_interaction())

    def test_is_allosteric_interaction(self):
        """Test Restraint.is_allosteric_interaction()"""
        def modify_atoms(atoms, allos):
            for a, n in zip(atoms, allos):
                a.isAS = n
        # AS-AS interaction
        r = self.make_restraint(modify_atoms, [True, True])
        self.assertTrue(r.is_allosteric_interaction())
        # AS-RS interaction
        r = self.make_restraint(modify_atoms, [True, False])
        self.assertFalse(r.is_allosteric_interaction())
        # RS-RS interaction
        r = self.make_restraint(modify_atoms, [False, False])
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
                              [Atom(i) for i in range(1,10)])
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
        s = BytesIO()
        r.write(s)
        self.assertEqual(s.getvalue(), 'R    3   1   9  12   2   2   1     '
                                       '3     2      10.0000   20.0000\n')
        # transform to multigaussian
        s = BytesIO()
        r.transform(None, 2, 30.0, truncated=False, fh=s)
        self.assertEqual(s.getvalue(),
                         'R    4   2   9  12   2   6   1     3     2       '
                         '0.5000    0.5000   10.0000   10.0000   '
                         '30.0000   30.0000\n')
        # transform to truncated gaussian
        s = BytesIO()
        tgparams = TruncatedGaussianParameters(delEmax=2.0, delEmaxNUC=3.0,
                                           slope=4.0, scl_delx=5.0, breaks={})
        r.transform(tgparams, 2, 30.0, truncated=True, nuc=True, fh=s)
        self.assertEqual(s.getvalue(),
               'R   50   2   9  12   2   9   1     3     2       '
               '3.0000    4.0000    5.0000    0.5000    0.5000   '
               '10.0000   10.0000   30.0000   30.0000\n')

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
                                   [Atom(i) for i in range(1,10)])
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
        s = BytesIO()
        r.write(s)
        self.assertEqual(s.getvalue(), 
                   'R    4   2   9  12   2   6   1     3     2       '
                   '0.8000    0.2000   10.0000   20.0000   30.0000   40.0000\n')
        # transform to multigaussian
        s = BytesIO()
        r.transform(None, 8, 70.0, truncated=False, fh=s)
        self.assertEqual(s.getvalue(),
                         'R    4   2   9  12   2   6   1     3     2       '
                         '0.5000    0.5000   10.0000   20.0000   '
                         '70.0000   70.0000\n')
        # transform to truncated gaussian
        s = BytesIO()
        tgparams = TruncatedGaussianParameters(delEmax=2.0, delEmaxNUC=3.0,
                                           slope=4.0, scl_delx=5.0, breaks={})
        r.transform(tgparams, 8, 70.0, truncated=True, nuc=True, fh=s)
        self.assertEqual(s.getvalue(),
               'R   50   2   9  12   2   9   1     3     2       '
               '3.0000    4.0000    5.0000    0.5000    0.5000   '
               '10.0000   20.0000   70.0000   70.0000\n')

    def test_cosine_restraint(self):
        """Test CosineRestraint class"""
        from allosmod.edit_restraints import CosineRestraint
        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = CosineRestraint("R 7 2 9 12 2 2 1 3 2 2.0 4.0",
                            [Atom(i) for i in range(1,10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        self.assertEqual(r.modal, 2)
        self.assertEqual(r.feat, 9)
        self.assertEqual(r.group, 12)
        self.assertAlmostEqual(r.phase, 2.0, places=1)
        self.assertAlmostEqual(r.force, 4.0, places=1)
        s = BytesIO()
        r.write(s)
        self.assertEqual(s.getvalue(),
                      'R    7   2   9  12   2   2   1     3     '
                      '2       2.0000    4.0000\n')

    def test_binormal_restraint(self):
        """Test BinormalRestraint class"""
        from allosmod.edit_restraints import BinormalRestraint
        class Atom(object):
            def __init__(self, ind):
                self.index = ind
                self.a = self
        r = BinormalRestraint("R 9 2 9 12 2 2 1 3 2 x y z",
                              [Atom(i) for i in range(1,10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        s = BytesIO()
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
                            [Atom(i) for i in range(1,10)])
        self.assertEqual([a.a.index for a in r.atoms], [3, 2])
        s = BytesIO()
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
        s = BytesIO()
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
        atoms = [Atom(i, True) for i in range(1,10)]
        # no filter
        s = BytesIO("R 3 1 9 12 2 2 1 3 2 10.00 20.00\n\n")
        rs = list(parse_restraints_file(s, atoms))
        self.assertEqual(len(rs), 1)
        self.assertEqual([a.a.index for a in rs[0].atoms], [3, 2])

        # rs-rs filter
        s = BytesIO("R 3 1 9 12 2 2 1 3 2 10.00 20.00\n\n")
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
        self.assertEqual(c.keys(), [])
        self.assertFalse(c[(1,4)])
        c[(1,4)] = True
        c[(5,2)] = True
        self.assertTrue(c[(1,4)])
        self.assertTrue(c[(4,1)])
        self.assertTrue(c[(2,5)])
        self.assertTrue(c[(5,2)])
        self.assertTrue(c[(Atom(ModellerAtom(5)),2)])
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
        with mock_method(allosmod.get_ss, 'get_ss', mock_get_ss):
            self.assertEqual(get_beta('empty'), {})
            self.assertEqual(get_beta('all_beta'), {1:True, 2:True, 3:True})
            self.assertEqual(get_beta('all_helix'), {})
            self.assertEqual(get_beta('some_beta'), {2: True})

    def test_get_nuc_restrained(self):
        """Test get_nuc_restrained()"""
        from allosmod.edit_restraints import get_nuc_restrained
        self.assertTrue(get_nuc_restrained('OP1', 'any residue'))
        self.assertTrue(get_nuc_restrained('N1', 'ADE'))
        self.assertTrue(get_nuc_restrained('C2', 'DT'))
        self.assertTrue(get_nuc_restrained('O2', 'U'))
        self.assertTrue(get_nuc_restrained('O6', 'G'))
        self.assertTrue(get_nuc_restrained('N1', 'CYT'))
        self.assertFalse(get_nuc_restrained('N2', 'CYT'))
        self.assertFalse(get_nuc_restrained('N2', 'URA'))

if __name__ == '__main__':
    unittest.main()
