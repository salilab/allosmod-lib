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

if __name__ == '__main__':
    unittest.main()
