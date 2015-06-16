import unittest
import subprocess
import os
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
        from allosmod.edit_restraints import RestraintEditor, Sigmas
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

if __name__ == '__main__':
    unittest.main()
