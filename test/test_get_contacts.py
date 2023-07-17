import unittest
import os
import sys
import subprocess
import collections
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)


MockResidue = collections.namedtuple('MockResidue', ['pdb_name', 'atoms'])
MockAtom = collections.namedtuple('MockAtom', ['name', 'x', 'y', 'z'])


class AtomList(object):
    def __init__(self, *atoms):
        self.atoms = atoms
        self.atom_map = {}
        for a in atoms:
            self.atom_map[a.name] = a

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.atoms[key]
        else:
            return self.atom_map[key]


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_contacts"""
        for args in ([], ['x'] * 3):
            check_output(['allosmod', 'get_contacts'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.get_contacts'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_get_contact_type(self):
        """Test get_contact_type()"""
        from allosmod.get_contacts import get_contact_type

        class Residue(object):
            def __init__(self, r):
                self.pdb_name = r
        self.assertEqual(get_contact_type(Residue('ALA'), Residue('GLY')), 2)
        self.assertEqual(get_contact_type(Residue('CYS'), Residue('TYR')), 1)
        self.assertEqual(get_contact_type(Residue('ASN'), Residue('TRP')), 3)
        self.assertEqual(get_contact_type(Residue('TRP'), Residue('ASN')), 3)
        self.assertEqual(get_contact_type(Residue('HEM'), Residue('ASN')), 0)
        self.assertEqual(get_contact_type(Residue('TRP'), Residue('HEM')), 0)

    def test_get_average_aa(self):
        """Test _get_average_aa()"""
        from allosmod.get_contacts import _get_average_aa
        # Should take CA in preference to O, in preference to N
        CA = MockAtom(name='CA', x=10, y=0, z=0)
        Ox = MockAtom(name='O', x=20, y=0, z=0)
        N = MockAtom(name='N', x=30, y=0, z=0)
        a = _get_average_aa(MockResidue(pdb_name='HIS',
                                        atoms=AtomList(CA, Ox, N)))
        self.assertEqual(a.average, ((10, 0, 0),))
        a = _get_average_aa(MockResidue(pdb_name='HIS', atoms=AtomList(Ox, N)))
        self.assertEqual(a.average, ((20, 0, 0),))
        a = _get_average_aa(MockResidue(pdb_name='HIS', atoms=AtomList(N)))
        self.assertEqual(a.average, ((30, 0, 0),))
        self.assertRaises(ValueError, _get_average_aa,
                          MockResidue(pdb_name='HIS', atoms=AtomList()))

    def test_get_contact_dist(self):
        """Test get_contact_dist()"""
        from allosmod.get_contacts import get_contact_dist

        class Residue(object):
            def __init__(self, av):
                self.average = av
        self.assertEqual(get_contact_dist(Residue(((0, 0, 0),)),
                                          Residue(((0, 5, 0),)), 16.0),
                         None)
        self.assertAlmostEqual(get_contact_dist(Residue(((0, 0, 0),)),
                                                Residue(((0, 5, 0),)), 100.),
                               5.0, places=1)

    def test_simple(self):
        """Simple complete run of get_contacts"""
        out = check_output(['allosmod', 'get_contacts',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'), '11.0'],
                           universal_newlines=True)
        lines = out.split('\n')
        self.assertEqual(len(lines), 12)
        self.assertEqual(lines[0][43:48], "9.955")
        self.assertEqual(lines[1][43:48], "9.194")


if __name__ == '__main__':
    unittest.main()
