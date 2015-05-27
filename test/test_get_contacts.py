import unittest
import modeller
import os
import sys
import subprocess
import collections
from test_pdb2ali import check_output

test_dir = os.path.dirname(sys.argv[0])

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_contacts"""
        for args in ([], ['x'] * 3):
            out = check_output(['allosmod', 'get_contacts'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
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

    def test_get_contact_dist(self):
        """Test get_contact_dist()"""
        from allosmod.get_contacts import get_contact_dist
        class Residue(object):
            def __init__(self, av):
                self.average = av
        self.assertEqual(get_contact_dist(Residue(((0,0,0),)),
                                          Residue(((0,5,0),)), 16.0),
                         None)
        self.assertAlmostEqual(get_contact_dist(Residue(((0,0,0),)),
                                                Residue(((0,5,0),)), 100.),
                               5.0, places=1)

    def test_simple(self):
        """Simple complete run of get_contacts"""
        out = check_output(['allosmod', 'get_contacts',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'), '11.0'])
        lines = out.split('\n')
        self.assertEqual(len(lines), 12)
        self.assertEqual(lines[0][43:48], "9.955")
        self.assertEqual(lines[1][43:48], "9.194")

if __name__ == '__main__':
    unittest.main()
