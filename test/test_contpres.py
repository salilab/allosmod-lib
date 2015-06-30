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
        """Test wrong arguments to contpres"""
        for args in ([], [''] * 4):
            out = check_output(['allosmod', 'contpres'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.contpres'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of contpres"""
        with open('break.dat', 'w') as fh:
            fh.write("foo\nbar\n")
        out = check_output(['allosmod', 'contpres',
                            os.path.join(test_dir, 'input',
                                         'test_contpres.rsr'),
                            os.path.join(test_dir, 'input',
                                         'test_contpres.pdb'), '50.0'])
        with open('break.dat') as fh:
            content = fh.read()
        self.assertEqual(content, "foo\nbar\n1 50.0\n3 50.0\n")
        with open('contpres.dat') as fh:
            content = fh.read()
        self.assertEqual(content, "1 4\n2 0\n3 4\n")
        os.unlink('break.dat')
        os.unlink('contpres.dat')

    def test_cdensity(self):
        """Simple contpres in cdensity mode"""
        with open('break.dat', 'w') as fh:
            fh.write("foo\nbar\n")
        out = check_output(['allosmod', 'contpres', '--cdensity_cutoff', '1.0',
                            os.path.join(test_dir, 'input',
                                         'test_contpres.rsr'),
                            os.path.join(test_dir, 'input',
                                         'test_contpres.pdb'), '50.0'])
        with open('break.dat') as fh:
            content = fh.read()
        self.assertEqual(content, "foo\nbar\n2 50.0\n")
        os.unlink('break.dat')
        os.unlink('contpres.dat')

if __name__ == '__main__':
    unittest.main()
