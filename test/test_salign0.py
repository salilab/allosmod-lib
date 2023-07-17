import unittest
import subprocess
import os
import sys
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  CA  CYS A   1       1.453   2.100   3.200  0.00  0.00           C
ATOM      2  CA  CYS A   2       1.453   0.000   0.000  0.00  0.00           C
ATOM      3  CA  CYS A   3       3.735   3.100   0.000  1.00  0.00           C
ATOM      4  CA  CYS A   4       2.735   4.100   8.000  1.00  0.00           C
"""


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to salign0"""
        for args in ([], ['1', '2', '3', '4']):
            check_output(['allosmod', 'salign0'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.salign0'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_determine_fit_atoms(self):
        """Test determine_fit_atoms() function"""
        import allosmod.salign0

        class MockAtom(object):
            def __init__(self, name):
                self.name = name

        class MockModel(object):
            def __init__(self, names):
                self.atoms = [MockAtom(n) for n in names]
        m = MockModel(['N', 'CA', 'CA', 'CA', 'P', 'P', 'P'])
        self.assertEqual(allosmod.salign0.determine_fit_atoms(m), 'CA')
        m = MockModel(['CA', 'CA', 'P', 'P', 'P'])
        self.assertEqual(allosmod.salign0.determine_fit_atoms(m), 'P')
        m = MockModel([])
        self.assertEqual(allosmod.salign0.determine_fit_atoms(m), 'CA')

    def test_simple(self):
        """Simple test of salign0"""
        with open('test1.pdb', 'w') as fh:
            fh.write(test_pdb)
        with open('test2.pdb', 'w') as fh:
            fh.write(test_pdb)
        check_output(['allosmod', 'salign0', '-v',
                      'test1.pdb', 'test2.pdb', 'out.aln'])
        os.unlink('out.aln')
        os.unlink('test1_fit.pdb')
        os.unlink('test2_fit.pdb')
        os.unlink('test1.pdb')
        os.unlink('test2.pdb')


if __name__ == '__main__':
    unittest.main()
