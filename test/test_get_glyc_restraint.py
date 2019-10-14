import os
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_glyc_restraint"""
        for args in ([], [''] * 3):
            out = check_output(['allosmod', 'get_glyc_restraint'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output([sys.executable, '-m',
                                'allosmod.get_glyc_restraint'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
    def test_simple(self):
        """Simple complete run of get_glyc_restraint"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test-allosmod.py'), 'w') as fh:
                fh.write("""line 1
def special_patches(self, aln):
    self.patch(residue_type='NGLB', residues=(self.residues['1:A'],
                                              self.residues['2:B']))
line 2
""")
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write("""
ATOM      1  N   ALA A   1      17.807  17.608   5.019  1.00 17.18           N
ATOM      2  CA  ALA A   1      17.121  17.162   6.197  1.00 15.60           C
ATOM      3  C   ALA A   1      18.085  17.018   7.343  1.00 14.54           C
ATOM      4  O   ALA A   1      19.244  16.654   7.119  1.00 15.42           O
HETATM    5  O1  BMA B   2      19.244  16.654   7.119  1.00 15.42           O
HETATM    6  C1  BMA B   2      19.244  16.654   7.119  1.00 15.42           O
""")
            out = check_output(['allosmod', 'get_glyc_restraint', 'test.pdb',
                                'test-allosmod.py'], stderr=subprocess.STDOUT,
                               retcode=0, universal_newlines=True, cwd=tmpdir)
            # One restraint between CA and C1:
            self.assertEqual(out, 'R    3   1   1   1   2   2   1        '
                                  '2     6    5.0000    0.0350\n')
            # special_patches should have been removed from .py file:
            with open(os.path.join(tmpdir, 'test-allosmod.py')) as fh:
                contents = fh.read()
            self.assertEqual(contents, 'line 1\nline 2\n')

    def test_bad_residue(self):
        """Simple get_glyc_restraint with invalid residue"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test-allosmod.py'), 'w') as fh:
                fh.write("""line 1
def special_patches(self, aln):
    self.patch(residue_type='NGLB', residues=(self.residues['1:A'],
                                              self.residues['2:B']))
line 2
""")
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write("""
ATOM      1  N   ALA A   1      17.807  17.608   5.019  1.00 17.18           N
ATOM      2  CA  ALA A   1      17.121  17.162   6.197  1.00 15.60           C
ATOM      3  C   ALA A   1      18.085  17.018   7.343  1.00 14.54           C
ATOM      4  O   ALA A   1      19.244  16.654   7.119  1.00 15.42           O
""")
            out = check_output(['allosmod', 'get_glyc_restraint', 'test.pdb',
                                'test-allosmod.py'], stderr=subprocess.STDOUT,
                               retcode=1, cwd=tmpdir)

if __name__ == '__main__':
    unittest.main()
