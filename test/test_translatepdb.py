import unittest
import subprocess
import os
import sys
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   2.100   3.200  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
"""


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to translatepdb"""
        for args in ([], ['foo', 'bar', 'baz', 'foo', 'bar']):
            check_output(['allosmod', 'translatepdb'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of translatepdb"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write(test_pdb)
            for out in (check_output(['allosmod', 'translatepdb', '--',
                                      'test.pdb', '1.2', '-2.3', '3.4'],
                                     universal_newlines=True, cwd=tmpdir),
                        check_output([sys.executable, '-m',
                                      'allosmod.translatepdb',
                                      '--', 'test.pdb', '1.2', '-2.3', '3.4'],
                                     universal_newlines=True, cwd=tmpdir)):
                lines = out.split('\n')
                self.assertEqual(len(lines), 4)
                x = float(lines[0][30:38])
                y = float(lines[0][38:46])
                z = float(lines[0][46:54])
                self.assertAlmostEqual(x, 2.653, places=3)
                self.assertAlmostEqual(y, -0.200, places=3)
                self.assertAlmostEqual(z, 6.600, places=3)
            os.unlink(os.path.join(tmpdir, 'test.pdb'))


if __name__ == '__main__':
    unittest.main()
