import unittest
import subprocess
import os
from allosmod.util import check_output

test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   2.100   3.200  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
"""

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to rotatepdb"""
        for args in ([], ['foo', 'bar', 'baz', 'foo', 'bar']):
            out = check_output(['allosmod', 'rotatepdb'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of rotatepdb"""
        with open('test.pdb', 'w') as fh:
            fh.write(test_pdb)
        for out in (check_output(['allosmod', 'rotatepdb', '--', 'test.pdb',
                                  '10', '-20', '30']),
                    check_output(['python', '-m', 'allosmod.rotatepdb',
                                  '--', 'test.pdb', '10', '-20', '30'])):
            lines = out.split('\n')
            self.assertEqual(len(lines), 4)
            x = float(lines[0][30:38])
            y = float(lines[0][38:46])
            z = float(lines[0][46:54])
            self.assertAlmostEqual(x, -0.615, places=3)
            self.assertAlmostEqual(y, 1.391, places=3)
            self.assertAlmostEqual(z, 3.801, places=3)
        os.unlink('test.pdb')

if __name__ == '__main__':
    unittest.main()
