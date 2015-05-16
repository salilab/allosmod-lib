import unittest
import subprocess
import os
from test_pdb2ali import check_output

test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   2.100   3.200  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
"""

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to getcofm"""
        for args in ([], ['foo', 'bar']):
            out = check_output(['allosmod', 'getcofm'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of getcofm"""
        with open('test.pdb', 'w') as fh:
            fh.write(test_pdb)
        for out in (check_output(['allosmod', 'getcofm', 'test.pdb']),
                    check_output(['python', '-m', 'allosmod.getcofm',
                                  'test.pdb'])):
            x = float(out[:8])
            y = float(out[9:17])
            z = float(out[18:26])
            self.assertAlmostEqual(x, 2.214, places=3)
            self.assertAlmostEqual(y, 1.733, places=3)
            self.assertAlmostEqual(z, 1.067, places=3)
        os.unlink('test.pdb')

if __name__ == '__main__':
    unittest.main()
