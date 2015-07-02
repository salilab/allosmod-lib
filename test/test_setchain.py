import unittest
import subprocess
import os
from allosmod.util import check_output

test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
"""

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to setchain"""
        for args in ([], ['foo', 'bar', 'baz']):
            out = check_output(['allosmod', 'setchain'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of setchain"""
        with open('test.pdb', 'w') as fh:
            fh.write(test_pdb)
        def check_inplace():
            check_output(['allosmod', 'setchain', '--in-place', 'test.pdb',
                          'X'])
            with open('test.pdb') as fh:
                return fh.read()
        for out in (check_output(['allosmod', 'setchain', 'test.pdb', 'X']),
                    check_output(['python', '-m', 'allosmod.setchain',
                                  'test.pdb', 'XYZ']),
                    check_inplace()):
            lines = out.split('\n')
            self.assertEqual(lines[2][17:25], 'CYS X   ')
        os.unlink('test.pdb')

if __name__ == '__main__':
    unittest.main()
