import unittest
import subprocess
import os
import sys
import utils
from subprocess import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
"""


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to setchain"""
        for args in ([], ['foo', 'bar', 'baz']):
            check_output(['allosmod', 'setchain'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of setchain"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write(test_pdb)

            def check_inplace():
                check_output(['allosmod', 'setchain', '--in-place', 'test.pdb',
                              'X'], cwd=tmpdir)
                with open(os.path.join(tmpdir, 'test.pdb')) as fh:
                    return fh.read()
            for out in (check_output(['allosmod', 'setchain', 'test.pdb', 'X'],
                                     universal_newlines=True, cwd=tmpdir),
                        check_output([sys.executable, '-m',
                                      'allosmod.setchain', 'test.pdb', 'XYZ'],
                                     universal_newlines=True, cwd=tmpdir),
                        check_inplace()):
                lines = out.split('\n')
                self.assertEqual(lines[2][17:25], 'CYS X   ')


if __name__ == '__main__':
    unittest.main()
