import unittest
import os
import sys
import subprocess
import utils
from subprocess import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to min_rmsd"""
        for args in ([], ['x'] * 3):
            check_output(['allosmod', 'min_rmsd'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.min_rmsd'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of min_rmsd"""
        with open('file1', 'w') as fh:
            fh.write("""
ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C  
ATOM     35  CA  VAL A   3      -0.964  -0.510   3.510  1.00  0.00           C  
ATOM     35  CA  VAL A   4      -4.964  -2.510   8.510  1.00  0.00           C  
""")   # noqa: W291
        with open('file2', 'w') as fh:
            fh.write("""
ATOM     14  CA  TYR A   2      -4.696  -0.319   4.300  1.00  0.00           C  
ATOM     35  CA  VAL A   3      -0.964  -0.510   3.510  1.00  0.00           C  
ATOM     35  CA  VAL A   4      -4.964  -2.510   8.510  1.00  0.00           C  
""")   # noqa: W291
        out = check_output(['allosmod', 'min_rmsd', 'file1', './file2'],
                           universal_newlines=True)
        self.assertEqual(out, "file1 file2 0.891\n")
        os.unlink("file1")
        os.unlink("file2")


if __name__ == '__main__':
    unittest.main()
