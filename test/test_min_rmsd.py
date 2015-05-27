import unittest
import os
import subprocess
from test_pdb2ali import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to min_rmsd"""
        for args in ([], ['x'] * 3):
            out = check_output(['allosmod', 'min_rmsd'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.min_rmsd'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of min_rmsd"""
        with open('file1', 'w') as fh:
            fh.write("""
ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C  
ATOM     35  CA  VAL A   3      -0.964  -0.510   3.510  1.00  0.00           C  
ATOM     35  CA  VAL A   4      -4.964  -2.510   8.510  1.00  0.00           C  
""")
        with open('file2', 'w') as fh:
            fh.write("""
ATOM     14  CA  TYR A   2      -4.696  -0.319   4.300  1.00  0.00           C  
ATOM     35  CA  VAL A   3      -0.964  -0.510   3.510  1.00  0.00           C  
ATOM     35  CA  VAL A   4      -4.964  -2.510   8.510  1.00  0.00           C  
""")
        out = check_output(['allosmod', 'min_rmsd', 'file1', './file2'])
        self.assertEqual(out, "file1 file2 0.891\n")
        os.unlink("file1")
        os.unlink("file2")

if __name__ == '__main__':
    unittest.main()
