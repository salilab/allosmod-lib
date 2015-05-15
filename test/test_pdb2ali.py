import unittest
import subprocess
import os

def check_output(args, stderr=None, *other):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr, *other)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise OSError("Process %s exited with code %d"
                      % (" ".join(args), p.returncode))
    return stdout

test_pdb = """EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25
ATOM      1  N   CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      2  CA  CYS A   1       1.453   0.000   0.000  0.00  0.00           C
ATOM      8  CA  MET A   2       3.735   3.100   0.000  1.00  0.00           C
ATOM     16  CA  TYR A   3       7.437   4.095   0.000  1.00  0.00           C
TER      28      TYR A   3
HETATM   29 FE   HEM A   4       8.785   5.457   0.370  1.00  0.00          FE
HETATM   72 ZN    ZN A   5       2.149   2.021   2.548  1.00  0.00          ZN
HETATM   73  O   HOH A   6       2.323   1.794   2.131  1.00  0.00           O
ATOM     75  CA  CYS B   7       3.453   2.000   2.000  0.00  0.00           C
ATOM     81  CA  MET B   8       5.735   5.100   2.000  1.00  0.00           C
ATOM     89  CA  TYR B   9       9.437   6.095   2.000  1.00  0.00           C
TER     101      TYR B   9
HETATM  102 FE   HEM B  10      10.639   8.288   1.766  1.00  0.00          FE
HETATM  145 ZN    ZN B  11       4.120   4.754   4.681  1.00  0.00          ZN
HETATM  146  O   HOH B  12       4.370   4.293   4.222  1.00  0.00           O
END
"""

class Tests(unittest.TestCase):
    def test_simple(self):
        """Simple test of pdb2ali"""
        with open('test.pdb', 'w') as fh:
            fh.write(test_pdb)
        out = check_output(['allosmod', 'pdb2ali', 'test.pdb'])
        self.assertEqual(out, """>P1;test.pdb
structureX:test.pdb:   1 :A:+10:B:::-1.00:-1.00
CMYh./CMYh.*
""")
        os.unlink('test.pdb')

    def test_rewrite_chain(self):
        """Make sure that pdb2ali rewrites empty chain IDs"""
        with open('test.pdb', 'w') as fh:
            fh.write(test_pdb.replace(' A ', '   '))
        out = check_output(['allosmod', 'pdb2ali', 'test.pdb'])
        self.assertEqual(out, """>P1;test.pdb
structureX:test.pdb:   1 :@:+10:B:::-1.00:-1.00
CMYh./CMYh.*
""")
        with open('test.pdb') as fh:
            lines = fh.readlines()
        # Empty chain should have been reassigned as '@'
        self.assertEqual(lines[1][21], '@')
        os.unlink('test.pdb')

if __name__ == '__main__':
    unittest.main()
