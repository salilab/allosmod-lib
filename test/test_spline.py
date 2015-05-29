import unittest
import subprocess
import os
import sys
from test_pdb2ali import check_output

test_dir = os.path.dirname(sys.argv[0])

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to spline"""
        for args in ([], ['x'] * 4):
            out = check_output(['allosmod', 'spline'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.spline'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of spline"""
        with open('test.pdb', 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA A   1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  VAL A   2      26.593  16.867   8.258  1.00120.51           C
""")
        check_output(['allosmod', 'spline', 'test.pdb',
                      os.path.join(test_dir, 'input', 'edited.rsr'), 'out.rsr'])
        rsr = open('out.rsr').read()
        expected = open(os.path.join(test_dir, 'input', 'converted.rsr')).read()
        self.assertEqual(rsr, expected)
        for f in ('test.pdb', 'out.rsr'):
            os.unlink(f)

if __name__ == '__main__':
    unittest.main()
