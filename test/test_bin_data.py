import unittest
import subprocess
import os
import sys
from test_pdb2ali import check_output

test_dir = os.path.dirname(sys.argv[0])

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to bin_data"""
        for args in ([], ['x'] * 6):
            out = check_output(['allosmod', 'bin_data'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.bin_data'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of bin_data"""
        with open('test.data', 'w') as fh:
            fh.write("""1.0 10
1.1 10
1.2 10
1.2 10
2.0\t20
3.0 30
4.0 40
5.0 50
6.0 60
""")
        out = check_output(['allosmod', 'bin_data', 'test.data', '0', 'm1',
                            '10', '20'])
        self.assertEqual(out, """0.92500 0.00000 4.00000
2.02500 0.00000 1.00000
3.12500 0.00000 1.00000
4.22500 0.00000 1.00000
4.77500 0.00000 1.00000
5.87500 0.00000 1.00000
""")
        out = check_output(['allosmod', 'bin_data', 'test.data', '1', 'm1',
                            '10', '20'])
        self.assertEqual(out, "9.72500 1.00000 4.00000\n")
        os.unlink('test.data')

if __name__ == '__main__':
    unittest.main()