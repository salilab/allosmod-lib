import unittest
import subprocess
import os
import sys
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to bin_data"""
        for args in ([], ['x'] * 6):
            check_output(['allosmod', 'bin_data'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.bin_data'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of bin_data"""
        with utils.temporary_directory() as tmpdir:
            testdata = os.path.join(tmpdir, 'test.data')
            with open(testdata, 'w') as fh:
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
                                '10', '20'],
                               universal_newlines=True, cwd=tmpdir)
            self.assertEqual(out, """0.92500 0.44444 4.00000
2.02500 0.11111 1.00000
3.12500 0.11111 1.00000
4.22500 0.11111 1.00000
4.77500 0.11111 1.00000
5.87500 0.11111 1.00000
""")
            out = check_output(['allosmod', 'bin_data', 'test.data', '1', 'm1',
                                '10', '20'],
                               universal_newlines=True, cwd=tmpdir)
            self.assertEqual(out, "9.72500 1.00000 4.00000\n")


if __name__ == '__main__':
    unittest.main()
