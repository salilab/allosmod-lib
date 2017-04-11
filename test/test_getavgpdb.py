import os
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to getavgpdb"""
        for args in ([], [''] * 5):
            out = check_output(['allosmod', 'getavgpdb'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.getavgpdb'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

if __name__ == '__main__':
    unittest.main()
