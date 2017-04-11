import unittest
import modeller
import os
import sys
import subprocess
import collections
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_ss"""
        for args in ([], ['x'] * 2):
            out = check_output(['allosmod', 'get_ss'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_ss'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_ss"""
        out = check_output(['allosmod', 'get_ss',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb')],
                           universal_newlines=True)
        self.assertEqual(out, '-\n-\n-\n-\n-\nS\nT\nT\n-\n-\n')

if __name__ == '__main__':
    unittest.main()
