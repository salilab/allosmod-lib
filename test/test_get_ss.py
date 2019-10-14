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
            out = check_output([sys.executable, '-m',
                                'allosmod.get_ss'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_ss"""
        out = check_output(['allosmod', 'get_ss',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb')],
                           universal_newlines=True)
        self.assertEqual(out, '-\n-\n-\n-\n-\nS\nT\nT\n-\n-\n')

    def test_rna(self):
        """Test get_ss with RNA-only file"""
        out = check_output(['allosmod', 'get_ss',
                           os.path.join(test_dir, 'input',
                                        'test_rna.pdb')],
                           universal_newlines=True)
        # No SS info for RNA, thus output should be empty
        self.assertEqual(out, '')

    def test_dssp_error(self):
        """Test handling of DSSP errors in get_ss"""
        out = check_output(['allosmod', 'get_ss',
                           os.path.join(test_dir, 'input',
                                        'not-exist.pdb')],
                           universal_newlines=True, retcode=1)
        self.assertEqual(out, '')

if __name__ == '__main__':
    unittest.main()
