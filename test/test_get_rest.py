import os
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_rest"""
        for args in ([], [''] * 2):
            out = check_output(['allosmod', 'get_rest'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output([sys.executable, '-m',
                                'allosmod.get_rest'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_rest"""
        with open('get_rest.in', 'w') as fh:
            pass
        out = check_output(['allosmod', 'get_rest',
                           os.path.join(test_dir, 'input',
                                        'asite_pdb1.pdb')],
                           universal_newlines=True)
        os.unlink('get_rest.in')
        # PDB file contains no sugars, so no restraints should be output
        self.assertEqual(out, '')


if __name__ == '__main__':
    unittest.main()
