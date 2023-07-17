import unittest
import os
import sys
import subprocess
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_qiavg_ca"""
        check_output(['allosmod', 'get_qiavg_ca'],
                     stderr=subprocess.STDOUT, retcode=2)
        check_output([sys.executable, '-m',
                      'allosmod.get_qiavg_ca'],
                     stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_qiavg_ca"""
        check_output(['allosmod', 'get_qiavg_ca',
                     os.path.join(test_dir, 'input',
                                  'test_qiavg.pdb'), '20.0',
                     os.path.join(test_dir, 'input',
                                  'test_qiavg_1.pdb'),
                     os.path.join(test_dir, 'input',
                                  'test_qiavg_2.pdb')])
        with open('qi_1.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents, """     1    0.4274
     2    1.1000
     3    0.8496
     4    0.8315
     5    0.8712
     6    0.8774
     7    0.8632
     8    0.9087
     9    1.1000
""")
        os.unlink('qi_1.dat')
        os.unlink('qi_2.dat')
        os.unlink('qi_avg.dat')

    def test_nres_mismatch(self):
        """Test get_qiavg_ca with different numbers of residues"""
        out = check_output(['allosmod', 'get_qiavg_ca',
                           os.path.join(test_dir, 'input',
                                        'test_qiavg.pdb'), '20.0',
                           os.path.join(test_dir, 'input',
                                        'test_rna.pdb')],
                           stderr=subprocess.STDOUT, retcode=1,
                           universal_newlines=True)
        self.assertTrue('different numbers of residues' in out)


if __name__ == '__main__':
    unittest.main()
