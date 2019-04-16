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
        """Test wrong arguments to get_qiavg_ca"""
        out = check_output(['allosmod', 'get_qiavg_ca'],
                           stderr=subprocess.STDOUT, retcode=2)
        out = check_output(['python', '-m',
                            'allosmod.get_qiavg_ca'],
                           stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_qiavg_ca"""
        out = check_output(['allosmod', 'get_qiavg_ca',
                           os.path.join(test_dir, 'input',
                                        'test_qiavg.pdb'), '20.0',
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_1.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_2.pdb')])
        with open('qi_1.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents,
"""     1    0.4274
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
                           stderr=subprocess.STDOUT, retcode=1)
        self.assertTrue('different numbers of residues' in out)

if __name__ == '__main__':
    unittest.main()
