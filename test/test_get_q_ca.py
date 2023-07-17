import unittest
import os
import sys
import subprocess
import utils
from subprocess import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_q_ca"""
        check_output(['allosmod', 'get_q_ca'],
                     stderr=subprocess.STDOUT, retcode=2)
        check_output([sys.executable, '-m',
                      'allosmod.get_q_ca'],
                     stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_q_ca"""
        check_output(['allosmod', 'get_q_ca',
                     os.path.join(test_dir, 'input',
                                  'test_qiavg.pdb'), '11.0',
                     os.path.join(test_dir, 'input',
                                  'test_qiavg_1.pdb'),
                     os.path.join(test_dir, 'input',
                                  'test_qiavg_2.pdb')])
        with open('qscore1to9.dat') as fh:
            contents = fh.read()
        self.assertEqual(
            contents,
            "     1    0.8314     0.8706     0.7474     0.0000\n"
            "     2    0.7676     0.8366     0.6433     0.0000\n")
        with open('qs_cut1to9.dat') as fh:
            contents = fh.read()
        self.assertEqual(
            contents,
            "     1    0.9248     0.9248     0.0000     0.0000\n"
            "     2    0.9169     0.9169     0.0000     0.0000\n")
        os.unlink('qscore1to9.dat')
        os.unlink('qs_cut1to9.dat')


if __name__ == '__main__':
    unittest.main()
