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
        """Test wrong arguments to get_qmatrix"""
        out = check_output(['allosmod', 'get_qmatrix'],
                           stderr=subprocess.STDOUT, retcode=2)
        out = check_output(['python', '-m',
                            'allosmod.get_qmatrix'],
                           stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_qmatrix"""
        out = check_output(['allosmod', 'get_qmatrix',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'), '11.0',
                           os.path.join(test_dir, 'input',
                                        'test_qmatrix_1.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_qmatrix_2.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_qmatrix_3.pdb')])
        # matrix is randomized, so sort it
        inds = {}
        mat = {}
        with open('qmatrix.dat') as fh:
            for n, line in enumerate(fh):
                fname, m1, m2, m3 = line.rstrip('\r\n').split()
                m1 = float(m1)
                m2 = float(m2)
                m3 = float(m3)
                mat_ind = int(fname[-5])-1
                inds[mat_ind] = n
                mat[mat_ind] = (m1,m2,m3)
        for i in range(3):
            self.assertAlmostEqual(mat[i][inds[i]], 1.0, places=3)
        self.assertAlmostEqual(mat[0][inds[1]], 0.9278, places=4)
        self.assertAlmostEqual(mat[1][inds[0]], 0.9278, places=4)
        self.assertAlmostEqual(mat[2][inds[0]], 0.8815, places=4)
       
        with open('cq_aq_qavg_qsd.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents, "Qa,b: 0.91\n")
        os.unlink('qmatrix.dat')
        os.unlink('cq_aq_qavg_qsd.dat')

if __name__ == '__main__':
    unittest.main()
