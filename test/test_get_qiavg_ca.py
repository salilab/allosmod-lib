import unittest
import modeller
import os
import sys
import subprocess
import collections
from test_pdb2ali import check_output

test_dir = os.path.dirname(sys.argv[0])

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
                                        'test_get_contacts.pdb'), '20.0',
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_1.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_2.pdb')])
        with open('qi_1.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents,
"""     1    0.4274
     2    1.0000
     3    0.8496
     4    0.8652
     5    0.8970
     6    0.9019
     7    0.8905
     8    0.9239
""")
        os.unlink('qi_1.dat')
        os.unlink('qi_2.dat')
        os.unlink('qi_avg.dat')

if __name__ == '__main__':
    unittest.main()
