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
        """Test wrong arguments to get_q_ca"""
        out = check_output(['allosmod', 'get_q_ca'],
                           stderr=subprocess.STDOUT, retcode=2)
        out = check_output(['python', '-m',
                            'allosmod.get_q_ca'],
                           stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_q_ca"""
        out = check_output(['allosmod', 'get_q_ca',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'), '11.0',
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_1.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_qiavg_2.pdb')])
        with open('qscore1to8.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents,
"""     1    0.8364     0.8706     0.7509     0.0000
     2    0.7251     0.8039     0.5280     0.0000
""")
        with open('qs_cut1to8.dat') as fh:
            contents = fh.read()
        self.assertEqual(contents,
"""     1    0.9373     0.9373     0.0000     0.0000
     2    0.9169     0.9169     0.0000     0.0000
""")
        os.unlink('qscore1to8.dat')
        os.unlink('qs_cut1to8.dat')

if __name__ == '__main__':
    unittest.main()
