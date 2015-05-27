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
        """Test wrong arguments to get_inter_contacts"""
        for args in ([], ['x'] * 4):
            out = check_output(['allosmod', 'get_inter_contacts'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_inter_contacts'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_inter_contacts"""
        out = check_output(['allosmod', 'get_inter_contacts',
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'),
                           os.path.join(test_dir, 'input',
                                        'test_get_contacts.pdb'), '8.0'])
        lines = out.split('\n')
        self.assertEqual(len(lines), 44)
        self.assertEqual(lines[0][43:48], "0.000")
        self.assertEqual(lines[1][43:48], "7.029")
        self.assertEqual(lines[-5], "      10   A      10   A  HEM  HEM  0      0.000  6  6")

if __name__ == '__main__':
    unittest.main()
