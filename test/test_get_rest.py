import os
import unittest
import subprocess
from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_rest"""
        for args in ([], [''] * 2):
            out = check_output(['allosmod', 'get_rest'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_rest'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

if __name__ == '__main__':
    unittest.main()
