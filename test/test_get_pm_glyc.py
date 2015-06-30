import os
import unittest
import subprocess
from test_pdb2ali import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_pm_glyc"""
        for args in ([], [''] * 7):
            out = check_output(['allosmod', 'get_pm_glyc'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_pm_glyc'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

if __name__ == '__main__':
    unittest.main()
