import unittest
import subprocess
from test_pdb2ali import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to edit_restraints"""
        for args in ([], ['' * 6]):
            out = check_output(['allosmod', 'edit_restraints'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.edit_restraints'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

if __name__ == '__main__':
    unittest.main()
