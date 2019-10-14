import os
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import allosmod.util
from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_loopadjres"""
        for args in ([''],):
            out = check_output(['allosmod', 'get_loopadjres'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output([sys.executable, '-m',
                                'allosmod.get_loopadjres'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_no_loop(self):
        """Test run with no loop"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'get_loopadjres'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=0)
            with open(os.path.join(tempdir, 'break.dat')) as fh:
                contents = fh.read()
            self.assertEqual(contents, '')

    def test_one_loop(self):
        """Test run with one loop"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "break.dat"), 'w') as fh:
                fh.write('line1\nline2\n')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFVV----------AFVV*
>P1;pm.pdb
sequence:::::::::
AFVVCCCCCCCCCCAFVV*""")
            out = check_output(['allosmod', 'get_loopadjres'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=0)
            with open(os.path.join(tempdir, 'break.dat')) as fh:
                contents = fh.read()
            self.assertEqual(contents, 'line1\nline2\n1 2.0\n2 2.0\n3 2.0\n'
                                '4 2.0\n15 2.0\n16 2.0\n17 2.0\n18 2.0\n')

if __name__ == '__main__':
    unittest.main()
