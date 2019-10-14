import os
import sys
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

ALIGN_ALI = """
>P1;pdb1
structureX:asite_pdb1:   1 :A:+30  :A:::-1.00:-1.00
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPN*

>P1;pdb2
structureX:asite_pdb2:   1 :A:+30  :A:::-1.00:-1.00
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPN*

>P1;pm.pdb
structureX:pm.pdb:   1 :A:+30  :A:::-1.00:-1.00
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPN*
"""

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to getavgpdb"""
        for args in ([], [''] * 5):
            out = check_output(['allosmod', 'getavgpdb'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output([sys.executable, '-m',
                                'allosmod.getavgpdb'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of getavgpdb"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'align.ali'), 'w') as fh:
                fh.write(ALIGN_ALI)
            out = check_output(['allosmod', 'getavgpdb',
                               os.path.join(test_dir, 'input',
                                            'asite_pdb1.pdb'),
                               os.path.join(test_dir, 'input',
                                            'asite_pdb2.pdb'),
                               'pdb1', 'pdb2'], cwd=tmpdir)
            self.assertEqual(sorted(os.listdir(tmpdir)),
                             ['align.ali', 'avgpdb.pdb', 'list', 'run.log'])

if __name__ == '__main__':
    unittest.main()
