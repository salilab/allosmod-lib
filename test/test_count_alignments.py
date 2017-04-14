import unittest
import modeller
import os
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to count_alignments"""
        for args in ([], ['x'] * 4):
            out = check_output(['allosmod', 'count_alignments'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.count_alignments'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def setup_inputs(self, tmpdir, seq='AW'):
        with open(os.path.join(tmpdir, 'align.ali'), 'w') as fh:
            fh.write("""C; Multiple alignment
>P1;5fd1
structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*
>P1;1bqx
structureN:1bqx:   1 :A: 77  :A:ferredoxin:Bacillus schlegelii:-1.00:-1.00
AYVITEPCIGTKCASCVEVCPVDCIHEGEDQYYIDPDVCIDCGACEAVCPVSAIYHEDFV-EEWKSYIQKNRDFF
KK-----------------------------*

>P1;pm.pdb
sequence:pm.pdb:1    : :54   : :ferredoxin:Peptococcus aerogenes: 2.00:-1.00
AY...VINDSC--IACGACKPECPVNIIQGS--IYAIDADSCIDCGSCASVCPVGAPNPED-----------------
-------------------------------*
""")

        with open(os.path.join(tmpdir, 'templates'), 'w') as fh:
            fh.write("5fd1\n1bqx\n")

    def test_simple(self):
        """Simple complete run of count_alignments"""
        with utils.temporary_directory() as tmpdir:
            self.setup_inputs(tmpdir)
            out = check_output(['allosmod', 'count_alignments',
                                'align.ali', 'templates', 'pm.pdb'],
                               universal_newlines=True, cwd=tmpdir)
            self.assertEqual(out.rstrip('\r\n'), "2")

    def test_count_alignments(self):
        """Test count_alignments function"""
        from allosmod.count_alignments import count_alignments
        with utils.temporary_directory() as tmpdir:
            self.setup_inputs(tmpdir)
            num_align, num_res = count_alignments(
                                      os.path.join(tmpdir, "align.ali"),
                                      "pm.pdb", ["5fd1", "1bqx"])
            self.assertEqual(num_align, 113)
            self.assertEqual(num_res, 57)

if __name__ == '__main__':
    unittest.main()
