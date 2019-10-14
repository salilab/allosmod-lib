import unittest
import modeller
import os
import sys
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_auto_align"""
        for args in ([], ['x'] * 5):
            out = check_output(['allosmod', 'get_auto_align'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output([sys.executable, '-m',
                                'allosmod.get_auto_align'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def setup_inputs(self, seq='AW'):
        with open('align.ali', 'w') as fh:
            fh.write("""C; Sample alignment
>P1;1fdx
sequence:1fdx:1    : :2   : ::: 2.00:-1.00
AY*
""")
        with open('templates', 'w') as fh:
            fh.write("5fd1\n")
        with open('5fd1', 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA A   1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  VAL A   2      26.593  16.867   8.258  1.00120.51           C
""")

    def test_simple(self):
        """Simple complete run of get_auto_align"""
        self.setup_inputs()
        check_output(['allosmod', 'get_auto_align', 'align.ali', '1fdx',
                      'templates', 'align_suggested.ali'])
        e = modeller.environ()
        a = modeller.alignment(e, file='align_suggested.ali')
        self.assertEqual(len(a), 2)
        self.assertEqual(a[0].code, '1fdx')
        self.assertEqual(a[1].code, '5fd1')
        for f in ('templates', '5fd1', 'align.ali', 'align_suggested.ali'):
            os.unlink(f)

if __name__ == '__main__':
    unittest.main()
