import unittest
import os
from io import StringIO
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import allosmod.util.align  # noqa: E402

TEST_ALIGNMENT = """
C; test alignment

>P1;pm.pdb
sequence:pm.pdb:  ::::::-1.00:-1.00
N/T--TVFQGVAGQSLQ*

>P1;TREMcopy.pdb
structureX:TREMcopy.pdb:  20 :A:+112 :A:::-1.00:-1.00
NT-TV/F-QGVAGQSLQ*
"""


class Tests(unittest.TestCase):
    def test_insert(self):
        """Test align.insert_gap()"""
        with utils.temporary_directory() as tmpdir:
            ali = os.path.join(tmpdir, 'test_insert.ali')
            with open(ali, 'w') as fh:
                fh.write(TEST_ALIGNMENT)
            s_out = StringIO()
            allosmod.util.align.insert_gap(ali, 1, 7, 9, s_out)
            out_lines = s_out.getvalue().split('\n')
            self.assertEqual(out_lines[2], 'N/T--TVFQG---VAGQSLQ')
            self.assertEqual(out_lines[6], 'NT-TV/F-QGVAG---QSLQ')

    def test_insert_fail(self):
        """Test align.insert_gap() failure"""
        with utils.temporary_directory() as tmpdir:
            ali = os.path.join(tmpdir, 'test_insert.ali')
            with open(ali, 'w') as fh:
                fh.write(TEST_ALIGNMENT)
            self.assertRaises(ValueError, allosmod.util.align.insert_gap,
                              ali, 1, 97, 99)


if __name__ == '__main__':
    unittest.main()
