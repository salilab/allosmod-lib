import unittest
import os
from io import BytesIO
import allosmod.util.align

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
        with open('test_insert.ali', 'w') as fh:
            fh.write(TEST_ALIGNMENT)
        s_out = BytesIO()
        allosmod.util.align.insert_gap('test_insert.ali', 1, 7, 9, s_out)
        out_lines = s_out.getvalue().split('\n')
        self.assertEqual(out_lines[2], 'N/T--TVFQG---VAGQSLQ')
        self.assertEqual(out_lines[6], 'NT-TV/F-QGVAG---QSLQ')
        os.unlink('test_insert.ali')

if __name__ == '__main__':
    unittest.main()
