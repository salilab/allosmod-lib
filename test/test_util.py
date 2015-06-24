import unittest
import os
from io import StringIO, BytesIO
import allosmod.util

class Tests(unittest.TestCase):
    def test_read_templates(self):
        """Test read_templates function"""
        with open('lst', 'w') as fh:
            fh.write('test1 A 1 10\ntest2\ntest3 B 4 5\n')
        self.assertEqual(allosmod.util.read_templates('lst'),
                         ['test1', 'test2', 'test3'])
        os.unlink('lst')

    def test_fix_newlines(self):
        """Test fix_newlines function"""
        def assert_content(exp_content):
            with open('testnl') as fh:
                content = fh.read()
            self.assertEqual(content, exp_content)
        with open('testnl', 'w') as fh:
            fh.write("test")
        allosmod.util.fix_newlines("testnl")
        assert_content("test")
        with open('testnl', 'w') as fh:
            fh.write("test\r\nbar\r\n")
        allosmod.util.fix_newlines("testnl")
        assert_content("test\nbar\n")
        os.unlink('testnl')

    def test_sequence(self):
        """Test Sequence class"""
        s = allosmod.util.Sequence()
        s.primary = 'ACG/VT--W'
        self.assertEqual(s.get_residues(), 'ACGVTW')

    def test_pir_file_empty(self):
        """Test read of empty PIR file"""
        sio = StringIO("C; comment\nR; comment\n\n")
        p = allosmod.util.PIRFile()
        seqs = list(p.read(sio))
        self.assertEqual(len(seqs), 0)

    def test_pir_file_unterminated(self):
        """Test read of PIR file with unterminated sequence"""
        sio = StringIO(""">P1;template
structureX:::::::::
AFVV
>P1;seq
sequence:::::::::
AFVV*
""")
        p = allosmod.util.PIRFile()
        self.assertRaises(allosmod.util.FileFormatError, list, p.read(sio))

    def test_pir_file_ok(self):
        """Test read of OK PIR file"""
        sio = StringIO(""">P1;template
structureX:::::::::
A-
FVV*
>P1;seq
:::::::::
AF/VV*
""")
        p = allosmod.util.PIRFile()
        seqs = list(p.read(sio))
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[0].primary, 'A-FVV')
        self.assertEqual(seqs[1].primary, 'AF/VV')
        self.assertEqual(seqs[0].prottyp, 'structureX')
        self.assertEqual(seqs[1].prottyp, 'sequence')

    def test_pir_file_bad_header(self):
        """Test read of PIR file with bad header"""
        sio = StringIO(""">P1;template
garbage
A-FVV*
""")
        p = allosmod.util.PIRFile()
        self.assertRaises(allosmod.util.FileFormatError, list, p.read(sio))

    def test_pir_file_no_header(self):
        """Test read of PIR file with no header"""
        sio = StringIO("AFVV*")
        p = allosmod.util.PIRFile()
        self.assertRaises(allosmod.util.FileFormatError, list, p.read(sio))

    def test_pir_file_write(self):
        """Test write of PIR file"""
        sio = BytesIO()
        s = allosmod.util.Sequence()
        s.code = "testcode"
        s.primary = 'AFV'
        p = allosmod.util.PIRFile()
        p.write(sio, s, width=2)
        self.assertEqual(sio.getvalue(),
                         ">P1;testcode\nsequence:::::::::\nAF\nV\n*\n")

    def test_pdb_filters(self):
        """Test PDB filters"""
        self.assertTrue(allosmod.util.atom_filter("ATOM foo"))
        self.assertFalse(allosmod.util.atom_filter("HETATM foo"))
        self.assertFalse(allosmod.util.atom_filter("foo"))

        self.assertTrue(allosmod.util.atom_hetatm_filter("ATOM foo"))
        self.assertTrue(allosmod.util.atom_hetatm_filter("HETATM foo"))
        self.assertFalse(allosmod.util.atom_hetatm_filter("foo"))

if __name__ == '__main__':
    unittest.main()
