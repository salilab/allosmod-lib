import unittest
import subprocess
import os
import modeller
from test_pdb2ali import check_output
from test_modeller import mock_method

def get_seq(code, pdb, seq):
    return """>P1;%s
structureX:%s:1::54::::2.00:-1.00
%s*
""" % (code, pdb, seq)

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_pm_initialstruct"""
        for args in ([], ['1', '2', '3', '4', '5', '6']):
            out = check_output(['allosmod', 'get_pm_initialstruct'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_pm_initialstruct'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_read_templates(self):
        """Test read_templates function"""
        import allosmod.get_pm_initialstruct
        with open('lst', 'w') as fh:
            fh.write('test1\ntest2\n')
        self.assertEqual(allosmod.get_pm_initialstruct.read_templates('lst'),
                         ['test1', 'test2'])
        os.unlink('lst')

    def test_get_target(self):
        """Test get_target function"""
        from allosmod.get_pm_initialstruct import get_target
        with open('test.aln', 'w') as fh:
            fh.write("""C; Sample alignment
>P1;1fdx
sequence:1fdx.pdb:1    : :54   : ::: 2.00:-1.00
AFVV*
>P1;5fd1
structureX:5fd1.pdb:1    :A:106  :A::: 1.90: 0.19
AFVV*
""")
        e = modeller.environ()
        self.assertEqual(get_target(e, 'foo', 'test.aln'), 'foo')
        self.assertEqual(get_target(e, None, 'test.aln'), '1fdx')
        os.unlink('test.aln')

    def test_find_het(self):
        """Test find_het function"""
        from allosmod.get_pm_initialstruct import find_het
        with open('test.aln', 'w') as fh:
            fh.write("C; Sample alignment\n")
            fh.write(get_seq('1fdx', '1fdx.pdb', 'AFVV'))
            fh.write(get_seq('5fd1', '5df1.pdb', 'AF.VV'))
            fh.write(get_seq('foo', 'foo.pdb', 'AF.VV'))
        self.assertEqual(find_het('test.aln', ['1fdx', 'foo']),
                         {'1fdx':None, 'foo': True})
        os.unlink('test.aln')

    def test_remove_hetatms(self):
        """Test remove_hetatms_from_aln_file function"""
        from allosmod.get_pm_initialstruct import remove_hetatms_from_aln_file
        orig_aln = "C; Sample alignment\n" \
                   + get_seq('1fdx', '1fdx.pdb', 'AFVV') \
                   + get_seq('5fd1', '5df1.pdb', 'AF.VV') \
                   + get_seq('foo', 'foo.pdb', 'AF.VV')
        nohet_aln = "C; Sample alignment\n" \
                    + get_seq('1fdx', '1fdx.pdb', 'AFVV') \
                    + get_seq('5fd1', '5df1.pdb', 'AFVV') \
                    + get_seq('foo', 'foo.pdb', 'AFVV')
        def make_aln():
            with open('test.aln', 'w') as fh:
                fh.write(orig_aln)
        def test_aln(expected):
            with open('test.aln') as fh:
                actual = fh.read()
            self.assertEqual(actual, expected)
        make_aln()
        remove_hetatms_from_aln_file('test.aln', 'foo', '1fdx')
        test_aln(nohet_aln)
        make_aln()
        remove_hetatms_from_aln_file('test.aln', '1fdx', 'foo')
        test_aln(orig_aln)
        os.unlink('test.aln')

    def test_simple(self):
        """Simple complete run of get_pm_initialstruct"""
        from allosmod.get_pm_initialstruct import get_pm_initialstruct
        with open('test.aln', 'w') as fh:
            fh.write(get_seq('1fdx', '1fdx', 'AY'))
            fh.write(get_seq('foo', 'foo', 'AW'))
        with open('templates', 'w') as fh:
            fh.write("1fdx\n")
        with open('1fdx', 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA     1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  TYR     2      26.593  16.867   8.258  1.00120.51           C
""")
        check_output(['allosmod', 'get_pm_initialstruct', '--target', 'foo',
                      '--keep-alignment', 'test.aln', 'templates',
                      '.', '1', 'slow'])
        e = modeller.environ()
        m = modeller.model(e, file='pred_1fdx/foo.B99990001.pdb')
        self.assertEqual([x.code for x in m.residues], ['A', 'W'])
        for f in ('1fdx', 'foo.B99990001.pdb', 'foo.ini', 'foo.sch',
                  'test.aln', 'foo.D00000001', 'foo.rsr',
                  'foo.V99990001'):
            os.unlink(os.path.join('pred_1fdx', f))
        os.rmdir('pred_1fdx')
        os.unlink('test.aln')
        os.unlink('templates')
        os.unlink('1fdx')

if __name__ == '__main__':
    unittest.main()
