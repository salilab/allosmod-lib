import unittest
import subprocess
import os
import shutil
import modeller
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output
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
            fh.write(get_seq('1fdx', '1fdx.pdb', 'AF\nVV'))
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

    def setup_inputs(self, seq='AW', subdir=''):
        with open(os.path.join(subdir, 'test.aln'), 'w') as fh:
            fh.write(get_seq('1fdx', '1fdx', 'AY'))
            fh.write(get_seq('foo', 'foo', seq))
        with open('templates', 'w') as fh:
            fh.write("1fdx\n")
        with open(os.path.join(subdir, '1fdx'), 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA     1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  TYR     2      26.593  16.867   8.258  1.00120.51           C
""")

    def test_simple(self):
        """Simple complete run of get_pm_initialstruct"""
        from allosmod.get_pm_initialstruct import get_pm_initialstruct
        self.setup_inputs()
        if os.path.exists('pred_1fdx'):
            shutil.rmtree('pred_1fdx')

        check_output(['allosmod', 'get_pm_initialstruct', '--target', 'foo',
                      '--keep-alignment', 'test.aln', 'templates',
                      '.', '1', 'slow'])
        e = modeller.environ()
        m = modeller.model(e, file='pred_1fdx/foo.B99990001.pdb')
        self.assertEqual([x.code for x in m.residues], ['A', 'W'])
        self.assertEqual(m.chains[0].name, 'A')
        for f in ('1fdx', 'foo.B99990001.pdb', 'foo.ini', 'foo.sch',
                  'test.aln', 'foo.D00000001', 'foo.rsr',
                  'foo.V99990001'):
            os.unlink(os.path.join('pred_1fdx', f))
        os.rmdir('pred_1fdx')
        os.unlink('test.aln')
        os.unlink('templates')
        os.unlink('1fdx')

    def test_opts(self):
        """Complete run of get_pm_initialstruct using different options"""
        from allosmod.get_pm_initialstruct import get_pm_initialstruct
        self.setup_inputs(seq='A/W')

        if os.path.exists('pred_1fdx'):
            shutil.rmtree('pred_1fdx')
        os.mkdir('pred_1fdx')
        check_output(['allosmod', 'get_pm_initialstruct', '--target', 'foo',
                      '--restraints-only', 'test.aln', 'templates',
                      '.', '1', 'fast'])
        for f in ('1fdx', 'family.mat', 'foo.ini', 'test.aln', 'test.aln.ali',
                  'foo.rsr'):
            os.unlink(os.path.join('pred_1fdx', f))
        os.rmdir('pred_1fdx')
        os.unlink('test.aln')
        os.unlink('templates')
        os.unlink('1fdx')

    def test_nochdir(self):
        """Complete run of get_pm_initialstruct using --no-chdir"""
        from allosmod.get_pm_initialstruct import get_pm_initialstruct
        if os.path.exists('pred_1fdx'):
            shutil.rmtree('pred_1fdx')
        os.mkdir('pred_1fdx')
        self.setup_inputs(seq='A/W', subdir='pred_1fdx')

        check_output(['allosmod', 'get_pm_initialstruct', '--target', 'foo',
                      '--restraints-only', '--no-chdir', '--csrfile',
                      'test.rsr', 'test.aln', '../templates', '.', '1', 'fast'],
                     cwd='pred_1fdx')
        for f in ('1fdx', 'family.mat', 'foo.ini', 'test.aln', 'test.aln.ali'):
            os.unlink(os.path.join('pred_1fdx', f))
        os.rmdir('pred_1fdx')

if __name__ == '__main__':
    unittest.main()
