import unittest
import subprocess
import os
import sys
import modeller
import utils
from subprocess import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


def get_seq(code, pdb, seq):
    return """>P1;%s
structureX:%s:1::54::::2.00:-1.00
%s*
""" % (code, pdb, seq)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_pm_initialstruct"""
        for args in ([], ['1', '2', '3', '4', '5', '6']):
            check_output(['allosmod', 'get_pm_initialstruct'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
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

    def setup_inputs(self, topdir, seq='AW', subdir=''):
        with open(os.path.join(topdir, subdir, 'test.aln'), 'w') as fh:
            fh.write(get_seq('1fdx', '1fdx', 'AY'))
            fh.write(get_seq('foo', 'foo', seq))
        with open(os.path.join(topdir, 'templates'), 'w') as fh:
            fh.write("1fdx\n")
        with open(os.path.join(topdir, subdir, '1fdx'), 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA     1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  TYR     2      26.593  16.867   8.258  1.00120.51           C
""")

    def test_simple(self):
        """Simple complete run of get_pm_initialstruct"""
        with utils.temporary_directory() as tmpdir:
            self.setup_inputs(tmpdir)

            check_output(['allosmod', 'get_pm_initialstruct',
                          '--target', 'foo',
                          '--keep-alignment', 'test.aln', 'templates',
                          '.', '1', 'slow'], cwd=tmpdir)
            e = modeller.environ()
            m = modeller.model(e, file=os.path.join(tmpdir, 'pred_1fdx',
                               'foo.B99990001.pdb'))
            self.assertEqual([x.code for x in m.residues], ['A', 'W'])
            self.assertEqual(m.chains[0].name, 'A')
            for f in ('1fdx', 'foo.B99990001.pdb', 'foo.ini', 'foo.sch',
                      'test.aln', 'foo.D00000001', 'foo.rsr',
                      'foo.V99990001'):
                os.unlink(os.path.join(tmpdir, 'pred_1fdx', f))

    def test_opts(self):
        """Complete run of get_pm_initialstruct using different options"""
        with utils.temporary_directory() as tmpdir:
            self.setup_inputs(tmpdir, seq='A/W')

            os.mkdir(os.path.join(tmpdir, 'pred_1fdx'))
            check_output(['allosmod', 'get_pm_initialstruct',
                          '--target', 'foo',
                          '--restraints-only', 'test.aln', 'templates',
                          '.', '1', 'fast'], cwd=tmpdir)
            for f in ('1fdx', 'family.mat', 'foo.ini', 'test.aln',
                      'test.aln.ali', 'foo.rsr'):
                os.unlink(os.path.join(tmpdir, 'pred_1fdx', f))

    def test_nochdir(self):
        """Complete run of get_pm_initialstruct using --no-chdir"""
        with utils.temporary_directory() as tmpdir:
            os.mkdir(os.path.join(tmpdir, 'pred_1fdx'))
            self.setup_inputs(tmpdir, seq='A/W', subdir='pred_1fdx')

            check_output(['allosmod', 'get_pm_initialstruct',
                          '--target', 'foo',
                          '--restraints-only', '--no-chdir', '--csrfile',
                          'test.rsr', 'test.aln', '../templates', '.', '1',
                          'fast'], cwd=os.path.join(tmpdir, 'pred_1fdx'))
            for f in ('1fdx', 'family.mat', 'foo.ini', 'test.aln',
                      'test.aln.ali'):
                os.unlink(os.path.join(tmpdir, 'pred_1fdx', f))


if __name__ == '__main__':
    unittest.main()
