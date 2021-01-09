import unittest
import modeller
import os
import sys
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import allosmod.util  # noqa: E402
from allosmod.util import check_output  # noqa: E402


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to make_mod_inputs"""
        for args in ([], ['x'] * 8):
            check_output(['allosmod', 'make_mod_inputs'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.make_mod_inputs'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def setup_inputs(self, seq='AW', dir='.'):
        def d(fname):
            return os.path.join(dir, fname)
        with open(d('align.ali'), 'w') as fh:
            fh.write("""C; Sample alignment
>P1;1fdx
sequence:1fdx:1    : :54   : ::: 2.00:-1.00
AY*
>P1;5fd1
structureX:5fd1:1    :A:2  :A::: 1.90: 0.19
AV*
""")
        with open(d('templates'), 'w') as fh:
            fh.write("5fd1\n")
        with open(d('5fd1'), 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA A   1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  VAL A   2      26.593  16.867   8.258  1.00120.51           C
""")
        with open(d('avgpdb.pdb'), 'w') as fh:
            fh.write("""
ATOM      2  CA  ALA A   1      27.449  14.935   5.140  1.00 29.87           C
ATOM      7  CA  TYR A   2      26.593  16.867   8.258  1.00120.51           C
""")

    def test_simple(self):
        """Simple complete run of make_mod_inputs"""
        with allosmod.util.temporary_directory() as tempdir:
            self.setup_inputs(dir=tempdir)
            check_output(['allosmod', 'make_mod_inputs', '--', '1fdx',
                          'templates', '-3333', '3', '3', '3', '4'],
                         cwd=tempdir)
            e = modeller.environ()
            for fname in ('random.ini', '1fdx.ini'):
                m = modeller.model(e, file=os.path.join(tempdir, fname))
                self.assertEqual([x.code for x in m.residues], ['A', 'Y'])
                # Should have converted CA-only in all-atom model
                self.assertEqual(len(m.atoms), 18)
            with open(os.path.join(tempdir, '1fdx.rsr')) as fh:
                self.assertEqual(len(fh.readlines()), 78)
            for f in ('templates', 'avgpdb.pdb', '5fd1', 'align.ali',
                      'random.ini', '1fdx.ini', '1fdx.rsr'):
                os.unlink(os.path.join(tempdir, f))

    def test_nucleic(self):
        """Test make_mod_inputs with nucleic acids"""
        import modeller.automodel
        from allosmod.make_mod_inputs import make_mod_inputs
        self.setup_inputs()

        def mock_make(cls, exit_stage):
            self.assertEqual(exit_stage, 1)
            self.assertAlmostEqual(cls.max_sc_sc_distance, 14.0, places=1)
        with utils.mock_method(modeller.automodel.automodel, 'make',
                               mock_make):
            make_mod_inputs('1fdx', ['5fd1'], -3333, [3, 3, 3], 4, True)
        for f in ('templates', 'avgpdb.pdb', '5fd1', 'align.ali',
                  'random.ini'):
            os.unlink(f)


if __name__ == '__main__':
    unittest.main()
