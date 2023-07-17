import unittest
import subprocess
import os
import sys
import utils
from subprocess import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

records = [['HETATM', ' CA', 'HID'],
           ['HETATM', ' CA', 'HIE'],
           ['HETATM', ' CA', 'HIP'],
           ['HETATM', ' CA', 'HSD'],
           ['HETATM', ' CA', 'HSE'],
           ['HETATM', ' CA', 'HSP'],
           ['ATOM',   ' CA', 'HID'],
           ['ATOM',   ' CA', 'HIE'],
           ['ATOM',   ' CA', 'HIP'],
           ['ATOM',   ' CA', 'HSD'],
           ['ATOM',   ' CA', 'HSE'],
           ['ATOM',   ' CA', 'HSP'],
           ['HETATM', ' CA', 'MSE'],
           ['ATOM',   'SE ', 'MSE'],
           ['HETATM', 'SE ', 'MSE'],
           ['ATOM',   'SE ', 'MET'],
           ['HETATM', 'SE ', 'MET'],
           ['ATOM',   ' SE', 'MSE'],
           ['HETATM', ' SE', 'MSE'],
           ['HETATM', ' SE', 'MET']]
test_pdb = "EXPDTA    THEORETICAL MODEL, MODELLER SVN 2015/05/15 09:37:25\n" \
    + "\n".join("%-6s   %2d %-3s  %s    %2d      19.244  16.654   7.119"
                "  1.00 15.42           S"
                % (x[0], n+1, x[1], x[2], n+1) for n, x in enumerate(records))

exp_records = [['ATOM', ' CA', 'HIS']] * 12 \
              + [['ATOM', ' CA', 'MET'],
                 ['ATOM', ' SD', 'MET'],
                 ['ATOM', ' SD', 'MET'],
                 ['ATOM', 'SE ', 'MET'],
                 ['HETATM', 'SE ', 'MET'],
                 ['ATOM', ' SE', 'MET'],
                 ['ATOM', ' SE', 'MET'],
                 ['HETATM', ' SE', 'MET']]


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to pdb_fix_res"""
        for args in ([], ['foo', 'bar']):
            check_output(['allosmod', 'pdb_fix_res'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple test of pdb_fix_res"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write(test_pdb)

            def check_inplace():
                check_output(['allosmod', 'pdb_fix_res', '--in-place',
                              'test.pdb'], cwd=tmpdir)
                with open(os.path.join(tmpdir, 'test.pdb')) as fh:
                    return fh.read()
            for out in (check_output(['allosmod', 'pdb_fix_res', 'test.pdb'],
                                     universal_newlines=True, cwd=tmpdir),
                        check_output([sys.executable, '-m',
                                      'allosmod.pdb_fix_res', 'test.pdb'],
                                     universal_newlines=True, cwd=tmpdir),
                        check_inplace()):
                lines = out.split('\n')
                del lines[-1]
                for n, line in enumerate(lines):
                    self.assertEqual(line[:6].rstrip(' '), exp_records[n][0])
                    self.assertEqual(line[12:15], exp_records[n][1])
                    self.assertEqual(line[17:20], exp_records[n][2])


if __name__ == '__main__':
    unittest.main()
