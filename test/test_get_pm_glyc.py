import os
import sys
import string
import shutil
import unittest
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)

import allosmod.util
import allosmod.get_pm_glyc
from allosmod.util import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_pm_glyc"""
        for args in ([], [''] * 7):
            out = check_output(['allosmod', 'get_pm_glyc'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_pm_glyc'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of get_pm_glyc script"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, 'list'), 'w') as fh:
                fh.write('test_pm_glyc.pdb')
            shutil.copy(os.path.join(test_dir, 'input', 'asite_pdb1.pdb'),
                        os.path.join(tempdir, 'test_pm_glyc.pdb'))
            shutil.copy(os.path.join(test_dir, 'input', 'test_pm_glyc.ali'),
                        os.path.join(tempdir, 'align.ali'))
            shutil.copy(os.path.join(test_dir, 'input', 'glyc.dat'),
                        tempdir)
            out = check_output(['allosmod', 'get_pm_glyc', 'pm.pdb', 'list',
                                '42', '1', 'true', 'test_pm_glyc.pdb'],
                               cwd=tempdir)
            # assert on generated files
            for f in ('align2.ali', 'model_ini.py', 'model_ini0.py',
                      'model_glyc.py', 'get_rest.in', 'allosmod.py'):
                os.unlink(os.path.join(tempdir, f))

    def test_bad_bond_type(self):
        """Test handling of invalid O1 bond type"""
        self.assertRaises(allosmod.get_pm_glyc.BondTypeError,
                          allosmod.get_pm_glyc.Sugar,
                          "NAG", '13aa', '54')

    def test_script_only(self):
        """Simple get_pm_glyc script, script generation only"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, 'list'), 'w') as fh:
                fh.write('test_pm_glyc.pdb')
            out = check_output(['allosmod', 'get_pm_glyc', 'pm.pdb', 'list',
                                '42', '1', 'true', 'script'],
                               cwd=tempdir)
            # assert on generated files
            for f in ('model_ini.py', 'model_ini0.py',
                      'model_glyc.py'):
                os.unlink(os.path.join(tempdir, f))

    def test_remove_overlaps(self):
        """Test remove_overlaps"""
        gaps = [(1,4), (10, 20), (15, 17)]
        self.assertEqual(allosmod.get_pm_glyc.remove_overlaps(gaps),
                         [(1,4), (10,20)])
        gaps = [(1,4), (3, 7)]
        self.assertEqual(allosmod.get_pm_glyc.remove_overlaps(gaps),
                         [(1,7)])

    def test_get_first_unused_chain(self):
        """Test get_first_unused_chain"""
        self.assertEqual(allosmod.get_pm_glyc.get_first_unused_chain({}), 'A')
        self.assertEqual(allosmod.get_pm_glyc.get_first_unused_chain(
                                                      {'1':'A', '2':'C'}), 'B')
        d = {}
        for n, chain in enumerate(string.ascii_uppercase):
            d[str(n+1)] = chain
        self.assertRaises(ValueError,
                          allosmod.get_pm_glyc.get_first_unused_chain, d)

    def test_get_residue_chains(self):
        """Test get_residue_chains"""
        fname = os.path.join(test_dir, 'input', 'asite_pdb1.pdb')
        self.assertEqual(allosmod.get_pm_glyc.get_residue_chains(fname),
                         dict.fromkeys([str(x) for x in range(1, 31)], 'A'))

    def test_count_residues(self):
        """Test count_residues()"""
        fname = os.path.join(test_dir, 'input', 'test_pm_glyc.ali')
        self.assertEqual(allosmod.get_pm_glyc.count_residues(fname, 'pm.pdb'),
                         30)

    def test_read_glyc_file(self):
        """Test read_glyc_file()"""
        fname = os.path.join(test_dir, 'input', 'glyc.dat')
        c = allosmod.get_pm_glyc.read_glyc_file(fname)
        self.assertEqual(len(c), 2)
        self.assertEqual(len(c[0]), 8)
        self.assertEqual(len(c[1]), 8)
        s = c[0][0]
        self.assertEqual(s.monomer, 'NAG')
        self.assertEqual(s.bond_type, 'NGLB')
        self.assertEqual(s.attach_res, 20)
        self.assertEqual(s.one_letter_code, '1')
        self.assertEqual(s.get_connect_atom(), 'ND2')

    def test_check_attachments(self):
        """Test _check_attachments()"""
        fname = os.path.join(test_dir, 'input', 'glyc.dat')
        c = allosmod.get_pm_glyc.read_glyc_file(fname)
        chain_for_res = {'1': 'A'}
        self.assertRaises(allosmod.get_pm_glyc.InvalidResidueError,
                          allosmod.get_pm_glyc._check_attachments,
                          c, chain_for_res)

if __name__ == '__main__':
    unittest.main()
