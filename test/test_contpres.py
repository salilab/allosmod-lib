import unittest
import os
import sys
import subprocess
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
test_dir = utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to contpres"""
        for args in ([], [''] * 4):
            check_output(['allosmod', 'contpres'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.contpres'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of contpres"""
        with utils.temporary_directory() as tmpdir:
            break_dat = os.path.join(tmpdir, 'break.dat')
            with open(break_dat, 'w') as fh:
                fh.write("foo\nbar\n")
            check_output(['allosmod', 'contpres',
                          os.path.join(test_dir, 'input',
                                       'test_contpres.rsr'),
                          os.path.join(test_dir, 'input',
                                       'test_contpres.pdb'), '50.0'],
                         cwd=tmpdir)
            with open(break_dat) as fh:
                content = fh.read()
            self.assertEqual(content, "foo\nbar\n1 50.0\n3 50.0\n")
            with open(os.path.join(tmpdir, 'contpres.dat')) as fh:
                content = fh.read()
            self.assertEqual(content, "1 4\n2 0\n3 4\n")

    def test_cdensity(self):
        """Simple contpres in cdensity mode"""
        with utils.temporary_directory() as tmpdir:
            break_dat = os.path.join(tmpdir, 'break.dat')
            with open(break_dat, 'w') as fh:
                fh.write("foo\nbar\n")
            check_output(['allosmod', 'contpres', '--cdensity_cutoff',
                          '1.0', os.path.join(test_dir, 'input',
                                              'test_contpres.rsr'),
                          os.path.join(test_dir, 'input',
                                       'test_contpres.pdb'), '50.0'],
                         cwd=tmpdir)
            with open(break_dat) as fh:
                content = fh.read()
            self.assertEqual(content, "foo\nbar\n2 50.0\n")
            os.unlink(break_dat)
            os.unlink(os.path.join(tmpdir, 'contpres.dat'))

    def test_cdensity_even_cutoff_negative(self):
        """Test contpres with even charge density, negative cutoff"""
        with utils.temporary_directory() as tmpdir:
            break_dat = os.path.join(tmpdir, 'break.dat')
            with open(break_dat, 'w') as fh:
                fh.write("foo\nbar\n")
            check_output(['allosmod', 'contpres', '--cdensity_cutoff',
                          '-1.0', os.path.join(test_dir, 'input',
                                               'contpres_evencd.rsr'),
                          os.path.join(test_dir, 'input',
                                       'contpres_evencd.pdb'), '50.0'],
                         cwd=tmpdir)
            with open(break_dat) as fh:
                content = fh.read()
            self.assertEqual(content, "foo\nbar\n1 50.0\n2 50.0\n")
            os.unlink(break_dat)
            os.unlink(os.path.join(tmpdir, 'contpres.dat'))

    def test_cdensity_even_cutoff_positive(self):
        """Test contpres with even charge density, positive cutoff"""
        with utils.temporary_directory() as tmpdir:
            break_dat = os.path.join(tmpdir, 'break.dat')
            with open(break_dat, 'w') as fh:
                fh.write("foo\nbar\n")
            check_output(['allosmod', 'contpres', '--cdensity_cutoff',
                          '1.0', os.path.join(test_dir, 'input',
                                              'contpres_evencd.rsr'),
                          os.path.join(test_dir, 'input',
                                       'contpres_evencd.pdb'), '50.0'],
                         cwd=tmpdir)
            with open(break_dat) as fh:
                content = fh.read()
            self.assertEqual(content, "foo\nbar\n")
            os.unlink(break_dat)
            os.unlink(os.path.join(tmpdir, 'contpres.dat'))


if __name__ == '__main__':
    unittest.main()
