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
        """Test wrong arguments to get_allosteric_site"""
        for args in ([], [''] * 5):
            check_output(['allosmod', 'get_allosteric_site'] + args,
                         stderr=subprocess.STDOUT, retcode=2)
            check_output([sys.executable, '-m',
                          'allosmod.get_allosteric_site'] + args,
                         stderr=subprocess.STDOUT, retcode=2)

    def test_get_fit_filename(self):
        """Test get_fit_filename()"""
        from allosmod.get_allosteric_site import get_fit_filename
        self.assertEqual(get_fit_filename('foo.pdb'), 'foo_fit.pdb')
        self.assertEqual(get_fit_filename('foo.ent'), 'foo.ent_fit.pdb')

    def test_simple(self):
        """Simple complete run of get_allosteric_site"""
        check_output(['allosmod', 'get_allosteric_site',
                      '--output_pdb', 'allostericsite.pdb',
                      '--atom_list', 'atomlistASRS',
                      os.path.join(test_dir, 'input',
                                   'asite_pdb1.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_ligand.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_pdb2.pdb'), '8.0'])
        with open('atomlistASRS') as fh:
            wc = len(fh.readlines())
        self.assertEqual(wc, 232)
        with open('allostericsite.pdb') as fh:
            wc = len(fh.readlines())
        self.assertEqual(wc, 119)
        os.unlink('atomlistASRS')
        os.unlink('allostericsite.pdb')

    def test_simple_no_outputs(self):
        """Simple complete run of get_allosteric_site, with no outputs"""
        check_output(['allosmod', 'get_allosteric_site',
                      os.path.join(test_dir, 'input',
                                   'asite_pdb1.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_ligand.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_pdb2.pdb'), '8.0'])

    def test_simple_no_site(self):
        """Simple complete run of get_allosteric_site, with no site"""
        check_output(['allosmod', 'get_allosteric_site',
                      '--output_pdb', 'allostericsite.pdb',
                      '--atom_list', 'atomlistASRS',
                      os.path.join(test_dir, 'input',
                                   'asite_pdb1.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_far_ligand.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_pdb2.pdb'), '8.0'], retcode=1)

    def test_simple_no_align(self):
        """Simple complete run of get_allosteric_site with no alignment"""
        check_output(['allosmod', 'get_allosteric_site',
                      '--output_pdb', 'allostericsite.pdb',
                      '--atom_list', 'atomlistASRS',
                      os.path.join(test_dir, 'input',
                                   'asite_pdb1.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_ligand.pdb'),
                      os.path.join(test_dir, 'input',
                                   'asite_ligand.pdb'), '8.0'], retcode=1)


if __name__ == '__main__':
    unittest.main()
