import unittest
import modeller
import os
import sys
import subprocess
import collections
from allosmod.util import check_output

test_dir = os.path.dirname(sys.argv[0])

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_allosteric_site"""
        for args in ([], [''] * 5):
            out = check_output(['allosmod', 'get_allosteric_site'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_allosteric_site'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_get_fit_filename(self):
        """Test get_fit_filename()"""
        from allosmod.get_allosteric_site import get_fit_filename
        self.assertEqual(get_fit_filename('foo.pdb'), 'foo_fit.pdb')
        self.assertEqual(get_fit_filename('foo.ent'), 'foo.ent_fit.pdb')

    def test_simple(self):
        """Simple complete run of get_allosteric_site"""
        out = check_output(['allosmod', 'get_allosteric_site',
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
        out = check_output(['allosmod', 'get_allosteric_site',
                            os.path.join(test_dir, 'input',
                                         'asite_pdb1.pdb'),
                            os.path.join(test_dir, 'input',
                                         'asite_ligand.pdb'),
                            os.path.join(test_dir, 'input',
                                         'asite_pdb2.pdb'), '8.0'])

if __name__ == '__main__':
    unittest.main()
