import unittest
import modeller
import os
import subprocess
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

from allosmod.util import check_output

def get_pdb_line(attyp, resnum, typ='ATOM'):
    return "%-6s    1  %-3s TYR A%4d      24.417  18.891   8.203  1.00  0.00           C\n" % (typ, attyp, resnum)

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to get_add_restraint"""
        for args in ([], ['x'] * 4):
            out = check_output(['allosmod', 'get_add_restraint'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.get_add_restraint'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_get_restraints(self):
        """Test get_restraints() function"""
        from allosmod.get_add_restraint import get_restraints
        with open('test.dat', 'w') as fh:
            fh.write("""HARM 10.0 0.5 1,2,3,4
HARM 4.0 0.1 2,1
UPBD 4.0 0.1 2,1,   
OTHEROPTION foo
""")
        r = list(get_restraints('test.dat', 'HARM'))
        self.assertEqual(len(r), 3)
        self.assertAlmostEqual(r[0].distance, 10.0, places=1)
        self.assertAlmostEqual(r[0].stddev, 0.5, places=1)
        self.assertEqual(r[0].resind1, '1')
        self.assertEqual(r[0].resind2, '2')
        self.assertAlmostEqual(r[1].distance, 10.0, places=1)
        self.assertAlmostEqual(r[1].stddev, 0.5, places=1)
        self.assertEqual(r[1].resind1, '3')
        self.assertEqual(r[1].resind2, '4')
        self.assertAlmostEqual(r[2].distance, 4.0, places=1)
        self.assertAlmostEqual(r[2].stddev, 0.1, places=1)
        self.assertEqual(r[2].resind1, '2')
        self.assertEqual(r[2].resind2, '1')

        r = list(get_restraints('test.dat', 'UPBD'))
        self.assertEqual(len(r), 1)
        self.assertAlmostEqual(r[0].distance, 4.0, places=1)
        self.assertAlmostEqual(r[0].stddev, 0.1, places=1)
        self.assertEqual(r[0].resind1, '2')
        self.assertEqual(r[0].resind2, '1')
        os.unlink('test.dat')

    def test_get_atom_indexes(self):
        """Test get_atom_indexes() function"""
        from allosmod.get_add_restraint import get_atom_indexes
        with open('test.dat', 'w') as fh:
            fh.write(get_pdb_line('DUM', 3))
            fh.write(get_pdb_line('DUM', 6, typ='HETATM'))
            fh.write('END\n')
        i = get_atom_indexes('test.dat')
        self.assertEqual(i, {'3': 1 ,'6': 2})
        with open('test.dat', 'w') as fh:
            fh.write(get_pdb_line('P', 4))
            fh.write(get_pdb_line('CA', 4))
        i = get_atom_indexes('test.dat')
        self.assertEqual(i, {'4': 2})
        with open('test.dat', 'w') as fh:
            pass
        i = get_atom_indexes('test.dat')
        self.assertEqual(i, {})
        os.unlink('test.dat')

    def test_simple(self):
        """Simple complete test of get_add_restraint command"""
        with utils.temporary_directory() as tmpdir:
            with open(os.path.join(tmpdir, 'test.pdb'), 'w') as fh:
                fh.write(get_pdb_line('P', 4))
                fh.write(get_pdb_line('CA', 4))
                fh.write(get_pdb_line('C', 4))
                fh.write(get_pdb_line('O', 6))
                fh.write(get_pdb_line('C', 6))
            with open(os.path.join(tmpdir, 'input.dat'), 'w') as fh:
                fh.write("""HARM 10.0 0.5 4,6
HARM 4.0 0.1 2,1
HARM 4.0 0.1 4,8
""")
            out = check_output(['allosmod', 'get_add_restraint', 'input.dat',
                                'test.pdb', 'HARM'],
                               universal_newlines=True, cwd=tmpdir)
            self.assertEqual(out, 'R    3   1   1  27   2   2   1     '
                                  '2     5      10.0000    0.5000\n')

if __name__ == '__main__':
    unittest.main()
