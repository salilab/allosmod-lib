import unittest
import os
import subprocess
from allosmod.util import check_output
import allosmod.setup

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to setup"""
        for args in ([''],):
            out = check_output(['allosmod', 'setup'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.setup'] + args,
                               stderr=subprocess.STDOUT, retcode=2)

    def test_make_config_file(self):
        """Test ConfigFile constructor"""
        c = allosmod.setup.ConfigFile(['foo', 'bar'], glyco=True, maxres=20)
        self.assertEqual(c.glyco, True)
        self.assertEqual(c.maxres, 20)
        self.assertEqual(c.num_templates, 2)
        c = allosmod.setup.ConfigFile(None)
        self.assertEqual(c.glyco, False)
        self.assertEqual(c.maxres, 0)
        self.assertEqual(c.num_templates, 0)

    def test_access_config_file(self):
        """Test ConfigFile access"""
        c = allosmod.setup.ConfigFile(None)
        c['deleMaX'] = 'foo'
        self.assertEqual(c['DELEMAX'], 'foo')
        self.assertTrue('DeLeMaX' in c)
        self.assertFalse('garbage' in c)
        self.assertEqual(c.get('DeLEMAX', 'bar'), 'foo')
        self.assertEqual(c.get('garbage', 'bar'), 'bar')

    def test_write_config_file(self):
        """Test ConfigFile.write()"""
        c = allosmod.setup.ConfigFile(None)
        c.write('test.dat')
        os.unlink('test.dat')

    def parse_config_file(self, contents, *args, **keys):
        with open('test.dat', 'w') as fh:
            fh.write(contents)
        c = allosmod.setup.ConfigFile(*args, **keys)
        try:
            return c, list(c.parse('test.dat'))
        finally:
            os.unlink('test.dat')

    def test_parse_config_file_invalid(self):
        """Test ConfigFile.parse() with invalid format"""
        c, errs = self.parse_config_file("DEVIATION=garbage", None)
        self.assertEqual(errs, ['Missing variable in test.dat: NRUNS',
                                'Invalid variable in test.dat: DEVIATION: '
                                'invalid literal for float(): garbage'])

    def test_parse_config_file_bad_delemax(self):
        """Test ConfigFile.parse() with bad delEmax"""
        c, errs = self.parse_config_file("NRUNS=1\ndelEmax=garbage", None)
        self.assertEqual(errs, ['Invalid variable in test.dat: DELEMAX: '
                                'invalid literal for float(): garbage'])

    def test_parse_config_file_bad_sampling(self):
        """Test ConfigFile.parse() with bad sampling"""
        c, errs = self.parse_config_file("NRUNS=1\nsampling=garbage", None)
        self.assertEqual(errs, ['Invalid variable in test.dat: SAMPLING: '
                           'not one of simulation, moderate_cm, moderate_am, '
                           'fast_cm'])

    def test_parse_config_file_bad_mdtemp(self):
        """Test ConfigFile.parse() with bad mdtemp"""
        c, errs = self.parse_config_file("NRUNS=1\nMDTEMP=garbage", None)
        self.assertEqual(errs, ['Invalid variable in test.dat: MDTEMP: '
                                'invalid literal for float(): garbage'])

    def test_parse_config_file_bad_boolean(self):
        """Test ConfigFile.parse() with bad boolean"""
        c, errs = self.parse_config_file("NRUNS=1\ncoarse=garbage", None)
        self.assertEqual(errs, ['Invalid variable in test.dat: COARSE: '
                            'not one of 1, yes, true, on, 0, no, false, off'])

    def test_parse_config_file_ok(self):
        """Test ConfigFile.parse() with ok file"""
        c, errs = self.parse_config_file("""
NRUNS=1
MDTEMP=SCAN
delEmax=100.0
rAS=40
SAMPLING=moderate_cm
COARSE=false
LOCALRIGID=yes
""", None)
        self.assertEqual(len(errs), 0)
        self.assertEqual(c['NRUNS'], 1)
        self.assertEqual(c['MDTEMP'], 'scan')
        self.assertEqual(c['RAS'], 40)
        self.assertAlmostEqual(c['delEmax'], 100.0, places=1)
        self.assertEqual(c['SAMPLING'], 'moderate_cm')
        self.assertEqual(c['COARSE'], False)
        self.assertEqual(c['LOCALRIGID'], True)

    def test_parse_config_file_defaults(self):
        """Test ConfigFile.parse() with defaults"""
        c, errs = self.parse_config_file("""
NRUNS=1
ASPDB=
""", ["temp1", "temp2"])
        self.assertEqual(len(errs), 0)
        self.assertEqual(c['NRUNS'], 1)
        self.assertAlmostEqual(c['MDTEMP'], 300.0, places=1)
        self.assertEqual(c['ASPDB'], 'temp1')

    def test_parse_config_file_glyco(self):
        """Test ConfigFile.parse() with glycocsylation"""
        c, errs = self.parse_config_file("NRUNS=1\nDEVIATION=4\n"
                                         "SAMPLING=moderate_am",
                                         None, glyco=True)
        self.assertEqual(len(errs), 0)
        self.assertEqual(c['deviation'], 0.)
        self.assertEqual(c['SAMPLING'], 'moderate_cm')

    def test_parse_config_file_one_struc(self):
        """Test ConfigFile.parse() with one structure"""
        c, errs = self.parse_config_file("NRUNS=1\nrAS=40.\n"
                                         "SAMPLING=moderate_am", ["template"])
        self.assertEqual(len(errs), 0)
        self.assertAlmostEqual(c['rAS'], 1000., places=1)

    def test_parse_config_file_big_system(self):
        """Test ConfigFile.parse() with big system"""
        c, errs = self.parse_config_file("NRUNS=1\nCOARSE=false", None,
                                         maxres=2000)
        self.assertEqual(len(errs), 0)
        self.assertEqual(c['COARSE'], True)

    def test_with_glyc2(self):
        """Test Setup.with_glyc2() method"""
        s = allosmod.setup.Setup()
        with open('allosmod.py', 'w') as fh:
            fh.write('line1\nline2\n')
        self.assertEqual(s.with_glyc2(), None)
        with open('allosmod.py', 'w') as fh:
            fh.write("self.patch(residue_type='foo')\n")
        self.assertEqual(s.with_glyc2(), True)
        os.unlink('allosmod.py')
        self.assertEqual(s.with_glyc2(), None)

    def test_get_other_pdb(self):
        """Test Setup.get_other_pdb() method"""
        s = allosmod.setup.Setup()
        s.config = {'ASPDB': 'test_aspdb'}
        s.templates = ['foo', 'bar', 'test_aspdb']
        self.assertEqual(s.get_other_pdb(), 'foo')
        s.templates = ['test_aspdb', 'bar']
        self.assertEqual(s.get_other_pdb(), 'bar')
        s.templates = []
        self.assertEqual(s.get_other_pdb(), 'test_aspdb')

    def test_simple_fail(self):
        """Simple complete failing run of setup"""
        with allosmod.util.temporary_directory() as tempdir:
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'Missing file: input.dat\n'
                         'Missing file: align.ali\nMissing file: list\n')

    def test_simple_ok(self):
        """Simple complete ok run of setup"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1\nHARM 10.0 0.5 1,2\n'
                          'LOBD 4.0 0.1 2,1\nUPBD 4.0 0.1 2,1')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "foo"), 'w') as fh:
                fh.write('ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C\n')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=0)
            os.unlink(os.path.join(tempdir, "lig.pdb"))
            os.unlink(os.path.join(tempdir, "qsub.sh"))
        self.assertEqual(out, '')

    def test_simple_bad_list(self):
        """Simple complete run of setup with bad list file"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('bar')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'Missing file: bar\nMissing sequence in '
                              'align.ali for file: bar\n'
                              'Missing LIGPDB in align.ali for file: bar\n'
                              'Missing ASPDB in align.ali for file: bar\n')

    def test_simple_bad_align(self):
        """Simple complete run of setup with bad alignment file"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFVAFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'Sequences in align.ali are not properly '
                              'aligned\nMissing file: foo\n')

    def test_simple_bad_aspdb(self):
        """Simple complete run of setup with bad ASPDB sequence"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1\nASPDB=bar')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'Missing file: foo\n'
                              'Missing ASPDB in align.ali for file: bar\n')

    def test_simple_bad_input_dat(self):
        """Simple complete run of setup with bad input.dat"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1\nMDTEMP=garbage\n')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "foo"), 'w') as fh:
                fh.write('ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C\n')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'Invalid variable in input.dat: MDTEMP: '
                              'invalid literal for float(): garbage\n')

    def test_simple_alignment_error(self):
        """Simple complete run of setup with alignment format error"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1\n')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "foo"), 'w') as fh:
                fh.write('ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C\n')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm.pdb
sequence:::::::::
AFV""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out, 'PIR sequence without terminating * at '
                              'end of file\n')

    def test_simple_no_target(self):
        """Simple complete run of setup with target not in alignment"""
        with allosmod.util.temporary_directory() as tempdir:
            with open(os.path.join(tempdir, "input.dat"), 'w') as fh:
                fh.write('NRUNS=1\n')
            with open(os.path.join(tempdir, "list"), 'w') as fh:
                fh.write('foo')
            with open(os.path.join(tempdir, "foo"), 'w') as fh:
                fh.write('ATOM     14  CA  TYR A   2      -6.696  -0.319   4.300  1.00  0.00           C\n')
            with open(os.path.join(tempdir, "align.ali"), 'w') as fh:
                fh.write(""">P1;foo
structureX:::::::::
AFV*
>P1;pm2.pdb
sequence:::::::::
AFV*""")
            out = check_output(['allosmod', 'setup'],
                               stderr=subprocess.STDOUT, cwd=tempdir, retcode=1)
        self.assertEqual(out,
                         'Missing sequence in align.ali for file: pm.pdb\n')

if __name__ == '__main__':
    unittest.main()
