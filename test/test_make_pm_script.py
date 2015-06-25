import os
import unittest
import subprocess
from test_pdb2ali import check_output

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments to make_pm_script"""
        for args in ([], ['' * 7]):
            out = check_output(['allosmod', 'make_pm_script'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
            out = check_output(['python', '-m',
                                'allosmod.make_pm_script'] + args,
                               stderr=subprocess.STDOUT, retcode=2)
        out = check_output(['allosmod', 'make_pm_script', '--', 'target',
                            'list', '-333', '1.0', 'garbage', '300.0'],
                            stderr=subprocess.STDOUT, retcode=2)

    def test_simple(self):
        """Simple complete run of make_pm_script"""
        with open('list', 'w') as fh:
            fh.write("foo\nbar\n")
        for t in ['script', 'moderate_cm', 'moderate_cm_simulation',
                  'fast_cm', 'fast_cm_simulation', 'moderate_am']:
            check_output(['allosmod', 'make_pm_script', '--', 'target', 'list',
                          '-333', '4.0', t, '300.0'])
        for f in ('model_run.py', 'list'):
            os.unlink(f)

if __name__ == '__main__':
    unittest.main()
