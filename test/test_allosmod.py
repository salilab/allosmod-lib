import unittest
import os
import utils
from utils import check_output
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):
    def test_no_args(self):
        """Check 'allosmod' with no arguments"""
        out = check_output(['allosmod'], universal_newlines=True)
        self.assertTrue("Use 'allosmod help' for help" in out)

    def test_help(self):
        """Check 'allosmod help'"""
        for args in (['help'], ['help', 'help']):
            out = check_output(['allosmod'] + args, universal_newlines=True)
            self.assertTrue('Get help on using' in out)

    def test_command_help(self):
        """Check 'allosmod help' for a subcommand"""
        out = check_output(['allosmod', 'help', 'pdb2ali'],
                           universal_newlines=True)
        self.assertTrue('Convert a PDB file' in out)

    def test_unknown_command(self):
        """Check 'allosmod' with an unknown command"""
        for args in (['bad-command'], ['help', 'bad-command']):
            check_output(['allosmod'] + args, retcode=1)


if __name__ == '__main__':
    unittest.main()
