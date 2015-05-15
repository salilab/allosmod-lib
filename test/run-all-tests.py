import unittest, sys, os, re
from optparse import OptionParser
import glob

# Only use coverage if it's new enough
try:
    import coverage
    if not hasattr(coverage.coverage, 'combine'):
        coverage = None
except ImportError:
    coverage = None

try:
    # Custom TestResult class that does not duplicate test name and docstring
    # when using unittest2
    from unittest import TextTestResult

    class _TestRunner(unittest.TextTestRunner):
        class _TestResult(TextTestResult):
            def getDescription(self, test):
                doc_first_line = test.shortDescription()
                if self.descriptions and doc_first_line:
                    return doc_first_line
                else:
                    return str(test)
        def _makeResult(self):
            return self._TestResult(self.stream, self.descriptions,
                                    self.verbosity)
except ImportError:
    # Older Pythons with unittest1 don't export TextTestResult and also
    # don't duplicate test names
    _TestRunner = unittest.TextTestRunner


class RunAllTests(unittest.TestProgram):
    """Custom main program that also displays a final coverage report"""
    def __init__(self, opts, *args, **keys):
        self.opts = opts
        if coverage:
            # Start coverage testing now before we import any modules
            cwd = os.path.dirname(sys.argv[0])
            self.topdir = os.path.abspath(os.path.join(cwd, '..', 'lib'))
            self.mods = glob.glob("%s/allosmod/*.py" % self.topdir)

            self.cov = coverage.coverage(branch=True, include=self.mods)
            self.cov.start()

        # Run the tests
        unittest.TestProgram.__init__(self, *args, **keys)

    def runTests(self):
        self.testRunner = _TestRunner(verbosity=self.verbosity)
        result = self.testRunner.run(self.test)

        if coverage:
            self.cov.stop()
            self.cov.combine()
            self.cov.use_cache(False)
            print >> sys.stderr, "\nPython coverage report\n"

            self.cov.file_locator.relative_dir = self.topdir + '/'
            self.cov.report(self.mods, file=sys.stderr)
            html = self.opts.html_coverage
            if html:
                self.cov.html_report(self.mods,
                                     directory=os.path.join(html, 'python'))
            for cov in glob.glob('.coverage.*'):
                os.unlink(cov)
        sys.exit(not result.wasSuccessful())


def regressionTest():
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    files = os.listdir(path)
    test = re.compile("^test_.*\.py$", re.IGNORECASE)
    files = filter(test.search, files)
    modnames = [os.path.splitext(f)[0] for f in files]

    modobjs = [__import__(m) for m in modnames]
    tests = [unittest.defaultTestLoader.loadTestsFromModule(o) for o in modobjs]
    return unittest.TestSuite(tests)

def parse_options():
    parser = OptionParser()
    parser.add_option("-v", dest="verbose", action='store_true',
                      help="verbose test output")
    parser.add_option("--html_coverage", dest="html_coverage", type="string",
                      default=None,
                      help="directory to write HTML coverage info into")
    return parser.parse_args()

if __name__ == "__main__":
    opts, args = parse_options()
    sys.argv = [sys.argv[0]] + args
    if opts.verbose:
        sys.argv.append('-v')
    RunAllTests(opts, defaultTest="regressionTest")
