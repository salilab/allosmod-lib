import unittest
import modeller.automodel.randomize
import contextlib
from io import BytesIO
import os
import utils
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import allosmod.modeller

@contextlib.contextmanager
def mock_method(cls, method_name, replacement=None):
    def mock(*args, **keys):
        mock.args = args
        mock.keys = keys
    if replacement is None:
        replacement = mock
    old_method = getattr(cls, method_name)
    setattr(cls, method_name, replacement)
    yield mock
    setattr(cls, method_name, old_method)

class WriteIntMock(allosmod.modeller.AllosModel):
    ints = None
    def __init__(self, *args, **keys):
        allosmod.modeller.AllosModel.__init__(self, *args, **keys)
        self.ints = []
    def write_int(self, ctr, suff):
        self.ints.append((ctr, suff))

test_pdb = """
ATOM      7  CA  HIS     7      17.121  17.162   6.197  1.00 15.60           C
"""

class Tests(unittest.TestCase):
    env = modeller.environ()

    def test_refine_none(self):
        """Test 'none' refinement"""
        m = allosmod.modeller.AllosModel(self.env, 'dummyaln', 'known', 'seq')
        m.read(file=BytesIO(test_pdb))
        atmsel = modeller.selection(m)
        with mock_method(modeller.optimizers.molecular_dynamics, 'optimize'):
            allosmod.modeller.none(atmsel, [])

    def test_refine_moderate(self):
        """Test 'moderate' refinement"""
        m = allosmod.modeller.AllosModel(self.env, 'dummyaln', 'known', 'seq')
        m.read(file=BytesIO(test_pdb))
        atmsel = modeller.selection(m)
        with mock_method(modeller.optimizers.molecular_dynamics, 'optimize'):
            allosmod.modeller.moderate(atmsel, [])

    def test_refine_moderate_am(self):
        """Test 'ModerateAM' refinement"""
        m = WriteIntMock(self.env, 'dummyaln', 'known', 'seq')
        m.read(file=BytesIO(test_pdb))
        atmsel = modeller.selection(m)
        with mock_method(modeller.optimizers.molecular_dynamics, 'optimize'):
            c = allosmod.modeller.ModerateAM(md_temp=300.0)
            c(atmsel, [])
        self.assertEqual(m.ints, [(x, 1) for x in range(501, 507)] +
                                 [(x, 1) for x in range(1001, 1101)])

    def test_refine_consttemp(self):
        """Test 'ConstTemp' refinement"""
        m = WriteIntMock(self.env, 'dummyaln', 'known', 'seq')
        m.read(file=BytesIO(test_pdb))
        atmsel = modeller.selection(m)
        with mock_method(modeller.optimizers.molecular_dynamics, 'optimize'):
            c = allosmod.modeller.ConstTemp(md_temp=300.0)
            c(atmsel, [])
        self.assertEqual(m.ints, [(x, 1) for x in range(501, 511)] +
                                 [(x, 1) for x in range(1001, 3001)])

    def test_allos_model(self):
        """Test AllodModel class"""
        m = allosmod.modeller.AllosModel(self.env, 'dummyaln', 'known', 'seq',
                                         deviation=0.0)
        self.assertEqual(m.rand_method, None)
        m = allosmod.modeller.AllosModel(self.env, 'dummyaln', 'known', 'seq',
                                         deviation=4.0)
        self.assertEqual(m.rand_method, modeller.automodel.randomize.xyz)

    def test_allos_model_refine(self):
        """Test AllodModel.refine method"""
        class MockAllosModel(allosmod.modeller.AllosModel):
            calls = None
            def __init__(self, *args, **keys):
                allosmod.modeller.AllosModel.__init__(self, *args, **keys)
                self.calls = []
            def initial_refine_hot(self, atmsel):
                self.calls.append('initial_refine_hot')
            def final_refine_hot(self, atmsel):
                self.calls.append('final_refine_hot')
            def fit_refined(self, fname):
                self.calls.append('fit_refined')
        def test_md_level(atmsel, actions):
            atmsel.get_model().calls.append('md_level')
        m = MockAllosModel(self.env, 'dummyaln', 'known', 'seq')
        m.fit_in_refine = 'NO_FIT'
        m.refine_hot_only = True
        m.md_level = test_md_level
        m.refine(modeller.selection(m), [])
        self.assertEqual(m.calls, ['initial_refine_hot', 'md_level',
                                   'final_refine_hot'])

        m = MockAllosModel(self.env, 'dummyaln', 'known', 'seq')
        m.fit_in_refine = 'NO_FIT'
        m.refine_hot_only = False
        m.md_level = None
        m.refine(modeller.selection(m), [])
        self.assertEqual(m.calls, [])

        m = MockAllosModel(self.env, 'dummyaln', 'known', 'seq')
        m.read(file=BytesIO(test_pdb))
        m.fit_in_refine = ''
        m.refine_hot_only = False
        m.md_level = None
        m.refine(modeller.selection(m), [])
        self.assertEqual(m.calls, ['fit_refined'])
        os.unlink('TO_BE_REFINED.TMP')

if __name__ == '__main__':
    unittest.main()
