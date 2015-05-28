from __future__ import print_function
import os
import unittest

class Tests(unittest.TestCase):
    def test_truncated_gaussian(self):
        """Test TruncatedGaussian math form"""
        import modeller
        from allosmod.modeller.forms import TruncatedGaussian
        env = modeller.environ()
        env.edat.dynamic_sphere= False
        with open('test.pdb', 'w') as fh:
            fh.write("ATOM      2  CA  ALA     1      27.449  14.935   5.140  1.00 29.87           C\n")
        m = modeller.model(env, file='test.pdb')
        os.unlink('test.pdb')
        feat = modeller.features.x_coordinate(m.atoms[0])
        f = TruncatedGaussian(group=modeller.physical.xy_distance, feature=feat,
                              dele_max=10, slope=4.0, scl_delx=0.7,
                              weights=(1,1), means=(14, 28), stdevs=(1, 2))
        m.restraints.add(f)
        sel = modeller.selection(m)
        e = []
        for i in range(30):
            m.atoms[0].x = 2.0 * i
            e.append(sel.objfunc())
        expected_e = [10.133, 10.133, 10.133, 10.133, 10.172, 5.263, 1.722, 0.542, 1.722, 5.260, 5.672, 3.607, 2.131, 1.246, 0.951, 1.246, 2.131, 3.607, 5.672, 8.314, 10.135, 10.133, 10.133, 10.133, 10.133, 10.133, 10.133, 10.133, 10.133, 10.133]
        for a, b in zip(e, expected_e):
            self.assertAlmostEqual(a, b, places=1)

if __name__ == '__main__':
    unittest.main()
