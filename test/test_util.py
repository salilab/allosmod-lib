import unittest
import os

class Tests(unittest.TestCase):
    def test_read_templates(self):
        """Test read_templates function"""
        import allosmod.util
        with open('lst', 'w') as fh:
            fh.write('test1\ntest2\n')
        self.assertEqual(allosmod.util.read_templates('lst'),
                         ['test1', 'test2'])
        os.unlink('lst')

if __name__ == '__main__':
    unittest.main()
