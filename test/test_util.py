import unittest
import os

class Tests(unittest.TestCase):
    def test_read_templates(self):
        """Test read_templates function"""
        import allosmod.util
        with open('lst', 'w') as fh:
            fh.write('test1 A 1 10\ntest2\ntest3 B 4 5\n')
        self.assertEqual(allosmod.util.read_templates('lst'),
                         ['test1', 'test2', 'test3'])
        os.unlink('lst')

if __name__ == '__main__':
    unittest.main()
