import unittest
from .processor import *

class TestProcessorMethods(unittest.TestCase):

    def test_build_matrix01(self):
        p = Processor()
        p.build_matrix(3, 3)
        self.assertEqual(p.encoded_matrix, [[0,0,0],[0,0,0],[0,0,0]])

    def test_build_matrix02(self):
        p = Processor()
        p.build_matrix(3, 4)
        self.assertNotEqual(p.encoded_matrix, [[0,0,0],[0,0,0],[0,0,0]])

if __name__ == '__main__':
    unittest.main()
