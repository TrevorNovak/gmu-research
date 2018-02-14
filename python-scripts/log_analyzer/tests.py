import unittest
from processor import *

class TestAnalyzerMethods(unittest.TestCase):

    def test_build_matrix(self):
        p = Processor()
        p.build_matrix(3, 3)
        self.assertEqual(p.encoded_matrix, [[0,0,0],[0,0,0],[0,0,0]])

if __name__ == '__main__':
    unittest.main()
