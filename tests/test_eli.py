#!/usr/bin/env python3

import unittest
from rimopy import eli

KERNELS_PATH = "./kernels"

class TestSum(unittest.TestCase):
    
    def test_getELI_Valladolid(self):
        res = eli.getELI(400, 41.652251, -4.7245321, "2022-01-17 02:30:00", KERNELS_PATH)
        self.assertGreater(res, 0, "Should be greater than 0")

if __name__ == '__main__':
    unittest.main()