#!/usr/bin/env python3

import unittest
from rimopy import eli

KERNELS_PATH = "./kernels"

VALL_LAT = 41.652251
VALL_LON = -4.7245321
VALL_ALT = 700
JAN_FULL_MOON_00 = "2022-01-17 00:00:00"
PROP_ERROR = 0.05 # up to 5% of error is allowed from AEMET's RimoApp

class TestSum(unittest.TestCase):
    
    def test_getELI_Valladolid(self):
        res = eli.getELI(400, VALL_LAT, VALL_LON, "2022-01-17 02:30:00", KERNELS_PATH)
        self.assertGreater(res, 0, "Should be greater than 0")

    def test_eli336_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(336, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 9.1239e-07
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

    def test_eli380_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(380, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 1.3348e-06
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

    def test_eli440_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(440, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 2.4528e-06
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

    def test_eli500_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(500, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 2.9900e-06
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

    def test_eli862_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(862, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 2.2911e-06
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

    def test_eli1662_uncorrected_Valladolid_20220117(self):
        res = eli.getELIPerNm(1662, VALL_LAT, VALL_LON, JAN_FULL_MOON_00, KERNELS_PATH, VALL_ALT, False)
        expected = 8.4489e-07
        self.assertAlmostEqual(res, expected, delta=expected*PROP_ERROR)

if __name__ == '__main__':
    unittest.main()