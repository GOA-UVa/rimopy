#!/usr/bin/env python3

import unittest
import sys
from rimopy import eli

KERNELS_PATH = "./kernels"

VALL_LAT = 41.652251
VALL_LON = -4.7245321
VALL_ALT = 700
JAN_FULL_MOON_00 = "2022-01-17 00:00:00"
JAN_FULL_MOON_17 = "2022-01-17 17:00:00"
FEB_NEW_MOON_00 = "2022-02-02 00:00:00"
DEFAULT_PROP_ERROR = 0.10 # up to 5% of error is allowed from AEMET's RimoApp by default

prop_error = DEFAULT_PROP_ERROR

def testValladolidNoCorr(ts: 'TestSum', wavelength, expected, date):
    res = eli.getELIPerNm(wavelength, VALL_LAT, VALL_LON, date, KERNELS_PATH, VALL_ALT, False)
    ts.assertAlmostEqual(res, expected, delta=expected*prop_error)

class TestSum(unittest.TestCase):
    
    def test_getELI_Valladolid(self):
        res = eli.getELI(400, VALL_LAT, VALL_LON, "2022-01-17 02:30:00", KERNELS_PATH)
        self.assertGreater(res, 0, "Should be greater than 0")

    def test_eli336_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 336, 9.1239e-07, JAN_FULL_MOON_00)

    def test_eli380_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 380, 1.3348e-06, JAN_FULL_MOON_00)

    def test_eli440_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 440, 2.4528e-06, JAN_FULL_MOON_00)

    def test_eli500_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 500, 2.9900e-06, JAN_FULL_MOON_00)

    def test_eli862_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 862, 2.2911e-06, JAN_FULL_MOON_00)

    def test_eli1011_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 1011, 1.8330e-06, JAN_FULL_MOON_00)

    def test_eli1662_uncorrected_Valladolid_20220117_00(self):
        testValladolidNoCorr(self, 1662, 8.4489e-07, JAN_FULL_MOON_00)

    def test_eli338_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 338, 1.2293e-06, JAN_FULL_MOON_17)

    def test_eli385_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 385, 1.5392e-06, JAN_FULL_MOON_17)

    def test_eli481_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 481, 4.0048e-06, JAN_FULL_MOON_17)

    def test_eli540_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 540, 3.8516e-06, JAN_FULL_MOON_17)

    def test_eli879_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 879, 2.7428e-06, JAN_FULL_MOON_17)

    def test_eli1020_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 1020, 2.2381e-06, JAN_FULL_MOON_17)

    def test_eli1654_uncorrected_Valladolid_20220117_17(self):
        testValladolidNoCorr(self, 1654, 1.0285e-06, JAN_FULL_MOON_17)

    def test_eli336_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 336, 5.4519e-10, FEB_NEW_MOON_00)

    def test_eli380_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 380, 1.0571e-09, FEB_NEW_MOON_00)

    def test_eli440_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 440, 1.8795e-09, FEB_NEW_MOON_00)

    def test_eli500_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 500, 2.7575e-09, FEB_NEW_MOON_00)

    def test_eli862_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 862, 2.1584e-09, FEB_NEW_MOON_00)

    def test_eli1011_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 1011, 7.5293e-10, FEB_NEW_MOON_00)

    def test_eli1662_uncorrected_Valladolid_20220202_00(self):
        testValladolidNoCorr(self, 1662, 7.9724e-10, FEB_NEW_MOON_00)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        prop_error = float(sys.argv[1])
    else:
        prop_error = DEFAULT_PROP_ERROR
    unittest.main(argv=['first-arg-is-ignored'], exit=False)