#!/usr/bin/env python3
# Script that tests against some AEMET RimoApp cases.

import unittest
import sys
import os

from rimopy import eli, esi


KERNELS_PATH = os.path.join(os.path.dirname(__file__), "./kernels")
VALL_LAT = 41.6636
VALL_LON = -4.70583
VALL_ALT = 705
JAN_FULL_MOON_00 = "2022-01-17 00:00:00"
JAN_FULL_MOON_17 = "2022-01-17 17:00:00"
FEB_NEW_MOON_00 = "2022-02-02 00:00:00"
DEFAULT_PROP_ERROR = 0  # 0% error, now the result is formatted as the one in AEMET's
ed_Vall = eli.EarthPoint(VALL_LAT, VALL_LON, JAN_FULL_MOON_00, VALL_ALT)
EXK = ["EarthStations.bsp", "EarthStations.tf"]
EXK_PATH = "./extra.temp.kernels"
VALL_NAME = "VALLADOLID"

prop_error = DEFAULT_PROP_ERROR
gfp = esi.GaussianFilterParams(1, 1)
calc = esi.ESICalculatorWehrli(
    esi.WehrliFile.ASC_WEHRLI, esi.ESIMethod.LINEAR_INTERPOLATION, gfp
)
eli_settings = eli.ELISettings(False, False, True, False)


def _testValladolidNoCorr(ts: "TestSum", wavelength, expected, date):
    ed_Vall.set_utc_times(date)
    eli_setts_2 = eli.ELISettings(False, False, True, True)
    res = eli.get_eli([wavelength], ed_Vall, KERNELS_PATH, calc, eli_setts_2)[0]
    # res = eli.get_eli_per_nm_from_extra_kernels(wavelength, ed_Vall.utc_times,
    #    KERNELS_PATH, EXK, EXK_PATH, VALL_NAME, calc, eli_settings)
    res = float("{:0.4e}".format(res))
    ts.assertAlmostEqual(res, expected, delta=expected * prop_error)


class TestSum(unittest.TestCase):
    def test_get_eli_Valladolid(self):
        ed_Vall_t = eli.EarthPoint(VALL_LAT, VALL_LON, "2022-01-17 02:30:00", VALL_ALT)
        res = eli.get_eli([400], ed_Vall_t, KERNELS_PATH)
        self.assertGreater(res, 0, "Should be greater than 0")

    def test_eli336_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 336, 9.1239e-07, JAN_FULL_MOON_00)

    def test_eli380_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 380, 1.3348e-06, JAN_FULL_MOON_00)

    def test_eli440_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 440, 2.4528e-06, JAN_FULL_MOON_00)

    def test_eli500_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 500, 2.9900e-06, JAN_FULL_MOON_00)

    def test_eli862_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 862, 2.2911e-06, JAN_FULL_MOON_00)

    def test_eli1011_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 1011, 1.8330e-06, JAN_FULL_MOON_00)

    def test_eli1662_uncorrected_Valladolid_20220117_00(self):
        _testValladolidNoCorr(self, 1662, 8.4489e-07, JAN_FULL_MOON_00)

    def test_eli338_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 338, 1.2293e-06, JAN_FULL_MOON_17)

    def test_eli385_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 385, 1.5392e-06, JAN_FULL_MOON_17)

    def test_eli481_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 481, 4.0048e-06, JAN_FULL_MOON_17)

    def test_eli540_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 540, 3.8516e-06, JAN_FULL_MOON_17)

    def test_eli879_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 879, 2.7428e-06, JAN_FULL_MOON_17)

    def test_eli1020_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 1020, 2.2381e-06, JAN_FULL_MOON_17)

    def test_eli1654_uncorrected_Valladolid_20220117_17(self):
        _testValladolidNoCorr(self, 1654, 1.0285e-06, JAN_FULL_MOON_17)

    def test_eli336_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 336, 5.4519e-10, FEB_NEW_MOON_00)

    def test_eli380_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 380, 1.0571e-09, FEB_NEW_MOON_00)

    def test_eli440_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 440, 1.8795e-09, FEB_NEW_MOON_00)

    def test_eli500_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 500, 2.7575e-09, FEB_NEW_MOON_00)

    def test_eli862_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 862, 2.1584e-09, FEB_NEW_MOON_00)

    def test_eli1011_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 1011, 7.5293e-10, FEB_NEW_MOON_00)

    def test_eli1662_uncorrected_Valladolid_20220202_00(self):
        _testValladolidNoCorr(self, 1662, 7.9724e-10, FEB_NEW_MOON_00)


def main():
    global prop_error
    if len(sys.argv) > 1:
        prop_error = float(sys.argv[1]) / 100
    else:
        prop_error = DEFAULT_PROP_ERROR
    unittest.main(argv=["first-arg-is-ignored"], exit=False)


if __name__ == "__main__":
    main()
