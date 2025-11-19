#!/usr/bin/env python3
# Script that tests against some CAELIS RIMO cases of the 1088 ph.

import unittest
import os
from datetime import datetime, timedelta as td

import numpy as np

from rimopy import eli, esi

KERNELS_PATH = os.path.join(os.path.dirname(__file__), "./kernels")

WLENS = [440.1, 500.1, 675.1, 870, 937.2, 1019.9, 1638.8]
#WLENS = [440, 500, 675, 870, 935, 1020, 1640]
E0_WEHRLI = [
    1.83425,
    1.9247,
    1.51323,
    0.957819,
    0.810338,
    0.713963,
    0.237624,
]
ESI_CALC = esi.ESICalculatorCustom(WLENS, E0_WEHRLI, True)
ELI_SETTS = eli.ELISettings(True, True, True, eli.MissingRCFBehavior.NEAREST)

CAELIS_PATH = os.path.join(os.path.dirname(__file__), "assets/caelis_sample.csv")


def _read_caelis_file():
    with open(CAELIS_PATH, encoding="utf-8") as fp:
        lines = fp.readlines()
    lines = [l.strip() for l in lines[1:]]
    lines = [l.split(",") for l in lines]
    lines = [[datetime.fromisoformat(l[0]), float(l[1]), float(l[2]), float(l[3])] for l in lines]
    return lines


class TestSum(unittest.TestCase):
    def _testTeide(self, wavelengths, expecteds, date):
        date = date + td(seconds=30)
        date = datetime.isoformat(date)
        ep = eli.EarthPoint(28.2725, -16.6427777, date, 3570)
        sims = eli.get_irradiance(
            wavelengths,
            earth_data=ep,
            kernels_path=KERNELS_PATH,
            esi_calc=ESI_CALC,
            eli_settings=ELI_SETTS,
        )
        np.testing.assert_allclose(sims, expecteds)

    def test_all(self):
        lines = np.array(_read_caelis_file())
        for i in range(len(lines) // 7):
            vals = np.array(lines[i * 7 : (i + 1) * 7])
            vals = vals[vals[:, 1].argsort()].T
            self.assertTrue(np.all(vals[0] == vals[0][0]))
            self._testTeide(vals[3].astype(float), vals[2].astype(float), vals[0][0])


def main():
    unittest.main(exit=False)


if __name__ == "__main__":
    main()
