# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

import os.path
import unittest

import astropy.units as u

from astro_metadata_translator.tests import MetadataAssertHelper

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class LcogtTestCase(unittest.TestCase, MetadataAssertHelper):
    """Test LCOGT translations."""

    datadir = os.path.join(TESTDIR, "data")

    def test_lcogt_translator(self):
        test_data = (
            (
                "fitsheader-lcogt-science.yaml",
                dict(
                    telescope="LCOGT 1m0-13 telescope",
                    instrument="fa01",
                    boresight_rotation_coord="sky",
                    dark_time=201.15662 * u.s,
                    detector_exposure_id=22938825,
                    detector_name="1",
                    detector_unique_name="S1",
                    detector_group="S",
                    detector_num=25,
                    detector_serial="DB-54",
                    exposure_id=229388,
                    exposure_group="229388",
                    exposure_time=20.068 * u.s,
                    focus_z=0 * u.mm,
                    group_counter_end=229388,
                    group_counter_start=229388,
                    has_simulated_content=False,
                    object="2023 PC",
                    observation_counter=229388,
                    observation_id="ct4m20130901t060255",
                    observation_type="science",
                    observation_reason="science",
                    observing_day=20230820,
                    physical_filter="w",
                    pressure=819.6660989 * u.hPa,
                    relative_humidity=31.1,
                    science_program="LCO2023B-005",
                    temperature=8.4 * u.deg_C,
                    visit_id=229388,
                    wcs_params=dict(max_sep=1.5),
                ),
            ),
            (
                "fitsheader-lcogt-2mcalib.yaml",
                dict(
                    telescope="LCOGT 1m0-13 telescope",
                    instrument="fa01",
                    boresight_rotation_coord="sky",
                    dark_time=201.15662 * u.s,
                    detector_exposure_id=22938825,
                    detector_name="1",
                    detector_unique_name="S1",
                    detector_group="S",
                    detector_num=25,
                    detector_serial="DB-54",
                    exposure_id=229388,
                    exposure_group="229388",
                    exposure_time=20.068 * u.s,
                    focus_z=0 * u.mm,
                    group_counter_end=229388,
                    group_counter_start=229388,
                    has_simulated_content=False,
                    object="2023 PC",
                    observation_counter=229388,
                    observation_id="ct4m20130901t060255",
                    observation_type="science",
                    observation_reason="science",
                    observing_day=20230820,
                    physical_filter="w",
                    pressure=819.6660989 * u.hPa,
                    relative_humidity=31.1,
                    science_program="LCO2023B-005",
                    temperature=8.4 * u.deg_C,
                    visit_id=229388,
                    wcs_params=dict(max_sep=1.5),
                ),
            ),
        )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(file, dir=self.datadir, **expected)

if __name__ == "__main__":
    unittest.main()
