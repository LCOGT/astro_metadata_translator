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
    maxDiff = None

    def test_lcogt_translator(self):
        test_data = (
            (
                "fitsheader-lcogt-science.yaml",
                dict(
                    telescope="LCOGT 1m0-13 Telescope",
                    instrument="fa01",
                    boresight_rotation_coord="sky",
                    dark_time=20.0* u.s,
                    detector_exposure_id=69679343269,
                    detector_name="DB-54",
                    detector_unique_name="DB-54",
                    detector_group=None,
                    detector_num=1,
                    detector_serial="DB-54",
                    exposure_id=69679343269,
                    exposure_group=None,
                    exposure_time=20.068 * u.s,
                    focus_z=0 * u.mm,
                    group_counter_end=69679343269,
                    group_counter_start=69679343269,
                    has_simulated_content=False,
                    object="2023 PC",
                    observation_counter=69679343269,
                    observation_id="69679343269",
                    observation_type="science",
                    observation_reason="science",
                    observing_day=20230820,
                    physical_filter="w",
                    pressure=819.6660989 * u.hPa,
                    relative_humidity=31.1,
                    science_program="LCO2023B-005",
                    temperature=8.4 * u.deg_C,
                    visit_id=541961607,
                    wcs_params=dict(max_sep=5.0),
                ),
            ),
            (
                "fitsheader-lcogt-2mcalib.yaml",
                dict(
                    telescope="LCOGT 2m0-01 Telescope",
                    instrument="ep04",
                    boresight_rotation_coord="sky",
                    dark_time=0 * u.s,
                    detector_exposure_id=69579051378,
                    detector_name="1234567",
                    detector_unique_name="1234567",
                    detector_group=None,
                    detector_num=1,
                    detector_serial="1234567",
                    exposure_id=69579051378,
                    exposure_group=None,
                    exposure_time=0.0 * u.s,
                    focus_z=0 * u.mm,
                    group_counter_end=69579051378,
                    group_counter_start=69579051378,
                    has_simulated_content=False,
                    object="N/A",
                    observation_counter=69579051378,
                    observation_id="69579051378",
                    observation_type="bias",
                    observation_reason="calibration",
                    observing_day=20230816,
                    physical_filter="gp",
                    pressure=712.5 * u.hPa,
                    relative_humidity=72.8,
                    science_program="calibrate",
                    temperature=6.8 * u.deg_C,
                    visit_id=541178539,
                    wcs_params=dict(max_sep=1.5),
                ),
            ),
        )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(file, dir=self.datadir, **expected)

if __name__ == "__main__":
    unittest.main()
