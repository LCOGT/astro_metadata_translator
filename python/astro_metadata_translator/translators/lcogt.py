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

"""Metadata translation code for LCOGT FITS headers."""

from __future__ import annotations

__all__ = ("LcogtTranslator",)

import logging
import posixpath
import re
from collections.abc import Iterator, Mapping, MutableMapping
from typing import TYPE_CHECKING, Any

import astropy.units as u
from astropy.coordinates import Angle, EarthLocation
from astropy.io import fits

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .fits import FitsTranslator
from .helpers import altaz_from_degree_headers, is_non_science, tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.time

log = logging.getLogger(__name__)


class LcogtTranslator(FitsTranslator):
    """Metadata translator for LCOGT standard headers."""

    name = "LCOGT"
    """Name of this translation class"""

    supported_instrument = "LCOGT"
    """Supports the LCOGT instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "LCOGT")
    """Default resource path root to use to locate header correction files."""

    # LCOGT has no rotator, and the instrument angle on sky is set to +Y=East,
    # +X=South which we define as a 90 degree rotation and an X-flip.
    _const_map = {
        "boresight_rotation_angle": Angle(90 * u.deg),
        "boresight_rotation_coord": "sky",
        "detector_num": 0,
        "detector_group": None,
        "exposure_group": None
    }

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "physical_filter": "FILTER",
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "focus_z": ("FOCOBOFF", dict(unit=u.mm)),
        "dark_time": ("REQTIME", dict(unit=u.s)),
        "boresight_airmass": ("AIRMASS", dict(checker=is_non_science)),
        "visit_id": "BLKUID",
        "object": "OBJECT",
        "science_program": "PROPID",
        "detector_serial": "DETECTID",
        "detector_name": "DETECTID",
        "detector_unique_name": "DETECTID",
        "instrument": ("INSTRUME", dict(default="Sinistro")),
        # Ensure that reasonable values are always available
        "relative_humidity": ("WMSHUMID", dict(default=40.0, minimum=0, maximum=100.0)),
        "temperature": ("WMSTEMP", dict(unit=u.deg_C, default=10.0, minimum=-10.0, maximum=40.0)),
        "pressure": ("WMSPRES", dict(unit=u.hPa, default=771.611, minimum=700.0, maximum=1150.0)),
    }

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
        supplied header.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        if (
            cls.is_keyword_defined(header, "ORIGIN")
            and cls.is_keyword_defined(header, "DATADICV")
            and cls.is_keyword_defined(header, "HDRVER")
            and "LCOGT-DIC" in header["DATADICV"]
            and "LCOGT" in header["ORIGIN"]
            and "LCOGT-HDR" in header["HDRVER"]
        ):
            return True
        return False

    @cache_translation
    def to_altaz_begin(self) -> astropy.coordinates.AltAz:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(
            self, (("ALTITUDE", "AZIMUTH"),), obstime=self.to_datetime_begin()
        )

    @cache_translation
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord:
        # Docstring will be inherited. Property defined in properties.py
        radecsys = ("RADECSYS",)
        radecpairs = (("RA", "DEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.hourangle, u.deg))

    @cache_translation
    def to_datetime_begin(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        value = self._from_fits_date_string(
            self._header["DATE-OBS"], scale=self._header["TIMESYS"].lower()
        )
        self._used_these_cards("DATE-OBS", "TIMESYS")
        return value

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py

        # Take a guess by adding on the exposure time
        value = self.to_datetime_begin() + self.to_exposure_time()
        return value

    # @cache_translation
    # def to_location(self) -> EarthLocation:
        # """Calculate the observatory location.

        # Returns
        # -------
        # location : `astropy.coordinates.EarthLocation`
            # An object representing the location of the telescope.
        # """

        # for long_key, lat_key, hgt_key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
            # if self.are_keys_ok([long_key, lat_key, hgt_key]):
                # print("Keys OK")
                # value = EarthLocation.from_geocentric(self._header[long_key], self._header[lat_key], self._header[hgt_key])
                # self._used_these_cards(long_key, lat_key)
                # break
        # return value

    @cache_translation
    def to_observation_type(self) -> str:
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        obstype = self._header["OBSTYPE"].strip().lower()
        self._used_these_cards("OBSTYPE")
        if obstype == "expose":
            return "science"
        return obstype

    @cache_translation
    def to_exposure_id(self) -> int:
        """Calculate exposure ID.

        Returns
        -------
        id : `int`
            ID of exposure.
        """
        mol_num = self._header["MOLUID"]
        frame_num = self._header["FRAMENUM"]
        value = int(str(mol_num) + str(frame_num))
        self._used_these_cards("MOLUID", "FRAMENUM")
        return value

    @cache_translation
    def to_detector_exposure_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id()

    @cache_translation
    def to_observation_id(self) -> int:
        """Return the observation id.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return str(self.to_exposure_id())

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the lifetime exposure number.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()

    @cache_translation
    def to_telescope(self) -> str:
        """Calculate the telescope name.

        Returns
        -------
        telname : `str`
            Telescope name. Combination of 'LCOGT' + TELESCOP + 'Telescope.
        """
        telname = self._header["TELESCOP"].strip().lower()
        self._used_these_cards("TELESCOP")
        telname = f"LCOGT {telname} Telescope"
        return telname

    @classmethod
    def determine_translatable_headers(
        cls, filename: str, primary: MutableMapping[str, Any] | None = None
    ) -> Iterator[MutableMapping[str, Any]]:
        """Given a file return all the headers usable for metadata translation.

        MegaPrime files are multi-extension FITS with a primary header and
        each detector stored in a subsequent extension.  MegaPrime uses
        ``INHERIT=F`` therefore the primary header will always be ignored
        if given.

        Parameters
        ----------
        filename : `str`
            Path to a file in a format understood by this translator.
        primary : `dict`-like, optional
            The primary header obtained by the caller. This is sometimes
            already known, for example if a system is trying to bootstrap
            without already knowing what data is in the file. Will be
            ignored.

        Yields
        ------
        headers : iterator of `dict`-like
            Each detector header in turn. The supplied header will never be
            included.

        Notes
        -----

        """
        # Since we want to scan many HDUs we use astropy directly to keep
        # the file open rather than continually opening and closing it
        # as we go to each HDU.
        with fits.open(filename) as fits_file:
            for hdu in fits_file:
                # Astropy <=4.2 strips the EXTNAME header but some CFHT data
                # have two EXTNAME headers and the CCD number is in the
                # second one.
                if hdu.name == "PRIMARY":
                    continue

                if hdu.name.startswith("SCI"):
                    # It may only be some data files that are broken so
                    # handle the expected form.
                    yield hdu.header
                    break

