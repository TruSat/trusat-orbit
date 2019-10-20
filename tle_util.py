#!/usr/bin/env python
""" tle_util.py - Utilities to import, export, validate and operate on Two-Line Element sets """
# Online TLE update portions from  from stvid/update_tle.py
# https://github.com/cbassa/stvid/blob/master/update_tle.py
# TODO May need facility to access archived TLEs (at SpaceTrack) 

from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version 3.6")
    sys.exit(1)

import re
import os
import copy
from hashlib import md5             # md5 finger printing
from shutil import copyfile
from io import BytesIO
from zipfile import ZipFile
from urllib.request import urlopen
from datetime import date, timedelta, datetime
from math import degrees, pi, pow, radians, sqrt  # Fast/precise math functions
from sgp4.ext import jday

import logging
log = logging.getLogger(__name__)

# import string
from getpass import getpass

from spacetrack import SpaceTrackClient

# Global variables
twopi = 2*pi
nocon = twopi/1440.0


""" TODOs:
 - Implement option for "strict" or "loose" processing of element fields (i.e., allow use of non-compliant TLE formats)

    # Todo 
    # - FILEOP Look at how Brandon provided file-load progress
    # - FILEOP Look at how people typically store TLEs in array
    # - IMPORT Update TLE, fingprint, schemas
    # - QUERY Elements nearest to an EPOCH (greater than, less than, some error bounds)
    # - FILEOP XF6 Insert 42 character names into line 0 of TLE from reference file (mcnames)
    # - FILEOP XF7 Update a tle file from another TLE file
    # - FILEOP XF9 Filter out elements with "DEB" with a space before that string. 
    #    But tles with names containing SPELDA, SYLDA, TANK, DPAF, and COVER
    #    The leading 0 and space are removed from "line 0" of the spacecom 3 line tles.
    # - FILEOP XF10 The leading 0 and space are removed from "line 0"
    # - FILEOP XF11 Sort TLEs into ascending NCat order
    # - FILEOP geteccen Filter mean motion less than 8.0 and eccentricity greater than 0.1
    # - FILEOP getleo Filter mean motion is greater than 5.0
    # - FILEOP getepoch (output epoch of each elset) only an .exe
"""


def extract_zip_to_memory(input_zip):
    """Return contents of zip file(s) in memory"""

    input_zip = ZipFile(input_zip)
    return {name: input_zip.read(name) for name in input_zip.namelist()}

class Error(Exception):
   """Base class for other exceptions"""
   pass

class TLEValueError(Error):
    """Raised when TLEs fail checksum and parameter validity checks."""
    pass

def checksum_tle_line(_line):
    """ Performs TLE-defined checksum on TLE line"""

    check = 0
    for char in _line[:-1]:
        if char.isdigit():
            check += int(char)
        if char == "-":
            check += 1

    _check_val = check % 10

    return(_check_val)

# TODO: Do something intelligent if num is larger than digits
def tle_fmt_int(num, digits=5):
    """ Return an integer right-aligned string with DIGITS of precision, all blank if num=0 
    Ignores sign.
    """
    if num:
        num = abs(num)
    else:
        return " "*digits
    string_int = "{:>{DIGITS}d}".format(num,DIGITS=digits)
    return string_int


def tle_fmt_epoch(EpochDateTime):
    """ Return an Epoch string in TLE format, with a total width of 14 characters

    YYDDD.dddddddd where Midnight Jan 1 2019 is 19001.00000000
    """
    # (year, month, day, hour, minute, second, wday, yday, dst) = EpochDateTime.timetuple()
    # Get just the variables we need - to avoid lint errors for unused variables
    (hour, minute, second, _, yday) = EpochDateTime.timetuple()[3:8]
    YY = EpochDateTime.strftime("%y")
    microsecond = EpochDateTime.microsecond

    frac_days = (yday) + (((microsecond/(3600*1e6)) + (second / 3600) + (minute/60) + hour)/24)
    return "{:2s}{:012.8f}".format(YY,frac_days)

def datetime_from_tle_fmt(tleformat):
    """ Return a datetime value from TLE epoch format

    YYDDD.dddddddd where Midnight Jan 1 2019 is 19001.00000000
    """
    _epoch_year = int(tleformat[0:2])
    if (_epoch_year >= 57):
        _epoch_year = 1900 + _epoch_year
    elif (_epoch_year < 57):
        _epoch_year = 2000 + _epoch_year

    _epoch_day = float(tleformat[2:])
    epoch = datetime(_epoch_year, 1, 1, 0, 0, 0, 0) + timedelta(days = (_epoch_day - 1) )
    return epoch

def fingerprint_file(file):
    """Open, read file and calculate MD5 on its contents"""
    with open(file,'rb') as fd:
        # read contents of the file
        _file_data = fd.read()    
        # pipe contents of the file through
        file_fingerprint = md5(_file_data).hexdigest()
    return file_fingerprint


def fingerprint_line(line):
    """ Creates a unique signature from a line."""
    return md5(line.encode('utf-8')).hexdigest()


def read_tle_decimal(pack):
    """Convert *pack* to decimal value.
    For example, the packed value -12345-6 corresponds to -0.12345 × 10-6
    Empty corresponds to a zero value.
    """
    if pack[0] in ["-", " ", "+"]:
        digits = pack[1:-2].strip()
        val = pack[0] + "." + digits + "e" + pack[-2:]
    else:
        digits = pack[:-2].strip()
        val = "." + digits + "e" + pack[-2:]
    return float(val)


def tle_fmt_decimal_pack(val):
    """Convert decimal value to TLE modified fortran notation:
    For example, the value -0.12345 × 10-6 corresponds to -12345-6 in TLE pack format    

    Written as to require no additional modules (performance is not critical).
    """
    if (val == 0):
        return " 00000-0"
    try:
        exponent_string = "{:e}".format(val)
        (mantissa,exponent) = exponent_string.split('e')
        mantissa = float(mantissa)
        exponent = int(exponent)
        exponent += 1
        mantissa *= 0.1
        significand = assumed_decimal_point(mantissa,5)
        if (val < 0):
            sign = "-"
        else:
            sign = " "
    except:
        log.error("Unable to parse exponent for packed decimal")
        return " "*5
    pack_string = "{SIGN}{VALUE}{EXPONENT:+1d}".format(
        SIGN=sign,
        VALUE=significand,
        EXPONENT=exponent
    )
    return pack_string


def launch_piece_letter_to_number(letters):
    """Convert 24 upper case letter 3-letter code to integer, omitting I (eye) and O (oh) from the alphabet.
    
    The TLE standard is 24 upper case letter 3-letter code, omitting I (eye) and O (oh) from the alphabet,
    with no representation for zero.
    
    The 1st piece of a launch is denoted by 'A', and subsequent pieces 'B', 'C', 'D'... 'Z'.
    The 25th (24*1 + 1) piece would be denoted by 'AA', and subsequent pieces 'AB', 'AC'... 'AZ', 'BA', BB', 'BC',... 'ZZ'.
    The 601st (24*24 + 24 + 1) piece would be denoted by 'AAA', and subsequent pieces, 'AAB', 'AAC'... AZZ', 'BAA', 'BAB'... 'ZZZ'
    This allows for a maximum of 24^3 + 24^2 + 24 pieces, or 14424 pieces for a single launch (ZZZ)

    Receives

    """

    letters = letters.strip()

    # Omit I (eye) and O (oh)
    dictionary = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8,
    'J':9,'K':10,'L':11,'M':12,'N':13,'P':14,'Q':15,'R':16,'S':17,
    'T':18,'U':19,'V':20,'W':21,'X':22,'Y':23,'Z':24}

    base = len(dictionary)
    strlen = len(letters)

    if (9 <= strlen <= 11): # Handle the case where we got passed the full international_designation
        letters = letters[8:].strip()
        strlen = len(letters)

    if strlen == 1:
        number = dictionary[letters]
    elif strlen == 2:
        first_number = dictionary[letters[0]]
        second_number = dictionary[letters[1]]
        number = (first_number * base) + second_number
    elif strlen == 3:
        first_number = dictionary[letters[0]]
        second_number = dictionary[letters[1]]
        third_number = dictionary[letters[2]]
        number = (first_number * base * base) + (second_number * base) + third_number
    else: # No valid data received
        return False
    return number


def delta_TLE(TLE1, TLE2):
    """ Difference of two TLEs - for evaluating parameter changes
    
    Return a TLE with sensible values of TLE2 - TLE1
    Non-numeric values are taken from TLE2
    Intended to allow a comparison of a pre/post-fit TLE pair
    """
    if (TLE1.satellite_number != TLE2.satellite_number):
        log.error("TLEs must be for same object")
        return False

    TLE = TruSatellite()

    TLE.satellite_number = TLE2.satellite_number
    TLE.classification = TLE2.classification
    TLE.designation = TLE2.designation
    TLE.ephemeris_type = TLE2.ephemeris_type
    TLE.element_num = TLE2.element_num
    TLE.epoch_datetime = TLE2.epoch_datetime

    TLE.jdsatepoch = TLE2.jdsatepoch - TLE1.jdsatepoch
    TLE.jdSGP4epoch = TLE2.jdSGP4epoch - TLE1.jdSGP4epoch

    TLE.line1_checksum = False
    TLE.line2_checksum = False

    TLE.mean_motion_derivative = TLE2.mean_motion_derivative - TLE1.mean_motion_derivative
    TLE.mean_motion_sec_derivative = TLE2.mean_motion_sec_derivative - TLE1.mean_motion_sec_derivative
    TLE.bstar = TLE2.bstar - TLE1.bstar

    TLE.inclination_degrees = TLE2.inclination_degrees - TLE1.inclination_degrees
    TLE.raan_degrees = TLE2.raan_degrees - TLE1.raan_degrees
    TLE.eccentricity = TLE2.eccentricity - TLE1.eccentricity
    TLE.arg_perigee_degrees = TLE2.arg_perigee_degrees - TLE1.arg_perigee_degrees
    TLE.mean_anomaly_degrees = TLE2.mean_anomaly_degrees - TLE1.mean_anomaly_degrees
    TLE.mean_motion_orbits_per_day = TLE2.mean_motion_orbits_per_day - TLE1.mean_motion_orbits_per_day
    TLE.orbit_number = TLE2.orbit_number - TLE1.orbit_number

    # Sourced from TLEs, but less useful directly
    TLE._id_launch_year = TLE2._id_launch_year
    TLE._id_launch_num = TLE2._id_launch_num
    TLE._id_launch_piece_letter = TLE2._id_launch_piece_letter
    TLE._epoch_year = TLE2._epoch_year
    TLE._epoch_day = TLE2._epoch_day

    # Derived quantities - MKS units
    TLE.inclination_radians = TLE2.inclination_radians - TLE1.inclination_radians
    TLE.raan_radians = TLE2.raan_radians - TLE1.raan_radians
    TLE.arg_perigee_radians = TLE2.arg_perigee_radians - TLE1.arg_perigee_radians
    TLE.mean_anomaly_radians = TLE2.mean_anomaly_radians - TLE1.mean_anomaly_radians
    TLE.mean_motion_radians_per_second = TLE2.mean_motion_radians_per_second - TLE1.mean_motion_radians_per_second    

    return TLE


# Note: Can't call it "Satellite" as that appears to interfere with Satellite class in python-sgp4
class TruSatellite(object):
    """ Class for TLE data objects and methods

    Do everything you would need to parse and validate an individual TLE to
    prepare it for next steps.

    Special things:
    - tle_file_signature is to reference the parent of this TLE, 
    default is "orphan" unless overridden (e.g., by script-generated TLE)
    - checksum defaults to "flag" to reject TLEs that don't arrive with a valid checksum
    Setting it to anything else overrides the received checksum with the correct checksum.
    """

    # Calculations and constants we want to compute exactly once
    # User can over-ride these with different values if desired
    _XKMPER = 6378.137       # WGS84 Earth Equatorial Radius
    _GE     = 398600.4418    # Earth gravitational constant km3/s2
    _GEsqrt = sqrt(_GE)
    _XKE    = sqrt((3600.0 * _GE) / (pow(_XKMPER,3)))


    def __init__(self, catalog=None, line0=None, line1=None, line2=None, tle_source_filename=None, tle_file_fingerprint=None, strict=True, checksum="flag"):
        self._tle_file = tle_source_filename
        self.line0     = line0
        self.line1     = line1
        self.line2     = line2
        self.tle_file_fingerprint = tle_file_fingerprint
        self._tle_source_filename = tle_source_filename
        self.strict    = strict
        self.checksum  = checksum

        # Variables users would likely want regular access to
        self.satellite_number = None
        self.classification = None  # Note: Set this to "O" (or some less ambiguous character) to indicate source?
        self.designation = None
        self.mean_motion_derivative = None          # TODO orbits_per day - and MKS versions
        self.mean_motion_sec_derivative = None      # TODO orbits_per day - and MKS versions
        self.bstar = None
        self.ephemeris_type = None
        self.element_num = None
        self.line1_checksum = None

        self.inclination_degrees = None
        self.raan_degrees = None
        self.eccentricity = None
        self.arg_perigee_degrees = None
        self.mean_anomaly_degrees = None    # revs / day
        self.mean_motion_orbits_per_day     = None
        self.orbit_number = None
        self.line2_checksum = None

        # Sourced from TLEs, but less useful directly
        self._id_launch_year = None
        self._id_launch_num = None
        self._id_launch_piece_letter = None
        self._epoch_year = None
        self._epoch_day = None

        # Derived quantities - MKS units
        self.inclination_radians = None
        self.raan_radians = None
        self.arg_perigee_radians = None
        self.mean_anomaly_radians = None
        self.mean_motion_radians_per_second = None

        # Derived quantities - other
        self.epoch_datetime = None
        self.epoch_string = None                            # Used for SQL INSERT, as it can't take a python datetime directly
        self.launch_piece_number = None
        self.tle_fingerprint = None                         # Created on import by DB
        self.analyst_object = None
        self.tle_good = None

        # Created on import by DB
        self.tle_id	= None
        self.import_timestamp = None

        # Convenience variables
        _GEsqrt = TruSatellite._GEsqrt
        _XKMPER = TruSatellite._XKMPER


        # If we are being asked to parse a TLE
        if (line1 and line2):
            self.line1 = line1.rstrip()
            self.line2 = line2.rstrip()
            if (line0):
                self.line0 = line0.rstrip()
                self.name  = re.sub("^0 ", "", self.line0)

            # Step through TLE import process
            try:
                self._checksum_tle()
                self._parse_tle()
                self._validity_check_tle()
                self._fingerprint_tle()
                self.derived_values()

            except TLEValueError:
                log.warning("{}: Encountered errors in processing the following TLE block:\t{}\n\t{}\n\t{}".format(self._tle_source_filename, self.line0,self.line1,self.line2))

    def derived_values(self):
        """ Calculate values which are determined from TLE parameters """
        self.epoch_string = self.epoch_datetime.isoformat(timespec='microseconds')

        (year, month, day, hour, minute, second) = self.epoch_datetime.timetuple()[:6]
        microseconds = int(self.epoch_datetime.strftime('%f'))
        sec_with_microseconds = second + microseconds/1.0E6

        self.jdsatepoch = jday(year, month, day, hour, minute, sec_with_microseconds)
        self.jdSGP4epoch = self.jdsatepoch - 2433281.5

        self.inclination_radians  = radians(self.inclination_degrees)
        self.raan_radians         = radians(self.raan_degrees)
        self.arg_perigee_radians  = radians(self.arg_perigee_degrees)
        self.mean_anomaly_radians = radians(self.mean_anomaly_degrees)
        self.mean_motion_radians_per_second = 2 * pi * self.mean_motion_orbits_per_day / 86400

        xpdotp   =  1440.0 / (2.0 *pi)  #  229.1831180523293
        self.mean_motion_radians_per_minute = self.mean_motion_orbits_per_day / xpdotp 

        if (self.designation and not self._id_launch_year):
            self._id_launch_year = int(self.designation[2:4])
        if (self.designation and not self._id_launch_num):
            self._id_launch_num = int(self.designation[5:8])
        if (self.designation and not self._id_launch_piece_letter):
            self._id_launch_piece_letter = self.designation[8:].strip()

        self.period = 2*pi/(self.mean_motion_radians_per_second)                            # In seconds
        self.semi_major_axis = pow(TruSatellite._GEsqrt / self.mean_motion_radians_per_second,2/3)  # in km
        self.perigee = self.semi_major_axis*(1 - self.eccentricity) - TruSatellite._XKMPER               # in km
        if(self.perigee < 0):
            log.warning("{}: Perigee {:0f} intersects the Earth.".format(self._tle_source_filename, self.perigee))

        self.apogee  = self.semi_major_axis*(1 + self.eccentricity) - TruSatellite._XKMPER               # in km

        # Add in the calculations for the non SGP4 things here...

        # TODO: Determine if we should create these if not defined from source information, although these are not stored in the DB
        # _epoch_day
        # _epoch_year

    def _parse_tle(self):
        """Parse fields in TLE data"""

        # Parse line 0
        self.name_long = self.line0[0:24].rstrip()
        self.name_long  = re.sub("^0 ", "", self.name_long)

        # Parse line 1
        try:
            self.satellite_number = int(self.line1[2:7])
        except ValueError:
            log.warning("{}: TLE Sat # NaN for {}".format(self._tle_source_filename, self.name_long))
            self.tle_good = False
 
        self.classification  = self.line1[7]

        self._id_launch_year         = self.line1[ 9:11]
        self._id_launch_num          = self.line1[11:14]
        self._id_launch_piece_letter = self.line1[14:17].strip()

        if (80000 <= self.satellite_number <= 89999):
            # Analyst object
            self.designation = self.line1[ 9:17].rstrip()
            self.analyst_object = True

            # It's probably the case that analyst objects don't have this detail defined
            # Doublecheck with T.S. Kelso
            self._id_launch_year = None
            self._id_launch_num = None
            self._id_launch_piece_letter = None
        elif (not self._id_launch_year.isspace() and not self._id_launch_num.isspace() and not self._id_launch_piece_letter.isspace()):
            try:
                self._id_launch_year  = int(self._id_launch_year)
                if (self._id_launch_year >= 57):
                    self._id_launch_year = 1900 + self._id_launch_year
                elif (self._id_launch_year < 57):
                    self._id_launch_year = 2000 + self._id_launch_year
            except ValueError:
                log.warning("{}: TLE Launch year NaN for {}".format(self._tle_source_filename, self.satellite_number))
                if(self.strict):
                    self.tle_good = False

            try:
                self._id_launch_num   = int(self._id_launch_num)
            except ValueError:
                log.warning("{}: TLE Launch number NaN for {}".format(self._tle_source_filename, self.satellite_number))
                if(self.strict):
                    self.tle_good = False

            else:
                try: 
                    assert self._id_launch_piece_letter.isupper() == True
                    assert self._id_launch_piece_letter.find('I') < 0
                    assert self._id_launch_piece_letter.find('O') < 0
                    self.launch_piece_number = launch_piece_letter_to_number(self._id_launch_piece_letter)                                                            
                except AssertionError: 
                    # This is a real error (and not an analyst object) if we're this far...
                    log.warning("{}: TLE invalid characters in launch piece field\n\t{}".format(self._tle_source_filename, self.line1))
                    if(self.strict):
                        self.tle_good = False

                try:
                    designation = "{:4d}-{:>03d}{:<3s}".format(self._id_launch_year,
                                                                    self._id_launch_num,
                                                                    self._id_launch_piece_letter)
                    self.designation = designation.rstrip()
                except (ValueError):
                    if(self.strict):
                        self.tle_good = False
                    else:
                        self.designation = self.line1[ 9:17].rstrip()
        else:
            log.warning("{}: Non-analyst object {} with no valid launch info.".format(self._tle_source_filename, self.satellite_number))

        self.epoch_datetime = datetime_from_tle_fmt(self.line1[18:32])

        self.mean_motion_derivative     = float(self.line1[33:43])
        self.mean_motion_sec_derivative = read_tle_decimal(self.line1[44:52])
        self.bstar                      = read_tle_decimal(self.line1[53:61])

        try:
            self.ephemeris_type = int(self.line1[62])
        except ValueError:
            self.ephemeris_type = 0
        self.element_num  = int(self.line1[64:68])

        # Parse line 2
        # Figure out where to do the error / type checking on these
        self.inclination_degrees  = float(self.line2[ 8:16])
        self.raan_degrees         = float(self.line2[17:25])
        self.eccentricity =   int(self.line2[26:33]) * 10 ** -7
        self.arg_perigee_degrees  = float(self.line2[34:42])
        self.mean_anomaly_degrees = float(self.line2[43:51])
        self.mean_motion_orbits_per_day  = float(self.line2[52:63])
        self.orbit_number =   int(self.line2[63:68])

        self.derived_values()
        
        # Create a tuple convenient for handing to python-SGP4 satrec variable
        # SGP expecting:
        #  Angles in radians
        #  Angle rates in radians per _minute_
        #  epoch time in days from jan 0, 1950. 0 hr
        self.satrec = [self.satellite_number, self.jdSGP4epoch, self.bstar, self.eccentricity, 
                       self.arg_perigee_radians, self.inclination_radians, self.mean_anomaly_radians, 
                       self.mean_motion_radians_per_minute, self.raan_radians]


    def _checksum_tle(self):
        """ Performs TLE-defined checksum on TLE """

        _check_val = [None, None, None]
        _linenum = 1

        for _line in [self.line1, self.line2]:
            _check_val[_linenum] = checksum_tle_line(_line)

            if (_check_val[_linenum] != int(_line[-1])) and (self.checksum!="fix"):
                log.warning("{}: TLEChecksumError: {}\n\tGot {} - should be: {}".format(self._tle_source_filename, _line,_line[-1],_check_val[_linenum]))
                self.tle_good = False
            elif (_line == 1):
                self.line1_checksum = _check_val[_line]
            elif (_line == 2):
                self.line2_checksum = _check_val[_line]
            _linenum += 1   


    def _validity_check_tle(self):
        """ Perform format and range-checking tests described at: https://celestrak.com/columns/v04n03/ """
        # TODO: Make sure the line1 object number matches the line2 object number 
        _thisyear = date.today().year
        _validity_errors = {}
        if (not (0 <= self.raan_degrees <= 360)):
             _validity_errors["RAAN_degrees"] = self.raan_degrees
        if (not (0 <= self.arg_perigee_degrees <= 360)):
             _validity_errors["arg_perigee_degrees"] = self.arg_perigee_degrees
        if (not (0 <= self.mean_anomaly_degrees <= 360)):
             _validity_errors["mean_anomaly_degrees"] = self.mean_anomaly_degrees
        if (not (0 <= self.inclination_degrees <= 180)):
             _validity_errors["inclination_degrees"] = self.inclination_degrees
        if (not (0 <= self.eccentricity < 1)):
             _validity_errors["eccentricity"] = self.eccentricity
            # Might need some allowance for future epochs here

        # TODO Determine a reasonable date in the future to check against epoch
        if (not (datetime(1957,10,4,0,0,0,0) <= self.epoch_datetime )):
             _validity_errors["epoch"] = self.epoch_datetime

        # Launch year might be None or string
        try:
            if (not (date(1957,1,1) <= date(self._id_launch_year,1,1) <= date.today() )):
                _validity_errors["launch_year"] = self._id_launch_year
        except TypeError:
            pass    # Should only be None types

        if (not (0 <= self.ephemeris_type <= 5)):
             _validity_errors["ephemeris_type"] = self.ephemeris_type
        if (len(_validity_errors) == 0):
            self.tle_good = True
        else: # 
            log.warning("{}: TLE failed validity error for sat_num {} ({})".format(self._tle_source_filename, self.satellite_number,_validity_errors))
            log.info("  {}\n  {}\n  {}".format(self.line0,self.line1,self.line2))
            self.tle_good = False


    def _fingerprint_tle(self):
        """ Creates a unique signature from a TLE.
        
        Incorporates those parts of the TLE that contribute to the orbit properties."""

        _line1fragment = self.line1[19:64] # Epoch year through Ephemeris type
        _line2fragment = self.line2[9:64]  # Inclination through Mean Motion

        _TLE_fingerprint_string = self.line1 + self.line2
                    
        self.tle_fingerprint = md5(_TLE_fingerprint_string.encode('utf-8')).hexdigest()

    def make_tle_lines(self):
        """ Creates TLE line1 and line2 from TLE class variables """
        # FIXME - Figure out why the roundtrip of tle_epoch dosen't precisely match an untouched TLE from the database

        tle_epoch = tle_fmt_epoch(self.epoch_datetime)
        eo_string = assumed_decimal_point(self.eccentricity,7)

        packed_designation = "{LAUNCH_YEAR:02d}{LAUNCH_NUM:03d}{LAUNCH_PIECE_LETTER:<3s}".format(
            LAUNCH_YEAR = self._id_launch_year,
            LAUNCH_NUM  = self._id_launch_num,
            LAUNCH_PIECE_LETTER = self._id_launch_piece_letter)

        # TODO: Deal with First Derivative xno, Second derivative xno, Bstar
        line1 = "1 {:5d}{:1} {:<8s} {:14s} {:<10.8f} {:8s} {:8s} {:1d}{:4s}00".format(
            self.satellite_number,
            self.classification,
            packed_designation,
            tle_epoch,
            self.mean_motion_derivative,
            tle_fmt_decimal_pack(self.mean_motion_sec_derivative),
            tle_fmt_decimal_pack(self.bstar),
            self.ephemeris_type,
            tle_fmt_int(self.element_num,digits=4)
            )

        line2 = "2 {:05d} {:8.4f} {:8.4f} {:7s} {:8.4f} {:8.4f} {:11.8f}{:5s}00".format(
            self.satellite_number, 
            self.inclination_degrees, 
            self.raan_degrees, 
            eo_string, 
            self.arg_perigee_degrees, 
            self.mean_anomaly_degrees, 
            self.mean_motion_orbits_per_day,
            tle_fmt_int(self.orbit_number,digits=5)
            )

        self.line1_checksum = checksum_tle_line(line1)
        self.line2_checksum = checksum_tle_line(line2)

        self.line0 = "0 {:<22s}".format(self.name) # Maintain the 24 Character Celestrak Standard, leading 0
        self.line1 = line1[:67] + "{:02d}".format(self.line1_checksum)
        self.line2 = line2[:67] + "{:02d}".format(self.line2_checksum)

class TLEFile(object):
    """TLEFile: Class for TLE file operations
    
        Returns an Dict of TruSatellite objects with NORAD id as KEY
        Note: This method of key-association will only store the first record for each NORAD-id it finds in the file.
    """

    def __init__(self, tle_file, strict=True, parse=True):
        self.tle_file = tle_file
        self._tle_fd = None                 # TLE file descriptor
        self.tle_file_fingerprint = None
        self._TLEs = []
        self.strict = strict
        self.parse = parse
        self._tle_basename = os.path.basename(self.tle_file)
        self.Satellites = {}

        if (strict):
            # FROM: https://www.orekit.org/static/jacoco/org.orekit.propagation.analytical.tle/TLE.java.html
            self._tle_line1_re = re.compile ('^1 [ 0-9]{5}[A-Z] [ 0-9]{5}[ A-Z]{3} [ 0-9]{5}[.][ 0-9]{8} (?:(?:[ 0+-][.][ 0-9]{8})|(?: [ +-][.][ 0-9]{7})) [ +-][ 0-9]{5}[+-][ 0-9] [ +-][ 0-9]{5}[+-][ 0-9] [ 0-9] [ 0-9]{4}[ 0-9]')
            self._tle_line2_re = re.compile ('^2 [ 0-9]{5} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{7} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{2}[.][ 0-9]{13}[ 0-9]')
        else:
            self._tle_line1_re = re.compile ('^1 ')
            self._tle_line2_re = re.compile ('^2 ')

        self.load_tles()

        if (parse):
            self.parse_tles()


    # FIXME: This can probably be eliminated in favor of general fingerprint_file()
    def fingerprint_file(self):
        """Open, read file and calculate MD5 on its contents"""
        with open(self.tle_file,'rb') as fd:
            # read contents of the file
            _file_data = fd.read()    
            # pipe contents of the file through
            self.tle_file_fingerprint = md5(_file_data).hexdigest()
        return self.tle_file_fingerprint


    def load_tles(self):
        """Read TLE data.
        
        Perform an MD5 checksum of the source file, and store with the Class variables.
        """

        def _line012(_l0, _l1, _l2):
            """ Mini routine to format the TLE set """
            _l0 = _l0.rstrip()
            if _l0.startswith('0 '): # Spacetrack 3-line format
                name = _l0[2:].rstrip()
            else:
                name = _l0.rstrip()
            return(name, _l1.rstrip(), _l2.rstrip())

        if (self.tle_file_fingerprint is None):
            self.fingerprint_file()

        l0 = l1 = ""
        tlecount = 0
        tlefileline = 0 
        with open(self.tle_file,'r') as tlefd:
            for l2 in tlefd:
                l2 = l2.rstrip()
                tlefileline += 1

                simplematch1 = (l1.startswith('1 ') and len(l1) >= 69)
                simplematch2 = (l2.startswith('2 ') and len(l2) >= 69)

                match1 = self._tle_line1_re.search(l1)                    
                match2 = self._tle_line2_re.search(l2)

                if (match1 and match2 and self.strict):
                    (name, line1, line2) = _line012(l0, l1, l2)
                    tlecount += 1
                    self._TLEs.append([name,line1,line2])
                elif ((simplematch1 and simplematch2) and not self.strict):
                    (name, line1, line2) = _line012(l0, l1, l2)
                    tlecount += 1
                    self._TLEs.append([name,line1,line2])
                elif ( (self.strict and simplematch1 and simplematch2) and (match1==None or match2==None)):
                    if not (match1 or match2):
                        log.warning("{}: Strict record structure checks failed for both TLE lines at file line: {}".format(self._tle_basename,tlefileline))
                        log.info("  {}\n  {}\n  {}\n".format(l0,l1,l2))
                    elif not match1:
                        log.warning("{}: Strict record structure checks failed for TLE line 1 at file line {}".format(self._tle_basename,tlefileline))
                        log.info("  {}\n  {}\n  {}\n".format(l0,l1,l2))
                    elif not match2:
                        log.warning("{}: Strict record structure checks failed for TLE line 2 at file line {}".format(self._tle_basename,tlefileline))
                        log.info("  {}\n  {}\n  {}\n".format(l0,l1,l2))
                    else:
                        log.warning("{}: Should not be here.".format(self._tle_basename))

                l0 = l1
                l1 = l2
        log.info("Read {} TLEs".format(tlecount))


    # TODO: Consider calling this parse_tles_unique, and storing only the most recent record provided
    # TODO: Create parse_tles_keep_dupes() which just creates an array with input order
    def parse_tles(self):
        log.info("Parsing...")
        for (line0, line1, line2) in self._TLEs:
            sat = TruSatellite(line0=line0, 
                            line1=line1, 
                            line2=line2, 
                            tle_file_fingerprint=self.tle_file_fingerprint,
                            tle_source_filename=self._tle_basename)
            # FIXME: Note that this method of storing satellites will not deal with multiple entries with the same satnum
            # Causes problems for satfit reading sat.txt files with multiple TLEs (will only show the last one)
            self.Satellites[sat.satellite_number] = sat
        return self.Satellites


def assumed_decimal_point(num_less_than_one, digits=7):
    """ Return a string with DIGITS of precision, with the decimal point removed 
    Ignores sign.
    """
    num = abs(num_less_than_one)
    string_num = "{0:.{DIGITS}f}".format(num,DIGITS=digits)
    return(string_num[2:])


def make_tle(*, name="None", ssn, desig="0000000", epoch_datetime, xincl, xnodeo, eo, omegao, xmo, xno, deg=True, quiet=False):
    """ write TLE to output file and to screen """

    # ssn               Spacecraft number
    # desig             International Designation
    # epoch_datetime    Epoch in pythonDATETIME format
    # xincl             inclination
    # xnodeo            RAAN
    # eo                eccentricity
    # omegao            argument of perigee
    # xmo               mean anomaly
    # xno               mean motion
    # TODO: Find standard variable names for First Derivative xno, Second derivative xno, Bstar

    line0 = None
    tle_epoch = tle_fmt_epoch(epoch_datetime)
    eo_string = assumed_decimal_point(eo,7)

    if (deg==False):
        xincl  = degrees(xincl)
        xnodeo = degrees(xnodeo)
        omegao = degrees(omegao)
        xmo    = degrees(xmo)

    # //   sprintf(bstar_string, "%12.4e", bstar*10);
    # //   bstar_fract[0] = bstar_string[0]; // sign
    # //   bstar_fract[1] = bstar_string[1];
    # //   bstar_fract[2] = bstar_string[3];
    # //   bstar_fract[3] = bstar_string[4];
    # //   bstar_fract[4] = bstar_string[5];
    # //   bstar_fract[5] = bstar_string[6];
    # //   bstar_fract[6] = '\0';
    # //   bstar_exp[0] = bstar_string[8];
    # //   bstar_exp[1] = bstar_string[11];
    # //   bstar_exp[2] = '\0';

    # //   double xns = 2160 * bstar * nn * c2;

    # //   sprintf(line1, "1 %05dU %-8s %014.8f %.8f  00000-0 %6s%2s 0    00"
    # //            ,ssn, desig, tle, xns, bstar_fract, bstar_exp);

    if name:
        line0  = "0 {:22s}".format(name) 
        if not quiet:
            print("{:s}".format(line0))

    line1 = "1 {:5d}U {:<8s} {:14s} 0.00000073  00000-0  50000-4 0    00".format(ssn,desig,tle_epoch)
    # TODO: Deal with First Derivative xno, Second derivative xno, Bstar
    line2 = "2 {:05d} {:8.4f} {:8.4f} {:7s} {:8.4f} {:8.4f} {:11.8f}    00".format(
            ssn, xincl, xnodeo, eo_string, omegao, xmo, xno)

    line1 = line1[:68] + str(checksum_tle_line(line1))
    line2 = line2[:68] + str(checksum_tle_line(line2))

    if not quiet:
        print("{:s}".format(line1))
        print("{:s}".format(line2))
    return(line0, line1, line2)


# TODO: Carefully consider if it makes sense to make classification="T" the default - if 3rd parties are calling this, do we want to "own it"?
def make_tle_from_SGP4_satrec(satrec, classification="T"):
    """ Make TLE record from python-SGP4 satrec variable

    Classification defaults to "T" (for TruSat) unless otherwise specified
    
    Input:
        satrec      python-SGP4 Satellite() Class variable

    Output:
        TLE         tle_util TruSatellite() Class variable (with TLE lines)
    """
    TLE = TruSatellite()

    TLE.line0 = satrec.line0

    TLE.sat_name = TLE.name	            = re.sub("^0 ", "", satrec.line0)
    TLE.satellite_number	            = satrec.satnum
    TLE.classification		            = classification  
    TLE.designation			            = satrec.intldesg
    TLE.epoch_datetime	                = satrec.epoch_datetime
    TLE.mean_motion_derivative		    = satrec.ndot
    TLE.mean_motion_sec_derivative	    = satrec.nddot
    TLE.bstar			                = satrec.bstar
    TLE.ephemeris_type	                = satrec.ephtype
    TLE.element_set_number	            = satrec.elnum
    TLE.inclination_degrees             = degrees(satrec.inclo)
    TLE.inclination_radians	            = satrec.inclo
    TLE.raan_degrees		            = degrees(satrec.nodeo)
    TLE.raan_radians		            = satrec.nodeo
    TLE.eccentricity		            = satrec.ecco
    TLE.arg_perigee_degrees	            = degrees(satrec.argpo)
    TLE.arg_perigee_radians	            = satrec.argpo
    TLE.mean_anomaly_degrees            = degrees(satrec.mo)
    TLE.mean_anomaly_radians            = satrec.mo
    TLE.mean_motion_orbits_per_day      = satrec.no_kozai / nocon
    TLE.mean_motion_radians_per_minute  = satrec.no_kozai
    TLE.mean_motion_radians_per_second  = satrec.no_kozai / 60
    TLE.orbit_number			        = satrec.revnum     # TODO: May need to calculate this based on period and time from epoch
    """ Further notes on Rev number:
    From: https://www.celestrak.com/columns/v04n03/#FAQ02
    The period from launch to the first ascending node is considered to be Rev 0 and Rev 1 begins when the first ascending node is reached. 
    Since many element sets are generated with epochs that place the satellite near its ascending node, 
    it is important to note whether the satellite has reached the ascending node when calculating subsequent rev numbers.
    """

    TLE.launch_piece_number	= launch_piece_letter_to_number(TLE.designation)

    TLE.derived_values()
    TLE.make_tle_lines()
    TLE._fingerprint_tle()

    return TLE


def append_tle_file(file_out, line0, line1, line2):
    try:
        with open(file_out, "a") as fp:
            fp.write("\n")
            if(line0):
                fp.write("{:s}\n".format(line0))
            fp.write("{:s}\n".format(line1))
            fp.write("{:s}\n".format(line2))
            fp.close()
    except IOError:
        log.warning("Can not open file: {}".format(file_out))    


def update_from_online(tle_path):
    log.info("Updataing TLEs in {}".format(tle_path))

    now = datetime.datetime.utcnow()
    time = now.strftime("%Y%m%d_%H%M%S")

    # Get Space Track TLEs
    catalog_tle = os.path.join(tle_path, 'catalog.tle')
    st = SpaceTrackClient(identity=cfg.get('Credentials', 'st-username'),
                        password=cfg.get('Credentials', 'st-password'))

    data = st.tle_latest(iter_lines=True, epoch='>now-30',
                        ordinal=1, format='3le')

    with open(catalog_tle, 'w') as fp:
        for line in data:
            # Fix missing leading zeros
            line = re.sub("^1     ", "1 0000", line)
            line = re.sub("^2     ", "2 0000", line)
            line = re.sub("^1    ", "1 000", line)
            line = re.sub("^2    ", "2 000", line)
            line = re.sub("^1   ", "1 00", line)
            line = re.sub("^2   ", "2 00", line)
            line = re.sub("^1  ", "1 0", line)
            line = re.sub("^2  ", "2 0", line)
            fp.write(line + '\n')

    copyfile(catalog_tle, os.path.join(tle_path, time + '_catalog.txt'))

    # Get classified TLEs
    resp = urlopen("http://www.prismnet.com/~mmccants/tles/classfd.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    classfd_tle = os.path.join(tle_path, 'classfd.tle')

    content = None
    outsize = 0
    with open(classfd_tle, 'rb') as infile:
        content = infile.read()
    with open(classfd_tle, 'wb') as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b'\n')

    copyfile(classfd_tle, os.path.join(tle_path, time + '_classfd.txt'))

    # Get int TLEs
    resp = urlopen("http://www.prismnet.com/~mmccants/tles/inttles.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    int_tle = os.path.join(tle_path, 'inttles.tle')

    content = None
    outsize = 0
    with open(int_tle, 'rb') as infile:
        content = infile.read()
    with open(int_tle, 'wb') as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b'\n')

    copyfile(int_tle, os.path.join(tle_path, time + '_inttles.txt'))

    # Create bulk catalog
    catalogs = [catalog_tle, classfd_tle]
    with open(os.path.join(tle_path, 'bulk.tle'), 'w') as outfile:
        for fname in catalogs:
            with open(fname) as infile:
                outfile.write(infile.read())