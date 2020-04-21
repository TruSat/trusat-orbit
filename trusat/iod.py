#!/usr/bin/env python3
""" iod.py - Utilities for importing, validating, and operating on IOD/RDE/UK positional formatting formats """
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
	print("This script requires Python version 3.6")
	sys.exit(1)

import os
import re
from datetime import datetime, timedelta

from math import asin, sin, cos, atan2

import logging
log = logging.getLogger(__name__)

from trusat.tle_util import launch_piece_letter_to_number, assumed_decimal_point

# REGEXP for valid angle string format with content
# Note that the columns containing the angle format are evaluated for validity with these REGEX
# Rules of the pattern match:
# 1. Force width of first field, column position of sign
# 2. Match only digits and spaces
# 3. No leading space, only trailing space
# 4. Trailing space must not have any other characters

""" Match IOD angle format, return Angle1 and Angle2 (with sign) groups
Regex101 description link - https://regex101.com/r/COxsWQ/6
"""
angle_content_IOD_re = re.compile(r"""
	 
	^ 				     # Start at the beginning of the string
	(?=[\d\s]{7}[+-])    # Lookahead for 7 decimal or spaces through to the sign
	(\d+\s*?) 	 	     # (Angle 1 Capture GROUP) Match 1+ digits followed by 0+ spaces
	(				     # Start Angle 2 Capture GROUP
		[+-]		     # Match sign
		(?=[\d\s]{1,6}$) # Lookahead for 1+ digits followed by 0+ spaces, up to 6 digits total
		(?!\d+\s+\d+$)   # Negative lookahead - don't match 1+ digits followed by 1+ spaces followed by 1+ more digits
		\d+\s*?			 # Match 1+ digits followed by 0+ spaces
	)					 # End Angle 2 Capture GROUP
	""",re.VERBOSE) # pylint: disable=anomalous-backslash-in-string

""" Match UK angle format, return Angle1 and Angle2 (with sign) groups
Regex101 description link - # https://regex101.com/r/COxsWQ/5
"""
angle_content_UK_re  = re.compile(r"""
	^ 				     # Start at the beginning of the string
	(?=[\d\s]{8}[+-])    # Lookahead for 8 decimal or spaces through to the sign
	(\d+\s*?) 	 	     # (Angle 1 Capture GROUP) Match 1+ digits followed by 0+ spaces
	(				     # Start Angle 2 Capture GROUP
		[+-]		     # Match sign
		(?=[\d\s]{1,7}$) # Lookahead for 1+ digits followed by 0+ spaces, up to 7 digits total
		(?!\d+\s+\d+$)   # Negative lookahead - match 1+ digits followed by 1+ spaces followed by 1+ more digits
		\d+\s*?			 # Match 1+ digits followed by 0+ spaces
	)					 # End Angle 2 Capture GROUP
	""",re.VERBOSE) # pylint: disable=anomalous-backslash-in-string

""" Match RDE angle format, return Angle1 and Angle2 (with sign) groups
Regex101 description link - https://regex101.com/r/COxsWQ/4
"""
angle_content_RDE_re  = re.compile(r"""
	^ 				     # Start at the beginning of the string
	(?=[\d\s]{6}[+-])    # Lookahead for 6 decimal or spaces through to the sign
	(\d+\s*?) 	 	     # (Angle 1 Capture GROUP) Match 1+ digits followed by 0+ spaces
	(				     # Start Angle 2 Capture GROUP
		[+-]		     # Match sign
		(?=[\d\s]{1,6}$) # Lookahead for 1+ digits followed by 0+ spaces, up to 6 digits total
		(?!\d+\s+\d+$)   # Negative lookahead - match 1+ digits followed by 1+ spaces followed by 1+ more digits
		\d+\s*?			 # Match 1+ digits followed by 0+ spaces
	)					 # End Angle 2 Capture GROUP
	""",re.VERBOSE) # pylint: disable=anomalous-backslash-in-string

def get_angle_content_IOD(AngleString):
	""" Validates and returns pair of angles from IOD-formatted AngleString, False for each angle if not valid. """
	match = angle_content_IOD_re.search(AngleString)	
	if (match):
		return ( match.group(1), match.group(2) )
	else:
		return (False, False)


def get_angle_content_UK(AngleString):
	""" Validates and returns pair of angles from UK-formatted AngleString, False for each angle if not valid. """
	match = angle_content_UK_re.search(AngleString)	
	if (match):
		return ( match.group(1), match.group(2) )
	else:
		return (False, False)


def get_angle_content_RDE(AngleString):
	""" Validates and returns pair of angles from RDE-formatted AngleString, False for each angle if not valid. """
	match = angle_content_RDE_re.search(AngleString)	
	if (match):
		return ( match.group(1), match.group(2) )
	else:
		return (False, False)


def radec_from_azel(az, el, latgd, lst):
	""" From David Vallado's astIOD.cpp
	/*------------------------------------------------------------------------------
	*
	*                           procedure radec_azel
	*
	* this procedure converts right ascension declination values with
	*   azimuth, and elevation.  notice the range is not defined because
	*   right ascension declination only allows a unit vector to be formed.
	*
	*  author        : david vallado                  719-573-2600   22 jun 2002
	*
	*  inputs          description                    range / units
	*    rtasc       - right ascension                0.0 to 2pi rad
	*    decl        - declination                    -pi/2 to pi/2 rad
	*    lst         - local sidedouble time            -2pi to 2pi rad
	*    latgd       - geodetic latitude              -pi/2 to pi/2 rad
	*
	*  outputs       :
	*    az          - azimuth                        0.0 to 2pi rad
	*    el          - elevation                      -pi/2 to pi/2 rad
	*
	*  locals        :
	*    lha         - local hour angle               -2pi to 2pi rad
	*    sinv        - sine value
	*    cosv        - cosine value
	*
	*  coupling      :
	*    arcsin      - arc sine function
	*    atan2       - arc tangent function that resolves quadrant ambiguites
	*
	*  references    :
	*    vallado       2013, 265, alg 27
	-----------------------------------------------------------------------------*/
	http://www.convertalot.com/celestial_horizon_co-ordinates_calculator.html
	"""
	decl = asin(sin(el) * sin(latgd) + cos(el) * cos(latgd) * cos(az))

	sinv = -(sin(az) * cos(el) * cos(latgd))  / (cos(latgd) * cos(decl))
	cosv = (sin(el) - sin(latgd) * sin(decl)) / (cos(latgd) * cos(decl))
	lha = atan2(sinv, cosv)
	rtasc = lst - lha
	return(rtasc, decl)


def azel_from_radec(rtasc, decl, latgd, lst):
	lha = lst - rtasc
	el = asin(sin(decl) * sin(latgd) + cos(decl) * cos(latgd) * cos(lha))
	sinv = -sin(lha) * cos(decl) * cos(latgd) / (cos(el) * cos(latgd))
	cosv = (sin(decl) - sin(el) * sin(latgd)) / (cos(el) * cos(latgd))
	az = atan2(sinv, cosv)
	return(az, el)


def launch_piece_number_to_letter(piece_num):
    """ Converts UK/RDE format 2-digit launch piece to three-letter TLE standard, left-justified, space-padded.

    The TLE standard is 24 upper case letter 3-letter code, omitting I (eye) and O (oh) from the alphabet,
    with no representation for zero.
    
    The 1st piece of a launch is denoted by 'A', and subsequent pieces 'B', 'C', 'D'... 'Z'.
    The 25th (24*1 + 1) piece would be denoted by 'AA', and subsequent pieces 'AB', 'AC'... 'AZ', 'BA', BB', 'BC',... 'ZZ'.
    The 601st (24*24 + 24 + 1) piece would be denoted by 'AAA', and subsequent pieces, 'AAB', 'AAC'... AZZ', 'BAA', 'BAB'... 'ZZZ'
    This allows for a maximum of 24^3 + 24^2 + 24 pieces, or 14424 pieces for a single launch (ZZZ)
    """

    # Zero is not a valid result, so we should never use the zero index
    dictionary='!ABCDEFGHJKLMNPQRSTUVWXYZ'
    piece_letters=''
    
    # FIXME: Figure out what the right thing to do is if we're passed 0 (Assign to 1=A?)
    x = int(piece_num)
    if (x == 0):
        x = 1

    # 14424 = 24*24*24 + 24*24 + 24
    # Just rail it high to preserve the format width
    if x > 14424:
        x = 14424
        log.warning("Exceeded maximum value (14424) for launch piece number")

    while x>0:
        x,idx = divmod(x,24)
        # With no zero value, zero remainder means we are at the maximum value in the alphabet
        # Not the first of the next!
        if(idx==0):
            idx = 24
            x-=1
        piece_letters = dictionary[idx] + piece_letters			

    return "{:<s}".format(piece_letters)


def right_zero_pad(val,length=8):
	""" Right-zero-pad short-form angle strings with zeros.
		This reduces amount of error checking required, and makes decimal conversions more consistent
		This will also make the IOD/UK/RDE angle lengths uniform (UK/RDE adds one more digit of precision on the right)
	"""
	val = val.rstrip()
	zpadval = val + "0"*(length-len(val))
	return zpadval


def angle_from_HHMMSSss(HHMMSSss):
	""" Returns a decimal angle provided a string in 'right ascension HHMMSSss format'.  Assumed positive value. """
	HHMMSSss = right_zero_pad(HHMMSSss,8)
	try:
		HH = int(HHMMSSss[0:2])
	except ValueError:
		HH = 0 # We should probably just declare it invalid...

	try:
		MM = int(HHMMSSss[2:4])
	except ValueError:
		MM = 0

	try:
		SS = int(HHMMSSss[4:6])
	except ValueError:
		SS = 0

	try:
		ss = int(HHMMSSss[6:8])
	except ValueError:
		ss = 0

	HHMMSSss_angle = 15.0*((HH) + MM/60.0 + SS / 3600.0 + ss*(0.01/3600.0))
	return HHMMSSss_angle

	
def angle_from_HHMMmmmm(HHMMmmmm):
	""" Returns a decimal angle provided a string in 'right ascension HHMMmmmm format'.  Assumed positive value. """
	HHMMmmmm = right_zero_pad(HHMMmmmm,8)
	try:
		HH = int(HHMMmmmm[0:2])
	except ValueError:
		HH = 0 # We should probably just declare it invalid...

	try:
		MM = int(HHMMmmmm[2:4])
	except ValueError:
		MM = 0

	try:
		mmmm = int(HHMMmmmm[4:8])
	except ValueError:
		mmmm = 0

	HHMMmmmm_angle = 15.0*((HH) + MM/60.0 + mmmm/600000.0)	
	return HHMMmmmm_angle
			

def angle_from_DDDMMSSs(DDDMMSSs):
	""" Returns a decimal angle provided a string in right ascension DDDMMSSs format.

	Input can be in the forms:
		DDDMMSSs (unsigned right ascension)
		+DDMMSSs (signed declination / elevation)
	"""
	DDDMMSSs = right_zero_pad(DDDMMSSs)
	if (DDDMMSSs[0] in ("+","-")):
		SIGN = DDDMMSSs[0]
		SIGN = int(SIGN + "1")
		DDDMMSSs = "0" + DDDMMSSs[1:]
	else:
		SIGN = 1.0

	try:
		DDD = int(DDDMMSSs[0:3])
	except ValueError:
		DDD = 0 # We should probably just declare it invalid...

	try:
		MM = int(DDDMMSSs[3:5])
	except ValueError:
		MM = 0

	try:
		SS = int(DDDMMSSs[5:7])
	except ValueError:
		SS = 0

	try:
		s = int(DDDMMSSs[7])
	except ValueError:
		s = 0

	DDDMMSSs_angle = SIGN*(DDD + MM/60.0 + SS/3600.0 + s*0.1/3600.0)
	return DDDMMSSs_angle


def angle_from_DDDMMmmm(DDDMMmmm):
	""" Returns a decimal angle provided a string in right ascension DDDMMmmm format

	Input can be in the forms:
		DDDMMmmm (unsigned right ascension)
		+DDMMmmm (signed declination / elevation)
	"""
	DDDMMmmm = right_zero_pad(DDDMMmmm)
	if (DDDMMmmm[0] in ("+","-")):
		SIGN = DDDMMmmm[0]
		SIGN = int(SIGN + "1")
		DDDMMmmm = "0" + DDDMMmmm[1:]
	else:
		SIGN = 1.0

	try:
		DDD = int(DDDMMmmm[0:3])
	except ValueError:
		DDD = 0 # We should probably just declare it invalid...

	try:
		MM = int(DDDMMmmm[3:5])
	except ValueError:
		MM = 0

	try:
		mmm = int(DDDMMmmm[5:8])
	except ValueError:
		mmm = 0

	DDDMMmmm_angle = SIGN*(DDD + MM/60.0 + mmm*(0.001/60.0))
	return DDDMMmmm_angle


def angle_from_DDDddddd(DDDddddd):
	""" Returns a decimal angle provided a string in right ascension DDDddddd format.
	
	Input can be in the forms:
		DDDddddd (unsigned right ascension)
		+DDddddd (signed declination / elevation)
	"""
	DDDddddd = right_zero_pad(DDDddddd)
	if (DDDddddd[0] in ("+","-")):
		SIGN = DDDddddd[0]
		SIGN = int(SIGN + "1")
		DDDddddd = "0" + DDDddddd[1:]
	else:
		SIGN = 1.0

	try:
		DDD = int(DDDddddd[0:3])
	except ValueError:
		DDD = 0 # We should probably just declare it invalid...

	try:
		ddddd = int(DDDddddd[3:8])
	except ValueError:
		ddddd = 0

	DDDddddd_angle =  SIGN*(DDD + (ddddd * 0.00001))
	return DDDddddd_angle


class Angle:
	"""Observed Angle

	# IOD - First four are OWTG system, minus one digit of precision
	#                    00000000001111
	#                    01234567890123
	# Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
	# 		 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
	# 		 3: RA/DEC = HHMMmmm+DDdddd MX   (MX in degrees of arc)
	# 		 4: AZ/EL  = DDDMMSS+DDMMSS MX   (MX in seconds of arc)
	# 		 5: AZ/EL  = DDDMMmm+DDMMmm MX   (MX in minutes of arc)
	# 		 6: AZ/EL  = DDDdddd+DDdddd MX   (MX in degrees of arc)
	# 		 7: RA/DEC = HHMMSSs+DDdddd MX   (MX in degrees of arc)

	# UK (RDE uses Code 1 with no sub-seconds)	
	#  		      Column 0000000000111111
	# 				     0123456789012345
	#   Code 1: RA/DEC = HHMMSSss+DDMMSSs
	# 	     2: RA/DEC = HHMMmmmm+DDMMmmm
	# 	     3: RA/DEC = HHMMmmmm+DDddddd
	# 	     4: AZ/EL  = DDDMMSSs DDMMSSs (elevation corrected for refraction)
	# 	     5: AZ/EL  = DDDMMmmm DDMMmmm (elevation corrected for refraction)
	# 	     6: AZ/EL  = DDDddddd DDddddd (elevation corrected for refraction)
	# 	     7: AZ/EL  = DDDMMSSs DDMMSSs (elevation not corrected for refraction)
	# 	     8: AZ/EL  = DDDMMmmm DDMMmmm (elevation not corrected for refraction)
	# 	     9: AZ/EL  = DDDddddd DDddddd (elevation not corrected for refraction)

	# TODO: Figure out what to do about elevation / refraction

	   Given a packed Angle string of two angle pairs matching the format above,
	   create one of the following formats.
	   RA =   HH.ddddd
	   DEC = DDD.ddddd

	   AZ =  DDD.ddddd
	   EL =   DD.ddddd

	   Where MX (Uncertainty) is not in degrees of arc, convert to degrees of arc

	The IOD and UK angle codes (1-6) are functionaly equivalent,
	noting the difference in precision.

	IOD Format 7 is unique to it.

	UK Formats 7-9 are identical in formatting to 4-6, noting the 
	considerations for atmospheric refraction.

	RDE Format 1 is the same as IOD/UK format without the sub-seconds
	"""
	EpochDict = {
		0 : 0,		# TODO: of Date (need to implement this "J-NOW" somewhere)
		1 : 1855,
		2 : 1875,
		3 : 1900,
		4 : 1950,
		5 : 2000,
		6 : 2050
	}

	def __init__(self,AngleFormatCode,EpochCode,AngleString,Uncertainty,UncertaintyString,RecordFormat="IOD"):
		self.AngleFormatCode   = AngleFormatCode
		self.EpochCode 		   = EpochCode
		self.AngleString 	   = AngleString
		self.Uncertainty 	   = Uncertainty
		self.UncertaintyString = UncertaintyString
		self.RecordFormat 	   = RecordFormat

		self.Epoch  = None
		self.Angle1 = None
		self.Angle2 = None
		self.RA     = 0.0
		self.DEC    = 0.0
		self.AZ     = 0.0
		self.EL     = 0.0
		self.ValidAngle = False

		self.Epoch = self.EpochDict[self.EpochCode]

		if (self.RecordFormat is "IOD"):
			(self.Angle1,self.Angle2) = get_angle_content_IOD(self.AngleString)
		elif (self.RecordFormat is "RDE"):
			(self.Angle1,self.Angle2) = get_angle_content_RDE(self.AngleString)
		elif (self.RecordFormat is "UK"):
			(self.Angle1,self.Angle2) = get_angle_content_UK(self.AngleString)
		else:
			log.error("Did not find a valid Record Format ({}) specified for angle unpacking.".format(self.RecordFormat))

		if (self.AngleFormatCode and self.Angle1 and self.Angle2):
		
			# 1:    IOD RA/DEC = HHMMSSs0+DDMMSS0 MX   (MX in seconds of arc)
			# 1:     UK RA/DEC = HHMMSSss+DDMMSSs SSSs (seconds of arc)
			# 1:    RDE RA/DEC = HHMMSS00+DDMMSS0 SSS  (seconds of arc)
			# Note, the "0" in these formats, throughout, are the ones we padded onto the IOD/RDE formats
			if (self.AngleFormatCode == 1):
				self.RA  = angle_from_HHMMSSss(self.Angle1)
				self.DEC = angle_from_DDDMMSSs(self.Angle2)

				if (self.RecordFormat is "IOD"):
					# Convert uncertainty in seconds of arc to fractional degrees
					self.Uncertainty *= (1/3600)
				elif (self.RecordFormat is "RDE"):
					self.Uncertainty *= (1/3600)
				elif (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0:3].lstrip() + "." + temp[3])/3600
					except ValueError:
						self.Uncertainty = 0

			# IOD 2: RA/DEC = HHMMmmm0+DDMMmm0 MX   (MX in minutes of arc)
			#  UK 2: RA/DEC = HHMMmmmm+DDMMmmm MMmm (minutes of arc)
			elif (self.AngleFormatCode == 2):
				self.RA  = angle_from_HHMMmmmm(self.Angle1)
				self.DEC = angle_from_DDDMMmmm(self.Angle2)

				if (self.RecordFormat is "IOD"):
					# Convert uncertainty in minutes of arc to fractional degrees
					self.Uncertainty *= (1/60)
				elif (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0:2].lstrip() + "." + temp[2:4])/60
					except ValueError:
						self.Uncertainty = 0

			# 3: IOD RA/DEC = HHMMmmm0+DDdddd0 MX   (MX in degrees of arc, no need to convert)
			# 3:  UK RA/DEC = HHMMmmmm+DDddddd Dddd (degrees of arc)
			elif (self.AngleFormatCode == 3):
				self.RA  = angle_from_HHMMmmmm(self.Angle1)
				self.DEC = angle_from_DDDddddd(self.Angle2)

				if (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0].lstrip() + "." + temp[1:4])
					except ValueError:
						self.Uncertainty = 0

			# 4: IOD AZ/EL  = DDDMMSS0+DDMMSS0 MX   (MX in seconds of arc)
			# 4:  UK AZ/EL  = DDDMMSSs DDMMSSs SSSs (seconds of arc) (elevation corrected for refraction)
			elif (self.AngleFormatCode == 4):
				self.AZ = angle_from_DDDMMSSs(self.Angle1)
				self.EL = angle_from_DDDMMSSs(self.Angle2)

				if (self.RecordFormat is "IOD"):
					# Convert uncertainty in seconds of arc to fractional degrees
					self.Uncertainty *= (1/3600)
				elif (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0:3].lstrip() + "." + temp[3])/3600
					except ValueError:
						self.Uncertainty = 0

			# 5: IOD AZ/EL  = DDDMMmm0+DDMMmm0 MX   (MX in minutes of arc)
			# 5:  UK AZ/EL  = DDDMMmmm+DDMMmmm MMmm (minutes of arc) (elevation corrected for refraction)
			elif (self.AngleFormatCode == 5):
				self.AZ = angle_from_DDDMMmmm(self.Angle1)
				self.EL = angle_from_DDDMMmmm(self.Angle2)

				if (self.RecordFormat is "IOD"):
					# Convert uncertainty in minutes of arc to fractional degrees
					self.Uncertainty *= (1/60)
				elif (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0:2].lstrip() + "." + temp[2:4])/60
					except ValueError:
						self.Uncertainty = 0

			# 6: IOD AZ/EL  = DDDdddd0+DDdddd0 MX   (MX in degrees of arc, no need to convert)
			# 6:  UK AZ/EL  = DDDddddd+DDddddd Dddd (degrees of arc) (elevation corrected for refraction)
			elif (self.AngleFormatCode == 6):
				self.AZ = angle_from_DDDddddd(self.Angle1)
				self.EL = angle_from_DDDddddd(self.Angle2)

				if (self.RecordFormat is "UK"):
					temp = right_zero_pad(self.UncertaintyString,length=4)
					try:
						self.Uncertainty = float(temp[0].lstrip() + "." + temp[1:4])
					except ValueError:
						self.Uncertainty = 0

			# 7:   RA/DECvv = HHMMSSs+DDdddd MX   (MX in degrees of arc, no need to convert)
			elif (self.AngleFormatCode == 7 and self.RecordFormat == "IOD"):
				self.RA  = angle_from_HHMMSSss(self.Angle1)
				self.DEC = angle_from_DDDddddd(self.Angle2)

			# TODO: Add in the remaining format codes for UK (7-9)
			else:
				log.error("Did not find an Angle case with AngleFormatCode '{}'".format(self.AngleFormatCode))
				self.ValidAngle = False
		if ( (self.RA and self.DEC) or (self.AZ and self.EL) ):
			self.ValidAngle = True
		else:
			self.ValidAngle = False
# End of the Angle parsing cases


def DateTime_frompacked(DateTimeString,format_type="IOD"):
	""" Unpack the datetime format from IOD strings 
	
	Of the string format:
		  00000000001111111111
	      01234567890123456789
	WANT: YYYYMMDDHHMMSSssssss 
	IOD:  YYYYMMDDHHMMSSsss
	UK:     YYMMDDHHMMSSssss
	RDE:    YYMMDDHHMMSS.ss

	Also, deal with cases where users provide rounded-up/over-the-max values for HH:MM:SS fields.
	"""
	# Parse YEAR ourselves, as datetime works on a 2000-year boundary, whereas TLE dates work on a [19,20]57 year boundary
	if (format_type in ["RDE","UK"]):
		YEAR = int(DateTimeString[0:2])
		if (YEAR < 57):		# Year of Sputnik launch, first cataloged object
			# Make the input string of consistent between formats for the rest of the parsing
			DateTimeString = '20' + DateTimeString 
		else:
			DateTimeString = '19' + DateTimeString 

	# Zero pad the datestring out to microsecond precision
	DateTimeString = right_zero_pad(DateTimeString,length=20)

	# Deal with decimal point in RDE format
	if (format_type == "RDE"):
		if (DateTimeString[14]=="."):
			SUBSEC = DateTimeString[15:].rstrip()
			MICROSECOND = right_zero_pad(SUBSEC,6)
			DateTimeString = DateTimeString[:14] + MICROSECOND
		else:
			DateTimeString = DateTimeString[:14] + "000000"

	YEAR   = int(DateTimeString[0:4])
	MONTH  = int(DateTimeString[4:6])
	DAY    = int(DateTimeString[6:8])
	HOUR   = int(DateTimeString[8:10])
	MINUTE = int(DateTimeString[10:12])
	SECOND = int(DateTimeString[12:14])
	MICROSECOND = int(DateTimeString[14:20])

	# Handle inputs with more than the maximum than datetime.strptime() can handle
	# This is fromt the users' software rounding up
	carry_over_time = 0
	if (SECOND > 59):
		carry_over_time += 60
		SECOND -= 60
	if (MINUTE > 59):
		carry_over_time += 60*60
		MINUTE -= 60
	if (HOUR > 23):
		carry_over_time += 24*60*60
		HOUR -= 24

	DateTimeString = f"{YEAR:04d}{MONTH:02d}{DAY:02d}{HOUR:02d}{MINUTE:02d}{SECOND:02d}{MICROSECOND:06d}"

	formatstring = '%Y%m%d%H%M%S%f'
	DateTimeUnpacked = datetime.strptime(DateTimeString,formatstring)
	DateTimeUnpacked = DateTimeUnpacked + timedelta(seconds=carry_over_time)

	return DateTimeUnpacked


class IOD:
	"""IOD/Visual Observation (includes UK and RDE encoded data)

	   NameCode, or NameString indicates packed data requiring further processing (to disentagle variable names)

       # Array of variable names for IOD parsing.
       # NameCode, or NameString indicates packed data requiring further processing (to disentagle variable names)
	"""

	# FIXME: There's some type of validation to be done making sure that ObjectNumber and InternationalDesignation are consistent with Celestrak SATCAT (or each other)
	def __init__(self, line=None):				# (Element #, start_col, end_col)
		self.line = line
		self.UserString = None 					# FIXME: A preMVP hack before we're parsing it fully.
		self.ObjectNumber = None				# (0, 0, 4)  	Object number (NORAD)
		self.LaunchYear = None					# (1, 6, 7)  	Launch year. Implied century.
		self.InternationalDesignation = None	# (2, 9, 14)  	International Designation (COSPAR)
		self.Station = 9999						# (3, 16, 19)  	Four digit station number
		self.StationStatusCode = None			# (4, 21, 21)  	Station status code
		self.StationStatus = None				# 				Only in the IOD format
		self.DateTimeString = None				# (5, 23, 39)  	DateTime string YYYYMMDDHHMMSSsss
		self.DateTime = None
		self.TimeUncertainty = None				# (6, 41, 42)  	Time Uncertainty, expressed as MX, Evaluated as M*10E(X-8).
		self.TimeStandardCode = None			#               In the UK, RDE formats
		self.AngleFormatCode = None				# (7, 44, 44)  	Angle format code
		self.EpochCode = None
		self.Epoch = None						# (8, 45, 45)  	Epoch code
		self.RA = None
		self.DEC = None
		self.AZ = None
		self.EL = None
		self.ValidPosition = 1					# Every IOD starts life with the potential to be valid - assume valid until proven otherwise
		self.PositionUncertainty = None			# (10, 62, 63) 	Positional uncertainty, expressed as MX, Evaluated as M*10E(X-8).
		self.OpticalCode = None					# (11, 65, 65) 	Optical behavior code
		self.VisualMagnitude = None				# (12, 66, 69) 	Visual magnitude. Implied decimal point.
		self.MagnitudeUncertainty = None		# (13, 71, 72) 	Magnitude uncertainty. Implied decimal point.
		self.FlashPeriod = None					# (14, 74, 79) 	Flash period in seconds. Implied decimal point.
		self.VisualMagnitude_high = None		#               In the UK, RDE formats
		self.VisualMagnitude_low = None			#               In the UK, RDE formats
		self.Remarks = None
		self.message_id = None
		self.IODType = None						# IOD, UK or RDE
		self.iod_string = line					# duplicate to self.line, but DB calls it iod_string FIXME: ? (rename it above, and throughout iod.py)

		self.obs_id = None						# Created automatically on DB import
		self.obsFingerPrint = None				# Created automatically on DB import
		self.import_timestamp = None			# Created automatically on DB import
		self.submitted = None					# Created automatically on DB import

	# TODO: The following functions are to-do
	def calc_AZ_EL_from_RA_DEC(self):
		# Create as a class function call, but put the actual workings as a module function
		# so it can be accessed separately
		pass

	def calc_RA_DEC_from_AZ_EL(self):
		# Create as a class function call, but put the actual workings as a module function
		# so it can be accessed separately
		pass
# end of class IOD


# Don't worry about the lint warnings on these: https://github.com/PyCQA/pylint/issues/1451
# Previous version for revert if necessary
#iod_format_re = re.compile ('^(\d{5}\s\d{2}\s\d{3}\D*)*\s{1,16}(\d{4}\s)(\D)\s(\d{8})(\d{0,9}$|\d{0,9}\s+\d{2}\s*\d*\s*\d*\s*[-+]*\d*\s*\d{0,2}\s[EFIRSXBHPADMNV]{0,1})') # pylint: disable=anomalous-backslash-in-string
# Simple format
#iod_format_re = re.compile ('^(\d{5}\s\d{2}\s\d{3}\D*)*\s{1,16}(\d{4}\s)(\D)\s(\d{8})(\d{0,9}$|\d{0,9}\s+\d\d.*)')

# This one catches all the good lines, and lets a ones with invalid VMAG/FLASH data through
# More description here: https://regex101.com/r/sevtPF/1
new_iod_format_re = re.compile(r"""
    ^	# Start at beginning of string
    (	# BEGIN NORAD/INTL DES GROUP 
		\d{5}\s
	    \d{2}\s
	    \d{3}
	    (?=[A-Z]+\s*)[\D\s]{3}(?<!\s\w)\s
	    |\s{16}	# Or match empty from beginning for 16 characters
	)	# END NORAD/INTL DES GROUP 
    [0-9A-HJ-NP-Z]{4}\s	# Station
    [EGFPBTCO ]\s	# Station status code
    [\d+]{8}			# YYYYMMDD
    (               # BEGIN optional remaining data GROUP
	\d*\s*$|		# End of the record, or
    (?=.{9})\d*\s*?\s # Match the next 9 characters of the time string, digits than optional space
    \d{2}\s			# Time uncertainty
    (				# BEGIN ANGLE DATA BLOCK
    [1-7]			# Valid angle format codes
    [\s0-6]\s		# Valid time epoch codes
    (?=[\d\s*]{7})\d+\s*?	# First group of Angle data
    [+-]					# Sign
    (?=[\d\s*]{6})\d+\s*?\s # Second group of Angle data
    \d{2}					# Position uncertainty
    |\s{20}					# 20 blank characters if no angle data
	) # END ANGLE DATA GROUP		
     (	# BEGIN Optical Behavior GROUP
	  \s[EFIRSXBHPADMNV]	# Valid optical behavior code
      ( # BEGIN visual magnitude GROUP
	   [+-]						# Visual Magnitude sign
       (?=[\d\s*?]{3})\d+\s*?\s	# Visual Magnitude
       (?=[\d\s*?]{2})\d+\s*?\s	# Visual Magnitude uncertainty
    	 (\s+\d+$)?				# Optional Flash period
      )? # END optional visual magnitude GROUP
     )? # END optional optical behavior GROUP
    ) # END optional remaining data GROUP
	""",re.VERBOSE)

def format_test_iod(line=False):
	match = new_iod_format_re.search(line)
	if match:
		return True
	else:
		return False


""" 
Notes on the UK format:

This format is deficient for the following reasons:
- Uses international designator as the identifier (no NORAD number)
- Uses a 2 digit numeric code for launch piece number, as opposed to a 3-letter code in the TLE standard.
- Column 33, "Time Standard Code" appears to be meaningless, as 0 is listed as possible by the example, which lists 1-3 being valid.
- Doesn't require leading zeros (i.e., on positional accuracy numbers)
However, it is superior for the following reason:
- One more digit of precision in reported angles

"""

# Previous version for revert if necessary
#uk_format_re = re.compile ('^\d{17}(?=\d+[ ]*)[ \d]{10}(?=\d+[ ]*)[\d ]{5}[1-3][1-9](?=\d*[ ]*)[\d ]{8}[-+ ](?=\d+[ ]*)[\d ]{7}[ \d]{4}[0-6$]') # pylint: disable=anomalous-backslash-in-string
# Weaker UK format test which matches the beginning of a UK line but NOT the beginning of an RDE line
# uk_format_re = re.compile('^(?=\d+[ ]*?)[\d ]{24}(?=\d)[ \d]{3}') # pylint: disable=anomalous-backslash-in-string
# uk_format_re = re.compile('^(\d{17})') # pylint: disable=anomalous-backslash-in-string

# More description here: https://regex101.com/r/YUKMjI/1
new_uk_format_re = re.compile (r"""
	^(\d{7}[0-9A-HJ-NP-Z]{4}\d{6})	# Intl Desig, Site Number, Date string
	(?=\d[\d\s]{9}).{10}	# Time of observation, must start with at least one digit, digits and spaces can follow
	(?=[\d\s]{5}).{5}		# Time accuracy, can be all blank
	([1-3]					# Time Standard code, 1-3 valid
	[1-9])					# Angle Format Code, 1-9 valid
	(?=\d[\d\s]{7}).{8}	    # AngleString1 - Must start with a digit, followed by 7 more digits or spaces
	([-+\s]					# Sign for ngleString2
	(?=\d[\d\s]{6}).{7})	# AngleString2 - Must start with a digit, followed by 6 more digits or spaces
	((?=[\d\s]{4}).{4})		# Position uncertainty, 4 digits or blanks
	[0-6]					# Epoch of start chart, 0-6 valid
	(\s*?$|					# End of minimum valid record, or continue matching
	(?=[\s]{13}).{13}		# Range and Range accuracy data (not used - accept any 13 characters)
	((?=[\s\d+-]{3}).{3})	# Brightest Visual Magnitude 3 digits, spaces or sign
	((INV|(?=[\s\d+-]{3}).{3}))	# Faintest Visual Magnitude 3 digits, spaces or sign - or the phrase INV if the satellite becomes invisible
	((?=[\d\s]{5}).{5})		# Flash period - 5 digits or spaces
	([SIRFXE ])				# Remarks code, valid characters
	)
	""",re.VERBOSE)

# uk_format_re = re.compile('^\d{17}(?=\d+[ ]*)[ \d]{10}(?=\d+[ ]*)[\d ]{5}[1-3][1-9](?=\d*[ ]*)[\d ]{8}[-+ ](?=\d+[ ]*)[\d ]{7}[ \d]{4}[0-6$][\d ]{8}[\d ]{5}(?=[+-]?\d+\s*)[-+\d ]{3}(?=[+-]?\d+\s*)[-+\d ]{3}(?=[ ]*\d+[ ]*)[\d ]{5}[ SIRFXE]{1}') # pylint: disable=anomalous-backslash-in-string
def format_test_uk(line=False):
	match = new_uk_format_re.search(line)
	if match:
		return True
	else:
		return False

# Previous version (currently in use due to https://github.com/consensys-space/trusat-orbit/issues/5) for revert if necessary
rde_format_re = re.compile('[0-9A-HJ-NP-Z]{4}\s\d{4}\s(?=\d.\d*[ ]*)[\d. ]{4}\d\s\d{3}[0-6]\n(^\d{2}\n(^\d{7}\s\d{6}.\d{2}\s\d{6}[-+ ]\d{6}(?=[- ]\d.\d)[- \d.]{4}(?=[- ]\d.\d)[- \d.]{4}.*\n){2,}){1,}999', re.MULTILINE) # pylint: disable=anomalous-backslash-in-string

# The following regex creates stalling-conditions on bulk import ref: https://github.com/consensys-space/trusat-orbit/issues/5
# The following detects a complete RDE record block
# https://regex101.com/r/AINrwf/1
new_rde_format_re_verbose = re.compile(r"""
[0-9A-HJ-NP-Z]{4}\s		 # Observing site number (4 digits)
\d{4}\s					 # UTC Year and Month of Observation, in YYMM format
(?=\d.\d*?\s*?)[\d. ]{3} # Time Accuracy, in seconds. Format is T.t
[1-3]					 # Time Standard Code, 1-3 valid
[1]\s					 # Position Format Code, Russell always uses 1
\d{3}					 # Position Accuracy, SSS
[0-6]\n					 # Epoch of star chart used to determine position

(?:						 # Beginning of optional repeat data block for daily data
 ^\d{2}.*\n				 # 2 digit day of month, followed by optional space then newline

 (?:					 # Beginning of optional repeat data block for observation lines 
  ^\d{7}\s				 # International Designation - 7 digit string
  \d{6}.\d{2}\s			 # UTC Time of Observation, in HHMMSS.ss
  \d{6}					 # Observed Right Ascension or Azimuth, in HHMMSS format
  [-+ ]\d{6}			 # Observed Declination or Elevation, in DDMMSS format (with sign)
  (?=[- ]\d.\d)[- \d.]{4}        # Brightest Visual Magnitude
  (?=[- ]\d.\d)[- \d.]{4}        # Faintest Visual magnitude
  (?=[- ]\s?\d.*?\d*?)[- \d.]{3} # Flash period
  [SIRFXE ].*\n 				 # Optical behavior code
 ){1,}					 # Require one data line per day
){1,}					 # Require one day per record
999						 # RDE record ends with 999
""",flags=re.MULTILINE|re.VERBOSE) # Need to OR flags to use multiple options

# With all extraneous white space removed, as to avoid the OR'ed flags in new_rde_format_re_verbose
new_rde_format_re_compact = re.compile(r"[0-9A-HJ-NP-Z]{4}\s\d{4}\s(?=\d.\d*?\s*?)[\d. ]{3}[1-3][1]\s\d{3}[0-6]\n(?:^\d{2}.*\n(?:^\d{7}\s\d{6}.\d{2}\s\d{6}[-+ ]\d{6}(?=[- ]\d.\d)[- \d.]{4}(?=[- ]\d.\d)[- \d.]{4}(?=[- ]\s?\d.*?\d*?)[- \d.]{3}[SIRFXE ].*\n){1,}){1,}999",re.MULTILINE)

# The following detects a valid data line within a RDE record block
# https://regex101.com/r/MDKDE6/2
new_rde_data_line_re = re.compile(r"""
  ^\d{7}\s				 # International Designation - 7 digit string
  \d{6}.\d{2}\s			 # UTC Time of Observation, in HHMMSS.ss
  \d{6}					 # Observed Right Ascension or Azimuth, in HHMMSS format
  [-+ ]\d{6}			 # Observed Declination or Elevation, in DDMMSS format (with sign)
  (?=[- ]\d.\d)[- \d.]{4}        # Brightest Visual Magnitude
  (?=[- ]\d.\d)[- \d.]{4}        # Faintest Visual magnitude
  (?=[- ]\s?\d.*?\d*?)[- \d.]{3} # Flash period
  [SIRFXE ].*$ 				 # Optical behavior code
""",re.VERBOSE)
rde_data_line_re = re.compile('\d{7}\s\d{6}.\d{2}\s\d{6}[-+ ]\d{6}(?=[- ]\d.\d)[- \d.]{4}(?=[- ]\d.\d)[- \d.]{4}')


def format_test_rde(block=False):
	match = new_rde_format_re.search(block)
	if match:
		return True
	else:
		return False


# Parse enter blocks of text instead of just lines
# This is the typical use-case, as IOD entries usually come in clusters
def parse_iod_lines(block=False):
	""" Parse a block of text containing IOD-formatted observation records, and return formatted records + count.

	Format documented at: http://www.satobs.org/position/IODformat.html
	"""
	lines = block.split('\n')
	IODentryList = []
	iod_count = 0
	iod_match_count = 0
	iod_block_line = 0
	earliest_time = datetime(1957,10,4,19,28,34,0)
	import_time = datetime.utcnow()
	for line in lines:
		iod_block_line +=1
		line = line.rstrip()
		match = format_test_iod(line)
		if (match):
			iod_match_count +=1
			log.debug("IOD match {} on line {}".format(iod_match_count,iod_block_line))
			IOD_line = IOD(line)
			IOD_line.IODType = "IOD"
			line_length = len(line)

			try:
				IOD_line.ObjectNumber = int(line[0:5])
			except ValueError:
				IOD_line.ObjectNumber = 0

			try:
				LaunchYearTwoDigit = int(line[6:8])
				if (LaunchYearTwoDigit < 57):		# Year of Sputnik launch, first cataloged object
					IOD_line.LaunchYear = 2000 + LaunchYearTwoDigit
				else:
					IOD_line.LaunchYear = 1900 + LaunchYearTwoDigit
				IOD_line.InternationalDesignation = str(IOD_line.LaunchYear) + '-' + line[9:15] # TODO: Figure out if we should strip() this or not

				# If the observer wasn't sure, let the DB trigger look it up from SATCAT
				if ("?" in IOD_line.InternationalDesignation):
					IOD_line.InternationalDesignation = "?"	
			except ValueError:
				IOD_line.InternationalDesignation=""

			try:
				IOD_line.Station = line[16:20]
			except ValueError: # No point in continuing if we don't have a station to report against
				continue

			IOD_line.StationStatusCode = line[21]	# Future TODO - Expand out the short and long description of the codes

			if (line_length >= 41):
				IOD_line.DateTimeString = line[23:40]
				IOD_line.DateTime = DateTime_frompacked(IOD_line.DateTimeString)
			else:
				IOD_line.DateTimeString = line[23:]
				IOD_line.DateTime = DateTime_frompacked(IOD_line.DateTimeString)

			# Flag records from the future and pre history
			if not (earliest_time < IOD_line.DateTime < import_time):
				IOD_line.ValidPosition = -1

			if (line_length >= 43):
				try:
					# Expressed as MX, where M = mantissa, and X = exponent input. Evaluated as M*10E(X-8).
					IOD_line.TimeUncertainty = float(line[41]) * 10 **(int(line[42]) - 8)
				except ValueError:
					IOD_line.TimeUncertainty = 0

			if (line_length >= 45):
				try:
					IOD_line.AngleFormatCode = int(line[44])
				except ValueError:
					IOD_line.AngleFormatCode = -1
					IOD_line.ValidPosition = -1	# Need format code to parse angle

				try:
					EpochCodeString = line[45]
					if (EpochCodeString == " "):
						IOD_line.EpochCode = 0
					else:
						IOD_line.EpochCode = int(line[45])
				except ValueError:
					IOD_line.EpochCode = -1
					IOD_line.ValidPosition = -1	# If EpochCode as read as Zero (not blank), we won't trigger this error
			else:
				IOD_line.AngleFormatCode = -1 
				IOD_line.EpochCode = -1
				IOD_line.ValidPosition = -1

			# TODO: Figure out which of the missing data parts to actually warn users about
			if (line_length>61 and IOD_line.AngleFormatCode >0):
				AngleString = line[47:61]

				try:
					IOD_line.PositionUncertainty = float(line[62]) * 10 **(int(line[63]) - 8)
				except ValueError:
					IOD_line.PositionUncertainty = None # Don't really need this (initialized state)

				if ( (IOD_line.AngleFormatCode>=0) and (IOD_line.EpochCode>=0) and AngleString):
					try:
						angle = Angle(
							IOD_line.AngleFormatCode, 
							IOD_line.EpochCode, 
							AngleString, 
							IOD_line.PositionUncertainty, 
							"",
							"IOD")
						if (angle.ValidAngle):
							IOD_line.AZ = angle.AZ
							IOD_line.EL = angle.EL
							IOD_line.RA = angle.RA
							IOD_line.DEC = angle.DEC
							IOD_line.Epoch = angle.Epoch
							IOD_line.PositionUncertainty = angle.Uncertainty
						else:
							IOD_line.ValidPosition = -1
					except:
						log.error("Problem angle: '{}' - Offending line:".format(AngleString))
						log.error(line)
						IOD_line.ValidPosition = -1
				elif ((IOD_line.AngleFormatCode>=0) and (IOD_line.EpochCode>=0)):
					# Note - Not an error, because the standards provide for "reporting for duty"
					# even if no observations are made.
					log.warning("Valid Angle ({}) and Epoch ({}) codes, but no valid position data. Offending line:".format(IOD_line.AngleFormatCode,IOD_line.EpochCode ))
					log.warning(line)
					IOD_line.ValidPosition = -1
				else:
					IOD_line.ValidPosition = -1
			else:
				IOD_line.ValidPosition = -1


			if (line_length >= 65):
				IOD_line.OpticalCode = line[65]

			if (line_length >= 70):
				# Cols 67-70 Visual Magnitude
				vmag = line[66:70]
				vmag_sign = vmag[0]
				vmag_whole = vmag[1:3]
				vmag_frac = vmag[3]

				try:
					IOD_line.VisualMagnitude = float(vmag_sign + vmag_whole + '.' + vmag_frac)
				except ValueError:
					IOD_line.VisualMagnitude = 99


			if (line_length >= 73):
				# Cols  72- 73: Magnitude uncertainty. Implied decimal point
				# Implied decimal point between cols 72 and 73. 
				vmag_unc = line[71:73]
				try:
					IOD_line.MagnitudeUncertainty = float(vmag_unc[0] + '.' + vmag_unc[1])
				except ValueError:
					IOD_line.MagnitudeUncertainty = 0

			if (line_length >= 80):
				# Cols  75- 80: Flash period in seconds. 
				# Implied decimal point between cols 77 and 78. BLANK IF NO DATA.
				flash_period = line[74:80]
				try:
					IOD_line.FlashPeriod = float(flash_period[0:3] + '.' + flash_period[3:6])
				except ValueError:
					IOD_line.FlashPeriod = 0

			if (line_length >= 80):
				IOD_line.Remarks = line[80:]
					
			iod_count += 1
			IODentryList.append(IOD_line)
	if (len(IODentryList)):
		return IODentryList
	else:
#		log.warning("Error: no data found. Call parse_iod_lines with a properly formatted IOD line comfirmed by format_test_iod")
		return(False)


# Parse enter blocks of text instead of just lines
# This is the typical use-case, as IOD entries usually come in clusters
def parse_uk_lines(block=False):
	""" Parse a block of text containing UK-formatted observation records, and return formatted records + count.

	Format documented at: http://www.satobs.org/position/UKformat.html
	"""
	lines = block.split('\n')
	IODentryList = []
	uk_count = 0
	earliest_time = datetime(1957,10,4,19,28,34,0)
	import_time = datetime.utcnow()
	for line in lines:
		line = line.rstrip()
		match = format_test_uk(line) # Handle this separately? I.E. assume its properly formatted?  Would require more exception handling below.
		if (match):
			UK_line = IOD(line)
			UK_line.IODType = "UK"
			line_length = len(line)

			# Note that UK Format does not specify Object Numbers
			# We're setting to 0, and relying on the SQL Trigger 'add_object_number' to update it on INSERT 
			# From the InternationalDesignation
			# FIXME: Figure out (from Mike McCants) how RDE/UK specify analyst objects
			# FIXME: Or figure out what do do with International Designators that don't map to the Celestrak SATCAT
			UK_line.ObjectNumber = 0

			LaunchYearTwoDigit = int(line[0:2])
			if (LaunchYearTwoDigit < 57):		# Year of Sputnik launch, first cataloged object
				UK_line.LaunchYear = 2000 + LaunchYearTwoDigit
			else:
				UK_line.LaunchYear = 1900 + LaunchYearTwoDigit

			piece_letters = launch_piece_number_to_letter(line[5:7])
			UK_line.InternationalDesignation = str(UK_line.LaunchYear) + '-' + line[2:5] + piece_letters # FIXME: Figure out if this should have trailing space like the IOD format does

			# 'Observing Site Number'
			UK_line.Station = line[7:11]

			UK_line.DateTimeString = line[11:27]
			UK_line.DateTime = DateTime_frompacked(UK_line.DateTimeString,"UK")
			# Flag records from the future and pre history
			if not (earliest_time < UK_line.DateTime < import_time):
				UK_line.ValidPosition = -1

			try:
				UK_line.TimeUncertainty = float(line[27] + '.' + line[28:32].rstrip())
			except ValueError:
				UK_line.TimeUncertainty = 0

			try:
				UK_line.TimeStandardCode = int(line[32])
			except ValueError:
				UK_line.TimeStandardCode = 0

			# Position Format Code
			try:
				UK_line.AngleFormatCode = int(line[33])
			except ValueError:
				log.error("Valid Angle Format Code not specified.")
				UK_line.AngleFormatCode = -1
				UK_line.ValidPosition = -1	# Need format code to parse angle

			AngleString = line[34:50]

			# Cols  51-54: Position Accuracy.
			try:
				PositionUncertaintyString = line[50:54]
			except ValueError:
				PositionUncertaintyString = None # Not required (initialized state)

			try:
				UK_line.EpochCode = int(line[54])
			except ValueError:
				log.error("Valid Epoch Code not specified.")
				UK_line.EpochCode = 0
				UK_line.ValidPosition = -1	# If EpochCode as read is Zero (not blank), we won't trigger this error

			if ( (UK_line.AngleFormatCode>=0) and (UK_line.EpochCode>=0) and AngleString ):
				try:
					angle = Angle(
						UK_line.AngleFormatCode, 
						UK_line.EpochCode, 
						AngleString, 
						0, 
						PositionUncertaintyString,
						"UK"
						)
					if (angle.ValidAngle):
						UK_line.AZ = angle.AZ
						UK_line.EL = angle.EL
						UK_line.RA = angle.RA
						UK_line.DEC = angle.DEC
						UK_line.Epoch = angle.Epoch
						UK_line.PositionUncertainty = angle.Uncertainty
					else:
						UK_line.ValidPosition = -1
				except:
					log.error("Problem angle: '{}'".format(AngleString))
					UK_line.ValidPosition = -1
			else:
				# Note - Not an error, because the standards provide for "reporting for duty"
				# even if no observations are made.
				# FIXME: Although this does not exist in the UK format (?
				log.warning("No Position Data")
				UK_line.ValidPosition = -1

			# This is the brightest stellar magnitude attained by the satellite during
			# the period of one minute centred on the time of the observation. It is
			# entered as a 3 digit number, in one of the following forms:

			# (1) if the satellite is brighter than magnitude +9.9
			# 	Column 69 is entered + or -
			# 	Columns 70 and 71 state the numerical value, formatted as Mm

			# (2) if the satellite is fainter than magnitude +9.9
			# 	The sign is omitted, and the numerical format is MMm

			if (line_length>=71):
				# Cols  69-71: Brightest Visual Magnitude
				vmag = line[68:71]
				if (vmag[0]=="+" or vmag[0]=="-"):
					UK_line.VisualMagnitude_high = float(vmag[0] + vmag[1] + '.' + vmag[2])
				else:
					UK_line.VisualMagnitude_high = float(vmag[0:2] + '.' + vmag[2])


			if (line_length>=74):
				# Cols  72-74: Faintest Visual Magnitude
				vmag = line[71:74]
				if (not vmag):
					pass
				elif (vmag == "INV"):
					UK_line.VisualMagnitude_low	= 99
				elif (vmag[0]=="+" or vmag[0]=="-"):
					UK_line.VisualMagnitude_low = float(vmag[0] + vmag[1] + '.' + vmag[2])
				elif(vmag):
					UK_line.VisualMagnitude_low = float(vmag[0:2] + '.' + vmag[2])

				if (UK_line.VisualMagnitude_high and not UK_line.VisualMagnitude_low):
					UK_line.VisualMagnitude = UK_line.VisualMagnitude_high	

			if (line_length>=79):
				# Cols  75-79: Flash Period
				# Time in seconds between successive maxima, formatted as SSSss
				FlashPeriodString = line[74:79]
				try:
					UK_line.FlashPeriod = float(FlashPeriodString[0:3] + '.' + FlashPeriodString[3:5])
				except ValueError:
					UK_line.FlashPeriod = None 

			if (line_length>=80):
			# UK format calls it "remarks" but it correctly maps to IOD-style Optical Behavior codes
				UK_line.OpticalCode = line[79]

			# Grab any trailing comments after the last valid field
			if (line_length>=81):
				UK_line.Remarks = line[80:].strip()

			uk_count += 1
			IODentryList.append(UK_line)
	if (len(IODentryList)):
		return IODentryList
	else:
#		log.warning("Error: No data found. Call parse_uk_lines with a properly formatted UK line comfirmed by format_test_uk")
		return(False)


# Parse enter blocks of text instead of just lines
# This is the typical use-case, as IOD entries usually come in clusters
# We are assuming there is only one full RDE record (with multiple observations) per message block
def parse_rde_record(block=False):
	""" Parse a block of text containing RDE-formatted observation records, and return formatted records + count.

	Format documented at: http://www.satobs.org/position/RDEformat.html
	"""
	match = rde_format_re.search(block) # FIXME: Using older version per https://github.com/consensys-space/trusat-orbit/issues/5
	IODentryList = []
	rde_count = 0
	earliest_time = datetime(1957,10,4,19,28,34,0)
	import_time = datetime.utcnow()
	if (match):
		lines = match.group(0).split('\n')
		
		# Start with valid position until demoted
		ValidPosition = 1

		# Note that RDE Format does not specify Object Numbers
		# We're setting to 0, and relying on the SQL Trigger 'add_object_number' to update it on INSERT 
		# From the InternationalDesignation
		ObjectNumber = 0

		# Match automatically returns the first line
		# Pop off the list and process the header which applies to all following records
		line = lines.pop(0)

		# 'Observing Site Number'
		Station = line[0:4]
		YYMM = line[5:9]

		try:
			TimeUncertainty = float(line[10:13])
		except ValueError:
			TimeUncertainty = 0

		try:
			TimeStandardCode = int(line[13])
		except ValueError:
			TimeStandardCode = 0

		# Col 15: Position Format Code.
		# Russell always uses Code 1, with whole arc seconds
		try:
			AngleFormatCode = int(line[14])
		except ValueError:
			log.error("Valid Angle Format Code not specified.")
			AngleFormatCode = 0
			ValidPosition = -1	# Need format code to parse angle

		# Cols  17-19: Position Accuracy.
		# Russell reports 120 arc seconds, as SSS
		try:
			PositionUncertaintyString = right_zero_pad(line[16:19],length=3)
			PositionUncertainty = float(PositionUncertaintyString.lstrip())

		except ValueError:
			PositionUncertainty = 0

		# Col 20: Epoch of star chart used to determine position
		# Russell always uses code 4.
		try:
			EpochCode = int(line[19])
		except ValueError:
			log.error("Valid Epoch Code not specified.")
			EpochCode = 0
			ValidPosition = -1	# If EpochCode as read is Zero (not blank), we won't trigger this error

		line = lines.pop(0)
		YYMMDD = YYMM + line[0:2]

		for line in lines:
			line = line.rstrip()
			if( len(line)>2 and line != "999"):
				RDE_line = IOD(line)

				# TODO: We'll need to do something else with this line.  Maybe UK-format it?
				RDE_line.IODType = "RDE"

				# Grab values from common lines
				RDE_line.Station = Station
				RDE_line.TimeUncertainty = TimeUncertainty
				RDE_line.TimeStandardCode = TimeStandardCode
				RDE_line.AngleFormatCode = AngleFormatCode
				RDE_line.PositionUncertainty = PositionUncertainty
				RDE_line.EpochCode = EpochCode
				RDE_line.ValidPosition = ValidPosition
				RDE_line.ObjectNumber = ObjectNumber

				LaunchYearTwoDigit = int(line[0:2])
				if (LaunchYearTwoDigit < 57):		# Year of Sputnik launch, first cataloged object
					RDE_line.LaunchYear = 2000 + LaunchYearTwoDigit
				else:
					RDE_line.LaunchYear = 1900 + LaunchYearTwoDigit

				piece_letters = launch_piece_number_to_letter(line[5:7])
				RDE_line.InternationalDesignation = str(RDE_line.LaunchYear) + '-' + line[2:5] + piece_letters

				# Doesn't exist in RDE format
				# RDE_line.StationStatusCode = line[21]  # TODO: - Expand out the short and long description of the codes

				RDE_line.DateTimeString = YYMMDD + line[8:17]
				RDE_line.DateTime = DateTime_frompacked(RDE_line.DateTimeString,"RDE")
				# Flag records from the future and pre history
				if not (earliest_time < RDE_line.DateTime < import_time):
					RDE_line.ValidPosition = -1

				AngleString = line[18:31]

				if ( (RDE_line.AngleFormatCode>=0) and (RDE_line.EpochCode>=0) and AngleString):
					try:
						angle = Angle(
							RDE_line.AngleFormatCode, 
							RDE_line.EpochCode, 
							AngleString, 
							RDE_line.PositionUncertainty, 
							"",
							"RDE"
							)
						if (angle.ValidAngle):
							RDE_line.AZ = angle.AZ
							RDE_line.EL = angle.EL
							RDE_line.RA = angle.RA
							RDE_line.DEC = angle.DEC
							RDE_line.Epoch = angle.Epoch
							RDE_line.PositionUncertainty = angle.Uncertainty
						else:
							RDE_line.ValidPosition = -1
					except:
						log.error("Problem angle: '{}'".format(AngleString))
						RDE_line.ValidPosition = -1
				else:
					# Note - Not an error, because the standards provide for "reporting for duty"
					# even if no observations are made.
					# FIXME: Although this does not exist in the RDE format (?
					log.warning("No Position Data")
					RDE_line.ValidPosition = -1

				# This is the brightest stellar magnitude attained by the satellite during
				# the period of one minute centred on the time of the observation. It is
				# entered as a 3 digit number, in one of the following forms:

				# (1) if the satellite is brighter than magnitude +9.9
				# 	Column 69 is entered + or -
				# 	Columns 70 and 71 state the numerical value, formatted as Mm

				# (2) if the satellite is fainter than magnitude +9.9
				# 	The sign is omitted, and the numerical format is MMm

				# Cols  32-35: Brightest Visual Magnitude
				#  This is the brightest stellar magnitude attained by the satellite during
				#  the period of one minute centred on the time of the observation, It is
				#  entered as a 3 digit number.

				#  If magnitude is negative, enter "-" in column 32
				#  Columns 32-35 state the numerical value, formatted as M.m
				vmag = line[31:35]
				RDE_line.VisualMagnitude_high = float(vmag.strip())

				# Cols  36-39: Faintest Visual Magnitude
				#  Format is the same as for Brightest Visual Magnitude.
				#  If the magnitude is constant, repeat the Brightest Visual Magnitude entry
				vmag = line[35:39]
				if (not vmag):
					RDE_line.VisualMagnitude_low	= None
				else:
					RDE_line.VisualMagnitude_low = float(vmag.strip())

				# Cols  40-42: Flash Period
				#  Time in seconds between successive maxima, formatted as SSS.sss
				#  If value is less than 1, show one leading zero.
				#  Show only significant trailing zeros.
				#  Omit decimal point for whole number values.

				FlashPeriodString = line[39:42]
				try:
					RDE_line.FlashPeriod = float(FlashPeriodString.strip())
				except ValueError:
					RDE_line.FlashPeriod = None 

				# Col      43: Remarks.
				RDE_line.OpticalCode = line[42]
				RDE_line.Remarks = line[43:].strip()
				RDE_line.line = RDE_line.line[0:20] + ' ' + YYMM + line.strip()
				rde_count += 1
				IODentryList.append(RDE_line)
			elif (line != "999"):
				YYMMDD = YYMM + line[0:2]
	if (len(IODentryList)):
		return IODentryList
	else:
#		log.warning("Error: no records found. Call parse_rde_record with a properly formatted RDE line comfirmed by format_test_rde")
		return(False)

""" The following three functions provide a simple interface to:
	1) Open a file
	2) Test for presence of requested record.
	3) Parse the records if they exist, and
	4) Return a list of parsed records, plus the count
"""

def get_iod_records_from_file(filename):
	""" Open a file and return a list with parsed IOD records, and record count """
	with open(filename) as file:
		IOD_file_lines = file.read()

	IOD_records = []
	IOD_records = parse_iod_lines(IOD_file_lines)

	if (len(IOD_records)):
		return IOD_records
	else:
		return False


def get_uk_records_from_file(filename):
	""" Open a file and return a list with parsed UK records, and record count """
	with open(filename) as file:
		UK_file_lines = file.read()

	UK_records = []
	UK_records = parse_uk_lines(UK_file_lines)

	if (len(UK_records)):
		return UK_records
	else:
		return False


def get_rde_records_from_file(filename):
	""" Open a file and return a list with parsed RDE records, and record count """
	with open(filename) as file:
		RDE_file_string = file.read()

	RDE_records = []
	RDE_records = parse_rde_record(RDE_file_string)

	if (len(RDE_records)):
		return RDE_records
	else:
		return False

def get_iod_records(block):
	""" Return a list with parsed IOD records """
	IOD_records = parse_iod_lines(block)

	if (len(IOD_records)):
		return IOD_records
	else:
		return False


def get_uk_records(block):
	""" Return a list with parsed UK records """
	UK_records = parse_uk_lines(block)

	if (len(UK_records)):
		return UK_records
	else:
		return False


def get_rde_records(block):
	""" Return a list with parsed RDE records """	
	RDE_records = parse_rde_record(block)

	if (len(RDE_records)):
		return RDE_records
	else:
		return False


def get_obs_from_text(block):
	""" Return a list of IOD Class variables from text block

	Returns data for the first format detected in the text block.
	For mixed data formats use get_obs_from_mixed_text(block)
	"""
	IOD_records = []
	# Start with most numerous records, move to least
	# Assume there's only one record type per file
	try:
		IOD_records = get_iod_records(block)
	except:
		pass

	if not (len(IOD_records)):
		try: 
			IOD_records = get_rde_records(block)
		except:
			pass

	if not (len(IOD_records)):
		try: 
			IOD_records = get_uk_records(block)
		except:
			pass

	return IOD_records


def get_obs_from_mixed_text(block):
	""" Return a list of IOD Class variables from text block - will detect IOD, UK and RDE formats

	For better performance with homogenous formats use get_obs_from_text(block)
	"""
	IOD_records = []
	UK_records  = []
	RDE_records = []
	try:
		IOD_records = get_iod_records(block)
	except:
		pass

	try: 
		UK_records = get_uk_records(block)
	except:
		pass

	try: 
		RDE_records = get_rde_records(block)
	except:
		pass

	if (UK_records):
		IOD_records.append(UK_records)	
	if (RDE_records):
		IOD_records.append(RDE_records)	

	return IOD_records


def get_IOD_record_counts(IOD_records):
	""" Simple interface to provide stats on a list of IOD Class variables

	Given a list of IOD_records, returns a dictionary of record counts by type
	"""
	TotalObsCount = len(IOD_records)

	TotalCount_IOD = 0
	TotalCount_RDE = 0
	TotalCount_UK  = 0

	for IOD in IOD_records:
		if   (IOD.IODType is "IOD"):
			TotalCount_IOD += 1
		elif (IOD.IODType is "RDE"):
			TotalCount_RDE += 1
		elif (IOD.IODType is "UK"):
			TotalCount_UK  += 1

	IOD_record_counts = {
		"total" : TotalObsCount,
		"IOD"	: TotalCount_IOD,
		"RDE"	: TotalCount_RDE,
		"UK"	: TotalCount_UK
	} 
	return IOD_record_counts


def write_uk_line(IOD):
	""" Writes a UK-formatted line given an IOD class variable as input """
	# TODO: Work in progress.  Currently stalled at converting decimal angle back to the string format.

	launch_year = int(IOD.InternationalDesignation[0:2])
	# launch_year = datetime.strftime(datetime(launch_year,1,1),'%y')

	launch_num = int(IOD.InternationalDesignation[5:8])
	launch_piece = launch_piece_letter_to_number(IOD.InternationalDesignation[8:11])

	obs_time_fmt = "%y%m%d%H%M%S%f"
	obs_time = datetime.strftime(IOD.DateTime,obs_time_fmt)[0:16]

	string_num = "{0:.4f}".format(IOD.TimeUncertainty)
	time_accuracy = string_num[0] + string_num[2:]

	

	uk_string = "{LAUNCHYR:02d}{LAUNCHNUM:03d}{LAUNCHPIECE:02d}{SITE:04d}{OBSTIME:16s}{TIMEACC:5s}".format(
		LAUNCHYR = launch_year,
		LAUNCHNUM = launch_num,
		LAUNCHPIECE = launch_piece,
		SITE = IOD.Station,
		OBSTIME = obs_time,
		TIMEACC = time_accuracy
	)

	print(uk_string)