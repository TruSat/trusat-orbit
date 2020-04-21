#!/usr/bin/env python
"""  Fit a TLE prediction to a reference TLE + new IOD observations

 Suite of utilities based on and extending Scott Campbell's C++ satfit code base
 https://github.com/interplanetarychris/scottcampbell-satfit) 
 for reading visual observations and updating TLEs 
 """

from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values
import sys

# As of 28 July 2019, python3.6 is the default "python3" in apt-get install python3 on Ubuntu
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version 3.6")
    sys.exit(1)

import configparser                 # config file parsing
import argparse                     # command line parsing
import os
from datetime import timedelta, datetime
from time import time                         # For performance timing
from math import (fabs, radians, sin, cos, pi, sqrt, fmod, acos, asin, atan, tan, degrees, modf, floor)    # Fast/precise math functions                      
import numpy as np

import logging
import string
import copy
import re
from getpass import getpass # For getting input without newline

import logging
log = logging.getLogger(__name__)

from spacetrack import SpaceTrackClient

# python SGP4 from git+https://github.com/interplanetarychris/python-sgp4@cython-7-dec-15-vallado
# Until the following pull request is approved
# https://github.com/brandon-rhodes/python-sgp4/pull/35

# Use local/dev version of python-sgp4

from sgp4.api import Satrec, SatrecArray, SGP4_ERRORS, jday
try:
    from caccelerated import *
    # from satfit_caccelerated import *
#    from sgp4.cpropagation import sgp4, sgp4init
#    from sgp4.cmodel import Satellite
except ImportError as e:
    print(e)
    from sgp4.propagation import sgp4, sgp4init
    from sgp4.model import Satellite
from sgp4 import earth_gravity

import trusat.iod

from trusat.tle_util import make_tle, append_tle_file, TLEFile, tle_fmt_epoch, datetime_from_tle_fmt, assumed_decimal_point, checksum_tle_line, TruSatellite, make_tle_from_SGP4_satrec

# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
backend_path = os.path.join(parentdir, "../trusat-backend")
sys.path.insert(1,backend_path) 
import database

from elfind import SGN, so2r

####### DECLARE GLOBAL CONSTANTS #######
TWOPI = 2*pi
NOCON = TWOPI/1440.0
DE2RA = pi/180.0

db = None
startDate = False

# FIXME: Legacy globals which are being worked out of the code
TLE_ref = None  # External TLE seed reference file
srch = "W"      # Search scope
epoch = None    # epoch of elements in tle format
ii = None   # inclination, degrees
om = None   # right ascension of ascending node, degrees
ec = None   # eccentricity
ww = None   # argument of the perigee, degrees
uu = None   # longitude
ma = None   # mean anomaly, degrees
nn = None   # mean motion, revolutions/day
c2 = None   # internal drag term
pgee = None # perigee
agee = None # apogee
bstar = None    # BSTAR drag term
whichconst = None  # Earth gravity constants

odata = None
ll = None
rd = None
iod_line = [] # Observation lines

# TODO: Document the meaning of these flags
# int
#    xi = 0, xe = 0, xw = 0, xn = 0;   // flags
xe = False
xi = False
xn = False
xw = False

# // please note that these variables are in TLE format, i.e. degrees, whereas
# // their counterparts inside SGP4 are in radians.

# double la, lo, hh;          // observer parameters
la = 0
lo = 0
hh = 0
# int first = 1, batch = 0, out = 0;
# int nobs, len;              // number of observations
nobs = 0             # number of observations
file = None             # input file string declared
# char buf[81];               // input buffer
# char name[81], tle1[81], tle2[81];    // original tle data lines
line1 = None  # read and write tle buffers
line2 = None  # read and write tle buffers
# double sum, nsum, xsum;     // sum of squares of residuals
sum = None
nsum = None
xsum = None
# double xrms;
# double zero, tleh;          // history parameters
# FILE *fp, *fpo;

# TODO: Define what these variables are for
estep = None
wstep = None
nstep = None

## Variable / Function reference from original Scott Campbell C++ code
# class Satellite(object):
#     """
#     /*************************** Class Satellite **********************************
#     *                                                                             *
#     *                            Two Constructors                                 *
#     *                                                                             *
#     *                                                                             *
#     *  Satellite sat(tle_date,                                                    *
#     *                inclination_degrees,                                         *
#     *                ascending_node_degrees,                                      *
#     *                eccentricity,                                                *
#     *                argument_of_perigee_degrees,                                 *
#     *                mean_anomaly_degrees,                                        *
#     *                mean_motion_revolutions_per_day,                             *
#     *                bstar)                                                       *
#     *     Eight argument constructor initializes satellite at epoch tle_date.     *
#     *                                                                             *
#     *                                                                             *
#     *  Satellite sat(tle_date, position_pointer, velocity_pointer, bstar)         *
#     *     Four argument constructor produces mean(SGP4)orbital elements           *
#     *                                                                             *
#     *                                                                             *
#     *                      Seven Member Functions                                 *
#     *                                                                             *
#     *  sat.delta_t(time);                                                         *
#     *                Updates satellite elements to new position at time.  Time    *
#     *                is either Julian day or TLE date. Program decides which      *
#     *                format based on relative magnitude.                          *
#     *                                                                             *
#     *  sat.delta_el(epoch_Julian date,                                            *
#     *               inclination_degrees,                                          *
#     *               ascending_node_degrees,                                       *
#     *               eccentricity,                                                 *
#     *               argument_of_perigee_degrees,                                  *
#     *               mean_anomaly_degrees,                                         *
#     *               mean_motion_degrees_per_day,                                  *
#     *               bstar)                                                        *
#     *         Change elements for an existing satellite object.                   *
#     *                                                                             *
#     *  sat.sxp4(tsince);   Analytic Orbit Propagator                              *
#     *                                                                             *
#     *  sat.sgp4(tsince);   Low Orbit Propagator                                   *
#     *                                                                             *
#     *  sat.sdp4(tsince);   Deep Orbit Propagator                                  *
#     *                                                                             *
#     *  sat.rv2el(rr, vv);  State vectors to mean elements                         *
#     *                                                                             *
#     *  print_el(sat);     Output mean elements to screen                         *
#     *                                                                             *
#     *  sat.print_rv();     Output state vectors to screen                         *
#     *                                                                             *
#     *                                                                             *
#     *                                                                             *
#     *                       Thirteen Output Values                                *
#     *                                                                             *
#     *  sat.rr;  Topocentric position vector to satellite at time, earth radii.    *
#     *                                                                             *
#     *  sat.vv;  Satellite velocity vector at time, earth radii / minute.          *
#     *                                                                             *
#     *  sat.jd;  Julian date epoch of satellite elements.                          *
#     *                                                                             *
#     *  sat.thetag;  Greenwich Hour Angle, degrees.                                *
#     *                                                                             *
#     *  sat.mjd;     Modified Julian Date.                                         *
#     *                                                                             *
#     *  sat.xincl;   Mean inclination, radians.                                    *
#     *                                                                             *
#     *  sat.xnodeo;  Mean longitude of the ascending node, radians.                *
#     *                                                                             *
#     *  sat.eo;      Mean eccentricity.                                            *
#     *                                                                             *
#     *  sat.omegao;  Mean argument of perigee, radians.                            *
#     *                                                                             *
#     *  sat.xmo;     Mean anomaly, radians.                                        *
#     *                                                                             *
#     *  sat.xno;     Mean motion, revolutions per day.                             *
#     *                                                                             *
#     *  sat.aodp;    Semi major axis, earth radii.                                 *
#     *                                                                             *
#     *  sat.c2;      internal drag term for ndot calculation                       *
#     *                                                                             *
#     *  sat.bstar;   Drag term.                                                    *
#     *                                                                             *                                                         *                                                                            **                                                                             **                                                                             *                                                                             *
#     *******************************************************************************/
#     """
#     def __init__(self, tle, ii, om, ec, ww, ma, nn, bstar):
#         self.tle   = tle
#         self.ii    = ii
#         self.om    = om
#         self.ec    = ec
#         self.ww    = ww
#         self.ma    = ma
#         self.nn    = nn
#         self.bstar = bstar

#     def sxp4(self, tsince):
#         # /* Period > 225 minutes is deep sIe */
#         # python-SGP4 determines this automatically, however
#         a1 = pow(xke / xno, 2/3.0)
#         r1 = cos(xincl)
#         temp = ck2 * 1.5 * (r1*r1 * 3.0 - 1.0) * pow( 1.0 - eo*eo, -1.5)
#         del1 = temp / (a1*a1)
#         ao = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (del1 * 1.654320987654321 + 1.0)))
#         delo = temp / (ao*ao)
#         xnodp = xno / (delo + 1.0)
        
#         # Select a deep-space/near-earth ephemeris
#         # If the object makes less than 6.4 revolutions around the earth...
#         if (TWOPI / (xnodp * 1440.0) >= (1.0 / 6.4)):
#             sdp4(tsince)    # yes,  it should be a deep-space (SDP4) ephemeris
#         else:
#             sgp4(tsince)

class SatrecMeta(object):
    def __init__(self, satnum=0):
        self.satnum = satnum

class Date(object):
    """ From date.h by Scott Campbell campbel7@the-i.net
    /****************************** Class Date ************************************
    *                                                                             *
    *                         Three Object Constructors                           *
    *                                                                             *
    *  Date t1;  Creates a date object, t1,  initialized to computer time at the  *
    *            instant of the creation of the object.                           *
    *                                                                             *
    *  Date t1(time);   Creates a date object, t1, initialized to the time,       *
    *                   which can be in either Julian date or TLE format.         *
    *                                                                             *
    *  Date t1(year, month, day, hour, min, sec);                                 *
    *            Creates a date object, t1, initialized to the calendar date and  *
    *            time passed by the six calendar variables                        *
    *                                                                             *
    *                                                                             *
    *                         Three Member Functions                              *
    *                                                                             *
    *  t1.now();   Re-initializes an existing date object, t1, to the computer    *
    *              time at the instant this command is executed by the computer.  *
    *                                                                             *
    *  t1.input();   Screen input of calendar date variables with current values  *
    *                presented as default. Just press ENTER to accept current     *
    *                value.  Useful to examine current values at run time.        *
    *                                                                             *
    *  t1.print();   Screen output of calendar date variables.                    *
    *                                                                             *
    *                                                                             *
    *                         Ten Output Values                                   *
    *                                                                             *
    *  t1.thetag;  Sidereal time                                                  *
    *  t1.jd;      Julian date                                                    *
    *  t1.mjd;     Modified Julian Date                                           *
    *  t1.tle;     Date in TLE format, only for years  2000 < year < 2100         *
    *  t1.doy;     day-of-year                                                    *
    *  t1.yy;      year                                                           *
    *  t1.mm;      month                                                          *
    *  t1.dd;      day of month, Greenwich                                        *
    *  t1.hr;      hour, Greenwich                                                *
    *  t1.mn;      minute                                                         *
    *  t1.ss;      seconds                                                        *
    *                                                                             *
    *  all output values are immediately available after any initialization       *
    *******************************************************************************/
    """
    # TODO: Make this also deal with an input of a Python datetime format
    def __init__(self, time=None, year=None, month=None, day=None, hour=None, min=None, sec=None, jd=None, jdF=None):
        self.time  = time               # Python datetime
        self.yy    = self.year  = year  # year
        self.mm    = self.month = month # month
        self.dd    = self.day   = day   # day of month, Greenwich
        self.hr    = self.hour  = hour  # hour, Greenwich
        self.mn    = self.min   = min   # minute
        self.ss    = self.sec   = sec   # seconds

        # Ten Output Values - above 5 and the below:
        self.thetag	= None  # Sidereal time in degrees
        self.jd		= jd    # Julian date (integer day)
        self.jdF	= jdF   # Julian date (day fraction)
        self.mjd	= None  # Modified Julian Date
        self.tle	= None  # Date in TLE format, only for years 2000 < year < 2100
        self.doy	= None  # day-of-year

        if (self.time is None and self.year is None):
            """ Creates a date object, t1,  initialized to (UTC) computer time at the 
            instant of the creation of the object."""
            self.time = datetime.utcnow()
            self.timevars_from_datetime()
        elif (self.time):
            """ Creates a date object, t1, initialized to the time,
            which can be in either Julian date or TLE format."""
            if (self.time < 2400000):    # this date is in tle format
                self.tle = self.time
                self.time = datetime_from_tle_fmt(self.tle)
                self.timevars_from_datetime()
                (self.jd, self.jdF) = jday(self.yy, self.mm, self.dd, self.hr, self.mm, self.ss)
                self.sidereal()
            else:                   # this date is julian
                # self.jd = self.time
                self.calcmjd()
                self.sidereal()
                self.time = jday_to_datetime(self.jd, self.jdF)
        elif (self.yy):
            """ Creates a date object, t1, initialized to the calendar date and
            time passed by the six calendar variables """
            # TODO: Deal with other variables potentially being "None"
            (_, subsec) = divmod(self.ss,1)
            subsec = int(subsec*1E6)
            intss = int(self.ss)
            self.time = datetime(self.yy,self.mm,self.dd,self.hr,self.mn,intss,subsec)
        else:
            # Shouldn't ever get here, default is to create the current time
            pass

        if (not self.time):
            self.time = datetime(self.yy,self.mm,self.dd,self.hr,self.mn,self.ss)

        # Fill out rest of internal variables
        if (not self.jd):
            (self.jd, self.jdF) = jday(self.yy, self.mm, self.dd, self.hr, self.mn, self.ss)
            self.calcmjd()

        if (not self.doy):
            self.doy = self.time.timetuple().tm_yday

        if (not self.tle):
            self.tle = tle_fmt_epoch(self.time)

        if (not self.thetag):
            self.sidereal()

    def now(self):
        """ Re-initializes an existing date object, t1, to the computer
        time at the instant this command is executed by the computer."""
        # TODO - only if needed

    def input(self):
        """Screen input of calendar date variables with current values
        presented as default. Just press ENTER to accept current
        value.  Useful to examine current values at run time."""
        # TODO - final polish

    def print(self):
        """Screen output of calendar date variables."""
        print("Year   {:d}".format(self.yy))
        print("Month  {:d}".format(self.mm))
        print("Day    {:d}".format(self.dd))
        print("Hour   {:d}".format(self.hr))
        print("Minute {:d}".format(self.mn))
        print("Second {:d}".format(self.ss))

    def sidereal(self):
        """calculate Greenwich sidereal hour angle"""
        t = ( (self.jd + self.jdF) - 2451545) / 36525
        thetag = 280.46061837 + 360.98564736629 * ((self.jd + self.jdF) - 2451545) \
                    + .000387933 * (t * t) - (t * t * t) / 38710000
        self.thetag = fmod(thetag, 360)     # degrees

    def timevars_from_datetime(self):
        self.yy = self.time.year
        self.mm = self.time.month
        self.dd = self.time.day
        self.hr = self.time.hour
        self.mn = self.time.minute
        self.ss = self.time.second

    def calcmjd(self):
        self.mjd = (self.jd + self.jdF) - 2400000.5


class Locate(object):
    """ An incomplete implementation of the Locate class by scott campbell

    Just enough to give the Locate.rre vector
    """
    def __init__(self, t1, lat, lon, hite):
        self.t1 = t1
        self.lat = lat
        self.lon = lon
        self.hite = hite

        self.hh = None
        self.la = None
        self.lo = None

        self.jd = self.t1.jd
        self.jdF = self.t1.jdF
        self.rre = [0, 0, 0]
    
        self.convert()
        self.rgeo()

    def convert(self):
        """ Convert units to km / rad / earth radii """
        xkmper = 6378.135               # 6378.1363 equatorial earth radius, km
        ##//   hh *= .3048              # feet to meters <- comment out for meters
        self.hh = self.hite * 0.001     # meters to kilometers
        self.hh /= xkmper               # kilometers to earth radii
        self.la = radians(self.lat)     # convert to radians
        self.lo = fmod(self.lon, 360)   # convert w long(-) to e long(+)

    def rgeo(self):
        """ vector to topocentric position """
        # t1 = Date(jd)               # calculate sidereal time
        theta0 = radians(fmod(self.t1.thetag + self.lo, 360))
        ff = 1 / 298.2572
        f = sqrt(1 - ff * (2 - ff) * sin(self.la)*sin(self.la))
        gc = 1 / f + self.hh
        gs = (1 - ff)*(1 - ff) / f + self.hh

        self.rre[0] = gc * cos(theta0) * cos(self.la)
        self.rre[1] = gc * sin(theta0) * cos(self.la)
        self.rre[2] = gs * sin(self.la)

def write_tle(file):
    pass

# // read site data, from file
def read_site(line):
    global la
    global lo
    global hh

    # sscanf(line + 16, "%4d", &site);    // read site number from iod_line
    site = int(line[16:20]) # read site number from iod_line

    # open stations.in file to read site data for each observation
    try:
        with open("data/stations.in", "r") as fp:
            for inp_str in fp:
                # global la, lo, hh
                # sscanf(inp_str, "%d %s %lf %lf %lf",&num, &fl, &la, &lo, &hh);
                arr = inp_str.split()
                num = arr[0]
                lat  = arr[2]
                lon  = arr[3]
                elev  = arr[4]
                try:
                    num = int(num)
                    la = float(lat)
                    lo = float(lon)
                    hh = float(elev)
                    if (num == site):
                        return
                except:
                    pass
            # end for
        print("\nno site data for {:d}".format(site))
        buf = input("[exit]")
        sys.exit(0)
    except:
        print("\nno stations.in file found")
        buf = input("[exit]")
        sys.exit(0)

def read_obssf(IOD_Records, Stations=None):
    """ decodes the iod_line data 
    
    This one is a more direct transcoding of Scott Campbell's, rather than the one in elfind
    which uses skyfield to get the julian data and topocentric position of the station
    """

    # for range in (nobs):
    #     # sscanf(iod_line[i] + 16, "%4d", &obscode);
    #     obscode = int(iod_line[14:18])
    #     year = iod_line
    #     sscanf(iod_line[i] + 23, "%4d %2d %2d %2d %2d %5s",
    #             &year, &month, &day, &hour, &min, &sbuf);
    #     sec = atof(sbuf);
    #     sec /= pow(10, strlen(sbuf) - 2);
    #     sscanf(iod_line[i] + 44, "%1d", &format);
    #     sscanf(iod_line[i] + 45, "%1d", &age);
    #     epoch = age;
    #     sscanf(iod_line[i] + 47, "%2lf %2lf %3lf %1s", &ra, &mm, &ss, &sn);
    #     sscanf(iod_line[i] + 55, "%2lf %2lf %2s", &dc, &dd, &dbuf);
    #     ds = atof(dbuf);
    #     if(strlen(dbuf) == 1) ds *= 10;
    #     sign = (sn[0] == '-') ? -1 : 1 ;

    # switch(format)      // change to HHhhhh, DDdddd
    # {
    #     /* Format 1: RA/DEC = HHMMSSs+DDMMSS */
    #     case 1 : ra += mm / 60 + ss / 36000;
    #             dc  = sign * (dc + dd / 60 + ds / 3600); break;
    #     /* Format 2: RA/DEC = HHMMmmm+DDMMmm */
    #     case 2 : ra += mm / 60 + ss / 60000;
    #             dc  = sign * (dc + dd / 60 + ds / 6000); break;
    #     /* Format 3: RA/DEC = HHMMmmm+DDdddd */
    #     case 3 : ra += mm / 60 + ss / 60000;
    #             dc  = sign * (dc + dd / 100 + ds / 10000); break;
    #     /* Format 7: RA/DEC = HHMMSSs+DDdddd */
    #     case 7 : ra += mm / 60 + ss / 36000;
    #             dc  = sign * (dc + dd / 100 + ds / 10000); break;
    #     default : printf("\nFormat not found\n");
    #             s_in("[exit]", buf);
    #             exit(0);
    # }
    # FIXME get rid of these globals

    # Sites = CosparSite("data/stations.in")
    csi = .0055878713278878
    zet = .0055888307019922
    the = .0048580335354883

    nobs = len(IOD_Records) # Number of iod-compliant formatted lines in the input file

    mtype = 'double,double,double,S4,int32'
    # mtype = np.dtype([('julian date','double'), ('ra radians', 'double'), 
    #                    ('dec radians', 'double'), ('Station', 'S4'),
    #                    ('Database observation id', 'int32')])
    ll    = np.zeros((nobs,3))
    # odata = np.zeros((nobs,5))
    odata = np.zeros((nobs,3))
    rd    = np.zeros((nobs,3))

    obs_meta = []

    i = 0
    for iod_line in IOD_Records:
        # Grab the most recent version of these variables for writing eventual TLE
        # FIXME Scott's original code grabs the 2nd or "middle" datetime for epoch
        ra = radians(iod_line.RA)
        dc = radians(iod_line.DEC)
        epoch_datetime = iod_line.obs_time
        ssn = iod_line.ObjectNumber
        station_num = iod_line.Station

        (year, month, day, hour, minute, second) = epoch_datetime.timetuple()[:6]
        microseconds = int(epoch_datetime.strftime('%f'))
        sec_with_microseconds = second + microseconds/1.0E6

        # self.jdsatepoch = jday(year, month, day, hour, minute, sec_with_microseconds)
        t1 = Date(year=year, month=month, day=day, hour=hour, min=minute, sec=sec_with_microseconds)

        # Skip the skyfield stuff while debugging against satfit.cpp
        # ts = load.timescale()
        # (year, month, day, hour, minute, second) = iod_line.DateTime.timetuple()[:6]
        # t_skyfield = ts.utc(year, month, day, hour, minute, second)
        # t1_jd = t_skyfield.tt

        if (iod_line.Epoch == 4): # precess from B1950 to J2000
            a = cos(dc) * sin(ra + csi)
            b = cos(the) * cos(dc) * cos(ra + csi) \
                - sin(the) * sin(dc)
            c = sin(the) * cos(dc) * cos(ra + csi) \
                + cos(the) * sin(dc)
            ra = atan(a / b) # ra - zet
            if (b < 0):
                ra += pi
                ra += zet # right ascension, radians
                ra += 1 / 30000
                dc = asin(c)
            if (abs(dc) > 1.4):
                dc = c / abs(c) * acos(sqrt(a*a + b*b))

        # precession from J2000
        t = ( (t1.jd+t1.jdF) - 2451545) / 36525  
        csi = radians((2306.2181 + .30188 * t + .017998 *t*t) * t) / 3600
        zet = radians((2306.2181 + 1.09468 * t + .018203 *t*t) * t ) / 3600
        the = radians((2004.3109 - .42665 * t - .041833 *t*t) * t ) / 3600
        a = cos(dc) * sin(ra + csi)
        b = cos(the) * cos(dc) * cos(ra + csi) \
            - sin(the) * sin(dc)
        c = sin(the) * cos(dc) * cos(ra + csi) \
            + cos(the) * sin(dc)
        ra = atan(a / b) # ra - zet
        if (b < 0):
            ra += pi
        ra += zet # right ascension, radians
        dc = asin(c)
        if (abs(dc) > 1.4):
            dc = c / abs(c) * acos(sqrt(a*a + b*b))

        # line-of-sight vectors
        ll[i][0] = cos(dc) * cos(ra)
        ll[i][1] = cos(dc) * sin(ra)
        ll[i][2] = sin(dc)

        odata[i][0] = (t1.jd+t1.jdF) # julian date
        odata[i][1] = ra # ra radians (observed)
        odata[i][2] = dc # dc radians (observed)
        # odata[i][3] = int(iod_line.Station) # station


        # odata[i][3] = iod_line.Station # station
        # odata[i][4] = iod_line.obs_id # Database observation id (ParsedIOD.obs_id). This inconveniently stores as a float, since its going into a numpy array.
        var = (iod_line.Station, iod_line.obs_id)
        obs_meta.append(var)

        # odata[i] = [(t1.jd+t1.jdF), ra, dc, iod_line.Station, iod_line.obs_id]

        # (la, lo, hh) = Sites.topos(iod_line.Station)
        # observer_location = Topos(latitude_degrees = la, longitude_degrees = lo, elevation_m = hh)

        # topocentric = observer_location.at(t_skyfield)
        # rd[i] = topocentric.position.km/_XKMPER   # elfind works in units of Earth radii
        if (Stations):
            la = Stations[station_num].latitude
            lo = Stations[station_num].longitude
            hh = Stations[station_num].elevation_m
        else:
            read_site(iod_line.iod_string)  # FIXME: This does a file open/scan for every call

        rre = Locate(t1, la, lo, hh)
        # Locate rre(t1.jd, la, lo, hh)
        rd[i] = rre.rre

        i+=1
    # end for
    return odata, ll, rd, obs_meta, t1


def print_satrec_coe(satrec):
    """
    Debug function for find the right satrec elements to access

    *    satn        - satellite number
    *    bstar       - sgp4 type drag coefficient              kg/m2er
    *    ecco        - eccentricity
    *    epoch       - epoch time in days from jan 0, 1950. 0 hr
    *    argpo       - argument of perigee (output if ds)
    *    inclo       - inclination
    *    mo          - mean anomaly (output if ds)
    *    no          - mean motion
    *    nodeo       - right ascension of ascending node

    """

    print("""
Epoch:             {EPOCH}
Semi-Major Axis:   {SEMI}
Eccentricity:      {ECC}
Inclination:       {INC}
RAAN (omega)       {OMEGA}
Arg P:             {ARGP}
True Anomaly (NU): {NU}
Mean Anomaly (mo): {MO}
Mean Motion  (no): {NO_KOZAI}
""".format(          
        EPOCH=satrec.jdsatepoch,
        SEMI=satrec.a,
        ECC=satrec.ecco,
        INC=degrees(satrec.inclo),
        OMEGA=degrees(satrec.nodeo),
        ARGP=degrees(satrec.argpo),
        NU="?",
        MO=degrees(satrec.mo),
        NO_KOZAI=satrec.no_kozai/NOCON))

    mean_elements="""
    # single averaged mean elements
    satrec.am {AM}  # Semi Major Axis
    satrec.em {EM}  # Eccentricity
    satrec.im {IM}  # Inclination
    satrec.Om {Om}  # RAAN
    satrec.mm {MM}  # Argument of Perigee
    satrec.nm {NM}  # Mean motion""".format(
                    AM=satrec.am,  
                    EM=satrec.em, 
                    IM=satrec.im, 
                    Om=satrec.Om, 
                    MM=satrec.mm,  
                    NM=satrec.nm)
    

def print_el(sat, satmeta, deg=False, quiet=False):
    """ print TLE to screen """

    newTLE = make_tle_from_SGP4_satrec(sat,satmeta,classification="T")

    # tle_epoch = tle_fmt_epoch(sat.epoch)

    # if (deg==False):
    #     xincl  = degrees(sat.inclo)
    #     xnodeo = degrees(sat.nodeo)
    #     omegao = degrees(sat.argpo)
    #     xmo    = degrees(sat.mo)
    #     xno    = sat.no_kozai/NOCON

    # line1 = "1 {:5d}U {:<8s} {:14s} 0.00000073  00000-0  50000-4 0    00".format(sat.satnum,sat.intldesg,tle_epoch)
    # # TODO: Deal with First Derivative xno, Second derivative xno, Bstar
    # line2 = "2 {:05d} {:8.4f} {:8.4f} {:7s} {:8.4f} {:8.4f} {:11.8f}    00".format(
    #         sat.satnum, xincl, xnodeo, eo_string, omegao, xmo, xno)

    # line1 = line1[:68] + str(checksum_tle_line(line1))
    # line2 = line2[:68] + str(checksum_tle_line(line2))

    if not quiet:
        print()
        print("{:s}".format(newTLE.line0))
        print("{:s}".format(newTLE.line1))
        print("{:s}".format(newTLE.line2))
        print("  {:.1f}km x {:.1f}km x {:.1f} degrees   {:.1f} min period".format(newTLE.apogee,newTLE.perigee,newTLE.inclination_degrees,newTLE.period/60))
    return(newTLE.line0, newTLE.line1, newTLE.line2)


def sort(iod_line, odata, obs_meta, ll, rd):
    """ Sort odata, ll, rd and iod_line by julian date in odata[i][0]

    Inputs:
        odata       Observation data 
        ll          Line of site vectors
        rd          Topocentric vectors
        iod_line

    Returns,
        Time-sorted versions of the same (same sort across all 4 variables)
    """
    nobs = len(odata)

    # FIXME: Accept this CPP code for now, but make this more pythonic.
    # Not needed to DB queries (which can be returned already sorted)
    unsorted = 1
    while (unsorted):
        for i in range(nobs-1):
            if (odata[i][0] > odata[i + 1][0]):    # switch low/high
                k = rd[i]
                rd[i] = rd[i + 1]
                rd[i + 1] = k

                k = ll[i]
                ll[i] = ll[i + 1]
                ll[i + 1] = k

                m = odata[i]
                odata[i] = odata[i + 1]
                odata[i + 1] = m

                buf = iod_line[i]
                iod_line[i] = iod_line[i + 1]
                iod_line[i + 1] = buf
                unsorted+=1
            # remove duplicates
            # Duplicates if the jd and ra positions are identical between rows
            if ( (odata[i][0] == odata[i + 1][0]) and (odata[i][1] == odata[i + 1][1]) ):
                for s in range (i, nobs-1):
                    rd[s] = rd[s + 1]
                    ll[s] = ll[s + 1]
                    odata[s] = odata[s + 1]
                    iod_line[s] = iod_line[s + 1]
                nobs-=1
        # This is replacement for GOTO logic in the original CPP
        if (unsorted == 1):
            unsorted = 0 # Didn't need to sort anything, clear the flag.
        else:
            unsorted = 1 # Set the flag to sort again.
    return [iod_line, odata, obs_meta, ll, rd]


# Version of print_fit intended to be non-interactive and store fit to variables
# New TruSat development in this version, to preserve original functionality of print_fit
# TODO: move to C-accelerated module
def calc_fit(sat, rd, ll, odata, obs_meta, last_rms, TLE_process):
    nobs = len(odata)
    if (nobs == 0):
        log.error("\nno obs")
        return False

    nrr = 0
    nvv = 0
    Perr = 0
    delt = 0
    xtrk = 0
    az = 0
    el = 0
    asp = 0
    alpha = 0
    sum = 0
    tempv = [0, 0, 1]
    temp =  [0, 0, 1]
    rr =    [0, 0, 1]
    nv =    [0, 0, 1]
    delr =  [0, 0, 1]
    zz =    [0, 0, 1]

    rr    = np.zeros(3)
    vv    = np.zeros(3)
    sat_rr    = np.zeros(3)
    sat_vv    = np.zeros(3)

    weight = None # TODO: Make this mean something

    for j in range (nobs):
        # advance satellite position
        (sat_rr, sat_vv) = delta_t(sat, odata[j][0])

        nrr = norm(sat_rr)
        nvv = norm(sat_vv)

        # computing elevation  # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        el = degrees( acose(np.dot(rd[j], ll[j]) / norm(rd[j])) )
        el = 90 - el
        el = round(el, 1)

        # computing aspect
        asp = degrees( acose(np.dot(ll[j], sat_vv) / nvv) )
        asp = 180 - asp
        asp = round(asp, 1)

        # computing azimuth # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        tempv = np.cross(rd[j], zz)
        nv = np.cross(tempv, rd[j])
        tempv = np.cross(rd[j], ll[j])
        temp = np.cross(tempv, rd[j])
        az = acose(np.dot(nv, temp) / (norm(nv)*norm(temp)))
        tempv = np.cross(temp, nv)
        if (np.dot(tempv, rd[j]) < 0):
            az = 2*pi - az
        if (norm(temp) == 0):
            az = 0.0
        az = degrees(az)
        az = round(az, 1)

        # observed satellite geocentric position vector, rr
        #? This is the one we want to anonymize the observation location
        rr = so2r(nrr, rd[j], ll[j])

        # geocentric position error angle in radians
        Perr = acose(np.dot(sat_rr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, delr, in e.r.
        delr = sat_rr - rr
        temp = np.cross(sat_rr, sat_vv)  # xtrk reference vector points left of track
        sign = SGN(np.dot(delr, temp))

        # observer velocity vector
        tempv = np.cross(zz, rd[j])
        temp = .004351409367 * tempv 
        # observed satellite velocity
        tempv = sat_vv - temp
        nvv = norm(tempv)

        # angle between delr vector and tempv vector, radians
        alpha = acose(np.dot(tempv, delr) / (nvv * norm(delr)))

        # magnitude of delr in direction of tempv, radians
        delt = atan(cos(alpha) * tan(Perr))   # geocentric range error

        # time error
        delt *= nrr / nvv                     # delta r in min
        delt *= 60                            # seconds
        delt  = round(delt, 2)

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # new use of delr variable, predicted line of sight vector
        delr = sat_rr - rd[j]
        nrr = norm(delr)

        # convert to unit vector
        delr = delr/nrr

        # topocentric position error angle in radians
        #? Transforming this to Geocentric would result in the error being reduced by
        #? the plane-projected angle
        #? Could perhaps account for this by storing the magnitude of that angle (and sign?)
        #? without exposing the precise location
        Perr = acose(np.dot(delr, ll[j]))

        # cross track error, as component of topocentric Perr
        #? This two would have an element of projected angle reduction
        #? Could probably correct for it?
        xtrk  = asin(sin(alpha) * sin(Perr))  # cross track magnitude, radians
        xtrk  = degrees(xtrk)                 # degrees
        xtrk *= sign                          # left of track is positive
        xtrk  = round(xtrk, 2)

        # sum position error in squared degrees
        Perr = degrees(Perr)
        sum += Perr*Perr

        (station_id, obs_id) = obs_meta[j]

        TLE_process.update(
            { obs_id : 
                {
                    "aspect"          : asp,
                    "cross_track_err" : xtrk,
                    "time_err"        : delt,
                    "position_err"    : Perr,
                    "obs_weight"      : weight                   
                }
            }
        )

    rms = sqrt(sum / nobs)
    return rms, TLE_process


def print_fit(satx, satxmeta, rd, ll, odata, obs_meta, last_rms):
    nrr = 0
    nvv = 0
    Perr = 0
    delt = 0
    xtrk = 0
    az = 0
    el = 0
    asp = 0
    alpha = 0
    sum = 0
    tempv = [0, 0, 1]
    temp =  [0, 0, 1]
    rr =    [0, 0, 1]
    nv =    [0, 0, 1]
    delr =  [0, 0, 1]
    zz =    [0, 0, 1]
    out = False

    rr    = np.zeros(3)
    vv    = np.zeros(3)
    satx_rr    = np.zeros(3)
    satx_vv    = np.zeros(3)

    nobs = len(odata)
    if (nobs == 0):
        print("\nno obs")
        return

    # copy sat
    # satx = copy.deepcopy(sat)

    fit_string = "\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr   Delta-Epoch"

    print(fit_string)

    # FIXME this (global?) variable
    if (out):
        with open(file, "a") as fp:
            fp.write("\n")
            fp.write(fit_string)

    for j in range (nobs):
        # advance satellite position
        (satx_rr, satx_vv) = delta_t(satx, odata[j][0])
        # jd = floor(odata[j][0])
        # fr = odata[j][0] - jd
        # (error, satx_rr_temp, satx_vv_temp) = satx.sgp4(jd,fr)

        # satx_rr = np.asarray(satx_rr_temp)
        # satx_vv = np.asarray(satx_vv_temp)

        # satx_rr = smult_rtn(1.0 / satx.radiusearthkm, satx_rr) # In Earth radii
        # satx_vv = smult_rtn(1.0 / (satx.radiusearthkm / 60.0), satx_vv)  # In Earth radii / min - seriously!

        # print(f"satx_rr: {satx_rr[0]} {satx_rr[1]} {satx_rr[2]}")
        # print(f"satx.t {satx.t}  epoch: {satx.jdsatepoch + satx.jdsatepochF}")

        nrr = norm(satx_rr)
        nvv = norm(satx_vv)

        # computing elevation  # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        el = degrees( acose(np.dot(rd[j], ll[j]) / norm(rd[j])) )
        el = 90 - el
        el = round(el, 1)

        # computing aspect
        asp = degrees( acose(np.dot(ll[j], satx_vv) / nvv) )
        asp = 180 - asp
        asp = round(asp, 1)

        # computing azimuth # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        tempv = np.cross(rd[j], zz)
        nv = np.cross(tempv, rd[j])
        tempv = np.cross(rd[j], ll[j])
        temp = np.cross(tempv, rd[j])
        az = acose(np.dot(nv, temp) / (norm(nv)*norm(temp)))
        tempv = np.cross(temp, nv)
        if (np.dot(tempv, rd[j]) < 0):
            az = 2*pi - az
        if (norm(temp) == 0):
            az = 0.0
        az = degrees(az)
        az = round(az, 1)

        # observed satellite geocentric position vector, rr
        rr = so2r(nrr, rd[j], ll[j])

        # geocentric position error angle in radians
        Perr = acose(np.dot(satx_rr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, delr, in e.r.
        delr = satx_rr - rr
        temp = np.cross(satx_rr, satx_vv)  # xtrk reference vector points left of track
        sign = SGN(np.dot(delr, temp))

        # observer velocity vector
        tempv = np.cross(zz, rd[j])
        temp = .004351409367 * tempv 
        # observed satellite velocity
        tempv = satx_vv - temp
        nvv = norm(tempv)

        # angle between delr vector and tempv vector, radians
        alpha = acose(np.dot(tempv, delr) / (nvv * norm(delr)))

        # magnitude of delr in direction of tempv, radians
        delt = atan(cos(alpha) * tan(Perr))   # geocentric range error

        # time error
        delt *= nrr / nvv                     # delta r in min
        delt *= 60                            # seconds
        delt  = round(delt, 2)

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # new use of delr variable, predicted line of sight vector
        delr = satx_rr - rd[j]
        nrr = norm(delr)

        # convert to unit vector
        delr = delr/nrr

        # topocentric position error angle in radians
        Perr = acose(np.dot(delr, ll[j]))

        # cross track error, as component of topocentric Perr
        xtrk  = asin(sin(alpha) * sin(Perr))  # cross track magnitude, radians
        xtrk  = degrees(xtrk)                 # degrees
        xtrk *= sign                          # left of track is positive
        xtrk  = round(xtrk, 2)

        # sum position error in squared degrees
        Perr = degrees(Perr)
        sum += Perr*Perr

        # Format time string
        # YYday HHMM:SSsss

        tsince_days = (odata[j][0] - (satx.jdsatepoch + satx.jdsatepochF) ) # time since epoch in days 
        tsince_minutes = tsince_days * 1440.0

        # Reconstruct the DateTime of the observations from the data we have here
        # TODO: Clean up this calculation
        obstime = jday_to_datetime(satx.jdsatepoch,satx.jdsatepochF) + timedelta(minutes=tsince_minutes)
        timestring = obstime.strftime('%y%j %H%M:%S')
        SSS = obstime.strftime('%f')
        SSS = int(1000*(int(SSS)/1E6))

        (station_id, obs_id) = obs_meta[j]

        # fit_string = "({:2d}) {:4s}  {}{:03d}  {:5.1f}  {:5.1f}  {:5.1f}  {:6.2f}   {:6.2f}  {:7.3f}  {:8.5f}".format(
        #     j + 1, odata[j][3].decode("utf-8"), timestring, SSS, az, el, asp, xtrk, delt, Perr, tsince_days)
        fit_string = "({:2d}) {:4s}  {}{:03d}  {:5.1f}  {:5.1f}  {:5.1f}  {:6.2f}   {:6.2f}  {:7.3f}  {:8.5f}".format(
            j + 1, station_id, timestring, SSS, az, el, asp, xtrk, delt, Perr, tsince_days)
        print(fit_string)

        # print fit to file
        # FIXME this (global?) variable
        if (out):
            with open(file, "a") as fp:
                fp.write(fit_string)

    rms = sqrt(sum / nobs)
    delta_rms = rms - last_rms
    print("\nrms{:12.5f} ({:.5f})".format(rms, delta_rms))
    return rms


def print_calc_fit(sat, rd, ll, odata, obs_meta, last_rms, TLE_process):
    fit_string = "\n      STA  YYday HHMM:SSsss   RA     DEC    ASP     XTRK    deltaT   Perr   Delta-Epoch"
    # fit_string = "\n      STA  YYday HHMM:SSsss   RA     DEC     XTRK    deltaT   Perr   Delta-Epoch"
    print(fit_string)
    for j in range(len(odata)):
        # # Format time string
        # # YYday HHMM:SSsss
        tsince_days = (odata[j][0] - (sat.jdsatepoch + sat.jdsatepoch) ) # time since epoch in days 
        tsince = tsince_days * 1440.0 # time since epoch in minutes # TODO: Clean up this calculation
        obstime = jday_to_datetime(sat.jdsatepoch, sat.jdsatepochF) + timedelta(minutes=tsince)
        timestring = obstime.strftime('%y%j %H%M:%S')
        SSS = obstime.strftime('%f')
        SSS = int(1000*(int(SSS)/1E6))

        (station_id, obs_id) = obs_meta[j]

        fit_string = "({:2d}) {:04s}  {}{:03d}  {:5.1f}  {:5.1f}  {:5.1f}  {:6.2f}   {:6.2f}  {:7.3f}  {:8.5f}".format(
            j + 1, station_id, timestring, SSS, 
            degrees(odata[j][1]), 
            degrees(odata[j][2]), 
            TLE_process[obs_id]["aspect"], 
            TLE_process[obs_id]["cross_track_err"], 
            TLE_process[obs_id]["time_err"], 
            TLE_process[obs_id]["position_err"],
            tsince_days)
        print(fit_string)
                 

# ////////////////// FUNCTIONS //////////////////////////////////////////////////


# TODO: move to C-accelerated module
def ww2ma(wx):
    """ find mean anomaly from true longitude and perigee """
    # FIXME uses uu, ec globals
    theta = radians(fmod(uu - wx, 360))
    e = acose((ec + cos(theta)) / (1 + ec * cos(theta)))
    if (theta > pi):
        e = 2 * pi - e
    ma = e - ec * sin(e)
    return degrees(ma)


def align(satx, rd, ll, odata):
    """ sets deltaT error at last obs epoch to zero """
    nobs = len(odata)

    delt = 1    # Give delt a starting value for the loop
    while(fabs(delt) > 1.0e-5): # FIXME global
        # advance satellite position
        (satx_rr, satx_vv) = delta_t(satx, odata[-1][0])
        nrr = norm(satx_rr)
        nvv = norm(satx_vv)

        # observed geocentric position vector, rr
        rr = so2r(nrr, rd[last], ll[last])

        # position error in radians
        Perr = acose(np.dot(satx_rr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, er
        delr = np.subtract(satx_rr, rr)

        # magnitude of delta r in direction of v, radians
        delt = Perr * np.dot(satx_vv, delr) / (nvv * norm(delr))

        # degrees
        delt = degrees(delt)

        # TODO not sure where first gets set
        if (first):
            first = 0
            ma = ma - 0.75*delt
        else:
            ma = ma - delt/5.0

        satx.delta_el(satx.jd, ii, om, ec, ww, ma, nn, bstar)   # FIXME python-SGP4

def fit_out(sat, satmeta, rd, ll, odata, obs_meta, sum):
    out = 1
    sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
    out = 0
    return sum

def elcor(sat, rd, ll, odata):
    # TODO Implement python function call to elcord here
    #    sprintf(file_in, "elcord %s", file);
    # FIXME Fix goto restart
    #    goto restart
    pass

def edit_file():
    # TODO: Figure out if it makes sense to implement this (in VIM, or ENV editor?)
    # sprintf(file_in, "%s", file);
    # //  system(file_in);
    # FIXME Fix goto restart
    # goto restart;
    pass

def edit_history():
    # TODO: Figure out if it makes sense to implement this (in VIM, or ENV editor?)
    # sprintf(file_in, "history%c%05d.txt", '\\', ssn);
    # //  system(file_in);
    # # FIXME Fix goto restart
    pass
    # accept_command()


def viewobs(iod_line):
    # FIXME Make the global variables local
    nobs = len(iod_line)
    print()
    for i in range (nobs):
        print("({:2d}) {}".format(i + 1, iod_line[i]))
    print()
    # accept_command()


def remove(rd, ll, odata, obs_meta, iod_line):
    buf = input("\nRemove Line Number : ")
    
    match = re.search("^([f])?(\d+)",buf)
    if (match):
        j = int(match.group(2).strip())
        if (match.group(1)== "f"):
            rd_del    = np.delete(rd,slice(0,j),0)
            ll_del    = np.delete(ll,slice(0,j),0)    
            odata_del = np.delete(odata,slice(0,j),0)
            obs_meta  = obs_meta[j:]
            iod_line  = iod_line[j:]
        else:
            rd_del    = np.delete(rd,j-1,0)
            ll_del    = np.delete(ll,j-1,0)    
            odata_del = np.delete(odata,j-1,0)
            obs_meta.pop(j-1) 
            iod_line.pop(j-1) 

    return rd_del, ll_del, odata_del, obs_meta, iod_line

def fit(sat, satmeta, rd, ll, odata, obs_meta, sum):
    out = 0
    sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
    print_el(sat, satmeta)
    return sum

def id():
    # requires SATID in same folder
    # TODO Implement python function call to elcord here
    # sprintf(buf, "satid %s", file);
    # system(buf);
    # FIXME Fix goto restart
    # goto restart
    pass

def history(sat, rd, ll, odata):

    while(True): # Forever loop
        print("\n(A)dd  (G)raphs  (E)dit  (Q)uit  ")

        buf = input(": ")
        buf = buf.strip()
        buf = buf.upper()

        ssn = sat.satnum
        file_in = "history/{:05d}.txt".format(ssn)

        try:
            fp =  open(file_in, "r")
        except: #FIXME: Trap appropriate error
            fp.close()
            # create history file
            with open(file_in, "w") as fpo:
                # write current TLE to history
                fpo.write("\n")
                fpo.write("{}".format(line1))
                fpo.write("{}".format(line2))
            print("\nHISTORY created with current TLE")
    
        # GRAPHS
        if (buf == 'G'):
            # TODO: Low priority features
            # i = 1
            # fp = fopen(file_in, "r")

            # file_out = "trendA.csv"

            # with fopen(file_out, "w") as fpo:
            #     fpo.write("time,inc,omega,ecc,perigee,motion,YYddd,,,{:05d}".format(ssn))

            # while(fgets(line1, 80, fp))
            # {
            # if(*line1 == '1' && strlen(line1) == 70)
            # {
            #     fgets(line2, 80, fp);
            #     if(*line2 == '2' && strlen(line2) == 70)
            #     {
            #     // first data line
            #     tleh = atof(line1 + 18);      // epoch in tle format at this point
            #     Date t1(tleh);
            #     tleh = t1.jd;
            #     if(i)
            #     {
            #         zero = tleh;
            #         i = 0;
            #     }
            #     tleh -= zero;

            #     // second data line
            #     ik = atof(line2 + 8);         // inclination, degrees
            #     ok = atof(line2 + 17);        // ascending node, degrees
            #         line2[25] = '.';            // add decimal point to eccentricity
            #     ek = atof(line2 + 25);        // eccentricity
            #         line2[25] = '.';            // add decimal point to eccentricity
            #     wk = atof(line2 + 34);        // perigee, degrees
            #     mk = atof(line2 + 43);        // mean anomaly, degrees
            #     // Make sure mean motion is null-terminated, since rev. no.
            #     // may immediately follow.
            #     line2[63] = '\0';
            #     nk = atof(line2 + 52);        // mean motion, revolutions/day
            #     }
            #     fprintf(fpo, "%f,%f,%f,%f,%f,%f,%.5s\n", tleh, ik, ok, ek, wk, nk, line1 + 18);
            #     theta = tleh;
            # }
            # }
            # Date t1(tle);
            # tleh = t1.jd;
            # tleh -= zero;
            # if(tleh > theta)
            # fprintf(fpo, "%f,%f,%f,%f,%f,%f,%5.0f\n", tleh, ii, om, ec, ww, nn, tle);
            # fclose(fpo);
            # sprintf(file_out, "trendA.xls");
            # // system(file_out);

            # fclose(fp);
            # // sprintf(file_out, "del trendA.csv");
            # // system(file_out);
            pass

        # EDIT
        elif (buf == 'E'):
            # TODO: Low priority features
            # file_in = "history%c%05d.txt".format(ssn)
            # sprintf(file_in, "history%c%05d.txt", '\\', ssn);
            # // system(file_in);
            pass

        # ADD
        elif (buf == 'A'):
            fp.close()
            # write current TLE to history
            write_tle(file_in)
            print("\nCurrent TLE added to HISTORY")

        # QUIT
        elif (buf == 'Q'):
            # FIXME create a local function to pass on TLE command
            # Looks like this is to write the TLE to screen
            # write_tle((char *)"")
            return True

def accept_command(db, sat, satmeta, rd, ll, odata, obs_meta, sum, uu, iod_line):
    """Accept a new command"""

    while(True):    # Forever loop
        print("\nEnter command")
        print(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee   (A)nomaly  (M)otion  (B)star")
        print(" (S)tep   I(D)   (T)ime  (F)it  (V)iew  (R)emove   (W)rite   (Q)uit")

        cmd = input(": ").strip()
        cmd = cmd.upper()

        start_rms = sum

        # Hidden functions
        if (cmd == "G"):    # Graph fit
            sum = fit_out(sat, satmeta, rd, ll, odata, obs_meta, sum)
        elif (cmd == "E"):  # Edit File
            edit_file()
        elif (cmd == "H"):  # History
            history(sat, rd, ll, odata)
        elif (cmd == "Y"):  # Edit History
            edit_history()
        elif (cmd == "C"):  # Elcor
            elcor(sat, rd, ll, odata)
        # elif (cmd == "O"):  # Discover
        #     discover(sat, rd, ll, odata)
        elif (cmd == "O"):  # Object search
            object_search(db,  object=sat.satnum)
        elif (cmd == "U"):  # Maneuver
            maneuver(sat, rd, ll, odata, obs_meta, sum, iod_line)

        # Visible functions
        elif (cmd == "S"):  # Step
            step(sat, rd, ll, odata, sum, uu, "S")
            print_el(sat, satmeta)       # print new elements
        elif (cmd == "L"):  # Loop-Step
            step(sat, rd, ll, odata, sum, uu, "L")
            print_el(sat, satmeta)       # print new elements
        elif (cmd == "Z"):
            step(sat, rd, ll, odata, sum, uu, "Z")
            print_el(sat, satmeta)       # print new elements
        elif (cmd == "I"):  # Incl
            print_el(sat, satmeta)
            incl(sat, satmeta, rd, ll, odata, obs_meta, sum)
        elif (cmd == "N"):  # Node
            print_el(sat, satmeta)
            node(sat, satmeta, rd, ll, odata, obs_meta, sum)
        elif (cmd == "X"):  # Eccentricity
            print_el(sat, satmeta)
            xntrcty(sat, satmeta, rd, ll, odata, obs_meta, sum)
        elif (cmd == "P"):  # Perigee
            print_el(sat, satmeta)
            perigee(sat, satmeta, rd, ll, odata, obs_meta, sum, uu)
        elif (cmd == "A"):  # Mean Anomaly
            print_el(sat, satmeta)
            anomaly(sat, satmeta, rd, ll, odata, obs_meta, sum, uu)
        elif (cmd == "M"):  # Mean Motion
            print_el(sat, satmeta)
            motion(sat, rd, ll, odata, obs_meta, sum)
        elif (cmd == "B"):  # Bstar
            print_el(sat, satmeta)
            bstar_func(sat, satmeta, rd, ll, odata, obs_meta, sum)

        elif (cmd == "D"):  # ID
            id()
        elif (cmd == "T"):  # Time
            time_func(sat, rd, ll, odata, obs_meta, sum)
        elif (cmd == "F"):  # Fit
            sum = fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
        elif (cmd == "V"):  # View observations
            viewobs(iod_line)
        elif (cmd == "R"):  # Remove observations
            (rd, ll, odata, obs_meta, iod_line) = remove(rd, ll, odata, obs_meta, iod_line)
        elif (cmd == "W"):  # Write elements
            success = write_el(db, sat, satmeta, rd, ll, odata, obs_meta, sum, start_rms)
            if (success):
                return True
        elif (cmd == "Q"):  # Quit
            main(db)
        elif (cmd == "."):
            move_epoch_to_previous_perigee(sat)
            print_el(sat, satmeta)
        elif (cmd == ","):
            # print(f"odata[-1][0]: {odata[-1][0]}")
            # (jd, fr) = divmod(odata[-1][0],1)
            # (error, rr, vv) = sat.sgp4(jd, fr) 
            # satrec = Satrec()
            # satrec.sgp4init(sat.satnum, 
            #                 (odata[-1][0]) - 2433281.5, # epoch time in days from jan 0, 1950. 0 hr
            #                 sat.bstar, sat.ndot, sat.nddot, sat.ecco, sat.argpo,
            #                 sat.inclo, sat.mo, sat.no_kozai, sat.nodeo)
            # satrec.jdsatepoch = jd
            # satrec.jdsatepochF = fr
            # sum = fit(satrec, satmeta, rd, ll, odata, sum)

            # print_el(satrec, satmeta)
            move_epoch_to_jd(sat,odata[-1][0])
            print_el(sat, satmeta)

def discover(sat, rd, ll, odata):       # partition search
    global srch   # FIXME global variable
    global uu
    global estep
    global wstep
    global nstep

    print("\n(W)ide  (N)arrow  [{}", srch)

    buf = input(": ").strip().upper()

    ec = sat.ecco
    ww = sat.argpo
    nn = sat.no_kozai/NOCON

    if (buf):
        srch = buf[0]

    # WIDE
    if (buf == 'W'):
        emin  = 0
        emax  = 0.2
        esize = 80
        wmin  = 0
        wmax  = 360
        wsize = 144
        nmin  = 11
        nmax  = 15
        nsize = 80

    # NARROW
    if (buf == 'N'):
        emin = ec * 0.9 # FIXME Access to this global variable
        emax = ec * 1.1 # FIXME Access to this global variable
        esize = 20
        wmin = ww - 2   # FIXME Access to this global variable
        wmax = ww + 2   # FIXME Access to this global variable
        wsize = 20
        nmin = nn - 0.1 # FIXME Access to this global variable
        nmax = nn + 0.1 # FIXME Access to this global variable
        nsize = 20

    if (buf == 'Z'):
        emin = ec - estep # FIXME Access to this global variable
        emax = ec + estep # FIXME Access to this global variable
        esize = 20
        wmin = ww - wstep # FIXME Access to this global variable
        wmax = ww + wstep # FIXME Access to this global variable
        wsize = 20
        nmin = nn - nstep # FIXME Access to this global variable
        nmax = nn + nstep # FIXME Access to this global variable
        nsize = 20

    print("\nemax [{:.7f}".format(emax))
    buf = input("]: ")
    emax = float(buf.strip())

    print("emin [{:.7f}".format(emin))
    buf = input("]: ")
    emin = float(buf.strip())

    estep = (emax - emin) / esize

    print("\nwmax [{:.4f}".format(wmax))
    buf = input("]: ")
    wmax = float(buf.strip())

    print("wmin [{:.4f}".format(wmin))
    buf = input("]: ")
    wmin = float(buf.strip())

    wstep = (wmax - wmin) / wsize

    print("\nnmax [{:.8f}".format(nmax))
    buf = input("]: ")
    nmax = float(buf.strip())

    print("nmin [{:.8f}".format(nmin))

    buf = input("]: ")
    nmin = float(buf.strip())

    nstep = (nmax - nmin) / nsize

    for wk in np.arange(wmin, wmax, wstep):
        theta = radians(uu - wk)
        for ek in np.arange(emin, emax, estep):
            e = acose((ek + cos(theta)) / (1 + ek * cos(theta)))
            if (theta > pi):
                e = 2 * pi - e
            mk = e - ek * sin(e)
            mk = degrees(mk)
            for nk in np.arange(nmin, nmax, nstep):
                # FIXME update class for this
                satx = diff_el(sat,ec=ek,ww=wk,ma=mk,nn=nk)
                # satx(sat.jd, ii, om, ek, wk, mk, nk, bstar)
                # establish the computed ra, dc, at jdo with no perturbations
                rms = find_rms(satx, rd, ll, odata)
                if (rms < sum):
                    sum = rms       # global
                    ec = ek
                    ww = wk
                    ma = mk
                    nn = nk
                # end for nk
            # end for ek
        # end for wk
    ww = fmod(ww, 360)
    ma = fmod(ma, 360)
    diff_el(sat,ec=ec,ww=ww,ma=ma,nn=nn)
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # update mean_anomaly
    anomaly_search(sat, rd, ll, odata, sum)
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # print new elements
    sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
    print_el(sat, satmeta)

    uu = longitude(sat)

    srch = 'Z' # FIXME, this is probably a global
    return sat

def incl(sat, satmeta, rd, ll, odata, obs_meta, sum):
    global srch

    ii = degrees(sat.inclo)
    om = degrees(sat.nodeo)

    while(True): # Forever loop
        print("\n(A)uto  (ii)  (Q)uit  ")
        buf = input(": ").strip()

        try: 
            ii = float(buf)
            xi = 1
            delta_el(sat, ii=ii)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()
            if (buf == 'A'):
                # partition search

                imax = ii + 2
                imin = ii - 2
                omax = om + 2
                omin = om - 2
                xi = 0

                ## This doesn't seem very "auto"
                # TODO: Add in a "seach" options which presents these
                # print("\nimax [{:.4f}".format(imax))
                # buf = input("]: ")
                # imax = float(buf.strip())

                # print("\nimin [{:.4f}".format(imin))
                # buf = input("]: ")
                # imin = float(buf.strip())

                # print("\nomax [{:.4f}".format(omax))
                # buf = input("]: ")
                # omax = float(buf.strip())
                
                # print("\nomin [{:.4f}".format(omin))
                # buf = input("]: ")
                # omin = float(buf.strip())

                node_search(sat, rd, ll, odata, obs_meta, sum, imax, imin, omax, omin)

                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

                srch = 'N' # FIXME: Fix this global
            elif (buf == 'Q'):
                return sat

def node(sat, rd, ll, odata, obs_meta, sum):
    global srch

    while(True): # Froever loop
        om = degrees(sat.nodeo)

        print("\n(A)uto  (om)  (Q)uit  ")
        buf = input(": ").strip()

        try:
            om = float(buf)
            delta_el(sat,om=om)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                satx = sat

                # partition search
                omax = om + 2
                omin = om - 2

                ## This doesn't seem very "auto"
                # print("\nomax [{:.4f}".format(omax))
                # buf = input("]: ")
                # omax = float(buf.strip())
                
                # print("\nomin [{:.4f}".format(omin))
                # buf = input("]: ")
                # omin = float(buf.strip())

                while((omax - omin) > 1e-5):
                    ostep = fabs(rtw(omax, omin) / 20)
                    for ok in np.arange (omin, omax, ostep):
                        satx = delta_el(sat,om=ok)
                        # satx.delta_el(sat.jd, ii, ok, ec, ww, ma, nn, bstar)
                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(satx, rd, ll, odata)
                        if (rms < sum):
                            sum = rms      # global
                            om = ok
                        # end for ok
                    # satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

                    omin = om - ostep
                    omax = om + ostep
                # om = fmod(om, 360)

                delta_el(satx,om=om)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

                srch = 'N' # FIXME: Fix this global variable
            elif (buf == 'Q'):
                return sat


def xntrcty(sat, rd, ll, odata, obs_meta, sum):
    global srch
    global xe

    while(True): # Forever loop
        print("\n(S)earch  (A)uto  (ec)  (Q)uit  ")
        buf = input(": ").strip()

        ec = sat.ecco

        emax = ec * 1.05
        emin = ec * 0.95

        try:
            ec = float(buf)
            xe = 1
            delta_el(sat,ec=ec)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                xe = 0

                ## Commenting this out, as it doesn't seem very "auto"
                # print("\nemax [{:.7f}".format(emax))
                # buf = input("]: ")
                # emax = float(buf.strip())

                # print("\nemin [{:.7f}".format(emin))
                # buf = input("]: ")
                # emin = float(buf.strip())

                while((emax - emin) > 1.e-8):
                    estep = (emax - emin) / 20
                    for ek in np.arange(emin, emax, estep):
                        delta_el(sat,ec=ek)
                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(sat, rd, ll, odata)
                        if (rms < sum):
                            sum = rms
                            ec = ek
                        # end for ek
                    delta_el(sat,ec=ec)
                    emin = ec - estep
                    emax = ec + estep

                delta_el(sat,ec=ec)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

                srch = 'N' # FIXME this global

            # SEARCH
            elif (buf == 'S'):
                print("\nemax [{:.7f}".format(emax))
                buf = input("]: ")
                emax = float(buf.strip())

                print("\nemin [{:.7f}".format(emin))
                buf = input("]: ")
                emin = float(buf.strip())

                print("\nestep [{:.0f}".format(estep))
                buf = input("]: ")
                estep = float(buf.strip())

                estep = (emax - emin) / estep

                ek = ec
                print("\neccentricity  sum")
                for ec in np.arange(emin,emax + estep, estep):
                    delta_el(sat,ec=ec)
                    sum = find_rms(sat, rd, ll, odata)
                    print("\n{:.7f}     {:7.4f}", ec, sum)

                print()
                ec = ek
                delta_el(sat,ec=ec)

            elif (buf == 'Q'):
                return sat

def perigee(sat, rd, ll, odata, obs_meta, sum, uu):
    global srch
    global xw

    while(True): # Forever loop
        pgee = int((sat.altp)*sat.radiusearthkm)
        agee = int( (sat.alta)*sat.radiusearthkm)

        print("\nPerigee = {:.6f} er".format(1+sat.altp), end="")
        print("     {:d} X {:d} km".format(pgee, agee))

        print("\n(S)earch  (A)uto  (ww)  (Q)uit  ")
        buf = input(": ").strip()

        try:
            ww = float(buf)
            uu = longitude(sat)
            ma = ww2ma(ww)
            xw = 1
            sat.delta_el(sat, ww=ww, ma=ma)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()
            # AUTO
            if (buf == 'A'):
                ww = degrees(sat.argpo)
                xw = 0
                wx = ww + 1
                while(fabs(wx - ww) > 1e-4):
                    wx = ww
                    print("\n{:8.4f}  {:.7f}".format(ww, sat.ecco))
                    wmax = ww + 1
                    wmin = ww - 1
                    emax = sat.ecco * 1.01
                    emin = sat.ecco * 0.99
                    perigee_search(sat, rd, ll, odata, obs_meta, sum, uu, wmax, wmin, emax, emin)
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

                srch = 'N' # FIXME Global

            # SEARCH
            elif (buf == 'S'):
                uu = longitude(sat)
                
                ww = degrees(sat.argpo)
                ma = degrees(sat.mo)
                wmax = ww + 2
                wmin = ww - 2
                wstep = 20

                print("\nwmax [{:.4f}".format(wmax))
                buf = input("]: ")
                wmax = float(buf.strip())

                print("\nwmin [{:.4f}".format(wmin))
                buf = input("]: ")
                wmin = float(buf.strip())

                print("\nwstep [{:.0f}".format(wstep))
                buf = input("]: ")
                wstep = float(buf.strip())
                wstep = (wmax - wmin) / wstep

                wk = ww
                mk = ma
                print("\nperigee      sum")
                for ww in np.arange(wmin, wmax + wstep, wstep):
                    ma = ww2ma(ww)
                    sat.delta_el(sat, ww=ww, ma=ma)
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                    sum = find_rms(sat, rd, ll, odata)
                    print("\n{:8.4f}  {:7.4f}".format(fmod(ww, 360), sum))
                print()
                ww = wk
                ma = mk
                sat.delta_el(sat, ww=ww, ma=ma)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

            # QUIT
            elif (buf == 'Q'):
                return sat

    # # FIXME: Looks like this might be a straggler from the beginning of anomaly() below
    # # Not sure how the following code ever gets used
    # # amax and amin are used as temporary variables
    # uu = longitude(degrees(sat.mo), degrees(sat.argpo))
    # amax = radians(uu)
    # amin = sin(degrees(sat.xincl)) * sin(amax)
    # amin *= amin
    # amin = sqrt(1 - amin)
    # amin = degrees((acose(cos(amax) / amin)))
    # if (fmod(uu, 360) > 180):
    #     amin = 360 - amin
    # if(degrees(sat.xincl) < 90):
    #     amax = fmod(degrees(sat.xnodeo) + amin - sat.thetag + 360, 360.0)
    # else:
    #     amax = fmod(degrees(sat.xnodeo) - amin - sat.thetag + 720, 360.0)
    # print("\nE Longitude = {:8.4f}".format(amax))


def anomaly(sat, satmeta, rd, ll, odata, obs_meta, sum, uu):
    global srch

    while(True): # Forever loop
        if ((sat.no_kozai/NOCON) < 1.5):
            # amax and amin are used as temporary variables
            uu = longitude(sat)
            amax = uu * DE2RA
            amin = sin(sat.inclo / DE2RA) * sin(amax)
            amin *= amin
            amin = sqrt(1 - amin)
            amin = ((acose(cos(amax) / amin))) / DE2RA
            if (fmod(uu, 360) > 180):
                amin = 360 - amin
            if((sat.inclo/DE2RA) < 90):
                amax = fmod((sat.nodeo/DE2RA) + amin - satmeta.thetag + 360, 360.0)
            else:
                amax = fmod((sat.nodeo/DE2RA) - amin - satmeta.thetag + 720, 360.0)
            print("\nE Longitude = {:8.4f}".format(amax))

        print("\n(S)earch  (A)uto  (ma)  (L)ast  (Q)uit  ")

        amax = 360
        amin = 0.0
        astep = 20

        buf = input(": ").strip()

        try:
            ma = float(buf)
            uu = longitude(sat)
            delta_el(sat, ma=ma)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                # anomaly_search(sat, rd, ll, odata, sum)
                (jd, fr, rr, vv, err) = get_sgp4_vec_vars(odata)
                anomaly_search(sat, rd, ll, jd, fr, rr, vv, err, sum)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements
                srch = 'N'   # FIXME Global
            # SEARCH
            elif (buf == 'S'):
                print("\namax [{:.7f}".format(amax))
                buf = input("]: ")
                amax = float(buf.strip())

                print("\namin [{:.7f}".format(amin))
                buf = input("]: ")
                amin = float(buf.strip())

                print("\nastep [{:.0f}".format(astep))
                buf = input("]: ")
                astep = float(buf.strip())
                astep = (amax - amin) / astep

                mk = ma
                print("\nanomaly        sum")
                for ma in np.arange(amin, amax + astep, astep):
                    delta_el(sat,ma=ma)
                    sum = find_rms(sat, rd, ll, odata)
                    print("\n{:8.4f}     {:7.4f}".format(ma, sum))
                print()
                ma = mk              # restore
                delta_el(sat,ma=ma)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

            # LAST
            elif (buf == 'L'):
                align(sat, rd, ll, odata)
                uu = longitude(sat)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

            # QUIT
            elif (buf == 'Q'):
                return True


def motion(sat, rd, ll, odata, obs_meta, sum):
    global xn

    while(True): # Forever loop
        print("\n(A)uto  (nn)  (Q)uit  ")
        buf = input(": ").strip()

        try:
            nn = float(buf)
            xn = 1
            delta_el(sat, nn=nn)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                xn = 0
                # update mean motion, no limits
                # motion_search(sat, rd, ll, odata)
                (jd, fr, rr, vv, err) = get_sgp4_vec_vars(odata)
                motion_search(sat, rd, ll, jd, fr, rr, vv, err, sum)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements
            elif (buf == 'Q'):
                return sat


def bstar_func(sat, satmeta, rd, ll, odata, obs_meta, sum):
    while(True): # Forever loop
        print("\n(A)uto  (b*)  (B)atch  (Q)uit  ")
        buf = input(": ").strip()

        bstar = sat.bstar

        try:
            bstar = float(buf)
            delta_el(sat,bstar=bstar)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)       # print new elements
        except:
            buf = buf.upper()
        
            # AUTO
            if (buf == 'A'):
                # update Bstar within limits
                bmax = bstar * 1.5
                bmin = bstar * 0.5
                if (bstar < 0):
                    bmin = bstar * 0.5
                    bmin = bstar * 1.5

                ## This doesn't seem very "auto"
                # TODO: make "Search" version
                # print("\nbmax [{:.8f}".format(bmax))
                # buf = input("]: ").strip()
                # bmax = float(buf)

                # print("\nbmin [{:.8f}".format(bmin))
                # buf = input("]: ").strip()
                # bmin = float(buf)

                while((bmax - bmin) > 1.e-12):
                    bstep = (bmax - bmin) / 20
                    for bk in np.arange(bmin, bmax, bstep):
                        delta_el(sat,bstar=bk)

                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(sat, rd, ll, odata)
                        if (rms < sum):
                            sum = rms
                            bstar = bk
                       # END for bk
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                    bmin = bstar - bstep
                    bmax = bstar + bstep

                delta_el(sat,bstar=bstar)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
                print_el(sat, satmeta)       # print new elements

            # BATCH
            elif (buf == 'B'):
                # create batch file
                with open("batch.inp", "w") as fp:
                    fp.write("{:s}\n\n\n\n\n\n\n".format("s"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n\n\n\n".format("s"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n\n\n\n".format("s"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n\n\n\n".format("s"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n\n\n\n".format("z"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n\n\n\n".format("z"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("b"))
                    fp.write("{:s}\n\n\n".format("a"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("w"))
                    fp.write("{:s}\n".format("u"))
                    fp.write("{:s}\n".format("q"))
                    fp.write("{:s}\n".format("q"))
                file_in = "satfit {:s} ba<batch.inp".format(file)
                ## TODO: Fix this batch stuff for UNIX
                # system(file_in);
                # sprintf(file_in, "DEL batch.inp");
                # // system(file_in);
                # sprintf(file_in, "satfit %s", file);
                # // system(file_in);
                # sys.exit(0)
            elif (buf == 'Q'):
                return None


def move_epoch_to_previous_perigee(sat):
    # if previously at perigee, back up one revolution
    if (sat.mo < radians(0.1)):  # Use SGP4 elements
        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) - NOCON/sat.no_kozai*(1 + sat.mo/TWOPI) # Use SGP4 elements
    # if not previously at perigee, back up to perigee
    else:
        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) - sat.mo*NOCON/(sat.no_kozai*TWOPI) # Use SGP4 elements
    delta_t(sat,t2_jd)
    sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
    # refine perigee
    for i in range(0, 30):
        # go forward
        if (sat.mm > radians(359.9)):
            t1_delta = NOCON/sat.no_kozai*(sat.mm/TWOPI - 1) # Use instantaneous mean elements except for no_kozai
        # back up
        else:
            t1_delta = sat.mm*NOCON/(sat.no_kozai*TWOPI) # Use instantaneous mean elements except for no_kozai
        t1_jd = t2_jd - t1_delta
        delta_t(sat,t1_jd)
        sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
    delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t1_jd)
    return sat

def maneuver(sat, rd, ll, odata, obs_meta, sum, iod_line):
    # Make a copy of original sat
    save_sat = copy.deepcopy(sat)
    nobs = len(odata)

    while(True): # Forever loop
        print("\nBoundary#  (P)erigee  (A)pogee  (O)b  (E)nergy  (R)estore  (Q)uit  ")

        buf = input(": ").strip()

        try:
            p = int(buf)
 
            # store old obs
            iod_linex = copy.copy(iod_line)
            llx       = copy.copy(ll)
            rdx       = copy.copy(rd)
            odatax    = copy.copy(odata)

            for i in range(p, nobs):
                iod_line[i-p] = iod_line[i-1]
                ll[i-p] = ll[i-1]
                rd[i-p] = rd[i-1]
                odata[i-p] = odata[i-1]
            nobsx = nobs
            nobs  = nobs - p + 1

            out = 0
            sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)
            print_el(sat, satmeta)
            print("\nperiod = {:f} days".format(NOCON/sat.no_kozai))
        except:
            buf = buf.upper()

            # Perigee
            if (buf == 'P'):
                # double time;      // days
                print("\n(P)revious  (N)ext  (Q)uit  ")

                buf = input(": ").strip()
                buf = buf.upper()
                if (buf  == 'Q'):
                    continue # Back to top of While loop
                if (buf == 'P'):
                    # if previously at perigee, back up one revolution
                    if (sat.mo < radians(0.1)):  # Use SGP4 elements
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) - NOCON/sat.no_kozai*(1 + sat.mo/TWOPI) # Use SGP4 elements
                    # if not previously at perigee, back up to perigee
                    else:
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) - sat.mo*NOCON/(sat.no_kozai*TWOPI) # Use SGP4 elements
                # NEXT: advance one revolution
                if (buf == 'N'):
                    # if previously at perigee, go forward one revolution
                    if (sat.mo > radians(359.9)):  # Use SGP4 elements
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) + NOCON/sat.no_kozai*(2 - sat.mo/TWOPI) # Use SGP4 elements
                    # if not previously at perigee, go forward to perigee
                    else:
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) + NOCON/sat.no_kozai*(1 - sat.mo/TWOPI) # Use SGP4 elements
                # t2 = Date(time=time)
                # move to time and ma at perigee
                delta_t(sat,t2_jd)
                sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
                # refine perigee
                for i in range(0, 30):
                    # go forward
                    if (sat.mm > radians(359.9)):
                        t1_delta = NOCON/sat.no_kozai*(sat.mm/TWOPI - 1) # Use instantaneous mean elements except for no_kozai
                    # back up
                    else:
                        t1_delta = sat.mm*NOCON/(sat.no_kozai*TWOPI) # Use instantaneous mean elements except for no_kozai
                    t1_jd = t2_jd - t1_delta
                    delta_t(sat,t1_jd)
                    sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
                print("\nPERIGEE")
                # Reinitialize satrec from current instantaneous mean elements and time (except mean motion)
                delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t1_jd)
                print_el(sat, satmeta)       # print new elements

                # perigee residual
                delta_t(sat,(sat.jdsatepoch + sat.jdsatepochF))
                save_sat = delta_t(save_sat,(sat.jdsatepoch + sat.jdsatepochF))
                delr = sat.rr - save_sat.rr           # compare sat and satm perigees
                print("\nperigee delta {:5.0f}".format(norm(delr)*sat.radiusearthkm))

            # Apogee
            elif (buf == 'A'):
                # time to travel from perigee to apogee, days
                print("\n(P)revious  (N)ext  (Q)uit  ")
                buf = input(": ").strip()
                buf = buf.upper()
                if (buf  == 'Q'):
                    continue # Back to top of While loop

                # Previous
                elif (buf == 'P'):
                    # if previously at or past apogee and past next perigee, back up to apogee
                    if (sat.mo < radians(180.1)):
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) - 0.5*NOCON/sat.no_kozai*(1 + sat.mo/pi) # Use SGP4 elements
                    # if previously past apogee and before next perigee, back up to apogee
                    else:
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) + 0.5*NOCON/sat.no_kozai*(1 - sat.mo/pi) # Use SGP4 elements
                # NEXT: advance to apogee
                elif (buf == 'N'):
                    # if previously at or past apogee and before perigee, go forward to apogee
                    if (sat.mo > radians(179.9)):
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) + 0.5*NOCON/sat.no_kozai*(3 - sat.mo/pi) # Use SGP4 elements
                    # if previously past apogee and past next perigee, go forward to apogee
                    else:
                        t2_jd = (sat.jdsatepoch + sat.jdsatepochF) + 0.5*NOCON/sat.no_kozai*(1 - sat.mo/pi) # Use SGP4 elements

                # move time and ma at apogee
                delta_t(sat,t2_jd)
                sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
                # loop to refine apogee, find when mean anomaly = pi
                for i in range(0, 30):
                    t1_jd = t2_jd + 0.5*NOCON/sat.no_kozai*(1 - sat.mm/pi) # Use instantaneous mean elements except for no_kozai
                    delta_t(sat,t1_jd)
                    sat.mm = posradang(sat.mm) # Use instantaneous mean element FIXME: cython SGP4 allows these to be negative sign
                print("\nAPOGEE")
                # Reinitialize satrec from current instantaneous mean elements and time (except mean motion)
                
                delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t1_jd)
                print_el(sat, satmeta)       # print new elements

                # apogee residual
                delta_t(sat,sat.jdsatepoch)
                save_sat = delta_t(save_sat,sat.jdsatepoch)
                delr = sat.rr - save_sat.rr           # compare sat and satm perigees
                print("\napogee delta {:5.0f}".format(norm(delr)*sat.radiusearthkm))

            # O(b) Pseudo observation?
            elif (buf == 'O'):
                # move up one to make room for pseudo ob
                # FIXME: Probably a pythonic way to do this 
                # for i in range(nobs-1, -1, -1):
                #     iod_line[i+1] = iod_line[i]
                #     ll[i+1] = ll[i]
                #     rd[i+1] = rd[i]
                #     odata[i+1] = odata[i]

                delta_t(sat,sat.jdsatepoch)
                # ll is unit vector in same direction as satm.rr
                ll = np.insert(ll,0,unit_vector(sat.rr),axis=0) # FIXME: This doesn't appear to check out with satfit.cpp
                
                # rd is unit vector (1 er) in same direction as satm.rr
                rd = np.insert(rd,0,unit_vector(sat.rr),axis=0)
                # odata holds epoch, ra, dc, obscode data

                # odata[0][0] = sat.jdsatepoch
                (ra, dc) = zrll(sat, rd[0]) # Get predicted position
                odata = np.insert(odata,0,np.array((sat.jdsatepoch,radians(ra),radians(dc))),axis=0)
                # odata[0][1] = radians(ra)
                # odata[0][2] = radians(dc)
                # odata[0][3] = 0000

                # print pseuso ob
                #           1         2         3         4         5         6
                # 0123456789012345678901234567890123456789012345678901234567890
                # 36868 10 039A   2018 G 20100909203108300 17 25 2216290+034350 57

                # norad = int(iod_line[1])
                # yy    = int(iod_line[1] + 6)
                # desig = "{:4s}".format(iod_line[1] + 9)
                # desig[4] = '\0'
                yy = int(sat.intldesg[2:4])
                desig = sat.intldesg[5:]
                t1 = Date(time=(sat.jdsatepoch))
                ra /= 15
                (ra, rm) = modf(ra)
                rm *= 60000
                # dm = fabs(fmod(dc, dm)*6000)
                (dc, dm) = modf(dc)
                dm *= 6000
                print("\nPSEUDO OB:")
                iod_line[0] = "\n{:05d} {:02d} {:4s}   0000 P {:4.0f}{:02.0f}{:02.0f}{:02.0f}{:02.0f}{:05.0f} 16 25 {:02.0f}{:05.0f}{:+03.0f}{:04.0f}\n".format(
                    sat.satnum, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss*1000, ra, rm, dc, dm)
                print(iod_line[0])
            # Energy 
            # FIXME: Need to finish porting this
            elif (buf == 'E'):
                buf = input("\ndelta specific energy(m^2/s^2): ")
                buf = buf.strip()

                try:
                    delta_t(sat,sat.jdsatepoch)
                    dE = float(buf)
                    dE /= 11300.168   # er^2 / min^2
                    vec = np.cross(sat.rr, sat.vv)
                    mu = norm(vec)
                    mu = mu*mu
                    mu = mu / sat.a # Was sat.aodp
                    mu = mu / (1 - sat.eo*sat.eo)
                    E2 = -0.5*mu / sat.a # Was sat.aodp
                    VV = sqrt(2*(E2 + dE + mu/norm(satm.rr)))  # new velocity magnitude
                    dV = unit_vector(sat.vv)   # velocity direction
                    dev = VV * dV               # new velocity vector
                    # sat = satm
                    sat.rv2el(sat.rr, vec)      # State vectors to mean elements #FIXME: python-SGP4
                    satE = delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t1_jd)
                    print_el(sat, satmeta)              # print new elements
                except:
                    pass # TODO: need better exception handling here
            # Restore 
            elif (buf == 'R'):
                print("\nRestore (O)bs  (E)lements  (Q)uit  ")
                buf = input(": ").strip()
                buf = buf.upper()

                # re-store old obs
                if (buf == 'O'):
                    nobs = nobsx
                    ll = llx
                    rd = rdx
                    odata = odatax
                # re-store old TLE
                elif (buf == 'E'):
                    # replace working elements with original
                    satm = save_sat        # original elements maneuverd to a node
                    sat  = save_sat        # new fit elements
                    print_el(sat, satmeta)         # print original elements
            # QUIT Maneuver        
            elif (buf == 'Q'):
                return sat


def write_el(db, sat, satmeta, rd, ll, odata, obs_meta, sum, start_rms):
    """
    Function to interact with source, and updated TLEs, including insert to TruSat database.

    TODO: Not sure what it should 'return' as it nominally does not update in the working variables.
    """
    global startDate

    # save_sat = copy.deepcopy(sat)
    while(True): # Forever loop
        # # When making a TLE, make sure to use the most recent obs_time available for setting epoch
        # epoch_jd = odata[-1][0]
        # print_satrec_coe(sat)

        # # tau = t.epoch - Mean Anomaly/Mean motion - Time of perigee crossing
        # tau = sat.jdsatepoch - (sat.mo/sat.no_kozai)
        # new_ma = sat.no_kozai*(epoch_jd-tau) # Solve for new Mean Anomaly at last observation time
        # sat = delta_t(sat,epoch_jd)

        # sat2 = re_initsat(sat, new_ma, epoch_jd)
        # sat2 = delta_t(sat2,epoch_jd)
        # print("After reinit")
        # print_satrec_coe(sat2)

        # print("sat1.rr {} sat1.vv {}".format(sat.rr,sat.vv))
        # print("sat2.rr {} sat2.vv {}".format(sat2.rr,sat2.vv))

        newTLE = make_tle_from_SGP4_satrec(sat,satmeta,classification="T")
        # tle = sat.epoch
        # ii = degrees(sat.inclo)
        # om = degrees(sat.nodeo)
        # ec = sat.ecco
        # ww = degrees(sat.argpo)
        # ma = degrees(sat.mo)
        # nn = sat.no_kozai / NOCON
        # c2 = sat.c2
        # bstar = sat.bstar

        print("\n(U)pdate (V)iew (A)ppend (O)riginal (R)estore", end='')
        if (db):
            print("\nDATABASE: (I)nsert new TLE (L)oop insert ra(W) insert (Q)uit", end='')
        else:
            print(" (Q)uit")

        buf = input(": ").strip().upper()

        # Append
        if (buf == 'A'):
            write_tle(file) # FIXME replace with tle_util version
        
        # View
        elif (buf == 'V'):
            print_el(sat, satmeta)

        # View Original
        elif (buf == 'O'):
            print_el(sat, satmeta)       # print new elements

        # Restore
        elif (buf == 'R'):
            # replace TLE in file with original
            with open(file, "w") as fp:
                fp.write("{:s}".format(TLE.name))
                fp.write("{:s}".format(TLE.line1))
                fp.write("{:s}".format(TLE.line2))
                for range in (nobs): # FIXME More pythonic way
                    fp.write("{:s}".format(iod_line[i]))
            # replace working elements with original
            sat = save_sat
            print_el(sat, satmeta)            # print original elements
        # Update
        elif (buf == 'U'):
            write_tle()     # updates lines and prints to screen # FIXME tle_util screen version
            #  print_file

            name = name.strip()
            pgee = (1-ec)*sat.aodp
            agee = (1+ec)*sat.aodp
            buf = "{:d} X {:d} km".format(int(((pgee-1)*6378.135)),
                                       int(((agee-1)*6378.135)))
            with open(file, "w") as fp:
                fp.write("{:<50s}{:19s}\n".format(name, buf))
                fp.write("{:s}\n".format(line1))
                fp.write("{:s}\n\n".format(line2))
                for i in range(nobs): # FIXME More pythonic way
                    fp.write("{:s}".format(iod_line[i]))
        # Insert TLE and TLE_process results to database
        elif (buf == 'I' or buf == 'L' or buf == 'W'):
            newTLE_process = {}
            (sum, newTLE_process) = calc_fit(sat, rd, ll, odata, obs_meta, sum, newTLE_process)

            remarks = input("TLE Remarks for TLE_process record: ").strip()

            result = db.addTruSatTLE(newTLE, newTLE_process, satmeta.parent_tle_id, start_rms, sum, remarks)
            if (result):
                startDate = newTLE.epoch_datetime
                if (buf == 'I'):
                    object_search(db, startDate=startDate, object=newTLE.satellite_number)
                elif (buf == 'L'):
                    object_process(db, startDate=startDate, object=newTLE.satellite_number)
                elif (buf == 'W'):
                    # raw_search(db)
                    return True
            else:
                log.error("DB Insert failed.")
        # QUIT write_el
        elif (buf == 'Q'):
            return True

def time_func(sat, rd, ll, odata, obs_meta, sum):
    return sat
    # while(True): # Forever loop
    #     nobs = len(odata)
    #     ec = sat.ecco
    #     ww = degrees(sat.argpo)
    #     nn = sat.no_kozai / NOCON
    #     xns = 2160 * sat.bstar * sat.c2 * sat.no_kozai / NOCON
    #     if (nobs > 0):
    #         time2 = odata[nobs - 1][0]
    #     else:
    #         t2.now()
    #         time2 = t2.jd

    #     delta_t(sat,time2) # FIXME: - this needs a delta time
    #     # sat.delta_t(time2) # FIXME python-SGP4
    #     z2 = sat.rr[2]

    #     while(z1 < 0 or z2 > 0):
    #         time1 = time2
    #         z1 = z2
    #         time2 -= .01
    #         delta_t(sat,time2) # FIXME: - this needs a delta time
    #         # sat.delta_t(time2)  # FIXME python-SGP4
    #         z2 = sat.rr[2]

    #     while(time1 - time2 > 1e-9):
    #         time3 = (time1 + time2) / 2
    #         delta_t(sat,time3) # FIXME: - this needs a delta time
    #         # sat.delta_t(time3)  # FIXME python-SGP4
    #         z3 = sat.rr[2]
    #         if (z3 < 0):
    #             time2 = time3
    #         else:
    #             time1 = time3

    #     t1 = Date(time=time2)
    #     t1.input()
    #     delta_t(sat,t1.jd) # FIXME: - this needs a delta time
    #     # sat.delta_t(t1.jd)             # advance to node # FIXME: python-SGP4
    #     # sat.rv2el(sat.rr, sat.vv)      # sgp4 elements # FIXME: python-SGP4
    #     tle = t1.tle
    #     sat.jd = t1.jd

    #     ii = degrees(sat.xincl)
    #     om = fmod(degrees(sat.xnodeo), 360)
    #     ec = sat.eo
    #     ww = fmod(degrees(sat.omegao), 360)
    #     ma = fmod(degrees(sat.xmo), 360)
    #     nn = sat.xno / NOCON
    #     bstar = sat.bstar
    #     print_el(sat)       # print new elements
    #     sum = find_rms(sat, rd, ll, odata, sum)
    #     print("\nrms{:12.5f}".format(sum))

    #     uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)


def initsat(TLE,gravconst="wgs72"):
    """Initialize SGP4 Satellte() from raw variables
    
    Could do this directly from the database, but treating it as a general function
    since we may get TLEs from other sources, or initialize other variables 
    other than SGP4-satrec.
    
    Inputs:
        TLE       TLE class variable (from DB or file)
        gravconst Which gravity constants to use in SGP4

    Output:
        satrec  Initialized SGP Satellite() class variable
        -or-    False on error
    """    
    # Initialize satrec variables, modeled after twoline2rv()
    satrec              = Satrec()

    # Line 1
    satnum              = TLE.satellite_number

    meta                = SatrecMeta(satnum=satrec.satnum)
    meta.line0          = TLE.line0
    meta.line1          = TLE.line1

    meta.classification = TLE.classification

    meta.intldesg       = TLE.designation
    meta.intldesg       = TLE.designation

    meta.epochyr        = 0 if not TLE._epoch_year else TLE._epoch_year
    meta.epochdays      = 0 if not TLE._epoch_day else TLE._epoch_day
    ndot                = TLE.mean_motion_derivative
    nddot               = TLE.mean_motion_sec_derivative
    bstar               = TLE.bstar
    meta.ephtype        = TLE.ephemeris_type
    meta.elnum          = 0 if not TLE.element_set_number else TLE.element_set_number

    # Line 2
    meta.line2          = TLE.line2
    inclo               = TLE.inclination_radians  # rad
    nodeo               = TLE.raan_radians         # rad
    ecco                = TLE.eccentricity   
    argpo               = TLE.arg_perigee_radians  # rad
    mo                  = TLE.mean_anomaly_radians # rad

    # TODO: Once mean_motion_radians_per_minute is added to the DB, use it directly here
    no_kozai            = TLE.mean_motion_orbits_per_day * NOCON # rad/min
    meta.revnum         = 0 if not TLE.orbit_number else TLE.orbit_number

    # Derived quantities
    # jdsatepoch     = floor(TLE.jdsatepoch)     # Julian date
    # jdsatepochF    = TLE.jdsatepoch % 1        # Julian date fraction
    meta.jdSGP4epoch      = TLE.jdsatepoch - 2433281.5
    meta.epoch_datetime   = TLE.epoch_datetime # Python datetime

    # Pass the source tle_id through the SGP4 class variable, for TLE genealogy
    meta.parent_tle_id     = TLE.tle_id

    # SGP4 mode variables
    operationmode  = u'i' # Unicode for cython
    # satrec.error          = 0

    if (gravconst == "wgs72old"):
        whichconst = earth_gravity.wgs72old
    elif (gravconst == "wgs84"):
        whichconst = earth_gravity.wgs84
    else:
        # Most popular const used by TLEs
        whichconst = earth_gravity.wgs72
    # satrec.whichconst     = whichconst  # Python extension: remembers its consts
    # satrec.whichconst = gravconst
    satrec.sgp4init(satnum, 
                            meta.jdSGP4epoch, # epoch time in days from jan 0, 1950. 0 hr
                            bstar, ndot, nddot, ecco, argpo,
                            inclo, mo, no_kozai, nodeo)
    
    # rtn_code = satrec.sgp4init("wgs72", operationmode, satrec.satnum, 
    #          satrec.jdsatepoch - 2433281.5, # epoch time in days from jan 0, 1950. 0 hr
    #          satrec.bstar, satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, 
    #          satrec.inclo, satrec.mo, satrec.no_kozai, satrec.nodeo, satrec)

    return satrec, meta

    # FIXME figure out how to deal with satrec.error codes
    print("Satrec.error {}".format(satrec.error))
    if (satrec.error == 0):
        return satrec
    elif (satrec.error == 1):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er")
        return False
    elif (satrec.error == 2):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("mean motion less than 0.0")
        return False
    elif (satrec.error == 3):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("pert elements, ecc < 0.0  or  ecc > 1.0")
        return False
    elif (satrec.error == 4):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("semi-latus rectum < 0.0")
        return False
    elif (satrec.error == 5):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("epoch elements are sub-orbital")
        return False
    elif (satrec.error == 6):
        log.error("sgp4init error {}".format(satrec.error))
        log.error("satellite has decayed")
        return False
    else:
        log.error("sgp4init error {}".format(satrec.error))
        log.error("Unknown error code")
        return False


def raw_search(db=False):
    global iod_line #FIXME: get rid of this global
    classification='U' # FIXME: Default to this unless overridden

    while(True):
        try:
            iod_obs_inp = input("\nEnter 1 or more IOD Obs IDs: ")
            iod_obs_inp = iod_obs_inp.strip()
            iod_obs = iod_obs_inp.split(' ')
        except:
            pass

        IOD_candidates = db.selectIODListat(iod_obs[0])
        if (len(iod_obs)==1):
            print("obs_id satnum STA  user              obs_time  RA   DEC")
            for r in IOD_candidates:
                print("{ID:7} {NUM:5} {STA:4} {USR} {TIME} ra:{RA:<8.4f} dec:{DEC:<8.4}".format(
                    ID=r.obs_id, 
                    NUM=r.ObjectNumber, 
                    STA=r.Station, 
                    USR=r.UserString, 
                    TIME=r.obs_time.isoformat(),
                    RA=r.RA,
                    DEC=r.DEC)
                )
            continue
        elif(len(iod_obs)==0):
            print("No observations found for obs_id {}".format(iod_obs))
            continue
        else:
            IODsq = db.selectIODlist(iod_obs)

        IODs = []
        iod_line = []
        satnum = None
        for i in IODsq:
            if (satnum and satnum != i.ObjectNumber ):
                # Send the user back to the sat selection list
                log.error("You selected observations for multiple satellites.  Try again")
                iod_obs = [iod_obs[0]]
                continue
            else:
                IODs.append(i)
                iod_line.append(i.iod_string)
                satnum = i.ObjectNumber

        try:
            try:
                tle_id_inp = input("\nEnter ref TLE ID: ")
                tle_id_inp = tle_id_inp.strip()
                tle_id = int(tle_id_inp)
            except:
                log.warning("NO TLE found for id '{}' ".format(tle_id))
                continue
            TLE = db.selectTLEid(tle_id)
            print("Found tle_id {} for reference:".format(TLE.tle_id))
            print("{}".format(TLE.name))
            print("{}".format(TLE.line1))
            print("{}".format(TLE.line2))
        except:
            log.warning("NO TLE found for id '{}' ".format(tle_id))
            continue

        Stations = db.getStationDictforIODs(IODs)
        if (not len(Stations)):
            log.warning("NO Station data found for observations.")
            continue
        break


    # // read all observations from input
    # // Observed topocentric vectors for each observer position.
    # // Equatorial coordinates.
    # TODO: Audit this against everything done in read_obs starting on line 1585
    # get line-of-sight vectors
    # (odata, ll, rd) = read_obs(IOD_Records)
    (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

    (sat, satmeta) = initsat(TLE)
    satmeta.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()
    satmeta.parent_tle_id = TLE.tle_id

    # Make a copy of original sat
    # save_sat = copy.deepcopy(sat)

    # // maneuvering sat
    # satm = copy.deepcopy(sat)

    # // calculate uu, degrees, for search
    uu = longitude(sat)

    # DB results already sorted, make sure the file ones are
    if(not db):
        [iod_line, odata, ll, rd, obs_meta] = sort(iod_line, odata, ll, rd, obs_meta)

    nobs = len(IODs)
    sum = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    print("\n{} Observations Found".format(nobs))
    print("\nStarting rms{:12.5f}".format(sum))

    # print dates
    # t2 = Date()
    # print("\TLE   : {:d}".format(t2.doy))
    # if(nobs > 0):
    #     t3 = Date(time=(odata[nobs - 1][0]))
    #     print("LAST OB : {:d}".format(t3.doy))
    age = odata[nobs - 1][0] - (sat.jdsatepoch + sat.jdsatepochF)
    print("TLE Age from last OBS: {:.2f} days".format(age))    

    # Accept a new command
    accept_command(db, sat, satmeta, rd, ll, odata, obs_meta, sum, uu, iod_line)


def iod_search(db=False):
    global iod_line #FIXME: get rid of this global
    classification='U' # FIXME: Default to this unless overridden

    while(True):
        try:
            iod_obs_inp = input("\nEnter 1 or more IOD Obs IDs: ")
            iod_obs_inp = iod_obs_inp.strip()
            iod_obs = iod_obs_inp.split(' ')
        except:
            pass

        IOD_candidates = db.selectIODListat(iod_obs[0])
        if (len(iod_obs)==1 and iod_obs_inp.find('.') < 0 ):
            print("obs_id satnum STA  user              obs_time  RA   DEC")
            for r in IOD_candidates:
                print("{ID:7} {NUM:5} {STA:4} {USR} {TIME} ra:{RA:<8.4f} dec:{DEC:<8.4}".format(
                    ID=r.obs_id, 
                    NUM=r.ObjectNumber, 
                    STA=r.Station, 
                    USR=r.UserString, 
                    TIME=r.obs_time.isoformat(),
                    RA=r.RA,
                    DEC=r.DEC)
                )
            continue
        elif(len(iod_obs)==0):
            print("No observations found for obs_id {}".format(iod_obs))
            continue
        elif(iod_obs_inp.find('.')>0):
            iod_obs = [int(iod_obs_inp.replace(".",""))]
            IODsq = db.selectIODlist(iod_obs)
        else:
            IODsq = db.selectIODlist(iod_obs)

        IODs = []
        iod_line = []
        satnum = None
        for i in IODsq:
            if (satnum and satnum != i.ObjectNumber ):
                # Send the user back to the sat selection list
                log.error("You selected observations for multiple satellites.  Try again")
                iod_obs = [iod_obs[0]]
                continue
            else:
                IODs.append(i)
                iod_line.append(i.iod_string)
                satnum = i.ObjectNumber

        try:
            if (classification=='T'):
                TLE = db.selectTLEEpochNearestDate(IODs[-1].obs_time, satnum,classification=classification)    # FIXME: Probably want to call elfind in a rebuild case
            else:
                TLE = db.selectTLEEpochNearestDate(IODs[-1].obs_time, satnum)    
            print("Using tle_id {} as reference:".format(TLE.tle_id))
            print("{}".format(TLE.name))
            print("{}".format(TLE.line1))
            print("{}".format(TLE.line2))
        except:
            log.warning("NO TLE found for object '{}' with EPOCH before {}".format(IODs[0].ObjectNumber,IODs[0].obs_time))
            continue

        Stations = db.getStationDictforIODs(IODs)
        if (not len(Stations)):
            log.warning("NO Station data found for observations.")
            continue
        break


    # // read all observations from input
    # // Observed topocentric vectors for each observer position.
    # // Equatorial coordinates.
    # TODO: Audit this against everything done in read_obs starting on line 1585
    # get line-of-sight vectors
    # (odata, ll, rd) = read_obs(IOD_Records)
    (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

    (sat, satmeta) = initsat(TLE)
    parent_tle_id = 0

    # FIXME: Figure out how this is used
    satmeta.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()

    # Make a copy of original sat
    # save_sat = copy.deepcopy(sat)

    # // maneuvering sat
    # satm = copy.deepcopy(sat)

    # // calculate uu, degrees, for search
    uu = longitude(sat)

    # DB results already sorted, make sure the file ones are
    if(not db):
        [iod_line, odata, ll, rd] = sort(iod_line, odata, ll, rd)

    nobs = len(IODs)
    sum = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    print("\n{} Observations Found".format(nobs))
    print("\nStarting rms{:12.5f}".format(sum))

    # print dates
    # t2 = Date()
    # print("\TLE   : {:d}".format(t2.doy))
    # if(nobs > 0):
    #     t3 = Date(time=(odata[nobs - 1][0]))
    #     print("LAST OB : {:d}".format(t3.doy))
    age = odata[nobs - 1][0] - (sat.jdsatepoch + sat.jdsatepochF)
    print("TLE Age from last OBS: {:.2f} days".format(age))    

    # Accept a new command
    accept_command(db, sat, satmeta, rd, ll, odata, obs_meta, sum, uu, iod_line)


def object_search(db=False,startDate=False,object=False):
    global iod_line #FIXME: get rid of this global
    global TLE_ref

    while(True):
        if not object:
            try:
                object_inp = input("\nEnter norad number of object: ")
                object = int(object_inp.strip())
            except:
                pass

        if not (startDate):
            result = db.findDateNewestTLE(object)
            if (result):
                startDate = result
                classification='T'
                print("Continuing processing at: {}".format(startDate))
            else:
                startDate = datetime(1957,10,4,19,28,34,0)
                classification='U'
                print("Processing since the beginning of time")
        else:
                # If we're entering the function with a start date, we've created a TTLE
                classification='T'

        # IOD_candidates = db.findObservationCluster(object,startDate=startDate,minObserverCount=1)
        IOD_candidates = db.findLastNIODs_noTLE(object, IOD_count=100)
        if (len(IOD_candidates)==0):
            print("No observations found for norad_number {}".format(object))
            main(db)

        IODs = db.selectIODlist(IOD_candidates)

        iod_line = []
        for i in IODs:
            iod_line.append(i.iod_string)

        try:
            tle_id = False
            while(True):
                if (not tle_id):
                    if (classification=='T'):
                        TLE = db.selectTLEEpochNearestDate(startDate, object)    # FIXME: Probably want to call elfind in a rebuild case
                    else:
                        TLE = db.selectTLEEpochNearestDate(IODs[-1].obs_time, object)    # FIXME: Probably want to call elfind in a rebuild case
                print("Using tle_id {} as reference:".format(TLE.tle_id))
                print("{}".format(TLE.name))
                print("{}".format(TLE.line1))
                print("{}".format(TLE.line2))

                print()
                print("[Return] for current TLE.  TLE [I]d#  Nearest to [F]irst/[L]ast Obs  [P]revious epoch\n[E]xternal  [Q]uit to main")
                TLE_pref = input(" :").strip()

                if (TLE_pref == ""):
                    break

                try: 
                    tle_id = int(TLE_pref)
                    TLE = db.selectTLEid(tle_id)
                except ValueError:
                    TLE_pref = TLE_pref.upper()
                    if (TLE_pref == "F"): # Nearest to first obs
                        TLE = db.selectTLEEpochNearestDate(IODs[0].obs_time, object)
                    elif (TLE_pref == "L"): # Nearest to last obs
                        TLE = db.selectTLEEpochNearestDate(IODs[-1].obs_time, object)
                    elif(TLE_pref == "P"): # Previous Epoch
                        TLE = db.selectTLEEpochBeforeDate(IODs[0].obs_time, object)
                    elif(TLE_pref == "E"): # External
                        try:
                            TLE = TLE_ref.Satellites[object]
                            TLE.tle_id = 0
                            tle_id = True
                        except:
                            if (TLE_ref):
                                log.warning("No external TLE found for object '{}'".format(object))
                    elif(TLE_pref == "I"): # TruSat TLE ID                        
                        TLE_id_inp = input("Enter TruSat TLE id: ")
                        TLE_id_inp = int(TLE_id_inp.strip())
                        TLE = db.selectTLEid(TLE_id_inp)
                    elif(TLE_pref == "Q"): # Quit to main menu
                        main(db)
        except:
            try:
                TLE = TLE_ref.Satellites[object]
                TLE.tle_id = 0
                tle_id = True
            except:
                if (TLE_ref):
                    log.warning("No external TLE found for object '{}'".format(object))
                continue

        Stations = db.getStationDictforIODs(IODs)
        if (not len(Stations)):
            log.warning("NO Station data found for observations.")
            continue
        break

    # TODO: Turn this, and the function from iod_search into a separate function

    # // read all observations from input
    # // Observed topocentric vectors for each observer position.
    # // Equatorial coordinates.
    # TODO: Audit this against everything done in read_obs starting on line 1585
    # get line-of-sight vectors
    # (odata, ll, rd) = read_obs(IOD_Records)
    (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

    (sat, satmeta) = initsat(TLE)
    satmeta.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()
    satmeta.parent_tle_id = TLE.tle_id

    # Make a copy of original sat
    # save_sat = copy.deepcopy(sat)

    # // maneuvering sat
    # satm = copy.deepcopy(sat)

    # // calculate uu, degrees, for search
    uu = longitude(sat)

    print(f"thetag: {satmeta.thetag}  uu: {uu}  gsto: {sat.gsto/DE2RA}")
    # DB results already sorted, make sure the file ones are
    if(not db):
        [iod_line, odata, ll, rd] = sort(iod_line, odata, ll, rd)

    nobs = len(IODs)
    sum = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    print("\n{} Observations Found".format(nobs))
    print("\nStarting rms{:12.5f}".format(sum))

    # print dates
    # t2 = Date()
    # print("\TLE   : {:d}".format(t2.doy))
    # if(nobs > 0):
    #     t3 = Date(time=(odata[nobs - 1][0]))
    #     print("LAST OB : {:d}".format(t3.doy))
    age = odata[nobs - 1][0] - (sat.jdsatepoch + sat.jdsatepochF)
    print("TLE Age from last OBS: {:.2f} days".format(age))    

    # Accept a new command
    accept_command(db, sat, satmeta, rd, ll, odata, obs_meta, sum, uu, iod_line)

def object_manual(db=False,startDate=False,object=False):
    global TLE_ref 

    """ Script-assisted manual processing """
    object = False
    NoTLEs = False

    objects = db.findObjectsWithIODsSubmittedAfterTLE()
    print("\n{} objects remaining with IODs newer than the latest TLE".format(len(objects)))
    if (len(objects) == 0):
        days = 30
        objects = db.findObjectsWithIODsButNoTLEs(days=days)
        num_objects = len(objects)
        print("Found {} objects needing TLEs in the last {} days".format(num_objects, days))
        NoTLEs = True

    for object in objects:
        result = True # Start the loop off
        while (result):
            if (NoTLEs):
                IOD_candidates = db.findLastNIODs_noTLE(object,IOD_count=20)
            else:
                IOD_candidates = db.findIODsSubmittedAfterPenultimateTLE(object)
            if (len(IOD_candidates)==0):
                print("No more observations found for norad_number {}".format(object))
                result = False
                continue

            IODs = db.selectIODlist(IOD_candidates)

            iod_line = []
            for i in IODs:
                iod_line.append(i.iod_string)

            try:
                TLE = db.selectTLEEpochBeforeDate(IODs[-1].submitted, object)    # FIXME: Probably want to call elfind in a rebuild case
                print("Using tle_id {} as reference:".format(TLE.tle_id))
                print("{}".format(TLE.name))
                print("{}".format(TLE.line1))
                print("{}".format(TLE.line2))
            except:
                log.warning("NO TLE found for object '{}' with EPOCH before {}".format(object,IODs[0].submitted))
                try:
                    TLE = TLE_ref.Satellites[object]
                    TLE.tle_id = 0
                    print("Using external TLE as reference:")
                    print("{}".format(TLE.name))
                    print("{}".format(TLE.line1))
                    print("{}".format(TLE.line2))
                except:
                    if (TLE_ref):
                        log.warning("NO external TLE found for object '{}'".format(object))
                    result = False
                    continue

            Stations = db.getStationDictforIODs(IODs)
            if (not len(Stations)):
                log.warning("NO Station data found for observations.")
                result = False
                continue

            (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

            (sat, satmeta) = initsat(TLE)
            satmeta.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()
            satmeta.parent_tle_id = TLE.tle_id

            # // calculate uu, degrees, for search
            uu = longitude(sat)

            nobs = len(IODs)

            start_rms = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
            print("\n{} Observations Found".format(nobs))
            print("\nStarting rms{:12.5f}".format(start_rms))

            # Accept a new command
            accept_command(db, sat, satmeta, rd, ll, odata, obs_meta, start_rms, uu, iod_line)
            result = False


def object_process(db=False,startDate=False,object=False):
    """ object_process: Loop through all observations for a particular object """
    result = True # Start the loop off
    while (result):
        IOD_candidates = db.findObservationCluster(object,startDate=startDate,minObserverCount=1)
        if (len(IOD_candidates)==0):
            print("No more observations found for norad_number {}".format(object))
            main(db)

        IODs = db.selectIODlist(IOD_candidates)

        iod_line = []
        for i in IODs:
            iod_line.append(i.iod_string)

        try:
            TLE = db.selectTLEEpochBeforeDate(IODs[-1].obs_time, object,classification='T')    # FIXME: Probably want to call elfind in a rebuild case
            print("Using tle_id {} as reference:".format(TLE.tle_id))
            print("{}".format(TLE.name))
            print("{}".format(TLE.line1))
            print("{}".format(TLE.line2))
        except:
            log.warning("NO TLE found for object '{}' with EPOCH before {}".format(object,IODs[0].obs_time))
            main(db)

        Stations = db.getStationDictforIODs(IODs)
        if (not len(Stations)):
            log.warning("NO Station data found for observations.")
            main(db)

        (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

        (sat, satmeta) = initsat(TLE)
        satmeta.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()
        satmeta.parent_tle_id = TLE.tle_id
        nobs = len(IODs)

        start_rms = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
        print("\n{} Observations Found".format(nobs))
        # print("\nStarting rms{:12.5f}".format(start_rms))

        age = odata[-1][0] - (sat.jdsatepoch + sat.jdsatepochF)
        print("TLE Age from last OBS: {:.2f} days".format(age))    

        sum = print_fit(sat, satmeta, rd, ll, odata, start_rms)

        move_epoch_to_jd(sat,odata[-1][0])
        # print_el(sat)
        # sum = print_fit(sat, rd, ll, odata, sum)

        # // calculate uu, degrees, for search
        uu = longitude(sat)
        step(sat, rd, ll, odata, obs_meta, sum, uu, "S")

        # # // calculate uu, degrees, for search
        # uu = longitude(sat)
        # sat = step(sat, rd, ll, odata, sum, uu, "Z")

        print_el(sat, satmeta)
        sum = print_fit(sat, satmeta, rd, ll, odata, obs_meta, sum)

        newTLE = make_tle_from_SGP4_satrec(sat,satmeta,classification="T")

        newTLE_process = {}
        (sum, newTLE_process) = calc_fit(sat, rd, ll, odata, obs_meta,sum, newTLE_process)

        remarks = "object_process automatic TLE"

        result = db.addTruSatTLE(newTLE, newTLE_process, satmeta.parent_tle_id, start_rms, sum, remarks)
        if (result):
            startDate = newTLE.epoch_datetime
        else:
            log.error("DB Insert failed.")
            main(db)

def object_needs_TLE(db=False):
    objects = db.findObjectsWithIODsButNoTLEs()
    num_objects = len(objects)
    print("Found {} objects needing TLEs".format(num_objects))
    for item in objects:
        print(item)
    main(db)

def object_tle_test(db=False):
    """ object_tle_test: Loop through all observations for a particular object, generate TLE_process RMS from existing TLEs without fit. """
    object = False
    currentTime = datetime.now()

    objects = db.findObjectsWithIODsNotUsedInTTLEs()
    for object in objects:
        first_try = True
        result = True # Start the loop off
        while (result):
            time_loop_start = time()
            if first_try:
                # Set up the beginning of the loop
                try:
                    # object_inp = input("\nEnter norad number of object: ")
                    # object = int(object_inp.strip())

                    TLE1 = False
                    startTime = 0
                    TLE2 = db.findFirstIODandTLE(object)
                    endTime = TLE2.epoch_datetime

                    first_try = False

                except:
                    print("No more unprocessed TLEs found for norad_number {}".format(object))
                    # main(db)
                    result = False
                    continue

            IODs = db.selectIODlistSubmitRange(object, startTime, endTime)

            if not IODs:
                print("No obs for object {} between {} and {}".format(object,startTime,endTime))
                TLE2 = db.findNextUnprocessedTLE(object,endTime)
                if not TLE2:
                    print("No more unprocessed TLEs found for norad_number {}".format(object))
                    print("Grabbing newest IODs after epoch {} ".format(endTime))
                    endTime = currentTime
                    IODs = db.selectIODlistSubmitRange(object, startTime, endTime)
                    TLE2 = db.selectTLEEpochNearestDate(endTime,object)
                    if not IODs:
                        result = False
                        continue
                else:
                    endTime = TLE2.epoch_datetime
                    continue

            iod_line = []
            for i in IODs:
                iod_line.append(i.iod_string)

            Stations = db.getStationDictforIODs(IODs)
            if (not len(Stations)):
                log.warning("NO Station data found for observations.")
                result = False
                continue

            (odata, ll, rd, obs_meta, t1) = read_obssf(IODs, Stations) # Use the scott-campbell way while debugging against satfit.cpp

            try:
                if (TLE1):
                    (sat, satmeta) = initsat(TLE1)
                    satmeta.parent_tle_id = TLE1.tle_id
                else:
                    (sat, satmeta) = initsat(TLE2)
                    satmeta.parent_tle_id = TLE2.tle_id

                print("\n\nUsing tle_id {} as reference".format(sat.parent_tle_id))
                print("{}".format(sat.line0))
                print("{}".format(sat.line1))
                print("{}".format(sat.line2))
            except:
                log.warning("Failed to initsat()")
                # main(db)
                result = False
                continue

            sat.thetag = t1.thetag  # FIXME: Find a different way to get this, or how its used in anomaly()

            # start_rms = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
            # sum = print_fit(sat, rd, ll, odata, start_rms)

            time_to_calc_fit = time()
            TLE_process = {}
            (sum, TLE_process) = calc_fit(sat, rd, ll, odata, obs_meta, 0, TLE_process)
            print_calc_fit(sat, rd, ll, odata, obs_meta, sum, TLE_process)
            remarks = "McCants TLE loop baseline"
            time_calc_fit = time()

            result = db.addTruSatTLE(False, TLE_process, sat.parent_tle_id, sum, sum, remarks, satellite_number=object, tle_result_id=TLE2.tle_id)
            time_to_insert = time()

            if (result):
                TLE1 = TLE2
                TLE2 = db.findNextUnprocessedTLE(object,TLE1.epoch_datetime)
                if not TLE2:
                    print("No more unprocessed TLEs found for norad_number {} after epoch {}".format(object, TLE1.epoch_datetime))
                    # main(db)
                    result = False
                    continue
                else:
                    startTime = TLE1.epoch_datetime
                    endTime = TLE2.epoch_datetime
                    time_end = time()
                    print("Performance: Loop {:3f} Initialize {:.3f}  calc_fit {:.3f}  insert {:3f}".format(
                        time_end - time_loop_start,
                        time_to_calc_fit - time_loop_start,
                        time_calc_fit - time_to_calc_fit,
                        time_to_insert - time_calc_fit
                    ))

            else:
                log.error("DB Insert failed.")
                # main(db)
                result = False
                continue

# /////////////////// MAIN //////////////////////////////////////////////////////

def main(db=False):
    """ satfit
    Fit a TLE prediction to a reference TLE + new IOD observations

    Inputs:
        iod_line    Array of IOD formatted observations
        file_in     source file (contains IODs plus TLEs)
        file_out    result file (contains file_in + new TLE predict)
    """
    
    global srch
    global whichconst
    global iod_line
    global TLE_ref 

    TLE_ref_file = "/Users/chris/Dropbox/Docker/satobs/tle/bulk.tle"
    # TLE_ref_file = "/Volumes/astronomy/TLE/bulk.tle"
    # TLE_ref_file = ""

    log = logging.getLogger()

    # make it print to the console.
    console = logging.StreamHandler()
    log.addHandler(console)

    quiet = False
    verbose = 1
    if (quiet == False):
        if verbose == 0:
            log.setLevel(logging.WARN) 
        elif verbose == 1:
            log.setLevel(logging.INFO) 
        elif verbose == 2:
            log.setLevel(logging.DEBUG) 
        log.debug("Log level set to {}".format(log.level))

    # if verbose:
    #     for arg in vars(args):
    #         log.debug("%s : %s",arg, getattr(args, arg))

    # TODO: Expand this code to deal with command line options like the original
    file = "sat.txt"
    # TODO: Look into batch capability

    if (db): # Make this switch on file
        pass
    elif (not db):
        # Read reference TLEs (once) for seeding new TruSat catalog observations
        try:
            TLE_ref = TLEFile(TLE_ref_file)
        except FileNotFoundError:
            TLE_ref = None

        # Set up database connection

        # Temporary database credentials hack
        try:
            CONFIG = os.path.abspath("../trusat-config.yaml")
            db = database.Database(CONFIG)
        except: 
            log.error("DB Login credentials not available.")

    srch = 'W' # initialize wide search

    # iod_obs_id = input("IOD db ID: ")
    while(True):
        print("\nEnter command")
        print("Database: (I)OD search (O)bject search (L)atest (M)anual Process Latest q(U)eue")
        print("          ra(W) search Mc(C)ants TLE baseline (N)eeds TLE")
        print("(Q)uit")

        cmd = input(": ").strip()
        cmd = cmd.upper()

        if (cmd == "I"):    # IOD Search
            iod_search(db)
        elif (cmd == "O"):  # Object search
            object_search(db)
        elif (cmd == "M"):  # Manual update
            object_manual(db)
        elif (cmd == "C"):  # McCants TLE baseline
            object_tle_test(db)
        elif (cmd == "N"):  # Objects needing TLEs
            object_needs_TLE(db)
        elif (cmd == "W"):  # Raw search
            raw_search(db)
        elif (cmd == "L"):  # Latest observations
            pass
        elif (cmd == "U"):  # View the queue
            pass
        elif (cmd == "Q"):  # Quit to command line
            sys.exit(0)

        # print("{} using station {} Observing Satellite {} : {:.3f} days away from epoch".format(observer_name,station_number,satellite.model.satnum,days))


    #   // find observation lines in input file
    #   // number of observations = nobs (global)
    #   // store observations to iod_line matrix (up to 200 observations)
    # IOD_Records = iod.get_iod_records_from_file(file)
    # for item in IOD_Records:
    #     iod_line.append(item.line)

    # Variables:
    #  ll[nobs+1][3];       // line of sight vectors; observer -> satellite
    #  rd[nobs+1][3];       // topocentric vectors to observer positions
    #  odata[nobs+1][4];    // observational data: jd, ra, dc, obscode

    #   // storage for maneuver function
    #   int nobsx;
    #   char iod_linex[200][81];
    #   double llx[nobs][3];
    #   double rdx[nobs][3];
    #   double odatax[nobs][4];


    #   // read first tle from input file
    #   // orbit parameters are passed from read_tle as global variables
    #   // echo tle to screen
    # TODO: Replace this with a call to DB query
    # TLEs = TLEFile(file, strict=False)
    FitSat = Satellite_cal(line0=name, line1=tle1, line2=tle2)
    # TODO: Need to fix a problem in tle_util.parse_tles()s which only retains the most-recently read TLE for a sat

if __name__ == '__main__':
    main()
