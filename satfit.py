#!/usr/bin/env python

from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values
import sys

# As of 28 July 2019, python3.6 is the default "python3" in apt-get install python3
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version 3.6")
    sys.exit(1)

import configparser                 # config file parsing
import argparse                     # command line parsing
import os
from datetime import date, timedelta, datetime
from time import time                         # For performance timing
from math import (fabs, radians, sin, cos, pi, sqrt, fmod, acos, asin, atan, tan, degrees)    # Fast/precise math functions                      
import numpy as np
import logging
import string
import copy
from getpass import getpass # For getting input without newline

from spacetrack import SpaceTrackClient

# These are necessary until Brandon Rhodes approves pull requests
# https://github.com/brandon-rhodes/python-sgp4/pull/35
sys.path.insert(1, '/Users/chris/Dropbox/code/MVP/python-sgp4')
# https://github.com/skyfielders/python-skyfield/pull/276
sys.path.insert(2, '/Users/chris/Dropbox/code/MVP/python-skyfield')

# from skyfield.iokit import Loader, download, parse_tle
# from skyfield import sgp4lib

try:
    from sgp4.cpropagation import sgp4, sgp4init
except ImportError as e:
    print(e)
    from sgp4.propagation import sgp4, sgp4init
from sgp4.ext import invjday, days2mdhms
from sgp4.io import twoline2rv
from sgp4.model import Satellite
from sgp4 import earth_gravity



# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
tle_path = os.path.join(parentdir, "sathunt-tle")
sys.path.insert(1,tle_path) 
from tle_util import make_tle, append_tle_file, TLEFile, tle_fmt_epoch, datetime_from_tle_fmt, assumed_decimal_point, checksum_tle_line, myjday

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
iod_path = os.path.join(parentdir, "sathunt-iod")
sys.path.insert(1,iod_path) 
import iod

from elfind import read_obs, rref, SGN, so2r

# ///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////
twopi = 2*pi
nocon = twopi/1440.0

# TODO: Make a class of these?
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

# TODO: Figure out what these flags mean
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
# // declare variables for search and history options
# double e, ik, ok, ek, wk, mk, nk, bk, theta, minsum, ra, dc;

# double
#   astep,   amin,    amax,
#   istep,   imin,    imax,
#   ostep,   omin,    omax,
#   estep,   emin,    emax,    esize,
#   wstep,   wmin,    wmax,    wsize,
#   nstep,   nmin,    nmax,    nsize,
#   bstep,   bmin,    bmax,
#            mmax,    mmin;

estep = None
wstep = None
nstep = None

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
#         # /* Period > 225 minutes is deep space */
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
#         if (twopi / (xnodp * 1440.0) >= (1.0 / 6.4)):
#             sdp4(tsince)    # yes,  it should be a deep-space (SDP4) ephemeris
#         else:
#             sgp4(tsince)


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
    def __init__(self, time=None, year=None, month=None, day=None, hour=None, min=None, sec=None):
        self.time  = time               # Python datetime
        self.yy    = self.year  = year  # year
        self.mm    = self.month = month # month
        self.dd    = self.day   = day   # day of month, Greenwich
        self.hr    = self.hour  = hour  # hour, Greenwich
        self.mn    = self.min   = min   # minute
        self.ss    = self.sec   = sec   # seconds

        # Ten Output Values - above 5 and the below:
        self.thetag	= None  # Sidereal time in degrees
        self.jd		= None  # Julian date
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
                self.jd = myjday(self.yy, self.mm, self.dd, self.hr, self.mm, self.ss)
                self.sidereal()
            else:                   # this date is julian
                self.jd = self.time
                self.calcmjd()
                (self.yy, self.mm, self.dd, self.hr, self.mn, self.ss) = invjday(self.jd)
                self.sidereal()
                (_, subsec) = divmod(self.ss,1)
                subsec = int(subsec*1E6)
                intss = int(self.ss)
                self.time = datetime(self.yy,self.mm,self.dd,self.hr,self.mn,intss,subsec)
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
            self.jd = myjday(self.yy, self.mm, self.dd, self.hr, self.mn, self.ss)
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
        t = (self.jd - 2451545) / 36525
        thetag = 280.46061837 + 360.98564736629 * (self.jd - 2451545) \
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
        self.mjd = self.jd - 2400000.5

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

def mag(v):
    """ Computes the magnitude of a vector ||v|| 

    Renamed from norm(v) in original Scott Campbell code
    to better correspond to function names in SGP4 code.
    """
    mag = np.sqrt(np.dot(v, v))
    return mag

def unit_vector(v):
    """ Returns the unit vector of the vector.  """
    u = v / mag(v)
    return u

def delta_t(sat,t):
    tsince = (t - sat.jdsatepoch) * 1440.0 # time since epoch in minutes
    (rr, vv) = sgp4(sat,tsince,sat.whichconst) 

    sat.rr = np.asarray(rr) / sat.radiusearthkm # In Earth radii
    sat.vv = np.asarray(vv) / (sat.radiusearthkm / 60.0)  # In Earth radii / min - seriously!

    return sat

def delta_el(sat, xincl=False, xnodeo=False, eo=False, omegao=False, xmo=False, xno=False,   bsr=False,
                    ii=False,      om=False, ec=False,     ww=False,  ma=False,  nn=False, bstar=False):
    """ delta_el - Reinitialize SGP4 satellite object with changed Keplerian elements
    
    Units can be passed in via radians or degrees, as is the convention for variable names
    """
    if (xincl):
        sat.inclo = xincl
    elif (ii):
        sat.inclo = radians(ii)

    if (xnodeo):
        sat.nodeo = posradang(xnodeo)
    elif (om):
        sat.nodeo = posradang(radians(om))

    if (eo):
        sat.ecco = eo   
    elif (ec):
        sat.ecco = ec   

    if (omegao):
        sat.argpo = posradang(omegao)
    elif (ww):
        sat.argpo = posradang(radians(ww))

    if (xmo):
        sat.mo = posradang(xmo)
    elif (ma):
        sat.mo = posradang(radians(ma))

    if (xno):
        sat.no_kozai = xno
    elif (nn):
        sat.no_kozai = nn * nocon

    if (bsr):
        sat.bstar = bsr
    elif (bstar):
        sat.bstar = bsr
    sgp4init(sat.whichconst, sat.operationmode, sat.satnum, sat.jdsatepoch, sat.bstar, 0, 0, sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo, sat)
    return sat

def posradang(a):
    if a < 0:
        return a + twopi
    elif a > twopi:
        return a - twopi
    else:
        return a

def acose(x):
    if ( x >= 1):
        rval = 0
    elif ( x <= -1):
        rval = pi
    else:
        rval = acos(x)
    return rval

# // read site data, subroutine of read_obssf
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

def read_obssf(IOD_Records):
    """ decodes the iod_line data 
    
    This one is a more direct transcoding of Scott Campbell's, rather than the one in elfind
    which uses skyfield to get the julian data and topocentric position of the station
    """

    # for range in (0, nobs):
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

    ll    = np.zeros((nobs,3))
    odata = np.zeros((nobs,4))
    rd    = np.zeros((nobs,3))

    i = 0
    for iod_line in IOD_Records:
        # Grab the most recent version of these variables for writing eventual TLE
        # FIXME Scott's original code grabs the 2nd or "middle" datetime for epoch
        epoch_datetime = iod_line.DateTime
        ssn = iod_line.ObjectNumber

        (year, month, day, hour, minute, second) = iod_line.DateTime.timetuple()[:6]
        microseconds = int(iod_line.DateTime.strftime('%f'))
        sec_with_microseconds = second + microseconds/1.0E6

        # self.jdsatepoch = myjday(year, month, day, hour, minute, sec_with_microseconds)
        t1 = Date(year=year, month=month, day=day, hour=hour, min=minute, sec=sec_with_microseconds)

        # Skip the skyfield stuff while debugging against satfit.cpp
        # ts = load.timescale()
        # (year, month, day, hour, minute, second) = iod_line.DateTime.timetuple()[:6]
        # t_skyfield = ts.utc(year, month, day, hour, minute, second)
        # t1_jd = t_skyfield.tt

        ra = radians(iod_line.RA)
        dc = radians(iod_line.DEC)

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
        t = (t1.jd - 2451545) / 36525  
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

        odata[i][0] = t1.jd # julian date
        odata[i][1] = ra # ra radians (observed)
        odata[i][2] = dc # dc radians (observed)
        odata[i][3] = iod_line.Station # station

        # (la, lo, hh) = Sites.topos(iod_line.Station)
        # observer_location = Topos(latitude_degrees = la, longitude_degrees = lo, elevation_m = hh)

        # topocentric = observer_location.at(t_skyfield)
        # rd[i] = topocentric.position.km/_XKMPER   # elfind works in units of Earth radii
        read_site(iod_line.line)  # FIXME: This does a file open/scan for every call
        rre = Locate(t1, la, lo, hh)
        # Locate rre(t1.jd, la, lo, hh)
        rd[i] = rre.rre

        i+=1
    # end for
    return odata, ll, rd

def print_el(sat, deg=False, quiet=False):
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
    tle_epoch = tle_fmt_epoch(sat.epoch)
    eo_string = assumed_decimal_point(sat.ecco,7)

    if (deg==False):
        xincl  = degrees(sat.inclo)
        xnodeo = degrees(sat.nodeo)
        omegao = degrees(sat.argpo)
        xmo    = degrees(sat.mo)
        xno    = degrees(sat.no_kozai) * 1440/360 # revs per day

    # if name:
    #     line0  = "{:24s}".format(name) 
    #     if not quiet:
    #         print("{:s}".format(line0))

    line1 = "1 {:5d}U {:<8s} {:14s} 0.00000073  00000-0  50000-4 0    00".format(sat.satnum,sat.intldesg,tle_epoch)
    # TODO: Deal with First Derivative xno, Second derivative xno, Bstar
    line2 = "2 {:05d} {:8.4f} {:8.4f} {:7s} {:8.4f} {:8.4f} {:11.8f}    00".format(
            sat.satnum, xincl, xnodeo, eo_string, omegao, xmo, xno)

    line1 = line1[:68] + str(checksum_tle_line(line1))
    line2 = line2[:68] + str(checksum_tle_line(line2))

    if not quiet:
        print()
        print("{:s}".format(line1))
        print("{:s}".format(line2))
    return(line0, line1, line2)


def longitude(ma, ww, ec):
    """Calculate true longitude from mean anomaly and argument of perigee
    Inputs: 
     ma     mean anomaly, degrees
     ww     argument of perigee, degrees

    Outputs:
     uu     True Longitude, degrees
    """

    ma = radians(ma)
    e = ma
    ek = 0
    while(fabs(e - ek) > 1e-6):
        ek = e
        e = ma + ec * sin(e)
    ma = degrees(ma) # FIXME: This appears to be here because of its prior global variable status
    theta = (ec - cos(e)) / (ec * cos(e) - 1)
    theta = acose(theta)
    if (e > pi):
        theta = 2 * pi - theta
    uu = ww + degrees(theta)
    return uu

def sort(iod_line, odata, ll, rd):
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
    # global nobs # FIXME shouldn't have to do this if we're not changing it

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
    return [iod_line, odata, ll, rd]

# TODO: This function (and those it calls) will benefit the best from accelerating
def find_rms(satx, rd, ll, odata): ## WORKS
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    # copy sat
    # TODO: evaluate if this is best done in main()
    # satx = copy.deepcopy(sat)
    nobs = len(odata)
    zum = 0
    for j in range(nobs):
        # advance satellite position
        satx = delta_t(satx,odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # delr = satx.rr - rd[j] # This was vmadd -1
        # nrr = mag(delr)       # range

        # convert to unit vector
        # delr = delr/nrr
        delr = unit_vector(satx.rr - rd[j])
        # topocentric position error in degrees
        Perr = degrees( acose( np.dot(delr, ll[j]) ) )

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)


def print_fit(sat, rd, ll, odata, last_rms): # WORKS
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

    nobs = len(odata)
    if (nobs == 0):
        print("\nno obs")
        return

    # copy sat
    satx = copy.deepcopy(sat)

    fit_string = "\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr"

    print(fit_string)

    # FIXME this (global?) variable
    if (out):
        with open(file, "a") as fp:
            fp.write("\n")
            fp.write(fit_string)

    for j in range (0, nobs):
        # advance satellite position
        delta_t(satx,odata[j][0])

        nrr = mag(satx.rr)
        nvv = mag(satx.vv)

        # computing elevation  # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        el = degrees( acose(np.dot(rd[j], ll[j]) / mag(rd[j])) )
        el = 90 - el
        el = round(el, 1)

        # computing aspect
        asp = degrees( acose(np.dot(ll[j], satx.vv) / nvv) )
        asp = 180 - asp
        asp = round(asp, 1)

        # computing azimuth # TODO: This doesn't quite jive with astropy conversion, however it looks like its just for display
        tempv = np.cross(rd[j], zz)
        nv = np.cross(tempv, rd[j])
        tempv = np.cross(rd[j], ll[j])
        temp = np.cross(tempv, rd[j])
        az = acose(np.dot(nv, temp) / (mag(nv)*mag(temp)))
        tempv = np.cross(temp, nv)
        if (np.dot(tempv, rd[j]) < 0):
            az = 2*pi - az
        if (mag(temp) == 0):
            az = 0.0
        az = degrees(az)
        az = round(az, 1)

        # observed satellite geocentric position vector, rr
        rr = so2r(nrr, rd[j], ll[j])

        # geocentric position error angle in radians
        Perr = acose(np.dot(satx.rr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, delr, in e.r.
        delr = satx.rr - rr
        temp = np.cross(satx.rr, satx.vv)  # xtrk reference vector points left of track
        sign = SGN(np.dot(delr, temp))

        # observer velocity vector
        tempv = np.cross(zz, rd[j])
        temp = .004351409367 * tempv 
        # observed satellite velocity
        tempv = satx.vv - temp
        nvv = mag(tempv)

        # angle between delr vector and tempv vector, radians
        alpha = acose(np.dot(tempv, delr) / (nvv * mag(delr)))

        # magnitude of delr in direction of tempv, radians
        delt = atan(cos(alpha) * tan(Perr))   # geocentric range error

        # time error
        delt *= nrr / nvv                     # delta r in min
        delt *= 60                            # seconds
        delt  = round(delt, 2)

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # new use of delr variable, predicted line of sight vector
        delr = satx.rr - rd[j]
        nrr = mag(delr)

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

        tsince = (odata[j][0] - satx.jdsatepoch)  * 1440.0 # time since epoch in minutes # TODO: Clean up this calculation
        obstime = satx.epoch + timedelta(minutes=tsince)
        timestring = obstime.strftime('%y%j %H%M:%S')
        SSS = obstime.strftime('%f')
        SSS = int(1000*(int(SSS)/1E6))
        fit_string = "({:2d}) {:04d}  {}{:03d}  {:5.1f}  {:5.1f}  {:5.1f}  {:6.2f}   {:6.2f}  {:7.3f}".format(
            j + 1, int(odata[j][3]), timestring, SSS, az, el, asp, xtrk, delt, Perr)
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

# ////////////////// FUNCTIONS //////////////////////////////////////////////////

def asym(a1, a2, a3):
    """ asymptotic extrapolation """
    if (fabs(a1 - a2) < 1.0e-5):
        b = 0.0
    else:
        b = (a2 - a3) / (a1 - a2)
    return (b * a2 - a3) / (b - 1.0)


def rtw(ao, ac):
    """ round the world """
    if (fabs(ao - ac) > 180):
        if (ao < 180):
            ao += 360
        else:
            ac += 360
    return ao - ac;


def ww2ma(wx):
    """ find mean anomaly from true longitude and perigee """
    # FIXME uses uu, ec globals
    theta = radians(fmod(uu - wx, 360))
    e = acose((ec + cos(theta)) / (1 + ec * cos(theta)))
    if (theta > pi):
        e = 2 * pi - e
    ma = e - ec * sin(e)
    return degrees(ma)

def zrll(satx, rd):
    """ Calculates predicted direction angles, right ascension(ra) and
    declination(dc), of the line of sight vector, rll, in equatorial coordinates
    from an observer position, rd, to the satellite position, sat.rr.
    A subroutine of diff_el.  Satellite object, sat, simplifies prediction

    Inputs:
        satx    perturbed Satellite object at TLE
        rd      input, topocentric vector to observer
    Outputs:
        ra      calculated right ascension, degrees
        dc      calculated declination, degrees

    """
    rll = satx.rr - rd # line of sign vector
    # return ra and dc, degrees
    dc = degrees(asin(rll[2] / mag(rll)))
    ra = degrees(atan(rll[1] / rll[0]))
    if (rll[0] < 0):
        ra += 180.0
    ra = fmod(ra, 360)
    return (ra, dc)


def diff_el(sat, rd, ll, odata, sum):
    """ differential correction of elements

    Variables:
        int i, j, k, c = 0;
        ra, dc;                 // right ascension, declination variables
        delta, el;              // small change amount
        mdata[150][7];          // define matrix for multiple regression
        dat[2][7];              // define output matrix of zpde utility
        b[6];                   // output deltas, sat position, velocity, b*
        rv[6][7];               // normal matrix
        rac, dcc, rms;


    """

    c = 0 
    dat = np.zeros((2,7))
    # mdata = np.zeros((150,7))
    # mdata = np.empty((0,7), int)
    mdata = []
    rv = np.zeros((6,7))
    b = np.zeros(6)
    nobs = len(odata)

    #loop:
    satx = copy.deepcopy(sat)
    while(True): # Forever loop (at least for 20 improvements)
        # begin the main loop of the program
        for i in range(nobs):
            satz = copy.deepcopy(sat)                      # differential sat

            delta_t(satx,odata[i][0])

            # first establish the computed ra, dc, at jdo with no perturbations
            (ra, dc) = zrll(satx, rd[i])     # output ra, dc, degrees
            rac = ra                        # store computed ra and dc
            dcc = dc

            # find the deltas and load into output matrix, dat
            dat[0][6] = rtw(degrees(odata[i][1]), rac)      # store delta_ra
            dat[1][6] = degrees(odata[i][2]) - dcc          # store delta_dc

            # 6 steps to build differential correction matrix
            j = 0
            if (xi):
                dat[0][j] = .001
                dat[1][j] = .001
            else:
                delta = 0.001                        # change
                el = copy.copy(satz.inclo)           # store reference
                # satz.inclo += delta                  # delta element
                satz = delta_el(satz,xincl=(satz.inclo + delta))
                satz = delta_t(satz, odata[i][0])    # recalculate with perturbed element # FIXME python-SGP4
                (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
                # satz.inclo = el                      # restore reference
                dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 1
            delta = 0.001                        # change
            el = copy.copy(satz.nodeo)          # store reference
            # satz.nodeo += delta                 # delta element
            satz = copy.deepcopy(sat)                      # differential sat
            satz = delta_el(satz,xnodeo=(satz.nodeo + delta))
            satz = delta_t(satz, odata[i][0])            # recalculate with perturbed element # FIXME: python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            # satz.nodeo = el                     # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

            # Results from this one are fairly different - dat[0][j] -453 vs -474
            j = 2
            if (xe):
                dat[0][j] = 0.00001
                dat[1][j] = 0.00001
            else:
                delta = 0.0001                       # change
                el = satz.ecco                         # store reference
                # satz.ecco += delta                     # delta element
                satz = copy.deepcopy(sat)                      # differential sat
                satz = delta_el(satz,eo=(satz.ecco + delta))
                satz = delta_t(satz, odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
                (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
                # satz.ecco = el                         # restore reference
                dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 3
            if (xw):
                dat[0][j] = 0.001
                dat[1][j] = 0.001
            else:
                delta = 0.001                        # change
                el = satz.argpo                     # store reference
                # satz.argpo += delta                 # delta element
                satz = copy.deepcopy(sat)                      # differential sat
                satz = delta_el(satz,omegao=(satz.argpo + delta))
                satz = delta_t(satz,odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
                (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
                # satz.argpo = el                     # restore reference
                dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 4
            delta = 0.001                        # change
            el = satz.mo                        # store reference
            # satz.mo += delta                    # delta element
            satz = copy.deepcopy(sat)                      # differential sat
            satz = delta_el(satz,xmo=(satz.mo + delta))
            satz = delta_t(satz,odata[i][0])            # recalculate with perturbed element # FIXME: python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            # satz.mo = el                        # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

            # TODO: Investigate whether we should perterb the no_unkozai or no_kozai (doesn't seem to make much difference)
            j = 5
            if (xn):
                dat[0][j] = 0.000001
                dat[1][j] = 0.000001
            else:
                delta = 0.00001                      # change
                el = satz.no_kozai                        # store reference
                # satz.no_kozai += delta                    # delta element
                satz = copy.deepcopy(sat)                      # differential sat
                satz = delta_el(satz,xno=(satz.no_kozai + delta))
                satz = delta_t(satz,odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
                (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
                # satz.no_kozai = el                        # restore reference
                dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            # mdata[2 * i]     = dat[0]   # numerical deltas transferred to
            # mdata[2 * i + 1] = dat[1]   # multiple regresssion matrix
            mdata.append(dat[0])
            mdata.append(dat[1])
            # END for i


        # multiple regression
        for j in range(6):
            for k in range (7):
                rv[j][k] = 0
                for i in range (nobs*2):
                    rv[j][k] = rv[j][k] + mdata[i][j] * mdata[i][k]
        b = rref(rv) # Returns false if singular matrix

        # Getting inf and -inf on some results
        if (np.isinf(b).any()):
            break

        saty = copy.deepcopy(sat) 
        # test update components with deltas
        saty.inclo    += b[0]*0.1
        saty.nodeo    += b[1]*0.1
        saty.ecco     += b[2]*0.1
        saty.argpo    += b[3]*0.1
        saty.mo       += b[4]*0.1
        saty.no_kozai += b[5]*0.1
        saty = delta_el(saty)
        rms = find_rms(saty, rd, ll, odata)
        if (rms < sum):
            sum = rms
            # update components with deltas
            sat.inclo    += b[0]*0.1
            sat.nodeo    += b[1]*0.1
            sat.ecco     += b[2]*0.1
            sat.argpo    += b[3]*0.1
            sat.mo       += b[4]*0.1
            sat.no_kozai += b[5]*0.1
            c+=1
            if (c < 20):
                continue # Back up to the top of the loop
            else:
                break
        else:
            break

    sat = delta_el(sat)
    print("                  {} diff_el iter".format(c))

    return sat
    # /*
    # // display computed deltas for all 6 components and an rms
    # for(i = 0; i < 6; i++)
    # {
    #     printf("\n");
    #     printf("   %9.6f", b[i]);
    # }
    # printf("\n\nrms%9.6f\n", rms);
    # s_in(": ", buf);
    # */
    # global elements updated
    # ii = degrees(sat.xincl)
    # om = degrees(sat.xnodeo)
    # ec = sat.eo
    # ww = degrees(sat.omegao)
    # ma = degrees(sat.xmo)
    # nn = sat.xno/nocon

def anomaly_search(sat, rd, ll, odata, sum):
    """ mean anomaly box search, no limits """
    # global ma   # Not supposed to need these unless we're assigning, but alas...
    # global sum    # Some loops (i.e. perigee_search) want to know the most recent sum to start
    # global nsum
    # global xsum

    step = 0.1
    # mk = ma # FIXME global
    mk = degrees(sat.mo)
    # sum = find_rms(sat, rd, ll, odata) # PRobably not needed

    max = 1 # CPP do loop evaluates while at the end
    min = 0
    while(fabs(max - min) > 1e-5):
        min = mk
        max = mk

        # nsum loop - until rms doesn't decrease since last loop
        while(True):
            min = mk - step
            sat = delta_el(sat, ma=min)
            # sat.delta_el(sat.jd, ii, om, ec, ww, min, nn, bstar) # FIXME python-SGP4
            nsum = find_rms(sat, rd, ll, odata)
            if (nsum < sum):
                mk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop - until rms doesn't decrease since last loop
        while(True): 
            max = mk + step
            sat = delta_el(sat, ma=max)
            # sat.delta_el(sat.jd, ii, om, ec, ww, max, nn, bstar)    # FIXME python-SGP4
            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                mk = max
                sum = xsum
                continue # Back to top of xsum loop
            break   
        step /= 2
    sat = delta_el(sat, ma=mk)
    # ma = fmod(mk, 360)
    return sat # Contains the ma at the end of the loop

def motion_search(sat, rd, ll, odata):
    """ mean motion box search, no limits """
    # nk = nn
    nk = sat.no_kozai/nocon
    sum = find_rms(sat, rd, ll, odata)

    # Start with this values to get through the loop once
    # CPP has the while evaluation at the end
    min = 0
    max = 1
    step = 0.1
    while(fabs(max - min) > 1e-10):
        min = nk
        max = nk

        # nsum loop - until rms doesn't decrease since last loop
        while(True):
            min = nk - step
            sat = delta_el(sat, nn=min)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, min, bstar)    # FIXME python-SGP4
            nsum = find_rms(sat, rd, ll, odata)
            if (nsum < sum):
                nk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop - until rms doesn't decrease since last loop
        while(True):
            max = nk + step
            sat = delta_el(sat, nn=max)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, max, bstar)    # FIXME python-SGP4
            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                nk = max
                sum = xsum
                continue
            break
        step /= 2
    nn = nk
    return sat # nn (mean motion) is incorporated in the last update for the sat variable.

def node_search(satx, rd, ll, odata, sum, imax, imin, omax, omin):
    """ partition search on node and inclination within set limits """
    global xi   # FIXME - Hold ii constant to user-specified value?  Gets set in incl()

    ii = degrees(satx.inclo)
    om = degrees(satx.nodeo)

    # satx = copy.deepcopy(sat)   # Working with the copy might not be necessary
    while((imax - imin) > 1e-5):
        istep = (imax - imin) / 20
        ostep = fabs(rtw(omax, omin) / 20)

        if (xi): # FIXME: to have the effect of running the for and while loops once for the fixed imin value
            imin  = ii
            imax  = ii + 1e-6
            istep = 10*imin
        for ik in np.arange(imin, imax, istep):
            for ok in np.arange(omin, omax, ostep):
                delta_el(satx, ii=ik, om=ok)
                # satx(tle, ik, ok, ec, ww, ma, nn, bstar); # FIXME SGP4

                # establish the computed ra, dc, at jdo with no perturbations
                rms = find_rms(satx, rd, ll, odata)
                if (rms < sum):
                    sum = rms
                    ii  = ik
                    om  = ok
            # END for ok
        # END for ik
        imin = ii - istep
        imax = ii + istep
        omin = om - ostep
        omax = om + ostep
    # om = fmod(om, 360)
    satx = delta_el(satx, ii=ii, om=om)
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4
    return satx


def perigee_search(sat, rd, ll, odata, sum, uu, wmax, wmin, emax, emin):
    """ partition search on perigee and eccentricity 
    
    Global variables used:
        uu      longitude
        xe      flag
        xn      flag
        xw      flag
    """
    # global xe
    # global xn
    # global xw
    # global nsum
    # global xsum
    # global estep
    # global wstep

    xe = 0
    xn = 0
    xw = 0

    # satx = copy.deepcopy(sat)

    # Grab the values we're searching for in the loop, in case we don't find a new optimal
    ec  = sat.ecco
    ww  = degrees(sat.argpo)
    ma  = degrees(sat.mo)

    if (sat.ecco > 0.1):
        wmax = degrees(sat.argpo) + 0.1
        wmin = degrees(sat.argpo) - 0.1
        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99

    # rms_time = 0
    # dl_time = 0
    # lc = 0
    while((wmax - wmin) > 1e-5):
        estep = (emax - emin) / 20
        wstep = (wmax - wmin) / 20
        for wk in np.arange(wmin, wmax, wstep):
            if (xw):
                wmin  = degrees(sat.argpo)
                wmax  = degrees(sat.argpo)
                wk    = degrees(sat.argpo)
                wstep = 0
            theta = radians(uu - wk)    # FIXME global (uu)
            for ek in np.arange(emin, emax, estep):
                if (xe):
                    emin  = sat.ecco
                    emax  = sat.ecco
                    ek    = sat.ecco
                    estep = 0
                e = acose((ek + cos(theta)) / (1 + ek * cos(theta)))
                if (theta > pi):
                    e = 2 * pi - e
                mk = e - ek * sin(e)
                mk = degrees(mk)
                # satx = copy.deepcopy(sat)
                # lc+=1
                # t_start = time()
                sat = delta_el(sat, ec=ek, ww=wk, ma=mk)
                # t_end = time()
                # dl_time += (t_end-t_start)

                # satx(sat.jd, ii, om, ek, wk, mk, nn, bstar) # FIXME python-SGP4
                # establish the computed ra, dc, at jdo with no perturbations
                # t_start = time()
                rms = find_rms(sat, rd, ll, odata)
                # t_end = time()
                # rms_time += (t_end-t_start)
                if (rms < sum):
                    sum = rms
                    ec  = ek
                    ww  = wk
                    ma  = mk
            # END for ek
        # END for wk

        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4
        # Could save a call here by checking for the existence of the variables, but what's one more time?
        sat = delta_el(sat, ec=ec, ww=ww, ma=ma)

        wmax = degrees(sat.argpo) + wstep
        wmin = degrees(sat.argpo) - wstep
        emax = sat.ecco + estep
        emin = sat.ecco - estep
               
    # print()
    # # print("PS: {} loop iterations".format(lc))
    # print("                  PS: find_rms: {:.2f}".format(rms_time))
    # print("                  PS: delta_el: {:.2f}".format(dl_time))

    # update mean_anomaly
    # 0.01685786247253418
    sat = anomaly_search(sat, rd, ll, odata, sum)
    # sat = delta_el(sat, xmo=radians(ma)) # Not needed from the last step of above
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4

    # update mean_motion
    if (not xn):
        # motion_search: 0.030543804168701172
        sat = motion_search(sat, rd, ll, odata)
        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4

    # calculate uu, degrees
    uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)

    # ww = fmod(degrees(sat.argpo), 360)
    # ma = fmod(degrees(sat.mo), 360)
    return sat

def align(sat, rd, ll, odata):
    """ sets deltaT error at last obs epoch to zero """
    nobs = len(odata)
    last = nobs - 1 # FIXME global
    satx = sat        # copy sat

    delt = 1    # Give delt a starting value for the loop
    while(fabs(delt) > 1.0e-5): # FIXME global
        # advance satellite position
        delta_t(satx,odata[last][0])
        nrr = mag(satx.rr)
        nvv = mag(satx.vv)

        # observed geocentric position vector, rr
        rr = so2r(nrr, rd[last], ll[last])

        # position error in radians
        Perr = acose(np.dot(satx.rr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, er
        delr = np.subtract(satx.rr, rr)

        # magnitude of delta r in direction of v, radians
        delt = Perr * np.dot(satx.vv, delr) / (nvv * mag(delr))

        # degrees
        delt = degrees(delt)

        # TODO not sure where first gets set
        if (first):
            first = 0;
            ma = ma - 0.75*delt
        else:
            ma = ma - delt/5.0

        satx.delta_el(satx.jd, ii, om, ec, ww, ma, nn, bstar)   # FIXME python-SGP4

def fit_out(sat, rd, ll, odata):
    out = 1
    sum = print_fit(sat, rd, ll, odata, sum)
    out = 0
    # accept_command()

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

def remove():
    # FIXME Make the global variables local

    buf = input("\nRemove Line Number : ")
    j = int(buf.strip())

    for i in range(j, nobs):
        rd[i - 1] = rd[i]
        ll[i - 1] = ll[i]
        odata[i - 1] = odata[i]
        iod_line[i - 1] = iod_line[i]
        nobs-=1

def fit(sat, rd, ll, odata, sum):
    out = 0
    sum = print_fit(sat, rd, ll, odata, sum)
    print_el(sat)
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


def accept_command(sat, rd, ll, odata, sum, uu):
    """Accept a new command"""

    while(True):    # Forever loop
        print("\nEnter command")
        print(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee   (A)nomaly  (M)otion  (B)star")
        print(" (S)tep   I(D)   (T)ime  (F)it  (V)iew  (R)emove   (W)rite   (Q)uit")

        cmd = input(": ").strip()
        cmd = cmd.upper()

        # Hidden functions
        if (cmd == "G"):    # Graph fit
            fit_out(sat, rd, ll, odata)
        elif (cmd == "E"):  # Edit File
            edit_file()
        elif (cmd == "H"):  # History
            sat = history(sat, rd, ll, odata)
        elif (cmd == "Y"):  # Edit History
            edit_history()
        elif (cmd == "C"):  # Elcor
            elcor(sat, rd, ll, odata)
        elif (cmd == "O"):  # Discover
            discover(sat, rd, ll, odata)
        elif (cmd == "U"):  # Maneuver
            sat = maneuver(sat, rd, ll, odata)

        # Visible functions
        elif (cmd == "S"):  # Step
            sat = step(sat, rd, ll, odata, sum, uu, "S")
        elif (cmd == "Z"):
            sat = step(sat, rd, ll, odata, sum, uu, "Z")
        elif (cmd == "I"):  # Incl
            print_el(sat)
            sat = incl(sat, rd, ll, odata, sum)
        elif (cmd == "N"):  # Node
            print_el(sat)
            sat = node(sat, rd, ll, odata, sum)
        elif (cmd == "X"):  # Eccentricity
            print_el(sat)
            sat = xntrcty(sat, rd, ll, odata, sum)
        elif (cmd == "P"):  # Perigee
            print_el(sat)
            sat = perigee(sat, rd, ll, odata, sum, uu)
        elif (cmd == "A"):  # Mean Anomaly
            print_el(sat)
            sat = anomaly(sat, rd, ll, odata, sum, uu)
        elif (cmd == "M"):  # Mean Motion
            print_el(sat)
            sat = motion(sat, rd, ll, odata, sum)
        elif (cmd == "B"):  # Bstar
            print_el(sat)
            sat = bstar(sat, rd, ll, odata, sum)

        elif (cmd == "D"):  # ID
            id()
        elif (cmd == "T"):  # Time
            time_func()
        elif (cmd == "F"):  # Fit
            sum = fit(sat, rd, ll, odata, sum)
        elif (cmd == "V"):  # View observations
            viewobs()
        elif (cmd == "R"):  # Remove observations
            remove()
        elif (cmd == "W"):  # Write elements
            write_el(sat)
        elif (cmd == "Q"):  # Quit
            sys.exit(0)
    
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
    nn = sat.no_kozai/nocon

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
    satx = diff_el(sat,ec=ec,ww=ww,ma=ma,nn=nn)
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # update mean_anomaly
    sat = anomaly_search(sat, rd, ll, odata, sum)
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # print new elements
    sum = print_fit(sat, rd, ll, odata, sum)
    print_el(sat)

    uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)

    srch = 'Z' # FIXME, this is probably a global
    return sat

def step(sat, rd, ll, odata, sum, uu, step_type):       # partition search within limits set below

    global srch

    nobs = len(odata)
    last_rms = sum
    # update mean_anomaly
    # anomaly_search: 0.015267133712768555
    sat = anomaly_search(sat, rd, ll, odata, sum)

    # sat = delta_el(sat, xmo=radians(ma)) # Not needed, as latest ma returned from above
    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
 
    # emax = ec * 1.1
    # emin = ec * 0.9
    # wmax = ww + 2
    # wmin = ww - 2
    # imax = ii + 2
    # imin = ii - 2
    # omax = om + 2
    # omin = om - 2

    emax = sat.ecco * 1.1
    emin = sat.ecco * 0.9
    wmax = degrees(sat.argpo) + 2
    wmin = degrees(sat.argpo) - 2
    imax = degrees(sat.inclo) + 2
    imin = degrees(sat.inclo) - 2
    omax = degrees(sat.nodeo) + 2
    omin = degrees(sat.nodeo) - 2

    print("\nPress Q to Quit    :\n")
    xsum = 0
    lc = 0
    stp_start = time()
    DE = 0
    while(fabs(sum-xsum)>1e-4):
        lc+=1
        xsum = sum
        ps_start = time()
        sat = perigee_search(sat, rd, ll, odata, sum, uu, wmax, wmin, emax, emin)
        # sat = delta_el(sat, omegao=radians(ww)) # Not needed as ww is returned in sat variable from above
        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        ns_start = time()
        sat = node_search(sat, rd, ll, odata, sum, imax, imin, omax, omin)
        # sat = delta_el(sat, xnodeo=radians(om))
        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

        de_start = time()
        if (nobs > 3 and step_type == 'Z'):
            sat = diff_el(sat, rd, ll, odata, sum)
            # sat = delta_el(sat)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            de_stop = time()
            DE = de_stop - de_start

        # set new limits
        # FIXME: Figure out globals
        # emax = 1.01 * ec
        # emin = 0.99 * ec
        # wmax = ww + 0.5
        # wmin = ww - 0.5
        # imax = ii + 0.5
        # imin = ii - 0.5
        # omax = om + 0.5
        # omin = om - 0.5

        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99
        wmax = degrees(sat.argpo) + 0.5
        wmin = degrees(sat.argpo) - 0.5
        imax = degrees(sat.inclo) + 0.5
        imin = degrees(sat.inclo) - 0.5
        omax = degrees(sat.nodeo) + 0.5
        omin = degrees(sat.nodeo) - 0.5

        # Shouldn't need this
        sum = find_rms(sat, rd, ll, odata)

        stp_lap = time()
        lap = stp_lap - ps_start
        PS = ns_start - ps_start
        NS = de_start - ns_start
        ELAPSED = stp_lap - stp_start
        print("rms{:12.5f}   Lap time: {:.2f}  PS {:.2f}  NS {:.2f}  DE {:.2f}   --  Elapsed {:.1f} / {}".format(sum, lap, PS, NS, DE, ELAPSED, lc))
        # buf = input("    : ")

        # try:
        #     buf = buf.strip().upper()
        #     if (buf == 'Q'):
        #         break
        # except NameError:
        #     continue

    print()
    # print("Step: {} loop iterations".format(lc))


    sum = print_fit(sat, rd, ll, odata, last_rms)
    print_el(sat)       # print new elements

    srch = 'N' # FIXME: Fix this global
    return sat


def incl(sat, rd, ll, odata, sum):
    global srch

    ii = degrees(sat.inclo)
    om = degrees(sat.nodeo)

    while(True): # Forever loop
        print("\n(A)uto  (ii)  (Q)uit  ")
        buf = input(": ").strip()

        try: 
            ii = float(buf)
            xi = 1
            sat = delta_el(sat, ii=ii)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
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

                sat = node_search(sat, rd, ll, odata, sum, imax, imin, omax, omin)

                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

                srch = 'N' # FIXME: Fix this global
            elif (buf == 'Q'):
                return sat

def node(sat, rd, ll, odata, sum):
    global srch

    while(True): # Froever loop
        om = degrees(sat.nodeo)

        print("\n(A)uto  (om)  (Q)uit  ")
        buf = input(": ").strip()

        try:
            om = float(buf)
            sat = delta_el(sat,om=om)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
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

                sat = delta_el(satx,om=om)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

                srch = 'N' # FIXME: Fix this global variable
            elif (buf == 'Q'):
                return sat


def xntrcty(sat, rd, ll, odata, sum):
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
            sat = delta_el(sat,ec=ec)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
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
                        sat = delta_el(sat,ec=ek)
                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(sat, rd, ll, odata)
                        if (rms < sum):
                            sum = rms
                            ec = ek
                        # end for ek
                    sat = delta_el(sat,ec=ec)
                    emin = ec - estep
                    emax = ec + estep

                sat = delta_el(sat,ec=ec)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

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
                    sat = delta_el(sat,ec=ec)
                    sum = find_rms(sat, rd, ll, odata)
                    print("\n{:.7f}     {:7.4f}", ec, sum)

                print()
                ec = ek
                sat = delta_el(sat,ec=ec)

            elif (buf == 'Q'):
                return sat

def perigee(sat, rd, ll, odata, sum, uu):
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
            uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)
            ma = ww2ma(ww)
            xw = 1
            sat.delta_el(sat, ww=ww, ma=ma)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
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
                    sat = perigee_search(sat, rd, ll, odata, sum, uu, wmax, wmin, emax, emin)
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

                srch = 'N' # FIXME Global

            # SEARCH
            elif (buf == 'S'):
                uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)
                
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


def anomaly(sat, rd, ll, odata, sum, uu):
    global srch

    while(True): # Forever loop
        if ((sat.no_kozai/nocon) < 1.5):
            # amax and amin are used as temporary variables
            uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)
            amax = radians(uu)
            amin = sin(degrees(sat.xincl)) * sin(amax)
            amin *= amin
            amin = sqrt(1 - amin)
            amin = degrees((acose(cos(amax) / amin)))
            if (fmod(uu, 360) > 180):
                amin = 360 - amin
            if(degrees(sat.xincl) < 90):
                amax = fmod(degrees(sat.xnodeo) + amin - sat.thetag + 360, 360.0)
            else:
                amax = fmod(degrees(sat.xnodeo) - amin - sat.thetag + 720, 360.0)
            print("\nE Longitude = {:8.4f}".format(amax))

        print("\n(S)earch  (A)uto  (ma)  (L)ast  (Q)uit  ")

        amax = 360
        amin = 0.0
        astep = 20

        buf = input(": ").strip()

        try:
            ma = float(buf)
            uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)
            sat = delta_el(sat, ma=ma)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                sat = anomaly_search(sat, rd, ll, odata, sum)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements
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
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                    sat = delta_el(sat,ma=ma)
                    sum = find_rms(sat, rd, ll, odata)
                    print("\n{:8.4f}     {:7.4f}".format(ma, sum))
                print()
                ma = mk              # restore
                sat = delta_el(sat,ma=ma)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

            # LAST
            elif (buf == 'L'):
                sat = align(sat, rd, ll, odata)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

            # QUIT
            elif (buf == 'Q'):
                return sat


def motion(sat, rd, ll, odata, sum):
    global xn

    while(True): # Forever loop
        print("\n(A)uto  (nn)  (Q)uit  ")
        buf = input(": ").strip()

        try:
            nn = float(buf)
            xn = 1
            sat = delta_el(sat, nn=nn)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                xn = 0
                # update mean motion, no limits
                sat = motion_search(sat, rd, ll, odata)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements
            elif (buf == 'Q'):
                return sat


def bstar(sat, rd, ll, odata, sum):
    while(True): # Forever loop
        print("\n(A)uto  (b*)  (B)atch  (Q)uit  ")
        buf = input(": ").strip()

        bstar = sat.bstar

        try:
            bstar = float(buf)
            sat = delta_el(sat,bstar=bstar)
            # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)       # print new elements
        except:
            buf = buf.upper()
        
            # AUTO
            if (buf == 'A'):
                # update Bstar within limits
                bmax = bstar * 1.1
                bmin = bstar * 0.9
                if (bstar < 0):
                    bmax = bmin
                    bmin = bstar * 1.1

                ## This doesn't seem very "auto"
                # TODO: make "Search" version
                # print("\nbmax [{:.8f}".format(bmax))
                # buf = input("]: ").strip()
                # bmax = float(buf)

                # print("\nbmin [{:.8f}".format(bmin))
                # buf = input("]: ").strip()
                # bmin = float(buf)

                while((bmax - bmin) > 1.e-9):
                    bstep = (bmax - bmin) / 20
                    for bk in np.arange(bmin, bmax, bstep):
                        sat = delta_el(sat,bstar=bk)
                        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bk)
                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(sat, rd, ll, odata)
                        if (rms < sum):
                            sum = rms
                            bstar = bk
                       # END for bk
                    # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                    bmin = bstar - bstep
                    bmax = bstar + bstep

                sat = delta_el(sat,bstar=bstar)
                # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = print_fit(sat, rd, ll, odata, sum)
                print_el(sat)       # print new elements

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
                return sat


def maneuver(sat, rd, ll, odata, sum):
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
            llx   = copy.copy(ll)
            rdx   = copy.copy(rd)
            odatax = copy.copy(odata)

            for i in range(p, nobs):
                iod_line[i-p] = iod_line[i-1]
                ll[i-p] = ll[i-1]
                rd[i-p] = rd[i-1]
                odata[i-p] = odata[i-1]
            nobsx = nobs
            nobs  = nobs - p + 1

            out = 0
            sum = print_fit(sat, rd, ll, odata, sum)
            print_el(sat)
            print("\nperiod = %f days".format(nocon/sat.xno))
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
                    if (ma < 0.1):
                        time = satm.jd - nocon/satm.xno*(1 + satm.xmo/twopi)
                    # if not previously at perigee, back up to perigee
                    else:
                        time = satm.jd - satm.xmo*nocon/(satm.xno*twopi)
                # NEXT: advance one revolution
                if (buf == 'N'):
                    # if previously at perigee, go forward one revolution
                    if (ma > 359.9):
                        time = satm.jd + nocon/satm.xno*(2 - satm.xmo/twopi)
                    # if not previously at perigee, go forward to perigee
                    else:
                        time = satm.jd + nocon/satm.xno*(1 - satm.xmo/twopi)
                # move to time and ma at perigee
                t2 = Date(time=time)
                delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                # satm.delta_t(t2.jd) # FIXME: Replace with python-SGP4 call
                # satm.rv2el(satm.rr, satm.vv) # FIXME: Replace with python-SGP4 call
                satm.jd = t2.jd
                ma = fmod(degrees(satm.xmo), 360)
                # refine perigee
                for i in range(0, 3):
                    # go forward
                    if (ma > 359.9):
                        time = nocon/satm.xno*(satm.xmo/twopi - 1)
                    # back up
                    else:
                        time = satm.xmo*nocon/(satm.xno*twopi)
                    t1 = Date(time=(satm.jd - time))
                    delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                    # satm.delta_t(t1.jd) # FIXME: Replace with python-SGP4 call 
                    satm.jd = t1.jd
                    # satm.rv2el(satm.rr, satm.vv)
                    ma = fmod(degrees(satm.xmo), 360)
                print("\nPERIGEE")
                satm.print_el()       # print new elements

                # perigee residual
                time = sat.jd                     # save sat epoch
                delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                # sat.delta_t(satm.jd)              # move sat to perigee FIXME python-SGP4
                delr = sat.rr - satm.rr           # compare sat and satm perigees
                delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                # sat.delta_t(time)                 # restore sat epoch  FIXME python-SGP4
                print("\nperigee delta {:5.0f}".format(mag(delr)*6378.135))

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
                    if (ma < 180.1):
                        time = satm.jd - 0.5*nocon/satm.xno*(1 + satm.xmo/pi)
                    # if previously past apogee and before next perigee, back up to apogee
                    else:
                        time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi)
                # NEXT: advance to apogee
                elif (buf == 'N'):
                    # if previously at or past apogee and before perigee, go forward to apogee
                    if (ma > 179.9):
                        time = satm.jd + 0.5*nocon/satm.xno*(3 - satm.xmo/pi)
                    # if previously past apogee and past next perigee, go forward to apogee
                    else:
                        time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi)

                # move time and satm.xmo to apogee
                t2 = Date(time=time)
                delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                # satm.delta_t(t2.jd) # FIXME python-SGP4
                satm.jd = t2.jd
                # satm.rv2el(satm.rr, satm.vv) # FIXME python-SGP4
                for i in range(0, 3):
                    # loop to refine apogee, find when satm.xmo = pi
                    t1 = Date(time=(satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi)))
                    delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                    # satm.delta_t(t1.jd) # FIXME python-SGP4
                    satm.jd = t1.jd
                    # satm.rv2el(satm.rr, satm.vv) # FIXME python-SGP4
                ma = fmod(degrees(satm.xmo), 360)
                print("\nAPOGEE")
                satm.print_el()       # print new elements

                # apogee residual
                time = sat.jd                     # save sat epoch
                delta_t(satm,t2.jd) # FIXME: - this needs a delta time
                # sat.delta_t(satm.jd)              # move sat to apogee # FIXME python-SGP4
                delr = sat.rr - satm.rr           # compare sat and satm apogees
                delta_t(sat,t2.jd) # FIXME: - this needs a delta time
                # sat.delta_t(time)                 # restore sat epoch # FIXME python-SGP4
                print("\napogee delta {:5.0f}".format(mag(delr)*6378.135))  # kilometers

            # O(b) Pseudo observation?
            elif (buf == 'O'):
                # move up one to make room for pseudo ob
                # FIXME: Probably a pythonic way to do this 
                for i in range(nobs-1, -1, -1):
                    iod_line[i+1] = iod_line[i]
                    ll[i+1] = ll[i]
                    rd[i+1] = rd[i]
                    odata[i+1] = odata[i]

                nobs+=1 # FIXME Global

                # ll is unit vector in same direction as satm.rr
                ll[0] = satm.rr/mag(satm.rr)
                
                # rd is unit vector (1 er) in same direction as satm.rr
                rd[0] = satm.rr/mag(satm.rr)
                # odata holds epoch, ra, dc, obscode data

                odata[0][0] = satm.jd
                (ra, dc) = zrll(satm, rd[0]) # Get predicted position
                odata[0][1] = radians(ra)
                odata[0][2] = radians(dc)
                odata[0][3] = 0000

                # print pseuso ob
                #           1         2         3         4         5         6
                # 0123456789012345678901234567890123456789012345678901234567890
                # 36868 10 039A   2018 G 20100909203108300 17 25 2216290+034350 57

                norad = int(iod_line[1])
                yy    = int(iod_line[1] + 6)
                desig = "{:4s}".format(iod_line[1] + 9)
                desig[4] = '\0'
                t1 = Date(time=(satm.jd))
                ra /= 15
                rm = fmod(ra, rm)*60000
                dm = fabs(fmod(dc, dm)*6000)
                print("\nPSEUDO OB:")
                iod_line[0] = "\n{:05d} {:02d} {:4s}   0000 P {:4.0f}{:02.0f}{:02.0f}{:02.0f}{:02.0f}{:05.0f} 16 25 {:02.0f}{:05.0f}{:+03.0f}{:04.0f}\n".form(
                    norad, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss*1000, ra, rm, dc, dm)
                print(iod_line[0])
            # Energy 
            elif (buf == 'E'):
                buf = input("\ndelta specific energy(m^2/s^2): ")
                buf = buf.strip()

                try:
                    dE = float(buf)
                    dE /= 11300.168   # er^2 / min^2
                    vec = np.cross(satm.rr, satm.vv)
                    mu = mag(vec)
                    mu = mu*mu
                    mu = mu / satm.aodp
                    mu = mu / (1 - satm.eo*satm.eo)
                    E2 = -0.5*mu / satm.aodp
                    VV = sqrt(2*(E2 + dE + mu/mag(satm.rr)))  # new velocity magnitude
                    dV = unit_vector(satm.vv)   # velocity direction
                    dev = VV * dV               # new velocity vector
                    sat = satm
                    sat.rv2el(sat.rr, vec)      # State vectors to mean elements #FIXME: python-SGP4
                    print_el(sat)              # print new elements
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
                    print_el(sat)         # print original elements
            # QUIT Maneuver        
            elif (buf == 'Q'):
                tle = sat.tle
                ii = degrees(sat.xincl)
                om = fmod(degrees(sat.xnodeo), 360)
                ec = sat.eo
                ww = fmod(degrees(sat.omegao), 360)
                ma = fmod(degrees(sat.xmo), 360)
                nn = sat.xno / nocon
                bstar = sat.bstar
                accept_command()

def write_el(sat):
    save_sat = copy.deepcopy(sat)
    while(True): # Forever loop
        tle = sat.epoch
        ii = degrees(sat.xincl)
        om = degrees(sat.xnodeo)
        ec = sat.ecco
        ww = degrees(sat.omegao)
        ma = degrees(sat.xmo)
        nn = sat.xno / nocon
        c2 = sat.c2
        bstar = sat.bstar

        print("\n(U)pdate (V)iew (A)ppend (O)riginal (R)estore (Q)uit")
        buf = input(": ").strip().upper()

        # Append
        if (buf == 'A'):
            write_tle(file) # FIXME replace with tle_util version
        
        # View
        elif (buf == 'V'):
            print_el(sat)

        # View Original
        elif (buf == 'O'):
            print_el(sat)       # print new elements

        # Restore
        elif (buf == 'R'):
            # replace TLE in file with original
            with open(file, "w") as fp:
                fp.write("{:s}".format(sat.name))
                fp.write("{:s}".format(sat.tle_line1))
                fp.write("{:s}".format(sat.tle_line2))
                for range in (0, nobs): # FIXME More pythonic way
                    fp.write("{:s}".format(iod_line[i]))
            # replace working elements with original
            sat = save_sat
            print_el(sat)            # print original elements
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
        # QUIT write_el
        elif (buf == 'Q'):
            return True

def time_func(sat, rd, ll, odata):
    pass
    # while(True): # Forever loop
    #     ec = sat.eo
    #     ww = degrees(sat.omegao)
    #     nn = sat.xno / nocon
    #     xns = 2160 * sat.bstar * sat.c2 * sat.xno / nocon
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
    #     nn = sat.xno / nocon
    #     bstar = sat.bstar
    #     print_el(sat)       # print new elements
    #     sum = find_rms(sat, rd, ll, odata)
    #     print("\nrms{:12.5f}".format(sum))

    #     uu = longitude(degrees(sat.mo), degrees(sat.argpo))

    #     accept_command()

# /////////////////// MAIN //////////////////////////////////////////////////////

def main():
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

    srch = 'W' # initialize wide search

    # TODO: Expand this code to deal with command line options like the original
    file = "sat.txt"
    # TODO: Look into batch capability

    # restart: # FIXME: get rid of this go-to block

    #   // find observation lines in input file
    #   // number of observations = nobs (global)
    #   // store observations to iod_line matrix (up to 200 observations)
    IOD_Records = iod.get_iod_records_from_file(file)
    nobs = len(IOD_Records)
    for item in IOD_Records:
        iod_line.append(item.line)

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
    TLEs = TLEFile(file, strict=False)

    # Grab just the first TLE instance from the file
    # TODO: Need to fix a problem in tle_util.parse_tles()s which only retains the most-recently read TLE for a sat
    for satnum in TLEs.Satellites:
        FitSat = TLEs.Satellites[satnum]
        name = FitSat.name
        tle1 = FitSat.line1
        tle2 = FitSat.line2
        print("{}".format(name))
        print("{}".format(tle1))
        print("{}".format(tle2))
        break

    # # Most popular const used by TLEs
    whichconst = earth_gravity.wgs72

    # # (rr,vv) = sgp4(id_satrec, tsince=0, whichconst=whichconst)
    # id_sat = sgp4.io.twoline2rv(UNIDsat.line1, UNIDsat.line2, whichconst, afspc_mode=False)
    # (year, monnth, day, hour, minute, second) = UNIDsat.epoch.timetuple()[:6]

    # // Create Satellite variable, sat, from TLE
    # TODO: - Replace with satellite Class variable
    # sat(tle, ii, om, ec, ww, ma, nn, bstar) # scott-campbell version
    (satn, epoch, bstar, ec, ww, ii, ma, nn, om) = FitSat.satrec
    # sat = twoline2rv(tle1,tle2,whichconst,afspc_mode)
    sat = Satellite()
    sat.whichconst = whichconst
    sat.error = 0
    sat.operationmode = False
    sat.satnum = satn
    sat.name = name
    sat.jdsatepoch = FitSat.jdsatepoch
    sat.epoch = FitSat.epoch
    sat.intldesg = FitSat.designation

    sgp4init(whichconst, False, satn, epoch, bstar, 0, 0, ec, ww, ii, ma, nn, om, sat)
#    sgp4init(satn, epoch, bstar, 0, 0, ec, ww, ii, ma, nn, om, sat)
#    sgp4init(whichconst, afspc_mode, satn, epoch, bstar, ec, ww, ii, ma, nn, om, sat)
 
    # Make a copy of original sat
    save_sat = copy.deepcopy(sat)

    # // maneuvering sat
    satm = copy.deepcopy(sat)

    # // calculate uu, degrees, for search
    uu = longitude(degrees(sat.mo), degrees(sat.argpo), sat.ecco)

    # // read all observations from input
    # // Observed topocentric vectors for each observer position.
    # // Equatorial coordinates.
    # TODO: Replace with DB query (combined with get_obs() above?)
    # TODO: Audit this against everything done in read_obs starting on line 1585
    # get line-of-sight vectors
    # (odata, ll, rd) = read_obs(IOD_Records)
    (odata, ll, rd) = read_obssf(IOD_Records) # Use the scott-campbell way while debugging against satfit.cpp

    # Should be able to have the query sort for us, although we'll need this for a file-op
    [iod_line, odata, ll, rd] = sort(iod_line, odata, ll, rd)

    sum = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    print("\n{} Observations Found".format(nobs))

    # print dates
    t2 = Date()
    print("\nTODAY   : {:d}".format(t2.doy))
    if(nobs > 0):
        t3 = Date(time=(odata[nobs - 1][0]))
        print("LAST OB : {:d}".format(t3.doy))

    # Accept a new command
    accept_command(sat, rd, ll, odata, sum, uu)

if __name__ == '__main__':
    main()