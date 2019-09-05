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

from spacetrack import SpaceTrackClient

# These are necessary until Brandon Rhodes approves pull requests
# https://github.com/brandon-rhodes/python-sgp4/pull/35
sys.path.insert(1, '/Users/chris/Dropbox/code/MVP/python-sgp4')
# https://github.com/skyfielders/python-skyfield/pull/276
sys.path.insert(2, '/Users/chris/Dropbox/code/MVP/python-skyfield')

# from skyfield.iokit import Loader, download, parse_tle
# from skyfield import sgp4lib
from sgp4.ext import jday, invjday, days2mdhms
from sgp4.io import twoline2rv
from sgp4.model import Satellite
from sgp4.propagation import sgp4init, sgp4
from sgp4 import earth_gravity

# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
tle_path = os.path.join(parentdir, "sathunt-tle")
sys.path.insert(1,tle_path) 
from tle_util import make_tle, append_tle_file, TLEFile, tle_fmt_epoch, datetime_from_tle_fmt, assumed_decimal_point, checksum_tle_line

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
iod_path = os.path.join(parentdir, "sathunt-iod")
sys.path.insert(1,iod_path) 
import iod

from elfind import read_obs, rref, SGN, so2r
from satid import unit_vector, mag

# ///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////
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

# // please note that these variables are in TLE format, i.e. degrees, whereas
# // their counterparts inside SGP4 are in radians.

# double la, lo, hh;          // observer parameters
# int first = 1, batch = 0, out = 0;
# int nobs, len;              // number of observations
nobs = 0             # number of observations
file = None             # input file string declared
# char buf[81];               // input buffer
# char name[81], tle1[81], tle2[81];    // original tle data lines
line1 = None  # read and write tle buffers
line2 = None  # read and write tle buffers
# double sum, nsum, xsum;     // sum of squares of residuals
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
                self.jd = jday(self.yy, self.mm, self.dd, self.hh, self.mm, self.ss)
                self.sidereal()
            else:                   # this date is julian
                self.jd = self.time
                self.calcmjd()
                (self.yy, self.mm, self.dd, self.hh, self.mn, self.ss) = invjday(self.jd)
                self.sidereal()
                (_, subsec) = divmod(self.ss,1)
                subsec = int(subsec*1E6)
                self.ss = int(self.ss)
                self.time = datetime(self.yy,self.mm,self.dd,self.hh,self.mn,self.ss, subsec)
        elif (self.yy):
            """ Creates a date object, t1, initialized to the calendar date and
            time passed by the six calendar variables """
            # TODO: Deal with other variables potentially being "None"
            self.time = datetime(self.yy,self.mm,self.dd,self.hh,self.mn,self.ss)
        else:
            # Shouldn't ever get here, default is to create the current time
            pass

        if (not self.time):
            self.time = datetime(self.yy,self.mm,self.dd,self.hh,self.mn,self.ss)

        # Fill out rest of internal variables
        if (not self.jd):
            self.jd = jday(self.yy, self.mm, self.dd, self.hh, self.mm, self.ss)
            self.calcmjd()

        if (not self.doy):
            self.doy = self.time.timetuple().tm_yday

        if (not self.tle):
            self.tle = tle_fmt_epoch(self.time)

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
        self.hh = self.time.hour
        self.mn = self.time.minute
        self.ss = self.time.second

    def calcmjd(self):
        self.mjd = self.jd - 2400000.5

def acose(x):
    if ( x >= 1):
        rval = 0
    elif ( x <= -1):
        rval = pi
    else:
        rval = acos(x)
    return rval

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
        xno    = degrees(sat.no_unkozai) * 1440/360 # revs per day

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
        print("{:s}".format(line1))
        print("{:s}".format(line2))
    return(line0, line1, line2)


def longitude(ma, ww):
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
    # nobs = len(odata)
    global nobs # FIXME shouldn't have to do this if we're not changing it

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

def find_rms(sat, rd, ll, odata):
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    zum = 0

    # copy sat
    satx = copy.deepcopy(sat)

    # nobs = len(odata)

    for j in range(nobs):
        # advance satellite position
        # TODO: replace this python-SGP4 function
        tsince = odata[j][0] - satx.jdsatepoch
        (satxrr, satxvv) = sgp4(satx,tsince,whichconst) 
        # satx.delta_t(odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # TODO: replace this with python-SGP4 code
        delr = satxrr - rd[j] # This was vmadd -1
        nrr = mag(delr)       # range

        # convert to unit vector
        delr = delr/nrr

        # topocentric position error in degrees
        Perr = degrees( acose( np.dot(delr, ll[j]) ) )

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)

def print_fit(sat, rd, ll, odata):
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
        # FIXME Update this to use python-SGP4
        tsince = odata[j][0] - satx.jdsatepoch
        (satxrr, satxvv) = sgp4(satx,tsince,whichconst) 
        # satx.delta_t(odata[j][0])
        nrr = mag(satxrr)
        nvv = mag(satxvv)

        # computing elevation
        el = degrees( acose(np.dot(rd[j], ll[j]) / mag(rd[j])) )
        el = 90 - el
        el = round(el, 1)

        # computing aspect
        asp = degrees( acose(np.dot(ll[j], satxvv) / nvv) )
        asp = 180 - asp
        asp = round(asp, 1)

        # computing azimuth
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
        Perr = acose(np.dot(satxrr, rr) / (nrr*nrr))

        # difference between computed and observed position vectors, delr, in e.r.
        delr = satxrr - rr
        temp = np.cross(satxrr, satxvv)  # xtrk reference vector points left of track
        sign = SGN(np.dot(delr, temp))

        # observer velocity vector
        tempv = np.cross(zz, rd[j])
        temp = .004351409367 * tempv 
        # observed satellite velocity
        tempv = satxvv - temp
        nvv = mag(tempv)

        # angle between delr vector and tempv vector, radians
        alpha = acose(np.dot(tempv, delr) / (nvv * mag(delr)))

        # magnitude of delr in direction of tempv, radians
        delt = atan(cos(alpha) * tan(Perr));   # geocentric range error

        # time error
        delt *= nrr / nvv                     # delta r in min
        delt *= 60                            # seconds
        delt  = round(delt, 2)

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # new use of delr variable, predicted line of sight vector
        delr = satxrr - rd[j]
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
        Perr = radians(Perr)
        sum += Perr*Perr

        # Format time string
        # YYday HHMM:SSsss
        timestring = satx.epoch.strftime('%y%j %H%M:%S')
        SSS = satx.epoch.strftime('%f')
        SSS = int(1000*(int(SSS)/1E6))
        fit_string = "({:2d}) {:04d}  {}{:03d}  {:5.1f}  {:5.1f}  {:5.1f}  {:6.2f}   {:6.2f}  {:7.3f}".format(
            j + 1, int(odata[j][3]), timestring, SSS, az, el, asp, xtrk, delt, Perr)
        print(fit_string)

        # print fit to file
        # FIXME this (global?) variable
        if (out):
            with open(file, "a") as fp:
                fp.write(fit_string)

    print("\nrms{:12.5f})".format(sqrt(sum / nobs)))

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


def diff_el(sat, rd, ll, odata):
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
    #loop:
    # begin the main loop of the program
    for i in range(0, nobs):
        satx = sat
        satz = sat                      # differential sat
        satx.delta_t(odata[i][0]);      # relocate satellite at new time # python-SGP4

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
            el = satz.xincl                      # store reference
            satz.xincl += delta                  # delta element
            satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            satz.xincl = el                      # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

        j = 1
        delta = 0.001                        # change
        el = satz.xnodeo                     # store reference
        satz.xnodeo += delta                 # delta element
        satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME: python-SGP4
        (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
        satz.xnodeo = el                     # restore reference
        dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
        dat[1][j] = (dc - dcc) / delta

        j = 2
        if (xe):
            dat[0][j] = 0.00001
            dat[1][j] = 0.00001
        else:
            delta = 0.0001                       # change
            el = satz.eo                         # store reference
            satz.eo += delta                     # delta element
            satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            satz.eo = el                         # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

        j = 3
        if (xw):
            dat[0][j] = 0.001
            dat[1][j] = 0.001
        else:
            delta = 0.001                        # change
            el = satz.omegao                     # store reference
            satz.omegao += delta                 # delta element
            satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            satz.omegao = el                     # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

        j = 4
        delta = 0.001                        # change
        el = satz.xmo                        # store reference
        satz.xmo += delta                    # delta element
        satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME: python-SGP4
        (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
        satz.xmo = el                        # restore reference
        dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
        dat[1][j] = (dc - dcc) / delta

        j = 5
        if (xn):
            dat[0][j] = 0.000001
            dat[1][j] = 0.000001
        else:
            delta = 0.00001                      # change
            el = satz.xno                        # store reference
            satz.xno += delta                    # delta element
            satz.delta_t(odata[i][0])            # recalculate with perturbed element # FIXME python-SGP4
            (ra, dc) = zrll(satz, rd[i])         # perturbed ra, dc
            satz.xno = el                        # restore reference
            dat[0][j] = rtw(ra, rac) / delta     # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

        mdata[2 * i]     = dat[0]   # numerical deltas transferred to
        mdata[2 * i + 1] = dat[1]   # multiple regresssion matrix
        # END for i

    # multiple regression
    for j in range(0, 6):
        for k in range (0, 7):
            rv[j][k] = 0
            for i in range (0, nobs*2):
                rv[j][k] = rv[j][k] + mdata[i][j] * mdata[i][k]
    rref(rv, b)

    saty = sat
    # test update components with deltas
    saty.xincl  += b[0]*0.1
    saty.xnodeo += b[1]*0.1
    saty.eo     += b[2]*0.1
    saty.omegao += b[3]*0.1
    saty.xmo    += b[4]*0.1
    saty.xno    += b[5]*0.1
    rms = find_rms(saty, rd, ll, odata)
    if (rms < sum):
        sum = rms
        # update components with deltas
        sat.xincl  += b[0]*0.1
        sat.xnodeo += b[1]*0.1
        sat.eo     += b[2]*0.1
        sat.omegao += b[3]*0.1
        sat.xmo    += b[4]*0.1
        sat.xno    += b[5]*0.1
        c+=1
        if (c < 20):
            next # Back up to the top of the loop
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
    ii = degrees(sat.xincl)
    om = degrees(sat.xnodeo)
    ec = sat.eo
    ww = degrees(sat.omegao)
    ma = degrees(sat.xmo)
    nn = sat.xno/nocon

def anomaly(sat, rd, ll, odata):
    """ box search, no limits """
    step = 0.1
    mk = ma # FIXME global
    sum = find_rms(sat, rd, ll, odata)

    while(fabs(max - min) > 1e-5):
        min = mk
        max = mk

        # nsum loop
        while(True):
            min = mk - step
            sat.delta_el(sat.jd, ii, om, ec, ww, min, nn, bstar) # FIXME python-SGP4
            nsum = find_rms(sat, rd, ll, odata)
            if (nsum < sum):
                mk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop
        while(True):
            max = mk + step
            sat.delta_el(sat.jd, ii, om, ec, ww, max, nn, bstar)    # FIXME python-SGP4
            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                mk = max
                sum = xsum
                continue # Back to top of xsum loop
        step /= 2
    ma = fmod(mk, 360)
    return ma

def motion_search(sat, rd, ll, odata):
    """ box search, no limits """
    step = 0.1
#    nk = nn
    nk = sat.no_unkozai
    sum = find_rms(sat, rd, ll, odata)

    # Start with this values to get through the loop once
    # CPP has the while evaluation at the end
    min = 0
    max = 1
    while(fabs(max - min) > 1e-10):
        min = nk
        max = nk

        # nsum loop
        while(True):
            min = nk - step
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, min, bstar)    # FIXME python-SGP4
            nsum = find_rms(sat, rd, ll, odata)
            if (nsum < sum):
                nk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop
        while(True):
            max = nk + step
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, max, bstar)    # FIXME python-SGP4
            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                nk = max
                sum = xsum
                continue
            step /= 2
    nn = nk

def node(sat, rd, ll, odata):
    """ partition search on node and inclination within set limits """
    while((imax - imin) > 1e-5):
        istep = (imax - imin) / 20
        ostep = fabs(rtw(omax, omin) / 20)
        for ik in range (imin, imax, istep):
            if (xi):
                imin  = ii
                imax  = ii
                ik    = ii
                istep = 0
            for ok in range(omin, omax, ostep):
                satx(tle, ik, ok, ec, ww, ma, nn, bstar)    # FIXME python-SGP4
                # establish the computed ra, dc, at jdo with no perturbations
                rms = find_rms(satx, rd, ll, odata)
                if (rms < sum):
                    sum = rms
                    ii  = ik
                    om  = ok
            # END for ok
        # END for ik
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4
        imin = ii - istep
        imax = ii + istep
        omin = om - ostep
        omax = om + ostep
    om = fmod(om, 360)

def perigee(sat, rd, ll, odata):
    """ partition search on perigee and eccentricity """
    if (ec > 0.1):
        wmax = ww + 0.1
        wmin = ww - 0.1
        emax = ec * 1.01
        emin = ec * 0.99
    while((wmax - wmin) > 1e-5):
        estep = (emax - emin) / 20
        wstep = (wmax - wmin) / 20
        for wk in range(wmin, wmax, wstep):
            if (xw):
                wmin  = ww
                wmax  = ww
                wk    = ww
                wstep = 0
            theta = radians(uu - wk)    # FIXME globals
            for ek in range(emin, emax, estep):
                if (xe):
                    emin  = ec
                    emax  = ec
                    ek    = ec
                    estep = 0
                e = acose((ek + cos(theta)) / (1 + ek * cos(theta)))
                if (theta > pi):
                    e = 2 * pi - e;
                mk = e - ek * sin(e)
                mk = degrees(mk)

                satx(sat.jd, ii, om, ek, wk, mk, nn, bstar) # FIXME python-SGP4
                # establish the computed ra, dc, at jdo with no perturbations
                rms = find_rms(satx, rd, ll, odata)
                if (rms < sum):
                    sum = rms
                    ec  = ek
                    ww  = wk
                    ma  = mk
            # END for ek
        # END for wk
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4

        emin = ec - estep
        emax = ec + estep
        wmin = ww - wstep
        wmax = ww + wstep
                
    # update mean_anomaly
    anomaly(sat, rd, ll, odata)
    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4

    # update mean_motion
    if (not xn):
        motion_search(sat, rd, ll, odata)
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar) # FIXME python-SGP4

    # calculate uu, degrees
    longitude()

    ww = fmod(ww, 360)
    ma = fmod(ma, 360)

def align(sat, rd, ll, odata):
    """ sets deltaT error at last obs epoch to zero """
    last = nobs - 1 # FIXME global
    satx = sat        # copy sat

    while(fabs(delt) > 1.0e-5): # FIXME global
        # advance satellite position
        satx.delta_t(odata[last][0]) # FIXME python-SGP4
        nrr = mag(satx.rr)
        nvv = mag(satx.vv)

        # observed geocentric position vector, rr
        rr = so2r(nrr, rd[last], ll[last])

        # position error in radians
        Perr = acose(dot(satx.rr, rr) / (nrr*nrr))

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

def fit_out():
   out = 1
   print_fit(sat, rd, ll, odata)
   out = 0
   accept_command()

def elcor():
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
    accept_command()

def viewobs():
    # FIXME Make the global variables local
    print()
    for i in range (0, nobs):
        print("({:2d}) {}".format(i + 1, iod_line[i]))
    print()
    accept_command()

def remove():
    # FIXME Make the global variables local
    global nobs
    global rd
    global ll
    global odata
    global iod_line

    buf = input("\nRemove Line Number : ")
    j = int(buf.strip())

    for i in range(j, nobs):
        rd[i - 1] = rd[i]
        ll[i - 1] = ll[i]
        odata[i - 1] = odata[i]
        iod_line[i - 1] = iod_line[i]
        nobs-=1

def fit(sat, rd, ll, odata):
   out = 0
   print_fit(sat, rd, ll, odata)
   print_el(sat)
   accept_command(sat, rd, ll, odata)

def id():
    # requires SATID in same folder
    # TODO Implement python function call to elcord here
    # sprintf(buf, "satid %s", file);
    # system(buf);
    # FIXME Fix goto restart
    # goto restart
    pass

def history():

    print("\n(A)dd  (G)raphs  (E)dit  (Q)uit  ")

    buf = input(": ")
    buf = buf.strip()
    buf = buf.upper()

    # QUIT
    if (buf == 'Q'):
        # FIXME create a local function to pass on TLE command
        # Looks like this is to write the TLE to screen
        # write_tle((char *)"")
        accept_command()

    # FIXME figure out how to get ssn from global variable
    file_in = "history/{:05d}.txt".format(ssn)

    try:
        fp =  open(file_in, "r")
    except: #FIXME: Trap appropriate error
        fp.close()
        # create history file
        with fopen(file_in, "w") as fpo:
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

    history()

def accept_command(sat, rd, ll, odata):
    """Accept a new command"""

    while(True):    # Forever loop
        print("\nEnter command")
        print(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee   (A)nomaly  (M)otion  (B)star")
        print(" (S)tep   I(D)   (T)ime  (F)it  (V)iew  (R)emove   (W)rite   (Q)uit")

        cmd = input(": ").strip()
        cmd = cmd.upper()

        # Hidden functions
        if (cmd == "G"):
            fit_out()
        elif (cmd == "E"):
            edit_file()
        elif (cmd == "H"):
            history()
        elif (cmd == "Y"):
            edit_history()
        elif (cmd == "C"):
            elcor()
        elif (cmd == "O"):
            discover()
        elif (cmd == "U"):
            maneuver()

        # Visible functions
        elif (cmd == "S"):
            step()
        elif (cmd == "Z"):
            step()
        elif (cmd == "I"):
            print_el(sat)
            incl()
        elif (cmd == "N"):
            print_el(sat)
            node()
        elif (cmd == "X"):
            print_el(sat)
            xntrcty()
        elif (cmd == "P"):
            print_el(sat)
            perigee()
        elif (cmd == "A"):
            print_el(sat)
            anomaly()
        elif (cmd == "M"):
            print_el(sat)
            motion(sat, rd, ll, odata)
        elif (cmd == "B"):
            print_el(sat)
            bstar()

        elif (cmd == "D"):
            id()
        elif (cmd == "T"):
            time_func()
        elif (cmd == "F"):
            fit(sat, rd, ll, odata)
        elif (cmd == "V"):
            viewobs()
        elif (cmd == "R"):
            remove()
        elif (cmd == "W"):
            write_el()
        elif (cmd == "Q"):
            sys.exit(0)
    
def discover():       # partition search
    srch = []   # FIXME global variable
    print("\n(W)ide  (N)arrow  [{}", srch[0])

    buf = input(": ").strip().upper()

    if (buf):
        srch[0] = buf[0]

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

    printf("nmin [{:.8f}".format(nmin))

    buf = input("]: ")
    nmin = float(buf.strip())

    nstep = (nmax - nmin) / nsize

    for wk in range (wmin, wmax, wstep):
        theta = radians(uu - wk)
        for ek in range(emin, emax, estep):
            e = acose((ek + cos(theta)) / (1 + ek * cos(theta)))
            if (theta > pi):
                e = 2 * pi - e
            mk = e - ek * sin(e)
            mk = degrees(mk)
            for nk in range(nmin, nmax, nstep):
                # FIXME update class for this
                satx(sat.jd, ii, om, ek, wk, mk, nk, bstar)
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
    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # update mean_anomaly
    anomaly(sat, rd, ll, odata)
    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # print new elements
    print_fit(sat, rd, ll, odata)
    print_el(sat)

    longitude()

    srch[0] = 'Z' # FIXME, this is probably a global
    accept_command()

def step(step_type):       # partition search within limits set below
    # update mean_anomaly
    anomaly(sat, rd, ll, odata)
    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    emax = ec * 1.1
    emin = ec * 0.9
    wmax = ww + 2
    wmin = ww - 2
    imax = ii + 2
    imin = ii - 2
    omax = om + 2
    omin = om - 2

    print("\nPress Q to Quit    :\n")
    while(True): # Forever loop
        perigee(sat, rd, ll, odata)
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        node(sat, rd, ll, odata)
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

        # FIXME: Figure out how to get nobs as a non-global
        # TODO: Figure out where buf = Z comes from
        if (nobs > 3 and step_type == 'Z'):
            diff_el(sat, rd, ll, odata)
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

        # set new limits
        # FIXME: Figure out globals
        emax = 1.01 * ec
        emin = 0.99 * ec
        wmax = ww + 0.5
        wmin = ww - 0.5
        imax = ii + 0.5
        imin = ii - 0.5
        omax = om + 0.5
        omin = om - 0.5

        print("rms{:12.5f}".format(sum))

        buf = input("    : ")
        buf = buf.strip()
        buf = buf.upper()

        if (buff == 'Q'):
            break

    print_fit(sat, rd, ll, odata)
    print_el(sat)       # print new elements

    srch[0] = 'N' # FIXME: Fix this global
    accept_command()

def incl():
    print("\n(A)uto  (ii)  (Q)uit  ")
    buf = input(": ").strip()

    try: 
        ii = float(buf)
        xi = 1
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        print_fit(sat, rd, ll, odata)
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

            print("\nimax [{:.4f}".format(imax))
            buf = input("]: ")
            imax = float(buf.strip())

            print("\nimin [{:.4f}".format(imin))
            buf = input("]: ")
            imin = float(buf.strip())

            print("\nomax [{:.4f}".format(omax))
            buf = input("]: ")
            omax = float(buf.strip())
            
            print("\nomin [{:.4f}".format(omin))
            buf = input("]: ")
            omin = float(buf.strip())

            node(sat, rd, ll, odata)

            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements

            srch[0] = 'N' # FIXME: Fix this global
    if (buf == 'Q'):
        accept_command()
    incl()

def node():
    print("\n(A)uto  (om)  (Q)uit  ")
    buf = input(": ").strip()

    try:
        om = float(buf)
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        print_fit(sat, rd, ll, odata)
        print_el(sat)       # print new elements
    except:
        buf = buf.upper()

        # AUTO
        if (buf == 'A'):
            satx = sat

            # partition search
            omax = om + 2
            omin = om - 2

            print("\nomax [{:.4f}".format(omax))
            buf = input("]: ")
            omax = float(buf.strip())
            
            print("\nomin [{:.4f}".format(omin))
            buf = input("]: ")
            omin = float(buf.strip())

            while((omax - omin) > 1e-5):
                ostep = fabs(rtw(omax, omin) / 20)
                for ok in range (omin, omax, ostep):
                    satx.delta_el(sat.jd, ii, ok, ec, ww, ma, nn, bstar)
                    # establish the computed ra, dc, at jdo with no perturbations
                    rms = find_rms(satx, rd, ll, odata)
                    if (rms < sum):
                        sum = rms      # global
                        om = ok
                    # end for ok
                satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

                omin = om - ostep
                omax = om + ostep
            om = fmod(om, 360)

            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements

            srch[0] = 'N' # FIXME: Fix this global variable
    if (buf == 'Q'):
        accept_command()
    else:
        node()


def xntrcty():
    print("\n(S)earch  (A)uto  (ec)  (Q)uit  ")
    buf = input(": ").strip()

    emax = ec * 1.05
    emin = ec * 0.95

    try:
        ec = float(buf)
        xe = 1
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        print_fit(sat, rd, ll, odata)
        print_el(sat)       # print new elements
    except:
        buf = buf.upper()

        # AUTO
        if (buf == 'A'):
            xe = 0

            print("\nemax [{:.7f}".format(emax))
            buf = input("]: ")
            emax = float(buf.strip())

            print("\nemin [{:.7f}".format(emin))
            buf = input("]: ")
            emin = float(buf.strip())

            while((emax - emin) > 1.e-8):
                estep = (emax - emin) / 20
                for ek in range(emin, emax, estep):
                    sat.delta_el(sat.jd, ii, om, ek, ww, ma, nn, bstar)
                    # establish the computed ra, dc, at jdo with no perturbations
                    rms = find_rms(sat, rd, ll, odata)
                    if (rms < sum):
                        sum = rms
                        ec = ek
                    # end for ek
                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                emin = ec - estep
                emax = ec + estep

            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements

            srch[0] = 'N' # FIXME this global

    # SEARCH
    if (buf == 'S'):
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
        for ec in range(emin,emax + estep, estep):
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            sum = find_rms(sat, rd, ll, odata)
            print("\n{:.7f}     {:7.4f}", ec, sum)

        print()
        ec = ek
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    if (buf == 'Q'):
        accept_command()
    else:
        xntrcty()

def perigee():
    pgee = (1-ec)*sat.aodp
    agee = (1+ec)*sat.aodp

    print("\nPerigee = {} er".format(pgee))
    print("     {:d} X {:d} km".format(int(((pgee-1)*6378.135)),
                                int(((agee-1)*6378.135))))

    print("\n(S)earch  (A)uto  (ww)  (Q)uit  ")

    buf = input(": ").strip()

    try:
        ww = float(buf)
        longitude()
        ma = ww2ma(ww)
        xw = 1
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        print_fit(sat, rd, ll, odata)
        print_el(sat)       # print new elements
    except:
        buf = buf.upper()
        # AUTO
        if (buf == 'A'):
            xw = 0
            wx = ww + 1
            while(fabs(wx - ww) > 1e-4):
                wx = ww
                print("\n{:8.4f}  {:.7f}".format(ww, ec))
                wmax = ww + 1
                wmin = ww - 1
                emax = ec * 1.01
                emin = ec * 0.99
                perigee(sat, rd, ll, odata)
                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements

            srch[0] = 'N' # FIXME Global

        # SEARCH
        if (buf == 'S'):
            longitude()
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
            for ww in range(wmin, wmax + wstep, wstep):
                ma = ww2ma(ww)
                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = find_rms(sat, rd, ll, odata)
                print("\n{:8.4f}  {:7.4f}".format(fmod(ww, 360), sum))
            print()
            ww = wk
            ma = mk
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

    # QUIT
    if (buf == 'Q'):
        accept_command()
    else:
        perigee()

    # FIXME: Looks like this might be a straggler from the beginning of anomaly() below
    # Not sure how the following code ever gets used
    # amax and amin are used as temporary variables
    longitude()
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


def anomaly():
    if (nn < 1.5):
        # amax and amin are used as temporary variables
        longitude()
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
        longitude()
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        print_fit(sat, rd, ll, odata)
        print_el(sat)       # print new elements
    except:
        buf = buf.upper()

        # AUTO
        if (buf == 'A'):
            anomaly(sat, rd, ll, odata)
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements
            srch[0] = 'N'   # FIXME Global
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
            astep = (amax - amin) / astep;

            mk = ma
            print("\nanomaly        sum")
            for ma in range(amin, amax + astep, astep):
                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                sum = find_rms(sat, rd, ll, odata)
                printf("\n{:8.4f}     {:7.4f}".format(ma, sum))
            print()
            ma = mk              # restore
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)

        # LAST
        elif (buf == 'L'):
            align(sat, rd, ll, odata)
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            longitude()
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements

        # QUIT
        elif (buf == 'Q'):
            accept_command()
    anomaly()

def motion(sat, rd, ll, odata):
    while(True): # Forever loop
        print("\n(A)uto  (nn)  (Q)uit  ")

        buf = input(": ").strip()

        try:
            nn = float(buf)
            xn = 1
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
            print_el(sat)       # print new elements
        except:
            buf = buf.upper()

            # AUTO
            if (buf == 'A'):
                xn = 0
                # update mean motion, no limits
                motion_search(sat, rd, ll, odata)
                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                print_fit(sat, rd, ll, odata)
                print_el(sat)       # print new elements
            elif (buf == 'Q'):
                accept_command()

def bstar():
    while(True): # Forever loop
        print("\n(A)uto  (b*)  (B)atch  (Q)uit  ");
        buf = input(": ").strip()

        try:
            bstar = float(buf)
            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
            print_fit(sat, rd, ll, odata)
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

                print("\nbmax [{:.8f}".format(bmax))
                buf = input("]: ").strip()
                bmax = float(buf)

                print("\nbmin [{:.8f}".format(bmin))
                buf = input("]: ").strip()
                bmin = float(buf)

                while((bmax - bmin) > 1.e-9):
                    bstep = (bmax - bmin) / 20
                    for bk in range(bmin, bmax, bstep):
                        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bk)
                        # establish the computed ra, dc, at jdo with no perturbations
                        rms = find_rms(sat, rd, ll, odata)
                        if (rms < sum):
                            sum = rms
                            bstar = bk
                       # END for bk
                    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                    bmin = bstar - bstep
                    bmax = bstar + bstep

                sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
                print_fit(sat, rd, ll, odata)
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
                accept_command()


def maneuver():
    while(True): # Forever loop
        print("\nBoundary#  (P)erigee  (A)pogee  (O)b  (E)nergy  (R)estore  (Q)uit  ")

        buf = input(": ").strip()

        try:
            p = int(buf)
 
            # store old obs
            iod_linex = iod_line
            llx   = ll
            rdx   = rd
            odata = odata

            for i in range(p, nobs):
                iod_line[i-p] = iod_line[i-1]
                ll[i-p] = ll[i-1]
                rd[i-p] = rd[i-1]
                odata[i-p] = odata[i-1]
            nobsx = nobs
            nobs  = nobs - p + 1

            out = 0
            print_fit(sat, rd, ll, odata)
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
                        time = satm.jd + nocon/satm.xno*(1 - satm.xmo/twopi);
                # move to time and ma at perigee
                t2 = Date(time)
                satm.delta_t(t2.jd) # FIXME: Replace with python-SGP4 call
                satm.rv2el(satm.rr, satm.vv) # FIXME: Replace with python-SGP4 call
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
                    t1 = Date(satm.jd - time)
                    satm.delta_t(t1.jd) # FIXME: Replace with python-SGP4 call 
                    satm.jd = t1.jd
                    satm.rv2el(satm.rr, satm.vv)
                    ma = fmod(degrees(satm.xmo), 360)
                print("\nPERIGEE")
                satm.print_el()       # print new elements

                # perigee residual
                time = sat.jd                     # save sat epoch
                sat.delta_t(satm.jd)              # move sat to perigee FIXME python-SGP4
                delr = sat.rr - satm.rr           # compare sat and satm perigees
                sat.delta_t(time)                 # restore sat epoch  FIXME python-SGP4
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
                t2 = Date(time)
                satm.delta_t(t2.jd) # FIXME python-SGP4
                satm.jd = t2.jd
                satm.rv2el(satm.rr, satm.vv) # FIXME python-SGP4
                for i in range(0, 3):
                    # loop to refine apogee, find when satm.xmo = pi
                    t1 = Date(satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi))
                    satm.delta_t(t1.jd) # FIXME python-SGP4
                    satm.jd = t1.jd
                    satm.rv2el(satm.rr, satm.vv) # FIXME python-SGP4
                ma = fmod(degrees(satm.xmo), 360)
                print("\nAPOGEE")
                satm.print_el()       # print new elements

                # apogee residual
                time = sat.jd                     # save sat epoch
                sat.delta_t(satm.jd)              # move sat to apogee # FIXME python-SGP4
                delr = sat.rr - satm.rr           # compare sat and satm apogees
                sat.delta_t(time)                 # restore sat epoch # FIXME python-SGP4
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
                t1 = Date(satm.jd)
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

def write_el():
    while(True): # Forever loop
        tle = sat.tle
        ii = degrees(sat.xincl)
        om = degrees(sat.xnodeo)
        ec = sat.eo
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
            save_print_el(sat)       # print new elements

        # Restore
        elif (buf == 'R'):
            # replace TLE in file with original
            with open(file, "w") as fp:
                fp.write("{:s}".format(name))
                fp.write("{:s}".format(tle1))
                fp.write("{:s}".format(tle2))
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
                for range in (0, nobs): # FIXME More pythonic way
                    fp.write("{:s}".format(iod_line[i]))
        # QUIT write_el
        elif (buf == 'Q'):
            accept_command()

def time_func():
    while(True): # Forever loop
        ec = sat.eo
        ww = degrees(sat.omegao)
        nn = sat.xno / nocon
        xns = 2160 * sat.bstar * sat.c2 * sat.xno / nocon
        if (nobs > 0):
            time2 = odata[nobs - 1][0]
        else:
            t2.now()
            time2 = t2.jd

        sat.delta_t(time2) # FIXME python-SGP4
        z2 = sat.rr[2]

        while(z1 < 0 or z2 > 0):
            time1 = time2
            z1 = z2
            time2 -= .01
            sat.delta_t(time2)  # FIXME python-SGP4
            z2 = sat.rr[2]

        while(time1 - time2 > 1e-9):
            time3 = (time1 + time2) / 2
            sat.delta_t(time3)  # FIXME python-SGP4
            z3 = sat.rr[2]
            if (z3 < 0):
                time2 = time3
            else:
                time1 = time3

        t1 = Date(time2)
        t1.input()
        sat.delta_t(t1.jd)             # advance to node # FIXME: python-SGP4
        sat.rv2el(sat.rr, sat.vv)      # sgp4 elements # FIXME: python-SGP4
        tle = t1.tle
        sat.jd = t1.jd

        ii = degrees(sat.xincl)
        om = fmod(degrees(sat.xnodeo), 360)
        ec = sat.eo
        ww = fmod(degrees(sat.omegao), 360)
        ma = fmod(degrees(sat.xmo), 360)
        nn = sat.xno / nocon
        bstar = sat.bstar
        print_el(sat)       # print new elements
        sum = find_rms(sat, rd, ll, odata)
        print("\nrms{:12.5f}".format(sum))

        longitude()

        accept_command()

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
    global satn
    global epoch
    global bstar
    global ec
    global ww
    global ii
    global ma
    global nn
    global om
    global file
    global nobs
    global whichconst
    global rd
    global ll
    global odata
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
    for satnum in TLEs.Satellites:
        FitSat = TLEs.Satellites[satnum]
        tle1 = FitSat.line1
        tle2 = FitSat.line2
        print("{}".format(tle1))
        print("{}".format(tle2))
        break

    # # Most popular const used by TLEs
    whichconst = earth_gravity.wgs72
    afspc_mode = False

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
    sat.jdsatepoch = FitSat.jdsatepoch
    sat.epoch = FitSat.epoch
    sat.intldesg = FitSat.designation

    sgp4init(whichconst, afspc_mode, satn, epoch, bstar, 0, 0, ec, ww, ii, ma, nn, om, sat)
#    sgp4init(satn, epoch, bstar, 0, 0, ec, ww, ii, ma, nn, om, sat)
#    sgp4init(whichconst, afspc_mode, satn, epoch, bstar, ec, ww, ii, ma, nn, om, sat)
 
    # Make a copy of original sat
    save_sat = copy.deepcopy(sat)

    # // maneuvering sat
    satm = copy.deepcopy(sat)

    # // calculate uu, degrees, for search
    uu = longitude(ma, ww)

    # // read all observations from input
    # // Observed topocentric vectors for each observer position.
    # // Equatorial coordinates.
    # TODO: Replace with DB query (combined with get_obs() above?)
    # TODO: Audit this against everything done in read_obs starting on line 1585
    # get line-of-sight vectors
    (odata, ll, rd) = read_obs(IOD_Records)

    # Should be able to have the query sort for us, although we'll need this for a file-op
    [iod_line, odata, ll, rd] = sort(iod_line, odata, ll, rd)

    sum = find_rms(sat, rd, ll, odata)
    print("{} Observations Found".format(nobs))

    # print dates
    t2 = Date()
    print("\nTODAY   : {:d}".format(t2.doy))
    if(nobs > 0):
        t3 = Date(time=odata[nobs - 1][0])
        print("\nLAST OB : {:d}".format(t3.doy))

    # Accept a new command
    accept_command(sat, rd, ll, odata)

if __name__ == '__main__':
    main()