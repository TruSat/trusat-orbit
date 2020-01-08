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

import os
from datetime import datetime
from time import time                         # For performance timing
from math import radians, pi                    
import pickle   # For bounty work quick setup

import logging
log = logging.getLogger(__name__)

# Scalar python SGP4 from git+https://github.com/interplanetarychris/python-sgp4@cython-7-dec-15-vallado
# Until the following pull request is approved
# https://github.com/brandon-rhodes/python-sgp4/pull/35

# Vector python SGP4 from git+https://github.com/interplanetarychris/python-sgp4@7-dec-15-vallado-tsince-vectorize

# Use local/dev version of python-sgp4
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sgp4_path = os.path.join(parentdir, "python-sgp4")
sys.path.insert(1,sgp4_path) 

try:
    from sgp4.cpropagation import sgp4, sgp4init
    from satfit_accelerated import *
    from sgp4.cmodel import Satellite
except ImportError as e:
    print(e)
    from sgp4.propagation import sgp4, sgp4init
    from sgp4.model import Satellite
from sgp4 import earth_gravity

# Our own functions
from tle_util import TruSatellite
# from satfit_accelerated import step, find_rms, move_epoch_to_jd, longitude

# ///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////
twopi = 2*pi
nocon = twopi/1440.0
de2ra = pi/180.0


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
                self.jd = jday(self.yy, self.mm, self.dd, self.hr, self.mm, self.ss)
                self.sidereal()
            else:                   # this date is julian
                self.jd = self.time
                self.calcmjd()
                self.sidereal()
                self.time = jday_to_datetime(self.jd)
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
            self.jd = jday(self.yy, self.mm, self.dd, self.hr, self.mn, self.ss)
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
    satrec                = Satellite(gravconst)

    satrec.line0          = TLE.line0

    # Line 1
    satrec.satnum         = TLE.satellite_number
    satrec.line1          = TLE.line1
    satrec.classification = TLE.classification
    satrec.intldesg       = TLE.designation
    satrec.epochyr        = 0 # TLE._epoch_year
    satrec.epochdays      = 0 # TLE._epoch_day
    satrec.ndot           = TLE.mean_motion_derivative
    satrec.nddot          = TLE.mean_motion_sec_derivative
    satrec.bstar          = TLE.bstar
    satrec.ephtype        = TLE.ephemeris_type
    satrec.elnum          = TLE.element_num

    # Line 2
    satrec.line2          = TLE.line2
    satrec.inclo          = TLE.inclination_radians  # rad
    satrec.nodeo          = TLE.raan_radians         # rad
    satrec.ecco           = TLE.eccentricity   
    satrec.argpo          = TLE.arg_perigee_radians  # rad
    satrec.mo             = TLE.mean_anomaly_radians # rad

    # TODO: Once mean_motion_radians_per_minute is added to the DB, use it directly here
    satrec.no_kozai       = TLE.mean_motion_orbits_per_day * nocon # rad/min
    satrec.revnum         = TLE.orbit_number

    # Derived quantities
    satrec.jdsatepoch     = TLE.jdsatepoch     # Julian date
    satrec.jdSGP4epoch    = satrec.jdsatepoch - 2433281.5
    # satrec.epoch_datetime = TLE.epoch_datetime # Python datetime

    # Pass the source tle_id through the SGP4 class variable, for TLE genealogy
    satrec.parent_tle_id = 0 # int(TLE.tle_id)

    # SGP4 mode variables
    satrec.operationmode  = u'i' # Unicode for cython
    satrec.error          = 0

    if (gravconst == "wgs72old"):
        whichconst = earth_gravity.wgs72old
    elif (gravconst == "wgs84"):
        whichconst = earth_gravity.wgs84
    else:
        # Most popular const used by TLEs
        whichconst = earth_gravity.wgs72
    # satrec.whichconst     = whichconst  # Python extension: remembers its consts
    # satrec.whichconst = gravconst
    rtn_code = sgp4init("wgs72", satrec.operationmode, satrec.satnum, 
             satrec.jdSGP4epoch, # epoch time in days from jan 0, 1950. 0 hr
             satrec.bstar, satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, 
             satrec.inclo, satrec.mo, satrec.no_kozai, satrec.nodeo, satrec)
    if (rtn_code is not True):
        if (satrec.error == 1):
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
    else:
        return satrec


# /////////////////// MAIN //////////////////////////////////////////////////////
def main():
    """ satfit_performance
    Fit a TLE prediction to a reference TLE + new IOD observations

    Facilitate standard test environment and performance metrics for assessment

    Inputs:
        iod_line    Array of IOD formatted observations
        file_in     source file (contains IODs plus TLEs)
        file_out    result file (contains file_in + new TLE predict)
    """
    
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

    # Convenience: Load IOD observation, station data pickle
    with open('satfit_performance.pickle', 'rb') as f:
        (odata, ll, rd, t1) = pickle.load(f)

    odata = [[2.45871561e+06, 1.55664956e+00, 1.24212334e+00, 4.17200000e+03, 4.08299000e+05],
 [2.45872254e+06, 1.54232514e+00, 9.74737447e-01, 4.17200000e+03, 4.11181000e+05],
 [2.45872254e+06, 1.59489793e+00, 9.67111755e-01, 4.17200000e+03, 4.11182000e+05],
 [2.45872254e+06, 1.66383389e+00, 9.55289479e-01, 4.17200000e+03, 4.11183000e+05],
 [2.45874841e+06, 3.47627697e+00, 9.32449088e-01, 4.17200000e+03, 4.16464000e+05],
 [2.45874841e+06, 3.41200155e+00, 9.48135077e-01, 4.17200000e+03, 4.16465000e+05],
 [2.45874841e+06, 3.34725501e+00, 9.61012729e-01, 4.17200000e+03, 4.16466000e+05],
 [2.45879826e+06, 4.25884112e+00, 1.15787821e+00, 4.17200000e+03, 4.21414000e+05],
 [2.45879826e+06, 4.35459538e+00, 1.16465128e+00, 4.17200000e+03, 4.21415000e+05],
 [2.45879826e+06, 4.43039711e+00, 1.16798418e+00, 4.17200000e+03, 4.21416000e+05]]

    ll = [[ 0.00456624,  0.32275503,  0.94647152],
 [ 0.01598115,  0.56115783,  0.82755452],
 [-0.01368066,  0.56751477,  0.82324955],
 [-0.05363986,  0.57487537,  0.81648091],
 [-0.5628066,  -0.19572568,  0.80308168],
 [-0.5620066,  -0.15578733,  0.8123293 ],
 [-0.56062117, -0.11695225,  0.81977197],
 [-0.17582553, -0.36071339,  0.91595373],
 [-0.13835716, -0.37005179,  0.91865063],
 [-0.10908354, -0.37652389,  0.91996225]]

    rd = [[ 0.59544982,  0.14061455,  0.78833855],
 [ 0.60865287, -0.06224648,  0.78833855],
 [ 0.60868465, -0.06193494,  0.78833855],
 [ 0.60872963, -0.06149125,  0.78833855],
 [ 0.55200408, -0.26386443,  0.78833855],
 [ 0.55219627, -0.26346199,  0.78833855],
 [ 0.55238828, -0.26305917,  0.78833855],
 [ 0.52008562, -0.32224818,  0.78833855],
 [ 0.52024938, -0.32198372,  0.78833855],
 [ 0.52036712, -0.32179341,  0.78833855]]


    # print(type(odata))
    # print(type(ll))
    # print(type(rd))
    # print(type(t1))

    # REFERENCE TLE
    line0 = 'SL-16 R/B'
    line1 = '1 22285U 92093B   19314.09990558 -.00000000  00000-0  24310-4 0  9990'
    line2 = '2 22285  71.0198 174.7928 0006190 244.6688 115.3794 14.15033886386443'

    TLE = TruSatellite(line0=line0, line1=line1, line2=line2)
    sat = initsat(TLE)
    sat.thetag = t1.thetag  # used in anomaly()

    # // calculate uu, degrees, for search
    uu = longitude(sat)

    rd = np.array([[0.59544982, 0.14061455, 0.78833855], [0.60865287, -0.06224648, 0.78833855], [0.60868465, -0.06193494, 0.78833855], [0.60872963, -0.06149125, 0.78833855], [0.55200408, -0.26386443, 0.78833855], [0.55219627, -0.26346199, 0.78833855], [0.55238828, -0.26305917, 0.78833855], [0.52008562, -0.32224818, 0.78833855], [0.52024938, -0.32198372, 0.78833855], [0.52036712, -0.32179341, 0.78833855]],dtype=np.double)
    ll = np.array([[0.00456624, 0.32275503, 0.94647152], [0.01598115, 0.56115783, 0.82755452], [-0.01368066, 0.56751477, 0.82324955], [-0.05363986, 0.57487537, 0.81648091], [-0.5628066, -0.19572568, 0.80308168], [-0.5620066, -0.15578733, 0.8123293], [-0.56062117, -0.11695225, 0.81977197], [-0.17582553, -0.36071339, 0.91595373], [-0.13835716, -0.37005179, 0.91865063], [-0.10908354, -0.37652389, 0.91996225]]
,dtype=np.double)
    odata = np.array([[2458715.61, 1.55664956, 1.24212334, 4172.0, 408299.0], [2458722.54, 1.54232514, 0.974737447, 4172.0, 411181.0], [2458722.54, 1.59489793, 0.967111755, 4172.0, 411182.0], [2458722.54, 1.66383389, 0.955289479, 4172.0, 411183.0], [2458748.41, 3.47627697, 0.932449088, 4172.0, 416464.0], [2458748.41, 3.41200155, 0.948135077, 4172.0, 416465.0], [2458748.41, 3.34725501, 0.961012729, 4172.0, 416466.0], [2458798.26, 4.25884112, 1.15787821, 4172.0, 421414.0], [2458798.26, 4.35459538, 1.16465128, 4172.0, 421415.0], [2458798.26, 4.43039711, 1.16798418, 4172.0, 421416.0]]
,dtype=np.double)

    # print(type(rd))
    # print(type(ll))
    # print(type(odata))
    print(type(sat))
    print(sat)
    start_rms = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    print(f"Start rms: {start_rms}")

    sat = move_epoch_to_jd(sat,t1.jd)
    sat = step(sat, rd, ll, odata, start_rms, uu, "L")
    stop_rms  = find_rms(sat, rd, ll, odata)

    print(f"Stop rms: {stop_rms}")

    # start_rms = find_rms(sat, rd, ll, odata) # Establish baseline rms for global
    # sat = step(sat, rd, ll, odata, start_rms, uu, "L")
    # stop_rms  = find_rms(sat, rd, ll, odata)

    # Approximate Fitted TLE
    line0 = 'SL-16 R/B'
    line1 = '1 22285T 92093B   19314.75942444 -0.00000000  00000-0  21616-4 0   03'
    line2 = '2 22285  71.0120 173.4241 0006190 243.6645 235.0266 14.15037779386403'


if __name__ == '__main__':
    main()