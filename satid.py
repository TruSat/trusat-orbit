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
from math import (acos, asin, atan, cos, sin, tan, degrees)    # Fast/precise math functions                      
import numpy as np
import logging
import string

from spacetrack import SpaceTrackClient

# These are necessary until Brandon Rhodes approves pull requests
# https://github.com/brandon-rhodes/python-sgp4/pull/35
sys.path.insert(1, '/Users/chris/Dropbox/code/preMVP/python-sgp4')
# https://github.com/skyfielders/python-skyfield/pull/276
sys.path.insert(2, '/Users/chris/Dropbox/code/preMVP/python-skyfield')

from skyfield.iokit import Loader, download, parse_tle
from skyfield import sgp4lib

# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
tle_path = os.path.join(parentdir, "sathunt-tle")
sys.path.insert(1,tle_path) 
from tle_util import make_tle, append_tle_file

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / mag(vector)

def proj(v2, v1):
    """ Returns the unit vector projection of v1 onto v2 """
    b = np.dot(v2, v1)/np.dot(v2, v2)
    temp = np.multiply(b, v2)

    # Make unit vector
    vp = unit_vector(temp)
    return vp

def flat_proj(v1, v2):
    """ Returns the flat projection of direction unit vector, v1 onto v2 """
    temp1 = np.cross(v1, v2)
    temp2 = np.cross(temp1, v1)
    return  proj(temp2, v2)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793

            Partially Ref: angle(vec1,vec2) in python-sgp4/ext.py
    """
    small     = 0.00000001
    undefined = 999999.1

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    magv1 = mag(v1)
    magv2 = mag(v1)

    if (magv1 * magv2 > small * small):
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    else:
        return undefined


def mag(v):
    """ Computes the magnitude of a vector ||v|| 

    Renamed from norm(v) in original Scott Campbell code
    to better correspond to function names in SGP4 code.
    """
    mag = np.sqrt(np.dot(v, v))
    return mag


def main():
    """ Interactive tool for finding an unknown TLE object within a library of TLEs 
    
    TODO:
    - Implment argv[1] = unid.txt, argv[2] = refer.tle
    - Make non-interactive callable version
    - Make stand-alone verison that uses python-SGP exclusively, not tle_util
    - Incorporate logic to read the (first?) TLE from the UNID file
    - Incorporate logic to warn/error the user if no TLEs found
    - Incorporate Perr/Alpha inputs as command line/config flags
    - Put in more compares for altitude, velocity, etc.
    """
    t0 = time()
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Utility to assist in ID of an unidentified (unid) satellite')
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify configuration file. [Default configuration.ini]",
                             dest='conf_file',
                             nargs='?',
                             const=1,
                             default='configuration.ini',
                             type=str,
                             metavar="FILE")
    conf_parser.add_argument("-d", "--datadir", 
                             help="data directory [default ./data]",
                             dest='datadir',
                             default='./data',                           
                             nargs='?',
                             const=1,                             
                             type=str,                             
                             metavar="DIR")
    conf_parser.add_argument("--tleref",
                             help="Specify TLE reference file. [Default refer.tle]",
                             dest='tle_ref',
                             nargs='?',
                             type=str,
                             metavar="REFTLE")
    conf_parser.add_argument("--tleunid",
                             help="Specify TLE unid file. [Default unid.tle]",
                             dest='tle_unid',
                             nargs='?',
                             type=str,
                             metavar="UNID")
    conf_parser.add_argument("--update", help="update TLEs from online sources",
                             action="store_true")
    conf_parser.add_argument("-dbname", "--database", 
                             help="database to USE",
                             dest='dbname',
                             default='opensatcat_dev',                           
                             nargs='?',
                             const=1,                             
                             type=str,                             
                             metavar="NAME")
    conf_parser.add_argument("-H", "--hostname", 
                             help="database hostname",
                             dest='dbhostname',
                             default='opensatcat.cvpypmmxjtv1.us-east-2.rds.amazonaws.com',
                             nargs='?',
                             const=1,
                             type=str,                             
                             metavar="HOSTNAME")
    conf_parser.add_argument("-u", "--user", 
                             help="database user name",
                             dest='dbusername',
                             nargs='?',
                             type=str,                             
                             metavar="USER")
    conf_parser.add_argument("-p", "--password", 
                             help="database user password",
                             dest='dbpassword',
                             nargs='?',
                             type=str,                             
                             metavar="PASSWD")
    conf_parser.add_argument("-t", "--dbtype", 
                             help="database type [INFILE, sqlserver, sqlite] \
                                   default: INFILE",
                             dest='dbtype',
                             nargs='?',
                             choices=['INFILE', 'sqlserver', 'sqlite'],
                             default='INFILE',
                             type=str,                             
                             metavar="TYPE")
    conf_parser.add_argument("-i", "--import", help="Import TLEs to database",
                             dest='importTLE',
                             action="store_true")
    conf_parser.add_argument("-q", "--quiet", help="Suppress console output",
                             dest='quiet',
                             action="store_true")
    conf_parser.add_argument("-V", "--verbose", 
                             help="increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.",
                             const=1, 
                             default=1, 
                             type=int, 
                             nargs="?")

    # Process commandline options and parse configuration
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    args = conf_parser.parse_args()
    log = logging.getLogger(__name__)

    # make it print to the console.
    console = logging.StreamHandler()
    log.addHandler(console)

    conf_file = args.conf_file
    tle_ref = args.tle_ref
    tle_unid = args.tle_unid
    
    update = args.update
    datadir = args.datadir
    dbname = args.dbname
    dbhostname = args.dbhostname
    dbusername = args.dbusername
    dbpassword = args.dbpassword
    dbtype = args.dbtype
    importTLE = args.importTLE
    verbose = args.verbose
    quiet = args.quiet

    # Set our python-skyfield data directory
    load = Loader(datadir)
    ts = load.timescale()

    if (quiet == False):
        if verbose == 0:
            log.setLevel(logging.WARN) 
        elif verbose == 1:
            log.setLevel(logging.INFO) 
        elif verbose == 2:
            log.setLevel(logging.DEBUG) 
        log.debug("Log level set to {}".format(log.level))

    if verbose:
        for arg in vars(args):
            log.debug("%s : %s",arg, getattr(args, arg))

    cfg.read([args.conf_file])
    log.info("Reading config from: {}".format(args.conf_file))

    # 1st arg in original version
    if not (tle_ref):
        try:
            tle_ref = cfg.get('Common', 'tle_ref')
        except KeyError:
            tle_ref = "refer.tle"

    # 2nd arg in original version
    if not (tle_unid):
        try:
            tle_unid = cfg.get('Common', 'tle_unid')
        except KeyError:
            tle_unid = "unid.txt"

    # # Read single (first?) TLE from UNIDentified TLE file
    # TLE_UNID = tle_util.TLEFile(tle_unid,strict=False)

    # for sat_num in TLE_UNID.Satellites: # Need learn a better method to get just the first/only record
    #         #// id_sat comparison variables
    #         #// Date t1(tle);
    #         #// Satellite id_sat(t1.jd, ii, om, ec, ww, ma, nn, bstar);
    #     UNIDsat = TLE_UNID.Satellites[sat_num]
    #     # echo tle to screen
    #     log.info("{LINE1}\n{LINE2}".format(LINE1=UNIDsat.line1, LINE2=UNIDsat.line2))

    # # Most popular const used by TLEs
    # whichconst = sgp4.earth_gravity.wgs72
    # afspc_mode = False
    # (satn, epoch, xbstar,  xecco, xargpo, xinclo, xmo, xno, xnodeo) = UNIDsat.satrec
    
    # # id_satrec = sgp4init(whichconst, afspc_mode, satn, epoch, xbstar,  xecco, xargpo, xinclo, xmo, xno, xnodeo)
    # # (rr,vv) = sgp4(id_satrec, tsince=0, whichconst=whichconst)

    # id_sat = sgp4.io.twoline2rv(UNIDsat.line1, UNIDsat.line2, whichconst, afspc_mode=False)
    # (year, monnth, day, hour, minute, second) = UNIDsat.epoch.timetuple()[:6]


    UNIDtle = load.tle(url=tle_unid,reload=False)
    # Make sure UNID satellite appears only once
    # UNIDtle = set(UNIDtle.values())
    if(not UNIDtle):
        log.error("No TLEs found in file: {}".format(tle_unid))
        log.error("Run elfind first?")
        sys.exit()
    # Get the single item out of the list
    # [UNID] = UNIDtle
    for satnum in UNIDtle: break
    UNID = UNIDtle[satnum]
    # t_unid = ts.ut1(jd=UNID.model.jdsatepoch)
    t_unid = UNID.epoch

    # Get r,v data at its own EPOCH
    # (rr, vv) = id_sat.propagate(year, monnth, day, hour, minute, second)
    (rr, vv, id_sat_err) = UNID._position_and_velocity_TEME_km(t_unid)

    id_sat_rr = np.array(rr)
    id_sat_vv = np.array(vv)

    # print(id_sat_rr)
    # print(id_sat_vv)

    # flat projection of UNID satellite direction unit vector, vp1
    vp1 = flat_proj(rr, vv)

    # Set Perr error bound
    err1 = input("   position error, degrees [2]: ") 
    err1 = err1 or 2
    err1 = float(err1)

    # Set alpha error bound
    err2 = input("   track angle error, degrees [20]: ")
    err2 = err2 or 20
    err2 = float(err2)

    # Read in REFERENCE element list, and loop through for potential solutions within error bounds
    REFtle = load.tle(url=tle_ref,reload=False)
    # Make sure REFtle satellites appears only once
    REFtle = set(REFtle.values())

    for ref_sat in REFtle:
        # log.debug("Comparing against {}".format(sat_num))
        # if(ref_sat.model.satnum == 26905):
        #     print("here")

        # Get r,v data at UNID epoch
        (rr, vv, ref_sat_err) = ref_sat._position_and_velocity_TEME_km(t_unid)

        ref_sat_rr = np.array(rr)
        ref_sat_vv = np.array(vv)


        # delr - satellite delta r vector
        delr = np.subtract(id_sat_rr, ref_sat_rr)

        # delr - flat projection of delta r unit vector
        delr = flat_proj(id_sat_rr, delr)

        # alpha - angle between delr and id_sat.vv, radians
        alpha = angle_between(delr, id_sat_vv)

        # Per - angle between position unit vectors, radians
        Perr = angle_between(ref_sat_rr, id_sat_rr)

        # delta - magnitude of Perr in direction of id_sat.vv (UNID velocity), radians
        delt = atan(tan(Perr) * cos(alpha))

        # delta_t - time of flight to Closest Point of Approach (cpa) seconds
        # rr, vv already in units of km, km/s. No need to convert.
        delta_t = delt * mag(id_sat_rr) / mag(id_sat_vv)

        # cpa - Closest Point of Approach (cpa), radians
        cpa = asin(sin(alpha) * sin(Perr))

        # vp2 - flat projection of REF satellite direction unit vector
        vp2 = flat_proj(ref_sat_rr, ref_sat_vv)

        # alpha - angle between direction unit vectors, radians
        alpha = acos(np.dot(vp1, vp2))

        # Calculate REF deltas from UNID
        try: 
            alpha = acos(cos(alpha)/cos(delt))
        except ValueError:
            alpha = float('nan')

        # Prepare for presentation to user
        alpha = degrees(alpha) # angle between direction unit vectors
        Perr = degrees(Perr)   # angle between position unit vectors

        # Compare UNID to REF using osculating elements (close enough)
        if((Perr < err1) and (alpha < err2)):

            # tle = None      # epoch of elements in tle format
            # ii = None       # inclination, degrees
            # om = None       # right ascension of ascending node, degrees
            # ec = None       # eccentricity
            # ww = None       # argument of the perigee, degrees
            # ma = None       # mean anomaly, degrees
            # nn = None       # mean motion, revolutions/day
            # uu = None       # true longitude
            # c2 = None       # bstar coefficient
            # bstar = None    # BSTAR drag term
            # name[81] = None

            # visually check match parameters using advanced mean elements
            # Write tle to screen
            (tle_line0, tle_line1, tle_line2) = make_tle(
                name=ref_sat.name,
                ssn=ref_sat.model.satnum, 
                epoch_datetime=ref_sat.epoch.utc_datetime(),
                xincl=ref_sat.model.inclo, 
                xnodeo=ref_sat.model.nodeo, 
                eo=ref_sat.model.ecco,
                omegao=ref_sat.model.argpo,
                xmo=ref_sat.model.mo, 
                xno=degrees(ref_sat.model.no_kozai*1440.0)/360.0,
                deg=False)

            log.info("   position error  {:4.1f}".format(Perr))
            log.info("track angle error  {:4.1f}\n".format(alpha))
            log.info("       time error  {:4.0f}".format(delta_t))
            log.info(" to closest point  {:4.1f}\n".format(degrees(cpa)))

            tle_file_path = os.path.join(datadir, tle_unid)
            append_tle_file(tle_file_path, tle_line0, tle_line1, tle_line2)

            get_continue = input("\n[Next]")

            # //        s_in("\n[Next]", buf);
            # //      }    // if match
        # //   }     // while
        # //   s_in("\n[Done]", buf);
    get_continue = input("\n[Done]")

        # //   system(id_file);
    # // }    // end main

if __name__ == '__main__':
    main()