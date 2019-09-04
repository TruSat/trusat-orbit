#!/usr/bin/env python3
import sys
import os # Only needed for paths while our repos are not public
import configargparse

if sys.version_info[0] != 3 or sys.version_info[1] < 7:
    print("This script requires Python version 3.7")
    sys.exit(1)

import os, re, math, csv
from skyfield.api import Topos
from skyfield.api import load
from skyfield.api import EarthSatellite

# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
database_path = os.path.join(parentdir, "sathunt-database")
sys.path.insert(1,database_path) 
import database

import logging
log = logging.getLogger(__name__)

# def getStation(station_id):
#     return True
    
# def topos(self,cospar=9999):
#     """Function to output Topos tuple for requested COSPAR station ID"""
#     try:
#         return ( [self.site[cospar]['latitude'],self.site[cospar]['longitude'], self.site[cospar]['elevation']] )

def main():
    # Read commandline options
    conf_parser = configargparse.ArgParser(default_config_files=['../trusat.config'], description='Utility to initalize IOD database from email file archive')
    # conf_parser.add_argument("-c", "--config",
    #                             help="config file path",
    #                         required=True, is_config_file=True, )
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
                            default='db.consensys.space',
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
                            help="database type [INFILE, sqlserver, sqlite] default: INFILE",
                            dest='dbtype',
                            nargs='?',
                            choices=['INFILE', 'sqlserver', 'sqlite'],
                            default='INFILE',
                            type=str,                             
                            metavar="TYPE")
    conf_parser.add_argument("-q", "--quiet", help="Suppress console output",
                            dest='quiet',
                            action="store_true")
    conf_parser.add_argument("-V", "--verbose", 
                            help="increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.",
                            const=1, 
                            default=0, 
                            type=int, 
                            nargs="?")

    args = conf_parser.parse_args()
    # Process commandline options and parse configuration
    dbname = args.dbname
    dbhostname = args.dbhostname
    dbusername = args.dbusername
    dbpassword = args.dbpassword
    dbtype = args.dbtype
    verbose = args.verbose
    quiet = args.quiet

    log = logging.getLogger()

    # make it print to the console.
    console = logging.StreamHandler()
    log.addHandler(console)

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

    if (dbtype == "sqlserver"):
        if dbusername == None:
            try: 
                dbusername = input("Username: ") 
            except Exception as error: 
                log.warning('ERROR: password must be specified {}'.format(error))
        if dbpassword == None:
            try: 
                dbpassword = getpass() 
            except Exception as error: 
                log.warning('ERROR: password must be specified {}'.format(error))

    # Set up database connection
    db = database.Database(dbname,dbtype,dbhostname,dbusername,dbpassword)

    ts = load.timescale()

    # 408821	"2019-08-24 13:38:39"	"leobarhorst@gmail.com"	37162	"2010-046A  "	4172	"E"	"20190824031442326"	"2019-08-24 03:14:42.3260"	0.1	NULL	"2"	"5"	2000	105.828	59.80349999999999	0	0	0.30000000000000004	"S"	NULL	NULL	NULL	NULL	NULL	NULL	"IOD"	"37162 10 046A   4172 E 20190824031442326 17 25 0703312+594821 37 S"	1	NULL	"370139ccab6591ba36c6fc2c33081a76"	"2019-08-25 17:44:27"
    # 408820	"2019-08-24 13:38:39"	"leobarhorst@gmail.com"	37162	"2010-046A  "	4172	"E"	"20190824031438341"	"2019-08-24 03:14:38.3410"	0.1	NULL	"2"	"5"	2000	106.752	58.856833333333334	0	0	0.30000000000000004	"S"	NULL	NULL	NULL	NULL	NULL	NULL	"IOD"	"37162 10 046A   4172 E 20190824031438341 17 25 0707008+585141 37 S"	1	NULL	"c24602b9fd0d66538baa8256813ec371"	"2019-08-25 17:44:26"
    # 408819	"2019-08-24 13:38:39"	"leobarhorst@gmail.com"	37162	"2010-046A  "	4172	"E"	"20190824031434357"	"2019-08-24 03:14:34.3570"	0.1	NULL	"2"	"5"	2000	107.57825	57.9475	0	0	0.30000000000000004	"S"	NULL	NULL	NULL	NULL	NULL	NULL	"IOD"	"37162 10 046A   4172 E 20190824031434357 17 25 0710313+575685 37 S"	1	NULL	"9a097ee7bb889d01b0eb43452018e510"	"2019-08-25 17:44:26"

    # iod_obs_id = input("IOD db ID: ")
    while(True):
        try:
            iod_obs_id = input("\nIOD db ID: ")
        except:
            break

        query_tmp = """SELECT object_number, station_number, obs_time, ra, declination, iod_string 
            FROM ParsedIOD
            WHERE obs_id={}""".format(iod_obs_id) 

        try:
            db.c.execute(query_tmp)
            [object_number, station_number, obs_time_sql, ra_obs, decl_obs, iod_string] = db.c.fetchone()
        except:
            log.warning("IOD ID '{}' not found.".format(iod_obs_id))
            continue

        try:
            TLE = db.selectTLEEpochBeforeDate(obs_time_sql, object_number)
            sat_name   = TLE[1]
            tle_line_1 = TLE[2]
            tle_line_2 = TLE[3]
            print("{}\n{}\n{}".format(sat_name,tle_line_1,tle_line_2))
        except:
            log.warning("NO TLE found for object '{}' with EPOCH before {}".format(object_number, obs_time_sql))
            continue

        query_tmp = """SELECT latitude, longitude, elevation_m, Observer.name as name FROM Station
            LEFT JOIN Observer ON Station.user = Observer.id
            WHERE station_num={}
            LIMIT 1;""".format(station_number) 
        try:
            db.c.execute(query_tmp)
            [lat, lon, elevation_m, observer_name] = db.c.fetchone()
            observer_location = Topos(latitude_degrees = lat, longitude_degrees = lon, elevation_m = elevation_m)
            print("{}: lat: {}  lon: {}  {}".format(observer_name, lat, lon, elevation_m))
        except:
            log.warning("Station data not found for '{}' not found.".format(station_number))
            continue

        tm = obs_time_sql
        tmseconds = tm.second + (tm.microsecond/1000000)
        t = ts.utc(tm.year,tm.month,tm.day,tm.hour,tm.minute,tmseconds)

        satellite = EarthSatellite(tle_line_1, tle_line_2, name=sat_name)

        #astrometric = observer_location.at(t).observe(mars)
        #alt, az, d = astrometric.apparent().altaz()

        days = t - satellite.epoch
        difference = satellite - observer_location
        topocentric = difference.at(t)
        # print(topocentric.speed().km_per_s)
        # print(topocentric.position.km)
        
        # (dx,dy,dz) = topocentric.position.km
        # print(dx,dy,dz)

        print("{} using station {} Observing Satellite {} : {:.3f} days away from epoch".format(observer_name,station_number,satellite.model.satnum,days))
        ra_predict, dec_predict, distance = topocentric.radec()  # ICRF ("J2000")
        print("RA:  {:7.4f} Obs: {:.4f} (delta: {:.3f})".format(ra_predict._degrees, ra_obs, (ra_obs-ra_predict._degrees) ))
        print("DEC: {:7.4f} Obs: {:.4f} (delta: {:.3f})".format(dec_predict._degrees, decl_obs, (decl_obs-dec_predict._degrees) ))
        print("Predicted range: {}".format(distance.km))

        #     # An attempt at computing residuals, based on sattools residuals.c
        #     # rx0 = ra._degrees
        #     # ry0 = dec._degrees

        #     # drx = IOD_line.RA  - rx0
        #     # dry = IOD_line.DEC - ry0

        #     # dt = -(rx0*drx+ry0*dry)/(drx*drx*dry*dry)
        #     # dr = math.sqrt( ((dry*rx0-drx*ry0)**2.0)/(drx*drx+dry*dry) )

        #     # if ((-rx0*drx-ry0*dry)<0.0):
        #     #  	dr *=-1
        #     # print("Residuals: ",dr,dt)
        # else:
        #     alt, az, distance = topocentric.altaz()
        #     print(az.degrees, IOD_line.AZ)
        #     print(alt.degrees, IOD_line.EL)

if (__name__ == '__main__'):
    main()