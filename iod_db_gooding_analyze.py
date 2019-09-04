#!/usr/bin/env python3
from __future__ import division         # Eliminate need for decimals on whole values
import sys
import os # Only needed for paths while our repos are not public
import configargparse

if sys.version_info[0] != 3 or sys.version_info[1] < 7:
    print("This script requires Python version 3.7")
    sys.exit(1)

import os, re, csv
from skyfield.api import Topos, load, EarthSatellite, Star
from skyfield.elementslib import osculating_elements_of
from skyfield.units import Angle

import json
import time
from datetime import datetime
from math import radians, sqrt

# # odpdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# # datdir = os.path.join(odpdir, "examples", "data")
# sys.path.append('/Users/chris/Dropbox/code/orbdetpy-master/')
# from orbdetpy import _Conversion, iodGooding

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

    # For Alt/Az conversions
    planets = load('de421.bsp')
    earth = planets['earth']
    # iod_obs_id = input("IOD db ID: ")
    while(True):
        try:
            iod_obs_id = input("\nEnter 3 IOD Obs IDs: ")
            iod_obs_id = iod_obs_id.strip()
            iod_obs = iod_obs_id.split(' ')
        except:
            break

        if (len(iod_obs)==1):
            query_tmp = """SELECT obs_id, object_number, station_number, user_string, obs_time, ra, declination FROM ParsedIOD 
                WHERE obs_id >= {OBS_ID}
                ORDER BY obs_id ASC
                LIMIT 10;""".format(OBS_ID=iod_obs[0])
            db.c.execute(query_tmp)
            for r in db.c.fetchall():
                print("{} {:5} {:4} {} {} ra:{:<8.4f} dec:{:<8.4}".format(
                    r[0], r[1], r[2], r[3], 
                    r[4].isoformat(),
                    r[5],r[6]
                ))
            continue
        elif(len(iod_obs)!=3):
            break

        query_tmp = """SELECT obs_id, object_number, ra, declination, station_number, Obs.latitude, Obs.longitude, Obs.elevation_m, Obs.obs_name,  obs_time, iod_string
            FROM ParsedIOD
			JOIN (SELECT latitude, longitude, elevation_m, Observer.name as obs_name, Station.station_num as station_num FROM Station,Observer WHERE Station.user = Observer.id) Obs ON ParsedIOD.station_number = Obs.station_num
            WHERE obs_id={} OR obs_id={} OR obs_id={}
            ORDER by obs_time ASC;""".format(iod_obs[0], iod_obs[1], iod_obs[2]) 

        try:
            db.cdict.execute(query_tmp)
            observation_list = db.cdict.fetchall()
        except:
            log.warning("IOD ID '{}' not found.".format(iod_obs_id))
            continue
        # [object_number, station_number, obs_time_sql, ra_obs, decl_obs, iod_string]
        obs_id = []
        object_number = []
        station_number = []
        ra_obs = []
        decl_obs = []
        ra_obs_rad = []
        decl_obs_rad = []
        gslat = []
        gslon = []
        gsalt = []
        obs_name = []
        observer_location = []
        obs_time_sql = []
        iod_string = []
        alt_predicti = []
        az_predicti = []
        rhoi = []
        sfazi = []
        sfalti = []

        for i in range(0,3):
            # (obs_id[i], object_number[i], ra_obs[i], decl_obs[i], station_number[i], lat[i], lon[i], elev[i], obs_name[i], obs_time_sql[i]) = observation_list[i]
            obs_id.append( observation_list[i]['obs_id'] )
            object_number.append( observation_list[i]['object_number'] )
            ra_obs.append( observation_list[i]['ra'] ) 
            decl_obs.append( observation_list[i]['declination'] ) 
            ra_obs_rad.append( radians(observation_list[i]['ra'] )) 
            decl_obs_rad.append( radians(observation_list[i]['declination'] )) 
            station_number.append( observation_list[i]['station_number'] )
            gslat.append( radians(observation_list[i]['latitude']) )
            gslon.append( radians(observation_list[i]['longitude']) )
            gsalt.append( observation_list[i]['elevation_m'] )
            obs_name.append( observation_list[i]['obs_name'] )
            obs_time_sql.append( observation_list[i]['obs_time'] )
            iod_string.append( observation_list[i]['iod_string'] )
            observer_location.append( Topos(
                latitude_degrees  = observation_list[i]['latitude'], 
                longitude_degrees = observation_list[i]['longitude'], 
                elevation_m       = observation_list[i]['elevation_m']) )

        try:
            epoch_time = obs_time_sql[0].strftime('%Y-%m-%d %H:%M:%S')

            TLE = db.selectTLEEpochBeforeDate(obs_time_sql[0], object_number[0])
            sat_name   = TLE[1]
            tle_line_1 = TLE[2]
            tle_line_2 = TLE[3]
            print("{}\n{}\n{}".format(sat_name,tle_line_1,tle_line_2))
        except:
            log.warning("NO TLE found for object '{}' with EPOCH before {}".format(object_number[0], obs_time_sql[0]))
            continue

        # query_tmp = """SELECT latitude, longitude, elevation_m, Observer.name as name FROM Station
        #     LEFT JOIN Observer ON Station.user = Observer.id
        #     WHERE station_num={}
        #     LIMIT 1;""".format(station_number) 
        # try:
        #     db.c.execute(query_tmp)
        #     [lat, lon, elevation_m, observer_name] = db.c.fetchone()
        #     observer_location = Topos(latitude_degrees = lat, longitude_degrees = lon, elevation_m = elevation_m)
        #     print("{}: lat: {}  lon: {}  {}".format(observer_name, lat, lon, elevation_m))
        # except:
        #     log.warning("Station data not found for '{}' not found.".format(station_number))
        #     continue

        satellite = EarthSatellite(tle_line_1, tle_line_2, name=sat_name)
        for i in range(0,3):
            print("\nObs ID: {}".format(obs_id[i]))
            tm = obs_time_sql[i]
            tmseconds = tm.second + (tm.microsecond/100000)
            t = ts.utc(tm.year,tm.month,tm.day,tm.hour,tm.minute,tmseconds)

            #astrometric = observer_location.at(t).observe(mars)
            #alt, az, d = astrometric.apparent().altaz()

            days = t - satellite.epoch
            difference = satellite - observer_location[i]
            topocentric = difference.at(t)
            oe = osculating_elements_of(satellite.at(t))
            # print(topocentric.speed().km_per_s)
            # print(topocentric.position.km)
            
            # (dx,dy,dz) = topocentric.position.km
            # print(dx,dy,dz)

            fake_sat = Star(ra=Angle(degrees=ra_obs[i]),
               dec=Angle(degrees=decl_obs[i]))
            earth_station = earth + observer_location[i]
            astrometric = earth_station.at(t).observe(fake_sat)
            (sfalt, sfaz, _) = astrometric.apparent().altaz(temperature_C=-273.25, pressure_mbar=0)
            sfalti.append(sfalt.radians)
            sfazi.append(sfaz.radians)

            if (i == 0):
                print("Time: {}".format(observation_list[i]['obs_time'].isoformat()))
            else:
                delta = observation_list[i]['obs_time'] - observation_list[i-1]['obs_time']
                print("Time delta: {:5.1f}s from prev".format(delta.total_seconds()))
            print("{} using station {} Observing Satellite {} : {:.3f} days away from epoch".format(obs_name[i],station_number[i],satellite.model.satnum,days))
            ra_predict, dec_predict, distance = topocentric.radec() # ICRF ("J2000")
            (alt_predict, az_predict, _) = topocentric.altaz()      # ICRF ("J2000")

            alt_predicti.append(alt_predict.radians)
            az_predicti.append(az_predict.radians)
            rhoi.append(distance.m)

            ra_error   =      ra_obs[i] - ra_predict._degrees
            decl_error =    decl_obs[i] - dec_predict._degrees
            az_error   =  sfaz._degrees - az_predict._degrees
            alt_error  = sfalt._degrees - alt_predict._degrees

            ra_dec_error = sqrt (ra_error*ra_error + decl_error*decl_error)
            alt_az_error = sqrt (az_error*az_error + alt_error*alt_error)

            x_radec_err = distance.km * ra_dec_error           
            x_altaz_err = distance.km * alt_az_error

            print("TLE RA:  {:8.4f} | Obs RA:  {:8.4f} | delta: {:6.3f} | total:   {:2.3f} deg".format(ra_predict._degrees, ra_obs[i], ra_error, ra_dec_error ))
            print("TLE DEC: {:8.4f} | Obs DEC: {:8.4f} | delta: {:6.3f} | pos err: {:2.3f} km".format(dec_predict._degrees, decl_obs[i], decl_error, x_radec_err ))
            print("TLE AZ:  {:8.4f} | Obs AZ:  {:8.4f} | delta: {:6.3f} | total:   {:2.3f} deg".format(az_predict._degrees, sfaz._degrees, az_error, alt_az_error ))
            print("TLE Alt: {:8.4f} | Obs Alt: {:8.4f} | delta: {:6.3f} | pos err: {:2.3f} km".format(alt_predict._degrees, sfalt._degrees, alt_error, x_altaz_err ))
            # print("TLE AZ:  {:8.4f} | Obs AZ:  {:8.4f} (radians)".format(az_predict.radians, sfaz.radians))
            # print("TLE Alt: {:8.4f} | Obs Alt: {:8.4f} (radians)".format(alt_predict.radians, sfalt.radians))
            print("TLE Range: {:8.4f} km".format(distance.km))
            # print("Predicted: {:8.4f} Az  {:8.4f} Alt (degrees)".format(az_predict._degrees, alt_predict._degrees))
            # print("Predicted: {:8.4f} Az  {:8.4f} Alt (radians)".format(az_predict.radians, alt_predict.radians))
            ([x, y, z], [xdot, ydot, zdot], _) = satellite._position_and_velocity_TEME_km(t)
            print(x*1000, y*1000, z*1000, xdot*1000, ydot*1000, zdot*1000)



# tmstr[i], ra[i], dec[i], gslat[i], gslon[i], gsalt[i])
# iodGooding(gslat, gslon, gsalt, tmstr, azi, eli, rho1init, rho3init)

        print("\n# {}".format(iod_obs_id))
        print("ra =    ", ra_obs_rad)
        print("dec =   ", decl_obs_rad)
        print("azi =   ", sfazi)
        print("eli =   ", sfalti)
        print("gslat = ", gslat)
        print("gslon = ", gslon)
        print("gsalt = ", gsalt)
        print("rho1init = ", rhoi[0])
        print("rho3init = ", rhoi[2])

        fmt = '%Y-%m-%dT%H:%M:%S.%fZ'
        print("tmstr =     ", [obs_time_sql[0].strftime(fmt), obs_time_sql[1].strftime(fmt), obs_time_sql[2].strftime(fmt) ])
        print()
        print("Az Predict  ", az_predicti)
        print("Alt Predict ", alt_predicti)

        print("\nSatellite:  {} {}".format(satellite.model.satnum, sat_name))
        print("a:        {:10.2f}".format(oe.semi_major_axis.m))
        print("e: ecc        {:6.2f}".format(oe.eccentricity))
        print("i: inc        {:6.2f}".format(oe.inclination._degrees))
        print("pa: AOP       {:6.2f}".format(oe.argument_of_periapsis._degrees))
        print("RAAN          {:6.2f}".format(oe.longitude_of_ascending_node._degrees))
        print("v: True Anom  {:6.2f}".format(oe.true_anomaly._degrees))
        print("period (min.) {:6.2f}".format(oe.period_in_days*24*60))

        for i in range(0,3):
            print(iod_string[i])

if (__name__ == '__main__'):
    main()