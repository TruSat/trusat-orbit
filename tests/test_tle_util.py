#!/usr/bin/env python

import os
import sys
from math import degrees, radians, pi, fmod, sqrt, pow
import copy
import json
from datetime import datetime

try:
    from unittest2 import TestCase, main, skipUnless
except:
    from unittest import TestCase, main, skipUnless

# python SGP4 from git+https://github.com/interplanetarychris/python-sgp4@cython-7-dec-15-vallado
# Until the following pull request is approved
# https://github.com/brandon-rhodes/python-sgp4/pull/35

# sys.path.insert(1, '../trusat-backend')

try:
    from trusat.caccelerated import *
except ImportError as e:
    print(e)
    from sgp4.propagation import sgp4, sgp4init
from sgp4.ext import invjday, days2mdhms, rv2coe
from sgp4.io import twoline2rv
from sgp4.model import Satellite
from sgp4 import earth_gravity

try:
    from trusat_backend import database
    DB = True
    print("Will attempt live database tests.")
except ImportError:
    DB = False
    print("DB module not available")


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) 

import trusat.tle_util as tle_util  
import trusat.satfit as satfit

# Global variables
TWOPI = 2*pi
NOCON = TWOPI/1440.0

line0 = "USA 186"
line1 = "1 28888U 05042A   14063.83505828 0.00032100  00000-0  29747-3 0    09"
line2 = "2 28888  96.8945 115.9842 0499287 308.6206  51.3792 14.85742003    02"
# Original TLE Epoch data points
# 2014-03-04T20:02:29.035392
# 23439.8350582798
# 2456721.33505828

def serialize(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, datetime):
        return obj.isoformat()

    return obj.__dict__

def compare_satrecs(sata,satb):
    # A most unholy quick hack until we can parse Class attributes more expertly
    if (sata.line0 != satb.line0):
        print("{} : A {} != B {}".format('line0',sata.line0,satb.line0))
    if (sata.satnum != satb.satnum):
        print("{} : A {} != B {}".format('satnum',sata.satnum,satb.satnum))
    if (sata.line1 != satb.line1):
        print("{} : A {} != B {}".format('line1',sata.line1,satb.line1))
    if (sata.classification != satb.classification):
        print("{} : A {} != B {}".format('classification',sata.classification,satb.classification))
    if (sata.intldesg != satb.intldesg):
        print("{} : A {} != B {}".format('intldesg',sata.intldesg,satb.intldesg))
    if (sata.epochyr != satb.epochyr):
        print("{} : A {} != B {}".format('epochyr',sata.epochyr,satb.epochyr))
    if (sata.epochdays != satb.epochdays):
        print("{} : A {} != B {}".format('epochdays',sata.epochdays,satb.epochdays))
    if (sata.ndot != satb.ndot):
        print("{} : A {} != B {}".format('ndot',sata.ndot,satb.ndot))
    if (sata.nddot != satb.nddot):
        print("{} : A {} != B {}".format('nddot',sata.nddot,satb.nddot))
    if (sata.bstar != satb.bstar):
        print("{} : A {} != B {}".format('bstar',sata.bstar,satb.bstar))
    if (sata.ephtype != satb.ephtype):
        print("{} : A {} != B {}".format('ephtype',sata.ephtype,satb.ephtype))
    if (sata.elnum != satb.elnum):
        print("{} : A {} != B {}".format('elnum',sata.elnum,satb.elnum))
    if (sata.line2 != satb.line2):
        print("{} : A {} != B {}".format('line2',sata.line2,satb.line2))
    if (sata.inclo != satb.inclo):
        print("{} : A {} != B {}".format('inclo',sata.inclo,satb.inclo))
    if (sata.nodeo != satb.nodeo):
        print("{} : A {} != B {}".format('nodeo',sata.nodeo,satb.nodeo))
    if (sata.ecco != satb.ecco):
        print("{} : A {} != B {}".format('ecco',sata.ecco,satb.ecco))
    if (sata.argpo != satb.argpo):
        print("{} : A {} != B {}".format('argpo',sata.argpo,satb.argpo))
    if (sata.mo != satb.mo):
        print("{} : A {} != B {}".format('mo',sata.mo,satb.mo))
    if (sata.no_kozai != satb.no_kozai):
        print("{} : A {} != B {}".format('no_kozai',sata.no_kozai,satb.no_kozai))
    if (sata.revnum != satb.revnum):
        print("{} : A {} != B {}".format('revnum',sata.revnum,satb.revnum))
    if (sata.jdsatepoch != satb.jdsatepoch):
        print("{} : A {} != B {}".format('jdsatepoch',sata.jdsatepoch,satb.jdsatepoch))
    if (sata.jdSGP4epoch != satb.jdSGP4epoch):
        print("{} : A {} != B {}".format('jdSGP4epoch',sata.jdSGP4epoch,satb.jdSGP4epoch))
    if (sata.epoch_datetime != satb.epoch_datetime):
        print("{} : A {} != B {}".format('epoch_datetime',sata.epoch_datetime,satb.epoch_datetime))
    if (sata.parent_tle_id != satb.parent_tle_id):
        print("{} : A {} != B {}".format('parent_tle_id',sata.parent_tle_id,satb.parent_tle_id))
    if (sata.operationmode != satb.operationmode):
        print("{} : A {} != B {}".format('operationmode',sata.operationmode,satb.operationmode))
    if (sata.error != satb.error):
        print("{} : A {} != B {}".format('error',sata.error,satb.error))
    if (sata.whichconst != satb.whichconst):
        print("{} : A {} != B {}".format('whichconst',sata.whichconst,satb.whichconst))
    if (sata.isimp != satb.isimp):
        print("{} : A {} != B {}".format('isimp',sata.isimp,satb.isimp))
    if (sata.method != satb.method):
        print("{} : A {} != B {}".format('method',sata.method,satb.method))
    if (sata.aycof != satb.aycof):
        print("{} : A {} != B {}".format('aycof',sata.aycof,satb.aycof))
    if (sata.con41 != satb.con41):
        print("{} : A {} != B {}".format('con41',sata.con41,satb.con41))
    if (sata.cc1 != satb.cc1):
        print("{} : A {} != B {}".format('cc1',sata.cc1,satb.cc1))
    if (sata.cc4 != satb.cc4):
        print("{} : A {} != B {}".format('cc4',sata.cc4,satb.cc4))
    if (sata.cc5 != satb.cc5):
        print("{} : A {} != B {}".format('cc5',sata.cc5,satb.cc5))
    if (sata.d2 != satb.d2):
        print("{} : A {} != B {}".format('d2',sata.d2,satb.d2))
    if (sata.d3 != satb.d3):
        print("{} : A {} != B {}".format('d3',sata.d3,satb.d3))
    if (sata.d4 != satb.d4):
        print("{} : A {} != B {}".format('d4',sata.d4,satb.d4))
    if (sata.delmo != satb.delmo):
        print("{} : A {} != B {}".format('delmo',sata.delmo,satb.delmo))
    if (sata.eta != satb.eta):
        print("{} : A {} != B {}".format('eta',sata.eta,satb.eta))
    if (sata.argpdot != satb.argpdot):
        print("{} : A {} != B {}".format('argpdot',sata.argpdot,satb.argpdot))
    if (sata.omgcof != satb.omgcof):
        print("{} : A {} != B {}".format('omgcof',sata.omgcof,satb.omgcof))
    if (sata.sinmao != satb.sinmao):
        print("{} : A {} != B {}".format('sinmao',sata.sinmao,satb.sinmao))
    if (sata.t != satb.t):
        print("{} : A {} != B {}".format('t',sata.t,satb.t))
    if (sata.t2cof != satb.t2cof):
        print("{} : A {} != B {}".format('t2cof',sata.t2cof,satb.t2cof))
    if (sata.t3cof != satb.t3cof):
        print("{} : A {} != B {}".format('t3cof',sata.t3cof,satb.t3cof))
    if (sata.t4cof != satb.t4cof):
        print("{} : A {} != B {}".format('t4cof',sata.t4cof,satb.t4cof))
    if (sata.t5cof != satb.t5cof):
        print("{} : A {} != B {}".format('t5cof',sata.t5cof,satb.t5cof))
    if (sata.x1mth2 != satb.x1mth2):
        print("{} : A {} != B {}".format('x1mth2',sata.x1mth2,satb.x1mth2))
    if (sata.x7thm1 != satb.x7thm1):
        print("{} : A {} != B {}".format('x7thm1',sata.x7thm1,satb.x7thm1))
    if (sata.mdot != satb.mdot):
        print("{} : A {} != B {}".format('mdot',sata.mdot,satb.mdot))
    if (sata.nodedot != satb.nodedot):
        print("{} : A {} != B {}".format('nodedot',sata.nodedot,satb.nodedot))
    if (sata.xlcof != satb.xlcof):
        print("{} : A {} != B {}".format('xlcof',sata.xlcof,satb.xlcof))
    if (sata.xmcof != satb.xmcof):
        print("{} : A {} != B {}".format('xmcof',sata.xmcof,satb.xmcof))
    if (sata.nodecf != satb.nodecf):
        print("{} : A {} != B {}".format('nodecf',sata.nodecf,satb.nodecf))
    if (sata.irez != satb.irez):
        print("{} : A {} != B {}".format('irez',sata.irez,satb.irez))
    if (sata.d2201 != satb.d2201):
        print("{} : A {} != B {}".format('d2201',sata.d2201,satb.d2201))
    if (sata.d2211 != satb.d2211):
        print("{} : A {} != B {}".format('d2211',sata.d2211,satb.d2211))
    if (sata.d3210 != satb.d3210):
        print("{} : A {} != B {}".format('d3210',sata.d3210,satb.d3210))
    if (sata.d3222 != satb.d3222):
        print("{} : A {} != B {}".format('d3222',sata.d3222,satb.d3222))
    if (sata.d4410 != satb.d4410):
        print("{} : A {} != B {}".format('d4410',sata.d4410,satb.d4410))
    if (sata.d4422 != satb.d4422):
        print("{} : A {} != B {}".format('d4422',sata.d4422,satb.d4422))
    if (sata.d5220 != satb.d5220):
        print("{} : A {} != B {}".format('d5220',sata.d5220,satb.d5220))
    if (sata.d5232 != satb.d5232):
        print("{} : A {} != B {}".format('d5232',sata.d5232,satb.d5232))
    if (sata.d5421 != satb.d5421):
        print("{} : A {} != B {}".format('d5421',sata.d5421,satb.d5421))
    if (sata.d5433 != satb.d5433):
        print("{} : A {} != B {}".format('d5433',sata.d5433,satb.d5433))
    if (sata.dedt != satb.dedt):
        print("{} : A {} != B {}".format('dedt',sata.dedt,satb.dedt))
    if (sata.del1 != satb.del1):
        print("{} : A {} != B {}".format('del1',sata.del1,satb.del1))
    if (sata.del2 != satb.del2):
        print("{} : A {} != B {}".format('del2',sata.del2,satb.del2))
    if (sata.del3 != satb.del3):
        print("{} : A {} != B {}".format('del3',sata.del3,satb.del3))
    if (sata.didt != satb.didt):
        print("{} : A {} != B {}".format('didt',sata.didt,satb.didt))
    if (sata.dmdt != satb.dmdt):
        print("{} : A {} != B {}".format('dmdt',sata.dmdt,satb.dmdt))
    if (sata.dnodt != satb.dnodt):
        print("{} : A {} != B {}".format('dnodt',sata.dnodt,satb.dnodt))
    if (sata.domdt != satb.domdt):
        print("{} : A {} != B {}".format('domdt',sata.domdt,satb.domdt))
    if (sata.e3 != satb.e3):
        print("{} : A {} != B {}".format('e3',sata.e3,satb.e3))
    if (sata.ee2 != satb.ee2):
        print("{} : A {} != B {}".format('ee2',sata.ee2,satb.ee2))
    if (sata.peo != satb.peo):
        print("{} : A {} != B {}".format('peo',sata.peo,satb.peo))
    if (sata.pgho != satb.pgho):
        print("{} : A {} != B {}".format('pgho',sata.pgho,satb.pgho))
    if (sata.pho != satb.pho):
        print("{} : A {} != B {}".format('pho',sata.pho,satb.pho))
    if (sata.pinco != satb.pinco):
        print("{} : A {} != B {}".format('pinco',sata.pinco,satb.pinco))
    if (sata.plo != satb.plo):
        print("{} : A {} != B {}".format('plo',sata.plo,satb.plo))
    if (sata.se2 != satb.se2):
        print("{} : A {} != B {}".format('se2',sata.se2,satb.se2))
    if (sata.se3 != satb.se3):
        print("{} : A {} != B {}".format('se3',sata.se3,satb.se3))
    if (sata.sgh2 != satb.sgh2):
        print("{} : A {} != B {}".format('sgh2',sata.sgh2,satb.sgh2))
    if (sata.sgh3 != satb.sgh3):
        print("{} : A {} != B {}".format('sgh3',sata.sgh3,satb.sgh3))
    if (sata.sgh4 != satb.sgh4):
        print("{} : A {} != B {}".format('sgh4',sata.sgh4,satb.sgh4))
    if (sata.sh2 != satb.sh2):
        print("{} : A {} != B {}".format('sh2',sata.sh2,satb.sh2))
    if (sata.sh3 != satb.sh3):
        print("{} : A {} != B {}".format('sh3',sata.sh3,satb.sh3))
    if (sata.si2 != satb.si2):
        print("{} : A {} != B {}".format('si2',sata.si2,satb.si2))
    if (sata.si3 != satb.si3):
        print("{} : A {} != B {}".format('si3',sata.si3,satb.si3))
    if (sata.sl2 != satb.sl2):
        print("{} : A {} != B {}".format('sl2',sata.sl2,satb.sl2))
    if (sata.sl3 != satb.sl3):
        print("{} : A {} != B {}".format('sl3',sata.sl3,satb.sl3))
    if (sata.sl4 != satb.sl4):
        print("{} : A {} != B {}".format('sl4',sata.sl4,satb.sl4))
    if (sata.gsto != satb.gsto):
        print("{} : A {} != B {}".format('gsto',sata.gsto,satb.gsto))
    if (sata.xfact != satb.xfact):
        print("{} : A {} != B {}".format('xfact',sata.xfact,satb.xfact))
    if (sata.xgh2 != satb.xgh2):
        print("{} : A {} != B {}".format('xgh2',sata.xgh2,satb.xgh2))
    if (sata.xgh3 != satb.xgh3):
        print("{} : A {} != B {}".format('xgh3',sata.xgh3,satb.xgh3))
    if (sata.xgh4 != satb.xgh4):
        print("{} : A {} != B {}".format('xgh4',sata.xgh4,satb.xgh4))
    if (sata.xh2 != satb.xh2):
        print("{} : A {} != B {}".format('xh2',sata.xh2,satb.xh2))
    if (sata.xh3 != satb.xh3):
        print("{} : A {} != B {}".format('xh3',sata.xh3,satb.xh3))
    if (sata.xi2 != satb.xi2):
        print("{} : A {} != B {}".format('xi2',sata.xi2,satb.xi2))
    if (sata.xi3 != satb.xi3):
        print("{} : A {} != B {}".format('xi3',sata.xi3,satb.xi3))
    if (sata.xl2 != satb.xl2):
        print("{} : A {} != B {}".format('xl2',sata.xl2,satb.xl2))
    if (sata.xl3 != satb.xl3):
        print("{} : A {} != B {}".format('xl3',sata.xl3,satb.xl3))
    if (sata.xl4 != satb.xl4):
        print("{} : A {} != B {}".format('xl4',sata.xl4,satb.xl4))
    if (sata.xlamo != satb.xlamo):
        print("{} : A {} != B {}".format('xlamo',sata.xlamo,satb.xlamo))
    if (sata.zmol != satb.zmol):
        print("{} : A {} != B {}".format('zmol',sata.zmol,satb.zmol))
    if (sata.zmos != satb.zmos):
        print("{} : A {} != B {}".format('zmos',sata.zmos,satb.zmos))
    if (sata.atime != satb.atime):
        print("{} : A {} != B {}".format('atime',sata.atime,satb.atime))
    if (sata.xli != satb.xli):
        print("{} : A {} != B {}".format('xli',sata.xli,satb.xli))
    if (sata.xni != satb.xni):
        print("{} : A {} != B {}".format('xni',sata.xni,satb.xni))
    if (sata.tumin != satb.tumin):
        print("{} : A {} != B {}".format('tumin',sata.tumin,satb.tumin))
    if (sata.mu != satb.mu):
        print("{} : A {} != B {}".format('mu',sata.mu,satb.mu))
    if (sata.radiusearthkm != satb.radiusearthkm):
        print("{} : A {} != B {}".format('radiusearthkm',sata.radiusearthkm,satb.radiusearthkm))
    if (sata.xke != satb.xke):
        print("{} : A {} != B {}".format('xke',sata.xke,satb.xke))
    if (sata.j2 != satb.j2):
        print("{} : A {} != B {}".format('j2',sata.j2,satb.j2))
    if (sata.j3 != satb.j3):
        print("{} : A {} != B {}".format('j3',sata.j3,satb.j3))
    if (sata.j4 != satb.j4):
        print("{} : A {} != B {}".format('j4',sata.j4,satb.j4))
    if (sata.j3oj2 != satb.j3oj2):
        print("{} : A {} != B {}".format('j3oj2',sata.j3oj2,satb.j3oj2))
    if (sata.am != satb.am):
        print("{} : A {} != B {}".format('am',sata.am,satb.am))
    if (sata.em != satb.em):
        print("{} : A {} != B {}".format('em',sata.em,satb.em))
    if (sata.im != satb.im):
        print("{} : A {} != B {}".format('im',sata.im,satb.im))
    if (sata.Om != satb.Om):
        print("{} : A {} != B {}".format('Om',sata.Om,satb.Om))
    if (sata.mm != satb.mm):
        print("{} : A {} != B {}".format('mm',sata.mm,satb.mm))
    if (sata.nm != satb.nm):
        print("{} : A {} != B {}".format('nm',sata.nm,satb.nm))
    if (sata.init != satb.init):
        print("{} : A {} != B {}".format('init',sata.init,satb.init))
    if (sata.no_unkozai != satb.no_unkozai):
        print("{} : A {} != B {}".format('no_unkozai',sata.no_unkozai,satb.no_unkozai))
    if (sata.a != satb.a):
        print("{} : A {} != B {}".format('a',sata.a,satb.a))
    if (sata.alta != satb.alta):
        print("{} : A {} != B {}".format('alta',sata.alta,satb.alta))
    if (sata.altp != satb.altp):
        print("{} : A {} != B {}".format('altp',sata.altp,satb.altp))
    if (sata.error_message != satb.error_message):
        print("{} : A {} != B {}".format('error_message',sata.error_message,satb.error_message))
    if (sata.om != satb.om):
        print("{} : A {} != B {}".format('om',sata.om,satb.om))
    # if (sata.rr_km != satb.rr_km):
    #     print("{} : A {} != B {}".format('rr_km',sata.rr_km,satb.rr_km))
    # if (sata.vv_kmpersec != satb.vv_kmpersec):
    #     print("{} : A {} != B {}".format('vv_kmpersec',sata.vv_kmpersec,satb.vv_kmpersec))
    # if (sata.rr != satb.rr):
    #     print("{} : A {} != B {}".format('rr',sata.rr,satb.rr))
    # if (sata.vv != satb.vv):
    #     print("{} : A {} != B {}".format('vv',sata.vv,satb.vv))



class Tests(TestCase):
    # Temporary database credentials hack
    if (DB):
        try:
            CONFIG = os.path.abspath("../trusat-config.yaml")
            db = database.Database(CONFIG)
        except: 
            log.error("DB Login credentials not available.")
            db = False
    else:
        db = False

    # TLE Epoch Test
    TestCase.epoch_string = '19236.07927356'
    TestCase.maxDiff = None
    
    def compare_satrecs_as_JSON(self, gotSat, expectedSat):
        # !TODO: there may be parts of the TruSatellite objects that are expected to differ
        # after a database round trip (e.g. tle_id, import_timestamp).
        # If so, these can be munged out here before comparison.
        gotSat.tle_id = None
        gotSat.import_timestamp = None
        self.assertEqual(
            json.dumps(gotSat, default=serialize, indent=4, sort_keys=True),
            json.dumps(expectedSat, default=serialize, indent=4, sort_keys=True)
        )

    def test_datetime_from_tle_fmt(self):
        TestCase.epoch_datetime = tle_util.datetime_from_tle_fmt(TestCase.epoch_string)

        self.assertEqual(TestCase.epoch_datetime.year,2019)
        self.assertEqual(TestCase.epoch_datetime.timetuple().tm_yday,236) 
        self.assertEqual(TestCase.epoch_datetime.month,8)
        self.assertEqual(TestCase.epoch_datetime.day,24)
        self.assertEqual(TestCase.epoch_datetime.hour,1)
        self.assertEqual(TestCase.epoch_datetime.minute,54)
        self.assertEqual(TestCase.epoch_datetime.second,9)
        self.assertEqual(TestCase.epoch_datetime.microsecond,235584)


    def test_tle_fmt_epoch(self):
        epoch_string_test = tle_util.tle_fmt_epoch(TestCase.epoch_datetime)
        self.assertEqual(epoch_string_test,TestCase.epoch_string)


    def assert_expected_TLE(self, TLE):
        """ Asserts that the TLE contents are as expected for the two tests below.

            Also calls satfit.initsat and checks that the result is as expected.
        """
        # Verify the data in the TLE made it exactly to the TLE object
        self.assertEqual(TLE.satellite_number,28888)
        self.assertEqual(TLE.classification,'U')
        self.assertEqual(TLE.designation,'2005-042A')
        self.assertEqual(TLE.mean_motion_derivative,0.00032100)
        self.assertEqual(TLE.mean_motion_sec_derivative,0.0)
        self.assertEqual(TLE.bstar,0.29747e-3)
        self.assertEqual(TLE.ephemeris_type,0)
        self.assertEqual(TLE.element_set_number,0)
        self.assertEqual(TLE.inclination_degrees,96.8945)
        self.assertEqual(TLE.raan_degrees,115.9842)
        self.assertEqual(TLE.eccentricity,0.0499287)
        self.assertEqual(TLE.arg_perigee_degrees,308.6206)
        self.assertEqual(TLE.mean_anomaly_degrees,51.3792)
        self.assertEqual(TLE.mean_motion_orbits_per_day,14.85742003)
        self.assertEqual(TLE.orbit_number,0)
        self.assertEqual(TLE.name,'USA 186')       

        # Other derived quantities
        self.assertEqual(tle_util.tle_fmt_epoch(TLE.epoch_datetime),'14063.83505828')
        # TODO: Add test for jdsatepoch and jdSGP4epoch
        # FIXME: TLE object currently has no "mean motion degrees" variable
        self.assertEqual(TLE.mean_motion_radians_per_minute,TLE.mean_motion_orbits_per_day*NOCON)
        self.assertEqual(TLE.mean_motion_radians_per_second,TLE.mean_motion_orbits_per_day*NOCON/60.0)
        self.assertEqual(TLE.inclination_degrees,degrees(TLE.inclination_radians))
        self.assertEqual(TLE.raan_degrees,degrees(TLE.raan_radians))
        self.assertEqual(TLE.arg_perigee_degrees,degrees(TLE.arg_perigee_radians))
        self.assertEqual(TLE.mean_anomaly_degrees,degrees(TLE.mean_anomaly_radians))
        self.assertEqual(TLE.launch_piece_number,tle_util.launch_piece_letter_to_number(TLE._id_launch_piece_letter))
        self.assertEqual(TLE._id_launch_year,2005)
        self.assertEqual(TLE._id_launch_num,42)
        self.assertEqual(TLE.line0,line0)
        self.assertEqual(TLE.line1,line1)
        self.assertEqual(TLE.line2,line2)

        # Add tests for perigee, period, apogee, semi_major_axis

        (sat1, sat1meta) = satfit.initsat(TLE)

        # Verify the TLE data is exactly translated to the initialized SGP4 satrec object
        self.assertEqual(tle_util.tle_fmt_epoch(sat1meta.epoch_datetime),'14063.83505828')
        self.assertEqual(TLE.epoch_datetime, sat1meta.epoch_datetime)
        self.assertEqual(TLE.eccentricity,sat1.ecco)
        self.assertEqual(TLE.inclination_radians, sat1.inclo)
        self.assertEqual(TLE.inclination_degrees, degrees(sat1.inclo))
        self.assertEqual(TLE.raan_radians, sat1.nodeo)
        self.assertEqual(TLE.raan_degrees, degrees(sat1.nodeo))
        self.assertEqual(TLE.arg_perigee_radians, sat1.argpo)
        self.assertEqual(TLE.arg_perigee_degrees, degrees(sat1.argpo))
        self.assertEqual(TLE.mean_anomaly_radians, sat1.mo)
        self.assertEqual(TLE.mean_anomaly_degrees, degrees(sat1.mo))
        self.assertEqual(TLE.mean_motion_radians_per_minute, sat1.no_kozai)
        self.assertEqual(TLE.mean_motion_orbits_per_day, sat1.no_kozai/NOCON)
        self.assertEqual(TLE.mean_motion_radians_per_second, sat1.no_kozai/60)  


    def test_satrec_to_TLE(self):
        TLE = tle_util.TruSatellite(line0=line0, line1=line1, line2=line2)

        # Verify the data in the TLE made it exactly to the TLE object
        self.assert_expected_TLE(TLE) 


    @skipUnless(db, "No database configured; skipping")
    def test_satrec_to_TLE_persisted(self):
        TLE = tle_util.TruSatellite(line0=line0, line1=line1, line2=line2)

        # Persist this to our database and get it back out, to make sure we do not lose or change data
        try:
            self.db.addTLE(TLE)
            TLE_ID = self.db.write_TLEs_to_db()  
            TLE_from_database = self.db.get_TLE(TLE_ID)

            # Verify the data in the TLE made it exactly to the TLE object
            self.compare_satrecs_as_JSON(TLE_from_database, TLE)
            self.assert_expected_TLE(TLE_from_database) 

        finally:
            self.db.clean()

        def test_tle_decimal(self):
            decm_assumed = "1234500"
            decm = -0.12345
            decm_pack = "-12345-6"

            dec_num = tle_util.read_tle_decimal(decm_pack)
            self.assertEqual(dec_num,decm_pack,msg="read_tle_decimal() failed")

            decm_pack_test = tle_util.tle_fmt_decimal_pack(dec_num)
            self.assertRegex(decm_pack_test,decm_pack,"tle_fmt_decimal_pack() failed")

            decm_assumed_test = tle_util.assumed_decimal_point(decm)
            self.assertEqual(decm_assumed_test,decm_assumed,msg="assumed_decimal_point() failed")


# Leftover test exploration - figure out if a unittest is useful here
def test_epoch_methods():
    line00 = line0 + " Apogee"
    line1  =  "1 28888U 05042A   14063.83505828 0.00032100  00000-0  29747-3 0    09"
    # 3 iters
    # line10 = '1 28888U 05042A   14063.82545230 0.00032100  00000-0  29747-3 0    09'
    # 30 iters
    line10  = '1 28888U 05042A   14063.82545230 0.00032096  00000-0  29747-3 0    03' # satfit.cpp
    line100 = '1 28888T 05042A   14063.82545230 0.00032100  00000-0  29747-3 0    09' # Mine
            
    line2  =  "2 28888  96.8945 115.9842 0499287 308.6206  51.3792 14.85742003    02"
    # 3 iters
    #line20 = '2 28888  96.8945 115.9758 0499291 308.6530   0.0000 14.85741387    00'
    # 30 iters
    line20  = '2 28888  96.8945 115.9758 0499282 308.6530   0.0000 14.85741387    00' # satfit.cpp
    line200 = '2 28888  96.8945 115.9758 0499293 308.6530 360.0000 14.85742003    07' # mine


    # sat0 : Elements regressed to perigee by satfit.cpp TLE0 epoch
    TLE0 = tle_util.TruSatellite(line0=line00, line1=line10, line2=line20)
    sat0 = satfit.initsat(TLE0)
    sat0 = satfit.delta_t(sat0,sat0.jdsatepoch)

    # sat2 : Elements regressed to perigee by satfit.py TLE0 epoch
    TLE2 = tle_util.TruSatellite(line0=line00, line1=line100, line2=line200)
    sat2 = satfit.initsat(TLE2)
    sat2 = satfit.delta_t(sat2,sat0.jdsatepoch)

    # "satx: Original elements back-propagated to perigee (epoch of sat0)"
    TLE = tle_util.TruSatellite(line0=line0, line1=line1, line2=line2)
    sat1 = satfit.initsat(TLE)
    sat1 = satfit.delta_t(sat1,sat0.jdsatepoch) # Changing tsince changes the mean elements

    # satA : Elements created from sat0 mean elements after propagating to TLE0 epoch
    satA = copy.deepcopy(sat0)    
    satA = satfit.delta_el(satA,xincl=satA.im,xnodeo=satA.Om,ec=satA.em,omegao=satA.om,xmo=satA.mm,xno=satA.no_kozai,bsr=sat1.bstar, jdsatepoch=sat0.jdsatepoch)
    satA = satfit.delta_t(satA,sat0.jdsatepoch)

    # satB : Elements created from sat1 mean elements after propagating to TLE0 epoch
    satB = copy.deepcopy(sat1)    
    satB = satfit.delta_el(satB,xincl=satB.im,xnodeo=satB.Om,ec=satB.em,omegao=satB.om,xmo=satB.mm,xno=satB.no_kozai,bsr=sat1.bstar, jdsatepoch=sat0.jdsatepoch)
    # satB.jdsatepoch = sat0.jdsatepoch 
    # satB.jdSGP4epoch = satB.jdsatepoch - 2433281.5
    # satB.epoch_datetime = satfit.jday_to_datetime(satB.jdsatepoch)
    satB = satfit.delta_t(satB,sat0.jdsatepoch)


    # satr : Elements created from rv2coe from sat1x r,v vectors TLE0 epoch
    satr = copy.deepcopy(sat1)    
    (p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper) = rv2coe(satr.rr_km, satr.vv_kmpersec, satr.whichconst.mu)
    # satx = satfit.delta_el(sat1,xincl=sat1.inclo,xnodeo=sat1.nodeo,ec=sat1.ecco,omegao=sat1.argpo,xmo=sat1.mo,xno=sat1.no_kozai,bsr=sat1.bstar)
    mean_motion = sqrt(satr.whichconst.mu/(pow(a,3)))*60.0 # radians per minute
    satr = satfit.delta_el(satr,xincl=incl,xnodeo=omega,ec=ecc,omegao=argp,xmo=m,xno=mean_motion,bsr=satr.bstar, jdsatepoch=sat0.jdsatepoch)
    # satr.jdsatepoch = sat0.jdsatepoch 
    # satr.jdSGP4epoch = satr.jdsatepoch - 2433281.5
    # satr.epoch_datetime = satfit.jday_to_datetime(satr.jdsatepoch)
    satr = satfit.delta_t(satr,sat0.jdsatepoch)


    print()
    print("sat1: Original elements back-propagated to perigee (epoch of sat0)")
    print("sat2 : Elements regressed to perigee by satfit.py TLE0 epoch")
    print("sat0 : Elements regressed to perigee by satfit.cpp TLE0 epoch")
    print("satA : Elements created from sat0 mean elements after propagating to TLE0 epoch")
    print("satB : Elements created from sat1 mean elements after propagating to TLE0 epoch")
    print("satr : Elements created from rv2coe from sat1x r,v vectors TLE0 epoch")
    print()
    print("At TLE epoch jd {}".format(sat0.jdsatepoch))
    print("sat1.rr {} sat1.vv {}".format(sat1.rr_km,sat1.vv_kmpersec))
    print("sat2.rr {} sat2.vv {}".format(sat2.rr_km,sat2.vv_kmpersec))
    print("sat0.rr {} sat0.vv {}".format(sat0.rr_km,sat0.vv_kmpersec))
    print("satA.rr {} satA.vv {}".format(satA.rr_km,satA.vv_kmpersec))
    print("satB.rr {} satB.vv {}".format(satB.rr_km,satB.vv_kmpersec))
    print("satr.rr {} satr.vv {}".format(satr.rr_km,satr.vv_kmpersec))

    dayahead = sat0.jdsatepoch + 1
    sat0 = satfit.delta_t(sat0,dayahead)
    sat1 = satfit.delta_t(sat1,dayahead) # Changing tsince changes the mean elements
    sat2 = satfit.delta_t(sat2,dayahead) # Changing tsince changes the mean elements
    satA = satfit.delta_t(satA,dayahead)
    satB = satfit.delta_t(satB,dayahead)
    print()
    print("At TLE epoch jd {} (Epoch + 1 day)".format(dayahead))
    print("sat1.rr {} sat1.vv {}".format(sat1.rr_km,sat1.vv_kmpersec))
    print("sat2.rr {} sat2.vv {}".format(sat2.rr_km,sat2.vv_kmpersec))
    print("sat0.rr {} sat0.vv {}".format(sat0.rr_km,sat0.vv_kmpersec))
    print("satA.rr {} satA.vv {}".format(satA.rr_km,satA.vv_kmpersec))
    print("satB.rr {} satB.vv {}".format(satB.rr_km,satB.vv_kmpersec))
    print("satr.rr {} satr.vv {}".format(satr.rr_km,satr.vv_kmpersec))

    dayahead = sat0.jdsatepoch + 10
    sat0 = satfit.delta_t(sat0,dayahead)
    sat1 = satfit.delta_t(sat1,dayahead) # Changing tsince changes the mean elements
    sat2 = satfit.delta_t(sat2,dayahead) # Changing tsince changes the mean elements
    satA = satfit.delta_t(satA,dayahead)
    satB = satfit.delta_t(satB,dayahead)
    print()
    print("At TLE epoch jd {} (Epoch + 10 day)".format(dayahead))
    print("sat1.rr {} sat1.vv {}".format(sat1.rr_km,sat1.vv_kmpersec))
    print("sat2.rr {} sat2.vv {}".format(sat2.rr_km,sat2.vv_kmpersec))
    print("sat0.rr {} sat0.vv {}".format(sat0.rr_km,sat0.vv_kmpersec))
    print("satA.rr {} satA.vv {}".format(satA.rr_km,satA.vv_kmpersec))
    print("satB.rr {} satB.vv {}".format(satB.rr_km,satB.vv_kmpersec))
    print("satr.rr {} satr.vv {}".format(satr.rr_km,satr.vv_kmpersec))


if __name__ == '__main__':
    main()