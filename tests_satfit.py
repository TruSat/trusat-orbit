#!/usr/bin/env python

try:
    from unittest2 import TestCase, main, expectedFailure
except:
    from unittest import TestCase, main, expectedFailure

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) 
  
import numpy as np  
from astropy.time import Time
from datetime import datetime

from satfit import jday_to_datetime, Date
from sgp4.ext import jday, invjday

class Tests(TestCase):
    def test_jday(self):
        print("jday...")
        jd = 2454115.05486 # Sunday 14 January 2007 at 13:18:59.9 

        # Reference Astropy as "answer"
        t_astropy = Time(jd, format='jd')

        jday_datetime = jday_to_datetime(jd)
        self.assertRegex(jday_datetime.isoformat(sep=' ',timespec='milliseconds'),t_astropy.iso,msg="jday_to_datetime() failed")

        (year, month, day, hour, minute, second) = invjday(jd)
        jday_jd = jday(year, month, day, hour, minute, second)
        self.assertEqual(jday_jd,jd,"jday() failed")

    def test_Date(self):
        print("Date()...")
        jd = 2454115.05486 # Sunday 14 January 2007 at 13:18:59.9 
        t1 = Date(jd)

        self.assertEqual(t1.jd,jd,"t1.jd failed")
        

        isodatestring = '2007-01-14T13:18:59.904002'





if __name__ == '__main__':
    main()