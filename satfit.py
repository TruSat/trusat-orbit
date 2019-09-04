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

from spacetrack import SpaceTrackClient

# These are necessary until Brandon Rhodes approves pull requests
# https://github.com/brandon-rhodes/python-sgp4/pull/35
sys.path.insert(1, '/Users/chris/Dropbox/code/preMVP/python-sgp4')
# https://github.com/skyfielders/python-skyfield/pull/276
sys.path.insert(2, '/Users/chris/Dropbox/code/preMVP/python-skyfield')

from skyfield.iokit import Loader, download, parse_tle
from skyfield import sgp4lib
from sgp4.ext import jday, invjday, days2mdhms

# The following 5 lines are necessary until our modules are public
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
tle_path = os.path.join(parentdir, "sathunt-tle")
sys.path.insert(1,tle_path) 
from tle_util import make_tle, append_tle_file

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
iod_path = os.path.join(parentdir, "sathunt-iod")
sys.path.insert(1,iod_path) 
import iod

from elfind import read_obs, rref
from satid import unit_vector, mag

# ///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////

# TODO: Make a class of these?
#    tle,      // epoch of elements in tle format
#    ii,       // inclination, degrees
#    om,       // right ascension of ascending node, degrees
#    ec,       // eccentricity
#    ww,       // argument of the perigee, degrees
#    uu,       // longitude
#    ma,       // mean anomaly, degrees
#    nn,       // mean motion, revolutions/day
#    c2,       // internal drag term
#    pgee,     // perigee
#    agee,     // apogee
#    bstar;    // BSTAR drag term

# TODO: Figure out what these flags mean
# int
#    xi = 0, xe = 0, xw = 0, xn = 0;   // flags

# // please note that these variables are in TLE format, i.e. degrees, whereas
# // their counterparts inside SGP4 are in radians.

# double la, lo, hh;          // observer parameters
# int first = 1, batch = 0, out = 0;
# int nobs, len;              // number of observations
# char file[81];              // input file string declared
# char buf[81];               // input buffer
# char name[81], tle1[81], tle2[81];    // original tle data lines
# char line1[81], line2[81];  // read and write tle buffers
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

# //////////// DECLARE FUNCTIONS ////////////////////////////////////////////////

# void get_obs(char *file, char iod_line[][81]);
# void read_obs(char iod_line[][81],
#             double odata[][4],
#             double ll[][3],
#             double rd[][3]);
# void sort(char iod_line[][81],
#             double odata[][4],
#             double ll[][3],
#             double rd[][3]);
# void read_tle(char *file);
# void print_fit(Satellite sat,  double rd[][3],
#               double ll[][3],  double odata[][4]);
# double find_rms(Satellite sat,  double rd[][3],
#                double ll[][3],  double odata[][4]);
# void anomaly(Satellite sat,  double rd[][3],
#             double ll[][3],  double odata[][4]);
# void motion(Satellite sat,  double rd[][3],
#            double ll[][3],  double odata[][4]);
# void perigee(Satellite sat,  double rd[][3],
#             double ll[][3],  double odata[][4]);
# void node(Satellite sat,  double rd[][3],
#          double ll[][3],  double odata[][4]);
# void align(Satellite sat,  double rd[][3],
#          double ll[][3],   double odata[][4]);
# void rref(double m[][7], double b[]);
# void zrll(Satellite sat, double* rd, double& ra, double& dc);
# void diff_el(Satellite sat,  double rd[][3],
#          double ll[][3],   double odata[][4]);
# void longitude(void);
# double ww2ma(double wx);
# void write_tle(char *file);
# void so2r(double r, double rd[], double ll[], double rr[]);
# inline double rtw(double ao, double ac);
# inline char *s_in(const char *prompt, char *buffer);

Class Date(object):
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
                self.time = tle_util.datetime_from_tle_fmt(self.tle)
                self.timevars_from_datetime()
                self.jd = jday(self.yy, self.mm, self.dd, self.hh, self.mm, self.ss)
                self.sidereal()
            else:                   # this date is julian
                self.jd = self.time
                self.calcmjd()
                (self.yy, self.mm, self.dd, self.hh, self.mn, self.ss) = invjday(self.jd)
                self.sidereal()
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
            self.tle = tle_util.tle_fmt_epoch(self.time)

    def self.now():
        """ Re-initializes an existing date object, t1, to the computer
        time at the instant this command is executed by the computer."""
        # TODO - only if needed

    def self.input():
        """Screen input of calendar date variables with current values
        presented as default. Just press ENTER to accept current
        value.  Useful to examine current values at run time."""
        # TODO - final polish

    def self.print():
        """Screen output of calendar date variables."""
        print("Year   {:d}".format(self.yy))
        print("Month  {:d}".format(self.mm))
        print("Day    {:d}".format(self.dd))
        print("Hour   {:d}".format(self.hr))
        print("Minute {:d}".format(self.mn))
        print("Second {:d}".format(self.ss))

    def self.sidereal():
        """calculate Greenwich sidereal hour angle"""
        t = (self.jd - 2451545) / 36525
        thetag = 280.46061837 + 360.98564736629 * (self.jd - 2451545) \
                    + .000387933 * (t * t) - (t * t * t) / 38710000
        self.thetag = fmod(thetag, 360)     # degrees

    def self.timevars_from_datetime():
        self.yy = self.time.year
        self.mm = self.time.month
        self.dd = self.time.day
        self.hh = self.time.hour
        self.mn = self.time.minute
        self.ss = self.time.seconds

    def self.calcmjd():
        self.mjd = self.jd - 2400000.5

def acose(x):
    if ( x >= 1):
        rval = 0
    elif ( x <= -1):
        rval = pi
    else:
        rval = acos(x)
    return rval

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
    nobs = len(odata)

    # FIXME: Accept this CPP code for now, but make this more pythonic.
    # Not needed to DB queries (which can be returned already sorted)
    unsorted = 1
    while (unsorted):
        for i in range(0, nobs - 1):
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
                for s in range (i, nobs - 1):
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
    # copy sat
    satx = sat

    nobs = len(odata)

    for j in range (0, nobs):
        # advance satellite position
        # TODO: replace this python-SGP4 function
        satx.delta_t(odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # TODO: replace this with python-SGP4 code
        delr = satx.rr - rd[j] # This was vmadd -1
        nrr = mag(delr)       # range

        # convert to unit vector
        delr = delr/nrr

        # topocentric position error in degrees
        Perr = degrees( acose( np.dot(delr, ll[j]) ) )

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)

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

def edit_file():
    # TODO: Figure out if it makes sense to implement this (in VIM, or ENV editor?)
    # sprintf(file_in, "%s", file);
    # //  system(file_in);
    # FIXME Fix goto restart
    # goto restart;

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
    buf = input("\nRemove Line Number : ")
    j = int(buf.strip())
    # FIXME Make the deletions pythonic
    for i in range (j, nobs):
        rd[i - 1] = rd[i]
        ll[i - 1] = ll[i]
        odata[i - 1] = odata[i]
        iod_line[i - 1] = iod_line[i])
    nobs-=1

def fit():
   out = 0
   print_fit(sat, rd, ll, odata)
   sat.print_el()
   accept_command()

def id():
    # requires SATID in same folder
    # TODO Implement python function call to elcord here
    # sprintf(buf, "satid %s", file);
    # system(buf);
    # FIXME Fix goto restart
    # goto restart

def accept_command():
    """Accept a new command"""

    while(True):    # Forever loop
        print("\nEnter command")
        print(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee   (A)nomaly  (M)otion  (B)star")
        print(" (S)tep   I(D)   (T)ime  (F)it  (V)iew  (R)emove   (W)rite   (Q)uit")


        cmd = input(": ") 
        cmd = cmd.strip()
        cmd = cmd.upper()

        # Hidden functions
        if (cmd == "G"):
            fit_out()
        if (cmd == "E"):
            edit_file()
        if (cmd == "H"):
            history()
        if (cmd == "Y"):
            edit_history()
        if (cmd == "C"):
            elcor()
        if (cmd == "O"):
            discover()
        if (cmd == "U"):
            maneuver()

        # Visible functions
        if (cmd == "S"):
            step()
        if (cmd == "Z"):
            step()
        if (cmd == "I"):
            sat.print_el()
            incl()
        if (cmd == "N"):
            sat.print_el()
            node()
        if (cmd == "X"):
            sat.print_el()
            xntrcty()
        if (cmd == "P"):
            sat.print_el()
            perigee()
        if (cmd == "A"):
            sat.print_el()
            anomaly()
        if (cmd == "M"):
            sat.print_el()
            motion()
        if (cmd == "B"):
            sat.print_el()
            bstar()

        if (cmd == "D"):
            id()
        if (cmd == "T"):
            time_func()
        if (cmd == "F"):
            fit()
        if (cmd == "V"):
            viewobs()
        if (cmd == "R"):
            remove()
        if (cmd == "W"):
            write_el()
        if (cmd == "Q"):
            sys.exit(0)
    

# /////////////////// MAIN //////////////////////////////////////////////////////

def main():
    """ satfit
    Fit a TLE prediction to a reference TLE + new IOD observations

    Inputs:
        iod_line    Array of IOD formatted observations
        file_in     source file (contains IODs plus TLEs)
        file_out    result file (contains file_in + new TLE predict)
    """
    # char srch[] = {'W'};                     // initialize wide search

    # TODO: Expand this code to deal with command line options like the original
    file = "sat.txt"
    # TODO: Look into batch capability

    # restart: # FIXME: get rid of this go-to block

    #   // find observation lines in input file
    #   // number of observations = nobs (global)
    #   // store observations to iod_line matrix (up to 200 observations)
    IOD_Records = iod.get_iod_records_from_file(file)
    num_file_obs = len(IOD_Records)

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
    TLEs = tle_util.TLEFile(file)

    # Grab just the first TLE instance from the file
    tle1 = TLEs[0].line1
    tle2 = TLEs[0].line2
    print({}.format(TLEs[0].line1))
    print({}.format(TLEs[0].line2))

    # // Create Satellite variable, sat, from TLE
    # TODO: - Replace with satellite Class variable
    Satellite sat(tle, ii, om, ec, ww, ma, nn, bstar);

    # Make a copy of original sat
    save_sat = sat

    # // maneuvering sat
    satm = sat 

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
    nobs = len(odata)
    print("{} Observations Found".format(nobs))

    # print dates
    t2 = Date()
    print("\nTODAY   : {:d}".format(t2.doy))
    if(nobs > 0):
        t3 = Date(time=odata[nobs - 1][0])
        print("\nLAST OB : {:d}".format(t3.doy))

    # Accept a new command
    accept()

if __name__ == '__main__':
    main()




# history:

#    printf("\n(A)dd  (G)raphs  (E)dit  (Q)uit  ");
#    s_in(": ", buf);

#    if ((*buf & 0x5f) == 'Q')
#    {
#      write_tle((char *)"");
#      goto accept;
#    }

#    sprintf(file_in, "history/%05d.txt", ssn);
#    if ((fp = fopen(file_in, "r")) == NULL)
#    {
#      fclose(fp);
#      // create history file
#      fpo = fopen(file_in, "w");
#      //  write current TLE to history
#      fprintf(fpo, "\n");
#      fprintf(fpo, "%s", line1);
#      fprintf(fpo, "%s", line2);
#      fclose(fpo);
#      printf("\nHISTORY created with current TLE\n");
#    }
#    fclose(fp);

#    if ((*buf & 0x5f) == 'G')
#    {
#      i = 1;
#      fp = fopen(file_in, "r");

#      sprintf(file_out, "trendA.csv");
#      fpo = fopen(file_out, "w");
#      fprintf(fpo, "time,inc,omega,ecc,perigee,motion,YYddd,,,%05d\n", ssn);

#      while(fgets(line1, 80, fp))
#      {
#        if(*line1 == '1' && strlen(line1) == 70)
#        {
#          fgets(line2, 80, fp);
#          if(*line2 == '2' && strlen(line2) == 70)
#          {
#            // first data line
#            tleh = atof(line1 + 18);      // epoch in tle format at this point
#            Date t1(tleh);
#            tleh = t1.jd;
#            if(i)
#            {
#              zero = tleh;
#              i = 0;
#            }
#            tleh -= zero;

#            // second data line
#            ik = atof(line2 + 8);         // inclination, degrees
#            ok = atof(line2 + 17);        // ascending node, degrees
#              line2[25] = '.';            // add decimal point to eccentricity
#            ek = atof(line2 + 25);        // eccentricity
#              line2[25] = '.';            // add decimal point to eccentricity
#            wk = atof(line2 + 34);        // perigee, degrees
#            mk = atof(line2 + 43);        // mean anomaly, degrees
#            // Make sure mean motion is null-terminated, since rev. no.
#            // may immediately follow.
#            line2[63] = '\0';
#            nk = atof(line2 + 52);        // mean motion, revolutions/day
#          }
#          fprintf(fpo, "%f,%f,%f,%f,%f,%f,%.5s\n", tleh, ik, ok, ek, wk, nk, line1 + 18);
#          theta = tleh;
#        }
#      }
#      Date t1(tle);
#      tleh = t1.jd;
#      tleh -= zero;
#      if(tleh > theta)
#         fprintf(fpo, "%f,%f,%f,%f,%f,%f,%5.0f\n", tleh, ii, om, ec, ww, nn, tle);
#      fclose(fpo);
#      sprintf(file_out, "trendA.xls");
#      // system(file_out);

#      fclose(fp);
#      // sprintf(file_out, "del trendA.csv");
#      // system(file_out);
#    }
#    if ((*buf & 0x5f) == 'E')
#    {
#      sprintf(file_in, "history%c%05d.txt", '\\', ssn);
#      // system(file_in);
#    }
#    if((*buf & 0x5f) == 'A')
#    {
#      fclose(fp);
#      //  write current TLE to history
#      write_tle(file_in);
#      printf("\nCurrent TLE added to HISTORY\n");
#    }

#    goto history;

# discover:       // partition search

#    printf("\n(W)ide  (N)arrow  [%c", srch[0]);
#    if (strlen(s_in("]: ", buf))) srch[0] = buf[0];
#    if ((*srch & 0x5f) == 'W')
#    {
#       emin = 0;             emax = .2;            esize = 80;
#       wmin = 0;             wmax = 360;           wsize = 144;
#       nmin = 11;            nmax = 15;            nsize = 80;
#    }
#    if ((*srch & 0x5f) == 'N')
#    {
#       emin = ec * .9;       emax = ec * 1.1;      esize = 20;
#       wmin = ww - 2;        wmax = ww + 2;        wsize = 20;
#       nmin = nn - .1;       nmax = nn + .1;       nsize = 20;
#    }
#    if ((*srch & 0x5f) == 'Z')
#    {
#       emin = ec - estep;    emax = ec + estep;    esize = 20;
#       wmin = ww - wstep;    wmax = ww + wstep;    wsize = 20;
#       nmin = nn - nstep;    nmax = nn + nstep;    nsize = 20;
#    }

#    printf("\nemax [%.7f", emax);
#    if (strlen(s_in("]: ", buf)))
#      emax = atof(buf);
#    printf("emin [%.7f", emin);
#    if (strlen(s_in("]: ", buf)))
#      emin = atof(buf);
#    estep = (emax - emin) / esize;

#    printf("\nwmax [%.4f", wmax);
#    if (strlen(s_in("]: ", buf)))
#      wmax = atof(buf);
#    printf("wmin [%.4f", wmin);
#    if (strlen(s_in("]: ", buf)))
#      wmin = atof(buf);
#    wstep = (wmax - wmin) / wsize;

#    printf("\nnmax [%.8f", nmax);
#    if (strlen(s_in("]: ", buf)))
#      nmax = atof(buf);
#    printf("nmin [%.8f", nmin);
#    if (strlen(s_in("]: ", buf)))
#      nmin = atof(buf);
#    nstep = (nmax - nmin) / nsize;

#    for(wk = wmin; wk < wmax; wk += wstep)
#    {
#       theta = (uu - wk) * de2ra;
#       for(ek = emin; ek < emax; ek += estep)
#       {
#          e = acose((ek + cos(theta)) / (1 + ek * cos(theta)));
#          if(theta > pi) e = 2 * pi - e;
#          mk = e - ek * sin(e);
#          mk = mk / de2ra;
#          for(nk = nmin; nk < nmax; nk += nstep)
#          {
#             Satellite satx(sat.jd, ii, om, ek, wk, mk, nk, bstar);
#             // establish the computed ra, dc, at jdo with no perturbations
#             rms = find_rms(satx, rd, ll, odata);
#             if(rms < sum)
#             {
#                sum = rms;       // global
#                ec = ek;
#                ww = wk;
#                ma = mk;
#                nn = nk;
#             }
#          }   // for nk
#       }   // for ek
#    }   // for wk
#    ww = mod(ww, 360);
#    ma = mod(ma, 360);
#    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#    // update mean_anomaly
#    anomaly(sat, rd, ll, odata);
#    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#    // print new elements
#    print_fit(sat, rd, ll, odata);
#    sat.print_el();

#    longitude();

#    srch[0] = 'Z';
#    goto accept;

# step:       // partition search within limits set below

#    // update mean_anomaly
#    anomaly(sat, rd, ll, odata);
#    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#    emax = ec * 1.1;
#    emin = ec * .9;
#    wmax = ww + 2;
#    wmin = ww - 2;
#    imax = ii + 2;
#    imin = ii - 2;
#    omax = om + 2;
#    omin = om - 2;

#    printf("\nPress Q to Quit    :\n\n");
#    while(1)
#    {
#       perigee(sat, rd, ll, odata);
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#       node(sat, rd, ll, odata);
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#       if(nobs > 3 && (*buf & 0x5f) == 'Z')
#       {
#         diff_el(sat, rd, ll, odata);
#         sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#       }

#       // set new limits
#       emax = 1.01 * ec;
#       emin = .99 * ec;
#       wmax = ww + .5;
#       wmin = ww - .5;
#       imax = ii + .5;
#       imin = ii - .5;
#       omax = om + .5;
#       omin = om - .5;

#       printf("rms%12.5f", sum);
#       s_in("    : ", buf);
#       if ((*buf & 0x5f) == 'Q') break;
#    }

#    print_fit(sat, rd, ll, odata);
#    sat.print_el();       // print new elements

#    srch[0] = 'N';
#    goto accept;

# incl:

#    printf("\n(A)uto  (ii)  (Q)uit  ");
#    s_in(": ", buf);

#    if (isdigit(buf[0]))
#    {
#      ii = atof(buf);
#      xi = 1;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      // partition search

#      imax = ii + 2;
#      imin = ii - 2;
#      omax = om + 2;
#      omin = om - 2;
#      xi = 0;

#      printf("\nimax [%.4f", imax);
#      if (strlen(s_in("]: ", buf)))
#        imax = atof(buf);
#      printf("imin [%.4f", imin);
#      if (strlen(s_in("]: ", buf)))
#        imin = atof(buf);

#      printf("\nomax [%.4f", omax);
#      if (strlen(s_in("]: ", buf)))
#        omax = atof(buf);
#      printf("omin [%.4f", omin);
#      if (strlen(s_in("]: ", buf)))
#        omin = atof(buf);

#      node(sat, rd, ll, odata);

#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements

#      srch[0] = 'N';
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto incl;

# node:

#    printf("\n(A)uto  (om)  (Q)uit  ");
#    s_in(": ", buf);

#    if (isdigit(buf[0]))
#    {
#      om = atof(buf);
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      Satellite satx = sat;

#      // partition search
#      omax = om + 2;
#      omin = om - 2;

#      printf("\nomax [%.4f", omax);
#      if (strlen(s_in("]: ", buf)))
#        omax = atof(buf);
#      printf("omin [%.4f", omin);
#      if (strlen(s_in("]: ", buf)))
#        omin = atof(buf);

#      while((omax - omin) > 1e-5)
#      {
#        ostep = fabs(rtw(omax, omin) / 20);
#        for(ok = omin; ok < omax; ok += ostep)
#        {
#           satx.delta_el(sat.jd, ii, ok, ec, ww, ma, nn, bstar);
#           // establish the computed ra, dc, at jdo with no perturbations
#           rms = find_rms(satx, rd, ll, odata);
#           if(rms < sum)
#           {
#              sum = rms;      // global
#              om = ok;
#           }
#        }   // for ok
#        satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#        omin = om - ostep;
#        omax = om + ostep;
#      }
#      om = mod(om, 360);

#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements

#      srch[0] = 'N';
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto node;

# xntrcty:

#    printf("\n(S)earch  (A)uto  (ec)  (Q)uit  ");
#    s_in(": ", buf);

#    emax = ec * 1.05;
#    emin = ec * .95;

#    if (isdigit(buf[0]))
#    {
#      ec = atof(buf);
#      xe = 1;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      xe = 0;
#      printf("\nemax [%.7f", emax);
#      if (strlen(s_in("]: ", buf)))
#        emax = atof(buf);
#      printf("emin [%.7f", emin);
#      if (strlen(s_in("]: ", buf)))
#        emin = atof(buf);

#      while((emax - emin) > 1.e-8)
#      {
#         estep = (emax - emin) / 20;
#         for(ek = emin; ek < emax; ek += estep)
#         {
#            sat.delta_el(sat.jd, ii, om, ek, ww, ma, nn, bstar);
#            // establish the computed ra, dc, at jdo with no perturbations
#            rms = find_rms(sat, rd, ll, odata);
#            if(rms < sum)
#            {
#               sum = rms;
#               ec = ek;
#            }
#         }   // for ek
#         sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#         emin = ec - estep;
#         emax = ec + estep;
#      }
 
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements

#      srch[0] = 'N';
#    }
#    if ((*buf & 0x5f) == 'S')
#    {
#      printf("\nemax [%.7f", emax);
#      if (strlen(s_in("]: ", buf)))
#        emax = atof(buf);
#      printf("emin [%.7f", emin);
#      if (strlen(s_in("]: ", buf)))
#        emin = atof(buf);
#      printf("estep [%.0f", estep);
#      if (strlen(s_in("]: ", buf)))
#        estep = atof(buf);
#      estep = (emax - emin) / estep;

#      ek = ec;
#      printf("\neccentricity  sum");
#      for(ec = emin; ec < emax + estep; ec += estep)
#      {
#        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#        sum = find_rms(sat, rd, ll, odata);
#        printf("\n%.7f     %7.4f", ec, sum);
#      }
#      printf("\n");
#      ec = ek;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto xntrcty;

# perigee:

#    pgee = (1-ec)*sat.aodp;
#    agee = (1+ec)*sat.aodp;

#    printf("\nPerigee = %f er", pgee);
#    printf("     %d X %d km\n", (int)((pgee-1)*6378.135),
#                                (int)((agee-1)*6378.135));

#    printf("\n(S)earch  (A)uto  (ww)  (Q)uit  ");
#    s_in(": ", buf);

#    if (isdigit(buf[0]))
#    {
#      longitude();
#      ww = atof(buf);
#      ma = ww2ma(ww);
#      xw = 1;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      xw = 0;
#      double wx = ww + 1;
#      while(fabs(wx - ww) > 1e-4)
#      {
#        wx = ww;
#        printf("\n%8.4f  %.7f", ww, ec);
#        wmax = ww + 1;
#        wmin = ww - 1;
#        emax = ec * 1.01;
#        emin = ec * .99;
#        perigee(sat, rd, ll, odata);
#        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      }
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements

#      srch[0] = 'N';
#    }
#    if ((*buf & 0x5f) == 'S')
#    {
#      longitude();
#      wmax = ww + 2;
#      wmin = ww - 2;
#      wstep = 20;

#      printf("\nwmax [%.4f", wmax);
#      if (strlen(s_in("]: ", buf)))
#        wmax = atof(buf);
#      printf("wmin [%.4f", wmin);
#      if (strlen(s_in("]: ", buf)))
#        wmin = atof(buf);
#      printf("wstep [%.0f", wstep);
#      if (strlen(s_in("]: ", buf)))
#        wstep = atof(buf);
#      wstep = (wmax - wmin) / wstep;

#      wk = ww;
#      mk = ma;
#      printf("\nperigee      sum");
#      for(ww = wmin; ww < wmax + wstep; ww += wstep)
#      {
#        ma = ww2ma(ww);
#        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#        sum = find_rms(sat, rd, ll, odata);
#        printf("\n%8.4f  %7.4f", mod(ww, 360), sum);
#      }
#      printf("\n");
#      ww = wk;
#      ma = mk;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto perigee;

#    // amax and amin are used as temporary variables
#    longitude();
#    amax = uu * de2ra;
#    amin = sin(sat.xincl/de2ra) * sin(amax);
#    amin *= amin;
#    amin = sqrt(1 - amin);
#    amin = (acose(cos(amax) / amin)) / de2ra;
#    if(fmod(uu, 360) > 180) amin = 360 - amin;
#    if(sat.xincl/de2ra < 90)
#         amax = fmod(sat.xnodeo/de2ra + amin - sat.thetag + 360, 360.);
#    else amax = fmod(sat.xnodeo/de2ra - amin - sat.thetag + 720, 360.);
#    printf("\nE Longitude = %8.4f\n", amax);

# anomaly:

#    if(nn < 1.5)
#    {
#      // amax and amin are used as temporary variables
#      longitude();
#      amax = uu * de2ra;
#      amin = sin(sat.xincl/de2ra) * sin(amax);
#      amin *= amin;
#      amin = sqrt(1 - amin);
#      amin = (acose(cos(amax) / amin)) / de2ra;
#      if(fmod(uu, 360) > 180) amin = 360 - amin;
#      if(sat.xincl/de2ra < 90)
#           amax = fmod(sat.xnodeo/de2ra + amin - sat.thetag + 360, 360.);
#      else amax = fmod(sat.xnodeo/de2ra - amin - sat.thetag + 720, 360.);
#      printf("\nE Longitude = %8.4f\n", amax);
#    }

#    printf("\n(S)earch  (A)uto  (ma)  (L)ast  (Q)uit  ");
#    s_in(": ", buf);

#    amax = 360.;
#    amin = 0.;
#    astep = 20;

#    if (isdigit(buf[0]))
#    {
#      ma = atof(buf);
#      longitude();
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      anomaly(sat, rd, ll, odata);
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#      srch[0] = 'N';
#    }
#    if ((*buf & 0x5f) == 'S')
#    {
#      printf("\namax [%.7f", amax);
#      if (strlen(s_in("]: ", buf)))
#        amax = atof(buf);
#      printf("amin [%.7f", amin);
#      if (strlen(s_in("]: ", buf)))
#        amin = atof(buf);
#      printf("astep [%.0f", astep);
#      if (strlen(s_in("]: ", buf)))
#        astep = atof(buf);
#      astep = (amax - amin) / astep;

#      mk = ma;
#      printf("\nanomaly        sum");
#      for(ma = amin; ma < amax + astep; ma += astep)
#      {
#        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#        sum = find_rms(sat, rd, ll, odata);
#        printf("\n%8.4f     %7.4f", ma, sum);
#      }
#      printf("\n");
#      ma = mk;              // restore
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#    }
#    if ((*buf & 0x5f) == 'L')
#    {
#      align(sat, rd, ll, odata);
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      longitude();
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto anomaly;

# motion:

#    printf("\n(A)uto  (nn)  (Q)uit  ");
#    s_in(": ", buf);

#    if (isdigit(buf[0]))
#    {
#      nn = atof(buf);
#      xn = 1;
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      xn = 0;
#      // update mean motion, no limits
#      motion(sat, rd, ll, odata);
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto motion;

# bstar:

#    printf("\n(A)uto  (b*)  (B)atch  (Q)uit  ");
#    s_in(": ", buf);

#    if (isdigit(buf[0]))
#    {
#      bstar = atof(buf);
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      // update Bstar within limits
#      bmax = bstar * 1.1;
#      bmin = bstar * .9;
#      if(bstar < 0)
#      {
#        bmax = bmin;
#        bmin = bstar * 1.1;
#      }
#      printf("\nbmax [%.8f", bmax);
#      if (strlen(s_in("]: ", buf)))
#        bmax = atof(buf);
#      printf("bmin [%.8f", bmin);
#      if (strlen(s_in("]: ", buf)))
#        bmin = atof(buf);

#      while((bmax - bmin) > 1.e-9)
#      {
#         bstep = (bmax - bmin) / 20;
#         for(bk = bmin; bk < bmax; bk += bstep)
#         {
#            sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bk);
#            // establish the computed ra, dc, at jdo with no perturbations
#            rms = find_rms(sat, rd, ll, odata);
#            if(rms < sum)
#            {
#               sum = rms;
#               bstar = bk;
#            }
#         }   // for bk
#         sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#         bmin = bstar - bstep;
#         bmax = bstar + bstep;
#      }
 
#      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'B')
#    {
#      // create batch file
#      fp = fopen("batch.inp", "w");
#      fprintf(fp, "%s\n\n\n\n\n\n\n", "s");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n\n\n\n", "s");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n\n\n\n", "s");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n\n\n\n", "s");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n\n\n\n", "z");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n\n\n\n", "z");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "b");
#      fprintf(fp, "%s\n\n\n", "a");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "w");
#      fprintf(fp, "%s\n", "u");
#      fprintf(fp, "%s\n", "q");
#      fprintf(fp, "%s\n", "q");
#      fclose(fp);
#      sprintf(file_in, "satfit %s ba<batch.inp", file);
#      // system(file_in);
#      sprintf(file_in, "DEL batch.inp");
#      // system(file_in);
#      sprintf(file_in, "satfit %s", file);
#      // system(file_in);
#      exit(0);
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;

#    goto bstar;

# maneuver:

#    printf("\nBoundary#  (P)erigee  (A)pogee  (O)b  (E)nergy  (R)estore  (Q)uit  ");
#    s_in(": ", buf);
#    if ((*buf & 0x5f) == 'Q')
#    {
#      tle = sat.tle;
#      ii = sat.xincl / de2ra;
#      om = mod(sat.xnodeo / de2ra, 360);
#      ec = sat.eo;
#      ww = mod(sat.omegao / de2ra, 360);
#      ma = mod(sat.xmo / de2ra, 360);
#      nn = sat.xno / nocon;
#      bstar = sat.bstar;
#      goto accept;
#    }
#    if (isdigit(buf[0]))
#    {
#      int p = atoi(buf);

#      // store old obs
#      for(i = 0; i < nobs; i++)
#      {
#        sprintf(iod_linex[i], "%s", iod_line[i]);
#        memcpy(llx[i], ll[i], sizeof(llx[i]));
#        memcpy(rdx[i], rd[i], sizeof(rdx[i]));
#        memcpy(odatax[i], odata[i], sizeof(odatax[i]));
#      }
#      for(i = p; i <= nobs; i++)
#      {
#        sprintf(iod_line[i-p], "%s", iod_line[i-1]);
#        memcpy(ll[i-p], ll[i-1], sizeof(ll[i-p]));
#        memcpy(rd[i-p], rd[i-1], sizeof(rd[i-p]));
#        memcpy(odata[i-p], odata[i-1], sizeof(odata[i-p]));
#      }
#      nobsx = nobs;
#      nobs  = nobs - p + 1;

#      out = 0;
#      print_fit(sat, rd, ll, odata);
#      sat.print_el();
#      printf("\nperiod = %f days\n", nocon/sat.xno);
#    }
#    if ((*buf & 0x5f) == 'P')
#    {
#      double time;      // days

#      printf("\n(P)revious  (N)ext  (Q)uit  ");
#      s_in(": ", buf);
#      if ((*buf & 0x5f) == 'Q') goto maneuver;
#      if ((*buf & 0x5f) == 'P')
#      {
#        // if previously at perigee, back up one revolution
#        if(ma < 0.1) time = satm.jd - nocon/satm.xno*(1 + satm.xmo/twopi);
#        // if not previously at perigee, back up to perigee
#        else         time = satm.jd - satm.xmo*nocon/(satm.xno*twopi);
#      }
#      // advance one revolution
#      if ((*buf & 0x5f) == 'N')
#      {
#        // if previously at perigee, go forward one revolution
#        if(ma > 359.9) time = satm.jd + nocon/satm.xno*(2 - satm.xmo/twopi);
#        // if not previously at perigee, go forward to perigee
#        else           time = satm.jd + nocon/satm.xno*(1 - satm.xmo/twopi);
#      }
#      // move to time and ma at perigee
#      Date t2(time);
#      satm.delta_t(t2.jd);
#      satm.rv2el(satm.rr, satm.vv);
#      satm.jd = t2.jd;
#      ma = mod(satm.xmo / de2ra, 360);
#      // refine perigee
#      for(i = 0; i < 3; i++)
#      {
#        // go forward
#        if(ma > 359.9) time = nocon/satm.xno*(satm.xmo/twopi - 1);
#        // back up
#        else time = satm.xmo*nocon/(satm.xno*twopi);
#        Date t1(satm.jd - time);
#        satm.delta_t(t1.jd);
#        satm.jd = t1.jd;
#        satm.rv2el(satm.rr, satm.vv);
#        ma = mod(satm.xmo / de2ra, 360);
#      }
#      printf("\nPERIGEE");
#      satm.print_el();       // print new elements

#      // perigee residual
#      double delr[3];
#      time = sat.jd;                     // save sat epoch
#      sat.delta_t(satm.jd);              // move sat to perigee
#      vmadd(sat.rr, satm.rr, delr, -1);  // compare sat and satm perigees
#      sat.delta_t(time);                 // restore sat epoch
#      printf("\nperigee delta %5.0f\n", norm(delr)*6378.135);
#    }
#    if ((*buf & 0x5f) == 'A')
#    {
#      double time;      // days

#      // time to travel from perigee to apogee, days
#      printf("\n(P)revious  (N)ext  (Q)uit  ");
#      s_in(": ", buf);
#      if ((*buf & 0x5f) == 'Q') goto maneuver;
#      if ((*buf & 0x5f) == 'P')
#      {
#        // if previously at or past apogee and past next perigee, back up to apogee
#        if(ma < 180.1) time = satm.jd - 0.5*nocon/satm.xno*(1 + satm.xmo/pi);
#        // if previously past apogee and before next perigee, back up to apogee
#        else           time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi);
#      }
#      // advance to apogee
#      if ((*buf & 0x5f) == 'N')
#      {
#        // if previously at or past apogee and before perigee, go forward to apogee
#        if(ma > 179.9) time = satm.jd + 0.5*nocon/satm.xno*(3 - satm.xmo/pi);
#        // if previously past apogee and past next perigee, go forward to apogee
#        else           time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi);
#      }
#      // move time and satm.xmo to apogee
#      Date t2(time);
#      satm.delta_t(t2.jd);
#      satm.jd = t2.jd;
#      satm.rv2el(satm.rr, satm.vv);
#      for(i = 0; i < 3; i++)
#      {
#        // loop to refine apogee, find when satm.xmo = pi
#        Date t1(satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi));
#        satm.delta_t(t1.jd);
#        satm.jd = t1.jd;
#        satm.rv2el(satm.rr, satm.vv);
#      }
#      ma = mod(satm.xmo / de2ra, 360);
#      printf("\nAPOGEE");
#      satm.print_el();       // print new elements

#      // apogee residual
#      double delr[3];
#      time = sat.jd;                     // save sat epoch
#      sat.delta_t(satm.jd);              // move sat to apogee
#      vmadd(sat.rr, satm.rr, delr, -1);  // compare sat and satm apogees
#      sat.delta_t(time);                 // restore sat epoch
#      printf("\napogee delta %5.0f\n", norm(delr)*6378.135);  // kilometers
#    }
#    if ((*buf & 0x5f) == 'O')
#    {
#      // move up one to make room for pseudo ob
#      for(i = nobs-1; i >= 0; i--)
#      {
#        sprintf(iod_line[i+1], "%s", iod_line[i]);
#        memcpy(ll[i+1], ll[i], sizeof(ll[i+1]));
#        memcpy(rd[i+1], rd[i], sizeof(rd[i+1]));
#        memcpy(odata[i+1], odata[i], sizeof(odata[i+1]));
#      }
#      nobs++;
#      // ll is unit vector in same direction as satm.rr
#      smult(1./norm(satm.rr), satm.rr, ll[0]);
#      // rd is unit vector (1 er) in same direction as satm.rr
#      smult(1./norm(satm.rr), satm.rr, rd[0]);
#      // odata holds epoch, ra, dc, obscode data
#      double ra, dc;
#      odata[0][0] = satm.jd;
#      zrll(satm, rd[0], ra, dc);
#      odata[0][1] = ra*de2ra;
#      odata[0][2] = dc*de2ra;
#      odata[0][3] = 0000;

#      // print pseuso ob
#      char desig[5];
#      int norad, yy;
#      double rm, dm;

# /*
#           1         2         3         4         5         6
# 0123456789012345678901234567890123456789012345678901234567890
# 36868 10 039A   2018 G 20100909203108300 17 25 2216290+034350 57
# */

#     norad = atoi(iod_line[1]);
#     yy    = atoi(iod_line[1] + 6);
#     sscanf(iod_line[1] + 9, "%4s", &desig);
#     desig[4] = '\0';
#     Date t1(satm.jd);
#     ra /= 15;
#     rm = modf(ra, &rm)*60000;
#     dm = fabs(modf(dc, &dm)*6000);
#     printf("\nPSEUDO OB:");
#     printf("\n%05d %02d %4s   0000 P %4.0f%02.0f%02.0f%02.0f%02.0f%05.0f 16 25 %02.0f%05.0f%+03.0f%04.0f\n",
#             norad, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr,
#             t1.mn, t1.ss*1000, ra, rm, dc, dm);
#     sprintf(iod_line[0], "%05d %02d %4s   0000 P %4.0f%02.0f%02.0f%02.0f%02.0f%05.0f 16 25 %02.0f%05.0f%+03.0f%04.0f\n",
#             norad, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr,
#             t1.mn, t1.ss*1000, ra, rm, dc, dm);


#    }
#    if ((*buf & 0x5f) == 'E')
#    {
#      s_in("\ndelta specific energy(m^2/s^2): ", buf);
#      if (isdigit(buf[0]))
#      {
#        double dE = atof(buf);
#        dE /= 11300.168;   // er^2 / min^2
#        double E2, VV, dV[3], vec[3];
#        cross(satm.rr, satm.vv, vec);
#        double mu = norm(vec);
#        mu = mu*mu;
#        mu = mu / satm.aodp;
#        mu = mu / (1 - satm.eo*satm.eo);
#        E2 = -0.5*mu / satm.aodp;
#        VV = sqrt(2*(E2 + dE + mu/norm(satm.rr)));  // new velocity magnitude
#        unitv(satm.vv, dV);         // velocity direction
#        smult(VV, dV, vec);         // new velocity vector
#        sat = satm;
#        sat.rv2el(sat.rr, vec);     // State vectors to mean elements
#        sat.print_el();             // print new elements
#      }
#    }
#    if ((*buf & 0x5f) == 'R')
#    {
#      printf("\nRestore (O)bs  (E)lements  (Q)uit  ");
#      s_in(": ", buf);
#      // re-store old obs
#      if ((*buf & 0x5f) == 'O')
#      {
#        nobs = nobsx;
#        for(i = 0; i < nobs; i++)
#        {
#          sprintf(iod_line[i], "%s", iod_linex[i]);
#          memcpy(ll[i], llx[i], sizeof(ll[i]));
#          memcpy(rd[i], rdx[i], sizeof(rd[i]));
#          memcpy(odata[i], odatax[i], sizeof(odata[i]));
#        }
#      }
#      // re-store old TLE
#      if ((*buf & 0x5f) == 'E')
#      {
#        // replace working elements with original
#        satm = save_sat;        // original elements maneuverd to a node
#        sat  = save_sat;        // new fit elements
#        sat.print_el();         // print original elements
#      }
#    }

#    goto maneuver;

# write_el:
#    tle = sat.tle;
#    ii = sat.xincl / de2ra;
#    om = sat.xnodeo / de2ra;
#    ec = sat.eo;
#    ww = sat.omegao / de2ra;
#    ma = sat.xmo / de2ra;
#    nn = sat.xno / nocon;
#    c2 = sat.c2;
#    bstar = sat.bstar;

#    printf("\n(U)pdate (V)iew (A)ppend (O)riginal (R)estore (Q)uit");
#    s_in(": ", buf);

#    if ((*buf & 0x5f) == 'A') write_tle(file);
#    if ((*buf & 0x5f) == 'V') sat.print_el();
#    if ((*buf & 0x5f) == 'O')
#    {
#       // view original
#       save_sat.print_el();       // print new elements
#    }
#    if ((*buf & 0x5f) == 'R')
#    {
#       // replace TLE in file with original
#       FILE *fp;
#       fp = fopen(file, "w");
#       fprintf(fp, "%s", name);
#       fprintf(fp, "%s", tle1);
#       fprintf(fp, "%s\n", tle2);
#       for(int i = 0; i < nobs; i++)
#       {
#          fprintf(fp, "%s", iod_line[i]);
#       }
#       fclose(fp);
#       // replace working elements with original
#       sat = save_sat;
#       sat.print_el();            // print original elements

#    }
#    if ((*buf & 0x5f) == 'U')
#    {

#       write_tle((char *)"");     // updates lines and prints to screen
#       //  print_file
#       FILE *fp;

#       len = strlen(name);
#       if (name[len-1] == '\n')
#           name[len-1] = '\0';
#       pgee = (1-ec)*sat.aodp;
#       agee = (1+ec)*sat.aodp;
#       sprintf(buf, "%d X %d km", (int)((pgee-1)*6378.135),
#                                  (int)((agee-1)*6378.135));
#       fp = fopen(file, "w");
#       fprintf(fp, "%-50s%19s\n", name, buf);
#       fprintf(fp, "%s\n", line1);
#       fprintf(fp, "%s\n\n", line2);
#       for(int i = 0; i < nobs; i++)
#       {
#          fprintf(fp, "%s", iod_line[i]);
#       }
#       fclose(fp);
#    }
#    if ((*buf & 0x5f) == 'Q') goto accept;
#    goto write_el;

# time_func:
#    double z1, z2, z3, xns, sinio, beta, time1, time2, time3;
#    ec = sat.eo;
#    ww = sat.omegao / de2ra;
#    nn = sat.xno / nocon;
#    xns = 2160 * sat.bstar * sat.c2 * sat.xno / nocon;
#    if(nobs > 0) time2 = odata[nobs - 1][0];
#    else
#    {
#       t2.now();
#       time2 = t2.jd;
#    }

#    sat.delta_t(time2);
#    z2 = sat.rr[2];
#    do
#    {
#       time1 = time2;
#       z1 = z2;
#       time2 -= .01;
#       sat.delta_t(time2);
#       z2 = sat.rr[2];
#    }while(z1 < 0 || z2 > 0);

#    while(time1 - time2 > 1e-9)
#    {
#       time3 = (time1 + time2) / 2;
#       sat.delta_t(time3);
#       z3 = sat.rr[2];
#       if(z3 < 0) time2 = time3;
#       else time1 = time3;
#    }
#    Date t1(time2);
#    t1.input();
#    sat.delta_t(t1.jd);             // advance to node
#    sat.rv2el(sat.rr, sat.vv);      // sgp4 elements
#    tle = t1.tle;
#    sat.jd = t1.jd;

#    ii = sat.xincl / de2ra;
#    om = mod(sat.xnodeo / de2ra, 360);
#    ec = sat.eo;
#    ww = mod(sat.omegao / de2ra, 360);
#    ma = mod(sat.xmo / de2ra, 360);
#    nn = sat.xno / nocon;
#    bstar = sat.bstar;
#    sat.print_el();       // print new elements
#    sum = find_rms(sat, rd, ll, odata);
#    printf("\nrms%12.5f\n", sum);

#    longitude();

#    goto accept;

# }    // end main

# ////////////////// FUNCTIONS //////////////////////////////////////////////////


# //  asymptotic extrapolation
# inline double asym(double a1, double a2, double a3)
# {
#   double b;
#   if(fabs(a1 - a2) < 1.e-5) b = 0.;
#   else b = (a2 - a3) / (a1 - a2);
#   return (b * a2 - a3) / (b - 1.);
# }

# /* Converts quasi scientific notation to double.  Format:  SnSn where S is
# either '-' or '\0' and n represents a sequence of 1 or more digits.  An
# implied decimal point exists to the left of the left n.  Total length 27
# chars max. Subroutine of read_tle*/
# double sci(char *string)
# {
#    char buf[30], *bufptr;

#    bufptr = buf;

#    if (string[1] == '\0')
#       return 0.;

#    /* get significand */
#    if (*string == '-')
#       *bufptr++ = '-';
#    *bufptr = '.';
#    while (isdigit(*++bufptr = *++string))
#       ;   /* copy significand */
#    *bufptr++ = 'E';

#    /* get exponent */
#    if (*string == '\0')   /* no sign on exponent */
#       ++string;

#    strcpy(bufptr, string);      /* copy exponent */
#    return atof(buf);
# }





# // read site data, subroutine of read_obs
# void read_site(char line[])
# {
#   int site, num;
#   char* fl;
#   char inp_str[81];
#   sscanf(line + 16, "%4d", &site);    // read site number from iod_line
#   // open stations.in file to read site data for each observation
#   FILE *fp;
#   if((fp = fopen("data/stations.in", "r")) != NULL)
#   {
#     while(fgets(inp_str, 80, fp))
#     {
#       // global la, lo, hh
#       sscanf(inp_str, "%d %s %lf %lf %lf",&num, &fl, &la, &lo, &hh);
#       if(num == site)
#       {
#         fclose(fp);
#         return;
#       }
#     } // end while
#     printf("\nno site data for %d\n", site);
#     s_in("[exit]", buf);
#     fclose(fp);
#     exit(0);
#   } // end if
#   else
#   {
#     printf("\nno stations.in file found\n");
#     s_in("[exit]", buf);
#     exit(0);;
#   }
# }





# // write TLE to output file and to screen
# void write_tle(char *file_out)
# {
#   // char ec_string[9];
#   char xns_string[13];
#   char bstar_string[13];
#   char bstar_fract[7];
#   char bstar_exp[3];

#   // sprintf(ec_string, "%.7f", ec);
#   // ec_string[0] = ec_string[2];
#   // ec_string[1] = ec_string[3];
#   // ec_string[2] = ec_string[4];
#   // ec_string[3] = ec_string[5];
#   // ec_string[4] = ec_string[6];
#   // ec_string[5] = ec_string[7];
#   // ec_string[6] = ec_string[8];

#   sprintf(bstar_string, "%12.4e", bstar*10);
#   bstar_fract[0] = bstar_string[0]; // sign
#   bstar_fract[1] = bstar_string[1];
#   bstar_fract[2] = bstar_string[3];
#   bstar_fract[3] = bstar_string[4];
#   bstar_fract[4] = bstar_string[5];
#   bstar_fract[5] = bstar_string[6];
#   bstar_fract[6] = '\0';
#   bstar_exp[0] = bstar_string[8];
#   bstar_exp[1] = bstar_string[11];
#   bstar_exp[2] = '\0';

#   double xns = 2160 * bstar * nn * c2;

#   sprintf(xns_string, "%.8f", xns);
#   if(xns_string[0] == '-')
#   {
#     xns_string[1] = xns_string[2];
#     xns_string[2] = xns_string[3];
#     xns_string[3] = xns_string[4];
#     xns_string[4] = xns_string[5];
#     xns_string[5] = xns_string[6];
#     xns_string[6] = xns_string[7];
#     xns_string[7] = xns_string[8];
#     xns_string[8] = xns_string[9];
#     xns_string[9] = xns_string[10];
#     xns_string[10] = xns_string[11];
#   }

#   sprintf(line1, "1 %05dU %-8s %014.8f %10s  00000-0 %6s%2s 0    00"
#            ,ssn, desig, tle, xns_string, bstar_fract, bstar_exp);
#   ccksum(line1);
#   // sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
#   //          ssn, ii, om, ec_string, ww, ma, nn);
#   sprintf(line2, "2 %05d %8.4lf %8.4lf %07ld %8.4lf %8.4lf %11.8lf    00",
#            ssn,
#            ii,
#            om,
#            (long) (ec * 10000000. + .5),
#            ww,
#            ma,
#            nn);
#   ccksum(line2);

#   if(strlen(file_out) > 0)
#   {
#     FILE *fp;
#     fp = fopen(file_out, "a");
#     fprintf(fp, "\n%s\n", line1);
#     fprintf(fp, "%s\n", line2);
#     fclose(fp);
#   }

#   printf("\n%s\n", line1);
#   printf("%s\n", line2);
# }

# inline char *s_in(const char *prompt, char *buffer)
# {
#   char *p;
 
#   printf("%s", prompt);
#   if (feof(stdin)) exit(0);
 
#   if (fgets(buffer,sizeof(buffer),stdin)) {
#     p = strchr(buffer, '\n');
#     if (p) {
#       *p = '\0';
#     }
#     return buffer;
#   } else {
#     return NULL;
#   }
# }

# // round the world
# inline double rtw(double ao, double ac)
# {
#   if(fabs(ao - ac) > 180)
#   {
#     if(ao < 180) ao += 360;
#     else ac += 360;
#   }
#   return ao - ac;
# }


# // find mean anomaly from true longitude and perigee
# double ww2ma(double wx)
# {
#    theta = mod(uu - wx, 360) * de2ra;
#    e = acose((ec + cos(theta)) / (1 + ec * cos(theta)));
#    if(theta > pi) e = 2 * pi - e;
#    ma = e - ec * sin(e);
#    ma /= de2ra;
#    return ma;
# }


# /* Calculates predicted direction angles, right ascension(ra) and
# declination(dc), of the line of sight vector, rll, in equatorial coordinates
# from an observer position, rd, to the satellite position, sat.rr.
# A subroutine of diff_el.  Satellite object, sat, simplifies prediction */
# void zrll(Satellite satx,  // perturbed Satellite object at TLE
#           double* rd,      // input, topocentric vector to observer
#           double& ra,      // calculated right ascension, degrees
#           double& dc)      // calculated declination, degrees
# {
#   double rll[3];                 // declare line of sight vector
#   vmadd(satx.rr, rd, rll, -1);   // rll = satx.rr - rd
#   // return ra and dc, degrees
#   dc = asin(rll[2] / norm(rll)) / de2ra;
#   ra = atan(rll[1] / rll[0]) / de2ra;
#   if(rll[0] < 0) ra += 180.;
#   ra = mod(ra, 360);
# }

# //////////////////////////////////////////////////////

# // differential correction of elements
# void diff_el(Satellite sat,  double rd[][3],
#             double ll[][3],  double odata[][4])
# {

#    int i, j, k, c = 0;
#    double ra, dc;                 // right ascension, declination variables
#    double delta, el;              // small change amount
#    double mdata[150][7];          // define matrix for multiple regression
#    double dat[2][7];              // define output matrix of zpde utility
#    double b[6];                   // output deltas, sat position, velocity, b*
#    double rv[6][7];               // normal matrix
#    double rac, dcc, rms;

#    loop:

#    // begin the main loop of the program
#    for(i = 0; i < nobs; i++)
#    {
#      Satellite satx = sat;
#      Satellite satz = sat;             // differential sat
#      satx.delta_t(odata[i][0]);        // relocate satellite at new time

#      // first establish the computed ra, dc, at jdo with no perturbations
#      zrll(satx, rd[i], ra, dc);        // output ra, dc, degrees
#      rac = ra;                         // store computed ra and dc
#      dcc = dc;

#      // find the deltas and load into output matrix, dat
#      dat[0][6] = rtw(odata[i][1]/de2ra, rac);      // store delta_ra
#      dat[1][6] = odata[i][2]/de2ra - dcc;          // store delta_dc

#      // 6 steps to build differential correction matrix
#      j = 0;
#      if(xi)
#      {
#         dat[0][j] = .001;
#         dat[1][j] = .001;
#      }
#      else
#      {
#         delta = .001;                         // change
#         el = satz.xincl;                      // store reference
#         satz.xincl += delta;                  // delta element
#         satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#         zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#         satz.xincl = el;                      // restore reference
#         dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#         dat[1][j] = (dc - dcc) / delta;
#      }

#      j = 1;
#      delta = .001;                         // change
#      el = satz.xnodeo;                     // store reference
#      satz.xnodeo += delta;                 // delta element
#      satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#      zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#      satz.xnodeo = el;                     // restore reference
#      dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#      dat[1][j] = (dc - dcc) / delta;

#      j = 2;
#      if(xe)
#      {
#         dat[0][j] = .00001;
#         dat[1][j] = .00001;
#      }
#      else
#      {
#         delta = .0001;                        // change
#         el = satz.eo;                         // store reference
#         satz.eo += delta;                     // delta element
#         satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#         zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#         satz.eo = el;                         // restore reference
#         dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#         dat[1][j] = (dc - dcc) / delta;
#      }

#      j = 3;
#      if(xw)
#      {
#         dat[0][j] = .001;
#         dat[1][j] = .001;
#      }
#      else
#      {
#         delta = .001;                         // change
#         el = satz.omegao;                     // store reference
#         satz.omegao += delta;                 // delta element
#         satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#         zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#         satz.omegao = el;                     // restore reference
#         dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#         dat[1][j] = (dc - dcc) / delta;
#      }

#      j = 4;
#      delta = .001;                         // change
#      el = satz.xmo;                        // store reference
#      satz.xmo += delta;                    // delta element
#      satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#      zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#      satz.xmo = el;                        // restore reference
#      dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#      dat[1][j] = (dc - dcc) / delta;

#      j = 5;
#      if(xn)
#      {
#         dat[0][j] = .000001;
#         dat[1][j] = .000001;
#      }
#      else
#      {
#         delta = .00001;                       // change
#         el = satz.xno;                        // store reference
#         satz.xno += delta;                    // delta element
#         satz.delta_t(odata[i][0]);            // recalculate with perturbed element
#         zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
#         satz.xno = el;                        // restore reference
#         dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
#         dat[1][j] = (dc - dcc) / delta;
#      }

#      memcpy(mdata[2 * i],dat[0],sizeof(mdata[2 * i]));        // numerical deltas transferred to
#      memcpy(mdata[2 * i + 1], dat[1], sizeof(mdata[2 * i + 1]));        // multiple regresssion matrix
#    }   // end for i

#    // multiple regression
#    for(j = 0; j < 6; j++)
#    {
#      for(k = 0; k < 7; k++)
#      {
#        rv[j][k] = 0;
#        for(i = 0; i < nobs*2; i++)
#        {
#          rv[j][k] = rv[j][k] + mdata[i][j] * mdata[i][k];
#        }
#      }
#    }

#    rref(rv, b);

#    Satellite saty = sat;
#    // test update components with deltas
#    saty.xincl  += b[0]*.1;
#    saty.xnodeo += b[1]*.1;
#    saty.eo     += b[2]*.1;
#    saty.omegao += b[3]*.1;
#    saty.xmo    += b[4]*.1;
#    saty.xno    += b[5]*.1;
#    rms = find_rms(saty, rd, ll, odata);
#    if(rms < sum)
#    {
#      sum = rms;
#      // update components with deltas
#      sat.xincl  += b[0]*.1;
#      sat.xnodeo += b[1]*.1;
#      sat.eo     += b[2]*.1;
#      sat.omegao += b[3]*.1;
#      sat.xmo    += b[4]*.1;
#      sat.xno    += b[5]*.1;
#      c++;
#      if(c < 20) goto loop;
#    }
# /*
#    // display computed deltas for all 6 components and an rms
#    for(i = 0; i < 6; i++)
#    {
#      printf("\n");
#      printf("   %9.6f", b[i]);
#    }
#    printf("\n\nrms%9.6f\n", rms);
#    s_in(": ", buf);
# */
#    // global elements updated
#    ii = sat.xincl/de2ra;
#    om = sat.xnodeo/de2ra;
#    ec = sat.eo;
#    ww = sat.omegao/de2ra;
#    ma = sat.xmo/de2ra;
#    nn = sat.xno/nocon;
# }

# //////////////////////////////////////////////////////

# // box search, no limits
# void anomaly(Satellite sat,  double rd[][3],
#             double ll[][3],  double odata[][4])
# {
#    double min, max, step;
#    step = .1;
#    mk = ma;
#    sum = find_rms(sat, rd, ll, odata);
#    do
#    {
#       min = mk;
#       max = mk;
#       nsum:
#       min = mk - step;
#       sat.delta_el(sat.jd, ii, om, ec, ww, min, nn, bstar);
#       nsum = find_rms(sat, rd, ll, odata);
#       if(nsum < sum)
#       {
#          mk = min;
#          sum = nsum;
#          goto nsum;
#       }
#       xsum:
#       max = mk + step;
#       sat.delta_el(sat.jd, ii, om, ec, ww, max, nn, bstar);
#       xsum = find_rms(sat, rd, ll, odata);
#       if(xsum < sum)
#       {
#          mk = max;
#          sum = xsum;
#          goto xsum;
#       }
#       step /= 2;
#    }while(fabs(max - min) > 1e-5);
#    ma = mod(mk, 360);
# }

# // box search, no limits
# void motion(Satellite sat,  double rd[][3],
#            double ll[][3],  double odata[][4])
# {
#    double min, max, step;
#    step = .1;
#    nk = nn;
#    sum = find_rms(sat, rd, ll, odata);
#    do
#    {
#       min = nk;
#       max = nk;
#       nsum:
#       min = nk - step;
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, min, bstar);
#       nsum = find_rms(sat, rd, ll, odata);
#       if(nsum < sum)
#       {
#          nk = min;
#          sum = nsum;
#          goto nsum;
#       }
#       xsum:
#       max = nk + step;
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, max, bstar);
#       xsum = find_rms(sat, rd, ll, odata);
#       if(xsum < sum)
#       {
#          nk = max;
#          sum = xsum;
#          goto xsum;
#       }
#       step /= 2;
#    }while(fabs(max - min) > 1e-10);
#    nn = nk;
# }

# // partition search on node and inclination within set limits
# void node(Satellite sat,  double rd[][3],
#          double ll[][3],  double odata[][4])
# {
#    double rms;
#    while((imax - imin) > 1e-5)
#    {
#       istep = (imax - imin) / 20;
#       ostep = fabs(rtw(omax, omin) / 20);
#       for(ik = imin; ik < imax; ik += istep)
#       {
#          if(xi)
#          {
#             imin = ii;
#             imax = ii;
#             ik = ii;
#             istep = 0;
#          }
#          for(ok = omin; ok < omax; ok += ostep)
#          {
#             Satellite satx(tle, ik, ok, ec, ww, ma, nn, bstar);
#             // establish the computed ra, dc, at jdo with no perturbations
#             rms = find_rms(satx, rd, ll, odata);
#             if(rms < sum)
#             {
#                sum = rms;
#                ii = ik;
#                om = ok;
#             }
#          }   // for ok
#       }   // for ik
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#       imin = ii - istep;
#       imax = ii + istep;
#       omin = om - ostep;
#       omax = om + ostep;
#    }
#    om = mod(om, 360);
# }

# // partition search on perigee and eccentricity
# void perigee(Satellite sat,  double rd[][3],
#             double ll[][3],  double odata[][4])
# {
#    double rms;
#    if(ec > .1)
#    {
#       wmax = ww + .1;
#       wmin = ww - .1;
#       emax = ec * 1.01;
#       emin = ec * .99;
#    }
#    while((wmax - wmin) > 1e-5)
#    {
#       estep = (emax - emin) / 20;
#       wstep = (wmax - wmin) / 20;
#       for(wk = wmin; wk < wmax; wk += wstep)
#       {
#          if(xw)
#          {
#             wmin = ww;
#             wmax = ww;
#             wk = ww;
#             wstep = 0;
#          }
#          theta = (uu - wk) * de2ra;
#          for(ek = emin; ek < emax; ek += estep)
#          {
#             if(xe)
#             {
#                emin = ec;
#                emax = ec;
#                ek = ec;
#                estep = 0;
#             }
#             e = acose((ek + cos(theta)) / (1 + ek * cos(theta)));
#             if(theta > pi) e = 2 * pi - e;
#             mk = e - ek * sin(e);
#             mk = mk / de2ra;

#             Satellite satx(sat.jd, ii, om, ek, wk, mk, nn, bstar);
#             // establish the computed ra, dc, at jdo with no perturbations
#             rms = find_rms(satx, rd, ll, odata);
#             if(rms < sum)
#             {
#                sum = rms;
#                ec = ek;
#                ww = wk;
#                ma = mk;
#             }
#          }   // for ek
#       }   // for wk
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#       emin = ec - estep;
#       emax = ec + estep;
#       wmin = ww - wstep;
#       wmax = ww + wstep;
#    }
                
#    // update mean_anomaly
#    anomaly(sat, rd, ll, odata);
#    sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

#    // update mean_motion
#    if(!xn)
#    {
#       motion(sat, rd, ll, odata);
#       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
#    }

#    // calculate uu, degrees
#    longitude();

#    ww = mod(ww, 360);
#    ma = mod(ma, 360);
# }

# /* find geocentric vector, rr, to satellite given the topocentric unit
# vector, ll, from the observer to the satellite, the geocentric vector, rd,
# to the observer, and the length, r, of the geocentric vector sat.rr. */
# void so2r(double r, double rd[], double ll[], double rr[])
# {
#   double ang1, ang2, nrd, rho;

#   nrd = norm(rd);
#   ang1 = acose(dot(rd, ll) / nrd);
#   if(ang1 < .001) rho = r - nrd;
#   else
#   {
#     ang2 = asin(nrd * sin(ang1) / r);
#     rho = r * sin(ang1 - ang2) / sin(ang1);
#   }
#   vmadd(rd, ll, rr, rho);
# }

# // sets deltaT error at last obs epoch to zero
# void align(Satellite sat,  double rd[][3],
#          double ll[][3],  double odata[][4])
# {
#    double mk, nrr, nvv, Perr, delt, rr[3], delr[3];
#    int last = nobs - 1;

#    Satellite satx = sat;        // copy sat

#    do
#    {
#      // advance satellite position
#      satx.delta_t(odata[last][0]);
#      nrr = norm(satx.rr);
#      nvv = norm(satx.vv);

#      // observed geocentric position vector, rr
#      so2r(nrr, rd[last], ll[last], rr);

#      // position error in radians
#      Perr = acose(dot(satx.rr, rr) / (nrr*nrr));

#      // difference between computed and observed position vectors, er
#      vmadd(satx.rr, rr, delr, -1);

#      // magnitude of delta r in direction of v, radians
#      delt = Perr * dot(satx.vv, delr) / (nvv * norm(delr));

#      // degrees
#      delt /= de2ra;

#      if(first)
#      {
#        first = 0;
#        ma = ma - 0.75*delt;
#      }
#      else ma = ma - delt/5.;
#      satx.delta_el(satx.jd, ii, om, ec, ww, ma, nn, bstar);
#    }while(fabs(delt) > 1.E-5);
# }

# void print_fit(Satellite sat,  double rd[][3],
#               double ll[][3],  double odata[][4])
# {
#   double nrr, nvv, Perr, delt, xtrk, az, el, asp, alpha, sum = 0,
#          tempv[3], temp[3], rr[3], nv[3], delr[3], zz[3] = {0, 0, 1};
#   int yy, day, hh, mm, ss, sign;
#   if(nobs == 0)
#   {
#      printf("\nno obs");
#      return;
#   }

#   // copy sat
#   Satellite satx = sat;

#   if(out)
#   {
#       FILE *fp;
#       fp = fopen(file, "a");
#       fprintf(fp, "\n\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr\n");
#       fclose(fp);
#   }

#   printf("\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr\n");
#   for(int j = 0; j < nobs; j++)
#   {
#      // advance satellite position
#      satx.delta_t(odata[j][0]);
#      nrr = norm(satx.rr);
#      nvv = norm(satx.vv);

#      //  computing elevation
#      el = acose(dot(rd[j], ll[j]) / norm(rd[j]));
#      el /= de2ra;
#      el = 90 - el;
#      el = round(el, 1);

#      // computing aspect
#      asp = acose(dot(ll[j], satx.vv) / nvv);
#      asp /= de2ra;
#      asp = 180 - asp;
#      asp = round(asp, 1);

#      // computing azimuth
#      cross(rd[j], zz, tempv);
#      cross(tempv, rd[j], nv);
#      cross(rd[j], ll[j], tempv);
#      cross(tempv, rd[j], temp);
#      az = acose(dot(nv, temp) / (norm(nv)*norm(temp)));
#      cross(temp, nv, tempv);
#      if(dot(tempv, rd[j]) < 0) az = 2*pi - az;
#      if(norm(temp) == 0) az = -0.0;
#      az /= de2ra;
#      az  = round(az, 1);

#      // observed satellite geocentric position vector, rr
#      so2r(nrr, rd[j], ll[j], rr);

#      // geocentric position error angle in radians
#      Perr = acose(dot(satx.rr, rr) / (nrr*nrr));

#      // difference between computed and observed position vectors, delr, in e.r.
#      vmadd(satx.rr, rr, delr, -1);
#      cross(satx.rr, satx.vv, temp);  // xtrk reference vector points left of track
#      sign = dot(delr, temp) < 0 ? -1 : 1;

#      // observer velocity vector
#      cross(zz, rd[j], tempv);
#      smult(.004351409367, tempv, temp);
#      // observed satellite velocity
#      vmadd(satx.vv, temp, tempv, -1);
#      nvv = norm(tempv);

#      // angle between delr vector and tempv vector, radians
#      alpha = acose(dot(tempv, delr) / (nvv * norm(delr)));

#      // magnitude of delr in direction of tempv, radians
#      delt = atan(cos(alpha) * tan(Perr));   // geocentric range error

#      // time error
#      delt *= nrr / nvv;                     // delta r in min
#      delt *= 60;                            // seconds
#      delt  = round(delt, 2);

#      // predicted topocentric coordinates (sat xyz - observer xyz)
#      // new use of delr variable, predicted line of sight vector
#      vmadd(satx.rr, rd[j], delr, -1);
#      nrr = norm(delr);

#      // convert to unit vector
#      smult(1/nrr, delr, delr);

#      // topocentric position error angle in radians
#      Perr = acose(dot(delr, ll[j]));

#      // cross track error, as component of topocentric Perr
#      xtrk  = asin(sin(alpha) * sin(Perr));  // cross track magnitude, radians
#      xtrk /= de2ra;                         // degrees
#      xtrk *= sign;                          // left of track is positive
#      xtrk  = round(xtrk, 2);

#      // sum position error in squared degrees
#      Perr /= de2ra;
#      sum += Perr*Perr;

#      yy = (int)satx.yy - 2000;
#      day = (int)satx.doy;
#      hh = (int)satx.hr;
#      mm = (int)satx.mn;
#      ss = (int)((satx.ss + .0001) * 1000);
#      if(ss >= 60000)
#      {
#         ss = (int)mod(ss, 60000);
#         mm += 1;
#      }
#      if(mm >= 60)
#      {
#         mm = (int)mod(mm, 60);
#         hh += 1;
#      }
#      if(hh >= 24)
#      {
#         hh = (int)mod(hh, 24);
#         day += 1;
#      }


#      printf("(%2d) %04d  %02d%03d %02d%02d:%05d  %5.1f  %5.1f  %5.1f  %6.2f   %6.2f  %7.3f\n",
#      j + 1, (int)odata[j][3], yy, day, hh, mm, ss, az, el, asp, xtrk, delt, Perr);

#      // print fit to file
#      if(out)
#      {
#         FILE *fp;
#         fp = fopen(file, "a");
#         fprintf(fp, "(%2d) %04d  %02d%03d %02d%02d:%05d  %5.1f  %5.1f  %5.1f  %6.2f   %6.2f  %7.3f\n",
#           j + 1, (int)odata[j][3], yy, day, hh, mm, ss, az, el, asp, xtrk, delt, Perr);
#         fclose(fp);
#      }
#   }
#   printf("\nrms%12.5f\n", sqrt(sum / nobs));
# }

