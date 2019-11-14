from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values
import sys

# As of 28 July 2019, python3.6 is the default "python3" in apt-get install python3 on Ubuntu
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version 3.6")
    sys.exit(1)

import os

# Use local/dev version of python-sgp4
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sgp4_path = os.path.join(parentdir, "python-sgp4")
sys.path.insert(1,sgp4_path) 

try:
    from sgp4.cpropagation import sgp4_scalar, sgp4_vector, sgp4init
except ImportError as e:
    print(e)
    from sgp4.propagation import sgp4_scalar, sgp4_vector, sgp4init

from sgp4.ext import jday, invjday
from datetime import datetime
from time import time
from math import fabs, cos, sin, pi, sqrt, fmod, acos, asin, atan
import copy
import numpy as np


import logging
log = logging.getLogger(__name__)

# ///////////// DECLARE GLOBAL VARIABLES ////////////////////////////////////////
twopi = 2*pi
nocon = twopi/1440.0
de2ra = pi/180.0

# Optimization: cache a method call instead of calling it on the object
mag = np.linalg.norm


def unit_vector(v):
    """ Returns the unit vector of the vector.  """
    u = v / mag(v)
    return u


def unit_vector_vec(v):
    """ Returns the unit vectors of the vector time series.  """
    u = np.apply_along_axis(unit_vector, 0, v)
    return u


def delta_t(sat,t):
    tsince = (t - sat.jdsatepoch) * 1440.0 # time since epoch in minutes

    (rr, vv) = sgp4_scalar(sat,np.array([tsince])) 

    sat.rr_km = rr
    sat.vv_kmpersec = vv

    sat.rr = np.asarray(rr) / sat.radiusearthkm # In Earth radii
    sat.vv = np.asarray(vv) / (sat.radiusearthkm / 60.0)  # In Earth radii / min - seriously!

    return sat


def delta_t_vec(sat,t):
    """ Perform propagation for a vector of jd time values
    Assume we're working with numpy array output from sgp4 
    Note we're changing to return r,v instead of sat
    """

    tsince = (t - sat.jdsatepoch) * 1440.0 # time since epoch in minutes

    (sat.rr_km, sat.vv_kmpersec) = sgp4_vector(sat,tsince) 
    if (sat.error):
        print(sat.error_message)
        return (False, False)
    else:
        sat.rr = sat.rr_km / sat.radiusearthkm # In Earth radii
        sat.vv = sat.vv_kmpersec / (sat.radiusearthkm / 60.0)  # In Earth radii / min - seriously!
        return sat.rr, sat.vv


def delta_el(sat, xincl=False, xnodeo=False,   eo=False, omegao=False, xmo=False,      xno=False,   bsr=False, 
                  inclo=False,  nodeo=False, ecco=False,  argpo=False,  mo=False, no_kozai=False,
                     ii=False,     om=False,   ec=False,     ww=False,  ma=False,       nn=False, bstar=False,  
            jdsatepoch=False, jdSGP4epoch=False, epoch_datetime=False):
    """ delta_el - Reinitialize SGP4 satellite object vi sgp4init() with changed Keplerian elements and/or epoch time
    
    Units can be passed in via radians or degrees, as is the convention for variable names
    """
    if (xincl):
        sat.inclo = xincl
    elif (inclo):
        sat.inclo = inclo
    elif (ii):
        sat.inclo = de2ra*ii

    if (xnodeo):
        sat.nodeo = posradang(xnodeo)
    elif (nodeo):
        sat.nodeo = posradang(nodeo)
    elif (om):
        sat.nodeo = posradang(de2ra*om)

    if (eo):
        sat.ecco = eo   
    elif (ecco):
        sat.ecco = ecco   
    elif (ec):
        sat.ecco = ec   

    if (omegao):
        sat.argpo = posradang(omegao)
    elif (argpo):
        sat.argpo = posradang(argpo)
    elif (ww):
        sat.argpo = posradang(de2ra*ww)

    if (xmo):
        sat.mo = posradang(xmo)
    elif (mo):
        sat.mo = posradang(mo)
    elif (ma):
        sat.mo = posradang(de2ra*ma)

    if (xno):
        sat.no_kozai = xno
    elif (no_kozai):
        sat.no_kozai = no_kozai
    elif (nn):
        sat.no_kozai = nn * nocon

    if (bsr):
        sat.bstar = bsr
    elif (bstar):
        sat.bstar = bstar

    if (jdsatepoch):
        sat.jdsatepoch = jdsatepoch
        sat.jdSGP4epoch = sat.jdsatepoch - 2433281.5
        sat.epoch_datetime = jday_to_datetime(sat.jdsatepoch)
    elif (jdSGP4epoch):
        sat.jdSGP4epoch = jdSGP4epoch
        sat.jdsatepoch = sat.jdSGP4epoch + 2433281.5
        sat.epoch_datetime = jday_to_datetime(sat.jdsatepoch)
    elif (epoch_datetime):
        sat.epoch_datetime = epoch_datetime

        (year, month, day, hour, minute, second) = sat.epoch_datetime.timetuple()[:6]
        microseconds = int(sat.epoch_datetime.strftime('%f'))
        sec_with_microseconds = second + microseconds/1.0E6

        sat.jdsatepoch = jday(year, month, day, hour, minute, sec_with_microseconds)
        sat.jdSGP4epoch = sat.jdsatepoch - 2433281.5
       
    sgp4init(sat.whichconst, sat.operationmode, sat.satnum, sat.jdSGP4epoch, sat.bstar, sat.ndot, sat.nddot, sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo, sat)
    return sat


def jday_to_datetime(jd):
    """ Returns a python datetime corresponding to a julian day """
    (yy, mm, dd, hr, mn, ss) = invjday(jd)
    (_, subsec) = divmod(ss,1)
    subsec = int(subsec*1E6)
    intss = int(ss)
    jday_datetime = datetime(yy,mm,dd,hr,mn,intss,subsec)    
    return jday_datetime


# TODO: accelerate
def posradang(a):
    """ Given a in radians, return equivalent circular value between 0 and 2 pi """
    if a < 0:
        return a + twopi
    elif a > twopi:
        return a - twopi
    else:
        return a


# TODO: accelerate
def rtw(ao, ac):
    """ round the world """
    if (fabs(ao - ac) > 180):
        if (ao < 180):
            ao += 360
        else:
            ac += 360
    return ao - ac


# TODO: accelerate
def acose(x):
    if ( x >= 1):
        rval = 0
    elif ( x <= -1):
        rval = pi
    else:
        rval = acos(x)
    return rval


def acose_vec(x):
    rval = np.where ( x >= 1, 0, np.where(x <= -1, pi, np.arccos(x)))
    return rval


def SGN(var):
    """ scott-campbell Legacy: Return sign of var """
    if (var<0):
        return -1
    else:
        return 1


# TODO: This function (and those it calls) will benefit the best from accelerating
# TODO: accelerate
def find_rms_scalar(satx, rd, ll, odata):
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    nobs = rd.shape[0] # Using this instead of len, as it might be better for Cython
    zum = 0

    # TODO Improve performance of vectorized version in branch: 
    # https://github.com/interplanetarychris/python-sgp4/tree/7-dec-15-vallado-tsince-vectorize
    for j in range(nobs):
        # advance satellite position
        satx = delta_t(satx,odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # converted to unit vector
        delr = unit_vector(satx.rr - rd[j])

        # topocentric position error in degrees
        Perr = ( acose( np.dot(delr, ll[j]) ) )/de2ra

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)


# TODO: accelerate
def find_rms_vector(satx, rd, ll, odata):
    """ find rms of observations against propagated TLE position
    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)
    Output:
        rms     RMS of observations against predict
    """
    # advance satellite position
    # Don't need the velocity component in this function
    (rrvec, _) = delta_t_vec(satx,odata[:,0])

    # predicted topocentric range (sat xyz - observer xyz)
    temp1 = rrvec - np.transpose(rd)
    delr = unit_vector_vec(temp1)

    temp2 = np.einsum('ij,ij->i',np.transpose(delr),ll)
    # topocentric position error in degrees
    Perr = acose_vec( temp2 ) / de2ra

    # sum position error in squared degrees
    rms = np.sqrt(Perr.dot(Perr)/Perr.size)

    return rms

# find_rms = find_rms_vector
find_rms = find_rms_scalar


# TODO: accelerate
def longitude(sat):
    """Calculate true longitude from mean anomaly and argument of perigee
    Inputs: 
      sat.ma     mean anomaly, radians
      sat.ecco   eccentricity
      sat.ww     argument of perigee, degrees

    Outputs:
      uu         True Longitude, degrees
    """

    ma = sat.mo
    ec = sat.ecco
    ww = sat.argpo / de2ra

    e = ma
    ek = 0
    while(fabs(e - ek) > 1e-6):
        ek = e
        e = ma + ec * sin(e)

    theta = (ec - cos(e)) / (ec * cos(e) - 1)
    theta = acose(theta)
    if (e > pi):
        theta = 2 * pi - theta
    uu = ww + theta / de2ra
    return uu


# TODO: accelerate
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
    dc = (asin(rll[2] / mag(rll)))/de2ra
    ra = (atan(rll[1] / rll[0]))/de2ra

    if (rll[0] < 0):
        ra += 180.0
    ra = fmod(ra, 360)

    return (ra, dc)


def rref(m): # scaled partial pivoting
    """ gaussian elimination """
    s = np.zeros(6)
    b = np.zeros(6)
    
    # calculate scale factors
    for i in range(6):
        s[i] = abs(m[i][0])
        for j in range(1, 6):
            if (s[i] < abs(m[i][j])):
                s[i] = abs(m[i][j])
    # end for i

    # swap rows according to scale
    for j in range (4):
        ROW = j
        for i in range (j + 1, 5):
            if (abs(m[ROW][j] / s[ROW]) < abs(m[i][j] / s[i])):
                ROW = i

        if (ROW != j):
            for k in range (j, 6):     # swap rows
                binv = m[j][k]
                m[j][k] = m[ROW][k]
                m[ROW][k] = binv

        binv = s[j]                 # swap scales
        s[j] = s[ROW]
        s[ROW] = binv
        # end if
    # end for j

    # # Alternate reference https://math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
    # # forward elimination 
    # for i in range(j + 1, 6):
    #   mult = m[i][j] / m[j][j]
    #   for k in range(j + 1, 7):
    #     m[i][k] = m[i][k] - mult * m[j][k]
    #   m[i][j] = 0
    # # end for i

    for j in range (6):
        mult = m[j][j]
        if (mult):
            for k in range (j, 7):
                m[j][k] = m[j][k] / mult
        for i in range(j + 1, 6):
            mult = m[i][j]
            for k in range (j + 1, 7):
                m[i][k] = m[i][k] - mult * m[j][k]
                m[i][j] = 0
        # end for i
    # end for j

    # test for singular matrix
    # Ref: https://stackoverflow.com/questions/13249108/efficient-pythonic-check-for-singular-matrix
    if np.linalg.cond(m) > 1/np.finfo(m.dtype).eps:
        log.error("Singular matrix")
        return [float("Inf")]

    # Alternate reference https://math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
    # back sustitution
    b[5] = m[5][6] / m[5][5]
    for i in range(4, 0, -1):
        bin = 0
        for k in range(i + 1, 5):
            bin = bin + m[i][k] * b[k]
        b[i] = (m[i][6] - bin) / m[i][i] # this line was giving an error at some points: "RuntimeWarning: divide by zero encountered in double_scalars"
    return b
# end rref


# TODO: accelerate
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
    xi = False
    xe = False
    xw = False
    xn = False

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
            dat[0][6] = rtw((odata[i][1])/de2ra, rac)      # store delta_ra
            dat[1][6] = ((odata[i][2])/de2ra) - dcc          # store delta_dc

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

    return sat


# TODO: accelerate
def anomaly_search(sat, rd, ll, odata, sum):
    """ mean anomaly box search, no limits """

    step = 0.1
    mk = sat.mo / de2ra
    min = 0
    max = 1 # Reference C++ source has the do loop evaluation at the end
    
    while (fabs(max - min) > 1e-5):
        min = mk
        max = mk

        # nsum loop - until rms doesn't decrease since last loop
        while (True):
            min = mk - step
            sat = delta_el(sat, ma=min)

            nsum = find_rms(sat, rd, ll, odata)
            if (nsum < sum):
                mk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop - until rms doesn't decrease since last loop
        while (True): 
            max = mk + step
            sat = delta_el(sat, ma=max)

            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                mk = max
                sum = xsum
                continue # Back to top of xsum loop
            break   
        step /= 2
    sat = delta_el(sat, ma=mk)

    return sat # Contains the ma at the end of the loop


# TODO: accelerate
def motion_search(sat, rd, ll, odata):
    """ mean motion box search, no limits """

    nk = sat.no_kozai/nocon
    sum = find_rms(sat, rd, ll, odata)

    # Start with this values to get through the loop once
    # Reference C++ source has the while evaluation at the end
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

            xsum = find_rms(sat, rd, ll, odata)
            if (xsum < sum):
                nk = max
                sum = xsum
                continue
            break
        step /= 2
    nn = nk
    return sat # nn (mean motion) is incorporated in the last update for the sat variable.


# TODO: accelerate
def node_search(satx, rd, ll, odata, sum, imax, imin, omax, omin):
    """ partition search on node and inclination within set limits """
    xi_set = False # unused here

    ii = satx.inclo / de2ra
    om = satx.nodeo / de2ra

    while((imax - imin) > 1e-5):
        istep = (imax - imin) / 20
        ostep = fabs(rtw(omax, omin) / 20)

        if (xi_set): # FIXME: to have the effect of running the for and while loops once for the fixed imin value
            imin  = ii
            imax  = ii + 1e-6
            istep = 10*imin

        for ik in np.arange(imin, imax, istep):
            for ok in np.arange(omin, omax, ostep):
                delta_el(satx, ii=ik, om=ok)

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

    satx = delta_el(satx, ii=ii, om=om)

    return satx


# TODO: accelerate
def perigee_search(sat, rd, ll, odata, sum, uu, wmax, wmin, emax, emin):
    """ partition search on perigee and eccentricity """

    xe = 0
    xn = 0
    xw = 0

    # Grab the values we're searching for in the loop, in case we don't find a new optimal
    ec  = sat.ecco
    ww  = sat.argpo/de2ra
    ma  = sat.mo/de2ra

    if (sat.ecco > 0.1):
        wmax = sat.argpo/de2ra + 0.1
        wmin = sat.argpo/de2ra - 0.1
        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99

    while((wmax - wmin) > 1e-5):
        estep = (emax - emin) / 20
        wstep = (wmax - wmin) / 20
        for wk in np.arange(wmin, wmax, wstep):
            if (xw):
                wmin  = sat.argpo/de2ra
                wmax  = sat.argpo/de2ra
                wk    = sat.argpo/de2ra
                wstep = 0
            theta = (uu - wk)*de2ra
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
                mk = (mk)/de2ra

                sat = delta_el(sat, ec=ek, ww=wk, ma=mk)

                rms = find_rms(sat, rd, ll, odata)

                if (rms < sum):
                    sum = rms
                    ec  = ek
                    ww  = wk
                    ma  = mk
            # END for ek
        # END for wk

        # Could save a call here by checking for the existence of the variables, but what's one more time?
        sat = delta_el(sat, ec=ec, ww=ww, ma=ma)

        wmax = sat.argpo/de2ra + wstep
        wmin = sat.argpo/de2ra - wstep
        emax = sat.ecco + estep
        emin = sat.ecco - estep
               
    # update mean_anomaly
    sat = anomaly_search(sat, rd, ll, odata, sum)

    # update mean_motion
    if (not xn):
        sat = motion_search(sat, rd, ll, odata)

    # calculate uu, degrees
    uu = longitude(sat)

    return sat


# TODO: accelerate
def step(sat, rd, ll, odata, sum, uu, step_type):       # partition search within limits set below
    nobs = odata.shape[0] # Using this instead of len, as it might be better for Cython
    last_rms = sum

    # first, update mean_anomaly 
    sat = anomaly_search(sat, rd, ll, odata, sum)

    emax = sat.ecco * 1.1
    emin = sat.ecco * 0.9
    wmax = sat.argpo / de2ra + 2
    wmin = sat.argpo / de2ra - 2
    imax = sat.inclo / de2ra + 2
    imin = sat.inclo / de2ra - 2
    omax = sat.nodeo / de2ra + 2
    omin = sat.nodeo / de2ra - 2

    if (step_type not in ["L","Z"]):
        print("\nPress Q to Quit    :\n")

    # Initialize loop progress variables
    xsum = 0  # Previous loop rms results
    lc   = 0  # Loop count
    DE   = 0  # Count of Differential Correction element loops
    stp_start = time() # Loop start time

    while( (fabs(sum-xsum)>1e-4) and lc <= 50 ):
        lc+=1
        xsum = sum
        ps_start = time() # Perigee Search start time
        sat = perigee_search(sat, rd, ll, odata, sum, uu, wmax, wmin, emax, emin)

        ns_start = time() # Node Search start time
        sat = node_search(sat, rd, ll, odata, sum, imax, imin, omax, omin)

        de_start = time() # Differential Correction of elements start time
        if (nobs > 3 and step_type == 'Z'):
            sat = diff_el(sat, rd, ll, odata, sum)
            de_stop = time()
            DE = de_stop - de_start

        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99
        wmax = sat.argpo / de2ra + 0.5
        wmin = sat.argpo / de2ra - 0.5
        imax = sat.inclo / de2ra + 0.5
        imin = sat.inclo / de2ra - 0.5
        omax = sat.nodeo / de2ra + 0.5
        omin = sat.nodeo / de2ra - 0.5

        # To save one call, would be nice if we just used the last calculated version from node_search or diff_el
        sum = find_rms(sat, rd, ll, odata)

        stp_lap = time()
        lap = stp_lap - ps_start
        PS = ns_start - ps_start
        NS = de_start - ns_start
        ELAPSED = stp_lap - stp_start
        print("rms{:12.5f}   Lap time: {:.2f}  PS {:.2f}  NS {:.2f}  DE {:.2f}  --  Elapsed {:.1f} / {}\t".format(sum, lap, PS, NS, DE, ELAPSED, lc),end='\r')
        sys.stdout.flush()

        if (step_type not in ["L","Z"]):
            print()
            buf = input("    : ")

            try:
                buf = buf.strip().upper()
                if (buf == 'Q'):
                    break
            except NameError:
                continue

    print()
    # sum = print_fit(sat, rd, ll, odata, last_rms)

    return sat


def move_epoch_to_jd(sat,t2_jd):
    sat = delta_t(sat,t2_jd)
    sat = delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t2_jd)
    return sat