# cython: boundscheck=False, wraparound=False, nonecheck=False, language_level=3, cdivision=True, profile=False, embedsignature=True
from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values


# python satfit_accelerated_setup.py build_ext --inplace

""" Optimization / Acceleration notes:
* Using cdef for functions that get called frequently, providing alias to them for non-Cython or performance testing
* Numpy is slow for 3-element vectors, raw python/cython operations are significantly (10x) faster
* Passing cython arrays by memory view is fastest (20x numpy ndarray)
* Passing scalars by reference appears to have no advantage
* Returning arrays by reference has *slight* (2-10%) advantage 
* A for loop gives (for generic length) has no penalty vs a one-line fixed function for length=3
* Accessing "Python" variables within the Satrec.* class comes at a small penalty, so where variables
*  are referenced within a loop, it is best to assign to a local variable (once) for reference
* OpenMP parallelization appears to offer no benefit (and is actually 40-100% slower)
* SGP4 (vectorized) is the primarily bottleneck in the current code-base, and appears to be optimized to its limit
* Memory view alloc/de-alloc occupies a non-trival amount (25%) of the compute time
"""
import sys

# As of 28 July 2019, python3.6 is the default "python3" in apt-get install python3 on Ubuntu
if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script requires Python version 3.6")
    sys.exit(1)


# # Use local/dev version of python-sgp4
# import os
# import inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sgp4_path = os.path.join(parentdir, "python-sgp4")
# sys.path.insert(1,sgp4_path) 

from sgp4.ext import invjday
from sgp4.api import SGP4_ERRORS, WGS72

from datetime import datetime
from time import time
import copy
from libc.math cimport fabs, cos, sin, M_PI, sqrt, fmod, acos, asin, atan, floor, modf

cimport cython
from cpython cimport array
import array
cimport numpy as np
import numpy as np

# Set up default array for (fast) array cloning
# https://cython.readthedocs.io/en/latest/src/tutorial/array.html
cdef array.array double_array_template = array.array('d')

np_double_array_template = np.zeros((3,), dtype=np.double)

import logging
log = logging.getLogger(__name__)

####### DECLARE GLOBAL CONSTANTS #######
cdef double TWOPI = 2*M_PI
cdef double NOCON = TWOPI/1440.0
cdef double DE2RA = M_PI/180.0
cdef Py_UNICODE OPSMODE = 'i'

def float_step(start, stop, step):
    """ Return evenly spaced values within a given interval. Does not include stop value. 
    """
    cdef num = start
    while num < stop:
        yield num
        num += step

cpdef dot(double[:] v1, double[:] v2):
    return dotc(v1,v2)

cdef double dotc(double[:] v1, double[:] v2):
    """ sum = v1 . v2 """
    cdef double sum = 0.0
    cdef Py_ssize_t i

    for i in range(3):
       sum += v1[i] * v2[i]
    return sum

cpdef norm(double[:] v):
    """ ||v|| """
    return sqrt(dotc(v, v))
    
cpdef unit_vector(double[:] v):
    return unit_vectorc(v)

cdef unit_vectorc(double[:] v):
    """ Returns the unit vector of the vector.  """
    cdef double magv
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v)

    u = array.clone(double_array_template, n, zero=False)

    magv = norm(v)
    for i in range(n):
        u[i] = v[i] / magv
    return u

cpdef smult_ref(double a, double[:] v, double[:] av):
    """ av = a * v """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v)

    for i in range(n):
        av[i] = a * v[i]

cpdef vmadd_ref(double[:] v1, double[:] v2, double[:] s, double a):
    """ s = v1 + a * v2   (used for subtraction) """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v1)

    for i in range(n):  
        s[i] = v1[i] + a * v2[i]

cpdef double posradang(double a):
    return posradangc(a)

cdef double posradangc(double a):
    """ Given a in radians, return equivalent circular value between 0 and 2 pi """
    # Since returns exit function, can eliminate overhead of if/elif/else
    if a < 0:
        return a + TWOPI
    if a > TWOPI:
        return a - TWOPI
    return a

cdef double rtw(double ao, double ac):
    """ round the world """
    if (fabs(ao - ac) > 180.0):
        if (ao < 180.0):
            ao += 360.0
        else:
            ac += 360.0
    return ao - ac

cpdef acose(double x):
    return acosec(x)

cdef double acosec(double x):
    # Since returns exit function, can eliminate overhead of if/elif/else
    if ( x >= 1.0 ):
        return 0.0
    if ( x <= -1.0 ):
        return  M_PI
    return acos(x)

# TODO: Switch over to C-math copysign()
cpdef SGN(double var):
    """ scott-campbell Legacy: Return sign of var """
    # Since returns exit function, can eliminate overhead of if/elif/else
    if (var < 0.0):
        return -1.0
    return 1.0

cpdef get_sgp4_vec_vars(double [:,:] odata):
    """ Convenience function to bridge "old variables" to new _sgp4 vectorized variables 
    
    Doesn't change any of the inputs, just reformats structure.
    """
    nobs = len(odata)
    # Format for _vectorized_sgp4
    jd = np.zeros((nobs,))
    fr = np.zeros((nobs,))

    for i in range(nobs):
        (jd[i], fr[i]) = divmod(odata[i][0],1)

    # Initialize result variables for vectorized sgp4
    rr  = np.zeros((1, nobs, 3), dtype='double')
    vv  = np.zeros((1, nobs, 3), dtype='double')
    err = np.zeros((1, nobs),    dtype='uint8')

    return jd, fr, rr, vv, err


def jday_to_datetime(jd, jdF):
    """ Returns a python datetime corresponding to a julian day """
    (yy, mm, dd, hr, mn, ss) = invjday(jd + jdF)
    (_, subsec) = divmod(ss,1)
    subsec = int(subsec*1E6)
    intss = int(ss)
    jday_datetime = datetime(yy,mm,dd,hr,mn,intss,subsec)    
    return jday_datetime


# 9x faster than delta_t_old
cpdef delta_t(sat, double t):
    cdef double tsince, jd, fr, error
    cdef Py_ssize_t i

    cdef double[:] rr_pyarr = array.clone(double_array_template, 3, zero=False)
    cdef double[:] vv_pyarr = array.clone(double_array_template, 3, zero=False)
    cdef double[:] rr_rtn = array.clone(double_array_template, 3, zero=False)
    cdef double[:] vv_rtn = array.clone(double_array_template, 3, zero=False)

    fr = modf(t, &jd)
    (error, rr, vv) = sat.sgp4(jd, fr) 

    for i in range(3):
        rr_pyarr[i] = rr[i]
        vv_pyarr[i] = vv[i]

    smult_ref(1.0 / sat.radiusearthkm, rr_pyarr, rr_rtn) # In Earth radii
    smult_ref(1.0 / (sat.radiusearthkm / 60.0), vv_pyarr, vv_rtn)  # In Earth radii / min - seriously!

    return rr_rtn, vv_rtn


cpdef delta_el(sat, xincl=False, xnodeo=False,   eo=False, omegao=False, xmo=False,      xno=False,   bsr=False, 
                  inclo=False,  nodeo=False, ecco=False,  argpo=False,  mo=False, no_kozai=False,
                     ii=False,     om=False,   ec=False,     ww=False,  ma=False,       nn=False, bstar=False,  
            jdsatepoch=False, jdSGP4epoch=False, epoch_datetime=False):
    """ delta_el - Reinitialize SGP4 satellite object vi sgp4init() with changed Keplerian elements and/or epoch time
    
    Units can be passed in via radians or degrees, as is the convention for variable names
    """
    cdef double jdsatepochF = 0

    # The older variables are currently not supported for delta_el
    sat_ndot  = sat.ndot
    sat_nddot = sat.nddot

    # Inclination
    if (xincl):
        sat_inclo = xincl
    elif (inclo):
        sat_inclo = inclo
    elif (ii):
        sat_inclo = DE2RA*ii
    else:
        sat_inclo = sat.inclo

    # Right Ascension of the Ascending Node (RAAN)
    if (xnodeo):
        sat_nodeo = posradangc(xnodeo)
    elif (nodeo):
        sat_nodeo = posradangc(nodeo)
    elif (om):
        sat_nodeo = posradangc(DE2RA*om)
    else:
        sat_nodeo = sat.nodeo

    # Eccentricity
    if (eo):
        sat_ecco = eo
    elif (ecco):
        sat_ecco = ecco
    elif (ec):
        sat_ecco = ec
    else:
        sat_ecco = sat.ecco

    # Argument of Perigee
    if (omegao):
        sat_argpo = posradangc(omegao)
    elif (argpo):
        sat_argpo = posradangc(argpo)
    elif (ww):
        sat_argpo = posradangc(DE2RA*ww)
    else:
        sat_argpo = sat.argpo

    # Mean anomaly
    if (xmo):
        sat_mo = posradangc(xmo)
    elif (mo):
        sat_mo = posradangc(mo)
    elif (ma):
        sat_mo = posradangc(DE2RA*ma)
    else:
        sat_mo = sat.mo

    # Mean motion
    if (xno):
        sat_no_kozai = copy.copy(xno)
    elif (no_kozai):
        sat_no_kozai = copy.copy(no_kozai)
    elif (nn):
        sat_no_kozai = nn * NOCON
    else:
        sat_no_kozai = sat.no_kozai

    if (bsr):
        sat_bstar = bsr
    elif (bstar):
        sat_bstar = bstar
    else:
        sat_bstar = sat.bstar

    # The following variables are not set by sgpinit, they need to be set separately
    if (jdsatepoch):
        jdsatepochF = jdsatepoch % 1
        jdsatepoch  = jdsatepoch - jdsatepochF
        jdSGP4epoch = (jdsatepoch + jdsatepochF) - 2433281.5

    elif (jdSGP4epoch):
        jdsatepoch  = floor(jdSGP4epoch + 2433281.5)
        jdsatepochF =      (jdSGP4epoch + 2433281.5) % 1
    else:
        jdsatepoch  = sat.jdsatepoch
        jdsatepochF = sat.jdsatepochF
       
    sat.sgp4init(WGS72, OPSMODE, sat.satnum, jdSGP4epoch, sat_bstar, sat_ndot, sat_nddot,
                sat_ecco, sat_argpo, sat_inclo, sat_mo, sat_no_kozai, sat_nodeo)

    ## Need to update this manually after initialization
    ## No longer needed - updated in sgp4init() https://github.com/brandon-rhodes/python-sgp4/pull/49#
    # sat.jdsatepoch = jdsatepoch
    # sat.jdsatepochF = jdsatepochF

    # sgp4init(sat.whichconst, OPSMODE, sat.whichconst, OPSMODE, sat.satnum, sat.jdSGP4epoch, sat.bstar, sat.ndot, sat.nddot, sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo, sat)
    # return sat

#TODO: Move this to satfit_caccelerated_profile.pyx when diff_el is ported
cdef delta_el_fast(sat):
    """ Update Satrec with new osculating elements, with variables set within Satrec
    """
    cdef double jdsatepoch  = sat.jdsatepoch
    cdef double jdsatepochF = sat.jdsatepochF
    cdef double jdSGP4epoch = (jdsatepoch + jdsatepochF) - 2433281.5

    sat.sgp4init(WGS72, OPSMODE, sat.satnum, jdSGP4epoch, sat.bstar, sat.ndot, sat.nddot,
                sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo)
    ## Need to update this manually after initialization
    ## No longer needed - updated in sgp4init() https://github.com/brandon-rhodes/python-sgp4/pull/49#
    # sat.jdsatepoch = jdsatepoch
    # sat.jdsatepochF = jdsatepochF


cdef delta_el_fast2(sat, double bstar, double ndot, double nddot,
                   double ecco, double argpo, double inclo, double mo, double no_kozai, double nodeo):
    """ Update Satrec with new osculating elements, with variables set within Satrec
    """
    cdef double jdsatepoch  = sat.jdsatepoch
    cdef double jdsatepochF = sat.jdsatepochF
    cdef double jdSGP4epoch = (jdsatepoch + jdsatepochF) - 2433281.5

    sat.sgp4init(WGS72, OPSMODE, sat.satnum, jdSGP4epoch, bstar, ndot, nddot,
                ecco, argpo, inclo, mo, no_kozai, nodeo)
    ## Need to update this manually after initialization
    ## No longer needed - updated in sgp4init() https://github.com/brandon-rhodes/python-sgp4/pull/49#
    # sat.jdsatepoch = jdsatepoch
    # sat.jdsatepochF = jdsatepochF


cpdef double find_rms(satx, double[:,:] rd, double[:,:] ll, double[:,:] odata):
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    cdef Py_ssize_t i, j
    cdef Py_ssize_t nobs = len(rd)
    cdef double zum = 0.0
    cdef double jd, fr, Perr, error, temp
    cdef tuple rr, vv # The scalar SGP4 function returns tuples

    cdef double[:] rr_pyarr = array.clone(double_array_template, 3, zero=False)
    cdef double[:] vv_pyarr = array.clone(double_array_template, 3, zero=False)
    cdef double[:] delr = array.clone(double_array_template, 3, zero=False)

    for j in range(nobs):
        # advance satellite position
        fr = modf(odata[j][0], &jd)
        (error, rr, vv) = satx.sgp4(jd, fr)

        for i in range(3):
            rr_pyarr[i] = rr[i]
            vv_pyarr[i] = vv[i]

        smult_ref(1.0 / satx.radiusearthkm, rr_pyarr, rr_pyarr) # In Earth radii
        smult_ref(1.0 / (satx.radiusearthkm / 60.0), vv_pyarr, vv_pyarr)  # In Earth radii / min - seriously!

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # converted to unit vector
        vmadd_ref(rr_pyarr, rd[j], delr, -1)
        delr = unit_vectorc(delr)

        # topocentric position error in degrees
        temp = dotc(delr, ll[j])
        Perr = ( acosec( temp ) )/DE2RA

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)


cpdef double find_rms_inline(satx, double[:,:] rd, double[:,:] ll, double[:] jd, double[:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    return find_rms_inlinec(satx, rd, ll, jd, fr, rr, vv, err)

cdef double find_rms_inlinec(satx, double[:,:] rd, double[:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    """ find rms of vectorized observations against propagated TLE position

    Ridiculously optimized as this is the most frequently called function

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    cdef Py_ssize_t i, j
    cdef Py_ssize_t nobs = len(rd)
    cdef double zum = 0.0
    cdef double Perr, error, temp, radiusearthkm

    cdef double[:] rr_pyarr = array.clone(double_array_template, 3, zero=False)

    # Local temp variables for inline version
    cdef double magv, sum

    # Run vectorized sgp4 against input arrays (once)
    satx._sgp4(jd, fr, err, rr, vv)

    radiusearthkm = satx.radiusearthkm # Touch python variable once, isolated, outside loop
    for j in range(nobs):
        # advance satellite position
        # (error, rr, vv) = satx._sgp4(jd, fr)
        # The following combines:
        # 1. Converting SGP4 tuple to array.array
        # 2. Changing units to Earth radii
        # 3. delr = vmadd_rtnc(rr_pyarr, rd[j], -1)
        # 4. norm(delr) (for later unit vector use)
        sum = 0.0
        for i in range(3):
            rr_pyarr[i] = (rr[0][j][i] / radiusearthkm) - rd[j][i] 
            sum += rr_pyarr[i]*rr_pyarr[i]
        magv = sqrt(sum)

        # topocentric position error in degrees
        # 5. dot(delr, ll[j])
        sum = 0.0
        for i in range(3):
            sum += (rr_pyarr[i] / magv) * ll[j][i]

        # This one is not worth in-lining
        Perr = acosec( sum ) / DE2RA

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)

cpdef double longitude(sat):
    return longitudec(sat)

cdef double longitudec(sat):
    """Calculate true longitude from mean anomaly and argument of perigee
    Inputs: 
      sat.ma     mean anomaly, radians
      sat.ecco   eccentricity
      sat.argpo  argument of perigee, radians

    Outputs:
      uu         True Longitude, degrees
    """
    cdef double uu, theta

    cdef double ma    = sat.mo
    cdef double ec    = sat.ecco
    cdef double argpo = sat.argpo 

    cdef double e = ma
    cdef double ek = 0
    while(fabs(e - ek) > 1e-6):
        ek = e
        e = ma + ec * sin(e)

    theta = (ec - cos(e)) / (ec * cos(e) - 1)
    theta = acosec(theta)
    if (e > M_PI):
        theta = TWOPI - theta

    uu = posradangc(argpo + theta)/DE2RA
    return uu


# TODO: accelerate
cdef zrll(satx, double[:] rd, double ra, double dc, double[:] satx_rr):
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
    cdef double[:] rll = np.empty_like(np_double_array_template)

    vmadd_ref(satx_rr, rd, rll, -1)

    # return ra and dc, degrees
    ra = ( atan(rll[1] / rll[0]) )/DE2RA
    dc = ( asin(rll[2] / norm(rll)) )/DE2RA

    ra = fmod(ra, 360)


cpdef rref(double [:,:] m): # scaled partial pivoting
    """ gaussian elimination """
    cdef double [:] s, b
    s = np.empty_like(np_double_array_template)
    b = np.empty_like(np_double_array_template)

    cdef Py_ssize_t i, j, k
    cdef double binv, mult, bin
    
    # calculate scale factors
    for i in range(6):
        s[i] = fabs(m[i][0])
        for j in range(1, 6):
            if (s[i] < fabs(m[i][j])):
                s[i] = fabs(m[i][j])
    # end for i

    # swap rows according to scale
    for j in range (4):
        ROW = j
        for i in range (j + 1, 5):
            if (fabs(m[ROW][j] / s[ROW]) < fabs(m[i][j] / s[i])):
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

    # # Alternate reference https:#math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
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
    if np.linalg.cond(m) > 1/sys.float_info.epsilon:
        log.error("Singular matrix")
        return [float("Inf")]

    # Alternate reference https:#math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
    # back sustitution
    b[5] = m[5][6] / m[5][5]
    for i in range(4, 0, -1):
        bin = 0
        for k in range(i + 1, 5):
            bin = bin + m[i][k] * b[k]
        b[i] = (m[i][6] - bin) / m[i][i] # this line was giving an error at some points: "RuntimeWarning: divide by zero encountered in double_scalars"
    return b
# end rref


# # TODO: accelerate
cpdef diff_el(satx, double [:,:] rd, double [:,:] ll, double [:,:] odata, double sum):
    """ differential correction of elements

    Variables:
        int i, j, k, c = 0;
        ra, dc;                 # right ascension, declination variables
        delta, el;              # small change amount
        mdata[150][7];          # define matrix for multiple regression
        dat[2][7];              # define output matrix of zpde utility
        b[6];                   # output deltas, sat position, velocity, b*
        rv[6][7];               # normal matrix
        rac, dcc, rms;
    """
    cdef Py_ssize_t i, j, k
    cdef double ra, dc, rac, dcc, deltar, rms, el

    # Set these to 0 instead of false?
    xi = False
    xe = False
    xw = False
    xn = False

    cdef int c = 0 
    dat = np.zeros((2,7))
    mdata = []
    b = array.clone(double_array_template, 6, zero=False)
    nobs = len(odata)

    # To make the compiler shut up for a bit
    ra = dc = 0
    # saty = satx
    # satz = satx
    # sat  = satx
    #loop:

    # delta_el_full template
    # delta_el_fast2(sat, bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo)

    cdef double bstar    = satx.bstar
    cdef double ndot     = satx.ndot
    cdef double nddot    = satx.nddot
    cdef double ecco     = satx.ecco
    cdef double argpo    = satx.argpo
    cdef double inclo    = satx.inclo
    cdef double mo       = satx.mo
    cdef double no_kozai = satx.no_kozai
    cdef double nodeo    = satx.nodeo

    (jd, fr, rr, vv, err) = get_sgp4_vec_vars(odata)


    while(True): # Forever loop (at least for 20 improvements)
        rv = np.zeros((6,7))
        # begin the main loop of the program
        for i in range(nobs):
            # satz = copy.deepcopy(sat)               # differential sat
            (satx_rr, _) = delta_t(satx,odata[i][0])

            # first establish the computed ra, dc, at jdo with no perturbations
            zrll(satx, rd[i], ra, dc, satx_rr)        # output ra, dc, degrees
            rac = ra                                  # store computed ra and dc
            dcc = dc

            # find the deltas and load into output matrix, dat
            dat[0][6] = rtw(odata[i][1]/DE2RA, rac)   # store delta_ra
            dat[1][6] = (odata[i][2]/DE2RA) - dcc     # store delta_dc

            # 6 steps to build differential correction matrix
            j = 0
            if (xi):
                dat[0][j] = .001
                dat[1][j] = .001
            else:
                delta = 0.001                         # change
                # el_copy = copy.copy(satx.inclo)     # store reference
                # satx.inclo += delta                 # delta element
                # delta_el_fast(satx)
                delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo, inclo + delta, mo, no_kozai, nodeo)
                (satx_rr, _) = delta_t(satx, odata[i][0]) # recalculate with perturbed element
                zrll(satx, rd[i], ra, dc, satx_rr)    # perturbed ra, dc
                # satx.inclo = el                     # restore reference
                dat[0][j] = rtw(ra, rac) / delta      # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 1
            delta = 0.001                             # change
            # el = satx.nodeo                         # store reference
            # satx.nodeo += delta                     # delta element
            # delta_el_fast(satx)
            delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo + delta)
            (satx_rr, _) = delta_t(satx, odata[i][0]) # recalculate with perturbed element # FIXME: python-SGP4
            zrll(satx, rd[i], ra, dc, satx_rr)        # perturbed ra, dc
            # satx.nodeo = el                         # restore reference
            dat[0][j] = rtw(ra, rac) / delta          # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

            # Results from this one are fairly different - dat[0][j] -453 vs -474
            j = 2
            if (xe):
                dat[0][j] = 0.00001
                dat[1][j] = 0.00001
            else:
                delta = 0.0001                        # change
                # el = satx.ecco                      # store reference
                # satx.ecco += delta                  # delta element
                # delta_el_fast(satx)
                delta_el_fast2(satx, bstar, ndot, nddot, ecco + delta, argpo, inclo, mo, no_kozai, nodeo)
                (satx_rr, _) = delta_t(satx, odata[i][0]) # recalculate with perturbed element # FIXME python-SGP4
                zrll(satx, rd[i], ra, dc, satx_rr)    # perturbed ra, dc
                # satx.ecco = el                      # restore reference
                dat[0][j] = rtw(ra, rac) / delta      # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 3
            if (xw):
                dat[0][j] = 0.001
                dat[1][j] = 0.001
            else:
                delta = 0.001                         # change
                # el = satx.argpo                     # store reference
                # satx.argpo += delta                 # delta element
                # delta_el_fast(satx)
                delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo + delta, inclo, mo, no_kozai, nodeo)
                (satx_rr, _) = delta_t(satx,odata[i][0])  # recalculate with perturbed element # FIXME python-SGP4
                zrll(satx, rd[i], ra, dc, satx_rr)    # perturbed ra, dc
                # satx.argpo = el                     # restore reference
                dat[0][j] = rtw(ra, rac) / delta      # perturbed - computed
                dat[1][j] = (dc - dcc) / delta

            j = 4
            delta = 0.001                             # change
            # el = satx.mo                            # store reference
            # satx.mo += delta                        # delta element
            # delta_el_fast(satx)
            delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo, inclo, mo + delta, no_kozai, nodeo)
            (satx_rr, _) = delta_t(satx,odata[i][0])  # recalculate with perturbed element
            zrll(satx, rd[i], ra, dc, satx_rr)        # perturbed ra, dc
            # satx.mo = el                            # restore reference
            dat[0][j] = rtw(ra, rac) / delta          # perturbed - computed
            dat[1][j] = (dc - dcc) / delta

            j = 5
            if (xn):
                dat[0][j] = 0.000001
                dat[1][j] = 0.000001
            else:
                delta = 0.00001                       # change
                # el = satx.no_kozai                  # store reference
                # satx.no_kozai += delta              # delta element
                # delta_el_fast(satx)
                delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai + delta, nodeo)
                (satx_rr, _) = delta_t(satx,odata[i][0])  # recalculate with perturbed element
                zrll(satx, rd[i], ra, dc, satx_rr)    # perturbed ra, dc
                # satx.no_kozai = el                  # restore reference
                dat[0][j] = rtw(ra, rac) / delta      # perturbed - computed
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

        # saty = copy.deepcopy(sat) 
        # test update components with deltas
        # satx.inclo    += b[0]*0.1
        # satx.nodeo    += b[1]*0.1
        # satx.ecco     += b[2]*0.1
        # satx.argpo    += b[3]*0.1
        # satx.mo       += b[4]*0.1
        # satx.no_kozai += b[5]*0.1
        # delta_el_fast(satx)

        delta_el_fast2(satx, bstar, ndot, nddot, 
                        ecco     + b[2]*0.1,
                        argpo    + b[3]*0.1,
                        inclo    + b[0]*0.1,
                        mo       + b[4]*0.1,
                        no_kozai + b[5]*0.1,
                        nodeo    + b[1]*0.1)
        # rms = find_rms_inline_scalar(satx, rd, ll, odata)
        rms = find_rms_inlinec(satx, rd, ll, jd, fr, rr, vv, err)
        if (rms < sum):
            sum = rms
            # update components with deltas
            inclo    += b[0]*0.1
            nodeo    += b[1]*0.1
            ecco     += b[2]*0.1
            argpo    += b[3]*0.1
            mo       += b[4]*0.1
            no_kozai += b[5]*0.1
            c+=1
            if (c < 20):
                continue # Back up to the top of the loop
            else:
                break
        else:
            break

    # delta_el_fast(satx)
    delta_el_fast2(satx, bstar, ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo)

cpdef anomaly_search(satx, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err, double sum):
    anomaly_searchc(satx, rd, ll, jd, fr, rr, vv, err, sum)

cdef anomaly_searchc(satx, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err, double sum):
    """ mean anomaly box search, no limits """
    cdef double nsum, xsum

    cdef double step = 0.1
    cdef double mk = satx.mo / DE2RA
    cdef double min = 0
    cdef double max = 1 # Reference C++ source has the do loop evaluation at the end

    while (fabs(max - min) > 1e-5):
        min = mk
        max = mk

        # nsum loop - until rms doesn't decrease since last loop
        while (True):
            min = mk - step
            # satx.mo=min*DE2RA
            delta_el_fast2(satx, satx.bstar, satx.ndot, satx.nddot,
                          satx.ecco, satx.argpo, satx.inclo, min*DE2RA, satx.no_kozai, satx.nodeo)
            # delta_el(satx, ma=min)

            nsum = find_rms_inlinec(satx, rd, ll, jd, fr, rr, vv, err)
            if (nsum < sum):
                mk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop - until rms doesn't decrease since last loop
        while (True): 
            max = mk + step
            # satx.mo=max*DE2RA
            delta_el_fast2(satx, satx.bstar, satx.ndot, satx.nddot,
                          satx.ecco, satx.argpo, satx.inclo, max*DE2RA, satx.no_kozai, satx.nodeo)
            # delta_el(satx, ma=max)

            xsum = find_rms_inlinec(satx, rd, ll, jd, fr, rr, vv, err)
            if (xsum < sum):
                mk = max
                sum = xsum
                continue # Back to top of xsum loop
            break   
        step /= 2
    # satx.mo = mk*DE2RA
    delta_el_fast2(satx, satx.bstar, satx.ndot, satx.nddot,
                  satx.ecco, satx.argpo, satx.inclo, mk*DE2RA, satx.no_kozai, satx.nodeo)
    # delta_el(satx, ma=mk)


cdef bstar_search(sat, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err, double sum): 
    cdef double bstep, bk, rms
    cdef double bstar = sat.bstar

    cdef double bmax = bstar * 1.1
    cdef double bmin = bstar * 0.9

    if (bstar < 0):
        bmax = bstar * 0.9
        bmin = bstar * 1.1

    while((bmax - bmin) > 1.e-10):
        bstep = (bmax - bmin) / 20
        for bk in float_step(bmin, bmax, bstep):
        
            # sat.bstar = bk
            # delta_el_fast(sat)
            delta_el_fast2(sat, bk, sat.ndot, sat.nddot,
                           sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo)
            # delta_el(sat,bstar=bk)

            # establish the computed ra, dc, at jdo with no perturbations
            rms = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)
            if (rms < sum):
                sum = rms
                bstar = bk
            # END for bk
        # sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar)
        bmin = bstar - bstep
        bmax = bstar + bstep
    # sat.bstar = bk
    # delta_el_fast(sat)
    delta_el_fast2(sat, bk, sat.ndot, sat.nddot,
                   sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo)


# TODO: accelerate
cdef motion_search(sat, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    """ mean motion box search, no limits """
    cdef double sum, nsum, xsum

    cdef double nk = sat.no_kozai/NOCON
    # FIXME: Maybe pass this in from calling function?
    sum = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)

    # Start with this values to get through the loop once
    # Reference C++ source has the while evaluation at the end
    cdef double min = 0
    cdef double max = 1
    cdef double step = 0.1
    while(fabs(max - min) > 1e-10):
        min = nk
        max = nk

        # nsum loop - until rms doesn't decrease since last loop
        while(True):
            min = nk - step
            # sat.no_kozai = min * NOCON
            # delta_el_fast(sat)
            delta_el_fast2(sat, sat.bstar, sat.ndot, sat.nddot,
                           sat.ecco, sat.argpo, sat.inclo, sat.mo, min * NOCON, sat.nodeo)
            # delta_el(sat, nn=min)

            nsum = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)
            if (nsum < sum):
                nk = min
                sum = nsum
                continue # back to top of nsum loop
            break # Go forward to xsum loop

        # xsum loop - until rms doesn't decrease since last loop
        while(True):
            max = nk + step
            # sat.no_kozai = max * NOCON
            # delta_el_fast(sat)
            delta_el_fast2(sat, sat.bstar, sat.ndot, sat.nddot,
                           sat.ecco, sat.argpo, sat.inclo, sat.mo, max * NOCON, sat.nodeo)
            # delta_el(sat, nn=max)

            xsum = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)
            if (xsum < sum):
                nk = max
                sum = xsum
                continue
            break
        step /= 2
    # return sat # nn (mean motion) is incorporated in the last update for the sat variable.


# # TODO: accelerate
cdef node_search(satx, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err, double sum, double imax, double imin, double omax, double omin):
    """ partition search on node and inclination within set limits """
    cdef bint xi_set = False # unused here
    cdef double istep, ik, ostep, ok, rms
    # cdef double rr[3] # Fix with real rr

    cdef double ii = satx.inclo / DE2RA
    cdef double om = satx.nodeo / DE2RA

    # Touch fixed python sat-variables outside loop, once
    # Saves about 0.2s per solution (not much, but multiplied by the whole catalog...)
    cdef double bstar    = satx.bstar
    cdef double ndot     = satx.ndot
    cdef double nddot    = satx.nddot
    cdef double ecco     = satx.ecco
    cdef double argpo    = satx.argpo
    cdef double mo       = satx.mo
    cdef double no_kozai = satx.no_kozai

    while((imax - imin) > 1e-5):
        istep = (imax - imin) / 20.0
        ostep = fabs(rtw(omax, omin) / 20.0)

        if (xi_set): # FIXME: to have the effect of running the for and while loops once for the fixed imin value
            imin  = ii
            imax  = ii + 1e-6
            istep = 10*imin

        for ik in float_step(imin, imax, istep):
            for ok in float_step(omin, omax, ostep):
                # satx.inclo = ik*DE2RA
                # satx.nodeo = ok*DE2RA
                # delta_el_fast(satx)
                delta_el_fast2(satx, bstar, ndot, nddot,
                               ecco, argpo, ik*DE2RA, mo, no_kozai, ok*DE2RA)
                # delta_el(satx, ii=ik, om=ok)

                # establish the computed ra, dc, at jdo with no perturbations
                rms = find_rms_inlinec(satx, rd, ll, jd, fr, rr, vv, err)
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

    # satx.inclo = ii*DE2RA
    # satx.nodeo = om*DE2RA
    # delta_el_fast(satx)
    delta_el_fast2(satx, bstar, ndot, nddot,
                   ecco, argpo, ii*DE2RA, mo, no_kozai, om*DE2RA)
    # delta_el(satx, ii=ii, om=om)


cdef perigee_search(sat, double [:,:] rd, double [:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err, double sum, double uu, double wmax, double wmin, double emax, double emin):
    """ partition search on perigee and eccentricity """
    cdef double estep, wstep, wk, ek, e, mk, rms

    cdef double xe = 0
    cdef double xn = 0
    cdef double xw = 0

    # Touch fixed python sat-variables outside loop, once
    cdef double bstar    = sat.bstar
    cdef double ndot     = sat.ndot
    cdef double nddot    = sat.nddot
    cdef double inclo    = sat.inclo
    cdef double no_kozai = sat.no_kozai
    cdef double nodeo    = sat.nodeo

    # Grab the values we're searching for in the loop, in case we don't find a new optimal
    cdef double ec  = sat.ecco
    cdef double ww  = sat.argpo/DE2RA
    cdef double ma  = sat.mo/DE2RA

    if (ec > 0.1):
        wmax = sat.argpo/DE2RA + 0.1
        wmin = sat.argpo/DE2RA - 0.1
        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99

    while((wmax - wmin) > 1e-5):
        estep = (emax - emin) / 20
        wstep = (wmax - wmin) / 20
        for wk in float_step(wmin, wmax, wstep):
            if (xw):
                wmin  = sat.argpo/DE2RA
                wmax  = sat.argpo/DE2RA
                wk    = sat.argpo/DE2RA
                wstep = 0
            theta = (uu - wk)*DE2RA
            # print(f"emin {emin}  emax {emax}  estep {estep}")
            for ek in float_step(emin, emax, estep):
                if (xe):
                    emin  = sat.ecco
                    emax  = sat.ecco
                    ek    = sat.ecco
                    estep = 0
                e = acosec((ek + cos(theta)) / (1 + ek * cos(theta)))
                if (theta > M_PI):
                    e = TWOPI - e
                mk = e - ek * sin(e)
                mk = (mk)/DE2RA

                # sat = delta_el(sat, ec=ek, ww=wk, ma=mk)
                # sat.ecco  = ek
                # sat.argpo = wk*DE2RA
                # sat.mo    = mk*DE2RA 
                delta_el_fast2(sat, bstar, ndot, nddot,
                              ek, wk*DE2RA, inclo, mk*DE2RA, no_kozai, nodeo)
                # delta_el_fast(sat)
                # delta_el(sat, ec=ek, ww=wk, ma=mk)

                rms = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)

                if (rms < sum):
                    sum = rms
                    ec  = ek
                    ww  = wk
                    ma  = mk
            # END for ek
        # END for wk

        # Could save a call here by checking for the existence of the variables, but what's one more time?
        # sat.ecco  = ec
        # sat.argpo = ww*DE2RA
        # sat.mo    = ma*DE2RA 
        # delta_el_fast(sat)
        delta_el_fast2(sat, bstar, ndot, nddot,
                       ec, ww*DE2RA, inclo, ma*DE2RA, no_kozai, nodeo)
        # delta_el(sat, ec=ec, ww=ww, ma=ma)

        wmax = sat.argpo/DE2RA + wstep
        wmin = sat.argpo/DE2RA - wstep
        emax = sat.ecco + estep
        emin = sat.ecco - estep
               
    # update mean_anomaly
    anomaly_searchc(sat, rd, ll, jd, fr, rr, vv, err, sum)

    # update mean_motion
    if (not xn):
        motion_search(sat, rd, ll, jd, fr, rr, vv, err)

    # calculate uu, degrees
    uu = longitudec(sat)


cpdef step(sat, double [:,:] rd, double [:,:] ll, double [:,:] odata, double sum, double uu, unicode step_type):       
    stepc(sat, rd, ll, odata, sum, uu, step_type)       


# # TODO: accelerate
cdef stepc(sat, double [:,:] rd, double [:,:] ll, double [:,:] odata, double sum, double uu, unicode step_type):       
    """ partition search within limits set below 
    """
    cdef int nobs = len(odata)
    cdef double last_rms = sum
    # cdef double rr[3] # Fix with real rr

    # first, update mean_anomaly 
    # anomaly_searchc(sat, rd, ll, odata, sum)
    # bstar_search(sat, rd, ll, odata, sum)

    cdef double emax = sat.ecco * 1.1
    cdef double emin = sat.ecco * 0.9
    cdef double wmax = sat.argpo / DE2RA + 2
    cdef double wmin = sat.argpo / DE2RA - 2
    cdef double imax = sat.inclo / DE2RA + 2
    cdef double imin = sat.inclo / DE2RA - 2
    cdef double omax = sat.nodeo / DE2RA + 2
    cdef double omin = sat.nodeo / DE2RA - 2
    # print(f"emin {emin}  emax {emax}")

    if (step_type not in ["L","Z"]):
        print("\nPress Q to Quit    :\n")

    # Initialize loop progress variables
    cdef double xsum = 0  # Previous loop rms results
    cdef int lc   = 0  # Loop count
    cdef int reversal = 0 # Dither detection
    cdef double DE   = 0  # Count of Differential Correction element loops
    cdef double stp_start = time() # Loop start time
    cdef double ps_start, ns_start, de_start, de_stop, stp_lap, PS, NS, ELAPSED
    cdef unicode buf

    cdef double[:] satx_rr = np.empty_like(np_double_array_template)

    # Parse out observation time data so we only need manipulate it once
    # Format for _vectorized_sgp4
    cdef double[:] jd = np.zeros((nobs,))
    cdef double[:] fr = np.zeros((nobs,))

    for i in range(nobs):
        jd[i] = odata[i][0] // 1.0
        fr[i] = odata[i][0]  % 1.0

    # Initialize result variables for vectorized sgp4
    cdef double[:,:,:] rr       = np.zeros((1, nobs, 3), dtype='double')
    cdef double[:,:,:] vv       = np.zeros((1, nobs, 3), dtype='double')
    cdef unsigned char[:,:] err = np.zeros((1, nobs),    dtype='uint8')

    # while( (fabs(sum-xsum)>1e-4) and lc <= 50 ):
    while( (fabs(sum-xsum)>1e-5) and lc < 1000 ):
        # Trying updating these in the loop, as it appears to further optimize the result (not sure about bstar)
        anomaly_searchc(sat, rd, ll, jd, fr, rr, vv, err, sum)
        # bstar_search(sat, rd, ll, odata, sum)

        lc+=1
        xsum = sum
        ps_start = time() # Perigee Search start time
        perigee_search(sat, rd, ll, jd, fr, rr, vv, err, sum, uu, wmax, wmin, emax, emin)

        ns_start = time() # Node Search start time
        node_search(sat, rd, ll, jd, fr, rr, vv, err, sum, imax, imin, omax, omin)

        if (sat.bstar != 0):
            bstar_search(sat, rd, ll, jd, fr, rr, vv, err, sum)

        de_start = time() # Differential Correction of elements start time
        if (nobs > 3 and step_type == 'Z'):
            diff_el(sat, rd, ll, odata, sum)
            de_stop = time()
            DE = de_stop - de_start

        emax = sat.ecco * 1.01
        emin = sat.ecco * 0.99
        wmax = sat.argpo / DE2RA + 0.5
        wmin = sat.argpo / DE2RA - 0.5
        imax = sat.inclo / DE2RA + 0.5
        imin = sat.inclo / DE2RA - 0.5
        omax = sat.nodeo / DE2RA + 0.5
        omin = sat.nodeo / DE2RA - 0.5

        # To save one call, would be nice if we just used the last calculated version from node_search or diff_el
        sum = find_rms_inlinec(sat, rd, ll, jd, fr, rr, vv, err)

        stp_lap = time()
        lap = stp_lap - ps_start
        PS = ns_start - ps_start
        NS = de_start - ns_start
        ELAPSED = stp_lap - stp_start

        # Figure out how to accelerate this piece?
        print_string = "rms{:12.5f}   Lap time: {:.2f}  PS {:.2f}  NS {:.2f}  DE {:.2f}  --  Elapsed {:.1f} / {}\t".format(sum, lap, PS, NS, DE, ELAPSED, lc)
        print(print_string,end='\r')
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
        else:
            if (sum > xsum):
                reversal += 1
                print("Dithering... {:2d}".format(reversal))
                if (reversal > 10):
                    break
    print(print_string)

def move_epoch_to_jd(sat,t2_jd):
    (satrr, satvv) = delta_t(sat,t2_jd)
    delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t2_jd)
