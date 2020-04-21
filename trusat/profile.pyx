# cython: boundscheck=False, wraparound=False, nonecheck=False, language_level=3, cdivision=True, profile=True, embedsignature=True
from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values

"""
python satfit_accelerated_setup.py build_ext --inplace

satfit_caccelerated_profile.pyx - file containing only performance-profiling comparison code
    not used in production

In general, functions under profile test should call the fastest-known supporting functions,
or most-compatible (where the fastest is not compatible). This allows to test only the contribution
of the top-function to performance differences.

See caccalerated.pyx for optimization notes.
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

from . import caccelerated as tsc

from sgp4.propagation import sgp4, sgp4init
from sgp4.api import SGP4_ERRORS

import copy
from libc.math cimport fabs, cos, sin, M_PI, sqrt, fmod, acos, asin, atan, floor, modf

cimport cython
from cpython cimport array
import array
from numpy cimport ndarray
cimport numpy as np
import numpy as np

# Set up default array for (fast) array cloning
# https://cython.readthedocs.io/en/latest/src/tutorial/array.html
cdef array.array double_array_template = array.array('d')

np_double_array_template = np.zeros((3,), dtype='double')

from cython.parallel cimport prange
cimport openmp

import logging
log = logging.getLogger(__name__)

####### DECLARE GLOBAL CONSTANTS #######
cdef double TWOPI = 2*M_PI
cdef double NOCON = TWOPI/1440.0
cdef double DE2RA = M_PI/180.0


def vectorized_init_variable(double [:] jd):
    # First element is "1" as we're only propagating 1 satellite
    cdef Py_ssize_t t_size = len(jd)
    r = np.zeros((1, t_size, 3), dtype='double')
    v = np.zeros((1, t_size, 3), dtype='double')
    e = np.zeros((1, t_size),    dtype='uint8')

def vectorized_init_fixed(double [:] jd):
    # First element is "1" as we're only propagating 1 satellite
    cdef Py_ssize_t t_size = 100
    r = np.zeros((1, t_size, 3), dtype='double')
    v = np.zeros((1, t_size, 3), dtype='double')
    e = np.zeros((1, t_size),    dtype='uint8')

cpdef vectorized_init_npview(double [:] jd, double [:,:,:] rv_template):
    r = np.empty_like(rv_template)
    v = np.empty_like(rv_template)
    e = np.empty_like(jd)

cpdef np_zeros_like(double[:,:,:] nparr):
    return np.zeros_like(nparr)

cpdef np_empty_like(double[:,:,:] nparr):
    return np.empty_like(nparr)

cpdef np_empty_like_c(double[:,:,::1] nparr):
    return np.empty_like(nparr)

cpdef np_zeros_like_c(double[:,:,::1] nparr):
    return np.zeros_like(nparr)

cpdef tuple_to_array(tuple tup):
    cdef double[:] pyarr = array.clone(double_array_template, 3, zero=False)
    cdef Py_ssize_t i

    for i in range(3):
        pyarr[i] = tup[i]

    return pyarr

cpdef np_asarray(tuple tup):
    cdef double[:] pyarr = array.clone(double_array_template, 3, zero=False)
    return np.asarray(tup)

cpdef func_len(double [:] arr):
    return len(arr)

cpdef func_np_shape(double [:] arr):
    return arr.shape[0]

cpdef floor_floor(double f):
    return floor(f)

cpdef floor_mod(double f):
    return int(f)

# This performs the best, use directly in code instead of function call
cpdef modf_sub(double param):
    cdef double i, frac
    frac = modf(param, &i)
    return i, frac

cpdef modf_divmod(double param):
    cdef double i, frac
    (i, frac) = divmod(param, 1)
    return i, frac

cpdef divmod_raw(double param):
    cdef double i, frac
    i = floor(param)
    frac = param % 1
    return i, frac

cpdef divmod_sub(double param):
    cdef double i, frac
    i = floor(param)
    frac = param - i
    return i, frac


# Test-only example to show that prange is not faster for small loops
# Result - faster than everything else, but 17% slower than non-parallel version
cpdef dot_prange(double[:] v1, double[:] v2):
    """ sum = v1 . v2 """
    cdef double sum = 0.0
    cdef Py_ssize_t i

    for i in prange(3, nogil=True):
        sum += v1[i] * v2[i]
    return sum

# Optimization: cache a method call instead of calling it on the object
mag = np.linalg.norm

cpdef norm_py(double[:] v):
    """ ||v|| """
    return sqrt(tsc.dot(v, v))

cpdef mag_raw(double[:] x):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(x)
    cdef double sum = 0

    for i in range(n):
        sum += x[i]*x[i]
    return sqrt(sum)

cdef double mag_rawc(double[:] x):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(x)
    cdef double sum = 0

    for i in range(n):
        sum += x[i]*x[i]
    return sqrt(sum)
    
cpdef unit_vector_np(ndarray[np.float64_t, ndim=1] v):
    """ Returns the unit vector of the vector.  """
    u = v / tsc.norm(v)
    return u

def unit_vector_def(double[:] v):
    """ Returns the unit vector of the vector.  """
    cdef double magv
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v)

    u = array.clone(double_array_template, n, zero=False)

    magv = mag_rawc(v)
    for i in range(n):
        u[i] = v[i] / magv
    return u

cpdef unit_vector_ref(double[:] v, double[:] u):
    """ Returns the unit vector of the vector.  """
    cdef double magv
    cdef Py_ssize_t i

    magv = mag_rawc(v)
    for i in range(3):
        u[i] = v[i] / magv

cpdef npdot(double[:] v1, double[:] v2):
    return np.dot(v1,v2)

cpdef npcross(double[:] v1, double[:] v2):
    return np.cross(v1,v2)

def smult_py(a, v):
    """ av = a * v """
    av = np.zeros(3)
    av = a * v

cpdef smult_rtn(double a, double[:] v):
    return smult_rtnc(a,v)

cdef double[:] smult_rtnc(double a, double[:] v):
    """ av = a * v """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v)
    cdef double[:] av = array.clone(double_array_template, n, zero=False)

    for i in range(n):
        av[i] = a * v[i]
    return av

# Test-only example to show that prange is not faster for small loops
# Result - faster than everything else, but 20% slower than non-parallel version
cpdef smult_ref_prange(double a, double[:] v, double[:] av):
    """ av = a * v """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v)

    for i in prange(n, nogil=True):
        av[i] = a * v[i]

cpdef vmadd_rtn(double[:] v1, double[:] v2, double a):
    return vmadd_rtnc(v1, v2, a)

cdef double[:] vmadd_rtnc(double[:] v1, double[:] v2, double a):
    """ s = v1 + v2 """
    cdef Py_ssize_t n = len(v1)
    cdef Py_ssize_t i
    cdef double[:] s = array.clone(double_array_template, n, zero=False)

    for i in range(n):
        s[i] = v1[i] + a * v2[i]
    return s

# Test-only example to show that prange is not faster for small loops
# Result - faster than everything else, but 13% slower than non-parallel version
cpdef vmadd_ref_prange(double[:] v1, double[:] v2, double[:] s, double a):
    """ s = v1 + a * v2   (used for subtraction) """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(v1)

    for i in prange(n, nogil=True):  
        s[i] = v1[i] + a * v2[i]

cpdef cross_rtn(double[:] v1, double[:] v2):
    """ b = v1 x v2  """
    b = array.clone(double_array_template, 3, zero=False)

    b[0] = v1[1] * v2[2] - v1[2] * v2[1]
    b[1] = v1[2] * v2[0] - v1[0] * v2[2]
    b[2] = v1[0] * v2[1] - v1[1] * v2[0]

    return b

cpdef cross_ref(double[:] v1, double[:] v2, double[:] b):
    """ b = v1 x v2  """
    b[0] = v1[1] * v2[2] - v1[2] * v2[1]
    b[1] = v1[2] * v2[0] - v1[0] * v2[2]
    b[2] = v1[0] * v2[1] - v1[1] * v2[0]

cpdef proj_rtn(double[:] v2, double[:] v1):
    """ vp = unit vector projection of v1 onto v2 """
    cdef double b
    temp = array.clone(double_array_template, 3, zero=False)

    b = tsc.dot(v2, v1)/tsc.dot(v2, v2)
    temp = smult_rtn(b, v2)
    return tsc.unit_vector(temp)

cpdef proj_ref(double[:] v2, double[:] v1, double[:] vp):
    """ vp = unit vector projection of v1 onto v2 """
    cdef double b
    
    cdef double[:] temp = array.clone(double_array_template, 3, zero=False)
    # cdef double[:] temp = np.empty_like(np_double_array_template)

    b = tsc.dot(v2, v1)/tsc.dot(v2, v2)
    tsc.smult_ref(b, v2, temp)
    unit_vector_ref(temp, vp)

cpdef pythonmodulo(double a):
    "Note that within cython, this does not fold to positive values"
    return a % 360.0

cdef double posradangc(double a):
    """ Given a in radians, return equivalent circular value between 0 and 2 pi """
    if a < 0:
        return a + TWOPI
    elif a > TWOPI:
        return a - TWOPI
    else:
        return a

cdef double acosec_nogil(double x) nogil:
    if ( x >= 1.0 ):
        return 0.0
    if ( x <= -1.0 ):
        return  M_PI
    return acos(x)

cpdef delta_t_old(sat, double t):
    cdef double tsince, jd, fr, error

    # rr = array.clone(double_array_template, 3, zero=False)
    # vv = array.clone(double_array_template, 3, zero=False)
    # rr = np.empty_like(np_double_array_template)
    # vv = np.empty_like(np_double_array_template)

    tsince = (t - sat.jdsatepoch) * 1440.0 # time since epoch in minutes

    (rr, vv) = sgp4(sat,tsince) 

    rr = smult_rtn(1.0 / sat.radiusearthkm, np.asarray(rr)) # In Earth radii
    vv = smult_rtn(1.0 / (sat.radiusearthkm / 60.0), np.asarray(vv))  # In Earth radii / min - seriously!

    return rr, vv

cpdef delta_el_old(sat, xincl=False, xnodeo=False,   eo=False, omegao=False, xmo=False,      xno=False,   bsr=False, 
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
        sat.inclo = DE2RA*ii

    if (xnodeo):
        sat.nodeo = posradangc(xnodeo)
    elif (nodeo):
        sat.nodeo = posradangc(nodeo)
    elif (om):
        sat.nodeo = posradangc(DE2RA*om)

    if (eo):
        sat.ecco = eo   
    elif (ecco):
        sat.ecco = ecco   
    elif (ec):
        sat.ecco = ec   

    if (omegao):
        sat.argpo = posradangc(omegao)
    elif (argpo):
        sat.argpo = posradangc(argpo)
    elif (ww):
        sat.argpo = posradangc(DE2RA*ww)

    if (xmo):
        sat.mo = posradangc(xmo)
    elif (mo):
        sat.mo = posradangc(mo)
    elif (ma):
        sat.mo = posradangc(DE2RA*ma)

    if (xno):
        sat.no_kozai = xno
    elif (no_kozai):
        sat.no_kozai = no_kozai
    elif (nn):
        sat.no_kozai = nn * NOCON

    if (bsr):
        sat.bstar = bsr
    elif (bstar):
        sat.bstar = bstar

    if (jdsatepoch):
        sat.jdsatepoch = jdsatepoch
        jdSGP4epoch = sat.jdsatepoch - 2433281.5
    elif (jdSGP4epoch):
        sat.jdsatepoch = jdSGP4epoch + 2433281.5
       
    # FIXME HACK
    jdSGP4epoch = sat.jdsatepoch - 2433281.5

    sgp4init('wgs72', 'i', sat.satnum, jdSGP4epoch, sat.bstar, sat.ndot, sat.nddot, 
             sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo, sat)


cpdef delta_el_new(sat, xincl=False, xnodeo=False,   eo=False, omegao=False, xmo=False,      xno=False,   bsr=False, 
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
       
    sat.sgp4init(sat.satnum, jdSGP4epoch, sat_bstar, sat_ndot, sat_nddot,
                sat_ecco, sat_argpo, sat_inclo, sat_mo, sat_no_kozai, sat_nodeo)

    # Need to update this manually after initialization
    sat.jdsatepoch = jdsatepoch
    sat.jdsatepochF = jdsatepochF

    # sgp4init(sat.whichconst, sat.operationmode, sat.satnum, sat.jdSGP4epoch, sat.bstar, sat.ndot, sat.nddot, sat.ecco, sat.argpo, sat.inclo, sat.mo, sat.no_kozai, sat.nodeo, sat)
    # return sat

cpdef find_rms_old(satx, double[:,:] rd, double[:,:] ll, double[:,:] odata):
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    cdef Py_ssize_t nobs = len(rd)
    cdef double Perr = 0.0
    cdef double zum = 0.0

    cdef double[:] delr, rr, vv
    delr = array.clone(double_array_template, 3, zero=False)

    # TODO Improve performance of vectorized version in branch: 
    # https:#github.com/interplanetarychris/python-sgp4/tree/7-dec-15-vallado-tsince-vectorize
    for j in range(nobs):
        # advance satellite position
        (rr, vv) = delta_t_old(satx,odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # converted to unit vector
        tsc.vmadd_ref(rr, rd[j], delr, -1)
        delr = tsc.unit_vector(delr)

        # topocentric position error in degrees
        Perr = ( tsc.acose( tsc.dot(delr, ll[j]) ) )/DE2RA

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)

cpdef find_rmspy(satx, double[:,:] rd, double[:,:] ll, double[:,:] odata):
    """ find rms of observations against propagated TLE position

    Inputs:
        sat     Class variable of satellite elements
        odata       Observation data (numpy array)
        ll          Line of site vectors (numpy array)
        rd          Topocentric vectors (numpy array)

    Output:
        rms     RMS of observations against predict
    """
    cdef int nobs = len(rd)
    cdef double Perr = 0.0
    cdef double zum = 0.0
    cdef double[:] rr, vv

    cdef double[:] delr
    delr = array.clone(double_array_template, 3, zero=False)

    # TODO Improve performance of vectorized version in branch: 
    # https:#github.com/interplanetarychris/python-sgp4/tree/7-dec-15-vallado-tsince-vectorize
    for j in range(nobs):
        # advance satellite position
        (rr, vv) = delta_t_old(satx,odata[j][0])

        # predicted topocentric coordinates (sat xyz - observer xyz)
        # converted to unit vector
        tsc.vmadd_ref(rr, rd[j], delr, -1)
        delr = tsc.unit_vector(delr)

        # topocentric position error in degrees
        Perr = ( tsc.acose( tsc.dot(delr, ll[j]) ) )/DE2RA

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)


cpdef double find_rms_inline_scalar(satx, double[:,:] rd, double[:,:] ll, double[:,:] odata):
    return find_rms_inlinec_scalar(satx, rd, ll, odata)


cdef double find_rms_inlinec_scalar(satx, double[:,:] rd, double[:,:] ll, double[:,:] odata):
    """ find rms of observations against propagated TLE position

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
    cdef double jd, fr, Perr, error, temp
    cdef tuple rr, vv # Fixme, would really rather not have to convert tuple to array.array

    cdef double[:] rr_pyarr = array.clone(double_array_template, 3, zero=False)

    # Local temp variables for inline version
    cdef double magv, sum

    for j in range(nobs):
        # Get day, day-fraction parts of observation date
        fr = modf(odata[j][0], &jd)

        # advance satellite position
        (error, rr, vv) = satx.sgp4(jd, fr)

        # Fixme, would really rather not have to convert tuple to array.array
        # The following combines:
        # 1. Converting tuple to array.array
        # 2. Changing units to Earth radii
        # 3. delr = tsc.vmadd_rtn(rr_pyarr, rd[j], -1)
        # 4. tsc.norm(delr) (for later unit vector use)
        sum = 0.0
        for i in range(3):
            rr_pyarr[i] = (rr[i] / satx.radiusearthkm) - rd[j][i] 
            sum += rr_pyarr[i]*rr_pyarr[i]
        magv = sqrt(sum)

        # topocentric position error in degrees
        # 5. tsc.dot(delr, ll[j])
        sum = 0.0
        for i in range(3):
            sum += (rr_pyarr[i] / magv) * ll[j][i]

        # This one is not worth in-lining
        Perr = tsc.acose( sum ) / DE2RA

        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)


cpdef double find_rms_inline_prange(satx, double[:,:] rd, double[:,:] ll, double[:] jd, double[:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    return find_rms_inlinec_prange(satx, rd, ll, jd, fr, rr, vv, err)

cpdef double func_satx(satx):
    return satx.radiusearthkm

cpdef double func_sgp4(satx, double[:,:] rd, double[:,:] ll, double[:] jd, double[:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    satx._sgp4(jd, fr, err, rr, vv)
    return satx.radiusearthkm

cpdef double func_arrs(double[:,:] rd, double[:,:] ll, double[:] jd, double[:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    return jd[0]

cdef double find_rms_inlinec_prange(satx, double[:,:] rd, double[:,:] ll, double [:] jd, double [:] fr, double [:,:,:] rr, double [:,:,:] vv, unsigned char [:,:] err):
    """ find rms of observations against propagated TLE position

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
    cdef double zum = 0.0
    cdef Py_ssize_t nobs = len(rd)
    cdef double Perr, error, temp, radiusearthkm
    # cdef tuple rr, vv # Fixme, would really rather not have to convert tuple to array.array
    cdef int num_threads, thread_number

    num_threads = openmp.omp_get_num_threads()

    cdef double[:,:] rr_pyarr = np.zeros((num_threads,3), dtype='double')

    # Local temp variables for inline version
    cdef double magv, sum

    # Run vectorized sgp4 against input arrays (once)
    satx._sgp4(jd, fr, err, rr, vv)

    radiusearthkm = satx.radiusearthkm # Touch python variable once, isolated, outside loop
    for j in prange(nobs, nogil=True):
        thread_number = openmp.omp_get_thread_num()
        # advance satellite position
        # (error, rr, vv) = satx._sgp4(jd, fr)
        # Fixme, would really rather not have to convert tuple to array.array
        # The following combines:
        # 1. Converting tuple to array.array
        # 2. Changing units to Earth radii
        # 3. delr = tsc.vmadd_rtn(rr_pyarr, rd[j], -1)
        # 4. tsc.norm(delr) (for later unit vector use)
        sum = 0.0
        for i in range(3):
            rr_pyarr[thread_number][i] = (rr[0][j][i] / radiusearthkm) - rd[j][i] 
            sum = sum + rr_pyarr[thread_number][i]*rr_pyarr[thread_number][i]
        magv = sqrt(sum)

        # topocentric position error in degrees
        # 5. tsc.dot(delr, ll[j])
        sum = 0.0
        for i in range(3):
            sum = sum + (rr_pyarr[thread_number][i] / magv) * ll[j][i]

        # This one is not worth in-lining
        Perr = acosec_nogil( sum ) / DE2RA     

        Perr = 1
        # sum position error in squared degrees
        zum += Perr*Perr
    return sqrt(zum / nobs)

def move_epoch_to_jd_old(sat,t2_jd):
    (rr, vv) = delta_t_old(sat,t2_jd)
    # sat = delta_el(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t2_jd)
    delta_el_old(sat,inclo=sat.im, nodeo=sat.Om, ecco=sat.em, argpo=sat.om, mo=sat.mm, no_kozai=sat.no_kozai, jdsatepoch=t2_jd)
    # return sat