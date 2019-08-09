# elfind.py
# Adapted from C++ version originally by Scott Campbell campbel7@the-i.net

from __future__ import print_function
from __future__ import division         # Eliminate need for decimals on whole values
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 7:
	print("This script requires Python version 3.7")
	sys.exit(1)

import configparser                 # config file parsing
import argparse                     # command line parsing
from math import (acos, asin, atan, atan2, cos, fmod,  # Fast/precise math functions
                  pi, pow, radians, sin, sqrt, tan)       
import numpy as np
import string

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

from skyfield.api import Topos, load
import iod
import tle_util

from iod_analyze import CosparSite
from satid_skyfield import mag, unit_vector

# Global variables (for now, work to eliminate these)
_XKMPER = 6378.137    # WGS84 Earth Equatorial Radius
CONVERGE = 1.0e-13    # Convergence tolerance - originally 1.0e-13 

# ///////////////// PHYSICAL CONSTANTS //////////////////////////////////////////

# /* dimensions & gravity of earth, World Geodetic System 1972 values */

de2ra = radians(1)
tu = .0093380977083333   # time unit in days ~ 13 minutes
twopi = 2*pi
xj3 = -2.53881E-6
ck2 = 5.413079E-4
xke = 7.43669161E-2 # (G*M)^(1/2)*(er/min)^(3/2) where G = 
                    # Newton's grav const, M = earth mass 

# osculating orbit elements
ssn = None
xincl = None    # Inclination
xnodeo = None
omegao = None
xmo = None
xno = None
epoch_datetime = None
# FIXME fix these global variables to a better container

def fmod2p(x):
  rval = fmod(x, twopi)

  if(rval < 0):
    rval += 2*pi
  return(rval)

def rv2el(rr2, vv2):
  """ vectors to mean elements """
  # classical osculating orbit elements calculated from vectors rr2, vv2
  # double xinck, xnodek, ek, mk, wk, xn, rk, uk, aodp, pl, rdotk, rfdotk, temp;

  global xincl
  global xnodeo
  global omegao
  global xno     # mean motion
  global xmo     # mean anomaly
  global eo      # eccentricity
  global omegao  # mean perigee
  # FIXME fix these global variables to a better container

  vk = (1 / xke)*vv2
  h = np.cross(rr2, vk)
  pl = np.dot(h, h)
  vz = np.array([0, 0, 1])
  n = np.cross(vz, h)
  rk = mag(rr2)
  rdotk = np.dot(rr2, vv2) / rk
  rfdotk = mag(h) * xke / rk
  temp = np.dot(rr2, n) / rk / mag(n)
  if (abs(temp) > 1):
    temp = SGN(temp)
  uk = acos(temp)
  if(rr2[2] < 0):
    uk = twopi - uk
  vz = np.cross(vk, h)
  vy = (-1 / rk) * rr2
  vec = np.add(vz,vy)
  ek = mag(vec)
  temp = h[2] / mag(h)
  if (abs(temp) > 1):
    temp = SGN(temp)
  xinck =  acos(temp)
  temp = n[0] / mag(n)
  if (abs(temp) > 1):
    temp = SGN(temp)
  xnodek = acos(temp)
  if (n[1] < 0):
    xnodek = fmod2p(twopi - xnodek)
  temp = np.dot(vec, n) / ek / mag(n)
  if (abs(temp) > 1):
    temp = SGN(temp)

  wk = acos(temp)
  if (vec[2] < 0):
    wk = fmod2p(twopi - wk)

  aodp = pl / (1 - ek*ek)

  xn = xke * pow(aodp, -1.5)

  """
  In the first loop the osculating elements rk, uk, xnodek, xinck, rdotk,
  and rfdotk are used as anchors to find the corresponding final SGP4
  mean elements r, u, xnodeo, xincl, rdot, and rfdot.  Several other final
  mean values based on these are also found: betal, cosio, sinio, theta2,
  cos2u, sin2u, x3thm1, x7thm1, x1mth2.  In addition, the osculating values
  initially held by aodp, pl, and xn are replaced by intermediate
  (not osculating and not mean) values used by SGP4.  The loop converges
  on the value of pl in about four iterations.
  """

  # seed value for first loop
  xincl = xinck
  u = uk

  for i in range(0, 99):
    a2 = pl
    betal = sqrt(pl / aodp)
    temp1 = ck2  / pl
    temp2 = temp1 / pl
    cosio = cos(xincl)
    sinio = sin(xincl)
    sin2u = sin(2*u)
    cos2u = cos(2*u)
    theta2 = cosio * cosio
    x3thm1 = 3 * theta2 - 1
    x1mth2 = 1 - theta2
    x7thm1 = 7 * theta2 - 1
    r = (rk - .5 * temp1 * x1mth2 * cos2u) \
         / (1. - 1.5 * temp2 * betal * x3thm1)
    u = uk + .25 * temp2 * x7thm1 * sin2u
    xnodeo = xnodek - 1.5 * temp2 * cosio * sin2u
    xincl = xinck - 1.5 * temp2 * cosio * sinio * cos2u
    rdot = rdotk + xn * temp1 * x1mth2 * sin2u
    rfdot = rfdotk - xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)
    temp = r * rfdot / xke
    pl = temp * temp

    # vis-viva equation
    temp = 2 / r - (rdot*rdot + rfdot*rfdot) / (xke*xke)
    aodp = 1 / temp

    xn = xke * pow(aodp, -1.5)
    if (abs(a2 - pl) < CONVERGE):
      log.debug("(pl) converged after {} iterations.".format(i))
      break

  """
  The next values are calculated from constants and a combination of mean
  and intermediate quantities from the first loop.  These values all remain
  fixed and are used in the second loop.
  """

  # preliminary values for the second loop
  ecose = 1 - r / aodp
  esine = r * rdot / (xke * sqrt(aodp))   # needed for Kepler's eqn.
  elsq = 1 - pl / aodp  # intermediate eccentricity squared
  a3ovk2 = -xj3 / ck2
  xlcof = .125 * a3ovk2 * sinio * (3 + 5 * cosio) \
        / (1 + cosio)
  aycof = .25 * a3ovk2 * sinio
  temp1 = esine / (1 + sqrt(1 - elsq))
  cosu = cos(u)
  sinu = sin(u)

  """
  The second loop normally converges in about six iterations to the final
  mean value for the eccentricity, eo.  The mean perigee, omegao, is also
  determined.  Cosepw and sinepw are found to twelve decimal places and
  are used to calculate an intermediate value for the eccentric anomaly,
  temp2.  Temp2 is then used in Kepler's equation to find an intermediate
  value for the true longitude, capu.
  """

  # seed values for loop 
  eo = sqrt(elsq)
  omegao = wk
  axn = eo * cos(omegao)

  for i in range(0, 99):
    a2 = eo
    beta = 1 - eo*eo
    temp = 1 / (aodp * beta)
    aynl = temp * aycof
    ayn = eo * sin(omegao) + aynl
    cosepw = r * cosu / aodp + axn - ayn * temp1
    sinepw = r * sinu / aodp + ayn + axn * temp1
    axn = cosepw * ecose + sinepw * esine
    ayn = sinepw * ecose - cosepw * esine
    omegao = fmod2p(atan2(ayn - aynl, axn))
    eo = axn / cos(omegao)
    if (abs(a2 - eo) < CONVERGE):
      log.debug("eccentricity (eo) and mean perigee (omegao) converged after {} iterations.".format(i))
      break   # eccentricity (eo) and mean perigee (omegao) exist at this point

  temp2 = atan2(sinepw, cosepw)
  capu = temp2 - esine             # Kepler's equation
  xll = temp * xlcof * axn

  # xll adjusts the intermediate true longitude,
  # capu, to the mean true longitude, xl          
  xl = capu - xll

  xmo = fmod2p(xl - omegao) # mean anomaly (xmo)

  """
  The third loop usually converges after three iterations to the
  mean semi-major axis, a1, which is then used to find the mean motion, xno.
  """

  a0 = aodp
  a1 = a0
  beta2 = sqrt(beta)
  temp = 1.5 * ck2 * x3thm1 / (beta * beta2)
  for i in range(0, 99):
    a2 = a1
    d0 = temp / (a0*a0)
    a0 = aodp * (1 - d0)
    d1 = temp / (a1*a1)
    a1 = a0 / (1 - d1 / 3 - d1*d1 - 134 * d1*d1*d1 / 81)
    if (abs(a2 - a1) < CONVERGE):
      log.debug("mean motion (xno), and semi-major axis (a1) converged after {} iterations.".format(i))
      break
  xno = pow(a1 , -1.5) / (tu * twopi) # Solution for mean motion (xno), and semi-major axis (a1) exist at this point
# end rv2el

def so2r(r, rd, ll):
  """ find topocentric vector, rr, to satellite 
  
  given the line of sight unit vector, ll, from the observer to the satellite, 
  the topocentric vector, rd, to the observer, and 
  the length of the vector, r.
  """

  nrd = mag(rd)
  ang1 = acos(np.dot(rd, ll) / nrd)

  if (ang1 < .001):
      rho = r - nrd
  else:
    ang2 = asin(nrd * sin(ang1) / r)
    rho = r * sin(ang1 - ang2) / sin(ang1)

  return np.add(rd, rho*ll)


def so2rv(od, rd, ll):
  """ fit circular orbit to two observations 
  
  Calculates rr2, vv2 and passes them onto rv2el
  """
  rro = 0
  rr = 1.1

  delt = od[1][0] - od[0][0]

  # use the first two positions in the array
  theta = acos(np.dot(ll[0], ll[1]))
  rrx = np.subtract(ll[1], ll[0])
  rry = np.cross(rd[1], rrx)
  sin_phi = mag(rry) / (mag(rd[1])*mag(rrx))
  sin_phi = sin_phi * sin_phi   # squared

  # initial estimate for rr
  # TODO: Figure out why this loop runs twice, and what 7 and 53.575 are.
  for _ in range (0, 2):
     rr = 1 + delt * sin_phi * 53.575 / (tan(theta / 2) * sqrt(rr))
     if (rr > 7):
        rr = 7

  while (abs(rr - rro) > 1.0E-8):
    rr1 = so2r(rr, rd[0], ll[0])
    rr2 = so2r(rr, rd[1], ll[1])
    theta = acos(np.dot(rr1, rr2) / (rr*rr))
    vv = theta * rr * tu / delt
    rro = 1 / (vv*vv)
    rr = .01 * (99 * rr + rro)  # weighted average 100:1

    if (rr > 8):
       rr = 8
       break
    if (rr < 1.002):
       rr = 1.002
       break

  rrx = np.cross(rr1, rr2)  # placeholder for intermediate calculation
  rry = np.cross(rrx, rr2)  # placeholder for intermediate calculation
  rrx = unit_vector(rry)    # placeholder for intermediate calculation
  vv2 = vv*rrx

  vv2 = xke * vv2
  rv2el(rr2, vv2)
  return(rr2, vv2)


def rref(m): # scaled partial pivoting
  """ gaussian elimination """
  b = []
  s = []
  # calculate scale factors
  for i in range (0, 5):
    s[i] = abs(m[i][0])
    for j in range(1, 5):
      if (s[i] < abs(m[i][j])):
        s[i] = abs(m[i][j])
  # end for i

  # swap rows according to scale
  for j in range (0, 4):
    ROW = j
    for i in range (j + 1, 5):
      if (abs(m[ROW][j] / s[ROW]) < abs(m[i][j] / s[i])):
        ROW = i

    if (ROW != j):
      for k in range (j, 6):     # swap rows
        bin = m[j][k]
        m[j][k] = m[ROW][k]
        m[ROW][k] = bin

      bin = s[j]                 # swap scales
      s[j] = s[ROW]
      s[ROW] = bin
    # end if

    # Alternate reference https://math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
    # forward elimination 
    for i in range(j + 1, 5):
      mult = m[i][j] / m[j][j]
      for k in range(j + 1, 6):
        m[i][k] = m[i][k] - mult * m[j][k]
      m[i][j] = 0
    # end for i
  # end for j

  # test for singular matrix
  # Ref: https://stackoverflow.com/questions/13249108/efficient-pythonic-check-for-singular-matrix
  if np.linalg.cond(m) > 1/np.finfo(m.dtype).eps:
    log.error("Singular matrix")
    sys.exit()
   
  # Alternate reference https://math.stackexchange.com/questions/2950727/gaussian-elimination-in-numerical
  # back sustitution
  b[5] = m[5][6] / m[5][5]
  for i in range(4, 0, -1):
    bin = 0
    for k in range(i + 1, 5):
      bin = bin + m[i][k] * b[k]
    b[i] = (m[i][6] - bin) / m[i][i]
  return b
# end rref


def f8g(rr1, vv1, delt):
  """ F and G series """
  f = []
  g = []

  r = mag(rr1)
  u = 1 / (r*r*r)
  p = np.dot(rr1, vv1) / (r*r)
  q = np.dot(vv1, vv1) / (r*r) - u
  p2 = p * p
  p4 = p2 * p2
  q2 = q * q
  u2 = u * u
  f[0] = 1
  f[1] = 0
  f[2] = -u / 2
  f[3] = p * u / 2
  f[4] = u * (u - 3 * (5 * p2 - q)) / 24
  f[5] = -p * u * (u - 7 * p2 + 3 * q) / 8
  f[6] = -u * (u2 - 6 * (35 * p2 - 4 * q) * u \
         + 45 * (21 * p4 - 14 * p2 * q + q2)) / 720
  f[7] = p * u * (u2 - 2 * (25 * p2 - 7 * q) * u \
         + 5 * (33 * p4 - 30 * p2 * q + 5 * q2)) / 80
  f[8] = u * (u2*u - 9 * (245 * p2 - 13 * q) * u2 \
         + 27 * (1925 * p4 - 910 * p2 * q + 41 * q2) * u \
         - 315 * (429 * p4*p2 - 495 * p4 * q \
         + 135 * p2 * q2 - 5 * q2*q)) / 40320
  g[0] = 0
  g[1] = 1
  g[2] = 0
  g[3] = -u / 6
  g[4] = p * u / 4
  g[5] = u * (u - 9 * (5 * p2 - q)) / 120
  g[6] = -p * u * (u - 2 * (7 * p2 - 3 * q)) / 24
  g[7] = -u * (u2 - 18 * (35 * p2 - 3 * q) * u \
         + 225 * (21 * p4 - 14 * p2 * q + q2)) / 5040
  g[8] = p * u * (u2 - 4 * (25 * p2 - 6 * q) * u \
         + 15 * (33 * p4 - 30 * p2 * q + 5 * q2)) / 320
  fr = 1 + delt*delt * (f[2] + delt * (f[3] + delt * (f[4] \
       + delt * (f[5] + delt * (f[6] + delt * (f[7] + delt * f[8]))))))
  gr = delt * (1 + delt*delt * (g[3] + delt * (g[4] + delt * (g[5] \
       + delt * (g[6] + delt * (g[7] + delt * g[8]))))))

  return fr, gr

def so3rv(od, rd, ll):
  """ fit orbit to three observations """

  # initial orbit
  (rr2, vv2) = so2rv(od, rd, ll)

  # refine initial orbit
  rr1 = rr2
  vv1 = vv2
  delt1 = (od[0][0] - od[1][0]) / tu
  delt3 = (od[2][0] - od[1][0]) / tu
  rr2[0] = 0
  rr2[1] = 0
  rr2[2] = 0

  rvx = np.subtract(rr1,rr2) # FIXME: Added this to seed the loop with a defined value
  while (mag(rvx) > 1.0E-12):
    (fr, gr) = f8g(rr1, vv1, delt1)
    f1 = fr
    g1 = gr
    (fr, gr) = f8g(rr1, vv1, delt3)
    f3 = fr
    g3 = gr
    rr2 = rr1
    vv2 = vv1

    mat = np.array([
     [f1 * ll[0][2],  0, -f1 * ll[0][0],  g1 * ll[0][2], 0,
      -g1 * ll[0][0], ll[0][2] * rd[0][0] - ll[0][0] * rd[0][2]],
      [0,             -f1 * ll[0][2],      f1 * ll[0][1], 0,
      -g1 * ll[0][2], g1 * ll[0][1],ll[0][1] * rd[0][2] - ll[0][2] * rd[0][1]],
      [ll[1][2],       0,     -ll[1][0],    0,       0,        0,
       ll[1][2] * rd[1][0] - ll[1][0] * rd[1][2]],
      [0,             -ll[1][2],           ll[1][1],      0,
       0,              0,    ll[1][1] * rd[1][2] - ll[1][2] * rd[1][1]],
      [f3 * ll[2][2],  0,                 -f3 * ll[2][0], g3 * ll[2][2],
       0,   -g3 * ll[2][0],   ll[2][2] * rd[2][0] - ll[2][0] * rd[2][2]],
      [0,             -f3 * ll[2][2],      f3 * ll[2][1], 0,   -g3 * ll[2][2],
       g3 * ll[2][1],   ll[2][1] * rd[2][2] - ll[2][2] * rd[2][1]]
      ])

    # pylint: disable=unbalanced-tuple-unpacking
    (rr1, vv1) = rref(mat)

    # averaging
    rr1 = (rr1 + 3*rr2) / 4                  #  rr1 = (3*rr2 + rr1) / 4
    vv1 = (vv1 + 3*vv2) / 4
    rvx = np.subtract(rr1,rr2)

  vv2 = xke * vv2
  rv2el(rr2, vv2)


def SGN(var):
  if (var<0):
    return -1
  else:
    return 1


def read_obs(iod_lines):
  """ decodes the iod_line data """

  global ssn
  global epoch_datetime
  # FIXME get rid of these globals

  Sites = CosparSite("data/stations.in")
  csi = .0055878713278878
  zet = .0055888307019922
  the = .0048580335354883

  nobs = len(iod_lines) # Number of iod-compliant formatted lines in the input file

  ll    = np.zeros((nobs,3))
  odata = np.zeros((nobs,3))
  rd    = np.zeros((nobs,3))

  i = 0
  for iod_line in iod_lines:
    # Grab the most recent version of these variables for writing eventual TLE
    # FIXME Scott's original code grabs the 2nd or "middle" datetime for epoch
    epoch_datetime = iod_line.DateTime
    ssn = iod_line.ObjectNumber

    ts = load.timescale()
    (year, month, day, hour, minute, second) = iod_line.DateTime.timetuple()[:6]
    t_skyfield = ts.utc(year, month, day, hour, minute, second)
    t1_jd = t_skyfield.tt

    ra = radians(iod_line.RA)
    dc = radians(iod_line.DEC)

    if (iod_line.Epoch == 4): # precess from B1950 to J2000
      a = cos(dc) * sin(ra + csi)
      b = cos(the) * cos(dc) * cos(ra + csi) \
          - sin(the) * sin(dc)
      c = sin(the) * cos(dc) * cos(ra + csi) \
          + cos(the) * sin(dc)
      ra = atan(a / b) # ra - zet
      if (b < 0):
        ra += pi
      ra += zet # right ascension, radians
      ra += 1 / 30000
      dc = asin(c)
      if (abs(dc) > 1.4):
        dc = c / abs(c) * acos(sqrt(a*a + b*b))

    # precession from J2000
    t = (t1_jd - 2451545) / 36525  
    csi = (2306.2181 + .30188 * t + .017998 *t*t) * t * de2ra / 3600
    zet = (2306.2181 + 1.09468 * t + .018203 *t*t) * t * de2ra / 3600
    the = (2004.3109 - .42665 * t - .041833 *t*t) * t * de2ra / 3600
    a = cos(dc) * sin(ra + csi)
    b = cos(the) * cos(dc) * cos(ra + csi) \
        - sin(the) * sin(dc)
    c = sin(the) * cos(dc) * cos(ra + csi) \
        + cos(the) * sin(dc)
    ra = atan(a / b) # ra - zet
    if (b < 0):
      ra += pi
    ra += zet # right ascension, radians
    dc = asin(c)
    if (abs(dc) > 1.4):
      dc = c / abs(c) * acos(sqrt(a*a + b*b))

    # line-of-sight vectors
    ll[i][0] = cos(dc) * cos(ra)
    ll[i][1] = cos(dc) * sin(ra)
    ll[i][2] = sin(dc)

    odata[i][0] = t1_jd # julian date
    odata[i][1] = ra # ra radians (observed)
    odata[i][2] = dc # dc radians (observed)

    (la, lo, hh) = Sites.topos(iod_line.Station)
    observer_location = Topos(latitude_degrees = la, longitude_degrees = lo, elevation_m = hh)

    topocentric = observer_location.at(t_skyfield)
    rd[i] = topocentric.position.km/_XKMPER   # elfind works in units of Earth radii
 
    i+=1
  # end for
  return odata, ll, rd


# MAIN
def main():
  """ Elfind finds and reads observations in the IOD format and computes a set of orbital elements matching the observations.  The observations must be made on a single pass from a single location.

  The default file for input is named "unid.txt".  However, any filename may be used as a command line argument to elfind. i.e.
    elfind input.txt
  would search for observations in the file "input.txt" and append the computed TLE at the end of the same file.  If two command line arguments are found, elfind interprets the first as the input file and the second  as the output file. If elfind is executed with no command argument i.e. 
    elfind 
  the file "unid.txt" is searched for properly formatted observations and the computed TLE is appended to the end of the file "unid.txt".

  A properly formatted observation is recorded as a single line in one of  the four IOD formats using right ascension and declination (1, 2, 3, 7)  Two observations gives the program enough information to find the orbital elements of the circular earth orbit through those observations. Three observations finds all parameters but the atmospheric drag.  A guess is made for the drag element.

  Adapted from C++ version originally by Scott Campbell campbel7@the-i.net
  """  
  log = logging.getLogger(__name__)

  # make it print to the console.
  console = logging.StreamHandler()
  log.addHandler(console)

  global ssn
  global epoch_datetime
  global xincl
  global xnodeo
  global omegao
  global xno     # mean motion
  global xmo     # mean anomaly
  global eo      # eccentricity
  global omegao  # argument perigee

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

  # if verbose:
  #     for arg in vars(args):
  #         log.debug("%s : %s",arg, getattr(args, arg))

  arguments = sys.argv[1:]
  count = len(arguments)

  if (count == 1):
    file_in=sys.argv[1]
    file_out=sys.argv[1]
  elif (count == 2):
    file_in=sys.argv[1]
    file_out=sys.argv[2]
  else: # default input/output files
    file_in = "unid.txt"
    file_out = "unid.txt"

  log.info("file_in: {}  file_out {}".format(file_in,file_out))

  IOD_Records = iod.iod_getrecords(file_in)
  num_file_obs = len(IOD_Records)

  if (num_file_obs < 2):
    log.warning("Not enough observations in file {%s}".format(file_in))
    sys.exit()

  # designate row numbers of obs to be used
  if(num_file_obs == 2):
    selected_lines = IOD_Records
  else:
    while (True):
      print("")
      for i in range(0,num_file_obs):
        print("({:2d}) {:s}".format(i+1, IOD_Records[i].line))

      try:
        userinput = input("Enter Q or the row numbers of 2 or 3 observations: ")
      except:
        sys.exit()

      if (userinput.upper() == "Q"):
        sys.exit()

      selected_lines = []
      if len(userinput)>1:
        for i in userinput.split():
          i = int(i)
          selected_lines.append(IOD_Records[i-1])
      else:
        log.error("Not enough lines specified.")
        continue

      num_selected_lines = len(selected_lines)

      # get line-of-sight vectors
      (odata, ll, rd) = read_obs(selected_lines)

      if (num_selected_lines == 2):
        so2rv(odata, rd, ll)
      else:
        so3rv(odata, rd, ll)
      # TODO: could grab this for the epoch day, but need to convert to TLE time
      # jd = odata[1][0]

      name = "UNID generated by elfind.py"
      (_, tle_line1, tle_line2) = tle_util.make_tle(ssn=ssn, name=name, epoch_datetime=epoch_datetime, xincl=xincl, xnodeo=xnodeo, eo=eo, omegao=omegao, xmo=xmo, xno=xno, deg=False)
      tle_util.append_tle_file(file_out, name, tle_line1, tle_line2)
# end main

if (__name__ == '__main__'):
  main()