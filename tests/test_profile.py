#!/usr/bin/env python
"""Utility to performance performance profiling on 
satfit and SGP4 routines.
"""
import sys
from timeit import Timer

from trusat.caccelerated import *
try:
    import trusat.profile as tsp
    profile=True
    print("Including performance profiling code.")
except ImportError: 
    profile=False
    print("Cython profile module not available, skipping some tests.")

from array import array
import numpy as np

from sgp4.propagation import sgp4, sgp4init
from sgp4.api import Satrec, SatrecArray, SGP4_ERRORS, WGS72
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

from trusat.satfit import Date

line0 = 'SL-16 R/B'
line1 = '1 22285U 92093B   19314.09990558 -.00000000  00000-0  24310-4 0  9990'
line2 = '2 22285  71.0198 174.7928 0006190 244.6688 115.3794 14.15033886386443'

satold = twoline2rv(line1,line2,wgs72)
satnew = Satrec.twoline2rv(line1,line2)

test_tuple           = (0.1, 0.2, 0.3)
test_list            = [0.1, 0.2, 0.3]
test_vec   = array('d',[0.1, 0.2, 0.3])
test_vec2  = array('d',[0.3, 0.2, 0.1])
rtn_vec    = array('d',[0, 0, 0])

test_vec_np = np.array([0.1,0.2,0.3],dtype=np.double)
rtn_vec_np  = np.zeros(3,dtype=np.double)
r = 0.0

b = 0.000012345
bmax = 1.1*b
bmin = 0.9*b
bstep = (bmax-bmin)/20.0

satnum   = 5
bstar    = 2.8098e-05
ndot     = 6.96919666594958e-13
nddot    = 0.0
ecco     = 0.1859667
argpo    = 5.790416027488515
inclo    = 0.5980929187319208
mo       = 0.3373093125574321
no_kozai = 0.04722944544077857
nodeo    = 6.08638547138321

sat2 = Satrec()
sat2_jdsatepoch  = 2451723.28495062

sat2.sgp4init(WGS72, 'i', satnum, sat2_jdsatepoch, bstar, ndot, nddot,
              ecco, argpo, inclo, mo, no_kozai, nodeo)


# Vectorized Python-SGP4 array initialization
jd = np.linspace(2451723, 2451723+900, 100)
fr = np.zeros((jd.shape[0],))
rv_template = np.zeros((1, len(jd), 3), dtype='double')

print()
np_zeros_template = np.zeros((1,10,3), dtype='double')
np_empty_template = np.empty((1,10,3), dtype='double')

def tac(baseline, time):
    comp = time / baseline
    return time, comp

num = 10000
baseline = Timer(lambda: np.zeros((3,), dtype='double')).timeit(number = num)
print("np.zeros((3,)                {:.6f} {:.2f}x (baseline)".format( *tac(baseline,Timer(lambda: np.zeros((3,), dtype='double')).timeit(number = num)) ))
print("np.empty((3,)                {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: np.empty((3,), dtype='double')).timeit(number = num)) ))
print("np.zeros((1,10,3)            {:.6f} {:.2f}x (baseline)".format( *tac(baseline,Timer(lambda: np.zeros((1,10,3), dtype='double')).timeit(number = num)) ))
print("np.empty((1,10,3)            {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: np.empty((1,10,3), dtype='double')).timeit(number = num)) ))
print()
print("np.empty_like(zeros)         {:.6f} {:.2f}x (baseline)".format( *tac(baseline,Timer(lambda: np.empty_like(np_zeros_template)).timeit(number = num)) ))
print("np.empty_like(empty)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: np.empty_like(np_empty_template)).timeit(number = num)) ))
if (profile):
    print("np_empty_like(zeros)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_empty_like(np_zeros_template)).timeit(number = num)) ))
    print("np_empty_like(empty)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_empty_like(np_empty_template)).timeit(number = num)) ))
    print("np_empty_like_c(zeros)       {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_empty_like_c(np_zeros_template)).timeit(number = num)) ))
    print("np_empty_like_c(empty)       {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_empty_like_c(np_empty_template)).timeit(number = num)) ))
print("np.zeros_like(empty)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: np.zeros_like(np_empty_template)).timeit(number = num)) ))
print("np.zeros_like(zeros)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: np.zeros_like(np_zeros_template)).timeit(number = num)) ))
if (profile):
    print("np_zeros_like(empty)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_zeros_like(np_empty_template)).timeit(number = num)) ))
    print("np_zeros_like(zeros)         {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_zeros_like(np_zeros_template)).timeit(number = num)) ))
    print("np_zeros_like_c(empty)       {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_zeros_like_c(np_empty_template)).timeit(number = num)) ))
    print("np_zeros_like_c(zeros)       {:.6f} {:.2f}x".format( *tac(baseline,Timer(lambda: tsp.np_zeros_like_c(np_zeros_template)).timeit(number = num)) ))
print()

if (profile):
    baseline = Timer(lambda: tsp.vectorized_init_variable(jd)).timeit(number = num)
    print("vectorized_init_variable()   {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vectorized_init_variable(jd)).timeit(number = num)) ))
    print("vectorized_init_fixed()      {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vectorized_init_fixed(jd)).timeit(number = num)) ))
    print("vectorized_init_npview()     {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vectorized_init_npview(jd, rv_template)).timeit(number = num)) ))
    print()

if (profile):
    baseline = Timer(lambda: tsp.tuple_to_array(test_tuple)).timeit(number = num)
    print("tuple_to_array()             {:.6f} {:.2f}x (baseline, inline loop)".format( *tac(baseline, Timer(lambda: tsp.tuple_to_array(test_tuple)).timeit(number = num)) ))
    print("np.asarray()                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.asarray(test_tuple)).timeit(number = num)) ))
    print("np_asarray()                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.np_asarray(test_tuple)).timeit(number = num)) ))
    print()

baseline = Timer(lambda: float_step(bmin,bmax,bstep)).timeit(number = num)
print("float_step():                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: float_step(bmin,bmax,bstep)).timeit(number = num)) ))
print("np.arrange():                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.arange(bmin,bmax,bstep)).timeit(number = num)) ))
print()

print("posradang():                 {:.6f}".format( Timer(lambda: posradang(1.87)).timeit(number = num)) )
if (profile):
    print("pythonmodulo():              {:.6f}".format( Timer(lambda: tsp.pythonmodulo(1.87)).timeit(number = num)) )
    print(f"340:{tsp.pythonmodulo(340)}  380:{tsp.pythonmodulo(380)}  -20:{tsp.pythonmodulo(-20)}")
print()

print("Only one part...")
print("param % 1:                   {:.6f}".format( Timer(lambda: sat2.jdsatepoch % 1).timeit(number = num)) )
if (profile):
    print("floor(f):                    {:.6f}".format( Timer(lambda: tsp.floor_floor(sat2.jdsatepoch)).timeit(number = num)) )
    print("int(f):                      {:.6f}".format( Timer(lambda: tsp.floor_mod(sat2.jdsatepoch)).timeit(number = num)) )
print()

if(profile):
    print("Both int and frac parts")
    baseline = Timer(lambda: tsp.modf_sub(sat2.jdsatepoch)).timeit(number = num)
    print("modf_sub(param,f):           {:.6f} {:.2f}x (baseline - use modf() directly)".format( *tac(baseline, Timer(lambda: tsp.modf_sub(sat2.jdsatepoch)).timeit(number = num)) ))
    print("divmod_sub() :               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.divmod_sub(sat2.jdsatepoch)).timeit(number = num)) ))
    print("divmod_raw() :               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.divmod_raw(sat2.jdsatepoch)).timeit(number = num)) ))
    print("divmod(param):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: divmod(sat2.jdsatepoch, 1)).timeit(number = num)) ))
    print("modf_divmod():               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.modf_divmod(sat2.jdsatepoch)).timeit(number = num)) ))
    print()

number = 100000
baseline = Timer(lambda: norm(test_vec)).timeit(number = num)
print("norm(pv):                    {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: norm(test_vec)).timeit(number = num)) ))
print("norm(np):                    {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: norm(test_vec_np)).timeit(number = num)) ))
if (profile):
    print("mag_raw(pv):                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.mag_raw(test_vec)).timeit(number = num)) ))
    print("mag_raw(np):                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.mag_raw(test_vec_np)).timeit(number = num)) ))
    print("norm_py(pv):                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.norm_py(test_vec)).timeit(number = num)) ))
    print("np.linalg.norm(np):          {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.mag(test_vec_np)).timeit(number = num)) ))
    print("np.linalg.norm(pv):          {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.mag(test_vec)).timeit(number = num)) ))
print()

number = 1000
baseline = Timer(lambda: unit_vector(test_vec)).timeit(number = num)
print("unit_vector(pv):             {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: unit_vector(test_vec)).timeit(number = num)) ))
print("unit_vector(np):             {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: unit_vector(test_vec_np)).timeit(number = num)) ))
if (profile):
    print("unit_vector_def(pv):         {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.unit_vector_def(test_vec)).timeit(number = num)) ))
    print("unit_vector_def(np):         {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.unit_vector_def(test_vec_np)).timeit(number = num)) ))
    print("unit_vector_ref(pv):         {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.unit_vector_ref(test_vec, rtn_vec)).timeit(number = num)) ))
    print("unit_vector_np(np):          {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.unit_vector_np(test_vec_np)).timeit(number = num)) ))
print()

baseline = Timer(lambda: dot(test_vec, test_vec)).timeit(number=num)
print("dot(pv):                     {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: dot(test_vec, test_vec)).timeit(number = num)) ))
print("dot(np):                     {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: dot(test_vec_np, test_vec_np)).timeit(number = num)) ))
print("np.dot(np):                  {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.dot(test_vec_np, test_vec_np)).timeit(number = num)) ))
print("np.dot(pv):                  {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.dot(test_vec, test_vec)).timeit(number = num)) ))
if (profile):
    print("npdot(pv):                   {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.npdot(test_vec, test_vec)).timeit(number = num)) ))
    print("dot_prange                   {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.dot_prange(test_vec, test_vec)).timeit(number = num)) ))
print()

if (profile):
    baseline = Timer(lambda: tsp.cross_rtn(test_vec, test_vec)).timeit(number=num)
    print("cross_rtn(pv):               {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: tsp.cross_rtn(test_vec, test_vec)).timeit(number = num)) ))
    print("cross_rtn(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.cross_rtn(test_vec_np, test_vec_np)).timeit(number = num)) ))
    print("cross_ref(pv):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.cross_ref(test_vec, test_vec, rtn_vec)).timeit(number = num)) ))
    print("cross_ref(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.cross_ref(test_vec_np, test_vec_np, rtn_vec_np)).timeit(number = num)) ))
    print("np.cross(pv) :               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.cross(test_vec, test_vec)).timeit(number = num)) ))
    print("npcross(pv):                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.npcross(test_vec, test_vec)).timeit(number = num)) ))
    print()

baseline = Timer(lambda: vmadd_ref(test_vec, test_vec2, rtn_vec, 1)).timeit(number = num)
print("vmadd_ref(pv):               {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: vmadd_ref(test_vec, test_vec2, rtn_vec, 1)).timeit(number = num)) ))
print("vmadd_ref(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: vmadd_ref(test_vec_np, test_vec_np, rtn_vec_np, 1)).timeit(number = num)) ))
if (profile):
    print("vmadd_rtn(pv):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vmadd_rtn(test_vec, test_vec2, 1)).timeit(number = num)) ))
    print("vmadd_rtn(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vmadd_rtn(test_vec_np, test_vec_np, 1)).timeit(number = num)) ))
    print("vmadd_ref_prange(pv):        {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.vmadd_ref_prange(test_vec, test_vec2, rtn_vec, 1)).timeit(number = num)) ))
print()

baseline = Timer(lambda: smult_ref(2, test_vec, rtn_vec)).timeit(number = num)
print("smult_ref(pv):               {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: smult_ref(2, test_vec, rtn_vec)).timeit(number = num)) ))
print("smult_ref(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: smult_ref(2, test_vec_np, rtn_vec_np)).timeit(number = num)) ))
if (profile):
    print("smult_rtn(pv):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.smult_rtn(2, test_vec)).timeit(number = num)) ))
    print("smult_rtn(np):               {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.smult_rtn(2, test_vec_np)).timeit(number = num)) ))
    print("smult_py(pv):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.smult_py(2, test_vec)).timeit(number = num)) ))
    print("smult_ref_prange(pv):        {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.smult_ref_prange(2, test_vec, rtn_vec)).timeit(number = num)) ))
print()

if (profile):
    print("Not currently using proj")
    baseline = Timer(lambda: tsp.proj_ref(test_vec, test_vec, rtn_vec)).timeit(number = num)
    print("proj_ref(pv):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.proj_ref(test_vec, test_vec, rtn_vec)).timeit(number = num)) ))
    print("proj_ref(np):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.proj_ref(test_vec_np, test_vec_np, rtn_vec_np)).timeit(number = num)) ))
    print("proj_rtn(pv):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.proj_rtn(test_vec, test_vec)).timeit(number = num)) ))
    print("proj_rtn(np):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.proj_rtn(test_vec_np, test_vec_np)).timeit(number = num)) ))
    print()

print("Not currently using np.asarray")
baseline = Timer(lambda: np.asarray(test_vec_np)).timeit(number = num)
print("np.asarray(np):              {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.asarray(test_vec_np)).timeit(number = num)) ))
print("np.asarray(pv):              {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.asarray(test_vec)).timeit(number = num)) ))
print("np.asarray(tuple):           {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.asarray(test_tuple)).timeit(number = num)) ))
print("np.asarray(list):            {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: np.asarray(test_list)).timeit(number = num)) ))
print()

baseline = Timer(lambda: delta_t(satnew,3600)).timeit(number = num)
print("delta_t:                     {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: delta_t(satnew,3600)).timeit(number = num)) ))
if (profile):
    print("delta_t_old:                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.delta_t_old(satold,3600)).timeit(number = num)) ))
print()

rd = np.array([[0.59544982, 0.14061455, 0.78833855], [0.60865287, -0.06224648, 0.78833855], [0.60868465, -0.06193494, 0.78833855], [0.60872963, -0.06149125, 0.78833855], [0.55200408, -0.26386443, 0.78833855], [0.55219627, -0.26346199, 0.78833855], [0.55238828, -0.26305917, 0.78833855], [0.52008562, -0.32224818, 0.78833855], [0.52024938, -0.32198372, 0.78833855], [0.52036712, -0.32179341, 0.78833855]],dtype=np.double)
ll = np.array([[0.00456624, 0.32275503, 0.94647152], [0.01598115, 0.56115783, 0.82755452], [-0.01368066, 0.56751477, 0.82324955], [-0.05363986, 0.57487537, 0.81648091], [-0.5628066, -0.19572568, 0.80308168], [-0.5620066, -0.15578733, 0.8123293], [-0.56062117, -0.11695225, 0.81977197], [-0.17582553, -0.36071339, 0.91595373], [-0.13835716, -0.37005179, 0.91865063], [-0.10908354, -0.37652389, 0.91996225]]
,dtype=np.double)
odata = np.array([[2458715.61, 1.55664956, 1.24212334, 4172.0, 408299.0], [2458722.54, 1.54232514, 0.974737447, 4172.0, 411181.0], [2458722.54, 1.59489793, 0.967111755, 4172.0, 411182.0], [2458722.54, 1.66383389, 0.955289479, 4172.0, 411183.0], [2458748.41, 3.47627697, 0.932449088, 4172.0, 416464.0], [2458748.41, 3.41200155, 0.948135077, 4172.0, 416465.0], [2458748.41, 3.34725501, 0.961012729, 4172.0, 416466.0], [2458798.26, 4.25884112, 1.15787821, 4172.0, 421414.0], [2458798.26, 4.35459538, 1.16465128, 4172.0, 421415.0], [2458798.26, 4.43039711, 1.16798418, 4172.0, 421416.0]]
,dtype=np.double)

if(profile):
    baseline = Timer(lambda: tsp.find_rms_inline_scalar(satnew, rd, ll, odata)).timeit(number=num)
    print("find_rms_inline_scalar:      {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.find_rms_inline_scalar(satnew, rd, ll, odata)).timeit(number = num)) ))
    print("find_rms:                    {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: find_rms(satnew, rd, ll, odata)).timeit(number = num)) ))
    print("find_rms_old:                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.find_rms_old(satold, rd, ll, odata)).timeit(number = num)) ))
    print("find_rmspy:                  {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.find_rmspy(satold, rd, ll, odata)).timeit(number = num)) ))
else:
    print("find_rms:                    {:.6f}".format( Timer(lambda: find_rms(satnew, rd, ll, odata)).timeit(number = num)) )

print()
find_rms_inline_baseline = baseline

tsince = (60 + satnew.jdsatepoch) * 1440.0 # time since epoch in minutes
(jd, fr) = divmod(satnew.jdsatepoch,1)
satnew.sgp4(jd, fr)

satrec3 = Satrec()
satrec3_jdsatepoch  = 2451723.28495062
satrec3.sgp4init(WGS72, 'i', satnum, satrec3_jdsatepoch-2433281.5, bstar, 
                ndot, nddot, ecco, argpo, inclo, mo, no_kozai, nodeo)

print("Single element SGP4")
baseline = Timer(lambda: satrec3.sgp4(jd,fr)).timeit(number = num)
print("sgp4_cpp(1):                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: satrec3.sgp4(jd,fr)).timeit(number = num)) ))
print("sgp4_py:                     {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: sgp4(satold, tsince)).timeit(number = num)) ))
print()

print("Vector element SGP4 vs scalar loop of same #")
sats = []
sats.append(satnew)
a = SatrecArray(sats)

for j in range(len(odata)):
    r = np.ndarray((len(a), j+1, 3))
    v = np.ndarray((len(a), j+1, 3))
    err = np.ndarray((len(a), j+1), 'uint8')
    jd = np.zeros((j+1,))
    fr = np.ndarray((j+1,))

    for i in range(j):
        jd[i] = (odata[i][0]) // 1.0
        fr[i] = (odata[i][0]) % 1.0
    _sgp4_cpp_arr_result = Timer(lambda: satrec3._sgp4(jd, fr, err, r, v)).timeit(number = num)
    print("_sgp4_cpp({:02d}):               {:.6f} {:.2f}x".format(j+1, _sgp4_cpp_arr_result, _sgp4_cpp_arr_result/(baseline*(j+1)) ))
print()

sats = []
sats.append(satnew)
a = SatrecArray(sats)

rr  = np.ndarray((len(a), len(odata), 3))
vv  = np.ndarray((len(a), len(odata), 3))
err = np.ndarray((len(a), len(odata)), 'uint8')
jd  = np.zeros((len(odata),))
fr  = np.ndarray((len(odata),))

for i in range(len(odata)):
    (jd[i], fr[i]) = divmod(odata[i][0],1)
print("_sgp4_cpp_arr(odata {:d}):     {:.6f}".format(len(odata), Timer(lambda: satrec3._sgp4(jd, fr, err, rr, vv)).timeit(number = num) ))

print()
print("Array of {} times, relative to find_rms_inline_scalar()".format(len(odata)))
print("find_rms_inline_arr:         {:.6f}  {:.2f}x".format( *tac(find_rms_inline_baseline,Timer(lambda: find_rms_inline(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))
if (profile):
    print("find_rms_inline_arr_prange:  {:.6f}  {:.2f}x".format( *tac(find_rms_inline_baseline,Timer(lambda: tsp.find_rms_inline_prange(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))

print()
large_arr_size = 100
# Larger, vectorized Python-SGP4 array initialization
jd  = np.linspace(odata[0][0], odata[0][0]+900, large_arr_size)
fr  = np.zeros((jd.shape[0],))
err = np.ndarray((len(a), len(jd)), 'uint8')
rr  = np.ndarray((len(a), len(jd), 3))
vv  = np.ndarray((len(a), len(jd), 3))
print("Array of {} times, relative to find_rms_inline_scalar()".format(jd.shape[0]))
find_rms_inline_baseline_large = Timer(lambda: find_rms_inline(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)
print("find_rms_inline_arr:         {:.6f}  {:.2f}x".format( *tac(find_rms_inline_baseline,Timer(lambda: find_rms_inline(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))
if (profile):
    print("find_rms_inline_arr_prange:  {:.6f}  {:.2f}x".format( *tac(find_rms_inline_baseline, Timer(lambda: tsp.find_rms_inline_prange(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))
print()

if (profile):
    print("Exploring overhead of passing varibles, _sgp4 vs find_rms contribution...")
    print("Relative to find_rms_inline_arr({:d})".format(large_arr_size))
    print("func_satx:                   {:.6f} {:.2f}x".format( *tac(find_rms_inline_baseline_large, Timer(lambda: tsp.func_satx(satnew)).timeit(number = num)) ))
    print("func_arrs:                   {:.6f} {:.2f}x".format( *tac(find_rms_inline_baseline_large, Timer(lambda: tsp.func_arrs(rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))
    print("func_sgp4:                   {:.6f} {:.2f}x".format( *tac(find_rms_inline_baseline_large, Timer(lambda: tsp.func_sgp4(satnew, rd, ll, jd, fr, rr, vv, err)).timeit(number = num)) ))
    print()

baseline = Timer(lambda: len(jd)).timeit(number = num)
print("len(jd):                     {:.6f} {:.2f}x (baseline)".format( *tac(baseline, Timer(lambda: len(jd)).timeit(number = num)) ))
print("jd.shape[0]:                 {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: jd.shape[0]).timeit(number = num)) ))
if (profile):
    print("func_len(jd):                {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.func_len(jd)).timeit(number = num)) ))
    print("func_np_shape(jd):           {:.6f} {:.2f}x".format( *tac(baseline, Timer(lambda: tsp.func_np_shape(jd)).timeit(number = num)) ))
print()

# Old Python
# satold = twoline2rv(line1,line2,wgs72)
# move_epoch_to_jd_old(satold,t1.jd)

# uu = longitude(satold)
# sum = find_rmspy(satold, rd, ll, odata)
# print(f"rms_old: {sum}")

# step(satold, rd, ll, odata, sum, uu, "L")

# rd = np.array([[0.59544982, 0.14061455, 0.78833855], [0.60865287, -0.06224648, 0.78833855], [0.60868465, -0.06193494, 0.78833855], [0.60872963, -0.06149125, 0.78833855], [0.55200408, -0.26386443, 0.78833855], [0.55219627, -0.26346199, 0.78833855], [0.55238828, -0.26305917, 0.78833855], [0.52008562, -0.32224818, 0.78833855], [0.52024938, -0.32198372, 0.78833855], [0.52036712, -0.32179341, 0.78833855]],dtype=np.double)
# ll = np.array([[0.00456624, 0.32275503, 0.94647152], [0.01598115, 0.56115783, 0.82755452], [-0.01368066, 0.56751477, 0.82324955], [-0.05363986, 0.57487537, 0.81648091], [-0.5628066, -0.19572568, 0.80308168], [-0.5620066, -0.15578733, 0.8123293], [-0.56062117, -0.11695225, 0.81977197], [-0.17582553, -0.36071339, 0.91595373], [-0.13835716, -0.37005179, 0.91865063], [-0.10908354, -0.37652389, 0.91996225]]
# ,dtype=np.double)
# odata = np.array([[2458715.61, 1.55664956, 1.24212334, 4172.0, 408299.0], [2458722.54, 1.54232514, 0.974737447, 4172.0, 411181.0], [2458722.54, 1.59489793, 0.967111755, 4172.0, 411182.0], [2458722.54, 1.66383389, 0.955289479, 4172.0, 411183.0], [2458748.41, 3.47627697, 0.932449088, 4172.0, 416464.0], [2458748.41, 3.41200155, 0.948135077, 4172.0, 416465.0], [2458748.41, 3.34725501, 0.961012729, 4172.0, 416466.0], [2458798.26, 4.25884112, 1.15787821, 4172.0, 421414.0], [2458798.26, 4.35459538, 1.16465128, 4172.0, 421415.0], [2458798.26, 4.43039711, 1.16798418, 4172.0, 421416.0]]
# ,dtype=np.double)

# # # New CPP hotness
# satnew = Satrec.twoline2rv(line1,line2)
# # move_epoch_to_jd(satnew,t1.jd)
# uu = longitude(satnew)
# # print(f"uu: {uu}")
# sum = find_rms(satnew, rd, ll, odata)
# # print(f"rms: {sum}")
# step(satnew, rd, ll, odata, sum, uu, "L")