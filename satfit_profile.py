#!/usr/bin/env python
"""Utility to performance performance profiling on 
satfit and SGP4 routines.
"""

from timeit import Timer
from satfit_caccelerated import *
from array import array
import numpy as np

test_vec = array('d',[0.1, 0.2, 0.3])
test_vec2 = array('d',[0.3, 0.2, 0.1])
rtn_vec = array('d',[0, 0, 0])

# np_arr_dbl_template = numpy.empty((3,), dtype='double')
# test_vec_np = numpy.empty_like(np_arr_dbl_template)
test_vec_np = np.array([0.1,0.2,0.3],dtype=np.double)
r = 0.0

print(test_vec)
print(unit_vector_raw(test_vec))
print(unit_vector_np(test_vec))
print("norm:          ", Timer(lambda: norm(test_vec)).timeit(number = 100000))
print("normc:         ", Timer(lambda: normc(test_vec)).timeit(number = 100000))
print("mag_raw:       ", Timer(lambda: mag_raw(test_vec)).timeit(number = 100000))
print("mag_rawnp:     ", Timer(lambda: mag_raw(test_vec_np)).timeit(number = 100000))
print()
print ("unit_vector_raw:       ",Timer(lambda: unit_vector_raw(test_vec)).timeit(number = 10000))
print ("unit_vector_rawnp:     ",Timer(lambda: unit_vector_raw(test_vec_np)).timeit(number = 10000))
print ("unit_vector_raw_ref:   ",Timer(lambda: unit_vector_raw_ref(test_vec, rtn_vec)).timeit(number = 10000))
print ("unit_vector_raw_def:   ",Timer(lambda: unit_vector_raw_def(test_vec)).timeit(number = 10000))
# print("posradang  0: ", Timer(lambda: posradang(0.0)).timeit(number = 10000))
# print("posradang  4: ", Timer(lambda: posradang(4.0)).timeit(number = 10000))
# print("posradang -1: ", Timer(lambda: posradang(-1.0)).timeit(number = 10000))
# print("posradang -1: ", Timer(lambda: posradang(-1.0)).timeit(number = 10000))
# print("rtw: ", Timer(lambda: rtw(361.0,-1.0)).timeit(number = 10000))
# print("acose 0.5: ", Timer(lambda: acose(0.5)).timeit(number = 10000))
# print("acose 1.5: ", Timer(lambda: acose(1.5)).timeit(number = 10000))
# print("SGN: ", Timer(lambda: SGN(1.5)).timeit(number = 10000))
print()
print("dot:      ", Timer(lambda: dot(test_vec, test_vec)).timeit(number = 10000))
print("dotnp:    ", Timer(lambda: dot(test_vec_np, test_vec_np)).timeit(number = 10000))
print("npdot:    ", Timer(lambda: npdot(test_vec, test_vec)).timeit(number = 10000))
print("np.dot:   ", Timer(lambda: np.dot(test_vec, test_vec)).timeit(number = 10000))
print()
print("cross_rtn:", Timer(lambda: cross_rtn(test_vec, test_vec)).timeit(number = 10000))
print("cross_np: ", Timer(lambda: cross_rtn(test_vec_np, test_vec_np)).timeit(number = 10000))
print("cross_ref:", Timer(lambda: cross_ref(test_vec, test_vec, rtn_vec)).timeit(number = 10000))
print("npcross:  ", Timer(lambda: npcross(test_vec, test_vec)).timeit(number = 10000))
print("np.cross: ", Timer(lambda: np.cross(test_vec, test_vec)).timeit(number = 10000))
print()
print("vmadd_rtn: ", Timer(lambda: vmadd_rtn(test_vec, test_vec2, 1)).timeit(number = 10000))
print("vmadd_np:  ", Timer(lambda: vmadd_rtn(test_vec_np, test_vec_np, 1)).timeit(number = 10000))
print("vmadd_ref: ", Timer(lambda: vmadd_ref(test_vec, test_vec2, rtn_vec, 1)).timeit(number = 10000))
print()
print("smult_rtn: ", Timer(lambda: smult_rtn(2, test_vec)).timeit(number = 10000))
print("smult_np:  ", Timer(lambda: smult_rtn(2, test_vec_np)).timeit(number = 10000))
print("smult_ref: ", Timer(lambda: smult_ref(2, test_vec, rtn_vec)).timeit(number = 10000))
print("smult_py:  ", Timer(lambda: smult_py(2, test_vec)).timeit(number = 10000))
print()
print("proj_rtn: ", Timer(lambda: proj_rtn(test_vec, test_vec)).timeit(number = 10000))
print("proj_np : ", Timer(lambda: proj_rtn(test_vec_np, test_vec_np)).timeit(number = 10000))
print("proj_ref: ", Timer(lambda: proj_ref(test_vec, test_vec, rtn_vec)).timeit(number = 10000))
print("proj_cdef:", Timer(lambda: proj_cdef(test_vec, test_vec)).timeit(number = 10000))
