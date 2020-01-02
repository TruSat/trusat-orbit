#!/usr/bin/env python
"""Utility to performance performance profiling on 
satfit and SGP4 routines.
"""

from timeit import Timer
from satfit_accelerated import *
from array import array
import numpy as np

test_vec = array('d',[0.1, 0.2, 0.3])
print(test_vec)
print(unit_vector_raw(test_vec))
print(unit_vector_np(test_vec))
print("norm:    ", Timer(lambda: norm(test_vec)).timeit(number = 100000))
print("mag_raw: ", Timer(lambda: mag_raw(test_vec)).timeit(number = 100000))
print ("unit_vector_raw: ",Timer(lambda: unit_vector_raw(test_vec)).timeit(number = 10000))
print ("unit_vector_np: " ,Timer(lambda: unit_vector_np(test_vec)).timeit(number = 10000))
print("posradang  0: ", Timer(lambda: posradang(0.0)).timeit(number = 10000))
print("posradang  4: ", Timer(lambda: posradang(4.0)).timeit(number = 10000))
print("posradang -1: ", Timer(lambda: posradang(-1.0)).timeit(number = 10000))
print("posradang -1: ", Timer(lambda: posradang(-1.0)).timeit(number = 10000))
print("rtw: ", Timer(lambda: rtw(361.0,-1.0)).timeit(number = 10000))
print("acose 0.5: ", Timer(lambda: acose(0.5)).timeit(number = 10000))
print("acose 1.5: ", Timer(lambda: acose(1.5)).timeit(number = 10000))
print("SGN: ", Timer(lambda: SGN(1.5)).timeit(number = 10000))
print("dot:      ", Timer(lambda: dot(test_vec, test_vec)).timeit(number = 10000))
print("npdot:    ", Timer(lambda: npdot(test_vec, test_vec)).timeit(number = 10000))
print("np.dot:   ", Timer(lambda: np.dot(test_vec, test_vec)).timeit(number = 10000))
print("cross:    ", Timer(lambda: cross(test_vec, test_vec)).timeit(number = 10000))
print("npcross:  ", Timer(lambda: npcross(test_vec, test_vec)).timeit(number = 10000))
print("np.cross: ", Timer(lambda: np.cross(test_vec, test_vec)).timeit(number = 10000))
