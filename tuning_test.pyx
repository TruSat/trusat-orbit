# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

import time
import sys

from cpython.array cimport array, clone
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free
import numpy as numpy
cimport numpy as numpy

cdef int loops

def timefunc(name):
    def timedecorator(f):
        cdef int L, i

        print("Running", name)
        for L in [1, 10, 100, 1000, 10000, 100000, 1000000]:
            start = time.clock()
            f(L)
            end = time.clock()
            print(format((end-start) / loops * 1e6, "2f"), end=" ")
            sys.stdout.flush()

        print("μs")
    return timedecorator

print()
print("INITIALISATIONS")
loops = 100000

@timefunc("cpython.array buffer")
def _(int L):
    cdef int i
    cdef array[double] arr, template = array('d')

    for i in range(loops):
        arr = clone(template, L, False)

    # Prevents dead code elimination
    str(arr[0])

@timefunc("cpython.array memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr
    cdef array template = array('d')

    for i in range(loops):
        arr = clone(template, L, False)

    # Prevents dead code elimination
    str(arr[0])

@timefunc("cpython.array raw C type")
def _(int L):
    cdef int i
    cdef array arr, template = array('d')

    for i in range(loops):
        arr = clone(template, L, False)

    # Prevents dead code elimination
    str(arr[0])

@timefunc("numpy.empty_like memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr
    template = numpy.empty((L,), dtype='double')

    for i in range(loops):
        arr = numpy.empty_like(template)

    # Prevents dead code elimination
    str(arr[0])

@timefunc("malloc")
def _(int L):
    cdef int i
    cdef double* arrptr

    for i in range(loops):
        arrptr = <double*> malloc(sizeof(double) * L)
        free(arrptr)

    # Prevents dead code elimination
    str(arrptr[0])

@timefunc("malloc memoryview")
def _(int L):
    cdef int i
    cdef double* arrptr
    cdef double[::1] arr

    for i in range(loops):
        arrptr = <double*> malloc(sizeof(double) * L)
        arr = <double[:L]>arrptr
        free(arrptr)

    # Prevents dead code elimination
    str(arr[0])

@timefunc("cvarray memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr

    for i in range(loops):
        arr = cvarray((L,),sizeof(double),'d')

    # Prevents dead code elimination
    str(arr[0])



print()
print("ITERATING")
loops = 1000

@timefunc("cpython.array buffer")
def _(int L):
    cdef int i
    cdef array[double] arr = clone(array('d'), L, False)

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arr[i]

    # Prevents dead-code elimination
    str(d)

@timefunc("cpython.array memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr = clone(array('d'), L, False)

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arr[i]

    # Prevents dead-code elimination
    str(d)

@timefunc("cpython.array raw C type")
def _(int L):
    cdef int i
    cdef array arr = clone(array('d'), L, False)

    cdef double d
    for i in range(loops):
        for i in range(L):
#            d = arr[i]
# cpython.array raw C type: Well damn, it's fast. And it's safe. Unfortunately it goes through Python to access its data fields. You can avoid that by using a wonderful trick:
            d = arr.data.as_doubles[i]
# which brings it up to the standard speed while removing safety! This makes this a wonderful replacement for malloc, being basically a pretty reference-counted version!

    # Prevents dead-code elimination
    str(d)

@timefunc("numpy.empty_like memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr = numpy.empty((L,), dtype='double')

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arr[i]

    # Prevents dead-code elimination
    str(d)

@timefunc("malloc")
def _(int L):
    cdef int i
    cdef double* arrptr = <double*> malloc(sizeof(double) * L)

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arrptr[i]

    free(arrptr)

    # Prevents dead-code elimination
    str(d)

@timefunc("malloc memoryview")
def _(int L):
    cdef int i
    cdef double* arrptr = <double*> malloc(sizeof(double) * L)
    cdef double[::1] arr = <double[:L]>arrptr

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arr[i]

    free(arrptr)

    # Prevents dead-code elimination
    str(d)

@timefunc("cvarray memoryview")
def _(int L):
    cdef int i
    cdef double[::1] arr = cvarray((L,),sizeof(double),'d')

    cdef double d
    for i in range(loops):
        for i in range(L):
            d = arr[i]

    # Prevents dead-code elimination
    str(d)

# >>> import tuning_test
# INITIALISATIONS
# Running cpython.array buffer
# 0.095690 0.121790 0.281630 0.195870 0.804040 0.830430 0.777540 μs
# Running cpython.array memoryview
# 0.452850 0.381720 0.647040 0.694140 1.140030 1.212530 1.369450 μs
# Running cpython.array raw C type
# 0.065290 0.053100 0.195280 0.153530 0.736020 0.698500 0.673310 μs
# Running numpy.empty_like memoryview
# 2.529230 2.716970 3.101880 3.095230 4.494870 4.146010 4.239880 μs
# Running malloc
# 0.080140 0.077730 0.180930 0.079910 0.633360 0.607400 26.902850 μs
# Running malloc memoryview
# 1.092770 0.994720 1.123100 1.052200 2.195600 1.998980 29.881840 μs
# Running cvarray memoryview
# 1.058430 0.939900 0.959430 1.318170 1.732960 1.662280 27.739710 μs

# ITERATING
# Running cpython.array buffer
# 0.025000 0.013000 0.026000 0.020000 0.023000 0.021000 0.072000 μs
# Running cpython.array memoryview
# 0.026000 0.017000 0.024000 0.017000 0.022000 0.021000 0.064000 μs
# Running cpython.array raw C type
# 0.028000 0.145000 1.364000 10.497000 101.763000 1001.862000 10567.756000 μs
# Running numpy.empty_like memoryview
# 0.034000 0.011000 0.019000 0.012000 0.013000 0.011000 0.299000 μs
# Running malloc
# 0.006000 0.005000 0.006000 0.007000 0.007000 0.005000 0.005000 μs
# Running malloc memoryview
# 0.013000 0.013000 0.014000 0.009000 0.010000 0.010000 0.195000 μs
# Running cvarray memoryview
# 0.009000 0.011000 0.013000 0.009000 0.009000 0.010000 0.175000 μs