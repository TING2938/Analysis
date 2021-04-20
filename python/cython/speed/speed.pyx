cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

cdef _speed(unsigned long long num):
    cdef unsigned long long i = num - 1
    while i > 1:
        if num % i == 0:
            return i
        i -= 1
    return num

def speed(num):
    return _speed(num)