import numpy as np

# Constants
B0 = 1.0
m = 1.0
qe = 1.0
c = 1.0

# Functions
def Bmod(x):
    return B0

def dBmod(x):
    return 0.0

def A(x):
    return B0*np.array([0.0, x[0], 0.0])

def dA(x):
    dA = np.zeros((3, 3))
    dA[0, 1] = B0
    return dA

def B(x):
    return np.array([0.0, 0.0, B0])

def dB(x):
    return np.zeros((3,3))

def e1(x):
    return np.array([0.0, 1.0, 0.0])

def e2(x):
    return np.array([1.0, 0.0, 0.0])

def de1(x):
    return np.zeros((3,3))

def de2(x):
    return np.zeros((3,3))
