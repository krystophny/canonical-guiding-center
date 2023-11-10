import numpy as np

# Constants
B0 = 1.0
m = 1.0
qe = 1.0
c = 1.0

# Functions
def Bmod(x):
    return B0 * (1.0 + x[2]**2)

def A(x):
    return Bmod(x)*np.array([-0.5 * x[1], 0.5 * x[0], 0])

def dA(x):
    dA = np.zeros((3, 3))
    dA[0, 1] = 0.5 * (1.0 + x[2]**2)
    dA[1, 0] = -0.5 * (1.0 + x[2]**2)
    dA[2, 0] = -x[1] * x[2]
    dA[2, 1] = x[0] * x[2]
    return dA*B0

def B(x):
    return B0 * np.array([-x[0] * x[2], -x[1] * x[2], 1 + x[2]**2])

def e1(x):
    e1 = np.array([0.0, 1.0+x[2]**2, x[1]*x[2]])
    return e1/np.sqrt(e1[0]**2 + e1[1]**2 + e1[2]**2)

def e2(x):
    e2 = np.array([
        1.0 + x[2]**2 * (2 + x[1]**2 + x[2]**2),
        -x[0] * x[1] * x[2]**2,
        x[0] * x[2] * (1 + x[2]**2)
    ])
    return e2/np.sqrt(e2[0]**2 + e2[1]**2 + e2[2]**2)
