import numpy as np

# Constants
B0 = 1.0
m = 1.0
qe = 1.0
c = 1.0

# Functions
def A(x):
    return B0 * (1.0 + x[2]**2)*np.array([-0.5 * x[1], 0.5 * x[0], 0])

def dA(x):
    dA = np.zeros((3, 3))
    dA[0, 1] = 0.5 * (1.0 + x[2]**2)
    dA[1, 0] = -0.5 * (1.0 + x[2]**2)
    dA[2, 0] = -x[1] * x[2]
    dA[2, 1] = x[0] * x[2]
    return dA*B0

def B(x):
    return B0 * np.array([-x[0] * x[2], -x[1] * x[2], 1 + x[2]**2])

def dB(x):
    dB = np.zeros((3, 3))
    dB[0, 0] = -x[2]
    dB[1, 1] = -x[2]
    dB[2, 0] = -x[0]
    dB[2, 1] = -x[1]
    dB[2, 2] = 2*x[2]
    return dB*B0

def Bmod(x):
    return B0 * np.sqrt(x[0]**2*x[2]**2 + x[1]**2*x[2]**2 + (1 + x[2]**2)**2)

def dBmod(x):
    Bmodval = Bmod(x)
    return B0 * np.array([
        x[0]**2 * x[2]**2 / Bmodval,
        x[1]**2 * x[2]**2 / Bmodval,
        x[2] * (x[0]**2 + x[1]**2 + 2.0*x[2]**2 + 2.0) / Bmodval
    ])

def e1(x):
    e1 = np.array([0.0, 1.0+x[2]**2, x[1]*x[2]])
    return e1/np.sqrt(e1[0]**2 + e1[1]**2 + e1[2]**2)

def de1(x):
    de1 = np.zeros((3, 3))
    de1[1, 1] = -((x[1] * x[2]**2 * (1 + x[2]**2)) / (1 + (2 + x[0]**2) * x[2]**2 + x[2]**4)**(3/2))
    de1[2, 1] = (x[2] * (1 + x[2]**2)**2) / (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)
    de1[1, 2] = (x[1]**2 * x[2] * (-1 + x[2]**2)) / (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)
    de1[2, 2] = (x[1] - x[1] * x[2]**4) / (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)
    return de1

def e2(x):
    e2 = np.array([
        1.0 + x[2]**2 * (2 + x[1]**2 + x[2]**2),
        -x[0] * x[1] * x[2]**2,
        x[0] * x[2] * (1 + x[2]**2)
    ])
    return e2/np.sqrt(e2[0]**2 + e2[1]**2 + e2[2]**2)

def de2(x):
    de2 = np.zeros((3, 3))
    de2[0, 0] = (x[0] * (x[2] + (2 + x[1]**2) * x[2]**3 + x[2]**5)**2) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2))
    de2[1, 0] = (x[1] * (x[2] + (2 + x[1]**2) * x[2]**3 + x[2]**5)**2) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2))
    de2[2, 0] = -((x[2] * (1 + x[2]**2) * (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)**2) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)))
    de2[0, 1] = -((x[0]**2 * x[1] * x[2]**4 * (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)))
    de2[1, 1] = (x[0] * x[2]**2 * (1 + x[2]**2 * (4 - (-6 + x[1]**4) * x[2]**2 + 4 * x[2]**4 + x[2]**6 + x[0]**2 * (1 + x[2]**2)**2))) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2))
    de2[2, 1] = (x[0] * x[1] * x[2]**3 * (1 + x[2]**2) * (2 + x[2]**2 * (x[0]**2 + 2 * (2 + x[1]**2 + x[2]**2)))) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2))
    de2[0, 2] = -((x[0]**2 * x[2] * (-1 + x[2]**4) * (1 + (2 + x[1]**2) * x[2]**2 + x[2]**4)) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)))
    de2[1, 2] = -((x[0] * x[1] * x[2] * (-1 + x[2]**4) * (2 + x[2]**2 * (x[0]**2 + 2 * (2 + x[1]**2 + x[2]**2)))) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2)))
    de2[2, 2] = (x[0] * (-1 + x[2]**2) * (1 + 4 * x[2]**2 - (-6 + x[0]**2 * x[1]**2 + x[1]**4) * x[2]**4 + 4 * x[2]**6 + x[2]**8)) / ((1 + (2 + x[1]**2) * x[2]**2 + x[2]**4) * (1 + (2 + x[0]**2 + x[1]**2) * x[2]**2 + x[2]**4)**(3/2))
    return de2


# Evaluate the expression with the given values of x[0], x[1], and x[2]
