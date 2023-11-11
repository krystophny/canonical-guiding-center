#%%
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, simplify, diff
from sympy.vector import CoordSys3D, curl, gradient

# Define symbols
z, p1, p2 = symbols('z p1 p2')

# Define the coordinate system
N = CoordSys3D('N')

# Define B0 and A as functions
Bmod = 1.0 + N.z**2
A = Bmod*(-0.5 * N.y * N.i + 0.5 * N.x * N.j)

# Calculate B as the curl of A and simplify
B = curl(A)
h = (B/Bmod)
# %%
