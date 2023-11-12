#%%
import numpy as np
from sympy import symbols, simplify, diff
from sympy.vector import CoordSys3D, divergence, curl, gradient
from sympy.plotting import plot3d
from spb import *

# Define symbols
z, p1, p2 = symbols('z p1 p2')

# Define the coordinate system
N = CoordSys3D('N')
i, j, k = N.base_vectors()
x, y, z = N.base_scalars()

# Magnetic bottle
Bscale = 1.0 + z**2

# Uniform field
#Bscale = 1.0

#%%
A = Bscale*(-0.5 * y * i + 0.5 * x * j)
A
#%%
B = curl(A)
B
#%%
Bmod = B.magnitude()
Bmod
#%%
h = B.normalize()
h
#%% Basis vector field
e1 = h.cross(A.normalize()).simplify()
e1

#%%
curl(e1).simplify()
#%%
divergence(e1).simplify()

# %%
