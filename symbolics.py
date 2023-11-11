#%%
import numpy as np
from sympy import symbols, simplify, diff
from sympy.vector import CoordSys3D, curl, gradient
from sympy.plotting import plot3d
from spb import *

# Define symbols
z, p1, p2 = symbols('z p1 p2')

# Define the coordinate system
N = CoordSys3D('N')
i, j, k = N.base_vectors()
x, y, z = N.base_scalars()

# Define B0 and A as functions
Bmod = 1.0 + z**2
A = Bmod*(-0.5 * y * i + 0.5 * x * j)

# Calculate B as the curl of A and simplify
B = curl(A)
Bmod = B.magnitude()
h = B.normalize()
e1 = h.cross(A.normalize()).simplify()
# %%

# Define the range for the variables (e.g., -1 to 1)
x_range = (-1, 1)
y_range = (-1, 1)
z_range = (-1, 1)

# Create a 3D vector field plot
plot = vector_field_3d(A, (x, *x_range), (y, *y_range), (z, *z_range), title='3D Vector Field', xlabel='X', ylabel='Y', zlabel='Z', show=False)

# %%
