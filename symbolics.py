from sympy import symbols, Eq, simplify, curl, cross
from sympy.vector import CoordSys3D, Vector
from sympy.plotting import plot3d

# Define symbols
z, x1, x2, x3 = symbols('z x1 x2 x3')
p1, p2 = symbols('p1 p2')

# Define the coordinate system
N = CoordSys3D('N')

# Define B0 and A as functions
B0 = 1 + z**2
A = ((-1/2) * B0 * x2, (1/2) * B0 * x1, 0)

# Calculate B as the curl of A and simplify
B = curl(A, N).simplify()

# Define e1 and e2 as the cross products of B with (1, 0, 0) and themselves, respectively
e1 = cross(B, N.i).simplify()
e2 = cross(B, e1).simplify()

# Create a VectorPlot3D of B, e1, and e2
plot3d(B.dot(N.i), B.dot(N.j), B.dot(N.k), (x1, -1, 1), (x2, -1, 1), (x3, -1, 1))
plot3d(e1.dot(N.i), e1.dot(N.j), e1.dot(N.k), (x1, -1, 1), (x2, -1, 1), (x3, -1, 1))
plot3d(e2.dot(N.i), e2.dot(N.j), e2.dot(N.k), (x1, -1, 1), (x2, -1, 1), (x3, -1, 1))

# Define p1 and p2 as expressions
p1_expr = ((p1 + A[0])**2 + (p2 + A[1])**2) / 2
p2_expr = ((p1 + A[0])**2 + (p2 + A[1])**2) / 2

# Simplify p1 and p2 expressions
p1_simplified = simplify(p1_expr)
p2_simplified = simplify(p2_expr)
