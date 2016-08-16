""" demo showing recursive prolonging and injecting """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *


# Create mesh hierarchy
mesh = UnitSquareMesh(5, 5)
L = 3
mesh_h = MeshHierarchy(mesh, L)

# Create function space / function hierarchies
V_h = FunctionSpaceHierarchy(mesh_h, 'DG', 1)
u_h = []
for i in range(len(V_h)):
    u_h.append(Function(V_h[i]))

# Create a function hierarchy to output to
u_h_ = []
for i in range(len(V_h)):
    u_h_.append(Function(V_h[i]))

# Interpolate an expression to the coarsest level of u_h
x = SpatialCoordinate(mesh_h[0])
u_h[0].interpolate(sin(x[0] * 2 * pi))

# Firstly prolong to the third coarsest level (skip one level)
prolong(u_h[0], u_h_[2])

# Print error
x = SpatialCoordinate(mesh_h[2])
Exact = Function(V_h[2]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[2]))

# Create another function hierarchy to output to
# one can build a function hierarchy of the generalised function space hierarchy
u_h_ = []
for i in range(len(V_h)):
    u_h_.append(Function(V_h[i]))

# Then prolong to the finest level
prolong(u_h[0], u_h_[-1])

# Print error
x = SpatialCoordinate(mesh_h[-1])
Exact = Function(V_h[-1]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[-1]))
