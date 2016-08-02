""" demo showing prolonging and injecting between levels on a generalised mesh hierarchy """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *


mesh = UnitIntervalMesh(5)

# number of refinement levels
L = 1

# refinement factor
M = 4

# create a generalised mesh hierarchy with refinement factor 4 instead of 2
GMH = GeneralisedMeshHierarchy(mesh, L, M)

# now create a DG0 function space hierarchy on this generalised mesh hierarchy
GFSH = GeneralisedFunctionSpaceHierarchy(GMH, 'DG', 1)

# note how the intermediate levels of the function space hierarchy are stored
OFSH = GFSH._full_hierarchy

# one can build a function hierarchy of the generalised function space hierarchy
u_h = []
for i in range(len(GFSH)):
    u_h.append(Function(GFSH[i]))

# create another one to output prolonged and injected solutions to
u_h_ = []
for i in range(len(GFSH)):
    u_h_.append(Function(GFSH[i]))

# Interpolate an expression to the coarsest level of u_h
x = SpatialCoordinate(GMH[0])
u_h[0].interpolate(sin(x[0] * 2 * pi))

# prolong to the finest output level
prolong(u_h[0], u_h_[1])

# Print error
x = SpatialCoordinate(GMH[1])
Exact = Function(GFSH[1]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[1]))

# Interpolate an expression to the finest level of u_h
x = SpatialCoordinate(GMH[1])
u_h[1].interpolate(sin(x[0] * 2 * pi))

# inject to the coarsest output level
inject(u_h[1], u_h_[0])

# Print error - this should be 0 as we're using the finer solution
x = SpatialCoordinate(GMH[0])
Exact = Function(GFSH[0]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[0]))
