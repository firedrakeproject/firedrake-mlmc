""" demo showing prolonging and injecting between levels on a generalised mesh hierarchy """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *


mesh = UnitIntervalMesh(5)

# number of refinement levels
L = 1

# refinements per levels
refinements_per_level = 2

# create a generalised mesh hierarchy with refinement factor 4 instead of 2
MH = MeshHierarchy(mesh, L, refinements_per_level=refinements_per_level)

# now create a DG0 function space hierarchy on this generalised mesh hierarchy
FSH = [FunctionSpace(m, 'DG', 1) for m in MH]

# note how the intermediate levels of the function space hierarchy are stored
OFSH = MH._unskipped_hierarchy

# one can build a function hierarchy of the generalised function space hierarchy
u_h = []
for i in range(len(FSH)):
    u_h.append(Function(FSH[i]))

# create another one to output prolonged and injected solutions to
u_h_ = []
for i in range(len(FSH)):
    u_h_.append(Function(FSH[i]))

# Interpolate an expression to the coarsest level of u_h
x = SpatialCoordinate(MH[0])
u_h[0].interpolate(sin(x[0] * 2 * pi))

# prolong to the finest output level
prolong(u_h[0], u_h_[1])

# Print error
x = SpatialCoordinate(MH[1])
Exact = Function(FSH[1]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[1]))

# Interpolate an expression to the finest level of u_h
x = SpatialCoordinate(MH[1])
u_h[1].interpolate(sin(x[0] * 2 * pi))

# inject to the coarsest output level
inject(u_h[1], u_h_[0])

# Print error - this should be 0 as we're using the finer solution
x = SpatialCoordinate(MH[0])
Exact = Function(FSH[0]).interpolate(sin(x[0] * 2 * pi))
print norm(assemble(Exact - u_h_[0]))
