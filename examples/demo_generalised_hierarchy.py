""" demo showing generalised mesh and function space hierarchies """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *

import numpy as np


mesh = UnitIntervalMesh(5)

# number of refinement levels
L = 3

# refinement factor
M = 4

# create a generalised mesh hierarchy with refinement factor 4 instead of 2
GMH = GeneralisedMeshHierarchy(mesh, L, M)

# note how the intermediate levels are stored
OMH = GMH._full_hierarchy

# note the decreasing resolution (min cell edge length) is now refined by factor of 4
min_cell_edge_length = np.zeros(L + 1)
for i in range(L + 1):
    min_cell_edge_length[i] = MinDx(GMH[i])

print 'min cell edge lengths refined by factor of 4: ', min_cell_edge_length

# now create a DG0 function space hierarchy on this generalised mesh hierarchy
GFSH = GeneralisedFunctionSpaceHierarchy(GMH, 'DG', 1)

# again note how the intermediate levels of the function space hierarchy are stored
OFSH = GFSH._full_hierarchy

# one can build a function hierarchy of the generalised function space hierarchy
u_h = []
for i in range(len(GFSH)):
    u_h.append(Function(GFSH[i]))

# as an example, let's interpolate an expression onto each one


def interp_exp(x):
    return sin(x[0])

# exact solution to see errors of evaluation of sin(pi / 6)
exact = np.sin(np.pi / 6)
error = np.zeros(L + 1)

for i in range(L + 1):
    u_h[i].interpolate(interp_exp(SpatialCoordinate(GMH[i])))
    error[i] = np.abs(u_h[i].at(np.pi / 6) - exact)

print('error of evaluation of (pi / 6) with decreasing resolution ' +
      'refined by factor of 4: ', error)
