""" demo of using the ensemble hierarchy object and updating ensemble members """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *

import matplotlib.pyplot as plt

import numpy as np


# Create mesh hierarchy
mesh = UnitIntervalMesh(5)
L = 2
refinements_per_level = 2
mesh_h = MeshHierarchy(mesh, L, refinements_per_level=refinements_per_level)

# Create function space / function hierarchies
V_h = tuple([FunctionSpace(m, 'DG', 0) for m in mesh_h])
u_h = []
for i in range(len(V_h)):
    u_h.append(Function(V_h[i]))

# Make an ensemble hierarchy
e_h = EnsembleHierarchy(V_h)

# Append several random states to the ensemble hierarchy on the bottom level
coarse = Function(V_h[0])
fine = Function(V_h[1])
for i in range(10):
    # Find the coordinates
    xc = SpatialCoordinate(mesh_h[0])
    xf = SpatialCoordinate(mesh_h[1])

    # Random shift
    x_shift = np.random.normal(0, 0.1, 1)[0]
    coarse.interpolate(sin(2 * pi * (xc[0] + x_shift)))
    fine.interpolate(sin(2 * pi * (xf[0] + x_shift)))
    S = State(coarse, fine)
    e_h.AppendToEnsemble(S)

# Next append several random states to the ensemble hierarchy on the finest level
coarse = Function(V_h[1])
fine = Function(V_h[2])
for i in range(10):
    # Find the coordinates
    xc = SpatialCoordinate(mesh_h[1])
    xf = SpatialCoordinate(mesh_h[2])

    # Random shift
    x_shift = np.random.normal(0, 0.1, 1)[0]
    coarse.interpolate(sin(2 * pi * (xc[0] + x_shift)))
    fine.interpolate(sin(2 * pi * (xf[0] + x_shift)))
    S = State(coarse, fine)
    e_h.AppendToEnsemble(S)

# Update statistics
e_h.UpdateStatistics()

# Plot the mlmc mean
plot(e_h.MultilevelExpectation)
plt.axis([0, 1, -3, 3])

# Update ensemble hierarchy with a random shift in the y axis
for i in range(10):
    x = np.random.normal(0.25, 0.1, 1)

    # Update state from state hierarchy
    e_h._state_hierarchy[0][i].state[0].dat.data[:] += x
    e_h._state_hierarchy[0][i].state[1].dat.data[:] += x
    e_h.UpdateEnsembleMember(e_h._state_hierarchy[0][i])

for i in range(10):
    x = np.random.normal(0.25, 0.1, 1)

    # Update state from state hierarchy
    e_h._state_hierarchy[1][i].state[0].dat.data[:] += x
    e_h._state_hierarchy[1][i].state[1].dat.data[:] += x
    e_h.UpdateEnsembleMember(e_h._state_hierarchy[1][i])

# Update statistics
e_h.UpdateStatistics()

# Plot the mlmc mean
plot(e_h.MultilevelExpectation)
plt.axis([0, 1, -3, 3])
plt.show()
