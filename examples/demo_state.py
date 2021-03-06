""" demo showing how to create a state for different levels in a mesh hierarchy """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *


# Create mesh hierarchy
mesh = UnitSquareMesh(5, 5)
L = 3
mesh_h = MeshHierarchy(mesh, L)

# Create function space / function hierarchies
V_h = [FunctionSpace(m, 'DG', 1) for m in mesh_h]

u_h = []
for i in range(len(V_h)):
    u_h.append(Function(V_h[i]))


# Loop through levels creating states and printing levels
for i in range(L):

    # Create a state
    S = State(u_h[i], u_h[i+1])

    # print levels
    print S.levels
