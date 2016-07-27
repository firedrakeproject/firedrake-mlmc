""" demo showing how to create a state for different levels in a mesh hierarchy """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *


# Create mesh hierarchy
mesh = UnitSquareMesh(5, 5)
L = 3
mesh_h = MeshHierarchy(mesh, L)

# Create function space / function hierarchies
V_h = FunctionSpaceHierarchy(mesh_h, 'DG', 1)
u_h = FunctionHierarchy(V_h)

# Loop through levels creating states and printing levels
for i in range(L):

    # Create a state
    S = State(u_h[i], u_h[i+1])

    # print levels
    print S.levels
