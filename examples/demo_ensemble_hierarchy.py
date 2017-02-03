""" demo of using the ensemble hierarchy objects """

from __future__ import division

from firedrake import *
from firedrake_mlmc import *

import numpy as np


# Create mesh hierarchy
mesh = UnitSquareMesh(5, 5)
L = 2
refinements_per_level = 2
mesh_h = MeshHierarchy(mesh, L, refinements_per_level=refinements_per_level)

# Create function space / function hierarchies
V_h = tuple([FunctionSpace(m, 'DG', 1) for m in mesh_h])
u_h = []
for i in range(len(V_h)):
    u_h.append(Function(V_h[i]))

# make an ensemble hierarchy
e_h = EnsembleHierarchy(V_h)

# check the refinement factor
print 'refinement factor of FSH = ', e_h.M

# Append several random states to the ensemble hierarchy on the bottom level
coarse = Function(V_h[0])
fine = Function(V_h[1])
for i in range(10):
    x = np.random.normal(0, 1, 1)[0]
    coarse.assign(sin(x))
    fine.assign(sin(x))
    S = State(coarse, fine)
    e_h.AppendToEnsemble(S)

# Next append several random states to the ensemble hierarchy on the finest level
coarse = Function(V_h[1])
fine = Function(V_h[2])
for i in range(10):
    x = np.random.normal(0, 1, 1)[0]
    coarse.assign(sin(x))
    fine.assign(sin(x))
    S = State(coarse, fine)
    e_h.AppendToEnsemble(S)

# print the ensemble hierarchy lengths
print e_h.nl

# print the ensemble hierarchy resolutions
print 'refinement of M =', e_h, ', ensemble hierarchy resolutions = ', e_h.dxl

# Update statistics

e_h.UpdateStatistics()

# print the data of the multilevel expectation
print 'Multilevel Exp cell basis coeffs: ', e_h.MultilevelExpectation.dat.data

# should be the same as just the first level mean
print 'First level mean cell basis coeffs: ', e_h.Mean[0].dat.data

# with both variances set to 0
print 'With variance in cell basis coeffs on first level: ', e_h.Variance[0].dat.data

print 'but none on the second: ', e_h.Variance[1].dat.data

print ('And thus the multilevel Monte Carlo approximation variance is the first level sample ' +
       'divided by the number of samples on that level: ', e_h.MultilevelVariance.dat.data)

# now clear the ensemble
e_h.UpdateStatistics(clear_ensemble=True)

print 'Even though the ensemble hierarchy has been cleared, ', e_h._hierarchy

print 'the online size of the ensembles used to compute the statistics still remain: ', e_h.nl
