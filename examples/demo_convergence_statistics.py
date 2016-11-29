""" Convergence demo of multilevel monte carlo estimates of statistics using
the ensemble hierarchy"""

from __future__ import division

from firedrake import *
from firedrake_mlmc import *

import matplotlib.pyplot as plot

import numpy as np


# define the mesh hierarchy
L = 3
mesh_hierarchy = MeshHierarchy(UnitSquareMesh(1, 1), L)

# define each level means
means = np.ones(L + 1)
means[0: L - 1] = 2 ** (-np.linspace(0, L - 2, L - 1))
means[-1] = 0

# deifne the function space hierarchy
fs_hierarchy = tuple([FunctionSpace(m, 'DG', 0) for m in mesh_hierarchy])


# define the function to calculate the multilevel monte carlo estimate of mean
def mlmc_estimate(N0, fs_hierarchy, means, state_type=Function):

    eh = EnsembleHierarchy(fs_hierarchy, state_type=state_type)
    L = len(fs_hierarchy) - 1
    ns = np.zeros(L)

    for i in range(L):
        n = int(N0 / (2 ** i))
        ns[i] = n
        for k in range(n):
            xc = np.random.normal(0, 1, 1)[0] + means[i]
            xf = np.random.normal(0, 1, 1)[0] + means[i + 1]
            if state_type == Function:
                coarse = Function(fs_hierarchy[i]).assign(xc)
                fine = Function(fs_hierarchy[i + 1]).assign(xf)
            if state_type == Constant:
                coarse = Constant(xc, domain=fs_hierarchy[i].mesh())
                fine = Constant(xf, domain=fs_hierarchy[i + 1].mesh())
            s = State(coarse, fine)
            eh.AppendToEnsemble(s)

    eh.UpdateStatistics()

    # check that all basis coeffs are same
    assert np.abs(np.mean(eh.MultilevelExpectation.dat.data[:]) -
                  eh.MultilevelExpectation.dat.data[0]) < 1e-4

    # find variance
    var = 0
    for i in range(L):
        var += np.mean(eh.Variance[i].dat.data[:]) / ns[i]
    mean = np.mean(eh.MultilevelExpectation.dat.data)

    return mean, var


# run the convergence loops with increasing n -> do once for constant and once
# for functions as state_type

s = 5
niter = 10

N0s = (12 * (2 ** np.linspace(0, s - 1, s))).astype(int)

rmse_const = np.zeros(s)
var_const = np.zeros(s)

rmse_func = np.zeros(s)
var_func = np.zeros(s)


for i in range(s):

    temp_rmse_const = np.zeros(niter)
    temp_var_const = np.zeros(niter)

    temp_rmse_func = np.zeros(niter)
    temp_var_func = np.zeros(niter)

    for j in range(niter):

        m_func, v_func = mlmc_estimate(N0s[i], fs_hierarchy, means, Function)
        m_const, v_const = mlmc_estimate(N0s[i], fs_hierarchy, means, Constant)

        temp_rmse_func[j] = np.square(m_func)
        temp_rmse_const[j] = np.square(m_const)

        temp_var_func[j] = v_func
        temp_var_const[j] = v_const

    rmse_func[i] = np.sqrt(np.mean(temp_rmse_func))
    rmse_const[i] = np.sqrt(np.mean(temp_rmse_const))

    var_func[i] = np.mean(temp_var_func)
    var_const[i] = np.mean(temp_var_const)

    print 'done sample size: ', N0s[i]

# plot results

plot.loglog(N0s, rmse_func, 'r*-')
plot.loglog(N0s, var_func, 'r*--')

plot.loglog(N0s, rmse_const, 'bo-')
plot.loglog(N0s, var_const, 'bo--')

plot.loglog(N0s, 1e1 * N0s.astype(float) ** (-1), 'k--')
plot.loglog(N0s, 5e0 * N0s.astype(float) ** (-0.5), 'k-')

plot.legend(['rmse (function)', 'variance (function)', 'rmse (constant)',
             'variance (constant)', 'linear decay', 'sqrt decay'])

plot.show()
