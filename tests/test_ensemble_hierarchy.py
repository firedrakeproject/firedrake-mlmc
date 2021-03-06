""" tests the ensemble hierarchy object """

from firedrake_mlmc import *

import pytest


def test_ensemble_hierarchy_same_fs_check():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]
    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    EH = EnsembleHierarchy(V_h)

    S = State(u_h_1, u_h_2)
    EH.AppendToEnsemble(S)

    assert EH[0][0][0].function_space() == EH[0][0][1].function_space()


def test_ensemble_hierarchy_constants_1():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]
    u_h_1 = Constant(0, domain=MH[0])
    u_h_2 = Constant(0, domain=MH[1])

    S = State(u_h_1, u_h_2)

    EH = EnsembleHierarchy(V_h)

    a = 0
    b = 1

    try:
        EH.AppendToEnsemble(S)

    except TypeError:
        a = 1
        b = 0

    assert a > b


def test_ensemble_hierarchy_constants_2():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]
    u_h_1 = Constant(0, domain=MH[0])
    u_h_2 = Constant(0, domain=MH[1])

    S = State(u_h_1, u_h_2)

    EH = EnsembleHierarchy(V_h, state_type=Constant)

    a = 0
    b = 1

    try:
        EH.AppendToEnsemble(S)

    except TypeError:
        a = 1
        b = 0

    assert a < b


def test_ensemble_hierarchy_attr():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]
    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    EH = EnsembleHierarchy(V_h)

    assert EH.M == 2

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    assert len(EH.nl) == 1

    assert EH.nl[0] == 10

    assert len(EH.dxl) == 1

    assert np.abs(EH.dxl[0] - np.sqrt((0.05 ** 2))) < 1e-8

    assert EH.L == 1


def test_ensemble_hierarchy_skip_levels():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[2])
    u_h_2 = Function(V_h[3])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    assert EH.L == 3

    assert len(EH.dxl) == 3

    assert len(EH.nl) == 3

    assert EH.nl[2] == 10

    assert EH.nl[0] == 0 and EH.nl[1] == 0


def test_ensemble_hierarchy_refinement():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1])
    u_h_2 = Function(V_h[2])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    assert EH.L == 2

    assert EH.L == len(EH)

    assert np.sum(np.abs((EH.dxl / EH.dxl[0]) - (2 ** - np.linspace(0, 1, 2)))) < 1e-8


def test_clear_ensemble():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1])
    u_h_2 = Function(V_h[2])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics(clear_ensemble=True)

    assert len(EH) == 2

    assert EH[0] == []

    assert EH[1] == []

    assert EH.nl == [10, 10]


def test_clear_identifier_hierarchy():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1])
    u_h_2 = Function(V_h[2])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics(clear_ensemble=True)

    assert len(EH._identifier_hierarchy) == 2

    assert EH._identifier_hierarchy[0] == []

    assert EH._identifier_hierarchy[1] == []


def test_keep_ensemble():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1])
    u_h_2 = Function(V_h[2])

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics()

    assert len(EH) == 2

    assert len(EH[0]) == 10

    assert len(EH[1]) == 10

    assert EH.nl == [10, 10]


def test_mean():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics()

    assert np.all(EH.MultilevelExpectation.dat.data == 3) == 1

    assert np.all(EH.Mean[0].dat.data == 2) == 1

    assert np.all(EH.Mean[1].dat.data == 1) == 1


def test_variance():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics()

    assert np.all(EH.MultilevelVariance.dat.data == 0) == 1

    assert np.all(EH.Variance[0].dat.data == 0) == 1

    assert np.all(EH.Variance[1].dat.data == 0) == 1


def test_mean_variance_compute():

    M = UnitSquareMesh(10, 10)
    L = 1
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0])
    u_h_2 = Function(V_h[1])

    for i in range(2):

        u_h_1.assign(i)
        u_h_2.assign(i)
        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    exact_mean = 0.5
    exact_variance = 0.5 - 0.25

    EH.UpdateStatistics()

    assert np.all(EH.Variance[0].dat.data == exact_variance) == 1

    assert np.all(EH.Mean[0].dat.data == exact_mean) == 1

    assert np.all(EH.MultilevelExpectation.dat.data == exact_mean) == 1


def test_identifier_hierarchy():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)
        assert EH._identifier_hierarchy[0][i] == S.identifier
        assert S.index == tuple([0, i])

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)
        assert EH._identifier_hierarchy[1][i] == S.identifier
        assert S.index == tuple([1, i])

    assert len(EH._identifier_hierarchy) == len(EH._hierarchy)
    assert len(EH._identifier_hierarchy[0]) == len(EH._hierarchy[0])
    assert len(EH._identifier_hierarchy[1]) == len(EH._hierarchy[1])


def test_state_hierarchy_1():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    ens1 = []
    ens2 = []

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)
        ens1.append(S)
        assert np.max(np.abs(EH._state_hierarchy[0][i].state[1].dat.data -
                      ens1[i].state[1].dat.data)) < 1e-5
        assert np.max(np.abs(EH._state_hierarchy[0][i].identifier - ens1[i].identifier)) < 1e-5

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)
        ens2.append(S)
        assert np.max(np.abs(EH._state_hierarchy[1][i].state[1].dat.data -
                      ens2[i].state[1].dat.data)) < 1e-5
        assert np.max(np.abs(EH._state_hierarchy[1][i].identifier - ens2[i].identifier)) < 1e-5


def test_updating_ensemble_members():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics()
    Mean = np.copy(EH.MultilevelExpectation.dat.data[:])

    # change the states
    for i in range(10):

        EH._state_hierarchy[1][i].state[1].assign(4)
        EH.UpdateEnsembleMember(EH._state_hierarchy[1][i])

    EH.UpdateStatistics()
    NewMean = EH.MultilevelExpectation.dat.data[:]
    assert np.abs(np.max((NewMean - Mean) - 1.0)) < 1e-5


def test_state_hierarchy_2():

    M = UnitSquareMesh(10, 10)
    L = 3
    MH = MeshHierarchy(M, L)

    V_h = [FunctionSpace(m, 'DG', 1) for m in MH]

    EH = EnsembleHierarchy(V_h)

    u_h_1 = Function(V_h[0]).assign(1)
    u_h_2 = Function(V_h[1]).assign(2)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    u_h_1 = Function(V_h[1]).assign(2)
    u_h_2 = Function(V_h[2]).assign(3)

    for i in range(10):

        S = State(u_h_1, u_h_2)
        EH.AppendToEnsemble(S)

    EH.UpdateStatistics()

    # change the states
    for i in range(10):

        EH._state_hierarchy[1][i].state[1].assign(4)
        EH.UpdateEnsembleMember(EH._state_hierarchy[1][i])
        assert np.max(np.abs(EH._state_hierarchy[1][i].state[1].dat.data - 4.0)) < 1e-5
        assert EH._state_hierarchy[1][i].identifier == EH._identifier_hierarchy[1][i]


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
