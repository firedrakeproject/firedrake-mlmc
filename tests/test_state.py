from firedrake_mlmc import *

import pytest


def test_state_levels_1():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = [FunctionSpace(m, 'DG', 0) for m in MH]

    F = Function(V[0])
    G = Function(V[1])

    b = 1
    a = 0

    try:
        State(F, G)

    except Warning:
        a = 1
        b = 0

    assert a < b


def test_constants():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    F = Constant(0, domain=MH[0])
    G = Constant(0, domain=MH[1])

    S = State(F, G)

    assert S.levels == tuple([0, 1])


def test_inputs():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = [FunctionSpace(m, 'DG', 0) for m in MH]

    F = Function(V[0])
    G = Constant(0, domain=MH[1])

    b = 1
    a = 0

    try:
        State(F, G)

    except TypeError:
        a = 1
        b = 0

    assert a > b


def test_state_levels_2():

    M = UnitSquareMesh(10, 10)

    V = FunctionSpace(M, 'DG', 0)

    F = Function(V)
    G = Function(V)

    b = 1
    a = 0

    try:
        State(F, G)

    except ValueError:
        a = 1
        b = 0

    assert a > b


def test_state_levels_3():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = [FunctionSpace(m, 'DG', 0) for m in MH]

    F = Function(V[0])
    G = Function(V[1])

    S = State(F, G)

    assert S.levels == tuple([0, 1])


def test_state_levels_4():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 2)
    V = [FunctionSpace(m, 'DG', 0) for m in M]

    F = Function(V[0])
    G = Function(V[2])

    b = 1
    a = 0

    try:
        State(F, G)

    except ValueError:
        a = 1
        b = 0

    assert a > b


def test_index_and_identifier():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = [FunctionSpace(m, 'DG', 0) for m in MH]

    F = Function(V[0])
    G = Function(V[1])

    S = State(F, G)

    assert S.identifier is None
    assert S.index is None

    EH = EnsembleHierarchy(V)

    EH.AppendToEnsemble(S)

    assert isinstance(S.identifier, int)
    assert isinstance(S.index, tuple)
    assert len(S.index) == 2


def test_reassigning_state():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = [FunctionSpace(m, 'DG', 0) for m in MH]

    F = Function(V[0])
    G = Function(V[1])

    S = State(F, G)

    EH = EnsembleHierarchy(V)

    EH.AppendToEnsemble(S)

    S.state[0].assign(1.0)
    S.state[1].assign(1.0)

    assert isinstance(S.identifier, int)
    assert isinstance(S.index, tuple)

    assert np.max(np.abs(S.state[0].dat.data - 1.0)) < 1e-5
    assert np.max(np.abs(S.state[1].dat.data - 1.0)) < 1e-5


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
