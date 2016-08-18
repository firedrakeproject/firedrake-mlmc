from firedrake_mlmc import *


def test_state_levels_1():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = FunctionSpaceHierarchy(MH, 'DG', 0)

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

    V = FunctionSpaceHierarchy(MH, 'DG', 0)

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

    except Warning:
        a = 1
        b = 0

    assert a > b


def test_state_levels_3():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = FunctionSpaceHierarchy(MH, 'DG', 0)

    F = Function(V[0])
    G = Function(V[1])

    S = State(F, G)

    assert S.levels == tuple([0, 1])


def test_state_levels_4():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 2)
    V = FunctionSpaceHierarchy(M, 'DG', 0)

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


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
