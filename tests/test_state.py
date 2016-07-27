from firedrake_mlmc import *


def test_state_levels_1():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = FunctionSpaceHierarchy(MH, 'DG', 0)

    F = FunctionHierarchy(V)

    b = 1
    a = 0

    try:
        State(F[0], F[1])

    except Warning:
        a = 1
        b = 0

    assert a < b


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

    F = FunctionHierarchy(V)

    S = State(F[0], F[1])

    assert S.levels == tuple([0, 1])


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
