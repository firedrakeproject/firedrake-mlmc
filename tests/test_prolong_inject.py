from firedrake_mlmc import *

import pytest


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_prolong_to_finest_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    G = Function(V[0]).interpolate(Expression("3"))
    F = Function(V[len(M) - 1]).interpolate(Expression("3"))

    H = Function(V[-1])

    prolong(G, H)

    assert get_level(H.function_space())[1] == len(M) - 1

    assert norm(assemble(H - F)) <= 0


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_prolong_same_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    F = Function(V[0]).interpolate(Expression("3"))
    G = Function(V[0]).interpolate(Expression("3"))

    H = Function(V[0])

    prolong(F, H)

    assert norm(assemble(H - G)) <= 0


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_inject_same_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    F = Function(V[0]).interpolate(Expression("3"))
    G = Function(V[0]).interpolate(Expression("3"))

    H = Function(V[0])

    inject(F, H)

    assert norm(assemble(H - G)) <= 0


def test_prolong_with_non_hierarchy_function():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)
    v = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    F = Function(V[0]).interpolate(Expression("3"))
    H = Function(v[-1])

    a = 0
    try:
        prolong(F, H)

    except ValueError:
        a = 1

    assert a == 1


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_prolong_to_any_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    F = Function(V[0]).interpolate(Expression("3"))
    G = Function(V[2]).interpolate(Expression("3"))

    H = Function(V[2])

    prolong(F, H)

    assert get_level(H.function_space())[1] == 2

    assert norm(assemble(H - G)) <= 0


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_inject_to_any_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)

    # check for prolonging to top level
    F = Function(V[-1]).interpolate(Expression("3"))
    G = Function(V[0]).interpolate(Expression("3"))

    H = Function(V[0])

    inject(F, H)

    assert get_level(H.function_space())[1] == 0

    assert norm(assemble(H - G)) <= 0


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
