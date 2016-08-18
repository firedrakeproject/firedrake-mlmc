from firedrake_mlmc import *

import pytest


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_generalised_prolong():

    M = 4

    L = 2

    GMH = GeneralisedMeshHierarchy(UnitSquareMesh(10, 10), L, M)

    V = GeneralisedFunctionSpaceHierarchy(GMH, 'DG', 0)

    # check by prolonging to finest level in generalised function hierarchy
    F_1 = Function(V[0]).interpolate(Expression("3"))
    F_2 = Function(V[2]).interpolate(Expression("3"))
    G_1 = Function(V._full_hierarchy[0]).interpolate(Expression("3"))
    G_2 = Function(V._full_hierarchy[4]).interpolate(Expression("3"))

    prolong(F_1, F_2)
    prolong(G_1, G_2)

    assert get_level(F_2.function_space())[1] == 2
    assert get_level(G_2.function_space())[1] == 4

    assert norm(assemble(F_2 - G_2)) <= 0


@pytest.mark.xfail(reason="recursive prolong not merged in master firedrake")
def test_generalised_inject():

    M = 4

    L = 2

    M = GeneralisedMeshHierarchy(UnitSquareMesh(10, 10), L, M)

    V = GeneralisedFunctionSpaceHierarchy(M, 'DG', 0)

    # check by injecting to coarsest level in generalised function hierarchy
    F_1 = Function(V[2]).interpolate(Expression("3"))
    F_2 = Function(V[0]).interpolate(Expression("3"))
    G_1 = Function(V._full_hierarchy[4]).interpolate(Expression("3"))
    G_2 = Function(V._full_hierarchy[0]).interpolate(Expression("3"))

    inject(F_1, F_2)
    inject(G_1, G_2)

    assert get_level(F_2.function_space())[1] == 0
    assert get_level(G_2.function_space())[1] == 0

    assert norm(assemble(G_2 - F_2)) <= 0


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
