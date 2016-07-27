from firedrake_mlmc import *


def test_prolong_to_finest_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)
    F = FunctionHierarchy(V)

    # check for prolonging to top level
    F[0].interpolate(Expression("3"))
    F[len(M) - 1].interpolate(Expression("3"))

    H = FunctionHierarchy(V)

    A = ProlongUpToFinestLevel(F[0], H)

    assert get_level(A)[1] == len(M) - 1

    assert norm(assemble(A - F[len(M) - 1])) <= 0


def test_prolong_with_non_hierarchy_function():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)
    F = Function(V[0])

    # check for prolonging to top level
    F.interpolate(Expression("3"))

    H = FunctionHierarchy(V)

    a = 0
    try:
        ProlongUpToFinestLevel(F, H)

    except IndexError:
        a = 1

    assert a == 1


def test_prolong_to_any_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)
    F = FunctionHierarchy(V)

    # check for prolonging to top level
    F[0].interpolate(Expression("3"))
    F[2].interpolate(Expression("3"))

    H = FunctionHierarchy(V)

    A = ProlongUpToAnyLevel(2, F[0], H)

    assert get_level(A)[1] == 2

    assert norm(assemble(A - F[2])) <= 0


def test_inject_to_any_level():

    M = MeshHierarchy(UnitSquareMesh(10, 10), 3)

    V = FunctionSpaceHierarchy(M, 'DG', 0)
    F = FunctionHierarchy(V)

    # check for prolonging to top level
    F[-1].interpolate(Expression("3"))
    F[0].interpolate(Expression("3"))

    H = FunctionHierarchy(V)

    A = InjectDownToAnyLevel(0, F[-1], H)

    assert get_level(A)[1] == 0

    assert norm(assemble(A - F[0])) <= 0


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
