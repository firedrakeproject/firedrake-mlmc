from firedrake_mlmc import *

import pytest


def test_refined_level():

    refinements_per_level = 2
    L = 3

    M = MeshHierarchy(UnitSquareMesh(4, 4), L, refinements_per_level=refinements_per_level)
    V = FunctionSpaceHierarchy(M, 'DG', 0)

    for i in range(L + 1):
        assert get_refined_level(V[i])[1] == i

    for i in range(int(refinements_per_level * L) + 1):
        if (i % refinements_per_level) == 0:
            assert get_refined_level(V._full_hierarchy[i])[1] == int(i / refinements_per_level)
        else:
            with pytest.raises(ValueError):
                get_refined_level(V._full_hierarchy[i])


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
