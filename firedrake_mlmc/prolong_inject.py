""" Prolonging and injecting multiple levels recursively in a FunctionHierarchy """

from __future__ import division  # Get proper divison
from __future__ import absolute_import

from firedrake import *  # noqa
from firedrake.mg.utils import *  # noqa


def TurnHierarchyIntoItemList(Hierarchy):

    """ Turns a :class:`FunctionHierarchy` into a list of individual :class:`Function`
        of the same :class:`FunctionSpace`

            :param Hierarchy: This is the :class:`FunctionHierarchy` one wants to turn into a list.
            :type Hierarchy: :class:`FunctionHierarchy`

    """

    L = len(Hierarchy)
    List = []

    for i in range(L):
        List.append(Function(Hierarchy[i].function_space()))

    return List


def ProlongUpToFinestLevel(Func, Hierarchy):

    """ Prolongs any :class:`Function` in a :class:`FunctionHierarchy` to the
        finest level in the empty output :class:`FunctionHierarchy`

            :param Function: The :class:`Function` one wants to prolong
            :type Function: :class:`Function`

            :param Hierarchy: To inset the prolonged :class:`Function` into
            :type Hierarchy: :class:`FunctionHierarchy`

    """

    CoarseIndex = get_level(Func)[1]

    if CoarseIndex == -1:
        raise IndexError('Coarse Index is not an index one can prolong from.' +
                         ' Function is not in a hierarchy')

    Hier = TurnHierarchyIntoItemList(Hierarchy)
    Hier[CoarseIndex] = Func

    if CoarseIndex == len(Hierarchy):
        Hier[-1] = Func
    for i in range(len(Hier) - 1 - CoarseIndex):
        prolong(Hier[i + CoarseIndex], Hier[i + 1 + CoarseIndex])

    Hierarchy[-1].assign(Hier[-1])

    return Hierarchy[-1]


def ProlongUpToAnyLevel(Level_to_prolong_to, Func, Hierarchy):

    """ Prolongs any :class:`Function` in a :class:`FunctionHierarchy` to the
        prescribed level in the empty output :class:`FunctionHierarchy`

            :param Level_to_prolong_to: The level to prolong :class`Function` to
            :type Level_to_prolong_to: int

            :param Function: The :class:`Function` one wishes to prolong
            :type Function: :class:`Function`

            :param Hierarchy: To inset the prolonged :class:`Function` into
            :type Hierarchy: :class:`FunctionHierarchy`

    """

    CoarseIndex = get_level(Func)[1]

    if CoarseIndex == -1:
        raise IndexError('Coarse Index is not an index one can prolong from.' +
                         ' Function is not in a hierarchy')

    if CoarseIndex > Level_to_prolong_to:
        raise ValueError('Cant prolong to level less than current')

    Hier = TurnHierarchyIntoItemList(Hierarchy)
    Hier[CoarseIndex] = Func

    if CoarseIndex == Level_to_prolong_to:  # this is the same as Level_to_prolong_to
        Hier[Level_to_prolong_to] = Func
    for i in range(Level_to_prolong_to - CoarseIndex):
        prolong(Hier[CoarseIndex + i], Hier[CoarseIndex + i + 1])

    Hierarchy[Level_to_prolong_to].assign(Hier[Level_to_prolong_to])

    return Hierarchy[Level_to_prolong_to]


def InjectDownToAnyLevel(Level_to_inject_to, Func, Hierarchy):

    """ Injects any :class:`Function` in a :class:`FunctionHierarchy` to the
        prescribed level in the empty output :class:`FunctionHierarchy`

            :param Level_to_inject_to: The level to inject :class:`Function` to
            :type Level_to_inject_to: int

            :param Function: The :class:`Function` one wishes to inject
            :type Function: :class:`Function`

            :param Hierarchy: To inset the injected :class:`Function` into
            :type Hierarchy: :class:`FunctionHierarchy`

    """

    FineIndex = get_level(Func)[1]

    if FineIndex == -1:
        raise IndexError('Fine Index is not an index one can inject from.' +
                         'Function is not in a hierarchy')

    if FineIndex < Level_to_inject_to:
        raise ValueError('Cant inject to level greater than current')

    Hier = TurnHierarchyIntoItemList(Hierarchy)
    Hier[FineIndex] = Func

    if FineIndex == Level_to_inject_to:
        Hier[Level_to_inject_to] = Func
    for i in range(FineIndex - Level_to_inject_to):
        inject(Hier[FineIndex - i], Hier[FineIndex - i - 1])

    Hierarchy[Level_to_inject_to].assign(Hier[Level_to_inject_to])

    return Hierarchy[Level_to_inject_to]
