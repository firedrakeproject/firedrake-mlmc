""" utilities for multigrid """

from __future__ import division  # Get proper divison
from __future__ import absolute_import

from firedrake import *  # noqa
from firedrake.mg.utils import *  # noqa


def get_refined_level(obj):
    """ this is effectively a wrapper for get_level which finds the level of the actual
        level of the refined hierarchy rather than the full hierarchy of the object

        If no level info is available, return ``None, -1``.
    """

    # get hierarchy and level information for object
    hierarchy, lvl = get_level(obj)

    # find refined level
    try:
        refinements_per_level = hierarchy.refinements_per_level
        # find out if part of a full hierarchy
        if (lvl / refinements_per_level) % 1 != 0:
            raise ValueError('object is not from a refined hierarchy, instead a full hierarchies')
        else:
            new_lvl = int(lvl / refinements_per_level)
        return hierarchy, new_lvl

    except AttributeError:
        return None, -1
