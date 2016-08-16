""" builds a object to store the fine and coarse state vectors """

from __future__ import division  # Get proper divison
from __future__ import absolute_import

from firedrake import *  # noqa
from firedrake.mg.utils import *  # noqa


class State(object):

    """ Tuple that stores the state :class:`Function` and properties of them for the coarse (index 0) and fine (index 1) discretizations for a single realisation of the user define system

        :param input_1: Coarse state :class:`Function` needed for discretization
        :type input_1: :class:`Function`

        :param input_2: Fine state :class:`Function` needed for discretization
        :type input_2: :class:`Function`

    """

    def __init__(self, input_1, input_2):

        # tuple of coarse and fine state :class:`Function`
        self.state = tuple([input_1, input_2])

        # give the state the attributes the levels of each fine / coarse
        # solution. be careful for states which are lists of multiple functions.
        self.levels = tuple([get_level(self.state[0].function_space())[1],
                             get_level(self.state[1].function_space())[1]])

        # add check that both inputs are Functions
        if type(self.state[0]) != Function:
            raise TypeError('state input is not a Function')

        if type(self.state[1]) != Function:
            raise TypeError('state input is not a Function')

        # add check for levels - non fatal -> perhaps change to actual error
        if self.levels[0] == -1 or self.levels[1] == -1:
            raise Warning('Levels of state may not be actual hierarchal levels.' +
                          ' Check if they belong to FunctionHierarchy! get_level has failed.')

        # check for non consecutive levels
        if self.levels[1] - 1 != self.levels[0]:
            raise ValueError('levels of inputs are not consecutive')

        super(State, self).__init__()
