""" builds a object to store the fine and coarse state vectors """

from __future__ import division  # Get proper divison
from __future__ import absolute_import

from firedrake import *  # noqa
from firedrake.mg.utils import *  # noqa


class State(object):

    """ Tuple that stores the state :class:`Function` or :class:`Constant` and properties
        of them for the coarse (index 0) and fine (index 1) discretizations for a single
        realisation of the user define system

        :param input_1: Coarse state :class:`Function` or :class:`Constant` needed for discretization
        :type input_1: :class:`Function` or :class:`Constant`

        :param input_2: Fine state :class:`Function` or :class:`Constant` needed for discretization
        :type input_2: :class:`Function` or :class:`Constant`

    """

    def __init__(self, input_1, input_2):

        # tuple of coarse and fine state :class:`Function` or :class:`Constant`
        self.state = tuple([input_1, input_2])

        # add check that both inputs are :class:`Function` or :class:`Constant`
        if not isinstance(self.state[0], Function) and not isinstance(self.state[0], Constant):
            raise TypeError('state input is not a Function or Constant')

        if not isinstance(self.state[1], Function) and not isinstance(self.state[1], Constant):
            raise TypeError('state input is not a Function or Constant')

        # finally check if both states are the same type
        if not isinstance(self.state[1], type(self.state[0])):
            raise TypeError('state inputs are not same type')

        # if they are a :class:`Constant` then they need domain
        if isinstance(self.state[1], Constant) or isinstance(self.state[0], Constant):
            if self.state[1].ufl_domain() is None or self.state[0].ufl_domain() is None:
                raise AttributeError('state Constants do not have an associated mesh')

        # give the state the attributes the levels of each fine / coarse solution
        self.levels = tuple([get_level(self.state[0].ufl_domain())[1],
                             get_level(self.state[1].ufl_domain())[1]])

        # add check for levels
        if self.levels[0] is None or self.levels[1] is None:
            raise ValueError('State not part of actual hierarchy. get_level has failed.')

        # check for non consecutive levels
        if self.levels[1] - 1 != self.levels[0]:
            raise ValueError('levels of inputs are not consecutive')

        super(State, self).__init__()
