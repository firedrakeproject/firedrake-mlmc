""" Generates an ensemble hierarchy object """

from __future__ import absolute_import

from firedrake_mlmc.min_dx import *

from firedrake import *  # noqa
from firedrake.mg.utils import *  # noqa


class EnsembleHierarchy(object):

    """ This class creates and builds a hierarchy of ensembles as well as containing properties of them, such as statistics.

        :param function_spaces: Tuple of function spaces corresponding to each level of a :class:`MeshHierarchy`.
        :type function_spaces: Tuple

        :param state_type: Type of states, either :class:`Function` or :class:`Constant`
        :type state_type: :class:`Function` or :class:`Constant`

    """

    def __init__(self, function_spaces, state_type=Function):

        # check that the state type is one of the required
        if state_type != Function and state_type != Constant:
            raise TypeError('state_type argument is not Constant or Function')

        # this is the hierarchy of function spaces
        self._fs_hierarchy = function_spaces

        if not isinstance(self._fs_hierarchy, (list, tuple)):
            raise TypeError('hierarchy of function spaces input is not a tuple or a list')

        # check that all levels of FS correspond to domains in the same hierarchy
        hierarchy = get_level(self._fs_hierarchy[0].mesh())[0]
        if hierarchy is None:
            raise ValueError('function spaces arent on meshes part of a hierarchy or the same hierarchy')
        else:
            self._mesh_hierarchy = hierarchy

        for fs in self._fs_hierarchy:
            assert get_level(fs.mesh())[0] == hierarchy

        # This is the actual hierarchy of ensembles
        self._hierarchy = []

        # Length of the hierarchy
        self.L = 0

        # Refinement factor
        self.M = self._mesh_hierarchy.refinements_per_level * 2

        # the min cell edge length in each level of ensemble hierarchy
        self.dxl = []

        # type of ensemble hierarchy input (default :class:`Function`)
        self._type = state_type

        # sample statistics
        if self._type == Function:
            self.MultilevelExpectation = Function(self._fs_hierarchy[-1])
            self.MultilevelExpectation.rename('MLMC Estimator')

        if self._type == Constant:
            self.MultilevelExpectation = Constant(0, domain=self._fs_hierarchy[-1].mesh())

        self.Mean = []  # : A list of the sample mean Functions of each of the different levels
        self.Variance = []  # : A list of the sample variance Functions of each of the different levels

        # running totals for sample statistics
        self.nl = []
        self.__sum_of_sq = []
        self.__sum = []

        super(EnsembleHierarchy, self).__init__()

    def __UpdateFunctionSpace(self, s, lvlc):

        """ Updates the :class:`FunctionSpace` of the appended State.state :class:`Function` to the
            finest level :class:`FunctionSpace` in the hierarchy

            :param S: The state
            :type S: :class:`State`.state

        """

        # Prolong both the state to this and put into a tuple
        u_1 = Function(self._fs_hierarchy[-1])
        u_2 = Function(self._fs_hierarchy[-1])
        if lvlc == len(self._fs_hierarchy) - 2:  # need this check to avoid error in prolonging to same level (finest)
            prolong(s[0], u_1)
            u_2.assign(s[1])
        else:
            prolong(s[0], u_1)
            prolong(s[1], u_2)
        Tuple_to_append = tuple([u_1, u_2])

        return Tuple_to_append

    def __UpdateDomain(self, s):

        """ Updates the :class:`Domain` of the appended State.state :class:`Constant` to the
            finest level :class:`Domain` in the hierarchy

            :param S: The state
            :type S: :class:`State`.state

        """

        # Prolong both the state to this and put into a tuple
        u_1 = Constant(s[0], domain=self._fs_hierarchy[-1].mesh())
        u_2 = Constant(s[1], domain=self._fs_hierarchy[-1].mesh())
        Tuple_to_append = tuple([u_1, u_2])

        return Tuple_to_append

    def AppendToEnsemble(self, S):

        """ Appends a state to a certain level of the ensemble hierarchy.

                :param S: The state
                :type S: :attr:`State'

        """

        if not isinstance(S.state, tuple):

            raise TypeError('input to append isnt a tuple')

        if not isinstance(S.state[0], self._type):

            raise TypeError('The state does not consist of correct types')

        if not isinstance(S.state[1], self._type):

            raise TypeError('The state does not consist of correct types')

        if get_level(S.state[0].ufl_domain())[1] is None:

            raise TypeError('state members arent part of a hierarchy')

        # Recover level of state (the coarser level)
        Level_to_append_to = S.levels[0]

        if Level_to_append_to < self.L:  # if there already exists that level in ensemble hierarchy

            # prolong the function space of state inputs to finest level of hierarchy
            # (:class:`Function` only!)
            if self._type == Function:
                T = self.__UpdateFunctionSpace(S.state, Level_to_append_to)

            if self._type == Constant:
                T = self.__UpdateDomain(S.state)

            # append state to required level
            self._hierarchy[Level_to_append_to].append(T)

            # update ensemble sizes
            self.nl[Level_to_append_to] += 1

            # update sums and squares
            if Level_to_append_to == 0:
                self.__sum[Level_to_append_to] = self.__sum[Level_to_append_to] + T[1]
                self.__sum_of_sq[Level_to_append_to] = (self.__sum_of_sq[Level_to_append_to] +
                                                        (T[1] ** 2))

            else:
                self.__sum[Level_to_append_to] = self.__sum[Level_to_append_to] + (T[1] - T[0])
                self.__sum_of_sq[Level_to_append_to] = (self.__sum_of_sq[Level_to_append_to] +
                                                        ((T[1] - T[0]) ** 2))

        if Level_to_append_to >= self.L:  # if the level doesn't exist in ensemble hierarchy

            # make sure previous level exists, otherwise make an empty ensemble on empty levels
            if Level_to_append_to - 1 > self.L:
                skipped_levels = Level_to_append_to - self.L

                for i in range(skipped_levels + 1):
                    self._hierarchy.append([])
                    self.nl.append(0)
                    self.dxl.append(0)

                    if self._type == Function:
                        self.Mean.append(Function(self._fs_hierarchy[-1]))
                        self.Mean[i].rename('Sample Mean')
                        self.Variance.append(Function(self._fs_hierarchy[-1]))
                        self.Variance[i].rename('Sample Variance')
                        self.__sum.append(Function(self._fs_hierarchy[-1]))
                        self.__sum_of_sq.append(Function(self._fs_hierarchy[-1]))

                    if self._type == Constant:
                        self.Mean.append(Constant(0,
                                                  domain=self._fs_hierarchy[-1].mesh()))
                        self.Variance.append(Constant(0,
                                                      domain=self._fs_hierarchy[-1].mesh()))
                        self.__sum.append(Constant(0,
                                                   domain=self._fs_hierarchy[-1].mesh()))
                        self.__sum_of_sq.append(Constant(0,
                                                         domain=self._fs_hierarchy[-1].mesh()))

            else:

                self._hierarchy.append([])
                self.nl.append(0)
                self.dxl.append(0)

                if self._type == Function:
                    self.Mean.append(Function(self._fs_hierarchy[-1]))
                    self.Mean[Level_to_append_to].rename('Sample Mean')
                    self.Variance.append(Function(self._fs_hierarchy[-1]))
                    self.Variance[Level_to_append_to].rename('Sample Variance')
                    self.__sum.append(Function(self._fs_hierarchy[-1]))
                    self.__sum_of_sq.append(Function(self._fs_hierarchy[-1]))

                if self._type == Constant:
                    self.Mean.append(Constant(0, domain=self._fs_hierarchy[-1].mesh()))
                    self.Variance.append(Constant(0, domain=self._fs_hierarchy[-1].mesh()))
                    self.__sum.append(Constant(0, domain=self._fs_hierarchy[-1].mesh()))
                    self.__sum_of_sq.append(Constant(0, domain=self._fs_hierarchy[-1].mesh()))

            # update length of ensemble hierarchy
            self.L = np.max([self.L, Level_to_append_to + 1])

            # prolong the function space of state inputs to finest level of hierarchy
            # (:class:`Function` only!)
            if self._type == Function:
                T = self.__UpdateFunctionSpace(S.state, Level_to_append_to)

            if self._type == Constant:
                T = self.__UpdateDomain(S.state)

            # append state to required level
            self._hierarchy[Level_to_append_to].append(T)

            # update ensemble sizes
            self.nl[Level_to_append_to] += 1

            # update sums and squares
            if Level_to_append_to == 0:
                self.__sum[Level_to_append_to] = self.__sum[Level_to_append_to] + T[1]
                self.__sum_of_sq[Level_to_append_to] = (self.__sum_of_sq[Level_to_append_to] +
                                                        (T[1] ** 2))

            else:
                self.__sum[Level_to_append_to] = self.__sum[Level_to_append_to] + (T[1] - T[0])
                self.__sum_of_sq[Level_to_append_to] = (self.__sum_of_sq[Level_to_append_to] +
                                                        ((T[1] - T[0]) ** 2))

            # compute min cell edge length for finer level of state
            self.dxl[Level_to_append_to] = MinDx(S.state[1].ufl_domain())

    def ClearEnsemble(self):
        """ Clears the ensemble hierarchy

        """

        for i in range(self.L):
            self._hierarchy[i] = []

    def UpdateStatistics(self, clear_ensemble=False):
        """ Updates statistics of ensemble hierarchy and has the option to clear the stored ensemble
            hierarchy afterwards

            :param clear_ensemble: Whether or not to clear the ensemble after updating stats (Fals)
            :type clear_ensemble: boolean

        """

        # update statistics on each level
        for i in range(self.L):
            self.Mean[i].assign(self.__sum[i] / self.nl[i])
            self.Variance[i].assign((self.__sum_of_sq[i] / self.nl[i]) - ((self.__sum[i] / self.nl[i]) ** 2))

        # update running multilevel monte carlo average
        if self._type == Function:
            mlmc_sum = Function(self._fs_hierarchy[-1])

        if self._type == Constant:
            mlmc_sum = Constant(0, domain=self._fs_hierarchy[-1].mesh())

        for i in range(self.L):
            mlmc_sum = mlmc_sum + self.Mean[i]

        self.MultilevelExpectation.assign(mlmc_sum)

        if clear_ensemble is True:

            # clear ensemble hierarchy
            self.ClearEnsemble()

    """ Iterative and Indexing functions """

    def __len__(self):
        """ Return the length of the ensemble hierarchy """
        return self.L

    def __iter__(self):
        """ Iterate over the ensembles in the hierarchy (from
        coarse to fine). """
        for en in self._hierarchy:
            yield en

    def __getitem__(self, idx):
        """ Return an ensemble in the hierarchy

            :arg idx: The index of the ensemble to return

        """
        return self._hierarchy[idx]
