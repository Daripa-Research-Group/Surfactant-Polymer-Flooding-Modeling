"""
This python script contains the class definition for polymers for the surfactant-flooding model

"""


class Polymer:
    def __init__(
        self,
        name,
        initial_concentration,
        e_coeff,
        n_coeff,
        viscosity_scalar=None,
        viscosity_matrix=None,
        vec_concentration=None,
    ):
        """
        Initializes a instance of the polymer class

        :param name: Name of the polymer
        :type name: enum 'PolymerList'

        :param initial_concentration: Initial concentration of polymer in the injected solution (scalar variable)
        :type concentration: float

        :param e_coeff: The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type e_coeff: List<int>

        :param n_coeff:  The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type n_coeff: List<int>

        :param viscosity_scalar: scalar quantitiey of the polymer viscosity
        :type viscosity_matrix: float, None

        :param viscosity_matrix: viscosity matrix of the polymer
        :type viscosity_matrix: np.array, None

        :param vec_concentration: vector representation of polymer concentration within resevoir
        :type vec_concentration: np.array, None
        """

        self.name = name
        self.initial_concentration = initial_concentration
        self.vec_concentration = vec_concentration
        self.viscosity_matrix = viscosity_matrix
        self.viscosity_scalar = viscosity_scalar
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff
