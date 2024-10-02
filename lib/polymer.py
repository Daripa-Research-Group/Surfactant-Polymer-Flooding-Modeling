"""
This python script contains the class definition for polymers for the surfactant-flooding model

"""

class Polymer:
    def __init__(self, name, concentration, e_coeff, n_coeff):
        """
        Initializes a instance of the polymer class

        :param name: Name of the polymer
        :type name: str

        :param concentration: concentration of polymer for surfactant-polymer flooding simulation in wppm
        :type concentration: float

        :param e_coeff: The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type e_coeff: List<int>

        :param n_coeff:  The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type n_coeff: List<int>
        """

        self.name = name
        self.concentration = concentration
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff

