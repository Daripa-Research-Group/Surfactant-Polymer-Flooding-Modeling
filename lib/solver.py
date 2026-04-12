"""
Will hold the classes relevant for solving the elliptic equations and the transport equations

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np

class TransportEquationSolver():
    """
    Will hold all the relevant functions and parameters necessary for solving the transport equations
    """
    def __init__(
            self,
            simulation_scenario,
            water_saturation,
            surfactant_concentration,
            polymer_concentration,
            aqueous_mobility,
            aqueous_viscosity,
            oleic_mobility,
            oleic_viscosity,
            
    ):
        """
        constructor for ``TransportEquationSolver``
        """
        pass

    def execute(self):
        pass

    def saturation_matrix_preprocessing(self):
        pass

    def surfactant_concentration_matrix_preprocessing(self):
        pass

    def polymer_concentration_matrix_preprocessing(self):
        pass
        

# class EllipticEquationSolver():
#     """
#     Will hold the relevant functions and parameters for solving the global pressure (u) and velocity (v)
#
#     (Not implemented for v1)
#     """
#     pass
