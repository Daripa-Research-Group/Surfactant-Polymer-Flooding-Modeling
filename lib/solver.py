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
            water,
            surfactant,
            polymer,
            pressure,
            velocity,
            source_flow_magnitude
    ):
        """
        constructor for ``TransportEquationSolver``
        """
        pass

    def execute(self):
        pass

    def _initialize(self):
        """
        initialize relevant matrices needed for computation
        """
        pass

    def _main_loop_computation(self):
        """
        primary loop that will run for computations
        """
        pass

    def _saturation_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the saturation matrix
        """
        pass

    def _surfactant_concentration_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the surfactant concentration matrix
        """
        pass

    def _polymer_concentration_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the polymer concentration matrix
        """
        pass

    def _bottom_grid_calculations(self):
        """
        will conduct calculations related to the bottom of the grid
        """
        pass

    def _top_grid_calculations(self):
        """
        will conduct calculatons related to the top of the grid
        """
        pass

    def _interior_grid_calculations(self):
        """
        will conduct calculations related to the interior of the grid
        """
        pass


# class EllipticEquationSolver():
#     """
#     Will hold the relevant functions and parameters for solving the global pressure (u) and velocity (v)
#
#     (Not implemented for v1)
#     """
#     pass
