"""
Will hold the classes relevant for solving the elliptic equations and the transport equations

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

from surfactant import Surfactant
from water import Water
from polymer import Polymer
from grid import Grid
from enumerations import *
import numpy as np

class TransportEquationSolver():
    """
    Will hold all the relevant functions and parameters necessary for solving the transport equations
    """
    def __init__(
            self,
            grid : Grid,
            water : Water,
            surfactant : Surfactant,
            polymer : Polymer,
            pressure : np.ndarray,
            velocity : np.ndarray,
            source_flow_magnitude: float
    ):
        """
        constructor for ``TransportEquationSolver``
        """
        pass

    #Dependent Properties
    ## Flows
    _total_flow = None
    @property
    def total_flow(self):
        return self._total_flow
    @total_flow.setter
    def total_flow(self, value):
        self._total_flow = value
    @property
    def polymer_flow(self):
        return self._total_flow*self._water.concetration_scalar
    @property
    def surfactant_flow(self):
        return self._total_flow*self._surfactant.concentration
    ## Grid properties
    @property
    def m(self):
        return self._grid.m
    @property
    def n(self):
        return self._grid.n
    @property
    def dx(self):
        return self._grid.dx
    @property
    def dy(self):
        return self._grid.dy
    @property
    def x(self):
        return self._grid.x
    @property
    def y(self):
        return self._grid.y
    @property
    def dt_array(self):
        pass
    ## Porosity
    @property
    def porosity(self):
        return 1
    ## Capillary pressure parameters
    @property
    def omega(self):
        return [0.1, 0.4] #[Ω_1, Ω_2]
    ## Capillary Number
    @property
    def critical_capillary_aqueous_initial(self): #Nca0
        return SimulationConstants.Aqueous_Phase_Critical_Capillary_Num.value 
    @property
    def critical_capillary_oleic_initial(self): #Nco0
        return SimulationConstants.Oleic_Phase_Critical_Capillary_Num.value
    
    # Functions for parameter definitions
    ## Residual saturations
    def _compute_residual_saturations(self):
        pass

    def _normalized_residual_saturations(self):
        pass

    def _derivative_residual_saturations(self): #FIXME: Will need to update function to work with autodiff (v2.0)
        pass
    
    ## Mobilities
    def _compute_lambda_a(self):
        pass

    def _compute_lambda_o(self):
        pass

    ## Capillary Number (ratio of viscous forces to surface tension forces)
    def _aqueous_capillary_number(self):
        pass
    def _oleic_capillary_number(self):

    






    def execute(self):
        pass

    def _initialize_matrices(self):
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
