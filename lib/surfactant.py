"""
This python script contains the class definition for surfactants for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""


from types import LambdaType
from lib.enumerations import SimulationConstants, SurfactantList
from lib.para import Box
import numpy as np

class Surfactant:
    """
    Contains property and calculations related to the Surfactant
    """
    def __init__(
        self,
        name : SurfactantList,
        initial_concentration : float,
        IFT_equation : LambdaType,
        derivative_IFT_equation : LambdaType,
        phi : np.ndarray,
        concentration_matrix : np.ndarray | None,
    ):
        """
        Creates instance of Surfactant class

        :param name: Name of the surfactant
        :type name: enum 'SurfactantList'

        :param concentration: Initial concentration in wppm of surfactant (scalar quantity)
        :type concentration: float

        :param IFT_conc_equ: expression that relates surfactant concentration to interfacial tension b/t oil and water
        :type IFT_conc_equ: lambda

        :param derivative_IFT_conc_equ: Deriviative of the equation relating IFT to surfactant concentration
        :type derivative_IFT_conc_equ: lambda

        :param vec_concentration: vector representation of surfactant concentration in resevoir
        :type vec_concentration: np.array, None
        """
        self.name = name
        self.concentration = initial_concentration
        self.concentration_matrix = concentration_matrix if (concentration_matrix is not None) else initial_concentration * np.zeros((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        self.IFT_conc_equ = IFT_equation
        self.derivative_IFT_conc_equ = derivative_IFT_equation
        self.is_surfactant = True if (initial_concentration > 0) else False
        self.phi = phi
    


    def initialize(self):
        """
        This function will initialize the surfactant object

        :return: Surfactant object
        :rtype: Surfactant
        """
        pass

    def compute_concentration(
            self,
            grid: tuple,
            mesh: Box,
            u: np.ndarray, 
            v: np.ndarray,
            dt: float,
            initial_water_saturation: float,
            water_saturation_matrix: np.ndarray,
            xmod : np.ndarray,
            ymod : np.ndarray,
            ):
        """
        This function will update the surfactant concentration matrix

        This function is derived from the section of the 'nmmoc_surf_mod_neumann'
        related to the polymer concentration matrix
        
        :param grid: The FEM grid used for simulation calculations (x and y variables from the MATLAB code)
        :type grid: tuple[NDArray[Any], ...]

        :param mesh: 'Box' object containing information for the FEM grid
        :type mesh: Box

        :param u: Matrix related to the global pressure
        :type u: np.ndarray

        :param v: Matrix related to the velocity matrix
        :type v: np.ndarray

        :param dt: time-step
        :type dt: float

        :param initial_water_saturation: the scalar quantity of the initial water saturation in sim
        :type initial_water_saturation: float

        :param water_saturation_matrix: the updated water saturation matrix
        :type water_saturation_matrix: np.ndarray

        :param xmod: x-dimension coordinate points for formulating the 'Cmod' matrix
        :type xmod: np.ndarray

        :param ymod: y-dimension coordinate points for formulating the 'Cmod' matrix
        :type ymod: np.ndarray

        :return: Surfactant concentration matrix
        :rtype: np.ndarray

        """
        pass
