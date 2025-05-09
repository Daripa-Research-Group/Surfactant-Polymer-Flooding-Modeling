"""
This python script contains the class definition for polymers for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
from enumerations import ModelType, PolymerList, SimulationConstants
from lib.Exceptions import SimulationCalcInputException

class Polymer:
    """
    Class definition for the polymer objects and the calculation that they can perform
    """
    def __init__(
        self,
        name,
        initial_concentration,
        e_coeff,
        n_coeff,
        rho,
        viscosity_scalar=None,
        viscosity_matrix=None,
        vec_concentration=None,
        shear_rate=None
    ):
        """
        Initializes a instance of the polymer class

        :param name: Name of the polymer
        :type name: enum 'PolymerList'

        :param concentration_scalar: holds the scalar quantity of the polymer concentration within the resevoir 
        :type concentration_scalar: float

        :param e_coeff: The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type e_coeff: List<int>

        :param n_coeff:  The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type n_coeff: List<int>

        :param rho: Density of polymer
        :type rho: float

        :param viscosity_scalar: scalar quantitiey of the polymer viscosity
        :type viscosity_matrix: float, None

        :param viscosity_matrix: viscosity matrix of the polymer
        :type viscosity_matrix: np.ndarray, None

        :param concentration_matrix: matrix representation of polymer concentration within resevoir
        :type concentration_matrix: np.ndarray, None

        :param shear_rate: Matrix that will hold the shear rate (the change in velocity normal to the direction of flow)
        :type shear_rate: np.ndarray, None
        """
        
        #PolymerList object
        self.name = name

        #properties related to the concentration (scalar concentration, matrix version of initial concentration, and current concentration matrix)
        self.concetration_scalar = initial_concentration
        self.init_concentration_matrix = initial_concentration * np.ones((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        self.concentration_matrix = vec_concentration
        
        #Properties related to the viscosity
        self.viscosity_matrix = viscosity_matrix
        self.viscosity_scalar = viscosity_scalar

        #Values required to formulate the numerical powerlaw function for viscosity calculations
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff
        
        #Polymer Density
        self.rho = rho 

        #shear rate matrix (needed when running 'shear thinning' model version)
        self.shear_rate = shear_rate

    def initialize(self):
        """
        Will initialize the viscosity and concentration matrices

        :return: a list of the initialized concentration, viscosity, and shear_rate matrices
        :rtype: List
        """
        pass


    def compute_viscosity(self, grid, u, v, aqueous_viscosity, model_type):
        """
        Compute polymer viscosity.
        This function is derived from 'compvis()' in the original MATLAB code (in file compvis.m).

        :param grid: The FEM grid used for simulation calculations (x and y variables from the MATLAB code)
        :type grid: np.ndarray

        :param u: Matrix related to the global pressure
        :type u: np.ndarray

        :param v: Matrix related to the velocity matrix
        :type v: np.ndarray

        :param aqueous_viscosity: Aqueous viscosity matrix
        :type aqueous_viscosity: np.ndarray

        :param model_type: Will state whether the model will include polymer shear thinning or not
        :type model_type: enum 'ModelType'

        :return: the viscosity_matrix & shear_rate matrix for the polymer within the grid
        :rtype: List
        """
        # x and y components from meshgrid
        x = grid[0]
        y = grid[1]

        if(self.concentration_matrix is None or self.init_concentration_matrix is None or self.concetration_scalar is None):
            raise SimulationCalcInputException("SimulationInputException: Polymer concentration matrix and/or scalar concentration value not initialized...") 

        #if model_type is NO SHEAR THINNING:
        if(model_type == ModelType.No_Shear_Thinning.value):
            ## the scalar viscosity is equal to the max within the aqueous viscosity matrix
            self.viscosity_scalar = np.max(aqueous_viscosity[0, :])
            self.viscosity_matrix = self.viscosity_scalar*np.ones((SimulationConstants.Grid_Size.value,SimulationConstants.Grid_Size.value))
        #Model Type is 'Sourav Implementation':
        elif(model_type == ModelType.Sourav_Implementation.value):
            #TODO: Will keep empty until properly understood how to implement
            pass
        #if polymer shear thinning is ON:
        elif(model_type == ModelType.Shear_Thinning_On.value):
            # Getting water density (Note: polymer density a property of class)
            rho_water = SimulationConstants.Water_Density.value

            # Formulating the numerically derived power law equation
            w1 = self.rho*self.concentration_matrix
            w2 = rho_water*(1-self.concentration_matrix)
            wppm = (w1/(w1+w2))*(10**6)
            w1_0 = self.rho*self.init_concentration_matrix #from the variable w10 in MATLAB code
            w2_0 = rho_water*(1-self.init_concentration_matrix) #from the variable w20 in MATLAB code
            wppm_0 = (w1_0/(w1_0+w2_0))*(10**6) #from the wppm0 variable in MATLAB code

            epsilon_0 = self.e_coeff[0]*(wppm_0**self.e_coeff[1])
            n_0 = np.min(self.n_coeff[0]*(wppm_0**self.n_coeff[1]))
            epsilon = self.e_coeff[0]*(wppm**self.e_coeff[1])
            n = np.min(self.n_coeff[0]*(wppm**self.n_coeff[1]))

            


    def compute_concentration(self, grid, u, v):
        """
        Update the polymer concentration matrix and the shear rate tensor

        This function is derived from the section of the 'nmmoc_surf_mod_neumann'
        related to the polymer concentration matrix

        :param grid: The FEM grid used for simulation calculations
        :type grid: np.ndarray

        :param u: matrix that holds the global pressure
        :type u: np.ndarray
        """
        pass
