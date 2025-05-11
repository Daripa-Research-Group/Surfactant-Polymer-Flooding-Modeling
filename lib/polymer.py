"""
This python script contains the class definition for polymers for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
from enumerations import ModelType, PolymerList, SimulationConstants
from lib.Exceptions import SimulationCalcInputException
from lib.para import Box

class Polymer:
    """
    Class definition for the polymer objects and their calculations
    """
    def __init__(
            self,
            name : PolymerList,
            e_coeff : list,
            n_coeff : list,
            rho : float,
            concentration_scalar : float,
            viscosity_scalar: float | None,
            viscosity_matrix: np.ndarray | None,
            concentration_matrix: np.ndarray | None,
            shear_rate: np.ndarray | None
            ):
        """
        Initializes a instance of the polymer class

        :param name: Name of the polymer
        :type name: enum 'PolymerList'

        :param e_coeff: The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type e_coeff: list<int>

        :param n_coeff:  The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        :type n_coeff: list<int>

        :param rho: Density of polymer
        :type rho: float

        :param concentration_scalar: Scalar quantity of concentration. When initializing, this param will equal the initial polymer concentration. 
        :type concentration_scalar: float
        
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
        self.concetration_scalar = concentration_scalar
        self.init_concentration_matrix = concentration_scalar * np.ones((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        self.concentration_matrix = concentration_matrix if(concentration_matrix is not None) else concentration_scalar * np.ones((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))

        
        #Properties related to the viscosity
        self.viscosity_matrix = viscosity_matrix if(viscosity_matrix is not None) else np.zeros((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        self.viscosity_scalar = viscosity_scalar if(viscosity_scalar is not None) else 0

        #Values required to formulate the numerical powerlaw function for viscosity calculations
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff
        
        #Polymer Density
        self.rho = rho 

        #shear rate matrix (needed when running 'shear thinning' model version)
        self.shear_rate = shear_rate if(shear_rate is not None) else np.zeros((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value)) 

    def initialize(self):
        """
        Will initialize the viscosity and concentration matrices

        :return: a list of the initialized concentration, viscosity, and shear_rate matrices
        :rtype: List
        """
        pass


    def compute_viscosity(
            self, 
            grid : tuple, 
            u : np.ndarray, 
            v : np.ndarray, 
            model_type : ModelType,
            aqueous_viscosity : np.ndarray | None, 
            ):
        """
        Compute polymer viscosity.
        This function is derived from 'compvis()' in the original MATLAB code (in file compvis.m).

        :param grid: The FEM grid used for simulation calculations (x and y variables from the MATLAB code)
        :type grid: tuple[NDArray[Any], ...]

        :param u: Matrix related to the global pressure
        :type u: np.ndarray

        :param v: Matrix related to the velocity matrix
        :type v: np.ndarray

        :param model_type: Will state whether the model will include polymer shear thinning or not
        :type model_type: enum 'ModelType'

        :param aqueous_viscosity: Aqueous viscosity matrix (will come from the 'Water' class). Only needed when shear thinning OFF
        :type aqueous_viscosity: np.ndarray, None
        
        :return: the viscosity_matrix (index 0) & shear_rate matrix (index 1) for the polymer within the grid
        :rtype: list
        """
        # x and y components from meshgrid
        x = grid[0]
        y = grid[1]

        if(self.concentration_matrix is None or self.init_concentration_matrix is None or self.concetration_scalar is None):
            raise SimulationCalcInputException("SimulationInputException: Polymer concentration matrix and/or scalar concentration value not initialized...") 

        #if model_type is NO SHEAR THINNING:
        if(model_type == ModelType.No_Shear_Thinning.value):
            if(aqueous_viscosity is None):
                raise SimulationCalcInputException("SimulationInputException: Aqueous viscosity matrix required but not provided. Please try again.")
            ## the scalar viscosity is equal to the max within the aqueous viscosity matrix
            self.viscosity_scalar = np.max(aqueous_viscosity[0, :])
            self.viscosity_matrix = self.viscosity_scalar*np.ones((SimulationConstants.Grid_Size.value,SimulationConstants.Grid_Size.value))
        #Model Type is 'Sourav Implementation':
        elif(model_type == ModelType.Sourav_Implementation.value):
            #TODO: Will keep empty until properly understood how to implement
            pass
        #if polymer shear thinning is ON:
        elif(model_type == ModelType.Shear_Thinning_On.value):
            if(aqueous_viscosity is not None):
                raise SimulationCalcInputException("SimulationInputException: Aqueous viscosity reliant on changing polymer viscosity. Update will be done within 'Water' Class")
            if(self.shear_rate is None or self.viscosity_matrix is None):
                raise SimulationCalcInputException("SimulationInputException: Either shear_matrix or viscosity_matrix are not initialized")
            # Getting water density and viscosity (Note: polymer density a property of class)
            rho_water = SimulationConstants.Water_Density.value
            viscosity_water = SimulationConstants.Water_Viscosity.value

            # Formulating the numerically derived power law equation
            w1 = self.rho*self.concentration_matrix
            w2 = rho_water*(1-self.concentration_matrix)
            wppm = (w1/(w1+w2))*(10**6)
            w1_0 = self.rho*self.init_concentration_matrix #from the variable w10 in MATLAB code
            w2_0 = rho_water*(1-self.init_concentration_matrix) #from the variable w20 in MATLAB code
            wppm_0 = (w1_0/(w1_0+w2_0))*(10**6) #from the wppm0 variable in MATLAB code
            
            ## Determining the epsilon and n coefficients for the power law equation 
            epsilon_0 = self.e_coeff[0]*(wppm_0**self.e_coeff[1])
            n_0 = np.min(self.n_coeff[0]*(wppm_0**self.n_coeff[1]))
            epsilon = self.e_coeff[0]*(wppm**self.e_coeff[1])
            n = np.min(self.n_coeff[0]*(wppm**self.n_coeff[1]))

            row = np.size(self.concentration_matrix, 0)
            col = np.size(self.concentration_matrix, 1)

            # Compute divergence terms
            a1 = np.gradient(v, axis=0)
            a2 = np.gradient(u, axis=1)
            a3 = np.gradient(u, axis=0)
            a4 = np.gradient(v, axis=1)

            pi_D = np.abs(-0.25 * (a1 + a2) ** 2 + a3 * a4)

            for i in range(row):
                for j in range(col):
                    if self.concentration_matrix[i, j] > 0:
                        self.shear_rate[i, j] = 2 * np.sqrt(pi_D[i, j])
                        if self.shear_rate[i, j] != 0:
                            self.viscosity_matrix[i, j] = epsilon_0[i, j] * (self.shear_rate[i, j] ** (n_0[i, j] - 1))
                            self.viscosity_matrix[i, j] = np.clip(self.viscosity_matrix[i, j], viscosity_water, 100)


        return [self.viscosity_matrix, self.shear_rate]

            


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
