"""
This python script contains the class definition for polymers for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np

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
        shear_rate = None
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

        #properties related to the concentration
        self.concetration_scalar = initial_concentration
        self.concentration_matrix = vec_concentration
        
        #Properties related to the viscosity
        self.viscosity_matrix = viscosity_matrix
        self.viscosity_scalar = viscosity_scalar

        #Values required to formulate the numerical powerlaw function for viscosity calculations
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff

        #shear rate matrix (needed when running 'shear thinning' model version)
        self.shear_rate = shear_rate

    def initialize(self):
        """
        Will initialize the viscosity and concentration matrices

        :return: a list of the initialized concentration and viscosity matrices
        :rtype: List
        """
        pass


    def compute_viscosity(self, grid, model_type):
        """
        Compute polymer viscosity.
        This function is derived from 'compvis()' in the original MATLAB code.

        :param grid: The FEM grid used for simulation calculations
        :type grid: np.ndarray

        :param model_type: Will state whether the model will include polymer shear thinning or not
        :type model_type: enum 'ModelType'

        :return: the viscosity_matrix for the polymer within the grid
        :rtype: np.ndarrray
        """
        pass


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

