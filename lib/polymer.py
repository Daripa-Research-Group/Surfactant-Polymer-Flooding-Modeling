"""
This python script contains the class definition for polymers for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
import scipy as sp
from scipy.sparse.linalg import bicgstab
from enumerations import ModelType, PolymerList, SimulationConstants
from Exceptions import SimulationCalcInputException
from grid import Grid

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
            phi: np.ndarray,
            viscosity_scalar: float | None = None,
            viscosity_matrix: np.ndarray | None = None,
            concentration_matrix: np.ndarray | None = None,
            shear_rate: np.ndarray | None = None
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

        :param phi: arrray used to initialize the concentration matrix (represents porosity of the resevoir)
        :type: np.ndarray
        
        :param viscosity_scalar: scalar quantity of the polymer viscosity
        :type viscosity_scalar: float, None

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
        self.concentration_matrix = concentration_matrix 
        
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

        #util param for initialization
        self.phi = phi # Will need to be created in the simulation class
        
    def initialize(
            self,
            grid_shape : tuple
            ):
        """
        Will initialize the viscosity, shear_rate, and concentration matrices

        :param grid_shape: contain the shape of the grid
        :type grid_shape: tuple

        :return: Initalized Polymer Object
        :rtype: Polymer
        """
        n = grid_shape[0]
        m = grid_shape[1]
        D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
        if(self.concentration_matrix is None):
            self.concentration_matrix = (~D)*self.concetration_scalar

        if(self.shear_rate is None):
            self.shear_rate = np.zeros((n+1,m+1))

        if(self.viscosity_scalar is None):
            beta1 = 15000 #constant that came from the MATLAB code
            self.viscosity_scalar = SimulationConstants.Water_Viscosity.value*(1+beta1+self.concetration_scalar)
        
        if(self.viscosity_matrix is None):
            self.viscosity_matrix = self.viscosity_scalar * np.ones((n+1,m+1))

        return self


    def compute_viscosity(
            self, 
            grid : tuple, 
            u : np.ndarray, 
            v : np.ndarray, 
            model_type : ModelType,
            aqueous_viscosity : np.ndarray | None = None, 
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

            


    def compute_concentration(
            self, 
            grid: tuple,
            mesh: Grid,
            u: np.ndarray, 
            v: np.ndarray,
            dt: float,
            f_lambda: np.ndarray,
            initial_water_saturation: float,
            water_saturation_matrix: np.ndarray,
            xmod : np.ndarray,
            ymod : np.ndarray,
            ):
        """
        Update the polymer concentration matrix and the shear rate tensor

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

        :param f_lambda: fraction of lambda_aqueous / lambda_total
        :type f_lambda: np.ndarray

        :param initial_water_saturation: the scalar quantity of the initial water saturation in sim
        :type initial_water_saturation: float

        :param water_saturation_matrix: the updated water saturation matrix
        :type water_saturation_matrix: np.ndarray

        :param xmod: x-dimension coordinate points for formulating the 'Cmod' matrix
        :type xmod: np.ndarray

        :param ymod: y-dimension coordinate points for formulating the 'Cmod' matrix
        :type ymod: np.ndarray

        :return: Polymer concentration matrix
        :rtype: np.ndarray
        """
        #initializing variables:
        if(self.concentration_matrix is None):
            raise SimulationCalcInputException("SimulationInputException: Polymer concentration matrix not initialized. Please initalize before running method.")
        x = grid[0]
        y = grid[1]
        m = mesh.m
        n = mesh.n
        dt_array = dt*np.ones((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        Qnew = water_saturation_matrix

        g1 = initial_water_saturation
        g2 = initial_water_saturation*self.concetration_scalar

        # Determining 'Cmod'
        x1d = x[0, :]
        y1d = y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)   
        
        # reorder vec_concentration if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            self.concentration_matrix = self.concentration_matrix[:, x_sort_idx]  # Sort columns of vec_concentration
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            self.concentration_matrix = self.concentration_matrix[y_sort_idx, :]  # Sort rows of vec_concentration
            
        interp = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d), self.concentration_matrix, method='linear', bounds_error=False, fill_value=None
        )
        Cmod = interp((xmod, ymod))
    
        # Using 'Cmod' and 'Qnew' to update the polymer concentration matrix
        idx = 1
        AAA = np.zeros((n * m, n * m))
        DDD = np.zeros((n * m, 1))

        while idx <= (m) * (n - 1) + 1:
            cnt = (idx - 1) // m  # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n, m))
            AA = BB
            CC = BB
            DD = np.zeros((m, 1))
            for i in range(m - 1):
                for j in range(n - 1):
                    if j == i:
                        if idx == 1:  # lowermost row of grid
                            if i == 1:  # leftmost point (source)
                                DD[i] = g2 / Qnew[cnt][i] + Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = 1 / dt_array[cnt][i] + g1 / Qnew[cnt][i]
                            else:
                                DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = 1 / dt_array[cnt][i]
                        elif idx == (m) * (n - 1) + 1:
                            if i == m - 1:
                                DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = (
                                    1 / dt_array[cnt][i] - g1 * f_lambda[cnt][i] / Qnew[cnt][i]
                                )
                            else:
                                DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = 1 / dt_array[cnt][i]
                        else:
                            DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                            BB[j][i] = 1 / dt_array[cnt][i]

            if cnt == 0:
                AAA[0:n, 0 : 2 * m] = np.hstack([BB, CC])
            elif cnt == n - 1:
                AAA[(m - 1) * n : m * n, (n - 2) * m : n * m] = np.hstack([AA, BB])
            else:
                AAA[cnt * n : (cnt + 1) * n, (cnt - 1) * m : (cnt + 2) * m] = np.hstack(
                    [AA, BB, CC]
                )

            DDD[cnt * m : (cnt + 1) * m] = DD

            idx += m

        Cnew_flat, info = bicgstab(AAA, DDD, rtol=10 ** (-10), maxiter=600)
        Cnew = Cnew_flat.reshape(m, n)
        self.concentration_matrix = Cnew

        return self.concentration_matrix
