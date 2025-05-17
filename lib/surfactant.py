"""
This python script contains the class definition for surfactants for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""


from types import LambdaType
from lib.enumerations import SimulationConstants, SurfactantList
from Exceptions import SimulationCalcInputException
from lib.para import Box
import numpy as np
import scipy as sp
from scipy.sparse.linalg import bicgstab

class Surfactant:
    """
    Contains property and calculations related to the surfactant object
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

        :param concentration_matrix: vector representation of surfactant concentration in resevoir
        :type concentration_matrix: np.array, None
        """
        self.name = name
        self.concentration = initial_concentration
        self.concentration_matrix = concentration_matrix 
        self.IFT_conc_equ = IFT_equation
        self.derivative_IFT_conc_equ = derivative_IFT_equation
        self.is_surfactant = True if (initial_concentration > 0) else False
        self.phi = phi

        #initializing the surfactant object
        self.initialize()
    


    def initialize(
            self
            ):
        """
        This function will initialize the surfactant object

        :return: Surfactant object
        :rtype: Surfactant
        """
        if(self.concentration_matrix is None):
            if(self.phi is None):
                raise SimulationCalcInputException("SimulationInputException: phi value not initalized...")
            D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
            self.concentration_matrix = (~D) * self.concentration
        return self

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
            param_coeff : dict
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

        :param param_coeff: contains the updated constants (calculated in the simulation class) for updating the surfactant concentration matrix
        :type param_coeff: dict

        :return: Surfactant concentration matrix
        :rtype: np.ndarray

        """
        # initializing constants
        if(self.concentration_matrix is None):
            raise SimulationCalcInputException("SimulationInputException: Surfactant concentration matrix not initialized. Please initizlize before running this function.")
        g1 = initial_water_saturation
        g3 = initial_water_saturation * self.concentration
        x = grid[0]
        y = grid[1]
        m = mesh.m
        n = mesh.n
        dx = mesh.dx
        dy = mesh.dy
        dt_array = dt*np.ones((SimulationConstants.Grid_Size.value, SimulationConstants.Grid_Size.value))
        Qnew = water_saturation_matrix

        x1d = x[0, :]
        y1d = y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)   
        
        # reorder surfactant.concentration_matrix if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            self.concentration_matrix = self.concentration_matrix[:, x_sort_idx]  # Sort columns of surfactant.concentration_matrix
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            self.concentration_matrix = self.concentration_matrix[y_sort_idx, :]  # Sort rows of surfactant.concentration_matrix
            
        interp = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d), self.concentration_matrix ,method='linear', bounds_error=False, fill_value=None
        )
        Gmod = interp((xmod, ymod))

        # Updating coefficients using interpolated surfactant concentration
        sigma_mod = self.IFT_conc_equ(Gmod)
        sigma_g_mod = self.derivative_IFT_conc_equ(Gmod)
        lambda_a = param_coeff['lambda_a']
        lambda_total = param_coeff['lambda_total']
        D = param_coeff['D']
        pc_g = param_coeff['pc_g']

        # intermediate parameters for code:
        F = D * pc_g / Qnew
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
                        if idx == 1:
                            if i == 0:
                                DD[i] = g3 / Qnew[cnt][i] + Gmod[cnt][i] / dt_array[cnt][i]
                                CC[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                    + g1 / Qnew[cnt][i]
                                )
                                BB[j][i + 1] = 2 * F[cnt][i] / (dx**2)
                            elif i == m - 1:  # Bottom right point
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                CC[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i - 1] = 2 * F[cnt][i] / (dx**2)
                            else:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                CC[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i - 1] = F[cnt][i] / (dx**2)
                                BB[j][i + 1] = 2 * F[cnt][i] / (dx**2)
                        elif idx == (m) * (n - 1) + 1:
                            if i == 0:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                    + g1 / Qnew[cnt][i]
                                )
                                BB[j][i + 1] = 2 * F[cnt][i] / (dx**2)
                            elif i == m - 1:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                    - ((g1 * lambda_a[cnt][i]) / [lambda_total[cnt][i]])
                                    / Qnew[cnt][i]
                                    + ((g3 * lambda_a[cnt][i]) / [lambda_total[cnt][i]])
                                    / (Qnew[cnt][i] * self.concentration)
                                )
                                BB[j][i - 1] = 2 * F[cnt][i] / (dx**2)
                            else:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = 2 * F[cnt][i] / (dy**2)
                                BB[j][i + 1] = F[cnt][i] / (dx**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i - 1] = F[cnt][i] / (dx**2)
                        else:
                            if i == 0:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i + 1] = 2 * F[cnt][i] / (dx**2)
                                CC[j][i] = F[cnt][i] / (dy**2)
                            elif i == m - 1:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i - 1] = 2 * F[cnt][i] / (dx**2)
                                CC[j][i] = F[cnt][i] / (dy**2)
                            else:
                                DD[i] = Gmod[cnt][i] / dt_array[cnt][i]
                                AA[j][i] = F[cnt][i] / (dy**2)
                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - ((2 / (dx**2)) + (2 / (dy**2))) * F[cnt][i]
                                )
                                BB[j][i - 1] = 2 * F[cnt][i] / (dx**2)
                                BB[j][i + 1] = 2 * F[cnt][i] / (dx**2)
                                CC[j][i] = F[cnt][i] / (dy**2)
            if cnt == 0:
                AAA[:n, : 2 * m] = np.hstack([BB, CC])
            elif cnt == n - 1:
                AAA[(m - 1) * n : m * n, (n - 2) * m : n * m] = np.hstack([AA, BB])
            else:
                AAA[cnt * n : (cnt + 1) * n, (cnt - 1) * m : (cnt + 2) * m] = np.hstack(
                    [AA, BB, CC]
                )

            DDD[cnt * m : (cnt + 1) * m] = DD
            idx += m

        Gnew_flat, info = bicgstab(AAA, DDD, rtol=10 ** (-10), maxiter=600)
        Gnew = Gnew_flat.reshape(m, n)
        
        self.concentration_matrix = Gnew

        return self.concentration_matrix
