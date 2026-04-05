"""
This python script contains the class definition for surfactants for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

from types import LambdaType
import numpy as np
from scipy.sparse.linalg import bicgstab
from enumerations import SurfactantList
from grid import Grid
from Exceptions import SimulationCalcInputException


class Surfactant:
    """
    Contains property and calculations related to the surfactant object
    """

    def __init__(
        self,
        name: SurfactantList,
        initial_concentration: float,
        phi: np.ndarray,
        IFT_equation: LambdaType | None = None,
        derivative_IFT_equation: (
            LambdaType | None
        ) = None,  # FIXME: can remove once implemented `autodiff` capabilities
        concentration_matrix: np.ndarray | None = None,
    ):
        """
        Creates instance of Surfactant class
        """
        self.name = name
        self.concentration = initial_concentration
        self.concentration_matrix = concentration_matrix
        self.IFT_conc_equ = IFT_equation
        self.derivative_IFT_conc_equ = derivative_IFT_equation  # FIXME: need to adjust when implementing 'autodiff'
        self.is_surfactant = True if (initial_concentration > 0) else False
        self.phi = phi

    #CLASS PROPERTIES
    @property
    def eval_IFT(self):
        """
        evaluate IFT at a given surfactant concentration_matrix
        """
        assert self.IFT_conc_equ is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownIFTEquation"
        )
        return self.IFT_conc_equ(self.concentration_matrix)

    @property
    def eval_dIFT_dGamma(self):  # FIXME: Need to adjust when implementing 'autodiff'
        """
        evaluate the dσ/dΓ at a particular surfactant concentration matrix
        """
        assert self.derivative_IFT_conc_equ is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownDerivativeIFTEquation"
        )
        return self.derivative_IFT_conc_equ(self.concentration_matrix)

    _name = None
    @property
    def name(self):
        """
        name (enum 'SurfactantList'): Name of the surfactant
        """
        return self._name
    @name.setter
    def name(self, value):
        self._name = value
    
    _concentration = None
    @property
    def concentration(self):
        """
        concentration (float): Initial concentration of surfactant (scalar quantity)
        """
        return self._concentration
    @concentration.setter
    def concentration(self, value):
        self._concentration = value

    _concentration_matrix = None
    @property
    def concentration_matrix(self):
        """
        concentration_matrix (np.ndarray, None): vector representation of surfactant concentration in resevoir
        """
        return self._concentration_matrix
    @concentration_matrix.setter
    def concentration_matrix(self, value):
        self._concentration_matrix = value

    _IFT_conc_equ = None
    @property
    def IFT_conc_equ(self):
        """
        IFT_conc_equ (lambda, None): expression that relates surfactant concentration to interfacial tension b/t oil and water
        """
        return self._IFT_conc_equ
    @IFT_conc_equ.setter
    def IFT_conc_equ(self, value):
        self._IFT_conc_equ = value

    _derivative_IFT_conc_equ = None
    @property
    def derivative_IFT_conc_equ(self):
        """
        derivative_IFT_conc_equ (lambda, None): Deriviative of the equation relating IFT to surfactant concentration
        """
        return self._derivative_IFT_conc_equ
    @derivative_IFT_conc_equ.setter
    def derivative_IFT_conc_equ(self, value):
        self._derivative_IFT_conc_equ = value

    _is_surfactant = None
    @property
    def is_surfactant(self):
        """
        flag for whether or not surfactant is in the simulation
        """
        return self._is_surfactant
    @is_surfactant.setter
    def is_surfactant(self, value):
        self._is_surfactant = value

    _phi = None
    @property
    def phi(self):
        """
        phi (np.ndarray): arrray used to initialize the concentration matrix (represents porosity of the resevoir)
        """
        return self._phi
    @phi.setter
    def phi(self, value):
        self._phi = value

    def initialize(
        self,
    ):
        """
        This function will initialize the surfactant object
        
        Returns: (Surfactant)
            Surfactant object
        """
        if self.concentration_matrix is None:
            if self.phi is None:
                raise SimulationCalcInputException(
                    "SimulationInputException: phi value not initalized..."
                )
            D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
            self.concentration_matrix = (~D) * self.concentration
        return self

    def compute_concentration(
        self,
        grid: Grid,
        water_sat: np.ndarray,
        const_parameters: dict,
        varying_parameters: dict,
        F: np.ndarray,
        Gmod: np.ndarray,
    ):
        """
        Computing the surfactant concentration matrix

        Raises:
            SimulationCalcInputException: Not all required inputs were provided

        Args:
            grid (Grid): FD mesh

            water_sat (np.ndarray): water saturation matrix

            const_parameters (dict): dictionary object with constant parameters used in calculation

            varying_parameters (dict): dictionary object with varying parameters used in calculation

            F (np.ndarray): intermediate matrix used in calcs

            Gmod (np.ndarray): bilinear interpolant for sur conc on redefined coordinates

        Returns: (dict)
            Returns the ``varying_parameters`` dict
        """
        # initializing constants
        assert self.concentration_matrix is not None, SimulationCalcInputException(
            "SimuationInputException: polymer concentration matrix not initialized. Please try again"
        )
        # Required constants:
        dx = const_parameters["FD_grid_constants"]["dx"]
        dy = const_parameters["FD_grid_constants"]["dy"]
        x = const_parameters["FD_grid_constants"]["x"]
        y = const_parameters["FD_grid_constants"]["y"]
        m = const_parameters["FD_grid_constants"]["m"]
        n = const_parameters["FD_grid_constants"]["n"]
        phi = self.phi
        omega1 = const_parameters["Pc_constants"]["omega1"]
        omega2 = const_parameters["Pc_constants"]["omega2"]
        Qnew = water_sat
        g1 = const_parameters["inlet_total_flow"]
        g3 = const_parameters["inlet_surfactant_flow"]
        KK = const_parameters["KK"]
        relative_permeability_formula = const_parameters[
            "relative_permeability_formula"
        ]

        # retrieving relevant parameters for updating the water saturation
        ## Time Step:
        dt = const_parameters["FD_grid_constants"]["dt"]
        dt_array = const_parameters["FD_grid_constants"]["dt_matrix"]

        pc_g = varying_parameters["capillary_pressure_and_derivatives"]["dpc_dg"]
        lambda_a = varying_parameters["mobility_parameters"]["lambda_a"]
        lambda_total = varying_parameters["mobility_parameters"]["lambda_total"]

        # intermediate parameters for code:
        idx = 1
        AAA = np.zeros((n * m, n * m))
        DDD = np.zeros((n * m, 1))

        while idx <= (m) * (n - 1) + 1:
            cnt = (idx - 1) // m  # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n, m))
            AA = np.copy(BB)
            CC = np.copy(BB)
            DD = np.zeros((m, 1))
            for i in range(m):
                for j in range(n):
                    if j == i:
                        if idx == 1:
                            if i == 0:
                                DD[i] = (
                                    g3 / Qnew[cnt][i] + Gmod[cnt][i] / dt_array[cnt][i]
                                )
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
                                BB[j][i + 1] = F[cnt][i] / (dx**2)
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
                                    - ((g1 * lambda_a[cnt][i]) / (lambda_total[cnt][i]))
                                    / Qnew[cnt][i]
                                    + ((g3 * lambda_a[cnt][i]) / (lambda_total[cnt][i]))
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
                                BB[j][i - 1] = F[cnt][i] / (dx**2)
                                BB[j][i + 1] = F[cnt][i] / (dx**2)
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

        return varying_parameters
