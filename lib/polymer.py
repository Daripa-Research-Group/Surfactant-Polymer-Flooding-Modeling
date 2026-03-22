"""
This python script contains the class definition for polymers for the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
import scipy as sp
from scipy.sparse.linalg import bicgstab
from .enumerations import ModelType, PolymerList, SimulationConstants
from .Exceptions import SimulationCalcInputException
from .grid import Grid


class Polymer:
    """
    Class definition for the polymer objects and their calculations
    """

    def __init__(
        self,
        name: PolymerList,
        e_coeff: np.ndarray,
        n_coeff: np.ndarray,
        rho: float,
        concentration_scalar: float,
        phi: np.ndarray,
        viscosity_scalar: float | None = None,
        viscosity_matrix: np.ndarray | None = None,
        concentration_matrix: np.ndarray | None = None,
        shear_rate: np.ndarray | None = None,
    ):
        """
        Initializes a instance of the polymer class
        """

        # PolymerList object
        self.name = name

        # properties related to the concentration (scalar concentration, matrix version of initial concentration, and current concentration matrix)
        self.concetration_scalar = concentration_scalar
        self.init_concentration_matrix = concentration_matrix
        self.concentration_matrix = concentration_matrix

        # Properties related to the viscosity
        self.viscosity_matrix = viscosity_matrix
        self.viscosity_scalar = viscosity_scalar

        # Values required to formulate the numerical powerlaw function for viscosity calculations
        self.e_coeff = e_coeff
        self.n_coeff = n_coeff

        # Polymer Density
        self.rho = rho

        # shear rate matrix (needed when running 'shear thinning' model version)
        self.shear_rate = shear_rate

        # util param for initialization
        self.phi = phi  # Will need to be created in the simulation class
 
    _name = None
    @property
    def name(self):
        """
        name (enum 'PolymerList'): Name of the polymer
        """
        return self._name
    @name.setter
    def name(self, value):
        self._name = value

    _concetration_scalar = None
    @property
    def concetration_scalar(self):
        """
        concentration_scalar (float): Scalar quantity of concentration. When initializing, this param will equal the initial polymer concentration.
        """
        return self._concetration_scalar
    @concetration_scalar.setter
    def concetration_scalar(self, value):
        self._concetration_scalar = value

    _init_concentration_matrix = None
    @property
    def init_concentration_matrix(self):
        """
        Initial matrix (at time t = 0) representation of polymer concentration within resevoir
        """
        return self._init_concentration_matrix
    @init_concentration_matrix.setter
    def init_concentration_matrix(self, value):
        self._init_concentration_matrix = value

    _concentration_matrix = None
    @property
    def concentration_matrix(self):
        """
        concentration_matrix (np.ndarray, None): matrix representation of polymer concentration within resevoir over time
        """
        return self._concentration_matrix
    @concentration_matrix.setter
    def concentration_matrix(self, value):
        self._concentration_matrix = value

    _viscosity_matrix = None
    @property
    def viscosity_matrix(self):
        """
        viscosity_matrix (np.ndarray, None): viscosity matrix of the polymer
        """
        return self._viscosity_matrix
    @viscosity_matrix.setter
    def viscosity_matrix(self, value):
        self._viscosity_matrix = value

    _viscosity_scalar = None
    @property
    def viscosity_scalar(self):
        """
        viscosity_scalar (float, None): scalar quantity of the polymer viscosity
        """
        return self._viscosity_scalar
    @viscosity_scalar.setter
    def viscosity_scalar(self, value):
        self._viscosity_scalar = value

    _e_coeff = None
    @property
    def e_coeff(self):
        """
        e_coeff (list[float]): The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        """
        return self._e_coeff
    @e_coeff.setter
    def e_coeff(self, value):
        self._e_coeff = value

    _n_coeff = None
    @property
    def n_coeff(self):
        """
        n_coeff (list[float]):  The coefficients used to determine epsilon for the empirical power law expression used to determine the viscosity of the aqueous phase
        """
        return self._n_coeff
    @n_coeff.setter
    def n_coeff(self, value):
        self._n_coeff = value

    _rho = None
    @property
    def rho(self):
        """
        rho (float): Density of polymer
        """
        return self._rho
    @rho.setter
    def rho(self, value):
        self._rho = value

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

    _shear_rate = None
    @property
    def shear_rate(self):
        """
        shear_rate (np.ndarray, None): Matrix that will hold the shear rate (the change in velocity normal to the direction of flow)
        """
        return self._shear_rate
    @shear_rate.setter
    def shear_rate(self, value):
        self._shear_rate = value




    def initialize(self, grid_shape: tuple):
        """
        Will initialize the viscosity, shear_rate, and concentration matrices
        
        Args:
        -----
            grid_shape (tuple): contain the shape of the grid
        
        Returns: (Polymer)
        -----------------
            Initalized Polymer Object
        """
        n = grid_shape[0]
        m = grid_shape[1]
        D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
        if self.concentration_matrix is None:
            self.concentration_matrix = (~D) * self.concetration_scalar

        if self.init_concentration_matrix is None:
            self.init_concentration_matrix = self.concetration_scalar * np.ones(
                (n + 1, m + 1)
            )

        if self.shear_rate is None:
            self.shear_rate = np.zeros((n + 1, m + 1))

        if self.viscosity_scalar is None:
            beta1 = 15000  # constant that came from the MATLAB code
            self.viscosity_scalar = SimulationConstants.Water_Viscosity.value * (
                1 + beta1 * self.concetration_scalar
            )

        if self.viscosity_matrix is None:
            self.viscosity_matrix = self.viscosity_scalar * np.ones((n + 1, m + 1))

        return self

    def compute_viscosity(
        self,
        grid: Grid,
        u: np.ndarray,
        v: np.ndarray,
        model_type: ModelType,
        aqueous_viscosity: np.ndarray | None = None,
    ):
        """
        Compute polymer viscosity.
        This function is derived from 'compvis()' in the original MATLAB code (in file compvis.m).
        
        Args:
        -----
            grid (Tuple[NDArray[Any], ...]): The FEM grid used for simulation calculations (x and y variables from the MATLAB code)

            u (np.ndarray): Matrix related to the global pressure

            v (np.ndarray): Matrix related to the velocity matrix

            model_type (enum 'ModelType'): Will state whether the model will include polymer shear thinning or not

            aqueous_viscosity (np.ndarray, None): Aqueous viscosity matrix (will come from the ``Water`` class). Only needed when shear thinning OFF
        
        Returns: (list)
        ---------------
            the viscosity_matrix (index 0) & shear_rate matrix (index 1) for the polymer within the grid
        """
        # x and y components from meshgrid
        x = grid.x
        y = grid.y

        if (
            self.concentration_matrix is None
            or self.init_concentration_matrix is None
            or self.concetration_scalar is None
        ):
            raise SimulationCalcInputException(
                "SimulationInputException: Polymer concentration matrix and/or scalar concentration value not initialized..."
            )

        # if model_type is NO SHEAR THINNING:
        if model_type.value == ModelType.No_Shear_Thinning.value:
            if aqueous_viscosity is None:
                raise SimulationCalcInputException(
                    "SimulationInputException: Aqueous viscosity matrix required but not provided. Please try again."
                )
            ## the scalar viscosity is equal to the max within the aqueous viscosity matrix
            self.viscosity_scalar = np.max(aqueous_viscosity[0, :])
            self.viscosity_matrix = self.viscosity_scalar * np.ones(
                (
                    SimulationConstants.Grid_Size.value,
                    SimulationConstants.Grid_Size.value,
                )
            )
        # Model Type is 'Sourav Implementation':
        elif model_type.value == ModelType.Sourav_Implementation.value:
            # TODO: Will keep empty until properly understood how to implement
            pass
        # if polymer shear thinning is ON:
        elif model_type.value == ModelType.Shear_Thinning_On.value:
            if aqueous_viscosity is not None:
                raise SimulationCalcInputException(
                    "SimulationInputException: Aqueous viscosity reliant on changing polymer viscosity. Update will be done within 'Water' Class"
                )
            if self.shear_rate is None or self.viscosity_matrix is None:
                raise SimulationCalcInputException(
                    "SimulationInputException: Either shear_matrix or viscosity_matrix are not initialized"
                )
            # Getting water density and viscosity (Note: polymer density a property of class)
            rho_water = SimulationConstants.Water_Density.value
            viscosity_water = SimulationConstants.Water_Viscosity.value

            # Formulating the numerically derived power law equation
            w1_0 = (
                self.rho * self.init_concentration_matrix
            )  # from the variable w10 in MATLAB code
            w2_0 = rho_water * (
                1 - self.init_concentration_matrix
            )  # from the variable w20 in MATLAB code
            wppm_0 = (w1_0 / (w1_0 + w2_0)) * (
                10**6
            )  # from the wppm0 variable in MATLAB code
            print(f"type w1_0: {np.shape(w1_0)}")
            print(f"type w2_0: {np.shape(w2_0)}")

            ## Determining the epsilon and n coefficients for the power law equation
            epsilon_0 = np.zeros(
                (
                    np.size(self.concentration_matrix, 0),
                    np.size(self.concentration_matrix, 1),
                )
            )
            n_0 = np.zeros(
                (
                    np.size(self.concentration_matrix, 0),
                    np.size(self.concentration_matrix, 1),
                )
            )
            print(f"type epsilon_0: {np.shape(n_0)}")
            print(f"type n_0: {np.shape(n_0)}")
            for r in range(np.size(self.concentration_matrix, 0)):
                for c in range(np.size(self.concentration_matrix, 1)):
                    epsilon_0[r, c] = self.e_coeff[0] * wppm_0[r, c] ** self.e_coeff[1]
                    n_0[r, c] = min(
                        self.n_coeff[0] * wppm_0[r, c] ** self.n_coeff[1], 1
                    )

            row = np.size(self.concentration_matrix, 0)
            col = np.size(self.concentration_matrix, 1)

            # Compute divergence terms
            a1 = self.divergence(x, v)
            a2 = self.divergence(y, u)
            a3 = self.divergence(x, u)
            a4 = self.divergence(y, v)

            pi_D = np.abs(-0.25 * ((a1 + a2) ** 2) + a3 * a4)
            for i in range(row):
                for j in range(col):
                    if self.concentration_matrix[i, j] > 0:
                        self.shear_rate[i, j] = 2 * np.sqrt(pi_D[i, j])
                        if not (self.shear_rate[i, j] == 0):
                            self.viscosity_matrix[i, j] = epsilon_0[i, j] * (
                                self.shear_rate[i, j] ** (n_0[i, j] - 1)
                            )
                            print(f"epsilon_0:{epsilon_0[i,j]}")
                            print(f"n_0:{n_0[i,j]}")
                            print(f"shear_rate:{self.shear_rate[i,j]}")
                            print(f"pi_D: {pi_D[i,j]}")
                            print("")
                            if self.viscosity_matrix[i, j] < viscosity_water:
                                self.viscosity_matrix[i, j] = viscosity_water
                            if self.viscosity_matrix[i, j] > 100:
                                self.viscosity_matrix[i, j] = 100

        return [self.viscosity_matrix, self.shear_rate]

    def compute_concentration(
        self,
        grid: Grid,
        water_sat: np.ndarray,
        u: np.ndarray,
        v: np.ndarray,
        xmod: np.ndarray,
        ymod: np.ndarray,
        const_parameters: dict,
        varying_parameters: dict,
    ):
        """
        Computes the polymer concentration

        Raises:
        -------
            SimulationCalcInputException: Not all required parameters were provided

        Args:
        -----
            grid (Grid): the FDMesh

            water_sat (np.ndarray): the water saturation matrix

            u (np.ndarray): Matrix related to the global pressure

            v (np.ndarray): Matrix related to the velocity matrix

            xmod (np.ndarray): x-dim characteristic coordinates based on Neumann boundary conditions

            ymod (np.ndarray): y-dim characteristic coordinates based on Neumann boundary conditions

            const_parameters (dict): constant parameters to help with calculations

            varying_parameters (dict): parameters that vary but assist with calculations for water saturation, polymer concentration, and surfactant concentration
        
        Returns: (dict)
        ---------------
            The varying parameters that were changed in this method

        """
        # initializing variables:
        # Assert statements to ensure that all parameters are property initialized:
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
        C = np.copy(self.concentration_matrix)
        g1 = const_parameters["inlet_total_flow"]
        g2 = const_parameters["inlet_polymer_flow"]
        KK = const_parameters["KK"]
        relative_permeability_formula = const_parameters[
            "relative_permeability_formula"
        ]

        # retrieving relevant parameters for updating the water saturation
        ## Time Step:
        dt = const_parameters["FD_grid_constants"]["dt"]
        dt_array = const_parameters["FD_grid_constants"]["dt_matrix"]

        # retrieving fractional flow variable
        f = varying_parameters["fractional_flow_parameters"]["f"]

        # Determining 'Cmod'
        x1d = x[0, :]
        y1d = y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)

        # reorder vec_concentration if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            C = C[:, x_sort_idx]  # Sort columns of vec_concentration
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            C = C[y_sort_idx, :]  # Sort rows of vec_concentration

        interp = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d),
            self.concentration_matrix,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )
        query_points = np.stack([ymod.ravel(), xmod.ravel()], axis=-1)
        Cmod = interp(query_points).reshape(xmod.shape)

        # Using 'Cmod' and 'Qnew' to update the polymer concentration matrix
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
                        if idx == 1:  # lowermost row of grid
                            if i == 1:  # leftmost point (source)
                                DD[i] = (
                                    g2 / Qnew[cnt][i] + Cmod[cnt][i] / dt_array[cnt][i]
                                )
                                BB[j][i] = 1 / dt_array[cnt][i] + g1 / Qnew[cnt][i]
                            else:
                                DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = 1 / dt_array[cnt][i]
                        elif idx == (m) * (n - 1) + 1:
                            if i == m - 1:
                                DD[i] = Cmod[cnt][i] / dt_array[cnt][i]
                                BB[j][i] = (
                                    1 / dt_array[cnt][i] - g1 * f[cnt][i] / Qnew[cnt][i]
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

        return varying_parameters

    def divergence(self, Fx, Fy, dx=1.0, dy=1.0):
        """
        Calculates Divergence

        Args:
        -----
            Fx (np.ndarray): Function #1
            
            Fy (np.ndarray): Function #2

            dx (float): change in the x-dimension

            dy (float): change in the y-dimension

        Raises: (np.ndarray)
        --------------------
            Div F = (δfx/δx) + (δfy/δy)
        """
        dFx_dx = np.gradient(Fx, dx, axis=1)
        dFy_dy = np.gradient(Fy, dy, axis=0)
        return dFx_dx + dFy_dy
