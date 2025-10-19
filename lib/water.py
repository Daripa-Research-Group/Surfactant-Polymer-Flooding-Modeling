"""
This python script contains the class definition for the water in the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

##EXTERNAL IMPORTS
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix
from scipy.linalg import fractional_matrix_power
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import bicgstab

##INTERNAL IMPORTS
from enumerations import ModelType, RelativePermeabilityFormula, SimulationConstants
from Exceptions import SimulationCalcInputException
from grid import Grid
from polymer import Polymer
from surfactant import Surfactant


class Water:
    """
    Contains the properties and methods for water in the SP-Flooding system
    """

    def __init__(
        self,
        init_water_saturation: float,
        init_aqueous_saturation: float,
        init_oleic_saturation: float,
        miuw: float,
        miuo: float,
        phi: np.ndarray,
    ):
        """
        Constructor for the 'Water' class

        :param init_water_saturation: initial water saturation
        :type init_water_saturation: float

        :param init_aqueous_saturation: initial aqueous phase saturation
        :type init_aqueous_saturation: float

        :param init_oleic_saturation: initial oil phase saturation
        :type init_oleic_saturation: float

        :param miuw: water viscosity
        :type miuw: float

        :param miuo: oil viscosity
        :type miuo: float

        :param phi: porosity matrix
        :type phi: np.ndarray
        """
        self.init_water_saturation = init_water_saturation
        self.init_aqueous_saturation = init_aqueous_saturation
        self.init_oleic_saturation = init_oleic_saturation
        self.miuw = miuw
        self.miuo = miuo
        self.water_saturation = None  # water saturation matrix
        self.viscosity_array = None  # aqueous viscosity matrix
        self.phi = phi

    def initialize(self, grid_shape: tuple):
        """
        Initializing 'Water' Object properties

        :param grid_shape: the n and m parameters from the 'Grid' class
        :type grid_shape: tuple

        :return: Updated 'Water' object
        :rtype: Water
        """
        # getting values for n and m from grid:
        n, m = grid_shape

        # initializing the water saturation matrix
        s0 = np.zeros((n + 1, m + 1))
        D = (self.phi > 1e-10) | (np.abs(self.phi) < 1e-10)
        s0 = np.logical_not(D).astype(float) + D.astype(float) * (
            1 - self.init_water_saturation
        )
        self.water_saturation = s0

        # initialize the aqueous viscosity matrix:
        self.viscosity_array = self.miuw * np.ones((n + 1, m + 1))

        return self

    def compute_viscosity(
        self,
        grid: Grid,
        model_type: ModelType,
        polymer: Polymer,
        u: np.ndarray | None = None,
        v: np.ndarray | None = None,
    ):
        """
        Compute aqueous viscosity.

        :param grid: Grid object for deterrmining matrix size
        :type grid: Grid

        :param model_type: Type of model we are running (Polymer shear thinning ON or OFF)
        :type model_type: enum 'ModelType'

        :param polymer: holds the information about the polymer in the sim
        :type polymer: Polymer

        :param u: global pressure matrix. Only needed when shear thinning ON.
        :type u: np.ndarray, None

        :param v: velocity matrix. Only needed when shear thinning ON.
        :type v: np.ndarray, None

        :return: updated aqueous viscosity matrix
        :rtype: np.ndarray
        """
        assert self.viscosity_array is not None, SimulationCalcInputException(
            "SimuationInputException: aqueous viscosity matrix not initialized. Please try again"
        )
        assert polymer is not None and isinstance(
            polymer, Polymer
        ), SimulationCalcInputException(
            "SimuationInputException: polymer object not initialized. \
                The polymer object must be initialized before updating aqueous viscosity. Please try again."
        )
        assert polymer.concentration_matrix is not None, SimulationCalcInputException(
            "SimuationInputException: Polymer Concentration matrix must be initialized. Please try again..."
        )
        assert polymer.shear_rate is not None, SimulationCalcInputException(
            "SimuationInputException: Polymer Shear Rate matrix must be initialized. Please try again..."
        )
        n = np.size(polymer.concentration_matrix, 0)
        m = np.size(polymer.concentration_matrix, 1)
        initial_polymer_concentration_scalar = polymer.concetration_scalar
        if (
            model_type.value == ModelType.No_Shear_Thinning.value
        ):  # no shear thinning polymer
            miuw = SimulationConstants.Water_Viscosity.value
            if initial_polymer_concentration_scalar == 0:
                self.viscosity_array = miuw * np.ones((n, m))
            else:
                beta1 = SimulationConstants.beta1.value
                self.viscosity_array = miuw * (1 + beta1 * polymer.concentration_matrix)
        elif (
            model_type.value == ModelType.Shear_Thinning_On.value
        ):  # shear thinning polymer
            # using the shear rate and polymer coefficients to understand how its viscosity changes
            assert u is not None, SimulationCalcInputException(
                "SimuationInputException: variables 'u' not initialized for shear-thinning-on model. Please try again"
            )
            assert v is not None, SimulationCalcInputException(
                "SimuationInputException: variables 'v' not initialized for shear-thinning-on model. Please try again"
            )

            # constants:
            rho_water = SimulationConstants.Water_Density.value
            viscosity_water = SimulationConstants.Water_Viscosity.value

            # relevant parameters for power law equation:
            w1 = polymer.rho * polymer.concentration_matrix
            w2 = rho_water * (1 - polymer.concentration_matrix)
            wppm = (w1 / (w1 + w2)) * (10**6)

            # determining ε and n for power law equation:
            # epsilon_val = polymer.e_coeff[0]*(wppm**polymer.e_coeff[1])
            # n_val = min(polymer.n_coeff[0]*(wppm**polymer.n_coeff[1]),1)

            # epsilon_val = np.zeros((np.size(polymer.concentration_matrix, 0), np.size(polymer.concentration_matrix, 1)))
            # n_val = np.zeros((np.size(polymer.concentration_matrix, 0), np.size(polymer.concentration_matrix, 1)))
            # print(f'type epsilon_0: {np.shape(epsilon_val)}')
            # print(f'type n_0: {np.shape(n_val)}')
            # for r in range(np.size(polymer.concentration_matrix, 0)):
            #     for c in range(np.size(polymer.concentration_matrix, 1)):
            #         # print(f'i: {r} | j: {c}')
            #         epsilon_val[r,c] = polymer.e_coeff[0] * wppm[r,c] ** polymer.e_coeff[1]
            #         n_val[r,c] = min(polymer.n_coeff[0] * wppm[r,c] ** polymer.n_coeff[1], 1)
            eps = 1e-12
            wppm_safe = np.maximum(wppm, eps)

            epsilon_val = polymer.e_coeff[0] * wppm_safe ** polymer.e_coeff[1]
            n_val = np.minimum(polymer.n_coeff[0] * wppm_safe ** polymer.n_coeff[1], 1)

            row = np.size(polymer.concentration_matrix, 0)
            col = np.size(polymer.concentration_matrix, 1)

            # Compute divergence terms
            a1 = np.gradient(v, axis=0)
            a2 = np.gradient(u, axis=1)
            a3 = np.gradient(u, axis=0)
            a4 = np.gradient(v, axis=1)

            pi_D = np.abs(-0.25 * (a1 + a2) ** 2 + a3 * a4)

            for ii in range(row):
                for jj in range(col):
                    # Applying constraints
                    self.viscosity_array[ii, jj] = epsilon_val[ii, jj] * (
                        polymer.shear_rate[ii, jj] ** (n_val[ii, jj] - 1)
                    )
                    if self.viscosity_array[ii, jj] < viscosity_water:
                        self.viscosity_array[ii, jj] = viscosity_water
                    if self.viscosity_array[ii, jj] > 100:
                        self.viscosity_array[ii, jj] = 100
        return self.viscosity_array

    def compute_residual_saturations(
        self, sigma: np.ndarray, u: np.ndarray, v: np.ndarray
    ):
        """
        Compute swr, sor based on capillary numbers (came from compres.m MATLAB file)

        :param sigma: interfacial tension (IFT)
        :type sigma: np.ndarray

        :param u: global pressure matrix.
        :type u: np.ndarray

        :param v: velocity matrix.
        :type v: np.ndarray

        :return residual saturation for oil (index 1) and water (index 0) phases
        :rtype: list
        """
        swr0 = self.init_aqueous_saturation
        sor0 = self.init_oleic_saturation

        Nco0 = 1.44e-4
        Nca0 = 1.44e-4

        vel_mag = np.sqrt(np.matmul(u, u) + np.matmul(v, v))
        nca = (vel_mag * self.viscosity_array) / sigma
        nco = (vel_mag * self.miuo) / sigma

        Nca = np.linalg.norm(nca)
        Nco = np.linalg.norm(nco)

        sor = sor0 * (Nco0 / Nco) ** 0.5213 if Nco >= Nco0 else sor0
        swr = swr0 * (Nca0 / Nca) ** 0.1534 if Nca >= Nca0 else swr0

        return [swr, sor]  # [residual water saturation, residual oil saturation]

    def compute_mobility(
        self,
        c: np.ndarray,
        sor: float,
        swr: float,
        aqueous: bool,
        rel_permeability_formula: RelativePermeabilityFormula,
        modified_water_saturation: np.ndarray | None = None,
    ):
        """
        Computing mobility (made using the compmob.m MATLAB file)

        :param c: polymer concentration matrix
        :type c: np.ndarray

        :param sor: residual saturation oil phase
        :type sor: float

        :param swr: residual saturation water phase
        :type swr: float

        :param aqueous: boolean for whether we are solving for aqoeous or oleic mobility
        :type aqueous: bool

        :param rel_permeability_formula: Select the type of relative Permeability formula from the ``RelativePermeabilityFormula`` Enum
        :type has_surfactant: enum ``RelativePermeabilityFormula``

        :param surfactant_conc: scalar quantity of the initial surfactant concentration
        :type surfactant_conc: float

        :return: aqueous or oleic mobility (depending on the 'aqueous' parameter)
        :rtype: np.ndarray
        """
        assert self.water_saturation is not None, SimulationCalcInputException(
            "SimuationInputException: water saturation matrix not initialized. Please try again"
        )
        assert self.viscosity_array is not None, SimulationCalcInputException(
            "SimuationInputException: viscosity matrix not initialized. Please try again"
        )
        s = (
            self.water_saturation
            if (modified_water_saturation is None)
            else modified_water_saturation
        )
        miua = self.viscosity_array
        if (
            rel_permeability_formula.value
            == RelativePermeabilityFormula.CoreyTypeEquation.value
        ):
            nsw0 = (s - self.init_aqueous_saturation) / (
                1 - self.init_aqueous_saturation
            )
            nso0 = (s - self.init_aqueous_saturation) / (
                1 - self.init_aqueous_saturation - self.init_oleic_saturation
            )
            krw0 = nsw0**3.5
            kro0 = ((1 - nso0) ** 2) * (1 - nso0**1.5)
        else:
            nsw = (s - swr) / (1 - swr)
            nso = (s - swr) / (1 - swr - sor)
            krw0 = nsw * (2.5 * swr * (nsw**2 - 1) + 1)
            kro0 = (1 - nso) * (1 - 5 * sor * nso)

        return krw0 / miua if aqueous else kro0 / self.miuo

    def compute_water_saturation(
        self,
        grid: Grid,
        surfactant: Surfactant,
        polymer: Polymer,
        u: np.ndarray,
        v: np.ndarray,
        xmod: np.ndarray,
        ymod: np.ndarray,
        const_parameters: dict,
        varying_parameters: dict,
    ):
        """
        Solving saturation equation (comes from part of the nmmoc_surf_mod_neumann.m file that
        is for calculating the water saturation)

        :raises SimulationCalcInputException: If water saturation matrix is None

        :param grid: the 'Grid' object
        :type grid: Grid

        :param surfactant: The surfactant object
        :type surfactant: Surfactant

        :param polymer: The polymer object
        :type polymer: Polymer

        :param u: global pressure matrix
        :type u: np.ndarray

        :param v: velocity matrix
        :type v: np.ndarray

        :param xmod: x-dimension coordinate points for formulating the 'Qmod' matrix
        :type xmod: np.ndarray

        :param ymod: y-dimension coordinate points for formulating the 'Qmod' matrix
        :type ymod: np.ndarray

        :param const_parameters: constant parameters used in the method
        :type const_parameters: dict

        :param varying_parameters: parameters whose values can change
        :type varying_parameters: dict

        :return: updated ``water_saturation`` matrix and ``varying_parameters`` dict
        :rtype: [np.ndarray, list]
        """
        # Assert statements to ensure that all parameters are property initialized:
        assert self.water_saturation is not None, SimulationCalcInputException(
            "SimuationInputException: water saturation matrix not initialized. Please try again"
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
        Q = self.water_saturation
        g1 = const_parameters["inlet_total_flow"]
        KK = const_parameters["KK"]
        relative_permeability_formula = const_parameters[
            "relative_permeability_formula"
        ]

        # retrieving relevant parameters for updating the water saturation
        ## Time Step:
        dt = const_parameters["FD_grid_constants"]["dt"]
        dt_array = const_parameters["FD_grid_constants"]["dt_matrix"]

        # Determining Qmod matrix
        x1d = x[0, :]
        y1d = y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)

        # reorder Q if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            Q = Q[:, x_sort_idx]  # Sort columns of S
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            Q = Q[y_sort_idx, :]  # Sort rows of Q

        interp_func = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d), Q, method="linear", bounds_error=False, fill_value=None
        )

        query_points = np.stack([ymod.ravel(), xmod.ravel()], axis=-1)
        Qmod = interp_func(query_points).reshape(xmod.shape)

        swr = varying_parameters["swr"]
        sor = varying_parameters["sor"]
        nsw = (Qmod - swr) / (1 - swr)
        nso = (Qmod - swr) / (1 - swr - sor)
        varying_parameters["nsw"] = nsw
        varying_parameters["nso"] = nso

        ## fractional flow and derivatives
        assert polymer.concentration_matrix is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPolymerConcentrationMatrix"
        )
        lambda_a = self.compute_mobility(
            c=polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=True,
            rel_permeability_formula=relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_o = self.compute_mobility(
            c=polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=False,
            rel_permeability_formula=relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_total = lambda_a + lambda_o
        varying_parameters["mobility_parameters"] = {
            "lambda_a": lambda_a,
            "lambda_o": lambda_o,
            "lambda_total": lambda_total,
        }
        ## fractional flow calculations
        f = lambda_a / lambda_total
        assert KK is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPermeabilityTensor"
        )
        D = KK * lambda_o * f
        assert self.viscosity_array is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnkownWaterViscosityMatrix"
        )
        varying_parameters["fractional_flow_parameters"] = {"f": f, "D": D}
        pc = (
            surfactant.eval_IFT
            * const_parameters["Pc_constants"]["omega2"]
            * np.sqrt(const_parameters["porosity"])
        ) / (
            np.matmul(
                KK ** (0.5),
                fractional_matrix_power(
                    1 - nso, 1 / const_parameters["Pc_constants"]["omega1"]
                ),
            )
        )
        pc_s = pc / (const_parameters["Pc_constants"]["omega1"] * (1 - nso))
        pc_g = (pc / surfactant.eval_IFT) * surfactant.eval_dIFT_dGamma + pc_s
        varying_parameters["capillary_pressure_and_derivatives"] = {
            "pc": pc,
            "dpc_ds": pc_s,
            "dpc_dg": pc_g,
        }
        f_c = (-1 * (lambda_o * lambda_a * self.miuo)) / (
            (lambda_total**2) * self.viscosity_array
        )
        f_g = (
            varying_parameters["relative_permeability_derivatives"]["dkra_dg"]
            * lambda_o
        ) / ((lambda_total**2) * self.viscosity_array)
        varying_parameters["fractional_flow_derivatives"]["df_dc"] = f_c
        varying_parameters["fractional_flow_derivatives"]["df_dg"] = f_g
        D_g = D * pc_g
        D_s = D * pc_s
        varying_parameters["fractional_flow_derivatives"]["dD_dg"] = D_g
        varying_parameters["fractional_flow_derivatives"]["dD_ds"] = D_s

        # Updating coefficients with interpolated saturations
        idx = 1
        AAA = np.zeros((n * m, n * m))
        DDD = np.zeros((n * m, 1))

        while (
            idx <= m * (n - 1) + 1
            and surfactant.concentration_matrix is not None
            and polymer.concentration_matrix is not None
        ):
            cnt = (idx - 1) // m  # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n, m))
            AA = np.copy(BB)
            CC = np.copy(BB)
            DD = np.zeros((m, 1))

            #'cnt+1' in matlab is 'cnt' in python as matlab indexes from 1 but python indexes from 0
            print(f"DT is of this type: {type(dt)} of value = {dt}")
            for i in range(m):
                for j in range(n):
                    if j == i:
                        if idx == 1:
                            if i == 0:  # first/left column
                                DD[i] = (
                                    (Qmod[cnt][i] / dt_array[cnt][i])
                                    + g1 * (1 - f[cnt][i])
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i + 1]) / (dx**2)
                                        + (D_g[cnt + 1][i] + D_g[cnt + 1][i]) / (dx**1)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i + 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i + 1]
                                    - (D_g[cnt][i] + D_g[cnt + 2][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt + 2][i]
                                )

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (dy**2)

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt + 1][i] + D_s[cnt][i + 1]) / (dx**2)
                                    - (D_s[cnt + 1][i] + D_s[cnt][i]) / (dy**2)
                                )

                                BB[j][i + 1] = (D_s[cnt][i] + D_s[cnt][i + 1]) * (dx**2)
                            elif i == m - 1:  # last/rightmost column
                                DD[i] = (
                                    Qmod[cnt][i] / dt_array[cnt][i]
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i - 1]) / (dx**2)
                                        + (D_g[cnt + 1][i] + D_g[cnt][i]) / (dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i - 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i - 1]
                                    - (D_g[cnt][i] + D_g[cnt + 1][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt + 1][i]
                                )

                                BB[j][i - 1] = (D_s[cnt][i] + D_s[cnt][i - 1]) / (dx**2)

                                BB[i][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt][i] + D_s[cnt][i - 1]) / (dx**2)
                                    - (D_s[cnt + 1][i] + D_s[cnt][i]) / (dy**2)
                                )

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (dy**2)
                            else:
                                DD[i] = (
                                    Qmod[cnt][i] / dt_array[cnt][i]
                                    - f_c[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt][i + 1]
                                            - polymer.concentration_matrix[cnt][i - 1]
                                        )
                                        / (2 * dx)
                                    )
                                    - f_g[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt][i + 1]
                                            - surfactant.concentration_matrix[cnt][
                                                i - 1
                                            ]
                                        )
                                        / (2 * dx)
                                    )
                                    + (
                                        (
                                            D_g[cnt][i + 1]
                                            + D_g[cnt][i - 1]
                                            + 2 * D_g[cnt][i]
                                        )
                                        / (2 * dx**2)
                                        + (D_g[cnt][i + 1] + D_g[cnt][i]) / (dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i + 1] + D_g[cnt][i])
                                    / (2 * dx**2)
                                    * surfactant.concentration_matrix[cnt][i + 1]
                                    - (D_g[cnt][i - 1] + D_g[cnt][i])
                                    / (2 * dx**2)
                                    * surfactant.concentration_matrix[cnt][i - 1]
                                    - (D_g[cnt][i] + D_g[cnt + 1][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt + 1][i]
                                )

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (
                                        D_s[cnt][i + 1]
                                        + D_s[cnt][i - 1]
                                        + 2 * D_s[cnt][i]
                                    )
                                    / (2 * dx**2)
                                    - (D_s[cnt + 1][i] + D_s[cnt][i]) / (dy**2)
                                )
                                BB[j][i - 1] = (D_s[cnt][i - 1] + D_s[cnt][i]) / (
                                    2 * dx**2
                                )
                                BB[j][i + 1] = (D_s[cnt][i + 1] + D_s[cnt][i]) / (
                                    2 * dx**2
                                )

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (dy**2)

                        elif idx == (m) * (n - 1) + 1:  # topmost row of grid
                            if i == 0:  # leftmost column
                                DD[i] = (
                                    (Qmod[cnt][i] / dt_array[cnt][i])
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i + 1]) / (dx**2)
                                        + (D_g[cnt - 1][i] + D_g[cnt][i]) / (dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i + 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i + 1]
                                    - (D_g[cnt][i] + D_g[cnt - 1][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt - 1][i]
                                )

                                AA[j][i] = (D_s[cnt][i] + D_s[cnt - 1][i]) / (dy**2)

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt][i] + D_s[cnt][i + 1]) / (dx**2)
                                    - (D_s[cnt - 1][i] + D_s[cnt][i]) / (dy**2)
                                )

                                BB[j][i + 1] = (D_s[cnt][i] + D_s[cnt][i + 1]) / (dx**2)
                            elif i == m - 1:  # rightmost column
                                DD[i] = (
                                    (Qmod[cnt][i] / dt_array[cnt][i])
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i - 1]) / (dx**2)
                                        + (D_g[cnt - 1][i] + D_g[cnt][i]) / (dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i - 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i - 1]
                                    - (D_g[cnt][i] + D_g[cnt - 1][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt - 1][i]
                                )

                                BB[j][i - 1] = (D_s[cnt][i] + D_s[cnt][i - 1]) / (dy**2)

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt][i] + D_s[cnt][i - 1]) / (dx**2)
                                    - (D_s[cnt - 1][i] + D_s[cnt][i]) / (dy**2)
                                )

                                AA[j][i] = (D_s[cnt][i] + D_s[cnt - 1][i]) / (dy**2)
                            else:
                                DD[i] = (
                                    Qmod[cnt][i] / dt_array[cnt][i]
                                    - f_c[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt][i + 1]
                                            - polymer.concentration_matrix[cnt][i - 1]
                                        )
                                        / (2 * dx)
                                    )
                                    - f_g[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt][i + 1]
                                            - surfactant.concentration_matrix[cnt][
                                                i - 1
                                            ]
                                        )
                                        / (2 * dx)
                                    )
                                    + (
                                        (
                                            D_g[cnt][i + 1]
                                            + D_g[cnt][i - 1]
                                            + 2 * D_g[cnt][i]
                                        )
                                        / (2 * dx**2)
                                        + (D_g[cnt][i + 1] + D_g[cnt][i]) / (dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i + 1] + D_g[cnt][i])
                                    / (2 * dx**2)
                                    * surfactant.concentration_matrix[cnt][i + 1]
                                    - (D_g[cnt][i - 1] + D_g[cnt][i])
                                    / (2 * dx**2)
                                    * surfactant.concentration_matrix[cnt][i - 1]
                                    - (D_g[cnt][i] + D_g[cnt - 1][i])
                                    / (dy**2)
                                    * surfactant.concentration_matrix[cnt - 1][i]
                                )

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (
                                        D_s[cnt][i + 1]
                                        + D_s[cnt][i - 1]
                                        + 2 * D_s[cnt][i]
                                    )
                                    / (2 * dx**2)
                                    - (D_s[cnt - 1][i] + D_s[cnt][i]) / (dy**2)
                                )
                                BB[j][i - 1] = (D_s[cnt][i - 1] + D_s[cnt][i]) / (
                                    2 * dx**2
                                )
                                BB[j][i + 1] = (D_s[cnt][i + 1] + D_s[cnt][i]) / (
                                    2 * dx**2
                                )

                                AA[j][i] = (D_s[cnt][i] + D_s[cnt - 1][i]) / (dy**2)

                        else:
                            if i == 0:
                                DD[i] = (
                                    Qmod[cnt][i] / dt_array[cnt][i]
                                    - f_c[cnt][i]
                                    * (
                                        v[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt + 1][i]
                                            - polymer.concentration_matrix[cnt][i]
                                        )
                                        / (2 * dy)
                                    )
                                    - f_g[cnt][i]
                                    * (
                                        v[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt + 1][i]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        / (2 * dy)
                                    )
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i + 1]) / (dx**2)
                                        + (
                                            D_g[cnt - 1][i]
                                            + 2 * D_g[cnt][i]
                                            + D_g[cnt + 1][i]
                                        )
                                        / (2 * dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i + 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i + 1]
                                    - (D_g[cnt][i] + D_g[cnt + 1][i])
                                    / (2 * dy**2)
                                    * surfactant.concentration_matrix[cnt + 1][i]
                                    - (D_g[cnt][i] + D_g[cnt - 1][i])
                                    / (2 * dy**2)
                                    * surfactant.concentration_matrix[cnt - 1][i]
                                )

                                AA[j][i] = (D_s[cnt][i] + D_s[cnt - 1][i]) / (2 * dy**2)

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt][i] + D_s[cnt][i + 1]) / (dx**2)
                                    - (
                                        D_s[cnt - 1][i]
                                        + 2 * D_s[cnt][i]
                                        + D_s[cnt + 1][i]
                                    )
                                    / (2 * dy**2)
                                )
                                BB[j][i + 1] = (D_s[cnt][i + 1] + D_s[cnt][i]) / (dx**2)

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (2 * dy**2)
                            elif i == m - 1:
                                DD[i] = (
                                    Qmod[cnt][i] / dt_array[cnt][i]
                                    - f_c[cnt][i]
                                    * (
                                        v[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt + 1][i]
                                            - polymer.concentration_matrix[cnt][i]
                                        )
                                        / (2 * dy)
                                    )
                                    - f_g[cnt][i]
                                    * (
                                        v[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt + 1][i]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        / (2 * dy)
                                    )
                                    + (
                                        (D_g[cnt][i] + D_g[cnt][i - 1]) / (dx**2)
                                        + (
                                            D_g[cnt - 1][i]
                                            + 2 * D_g[cnt][i]
                                            + D_g[cnt + 1][i]
                                        )
                                        / (2 * dy**2)
                                    )
                                    * surfactant.concentration_matrix[cnt][i]
                                    - (D_g[cnt][i] + D_g[cnt][i - 1])
                                    / (dx**2)
                                    * surfactant.concentration_matrix[cnt][i - 1]
                                    - (D_g[cnt][i] + D_g[cnt + 1][i])
                                    / (2 * dy**2)
                                    * surfactant.concentration_matrix[cnt + 1][i]
                                    - (D_g[cnt][i] + D_g[cnt - 1][i])
                                    / (2 * dy**2)
                                    * surfactant.concentration_matrix[cnt - 1][i]
                                )

                                AA[j][i] = (D_s[cnt][i] + D_s[cnt - 1][i]) / (2 * dy**2)

                                BB[j][i] = (
                                    1 / dt_array[cnt][i]
                                    - (D_s[cnt][i] + D_s[cnt][i - 1]) / (dx**2)
                                    - (
                                        D_s[cnt - 1][i]
                                        + 2 * D_s[cnt][i]
                                        + D_s[cnt + 1][i]
                                    )
                                    / (2 * dy**2)
                                )
                                BB[j][i - 1] = (D_s[cnt][i] + D_s[cnt][i - 1]) / (dx**2)

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (2 * dy**2)
                            else:
                                DD[i] = (
                                    (Qmod[cnt][i] / dt_array[cnt][i])
                                    - f_c[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt][i + 1]
                                            - polymer.concentration_matrix[cnt][i - 1]
                                        )
                                        / (2 * dx)
                                        + v[cnt][i]
                                        * (
                                            polymer.concentration_matrix[cnt + 1][i]
                                            - polymer.concentration_matrix[cnt - 1][i]
                                        )
                                        / (2 * dy)
                                    )
                                    - f_g[cnt][i]
                                    * (
                                        u[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt][i + 1]
                                            - surfactant.concentration_matrix[cnt][
                                                i - 1
                                            ]
                                        )
                                        / (2 * dx)
                                        + v[cnt][i]
                                        * (
                                            surfactant.concentration_matrix[cnt + 1][i]
                                            - surfactant.concentration_matrix[cnt - 1][
                                                i
                                            ]
                                        )
                                        / (2 * dy)
                                    )
                                    - (
                                        D_g[cnt][i + 1]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt][i + 1]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        - D_g[cnt][i - 1]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt][i - 1]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        + D_g[cnt][i - 1]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt][i - 1]
                                            - surfactant.concentration_matrix[cnt][
                                                i + 1
                                            ]
                                        )
                                        + D_g[cnt + 1][i]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt + 1][i]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        + D_g[cnt - 1][i]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt - 1][i]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                        + D_g[cnt][i]
                                        / (2 * dx**2)
                                        * (
                                            surfactant.concentration_matrix[cnt + 1][i]
                                            - surfactant.concentration_matrix[cnt][i]
                                        )
                                    )
                                )
                                AA[j][i] = (D_s[cnt - 1][i] + D_s[cnt][i]) / (2 * dy**2)

                                CC[j][i] = (D_s[cnt][i] + D_s[cnt + 1][i]) / (2 * dy**2)

                                BB[j][i] = 1 / dt_array[cnt][i] - (
                                    (1 / (2 * dx**2))
                                    * (D_s[cnt][i] + 2 * D_s[cnt][i] + D_s[cnt][i + 1])
                                    + (1 / (2 * dy**2))
                                    * (
                                        D_s[cnt - 1][i]
                                        + 2 * D_s[cnt][i]
                                        + D_s[cnt + 1][i]
                                    )
                                )
                                BB[j][i + 1] = (D_s[cnt][i] + D_s[cnt][i + 1]) / (
                                    2 * dx**2
                                )
                                BB[j][i - 1] = (D_s[cnt][i - 1] + D_s[cnt][i]) / (
                                    2 * dx**2
                                )

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

        # Entering the bicgstab for the saturation calculations
        # bicgstab (Biconjugate Gradient Stabilized) - Iterative algorithm to solve large, sparse, and non-symmetric linear systems of the form Ax = b

        Qnew_flat, info = bicgstab(AAA, DDD, rtol=10 ** (-10), maxiter=600)
        Qnew = Qnew_flat = Qnew_flat.reshape(m, n)

        Qnew[Qnew > 1] = 1

        self.water_saturation = Qnew

        return [Qnew, varying_parameters]
