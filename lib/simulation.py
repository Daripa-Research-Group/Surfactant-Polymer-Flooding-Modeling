"""
This python script contains the class definition for running simulations

This Python code has been derived from the MATLAB Surfactant-Polymer Flooding Simulation 
developed by Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar and Carlos Acosta Caripo

"""

import os

from .grid import Grid, FEMesh
from .enumerations import (
    ModelType,
    PolymerList,
    RelativePermeabilityFormula,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from .Exceptions import SimulationCalcInputException, UserInputException
from .polymer import Polymer
from .surfactant import Surfactant
import numpy as np
import scipy as sp
from .water import Water
from scipy.io import loadmat
from scipy.sparse.linalg import bicgstab
from scipy.linalg import fractional_matrix_power
from scipy.interpolate import RegularGridInterpolator

os.makedirs(
    "memmaps", exist_ok=True
)  # ensures that the program works on computers with RAM constraints


class Simulation:
    """
    Simulation class to run SP-flooding simulations based on MATLAB translation.
    """

    def __init__(self, user_input_dict: dict):
        """
        This method will check the ``user_input_dict`` and initialize the simulation
        
        Raises:
        ------
            UserInputException: If there is a issue with the user inputs in ``user_input_dict``
            SimulationCalcInputException: If there is an issue with the execution of a calculation during runtime
        
        Args:
        -----
            user_input_dict (dict): dictionary containing the information from the GUI
        """
        ## Performs checks on the user input dictionary passed in:
        try:
            model_type = ModelType(user_input_dict["model_type"])
            if model_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadModelTypeAssignment", user_input_dict
            )

        try:
            reservoir_geometry = ResevoirGeometry(user_input_dict["reservoir_geometry"])
            if reservoir_geometry is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadReservoirGeometryAssignment", user_input_dict
            )

        try:
            permeability_flag = PermeabilityType(user_input_dict["permeability"])
            if permeability_flag is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadPermeabilityAssignment", user_input_dict
            )

        try:
            polymer_type = PolymerList.get_by_value(user_input_dict["polymer_type"])
            if polymer_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadPolymerTypeAssignment", user_input_dict
            )

        try:
            polymer_concentration = user_input_dict["polymer_concentration"]
            if polymer_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadPolymerConcentrationAssignment", user_input_dict
            )

        try:
            surfactant_type = SurfactantList.get_by_value(
                user_input_dict["surfactant_type"]
            )
            if surfactant_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadSurfactantTypeAssignment", user_input_dict
            )

        try:
            surfactant_concentration = user_input_dict["surfactant_concentration"]
            if surfactant_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadSurfactantConcentrationAssignment", user_input_dict
            )

        ## Instantiates the required simulation properties:

        # initializing properties that hold simulation constants
        self.grid_size = SimulationConstants.Grid_Size.value
        self.source_flow_magnitude = SimulationConstants.Source_Flow_Magnitude.value

        # Initializing Simulation Flags:
        self.permeability_flag = permeability_flag
        self.reservoir_geometry = reservoir_geometry
        self.model_type = model_type
        self.relative_permeability_formula = (
            RelativePermeabilityFormula.AmaefuleHandEquation
        )  # WILL NEED TO UPDATE TO INCLUDE IN GUI

        # Initializing sim properties
        self.mesh, self.FE_mesh = self._create_mesh()
        grid_shape = (self.mesh.n, self.mesh.m)
        self.x, self.y = self.mesh.get_meshgrid
        self.phi = None  # FIXME: Should check if this functionality code
        self.KK = None  # Permeability tensor
        self.time_step = None
        self.u, self.v = self._initialize_pressure_and_velocity()
        self._initialize_simulation()  # will initialize phi, KK, and time_step
        if (
            self.phi is None or self.KK is None or self.time_step is None
        ):  # Raise Exception if not properly initialized...
            raise SimulationCalcInputException(
                "SimulationCalcInputError:BadInitialSimulationPropertiesCalculation"
            )

        # Initalizing Polymer Object
        self.polymer = Polymer(
            name=polymer_type,
            e_coeff=polymer_type.e_coeff,
            n_coeff=polymer_type.n_coeff,
            rho=polymer_type.Density,
            concentration_scalar=polymer_concentration,
            phi=self.phi,
        )
        self.polymer.initialize(grid_shape=grid_shape)

        # Initializing Surfactant Object
        self.surfactant = Surfactant(
            name=surfactant_type,
            initial_concentration=surfactant_concentration,
            IFT_equation=surfactant_type.IFT_equation,
            derivative_IFT_equation=surfactant_type.derivative_IFT_equation,
            phi=self.phi,
        )
        self.surfactant.initialize()

        # Initializing Water Object
        self.water = Water(
            init_water_saturation=SimulationConstants.Initial_Residual_Water_Saturation.value,
            init_aqueous_saturation=SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value,
            init_oleic_saturation=SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value,
            miuw=SimulationConstants.Water_Viscosity.value,
            miuo=SimulationConstants.Oil_Viscosity.value,
            phi=self.phi,
        )
        self.water.initialize(grid_shape=grid_shape)

        # Properties for Exporting Simulation Results
        self.COC = np.zeros((1, 2000))  # Cumulative Oil Recovered
        self.miuaTcal = np.zeros((1, 2000))  # Total Aqueous Viscosity
        self.lambdaTcal = np.zeros((1, 2000))  # Total Mobility (λ_o + λ_a)

        ##The following properties require memmaps:
        self.ProdRate, self.CROIP = (
            self._initialize_memmap_properties()
        )  # ProdRate (Production Rate) / CROIP (Cummulative Remaining Oil In Place)
        self.MFW = []
        self.integrated_inlet_flow = 0  # "src_total" in the MATLAB version of the code

    # Dependent Property of Simulation Class
    _source_prod_flow = None

    @property
    def source_prod_flow(self) -> np.ndarray:
        """
        Returns: (np.ndarray)
        ---------------------
            The matrix with the source & and production well flow rates
        """
        # setting permeability state
        if self._source_prod_flow is None:
            self._source_prod_flow = np.zeros((self.mesh.n + 1, self.mesh.m + 1))
            bool_Homogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Homogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Quarter_Five_Spot = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (
                self.reservoir_geometry.value
                == ResevoirGeometry.Quarter_Five_Spot.value
            )
            if bool_Homogenous_and_Rectilinear or bool_Heterogenous_and_Rectilinear:
                self._source_prod_flow[:, 0] = (
                    self.source_flow_magnitude
                )  # intensity of injection well = src
                self._source_prod_flow[:, -1] = (
                    -1 * self.source_flow_magnitude
                )  # intensity of production well = -src
            elif bool_Heterogenous_and_Quarter_Five_Spot:  # Quarter-Five Spot
                self._source_prod_flow[0, 0] = (
                    self.source_flow_magnitude
                )  # Intensity of injection well = src
                self._source_prod_flow[-1, -1] = (
                    -1 * self.source_flow_magnitude
                )  # Intensity of production well = -src

        return self._source_prod_flow

    _scenario_flag = None

    @property
    def scenario_flag(self) -> int:
        """
        Determines the scenario based on the chosen reservoir geometry and permeability

        Returns: (int)
        --------------
            Integer value that represents a type of scenario run
        """
        if self._scenario_flag is None:
            bool_Homogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Homogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Quarter_Five_Spot = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (
                self.reservoir_geometry.value
                == ResevoirGeometry.Quarter_Five_Spot.value
            )

            if bool_Homogenous_and_Rectilinear:
                self._scenario_flag = 1
            elif bool_Heterogenous_and_Rectilinear:
                self._scenario_flag = 2
            elif bool_Heterogenous_and_Quarter_Five_Spot:
                self._scenario_flag = 3
            else:
                raise SimulationCalcInputException(
                    "SimulationCalcInputError:InvalidSimulationCase"
                )

        return self._scenario_flag

    # Private Methods of the Simulation Class
    def _initialize_pressure_and_velocity(self):
        """
        (private method)

        Initializing global pressure ('u') and velocity matrices ('v')
        Will use the ``n`` and ``m`` properties from ``Grid`` Class for initialization

        Returns: (tuple[np.ndarray, np.ndarray])
        ----------------------------------------
            The global pressure matrix (index 0) and velocity matrix (index 1)
        """
        u = np.zeros((self.mesh.n + 1, self.mesh.m + 1), dtype=np.complex128)
        v = np.zeros((self.mesh.n + 1, self.mesh.m + 1), dtype=np.complex128)

        return u, v

    def _compute_pressure_and_velocity_matrices(self, sparsed_A, B, beta):
        """
        (private method)

        dependent property to calculate the pressure matrix (``u``)
        and velocity matrix (``v``). Will rely on functions in the ``FEMesh`` class.
        
        Returns: (list[np.ndarray])
        ---------------------------
            List of update pressure matrix and velocity matrix => [u,v]
        """
        max_iterations = 1000
        new_u, convergence_flag = bicgstab(
            sparsed_A, B, rtol=1e-10, atol=0, maxiter=max_iterations
        )  # new_u is of shape (900, )
        if convergence_flag != 0:
            import warnings
            warnings.warn(f"BiCGSTAB: convergence issue (info={convergence_flag}) in pressure solver")

        new_v = np.zeros((self.FE_mesh.n + 1, self.FE_mesh.m + 1), dtype=object)
        for i in range(self.FE_mesh.m + 1):
            for j in range(self.FE_mesh.n + 1):
                new_v[j, i] = new_u[j * (self.FE_mesh.m + 1) + i]

        px, py = self._get_gradient(new_v)
        new_u = (-1 * beta) * px
        new_v = (-1 * beta) * py

        return new_u, new_v

    def _get_gradient(self, vn):
        """
        (private method)

        Helper function to determine the gradients with respect to x and y dimensions
        
        Returns: (tuple[_Array[tuple[int, int], float64], NDArray[float64]])
        --------------------------------------------------------------------
            Tuple with px py which are numpy matrices that hold the gradient wrt to x and y dimensions
        """
        m = self.mesh.m
        n = self.mesh.n

        dx = self.mesh.dx
        dy = self.mesh.dy

        px = np.zeros((n + 1, m + 1))
        py = np.copy(px)

        for i in range(m + 1):
            for j in range(n + 1):
                if i != 0:
                    px[j, i] = (vn[j, i] - vn[j, i - 1]) / dx
                if i != m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i]) / dx
                if i != 0 and i != m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i - 1]) / (2 * dx)
                if j != 0:
                    py[j, i] = (vn[j, i] - vn[j - 1, i]) / dy
                if j != n:
                    py[j, i] = (vn[j + 1, i] - vn[j, i]) / dy
                if j != 0 and j != n:
                    py[j, i] = (vn[j + 1, i] - vn[j - 1, i]) / (2 * dy)

        return px, py

    def _initialize_memmap_properties(self):
        """
        (private method)

        Will initialize the ``ProdRate`` and ``CROIP`` properties.
        Using memmaps to allow window's users to run program.
        
        Returns: (tuple[np.ndarray, np.ndarray])
        ----------------------------------------
            Initialized ``ProdRate`` and ``CROIP`` properties
        """
        os.makedirs("memmaps", exist_ok=True)
        tf = 500
        dt = self.mesh.dx / self.source_flow_magnitude
        timestamps = int(np.floor(tf / dt))
        CROIP = np.memmap(
            "memmaps/CROIP.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )
        ProdRate = np.memmap(
            "memmaps/ProdRate.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )

        return ProdRate, CROIP

    def _create_mesh(self):
        """
        (private method)

        Returns: (Tuple[Grid, FEMesh])
        -----------------------------
            Initialized FD and FE mesh
        """
        FD_mesh = Grid(self.grid_size, self.grid_size)

        FE_mesh = FEMesh(self.grid_size, self.grid_size)

        return FD_mesh, FE_mesh

    def _initialize_simulation(self):
        """
        (private method)

        Sets up initial reservoir fields, permeability, and time step.

        Returns: (None)
        ---------------
            Initialized properties of simulation. Required in ``__init__`` function
        """
        self.phi = self._compute_phi()

        # Initialize permeability matrix
        self.KK = self._compute_permeability()

        # Calculate time step
        self.time_step = self.mesh.dx / self.source_flow_magnitude
        if (
            self.permeability_flag == PermeabilityType.Heterogenous
            and self.reservoir_geometry == ResevoirGeometry.Quarter_Five_Spot
        ):
            self.time_step *= 100

    def _compute_phi(self):
        """
        (private method)

        Returns: (np.ndarray)
        ---------------------
            Computes and returns the level set function phi at each grid point. 
            (Equivalent to MATLAB get_phi_test function.)
        """
        m = self.mesh.m
        n = self.mesh.n
        dx = self.mesh.dx
        dy = self.mesh.dy
        left = self.mesh.left
        bottom = self.mesh.bottom

        jj, ii = np.meshgrid(np.arange(1, n + 2), np.arange(1, m + 2))
        x_coords = left + (ii - 1) * dx
        y_coords = bottom + (jj - 1) * dy

        phi = self._z_func_test(x_coords, y_coords)
        return phi

    def _z_func_test(self, x, y):
        """
        (private method)

        Args:
        -----
            x (np.ndarry): x-dimension coordinates 

            y (np.ndarray): y-dimension coordinate points
        
        Returns: 
        --------
            Compute and returns the initial position of the water front.
            (Equivalent to MATLAB z_func_test.)

        """
        init_front_hs = 0.1
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        if bool_Homogenous_and_Rectilinear:
            # Homogeneous
            out = y - init_front_hs + 0.01 * np.cos(80 * np.pi * x)
        elif bool_Heterogenous_and_Rectilinear:
            # Rectilinear Heterogeneous
            out = y - init_front_hs
        elif bool_Heterogenous_and_Quarter_Five_Spot:
            # Quarter five spot
            out = np.square(x) + np.square(y) - 0.015
        else:
            out = np.zeros_like(x)

        return out

    def _compute_permeability(self):
        """
        (private method)
        
        Returns:
        -------
            Compute and returns the permeability matrix KK based on the ``scenario_flag``.
            (Equivalent to MATLAB KKdef function.)
        """
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        m = self.mesh.m
        KK = None
        if bool_Homogenous_and_Rectilinear:
            # Homogeneous
            Kmax = 100
            KK = Kmax * np.ones((m + 1, m + 1))
        elif bool_Heterogenous_and_Rectilinear:
            # Heterogeneous
            Kmax = 100
            KK = Kmax * (
                0.5
                * (1 - 10 ** (-7))
                * (
                    np.sin(6 * np.pi * np.cos(self.x))
                    * np.cos(4 * np.pi * np.sin(3 * self.y))
                    - 1
                )
                + 1
            )
        elif bool_Heterogenous_and_Quarter_Five_Spot:
            # Load Upper Ness formation (SPE10)
            mat_data = loadmat(
                "./lib/Resources/KK30Ness.mat"
            )  # FIXME: when using master_surf_grid need to change this path
            if "KK" not in mat_data:
                raise SimulationCalcInputException(
                    "SimulationInputException: KK matrix not found in KK30Ness.mat file."
                )
            KK = mat_data["KK"]
        return KK

    def _transport_equation_solver(self, dt):
        """
        (Private method)

        This method of the ``Simulation`` class will solve the transport equations
        for the surfactant concentration, polymer concentration, and water saturation

        Based on the 'nmmoc_surf_mod_neumann.m' function within the MATLAB version

        Note:
            dg - derivative with respect to surfactant concentration
            ds - derivative with respect to water saturation
            dc - derivative with respect to polymer concentration

        Returns: (tuple[np.ndarray, np.ndarray, np.ndarray])
        ----------------------------------------------------
            tuple[Water, Polymer, Surfactant]
        """
        # Initialize constant parameters
        const_parameters = {}
        const_parameters["inlet_total_flow"] = self.source_flow_magnitude  # G1
        const_parameters["inlet_polymer_flow"] = (
            self.polymer.concetration_scalar * self.source_flow_magnitude
        )  # G2
        const_parameters["inlet_surfactant_flow"] = (
            self.surfactant.concentration * self.source_flow_magnitude
        )  # G3
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UndefinedWaterSaturationMatrix"
        )
        n = np.shape(self.water.water_saturation)[0]
        m = np.shape(self.water.water_saturation)[1]
        const_parameters["FD_grid_constants"] = {
            "n": n,
            "m": m,
            "dx": self.mesh.dx,
            "dy": self.mesh.dy,
            "dt": dt,
            "dt_matrix": dt * np.ones((n, m)),
            "x": self.mesh.x,
            "y": self.mesh.y,
        }
        ## parameters around Pc calculations
        const_parameters["Pc_constants"] = {
            "omega1": SimulationConstants.Capillary_Pressure_Param_1.value,
            "omega2": SimulationConstants.Capillary_Pressure_Param_2.value,
        }
        ## parameters around capillary number for residual saturation calcs
        const_parameters["resid_saturation_constants"] = {
            "Nco0": 10 ** (-5),
            "Nca0": 10 ** (-5),
        }
        ## porosity parameter
        const_parameters["porosity"] = 1

        ## Permeability
        const_parameters["KK"] = self.KK
        const_parameters["relative_permeability_formula"] = (
            self.relative_permeability_formula
        )

        # Initialize variable parameters
        varying_parameters = {}
        ## residual water and oil saturations
        [swr, sor] = self.water.compute_residual_saturations(
            sigma=self.surfactant.eval_IFT, u=self.u, v=self.v
        )
        varying_parameters["swr"] = swr
        varying_parameters["sor"] = sor
        nsw = (self.water.water_saturation - swr) / (1 - swr)
        nso = (self.water.water_saturation - swr) / (1 - swr - sor)
        varying_parameters["nsw"] = nsw
        varying_parameters["nso"] = nso
        ## recomputing mobilities
        assert (
            self.polymer.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPolymerConcentrationMatrix"
        )
        lambda_a = self.water.compute_mobility(
            c=self.polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=True,
            rel_permeability_formula=self.relative_permeability_formula,
        )
        lambda_o = self.water.compute_mobility(
            c=self.polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=False,
            rel_permeability_formula=self.relative_permeability_formula,
        )
        lambda_total = lambda_a + lambda_o
        varying_parameters["mobility_parameters"] = {
            "lambda_a": lambda_a,
            "lambda_o": lambda_o,
            "lambda_total": lambda_total,
        }
        ## fractional flow calculations
        f = lambda_a / lambda_total
        assert self.KK is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPermeabilityTensor"
        )
        D = self.KK * lambda_o * f
        varying_parameters["fractional_flow_parameters"] = {"f": f, "D": D}
        ## calculating dσ/dΓ
        varying_parameters["IFT_and_derivative"] = {
            "IFT": self.surfactant.eval_IFT,
            "dIFT_dg": self.surfactant.eval_dIFT_dGamma,
        }
        ## compute capillary number
        assert self.water.viscosity_array is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownAqueousViscosityMatrix"
        )
        nca = (
            np.sqrt(np.matmul(self.u, self.u) + np.matmul(self.v, self.v), dtype=np.complex128)
            * self.water.viscosity_array.astype(np.complex128)
            / self.surfactant.eval_IFT.astype(np.complex128)
        )
        nco = (
            np.sqrt(np.matmul(self.u, self.u) + np.matmul(self.v, self.v), dtype=np.complex128)
            * SimulationConstants.Oil_Viscosity.value
            / self.surfactant.eval_IFT
        )
        norm_nca = np.linalg.norm(nca)  # 2-norm of nca
        norm_nco = np.linalg.norm(nco)  # 2-norm of nco
        ## compute derivatives of residual saturations with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        dswr_dg = np.zeros((n, m))
        dsor_dg = np.zeros((n, m))
        swr_0 = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
        sor_0 = SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value
        assert (
            self.surfactant.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownSurfactantConcentrationMatrix"
        )
        for j in range(n):
            for i in range(m):
                if norm_nca >= const_parameters["resid_saturation_constants"]["Nca0"]:
                    dswr_dg[j, i] = -(
                        swr_0
                        * 0.1534
                        * 10.001
                        * (
                            const_parameters["resid_saturation_constants"]["Nca0"]
                            ** 0.1534
                        )
                    ) / (
                        (
                            np.sqrt((self.u[j, i] ** 2) + (self.v[j, i] ** 2))
                            * self.water.viscosity_array[j, i]
                        )
                        ** (0.1534)
                        * self.surfactant.eval_IFT[j, i] ** (0.8466)
                        * (self.surfactant.concentration_matrix[j, i] + 1) ** 2
                    )
                if norm_nco >= const_parameters["resid_saturation_constants"]["Nco0"]:
                    dsor_dg[j, i] = -(
                        sor_0
                        * 0.5213
                        * 10.001
                        * (
                            const_parameters["resid_saturation_constants"]["Nca0"]
                            ** 0.5213
                        )
                    ) / (
                        (
                            np.sqrt((self.u[j, i] ** 2) + (self.v[j, i] ** 2))
                            * self.water.viscosity_array[j, i]
                        )
                        ** (0.5213)
                        * self.surfactant.eval_IFT[j, i] ** (0.4787)
                        * (self.surfactant.concentration_matrix[j, i] + 1) ** 2
                    )
        varying_parameters["resid_saturation_derivatives"] = {
            "dswr_dg": dswr_dg,
            "dsor_dg": dsor_dg,
        }
        ## compute derivatives of normalized saturation with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["normalized_saturation_derivatives"] = {
            "dnsw_dg": dswr_dg * (self.water.water_saturation - 1) / (1 - swr) ** 2,
            "dnso_dg": (
                dswr_dg * (self.water.water_saturation + sor - 1)
                + dsor_dg * (self.water.water_saturation - swr)
            )
            / (1 - swr - sor) ** 2,
        }
        ## computing relative permeability with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["relative_permeability_derivatives"] = {
            "dkra_dg": 2.5 * dswr_dg * (nsw**3 - nsw)
            + (self.water.water_saturation - 1)
            * (2.5 * swr * (3 * nsw**2 - 1) + 1)
            * varying_parameters["normalized_saturation_derivatives"]["dnsw_dg"]
            / (1 - swr) ** 2,
            "dkro_dg": 1
            - 5 * sor * nso
            + (1 - nso) * (1 - 5 * nso * dsor_dg)
            - (1 + 5 * sor - 10 * sor * nso)
            * varying_parameters["normalized_saturation_derivatives"]["dnso_dg"],
            "dkra_ds": 2.5 * swr * (3 * (nsw) ** 2 - 1) + 1,
            "dkro_ds": 10 * sor * nso - 5 * sor - 1,
        }
        ## computing capillary pressure derivatives with respect to concentrations and saturations (FIXME: Need to update when inplementing autodiff!)
        pc = (
            self.surfactant.eval_IFT
            * const_parameters["Pc_constants"]["omega2"]
            * np.sqrt(const_parameters["porosity"])
        ) / (
            np.matmul(
                self.KK ** (0.5),
                fractional_matrix_power(
                    1 - nso, 1 / const_parameters["Pc_constants"]["omega1"]
                ),
            )
        )
        dpc_ds = pc / (const_parameters["Pc_constants"]["omega1"] * (1 - nso))
        dpc_dg = (
            pc / self.surfactant.eval_IFT
        ) * self.surfactant.eval_dIFT_dGamma + dpc_ds
        varying_parameters["capillary_pressure_and_derivatives"] = {
            "pc": pc,
            "dpc_ds": dpc_ds,
            "dpc_dg": dpc_dg,
        }
        ## computing fractional flow derivatives with respect to concentrations and saturations (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["fractional_flow_derivatives"] = {
            "df_ds": varying_parameters["relative_permeability_derivatives"]["dkra_ds"]
            * lambda_o
            / (lambda_total**2 * self.water.viscosity_array)
            - varying_parameters["relative_permeability_derivatives"]["dkro_ds"]
            * lambda_a
            / (lambda_total**2 * self.water.miuo),
            "df_dc": (-1 * (lambda_o * lambda_a * self.water.miuo))
            / ((lambda_total**2) * self.water.viscosity_array),
            "df_dg": (
                (
                    varying_parameters["relative_permeability_derivatives"]["dkra_dg"]
                    * lambda_o
                )
                / ((lambda_total**2) * self.water.viscosity_array)
            )
            - (
                (
                    varying_parameters["relative_permeability_derivatives"]["dkro_dg"]
                    * lambda_a
                )
                / ((lambda_total**2) * self.water.viscosity_array)
            ),
            "dD_dg": D
            * varying_parameters["capillary_pressure_and_derivatives"]["dpc_dg"],
            "dD_ds": D
            * varying_parameters["capillary_pressure_and_derivatives"]["dpc_ds"],
        }
        # Update Water Saturation Matrix
        ## Calculate ``xmod`` and ``ymod``
        [xmod, ymod] = self._characteristic_coordinates(
            1,
            self.water.water_saturation,
            self.water.water_saturation,
            const_parameters,
            varying_parameters,
        )

        ## Pass in parameters into ``compute_water_saturation`` method of the ``Water`` class
        Q_old = np.copy(self.water.water_saturation)
        Qmod, varying_parameters = self.water.compute_water_saturation(
            grid=self.mesh,
            surfactant=self.surfactant,
            polymer=self.polymer,
            u=self.u,
            v=self.v,
            xmod=xmod,
            ymod=ymod,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
        )

        # Update the Polymer Concentration Matrix
        [xmod, ymod] = self._characteristic_coordinates(
            2, Q_old, self.water.water_saturation, const_parameters, varying_parameters
        )
        C_old = np.copy(self.polymer.concentration_matrix)
        varying_parameters = self.polymer.compute_concentration(
            grid=self.mesh,
            water_sat=self.water.water_saturation,
            u=self.u,
            v=self.v,
            xmod=xmod,
            ymod=ymod,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
        )

        # Update the Surfactant Concentration matrix
        [xmod, ymod] = self._characteristic_coordinates(
            3, Q_old, self.water.water_saturation, const_parameters, varying_parameters
        )
        G_old = np.copy(self.surfactant.concentration_matrix)
        G = np.copy(self.surfactant.concentration_matrix)
        x1d = self.mesh.x[0, :]
        y1d = self.mesh.y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)

        # reorder surfactant.concentration_matrix if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            G = G[:, x_sort_idx]  # Sort columns of surfactant.concentration_matrix
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            G = G[y_sort_idx, :]  # Sort rows of surfactant.concentration_matrix

        interp = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d),
            G,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )
        query_points = np.stack([ymod.ravel(), xmod.ravel()], axis=-1)
        Gmod = interp(query_points).reshape(xmod.shape)

        # Updating coefficients using interpolated surfactant concentration
        assert self.surfactant.IFT_conc_equ is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownIFTEquation"
        )
        assert (
            self.surfactant.derivative_IFT_conc_equ is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownDerivativeIFTEquation"
        )
        sigma_mod = self.surfactant.IFT_conc_equ(Gmod)
        sigma_g_mod = self.surfactant.derivative_IFT_conc_equ(Gmod)
        [swr, sor] = self.water.compute_residual_saturations(
            sigma=sigma_mod, u=self.u, v=self.v
        )
        lambda_a = self.water.compute_mobility(
            c=C_old,
            sor=float(sor),
            swr=float(swr),
            aqueous=True,
            rel_permeability_formula=self.relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_o = self.water.compute_mobility(
            c=C_old,
            sor=float(sor),
            swr=float(swr),
            aqueous=False,
            rel_permeability_formula=self.relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_total = lambda_a + lambda_o
        varying_parameters["mobility_parameters"] = {
            "lambda_a": lambda_a,
            "lambda_o": lambda_o,
            "lambda_total": lambda_total,
        }
        f = lambda_a / lambda_total
        assert self.KK is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPermeabilityTensor"
        )
        D = self.KK * lambda_o * f
        varying_parameters["fractional_flow_parameters"] = {"f": f, "D": D}
        dpc_ds = pc / (const_parameters["Pc_constants"]["omega1"] * (1 - nso))
        dpc_dg = (pc / sigma_mod) * sigma_g_mod + dpc_ds
        varying_parameters["capillary_pressure_and_derivatives"] = {
            "pc": pc,
            "dpc_ds": dpc_ds,
            "dpc_dg": dpc_dg,
        }
        F = D * dpc_dg / self.water.water_saturation
        varying_parameters = self.surfactant.compute_concentration(
            grid=self.mesh,
            water_sat=self.water.water_saturation,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
            F=F,
            Gmod=Gmod,
        )

        # Returning updated Water, Polymer, and Surfactant objects
        return self.water, self.polymer, self.surfactant

    def _characteristic_coordinates(
        self,
        flag,
        old_water_saturation_matrix,
        new_water_saturation_matrix,
        const_parameters,
        varying_parameters,
    ):
        """
        (private method)

        Compute redefined characteristic coordinates (xmod, ymod) according to Neumann boundary conditions.

        will be a helper function to the ``self._transport_equation_solver()`` method.

        Args:
        ------
            flag (int): scenario flag for the simulation the user wants to run

            old_water_saturation_matrix (np.ndarray): water saturation matrix from previous iteration

            new_water_saturation_matrix (np.ndarray): water saturation from current iteration

            const_parameters (dict): constant parameters to help with calculations

            varying_parameters (dict): parameters that vary but assist with calculations for water saturation, polymer concentration, and surfactant concentration

        Returns: (tuple[np.ndarray, np.ndarray])
        ---------------------------------------
            This method returns ``xmod`` and ``ymod``, which are the modified characteristic coordinates according to the Neumann boundary conditions
        """
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownWaterSaturationMatrix"
        )
        assert (
            self.surfactant.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownSurfactantConcentrationMatrix"
        )
        x, y = (
            const_parameters["FD_grid_constants"]["x"],
            const_parameters["FD_grid_constants"]["y"],
        )
        dt_matrix = const_parameters["FD_grid_constants"]["dt_matrix"]
        xjump = None
        yjump = None
        f = varying_parameters["fractional_flow_parameters"]["f"]
        f_s = varying_parameters["fractional_flow_derivatives"]["df_ds"]
        D = varying_parameters["fractional_flow_parameters"]["D"]
        pc_s = varying_parameters["capillary_pressure_and_derivatives"]["dpc_ds"]
        pc_g = varying_parameters["capillary_pressure_and_derivatives"]["dpc_dg"]
        sold = old_water_saturation_matrix
        snew = new_water_saturation_matrix

        if flag == 1:
            xjump = x - f_s * self.u * dt_matrix
            yjump = y - f_s * self.v * dt_matrix
        elif flag == 2:
            # Calculate gradients
            sx, sy = self._get_gradient(sold)
            gx, gy = self._get_gradient(self.surfactant.concentration_matrix)

            xjump = (
                x
                - (
                    (f / snew) * self.u
                    + (D * pc_s / snew) * sx
                    + (D * pc_g / snew) * gx
                )
                * dt_matrix
            )
            yjump = (
                y
                - (
                    (f / snew) * self.v
                    + (D * pc_s / snew) * sy
                    + (D * pc_g / snew) * gy
                )
                * dt_matrix
            )
        elif flag == 3:
            sx, sy = self._get_gradient(sold)

            xjump = x - ((f / snew) * self.u + (D * pc_s / snew) * sx) * dt_matrix
            yjump = y - ((f / snew) * self.v + (D * pc_s / snew) * sy) * dt_matrix

        # Apply Neumann reflection conditions
        if xjump is None or yjump is None:
            raise SimulationCalcInputException(
                "SimulationInputException:UnknownXJumpYJumpMatrices"
            )

        xmod = np.copy(x)
        ymod = np.copy(y)

        for j in range(np.shape(y)[0]):
            for i in range(np.shape(x)[1]):
                if xjump[j, i] <= 1 and yjump[j, i] <= 1:
                    xmod[j, i] = np.abs(xjump[j, i])
                    ymod[j, i] = np.abs(yjump[j, i])
                elif xjump[j, i] > 1 and yjump[j, i] <= 1:
                    xmod[j, i] = 2 - xjump[j, i]
                    ymod[j, i] = np.abs(yjump[j, i])
                elif xjump[j, i] <= 1 and yjump[j, i] > 1:
                    xmod[j, i] = np.abs(xjump[j, i])
                    ymod[j, i] = 2 - yjump[j, i]
                elif xjump[j, i] > 1 and yjump[j, i] > 1:
                    xmod[j, i] = 2 - xjump[j, i]
                    ymod[j, i] = 2 - yjump[j, i]

        return xmod, ymod

    def _export_results(self):
        """
        (private method)

        Saves simulation sim_results to CSV files.
        """
        import os

        os.makedirs("sim_results", exist_ok=True)

        np.savetxt("sim_results/COC.csv", self.COC, delimiter=",")
        if hasattr(self, "lambdaTcal"):
            np.savetxt(
                "sim_results/lambdaTcal.csv", np.array(self.lambdaTcal), delimiter=","
            )
        if hasattr(self, "miuaTcal"):
            np.savetxt(
                "sim_results/miuaTcal.csv", np.array(self.miuaTcal), delimiter=","
            )

        print("Simulation sim_results exported to /sim_results/ folder.")
        
    def _compute_MFW(self, UU):
        """
        (private function)

        Args:
        -----
            UU (np.ndarray): water saturation matrix

        Returns: (np.array)
        -------
            list of values which are the mean finger width during each iteration.
            Note: Only will run under the Rectilinear Homogenous and Rectilinear Heterogeneous simulation scenarios
        """
        # post processing of finger width 
        interface = np.zeros((29, 1))
        mean_UU_save = np.zeros((29, 29))
        store_UU = 0
        
        # find average concentration for each Y level
        iter_x = 0
        iter_y = 0
        iter_counter = 0
        
        while iter_y < 29:
            while iter_x < 29:
                
                if UU[iter_y, iter_x] > 0.21 and UU[iter_y, iter_x] < 0.99:
                    iter_counter += 1
                    store_UU += UU[iter_y, iter_x]
                    
                iter_x += 1
                
            mean_UU_save.T.flat[iter_y] = store_UU / iter_counter
            iter_counter = 0
            store_UU = 0
            iter_y += 1
            iter_x = 0
        
        iter_y = 0
        iter_x = 0
        while iter_y < 29:
            while iter_x < 29:
                if UU[iter_y, iter_x] < mean_UU_save.T.flat[iter_y]:
                    interface[iter_y, 0] = iter_x
                    break
                iter_x += 1
                
            iter_x = 0
            iter_y += 1
        
        # find location of interface front
        iter_x = 0
        iter_y = 0
        while iter_x < 29:
            if np.mean(UU[:, iter_x]) < np.mean(mean_UU_save[:, 0]):
                break
            iter_x += 1
        
        # find presence of saturation along the mixing layer
        rows = UU.shape[0]
        cols = UU.shape[1]
        check_concentration = np.zeros((rows, cols))
        for ii in range(rows):
            for jj in range(cols):
                if UU[ii, jj] < mean_UU_save[0, iter_y]:
                    check_concentration[ii, jj] = 1
                else:
                    check_concentration[ii, jj] = 0
        
        last = 0
        counter = 0
        mean_finger_width = 0
        total_concentration = 0
        jj = iter_x - 1
        for i in range(1):
            for ii in range(cols):
                if check_concentration[ii, jj] == 1:
                    new_last = 1
                    total_concentration += 1
                else:
                    new_last = 0
                
                if (new_last == 1 and last == 0) or (new_last == 0 and last == 1):
                    counter += 1
                last = new_last
        
            old_MFW = mean_finger_width
            mean_finger_width = 2 * total_concentration / (counter + 1)
            mean_finger_width = np.max(mean_finger_width, old_MFW)
            jj += 1
            
        return interface, mean_finger_width, iter_x

    # Public Method of Simulation Class
    def run(self):
        """
        Executes simulation loop.
        
        Raises:
        ------
            SimulationCalcInputException: If there is an issue with the execution of a calculation during runtime
        """
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:WaterSaturationMatrixUnavailable"
        )
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        try:
            ## STEP 1: initializing start time, end time, and time step
            t = 0
            t_cal = 0
            t_stop = 500
            dt = self.mesh.dx / self.source_flow_magnitude
            # dt value if running under Quarter Five Spot & Heterogeneous scenario:
            if bool_Heterogenous_and_Quarter_Five_Spot:
                dt *= 100

            ## STEP 2: Initiating the primary 'while' loop that will keep running until water shows up in production well
            while(t < t_stop and self.water.water_saturation[self.mesh.n, self.mesh.m] <= 0.70):
                print(f'{t},{self.water.water_saturation[self.mesh.n, self.mesh.m]}')
                # while t < 1:
                ## STEP 2.1: Increment time and amount of feed used:
                self.integrated_inlet_flow += self.source_flow_magnitude
                t += dt
                ## STEP 2.2: Compute viscosities:
                if (
                    self.model_type.value == ModelType.No_Shear_Thinning.value
                ):  # if No Polymer Shear Thinning
                    self.water.compute_viscosity(
                        grid=self.mesh,
                        model_type=self.model_type,
                        polymer=self.polymer,
                        u=self.u,
                        v=self.v,
                    )
                    self.polymer.compute_viscosity(
                        grid=self.mesh,
                        u=self.u,
                        v=self.v,
                        model_type=self.model_type,
                        aqueous_viscosity=self.water.viscosity_array,
                    )
                elif self.model_type.value == ModelType.Shear_Thinning_On.value:
                    self.polymer.compute_viscosity(
                        grid=self.mesh,
                        u=self.u,
                        v=self.v,
                        model_type=self.model_type,
                        aqueous_viscosity=None,
                    )
                    self.water.compute_viscosity(
                        grid=self.mesh,
                        model_type=self.model_type,
                        polymer=self.polymer,
                        u=self.u,
                        v=self.v,
                    )
                ## STEP 2.2: Computing Residual Saturation:
                assert (
                    self.surfactant.IFT_conc_equ is not None
                ), SimulationCalcInputException(
                    "SimulationCalcInputError:SurfactantIFTEquationUnavailable"
                )
                interfacial_tension_matrix = self.surfactant.IFT_conc_equ(
                    self.surfactant.concentration_matrix
                )
                [resid_water_saturation, resid_oleic_saturation] = (
                    self.water.compute_residual_saturations(
                        sigma=interfacial_tension_matrix, u=self.u, v=self.v
                    )
                )
                ## STEP 2.3: Compute mobilities:
                assert (
                    self.polymer.concentration_matrix is not None
                ), SimulationCalcInputException(
                    "SimulationCalcInputError:PolymerConcentrationMatrixUnavailable"
                )
                aqueous_mobility = self.water.compute_mobility(
                    c=self.polymer.concentration_matrix,
                    sor=float(resid_oleic_saturation),
                    swr=float(resid_water_saturation),
                    aqueous=True,
                    rel_permeability_formula=self.relative_permeability_formula,
                )
                oleic_mobility = self.water.compute_mobility(
                    c=self.polymer.concentration_matrix,
                    sor=float(resid_oleic_saturation),
                    swr=float(resid_water_saturation),
                    aqueous=False,
                    rel_permeability_formula=self.relative_permeability_formula,
                )
                total_mobility = aqueous_mobility + oleic_mobility
                assert self.KK is not None, SimulationCalcInputException(
                    "SimulationCalcInputError:PermeabilityTensorUnavailable"
                )
                beta = self.KK * total_mobility

                ## STEP 2.4: Calculating Global Pressure and velocity
                ### STEP 2.4.1: setting FEM Mesh
                self.FE_mesh.set_triangulation()
                self.FE_mesh.set_FE_meshgrid(beta)
                self.FE_mesh.set_right_hand(self.source_prod_flow)
                self.FE_mesh.get_A_B_matrices()
                ### STEP 2.4.2: updating the pressure & velocity matrices
                u_old = self.u  # storing old pressure matrix
                v_old = self.v  # storing old velocity matrix
                self.u, self.v = self._compute_pressure_and_velocity_matrices(
                    self.FE_mesh.sparsed_A, self.FE_mesh.B, beta
                )

                ## STEP 2.5: Solving Transport Equations
                self._transport_equation_solver(dt)

                ## Step 2.6: MFW post processing (excluding QFS)
                if (self.scenario_flag != 3): # FIXME: compute_MFW currently operates for rectilinear geometries. Implement MFW computation for QFS
                    interface, MFW_val, _ = self._compute_MFW(self.water.water_saturation)
                    self.MFW.append(MFW_val)
        except Exception as e:
            print(e)
