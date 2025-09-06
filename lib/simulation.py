"""
This python script contains the class definition for running simulations

This Python code has been derived from the MATLAB Surfactant-Polymer Flooding Simulation 
developed by Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar and Carlos Acosta Caripo

"""

import os

import numpy as np
from grid import Grid, FEMesh
from enumerations import (
    ModelType,
    PolymerList,
    RelativePermeabilityFormula,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from Exceptions import SimulationCalcInputException, UserInputException
from polymer import Polymer
from surfactant import Surfactant
from water import Water
from scipy.io import loadmat
from scipy.sparse.linalg import bicgstab

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

        :raises UserInputException: If there is a issue with the user inputs in ``user_input_dict``
        :raises SimulationCalcInputException: If there is an issue with the execution of a calculation during runtime

        :param user_input_dict: dictionary containing the information from the GUI
        :type user_input_dict: dict
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
        self.MFW = (
            []
        )  # Mean Finger Width (will be converted into a numpy array when reporting)
        self.integrated_inlet_flow = 0  # "src_total" in the MATLAB version of the code

    # Dependent Property of Simulation Class
    _source_prod_flow = None

    @property
    def source_prod_flow(self):
        """
        Return
        ------
        :return: returns the matrix with the source & and production well flow rates
        :rtype: np.ndarray
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

    def _initialize_pressure_and_velocity(self):
        """
        (private method)

        Initializing global pressure ('u') and velocity matrices ('v')
        Will use the ``n`` and ``m`` properties from ``Grid`` Class for initialization

        Return
        ------
        :return: the global pressure matrix (index 0) and velocity matrix (index 1)
        :rtype: list[np.ndarray]
        """
        u = np.zeros((self.mesh.n + 1, self.mesh.m + 1))
        v = np.zeros((self.mesh.n + 1, self.mesh.m + 1))

        return u, v

    def _compute_pressure_and_velocity_matrices(self, sparsed_A, B, beta):
        """
        (private method)

        dependent property to calculate the pressure matrix (``u``)
        and velocity matrix (``v``). Will rely on functions in the ``FEMesh`` class.

        :return: list of update pressure matrix and velocity matrix => [u,v]
        :rtype: list[np.ndarray]
        """
        max_iterations = 1000
        new_u, convergence_flag = bicgstab(
            sparsed_A, B, maxiter=max_iterations
        )  # new_u is of shape (900, )
        assert convergence_flag == 0, SimulationCalcInputException("ConvergenceFailure")

        new_v = np.zeros((self.FE_mesh.n + 1, self.FE_mesh.m + 1), dtype=object)
        for i in range(self.FE_mesh.m + 1):
            for j in range(self.FE_mesh.n + 1):
                new_v[j, i] = new_u[j * (self.FE_mesh.m + 1) + i]

        px, py = self._get_gradient(new_v)
        new_u = (-1 * beta) * px
        new_v = (-1 * beta) * py

        return new_u, new_v

    def _get_gradient(self, vn):
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

        :return: Initialized ``ProdRate`` and ``CROIP`` properties
        :rtype: list
        """
        os.makedirs("memmaps", exist_ok=True)
        tf = 500
        dt = self.mesh.dx / self.source_flow_magnitude
        timestamps = int(np.floor(tf / dt))
        CROIP = np.memmap(
            "memmaps/ProdRate.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )
        ProdRate = np.memmap(
            "memmaps/ProdRate.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )

        return ProdRate, CROIP

    def _create_mesh(self):
        """
        (private method)

        :return: Initialized FD and FE mesh
        :rtype: Tuple[Grid, FEMesh]
        """
        FD_mesh = Grid(self.grid_size, self.grid_size)

        FE_mesh = FEMesh(self.grid_size, self.grid_size)

        return FD_mesh, FE_mesh

    # def _generate_grid(self):
    #     x = np.arange(self.mesh.left, self.mesh.right + self.mesh.dx, self.mesh.dx)
    #     y = np.arange(self.mesh.bottom, self.mesh.top + self.mesh.dy, self.mesh.dy)
    #     return np.meshgrid(x, y)

    def _initialize_simulation(self):
        """
        (private method)

        Sets up initial reservoir fields, permeability, and time step.

        :return: initialized properties of simulation. Required in ``__init__`` function
        :rtype: None
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

        Compute level set function phi at each grid point.
        Equivalent to MATLAB get_phi_test function.
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

        Compute the initial position of the water front.
        Equivalent to MATLAB z_func_test.

        Takes array user_input_dict x, y.
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

        Compute permeability matrix KK based on the flag.
        Equivalent to MATLAB KKdef function.
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
        # elif flag == 3:
        #     # Impermeable block at center
        #     KK = 3000 * np.ones((m + 1, m + 1))
        #     center = m // 2
        #     delta = m // 8
        #     KK[
        #         center - delta : center + delta + 1,
        #         center - delta : center + delta + 1,
        #     ] = 3
        # elif flag == 4:
        #     # Impermeable blocks off-center
        #     KK = 3000 * np.ones((m + 1, m + 1))
        #     KK[
        #         (3*m)//4 - m//12 : (3*m)//4 + m//12 + 1,
        #         (2*m)//3 - m//12 : (2*m)//3 + m//12 + 1
        #     ] = 3
        #     KK[
        #         m//3 - m//10 : m//3 + m//10 + 1,
        #         m//3 - m//10 : m//3 + m//10 + 1
        #     ] = 3
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
        # elif flag == 6:
        #     # Load Tarbert formation (SPE10)
        #     mat_data = loadmat('/Resources/KK30Tabert.mat')
        #     if 'KK' not in mat_data:
        #         raise SimulationCalcInputException('SimulationInputException: KK matrix not found in KK30Tabert.mat file.')
        #     KK = mat_data['KK']
        # else:
        #     raise SimulationCalcInputException("SimulationInputException: Unknown permeability flag.")
        return KK

    def _characteristic_coordinates(
        self, x, y, s, snew, g, f, f_s, D, pc_s, pc_g, u, v, dt, para, flag
    ):
        """
        (private method)

        Compute redefined characteristic coordinates (xmod, ymod) according to Neumann boundary conditions.
        """
        dx, dy = para.box.dx, para.box.dy
        m, n = para.box.m, para.box.n
        xjump = None
        yjump = None
        if flag == 1:
            xjump = x - f_s * u * dt
            yjump = y - f_s * v * dt
        elif flag == 2:
            # Calculate gradients
            sx, sy = np.gradient(s, dx, dy, edge_order=2)
            gx, gy = np.gradient(g, dx, dy, edge_order=2)

            xjump = (
                x
                - ((f / snew) * u + (D * pc_s / snew) * sx + (D * pc_g / snew) * gx)
                * dt
            )
            yjump = (
                y
                - ((f / snew) * v + (D * pc_s / snew) * sy + (D * pc_g / snew) * gy)
                * dt
            )
        elif flag == 3:
            sx, sy = np.gradient(s, dx, dy, edge_order=2)

            xjump = x - ((f / snew) * u + (D * pc_s / snew) * sx) * dt
            yjump = y - ((f / snew) * v + (D * pc_s / snew) * sy) * dt

        # Apply Neumann reflection conditions
        if xjump is None or yjump is None:
            raise SimulationCalcInputException(
                "SimulationInputException: xjump or yjump not initialized..."
            )

        xmod = np.where(xjump <= 1, np.abs(xjump), 2 - xjump)
        ymod = np.where(yjump <= 1, np.abs(yjump), 2 - yjump)

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

    def run(self):
        """
        Executes simulation loop.

        :raises SimulationCalcInputException: if relevant inputs for calculation not provided or not initialized

        :return: Dictionary with relevant results for plotting and data analysis
        :rtype: dict
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
            # while(t < t_stop and self.water.water_saturation[self.mesh.n, self.mesh.m] <= 0.70):
            while t < 1:
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
                    surfactant_conc=self.surfactant.concentration,
                )
                oleic_mobility = self.water.compute_mobility(
                    c=self.polymer.concentration_matrix,
                    sor=float(resid_oleic_saturation),
                    swr=float(resid_water_saturation),
                    aqueous=False,
                    rel_permeability_formula=self.relative_permeability_formula,
                    surfactant_conc=self.surfactant.concentration,
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

                u_old = self.u
                v_old = self.v
                self.u, self.v = self._compute_pressure_and_velocity_matrices(
                    self.FE_mesh.sparsed_A, self.FE_mesh.B, beta
                )

                break

        except Exception as e:
            print(e)
