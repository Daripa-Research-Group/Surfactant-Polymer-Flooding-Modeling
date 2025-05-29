"""
This python script contains the class definition for running simulations
"""
import os

import numpy as np
from para import Box
from enumerations import (
    ModelType,
    PolymerList,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from Exceptions import SimulationCalcInputException, UserInputException
from polymer import Polymer
from surfactant import Surfactant
from scipy.io import loadmat
os.makedirs("memmaps", exist_ok=True)


class Simulation:
    """
    Simulation class to run SP-flooding simulations based on MATLAB translation.
    """

    def __init__(
            self,
            user_input_dict : dict
            ):
        """
        This method will check the user user_input_dict and initialize the simulation
        """
        ## Performs checks on the user input dictionary passed in:
        try:
            model_type = ModelType(user_input_dict["model_type"])
            if model_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Model Type not selected. Please try again.", user_input_dict)

        try:
            reservoir_geometry = ResevoirGeometry(user_input_dict["reservoir_geometry"])
            if reservoir_geometry is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Reservoir Geometry not selected. Please try again.", user_input_dict)

        try:
            permeability_flag = PermeabilityType(user_input_dict["permeability"])
            if permeability_flag is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Permeability not properly selected. Please try again.", user_input_dict)

        try:
            polymer_type = PolymerList.get_by_value(user_input_dict["polymer_type"])
            if polymer_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Polymer not properly selected. Please try again.", user_input_dict)

        try:
            polymer_concentration = user_input_dict["polymer_concentration"]
            if polymer_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Polymer concentration not given. Please try again.", user_input_dict)

        try:
            surfactant_type = SurfactantList(user_input_dict["surfactant_type"])
            if surfactant_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Surfactant not properly selected. Please try again.", user_input_dict)

        try:
            surfactant_concentration = user_input_dict["surfactant_concentration"]
            if surfactant_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException("Surfactant concentration not given. Please try again.", user_input_dict)
        
        ## Instantiates the required simulation properties:
        
        # Initializing Simulation Flags:
        self.permeability_flag = permeability_flag
        self.reservoir_geometry = reservoir_geometry
        self.model_type = model_type
        
        #initializing grid size
        self.grid_size = SimulationConstants.Grid_Size.value

        # Initializing sim properties
        self.mesh = self._create_mesh()
        self.x, self.y = self._generate_grid()
        self.phi = None  # Level set function (relates to porosity)
        self.KK = None  # Permeability tensor
        self.time_step = None
        self.initialize_simulation() #will initialize phi, KK, and time_step
        if(self.phi is None or self.KK is None or self.time_step is None): #Raise Exception if not properly initialized...
            raise SimulationCalcInputException("SimulationInputException: phi, KK, or time_step not properly initialized. Please try again...") 
        
        # Initalizing Polymer Object
        self.polymer = Polymer(
            name = polymer_type,
            e_coeff = polymer_type.n_coeff,
            n_coeff = polymer_type.e_coeff,
            rho = polymer_type.Density,
            concentration_scalar = polymer_concentration,
            phi = self.phi,
        )
        # Initializing Surfactant Object
        self.surfactant = Surfactant(    
            name = surfactant_type,
            initial_concentration = surfactant_concentration,
            IFT_equation = surfactant_type.IFT_equation,
            derivative_IFT_equation =  surfactant_type.derivative_IFT_equation,
            phi = self.phi,
        )
        # TODO: Initializing Water Object








    # def __init__(
    #     self,
    #     sim_id: int,
    #     grid_size: int,
    #     polymer: Polymer,
    #     surfactant: Surfactant,
    #     reservoir_geometry: ResevoirGeometry,
    #     permeability_flag: PermeabilityType,
    #     model_type: ModelType,
    #     plot_type,
    #     init_water_saturation=SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value,
    #     init_oleic_saturation=SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value,
    #     source_flow_magnitude=SimulationConstants.Source_Flow_Magnitude.value,
    # ):
    #     self.sim_id = sim_id
    #     self.grid_size = grid_size
    #     self.polymer = polymer
    #     self.surfactant = surfactant
    #
    #     self.reservoir_geometry = reservoir_geometry
    #     self.permeability_flag = permeability_flag
    #     self.model_type = model_type
    #     self.plot_type = plot_type
    #
    #     self.init_water_saturation = init_water_saturation
    #     self.init_oleic_saturation = init_oleic_saturation
    #     self.source_flow_magnitude = source_flow_magnitude
    #
    #     # Internal initialization
    #     self.mesh = self._create_mesh()
    #     self.x, self.y = self._generate_grid()
    #     self.phi = None  # Level set function
    #     self.KK = None  # Permeability field
    #     self.time_step = None
    #
    #     # Initial field variables
    #     self.s0 = 0.79  # Initial oil saturation
    #     self.c0 = self.polymer.initial_concentration
    #     self.g0 = self.surfactant.concentration if surfactant else 0
    #
    #     self.U = None  # Water saturation
    #     self.C = None  # Polymer concentration
    #     self.G = None  # Surfactant concentration
    #     self.c0_array = None
    #     self.miuw = SimulationConstants.Water_Viscosity.value
    #     self.miuo = SimulationConstants.Oil_Viscosity.value
    #     self.miup = None
    #     self.miup_array = None

    def _create_mesh(self):
        mesh = Box()
        mesh.m = self.grid_size
        mesh.n = self.grid_size
        mesh.calculate_spacing
        return mesh

    def _generate_grid(self):
        x = np.arange(self.mesh.left, self.mesh.right + self.mesh.dx, self.mesh.dx)
        y = np.arange(self.mesh.bottom, self.mesh.top + self.mesh.dy, self.mesh.dy)
        return np.meshgrid(x, y)

    def initialize_simulation(self):
        """
        Sets up initial reservoir fields, permeability, and time step.
        """
        self.phi = self._compute_phi()

        # Initialize permeability matrix
        self.KK = self._compute_permeability()

        # Calculate time step
        self.time_step = self.mesh.dx / self.source_flow_magnitude
        if self.permeability_flag == PermeabilityType.Heterogenous and self.reservoir_geometry == ResevoirGeometry.Quarter_Five_Spot:
            self.time_step *= 100

    def _compute_phi(self):
        """
        Compute level set function phi at each grid point.
        Equivalent to MATLAB get_phi_test function.
        """
        m = self.mesh.m
        n = self.mesh.n
        dx = self.mesh.dx
        dy = self.mesh.dy
        left = self.mesh.left
        bottom = self.mesh.bottom

        jj, ii = np.meshgrid(np.arange(1, n+2), np.arange(1, m+2))
        x_coords = left + (ii - 1) * dx
        y_coords = bottom + (jj - 1) * dy

        phi = self._z_func_test(x_coords, y_coords)
        return phi
    
    def _z_func_test(self, x, y):
        """
        Compute the initial position of the water front.
        Equivalent to MATLAB z_func_test.
        
        Takes array user_input_dict x, y.
        """
        init_front_hs = 0.1
        permeability_input = self.permeability_flag.value

        if permeability_input == 1:
            # Homogeneous
            out = y - init_front_hs + 0.01 * np.cos(80 * np.pi * x)
        elif permeability_input == 2:
            # Rectilinear Heterogeneous
            out = y - init_front_hs
        elif permeability_input == 5:
            # Quarter five spot 
            out = np.square(x) + np.square(y) - 0.015
        else:
            out = np.zeros_like(x)
        
        return out

    def _compute_permeability(self):
        """
        Compute permeability matrix KK based on the flag.
        Equivalent to MATLAB KKdef function.
        """
        flag = self.permeability_flag.value
        m = self.mesh.m

        if flag == 1:
            # Homogeneous
            Kmax = 100
            KK = Kmax * np.ones((m + 1, m + 1))

        elif flag == 2:
            # Heterogeneous
            Kmax = 100
            KK = Kmax * (
                0.5 * (1 - 10 ** (-7)) * (np.sin(6 * np.pi * np.cos(self.x)) * np.cos(4 * np.pi * np.sin(3 * self.y)) - 1) + 1
            )

        elif flag == 3:
            # Impermeable block at center
            KK = 3000 * np.ones((m + 1, m + 1))
            center = m // 2
            delta = m // 8
            KK[
                center - delta : center + delta + 1,
                center - delta : center + delta + 1,
            ] = 3

        elif flag == 4:
            # Impermeable blocks off-center
            KK = 3000 * np.ones((m + 1, m + 1))
            KK[
                (3*m)//4 - m//12 : (3*m)//4 + m//12 + 1,
                (2*m)//3 - m//12 : (2*m)//3 + m//12 + 1
            ] = 3
            KK[
                m//3 - m//10 : m//3 + m//10 + 1,
                m//3 - m//10 : m//3 + m//10 + 1
            ] = 3

        elif flag == 5:
            # Load Upper Ness formation (SPE10)
            mat_data = loadmat('/Resources/KK30Ness.mat')
            if 'KK' not in mat_data:
                raise Exception('KK matrix not found in KK30Ness.mat file.')
            KK = mat_data['KK']
        elif flag == 6:
            # Load Tarbert formation (SPE10)
            mat_data = loadmat('/Resources/KK30Tabert.mat')
            if 'KK' not in mat_data:
                raise Exception('KK matrix not found in KK30Tabert.mat file.')
            KK = mat_data['KK']

        else:
            raise Exception("Unknown permeability flag.")

        return KK
    
    def _compvis(self, c, u, v):
        """
        Computes aqueous viscosity and shear rate based on polymer concentration and velocity fields.
        """
        viscosity_flag = self.model_type.value
        miuw = self.miuw
        miuo = self.miuo
        beta1 = 15000  

        gamma_dot = np.zeros_like(c)
        n, m = c.shape

        if viscosity_flag == 1:
            # No-shear-thinning model
            if self.c0 == 0:
                miua = miuw * np.ones((n, m))
            else:
                miua = miuw * (1 + beta1 * c)

        elif viscosity_flag == 2:
            # Sourav's implementation
            if self.c0 == 0:
                miua = miuw * np.ones((n, m))
            else:
                miua = miuo * (0.5 + c)

        elif viscosity_flag == 3:
            # Shear-thinning (Dynamic viscosity model)
            rho_water = 1000
            rho_polymer = self.polymer.name.Density

            w1 = rho_polymer * c
            w2 = rho_water * (1 - c)
            wppm = (w1 / (w1 + w2)) * 1e6

            w10 = rho_polymer * self.c0_array
            w20 = rho_water * (1 - self.c0_array)
            wppm0 = (w10 / (w10 + w20)) * 1e6

            eps_coeff = self.polymer.e_coeff
            n_coeff = self.polymer.n_coeff

            epsilon0 = eps_coeff[0] * wppm0 ** eps_coeff[1]
            power_n0 = np.minimum(n_coeff[0] * wppm0 ** n_coeff[1], 1)

            epsilon = eps_coeff[0] * wppm ** eps_coeff[1]
            power_n = np.minimum(n_coeff[0] * wppm ** n_coeff[1], 1)

            miua = miuw * np.ones((n, m))

            # Compute divergence terms
            a1 = np.gradient(v, axis=0)
            a2 = np.gradient(u, axis=1)
            a3 = np.gradient(u, axis=0)
            a4 = np.gradient(v, axis=1)

            pi_D = np.abs(-0.25 * (a1 + a2) ** 2 + a3 * a4)

            for i in range(n):
                for j in range(m):
                    if c[i, j] > 0:
                        gamma_dot[i, j] = 2 * np.sqrt(pi_D[i, j])

                        if gamma_dot[i, j] != 0:
                            self.miup_array[i, j] = epsilon0[i, j] * (gamma_dot[i, j] ** (power_n0[i, j] - 1))
                            miua[i, j] = epsilon[i, j] * (gamma_dot[i, j] ** (power_n[i, j] - 1))

                            miua[i, j] = np.clip(miua[i, j], miuw, 100)
                            self.miup_array[i, j] = np.clip(self.miup_array[i, j], miuw, 100)

        else:
            raise Exception("Invalid viscosityFlag value.")

        return miua, gamma_dot
    
    def _compres(self, sigma, u, v, miua):
        """
        Computes residual saturations swr and sor as functions of capillary numbers.
        """
        miuo = self.miuo
        swr0 = self.init_water_saturation
        sor0 = self.init_oleic_saturation

        # Critical capillary numbers
        Nco0 = 1.44e-4
        Nca0 = 1.44e-4

        # Compute capillary numbers
        velocity_mag = np.sqrt(u ** 2 + v ** 2)
        nca = (velocity_mag * miua) / sigma
        nco = (velocity_mag * miuo) / sigma

        Nca = np.linalg.norm(nca)
        Nco = np.linalg.norm(nco)

        # Update residual saturations
        if Nco < Nco0:
            sor = sor0
        else:
            sor = sor0 * (Nco0 / Nco) ** 0.5213

        if Nca < Nca0:
            swr = swr0
        else:
            swr = swr0 * (Nca0 / Nca) ** 0.1534

        return swr, sor
    
    def _save_lambda_miua(self, tcal, lambda_a, miua):
        if (tcal % 200) == 0:
            if not hasattr(self, "lambdaTcal"):
                self.lambdaTcal = []
                self.miuaTcal = []

            mid_index = self.mesh.n // 2
            self.miuaTcal.append(miua[mid_index, :].copy())
            self.lambdaTcal.append((lambda_a[:, mid_index] / (1e-12 + (1 - lambda_a[:, mid_index]))).copy())

    def _compmob(self, s, miua, c, sor, swr, aqueous=True):
        """
        Computes phase mobilities (aqueous or oleic).
        """
        miuo = self.miuo

        surf = 1 if self.surfactant is not None and self.surfactant.concentration > 0 else 0

        if surf == 0:
            # No surfactant case
            nsw0 = (s - self.init_water_saturation) / (1 - self.init_water_saturation)
            nso0 = (s - self.init_water_saturation) / (1 - self.init_water_saturation - self.init_oleic_saturation)

            krw0 = nsw0 ** 3.5
            kro0 = ((1 - nso0) ** 2) * (1 - nso0 ** 1.5)

            if aqueous:
                lambda_val = krw0 / miua
            else:
                lambda_val = kro0 / miuo

        else:
            # With surfactant case
            nsw = (s - swr) / (1 - swr)
            nso = (s - swr) / (1 - swr - sor)

            krw = nsw * (2.5 * swr * (nsw ** 2 - 1) + 1)
            kro = (1 - nso) * (1 - 5 * sor * nso)

            if aqueous:
                lambda_val = krw / miua
            else:
                lambda_val = kro / miuo

        return lambda_val

    def _initialize_fields(self):
        """
        Initialize water saturation (s0), polymer concentration (c0), and surfactant concentration (g0).
        Equivalent to MATLAB s0c0 function.
        """
        m = self.mesh.m
        n = self.mesh.n
        phi = self.phi

        s0 = np.zeros((n+1, m+1))
        c0 = np.zeros_like(s0)
        g0 = np.zeros_like(s0)

        D = (phi > 1e-10) | (np.abs(phi) < 1e-10)

        # NOTE: ~D means "outside domain", D means "inside domain"
        logical_not_D = np.logical_not(D)
        s0 = logical_not_D.astype(np.float64) + D.astype(np.float64) * (1 - self.init_water_saturation)
        c0 = logical_not_D.astype(np.float64) * self.c0
        g0 = logical_not_D.astype(np.float64) * self.g0

        return s0, c0, g0
    
    def _solve_flow(self, beta):
        """
        Solves the elliptic system for pressure and computes velocity fields (u, v).
        """
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve

        # Set up dimensions
        n, m = self.mesh.n, self.mesh.m
        dx, dy = self.mesh.dx, self.mesh.dy

        # Build system matrices (discrete Laplacian)
        num_nodes = (n + 1) * (m + 1)

        rows = []
        cols = []
        data = []
        b = np.zeros(num_nodes)

        def idx(i, j):
            return i * (m + 1) + j

        for i in range(n + 1):
            for j in range(m + 1):
                center = idx(i, j)

                if i > 0:
                    up = idx(i - 1, j)
                    rows.append(center)
                    cols.append(up)
                    data.append(-beta[i, j] / dy ** 2)

                if i < n:
                    down = idx(i + 1, j)
                    rows.append(center)
                    cols.append(down)
                    data.append(-beta[i, j] / dy ** 2)

                if j > 0:
                    left = idx(i, j - 1)
                    rows.append(center)
                    cols.append(left)
                    data.append(-beta[i, j] / dx ** 2)

                if j < m:
                    right = idx(i, j + 1)
                    rows.append(center)
                    cols.append(right)
                    data.append(-beta[i, j] / dx ** 2)

                # Center coefficient
                rows.append(center)
                cols.append(center)
                data.append(2 * beta[i, j] * (1/dx**2 + 1/dy**2))

        # Assemble sparse matrix A
        A = coo_matrix((data, (rows, cols)), shape=(num_nodes, num_nodes)).tocsc()

        # Right-hand side is zeros (steady-state problem with fixed flux at wells)
        B = b

        # Solve the linear system
        p_flat = spsolve(A, B)

        # Reshape back to 2D pressure field
        p = p_flat.reshape((n + 1, m + 1))

        # Compute velocities u, v = -beta * gradient(p)
        u = np.zeros_like(p)
        v = np.zeros_like(p)

        # Compute gradients
        u[1:-1, :] = -(p[2:, :] - p[:-2, :]) / (2 * dy)
        v[:, 1:-1] = -(p[:, 2:] - p[:, :-2]) / (2 * dx)

        # Multiply by beta
        u = beta * u
        v = beta * v

        return u, v
    
    def _characteristic_coordinates(self, x, y, s, snew, g, f, f_s, D, pc_s, pc_g, u, v, dt, para, flag):
        """
        Compute redefined characteristic coordinates (xmod, ymod) according to Neumann boundary conditions.
        """
        dx, dy = para.box.dx, para.box.dy
        m, n = para.box.m, para.box.n

        if flag == 1:
            xjump = x - f_s * u * dt
            yjump = y - f_s * v * dt
        elif flag == 2:
            # Calculate gradients
            sx, sy = np.gradient(s, dx, dy, edge_order=2)
            gx, gy = np.gradient(g, dx, dy, edge_order=2)

            xjump = x - ((f / snew) * u + (D * pc_s / snew) * sx + (D * pc_g / snew) * gx) * dt
            yjump = y - ((f / snew) * v + (D * pc_s / snew) * sy + (D * pc_g / snew) * gy) * dt
        elif flag == 3:
            sx, sy = np.gradient(s, dx, dy, edge_order=2)

            xjump = x - ((f / snew) * u + (D * pc_s / snew) * sx) * dt
            yjump = y - ((f / snew) * v + (D * pc_s / snew) * sy) * dt

        # Apply Neumann reflection conditions
        xmod = np.where(xjump <= 1, np.abs(xjump), 2 - xjump)
        ymod = np.where(yjump <= 1, np.abs(yjump), 2 - yjump)

        return xmod, ymod

    def _solve_transport(self, u, v, sigma):
        from scipy.sparse import csc_matrix
        from scipy.interpolate import RegularGridInterpolator

        dt = self.time_step
        dx, dy = self.mesh.dx, self.mesh.dy
        n, m = self.mesh.n, self.mesh.m

        x, y = np.meshgrid(
            np.arange(self.mesh.left, self.mesh.right + dx, dx),
            np.arange(self.mesh.bottom, self.mesh.top + dy, dy)
        )

        S = self.water_saturation.copy()
        C = self.polymer.vec_concentration.copy()
        G = self.surfactant.vec_concentration.copy()

        phi = 1
        omega1 = SimulationConstants.Capillary_Pressure_Param_1.value
        omega2 = SimulationConstants.Capillary_Pressure_Param_2.value

        [swr, sor] = self.compres(sigma, u, v, self.aqueous_viscosity)
        nsw = (S - swr) / (1 - swr)
        nso = (S - swr) / (1 - swr - sor)

        lambda_a = self.compmob(S, self.aqueous_viscosity, C, sor, swr, 1, 1)
        lambda_o = self.compmob(S, self.aqueous_viscosity, C, sor, swr, 0, 1)
        lambda_total = lambda_a + lambda_o

        f = lambda_a / lambda_total
        D = self.KK * lambda_o * f

        sigma_g = -10.001 / (G + 1) ** 2

        pc = (sigma * omega2 * phi ** 0.5) / (self.KK ** 0.5 * (1 - nso) ** (1/omega1))
        pc_s = pc / (omega1 * (1 - nso))
        pc_g = (pc / sigma) * sigma_g + pc_s

        xmod, ymod = self._characteristic_coordinates(x, y, S, S, G, f, f, D, pc_s, pc_g, u, v, dt, self, flag=1)

        interp_S = RegularGridInterpolator((y[:, 0], x[0, :]), S)
        coords = np.array([ymod.flatten(), xmod.flatten()]).T
        Snew = interp_S(coords).reshape(n+1, m+1)
        Snew = np.clip(Snew, 0, 1)

        self.water_saturation = Snew

        # Now solve polymer concentration
        xmod_c, ymod_c = self._characteristic_coordinates(x, y, S, Snew, G, f, f, D, pc_s, pc_g, u, v, dt, self, flag=2)
        interp_C = RegularGridInterpolator((y[:, 0], x[0, :]), C)
        coords_c = np.array([ymod_c.flatten(), xmod_c.flatten()]).T
        Cmod = interp_C(coords_c).reshape(n+1, m+1)
        Cnew = np.clip(Cmod, 0, self.polymer.initial_concentration)

        self.polymer.vec_concentration = Cnew

        # Now solve surfactant concentration
        xmod_g, ymod_g = self._characteristic_coordinates(x, y, S, Snew, G, f, f, D, pc_s, pc_g, u, v, dt, self, flag=3)
        interp_G = RegularGridInterpolator((y[:, 0], x[0, :]), G)
        coords_g = np.array([ymod_g.flatten(), xmod_g.flatten()]).T
        Gmod = interp_G(coords_g).reshape(n+1, m+1)
        Gnew = np.clip(Gmod, 0, self.surfactant.concentration)

        self.surfactant.vec_concentration = Gnew

        # Return oil/water production and ROIP
        lambda_a = self.compmob(Snew, self.aqueous_viscosity, Cnew, sor, swr, 1, 1)
        lambda_o = self.compmob(Snew, self.aqueous_viscosity, Cnew, sor, swr, 0, 1)
        lambda_total = lambda_a + lambda_o

        prod_oil_vol = lambda_o[-1, -1] * self.src / lambda_total[-1, -1]
        prod_water_vol = lambda_a[-1, -1] * self.src / lambda_total[-1, -1]
        ROIP = 100 * np.sum(1 - Snew) / (n * m)

        return {
            "S": Snew,
            "C": Cnew,
            "G": Gnew,
            "prod_oil_vol": prod_oil_vol,
            "prod_water_vol": prod_water_vol,
            "ROIP": ROIP
        }


    def _export_results(self):
        """
        Saves simulation sim_results to CSV files.
        """
        import os
        os.makedirs("sim_results", exist_ok=True)

        np.savetxt("sim_results/COC.csv", self.COC, delimiter=",")
        if hasattr(self, 'lambdaTcal'):
            np.savetxt("sim_results/lambdaTcal.csv", np.array(self.lambdaTcal), delimiter=",")
        if hasattr(self, 'miuaTcal'):
            np.savetxt("sim_results/miuaTcal.csv", np.array(self.miuaTcal), delimiter=",")
        
        print("Simulation sim_results exported to /sim_results/ folder.")


    def run(self):
        """
        Executes simulation loop.
        """
        self.initialize_simulation()

        # --- Initialize simulation tracking variables ---
        self.tstop = 500  # seconds
        self.src_total = 0
        self.sumUU = 0

        self.COC = np.zeros((1, 2000), dtype=np.float64)
        
        num_time_steps = int(np.floor(self.tstop / self.time_step))
        self.ProdRate = np.memmap(
            "memmaps/ProdRate.dat",
            dtype="float64",
            mode="w+",
            shape=(1, num_time_steps)
        )
        self.CROIP = np.memmap(
            "memmaps/CROIP.dat",
            dtype="float64",
            mode="w+",
            shape=(1, num_time_steps)
        )

        self.miuaSave = 0
        self.shearSave = 0
        self.concSave = 0

        # Initialize velocity fields
        u = np.zeros((self.mesh.n + 1, self.mesh.m + 1), dtype=np.float64)
        v = np.zeros_like(u)

        t = 0
        tcal = 0

        print(f"Running simulation {self.sim_id} with grid size {self.grid_size}x{self.grid_size}")
        print(f"Model Type: {self.model_type.name}")
        print(f"Initial polymer concentration: {self.c0}")
        print(f"Initial surfactant concentration: {self.g0}")
        print(f"Calculated dt: {self.time_step:.5e}")

        breakthrough = False
        # --- Main time-stepping loop ---
        while t < self.tstop and self.U[-1, -1] <= 0.70:
            print("I am inside nowww")
            self.src_total += self.source_flow_magnitude
            t += self.time_step

            inner_iter = 0
            epsilon = 10

            # Step 1: Surface tension
            sigma = 10.001 / (self.G + 1) - 0.001

            # Step 2: Aqueous viscosity and shear
            miua, shear = self._compvis(self.C, u, v)

            # Step 3: Apply no shear-thinning if necessary
            if self.model_type == ModelType.No_Shear_Thinning:
                self.miup = np.max(miua[0, :])
                self.miup_array = self.miup * np.ones_like(self.U)

            # Step 4: Residual saturations
            swr, sor = self._compres(sigma, u, v, miua)

            # Step 5: Mobilities
            lambda_a = self._compmob(self.U, miua, self.C, sor, swr, aqueous=True)
            lambda_o = self._compmob(self.U, miua, self.C, sor, swr, aqueous=False)
            lambda_total = lambda_a + lambda_o

            # Step 6: Beta
            beta = self.KK * lambda_total

            # Step 7: Solve flow field
            u, v = self._solve_flow(beta)

            # Step 8: Shear effects if needed
            if self.shear_flag:
                miua = self._apply_shear_effects(u, v, miua)

            # Step 9: Solve transport equations
            # self.U, self.C, self.G, OC, WC, ROIP = self._solve_transport(u, v, sigma)
            transport_results = self._solve_transport(u, v, sigma)
            self.U = self.water_saturation = transport_results["S"]
            self.C = self.polymer.vec_concentration = transport_results["C"]
            self.G = self.surfactant.vec_concentration = transport_results["G"]
            OC = transport_results["prod_oil_vol"]
            WC = transport_results["prod_water_vol"]
            ROIP = transport_results["ROIP"]            

            # Step 10: Update oil recovery tracking
            if tcal == 0:
                self.COC[0, tcal] = OC
            else:
                self.COC[0, tcal] = self.COC[0, tcal - 1] + OC

            self.ProdRate[0, tcal] = OC / self.time_step
            self.CROIP[0, tcal] = ROIP

            tcal += 1

            # Step 11: Save intermediate fields
            if tcal % 200 == 0:
                self._save_lambda_miua(tcal, lambda_a, miua)
                
        # --- Flush memmaps to disk ---
        self.ProdRate.flush()
        self.CROIP.flush()

        # After time-stepping
        self._export_results()

        print(f"Simulation {self.sim_id} finished. Results saved to 'Results/' and 'memmaps/' folders.")
