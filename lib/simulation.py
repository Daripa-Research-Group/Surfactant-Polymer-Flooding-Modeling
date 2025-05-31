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
from water import Water #TODO: Merge into water refactor branch
from scipy.io import loadmat
os.makedirs("memmaps", exist_ok=True) #ensures that the program works on computers with RAM constraints


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

        :param user_input_dict: dictionary containing the information from the GUI
        :type user_input_dict: dict
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
        
        #initializing properties that hold simulation constants
        self.grid_size = SimulationConstants.Grid_Size.value
        self.source_flow_magnitude = SimulationConstants.Source_Flow_Magnitude.value

        # Initializing sim properties
        self.mesh = self._create_mesh() # TODO: This needs to be replaced with 'Grid' object
        self.x, self.y = self._generate_grid()
        self.phi = None  # Level set function (relates to porosity)
        self.KK = None  # Permeability tensor
        self.time_step = None
        self.initialize_simulation() #will initialize phi, KK, and time_step
        if(self.phi is None or self.KK is None or self.time_step is None): #Raise Exception if not properly initialized...
            raise SimulationCalcInputException("SimulationInputException: phi, KK, or time_step not properly initialized. Please try again...") 
        
        # TODO: Initializing global pressure ('u') and velocity matrices ('v'):
            # Need to be initialized using the Grid class
        self.u = None
        self.v = None

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

        # TODO:Initializing Water Object

        # TODO:Properties for Exporting Simulation Results
        self.COC = np.zeros((1,2000))
        self.miuaTcal = np.zeros((1,2000))
        self.lambdaTcal = np.zeros((1,2000))

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
                raise SimulationCalcInputException('SimulationInputException: KK matrix not found in KK30Ness.mat file.')
            KK = mat_data['KK']
        elif flag == 6:
            # Load Tarbert formation (SPE10)
            mat_data = loadmat('/Resources/KK30Tabert.mat')
            if 'KK' not in mat_data:
                raise SimulationCalcInputException('SimulationInputException: KK matrix not found in KK30Tabert.mat file.')
            KK = mat_data['KK']
        else:
            raise SimulationCalcInputException("SimulationInputException: Unknown permeability flag.")
        return KK
    
    def _characteristic_coordinates(self, x, y, s, snew, g, f, f_s, D, pc_s, pc_g, u, v, dt, para, flag):
        """
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

            xjump = x - ((f / snew) * u + (D * pc_s / snew) * sx + (D * pc_g / snew) * gx) * dt
            yjump = y - ((f / snew) * v + (D * pc_s / snew) * sy + (D * pc_g / snew) * gy) * dt
        elif flag == 3:
            sx, sy = np.gradient(s, dx, dy, edge_order=2)

            xjump = x - ((f / snew) * u + (D * pc_s / snew) * sx) * dt
            yjump = y - ((f / snew) * v + (D * pc_s / snew) * sy) * dt

        # Apply Neumann reflection conditions
        if(xjump is None or yjump is None):
            raise SimulationCalcInputException("SimulationInputException: xjump or yjump not initialized...")
        
        xmod = np.where(xjump <= 1, np.abs(xjump), 2 - xjump)
        ymod = np.where(yjump <= 1, np.abs(yjump), 2 - yjump)

        return xmod, ymod

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

        :return: Dictionary with relevant results for plotting and data analysis
        :rtype: dict
        """
        try:
            # TODO: Updating the Grid with the initial source and production well flowrates (actual steps will be done within the 'Grid' Class)


            # TODO: Running primary while loop to iterate through time-steps
            t = 0
            tf = 500
            dt = self.mesh.dx / self.source_flow_magnitude
            tcal = 0
            tsave = 0
            viscosity_aqueous_save = 0
            shear_force_save = 0
            concentration_save = 0
            source_flow_magnitude_total = 0
            sum_of_saturation_matrix = 0

            # This while loop will run as long as little to no water in production well and current timestep < final timestep
            while(t < tf and self.water.water_saturation[self.mesh.n, self.mesh.m] <= 0.70):
                #update the total source flow magnitude:
                source_flow_magnitude_total += self.source_flow_magnitude

                #update time:
                t += dt
                innerIter = 0
                
                #compute viscosity:
                if(self.model_type == ModelType.No_Shear_Thinning):
                    # Will need the aqueous viscosity for no shear thinning
                    if(self.x is None or self.y is None or self.u is None or self.v is None):
                        raise SimulationCalcInputException("SimulationInputException: Grid, global pressure, or velocity not initialized. Please try again.")
                    grid =(self.x, self.y)
                    [polymer_viscosity_matrix, shear_rate_matrix] = \
                            self.polymer.compute_viscosity(grid=grid, u=self.u, v=self.v, model_type=self.model_type, aqueous_viscosity=self.water.viscosity_array)
                elif(self.model_type == ModelType.Shear_Thinning_On):
                    if(self.x is None or self.y is None or self.u is None or self.v is None):
                        raise SimulationCalcInputException("SimulationInputException: Grid, global pressure, or velocity not initialized. Please try again.")
                    grid =(self.x, self.y)
                    [polymer_viscosity_matrix, shear_rate_matrix] = \
                            self.polymer.compute_viscosity(grid=grid, u=self.u, v=self.v, model_type=self.model_type)
                else:
                    raise SimulationCalcInputException("SimulationInputException: Improper model type provided... please try again...")

            #TODO:Update 'tsave':
            tsave += tcal

            #TODO:calculating mobilities of wetting and non-wettting phases
                #need to invoke function within 'Water' class

            #TODO:Update beta value

            #TODO:Update grid & update global pressure (u) and velocity(v):

            #TODO: Solve transport equations to update concentration & saturation matrices:



        except Exception as e:
            print(e)

