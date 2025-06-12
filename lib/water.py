"""
This python script contains the class definition for the water in the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from enumerations import ModelType, SimulationConstants
from Exceptions import SimulationCalcInputException
from grid import Grid
from polymer import Polymer

class Water:
    def __init__(
            self, 
            init_water_saturation: float, 
            init_oleic_saturation: float, 
            miuw: float, 
            miuo: float,
            phi: np.ndarray
            ):
        """
        Constructor for the 'Water' class

        :param init_water_saturation: initial water saturation
        :type init_water_saturation: float

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
        self.init_oleic_saturation = init_oleic_saturation
        self.miuw = miuw
        self.miuo = miuo
        self.water_saturation = None
        self.viscosity_array = None
        self.phi = phi

    def initialize(
            self, 
            grid_shape: tuple
            ):
        """
        Initializing 'Water' Object properties
        
        :param grid_shape: the n and m parameters from the 'Grid' class
        :type grid_shape: tuple

        :return: Updated 'Water' object
        :rtype: Water
        """
        #getting values for n and m from grid:
        n, m = grid_shape

        #initializing the water saturation matrix
        s0 = np.zeros((n + 1, m + 1))
        D = (self.phi > 1e-10) | (np.abs(self.phi) < 1e-10)
        s0 = np.logical_not(D).astype(float) + D.astype(float) * (1 - self.init_water_saturation)
        self.water_saturation = s0

        #initialize the aqueous viscosity matrix:
        self.viscosity_array = self.miuw*np.ones((n,m))

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
        Compute aqueous viscosity (NO shear-thinning version).

        :param grid: Grid object for deterrmining matrix size
        :type grid: Grid

        :param polymer: holds the information about the polymer in the sim
        :type polymer: Polymer

        :param u: global pressure matrix. Only needed when shear thinning ON.
        :type u: np.ndarray, None

        :param v: velocity matrix. Only needed when shear thinning ON.
        :type v: np.ndarray, None

        :return: updated aqueous viscosity matrix
        :rtype: np.ndarray
        """
        assert self.viscosity_array is not None, SimulationCalcInputException("SimuationInputException: aqueous viscosity matrix not initialized. Please try again")
        assert polymer is not None, SimulationCalcInputException("SimuationInputException: polymer object not initialized. \
                The polymer object must be initialized before updating aqueous viscosity. Please try again.")
        n = np.size(polymer.concentration_matrix, 0)
        m = np.size(polymer.concentration_matrix, 0)
        initial_polymer_concentration_scalar = polymer.concetration_scalar
        if(model_type.value == ModelType.No_Shear_Thinning.value): #no shear thinning polymer
            miuw = SimulationConstants.Water_Viscosity.value
            if(initial_polymer_concentration_scalar == 0):
                self.viscosity_array = miuw*np.ones((n,m))
            else:
                beta1 = SimulationConstants.beta1.value
                self.viscosity_array = miuw*(1+beta1*polymer.concentration_matrix)
        elif(model_type.value == ModelType.Shear_Thinning_On.value): #shear thinning polymer
            # using the shear rate and polymer coefficients to understand how its viscosity changes
            assert u is not None, SimulationCalcInputException("SimuationInputException: variables 'u' not initialized for shear-thinning-on model. Please try again")
            assert v is not None, SimulationCalcInputException("SimuationInputException: variables 'v' not initialized for shear-thinning-on model. Please try again")
            
            #constants:
            rho_water = SimulationConstants.Water_Density.value
            viscosity_water = SimulationConstants.Water_Viscosity.value

            #relevant parameters for power law equation:
            w1 = polymer.rho*polymer.concentration_matrix
            w2 = rho_water*(1-polymer.concentration_matrix)
            wppm = (w1/(w1+w2))*(10**6)
            w1_0 = polymer.rho*polymer.init_concentration_matrix #from the variable w10 in MATLAB code
            w2_0 = rho_water*(1-polymer.init_concentration_matrix) #from the variable w20 in MATLAB code
            wppm_0 = (w1_0/(w1_0+w2_0))*(10**6) #from the wppm0 variable in MATLAB code
            
            #determining ε and n for power law equation:
            epsilon_0 = polymer.e_coeff[0]*(wppm_0**polymer.e_coeff[1])
            n_0 = np.min(polymer.n_coeff[0]*(wppm_0**polymer.n_coeff[1]))
            epsilon_val = polymer.e_coeff[0]*(wppm**polymer.e_coeff[1])
            n_val = np.min(polymer.n_coeff[0]*(wppm**polymer.n_coeff[1]))

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
                    self.viscosity_array[ii,jj] = epsilon_val[ii,jj]*(polymer.shear_rate[ii,jj]**(n_val[ii,jj]-1))
                    if self.viscosity_array[ii, jj] < viscosity_water:
                        self.viscosity_array[ii, jj] = viscosity_water
                    if self.viscosity_array[ii, jj] > 100:
                        self.viscosity_array[ii, jj] = 100
        return self.viscosity_array 
        

    def compute_residual_saturations(
            self, 
            sigma: np.ndarray, 
            u: np.ndarray, 
            v: np.ndarray
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
        swr0 = self.init_water_saturation
        sor0 = self.init_oleic_saturation

        Nco0 = 1.44e-4
        Nca0 = 1.44e-4

        vel_mag = np.sqrt(u ** 2 + v ** 2)
        nca = (vel_mag * self.viscosity_array) / sigma
        nco = (vel_mag * self.miuo) / sigma

        Nca = np.linalg.norm(nca)
        Nco = np.linalg.norm(nco)

        sor = sor0 * (Nco0 / Nco) ** 0.5213 if Nco >= Nco0 else sor0
        swr = swr0 * (Nca0 / Nca) ** 0.1534 if Nca >= Nca0 else swr0

        return [swr, sor] #[residual water saturation, residual oil saturation]

    def compute_mobility(
            self, 
            c: np.ndarray, 
            sor: np.ndarray, 
            swr: np.ndarray, 
            aqueous: bool, 
            has_surfactant: bool, 
            surfactant_conc: float
            ):
        """
        Computing mobility (made using the compmob.m MATLAB file)
        
        :param c: polymer concentration matrix
        :type c: np.ndarray

        :param sor: residual saturation oil phase
        :type sor: np.ndarray

        :param swr: residual saturation water phase
        :type swr: np.ndarray

        :param aqueous: boolean for whether we are solving for aqoeous or oleic mobility
        :type aqueous: bool

        :param has_surfactant: whether or not there is surfactant in the system
        :type has_surfactant: bool

        :param surfactant_conc: scalar quantity of the initial surfactant concentration
        :type surfactant_conc: float

        :return: aqueous or oleic mobility (depending on the 'aqueous' parameter)
        :rtype: np.ndarray
        """
        assert self.water_saturation is not None, SimulationCalcInputException("SimuationInputException: water saturation matrix not initialized. Please try again")
        assert self.viscosity_array is not None, SimulationCalcInputException("SimuationInputException: viscosity matrix not initialized. Please try again")
        s = self.water_saturation
        miua = self.viscosity_array
        if not has_surfactant or surfactant_conc == 0:
            nsw0 = (s - self.init_water_saturation) / (1 - self.init_water_saturation)
            nso0 = (s - self.init_water_saturation) / (1 - self.init_water_saturation - self.init_oleic_saturation)
            krw0 = nsw0 ** 3.5
            kro0 = ((1 - nso0) ** 2) * (1 - nso0 ** 1.5)
        else:
            nsw = (s - swr) / (1 - swr)
            nso = (s - swr) / (1 - swr - sor)
            krw0 = nsw * (2.5 * swr * (nsw ** 2 - 1) + 1)
            kro0 = (1 - nso) * (1 - 5 * sor * nso)

        return krw0 / miua if aqueous else kro0 / self.miuo

    def compute_water_saturation(
            self,
            grid: Grid,
            dt: float, 
            KK: np.ndarray, 
            lambda_a: np.ndarray, 
            lambda_o: np.ndarray, 
            sigma: np.ndarray, 
            surfactant: Surfactant, 
            u: np.ndarray,
            v: np.ndarray,
            ):
        """
        Solving saturation equation (comes from part of the nmmoc_surf_mod_neumann.m file that 
        is for calculating the water saturation)

        :param grid: the 'Grid' object
        :type grid: Grid

        :param dt: time step
        :type dt: float

        :param KK: the permeability tensor
        :type KK: np.ndarray

        :param lambda_a: aqueous mobility
        :type lambda_a: np.ndarray

        :param lambda_o: oleic mobility
        :type lambda_o: np.ndarray

        :param sigma: IFT (interfacial tension)
        :type sigma: np.ndarray

        :param G: surfactant concentration (matrix form)
        :type G: np.ndarray

        :param u: global pressure matrix
        :type u: np.ndarray

        :param v: velocity matrix
        :type v: np.ndarray

        :return: updated water saturation matrix
        :rtype: np.ndarray
        """
        assert self.water_saturation is not None, SimulationCalcInputException("SimuationInputException: water saturation matrix not initialized. Please try again")
        dx, dy = grid.dx, grid.dy
        phi = 1  # Constant
        omega1 = SimulationConstants.Capillary_Pressure_Param_1.value
        omega2 = SimulationConstants.Capillary_Pressure_Param_2.value

        S = self.water_saturation
        nsw = (S - self.init_water_saturation) / (1 - self.init_water_saturation)
        nso = (S - self.init_water_saturation) / (1 - self.init_water_saturation - self.init_oleic_saturation)

        lambda_total = lambda_a + lambda_o
        f = lambda_a / lambda_total
        D = KK * lambda_o * f
        sigma_g = surfactant.derivative_IFT_equation(surfactant.concentration_matrix) #-10.001 / (G + 1) ** 2

        pc = (sigma * omega2 * phi ** 0.5) / (KK ** 0.5 * (1 - nso) ** (1 / omega1))
        pc_s = pc / (omega1 * (1 - nso))
        pc_g = (pc / sigma) * sigma_g + pc_s
        
        x = grid.x
        y = grid.y
        xjump = x - f * u * dt
        yjump = y - f * v * dt
        xmod = np.where(xjump <= 1, np.abs(xjump), 2 - xjump)
        ymod = np.where(yjump <= 1, np.abs(yjump), 2 - yjump)

        interp = RegularGridInterpolator((y[:, 0], x[0, :]), S)
        coords = np.array([ymod.flatten(), xmod.flatten()]).T
        Snew = interp(coords).reshape(S.shape)
        Snew = np.clip(Snew, 0, 1)

        self.water_saturation = Snew

        return Snew
