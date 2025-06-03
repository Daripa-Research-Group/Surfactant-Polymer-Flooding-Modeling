"""
This python script contains the class definition for the water in the surfactant-flooding model

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from enumerations import SimulationConstants


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
        n, m = grid_shape
        s0 = np.zeros((n + 1, m + 1))
        D = (self.phi > 1e-10) | (np.abs(self.phi) < 1e-10)
        s0 = np.logical_not(D).astype(float) + D.astype(float) * (1 - self.init_water_saturation)
        self.water_saturation = s0
        return s0

    def compute_viscosity(
            self, 
            c: np.ndarray, 
            u: np.ndarray, 
            v: np.ndarray, 
            c0: float
            ):
        """
        Compute aqueous viscosity (NO shear-thinning version).
        """
        beta1 = 15000
        if c0 == 0:
            miua = self.miuw * np.ones_like(c)
        else:
            miua = self.miuw * (1 + beta1 * c)
        self.viscosity_array = miua
        return miua, np.zeros_like(c) 

    def compute_residual_saturations(
            self, 
            sigma: np.ndarray, 
            u: np.ndarray, 
            v: np.ndarray
            ):
        """
        Compute swr, sor based on capillary numbers.
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

        return swr, sor

    def compute_mobility(
            self, 
            s: np.ndarray, 
            c: np.ndarray, 
            miua: np.ndarray, 
            sor, 
            swr, 
            aqueous: bool, 
            has_surfactant: bool, 
            surfactant_conc: float
            ):
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

    def solve_saturation_transport(
            self,
            x, 
            y, 
            dt, 
            KK, 
            lambda_a, 
            lambda_o, 
            sigma, 
            G, 
            para, 
            u, 
            v
            ):
        dx, dy = para.grid.dx, para.grid.dy
        phi = 1  # Constant
        omega1 = SimulationConstants.Capillary_Pressure_Param_1.value
        omega2 = SimulationConstants.Capillary_Pressure_Param_2.value

        S = self.water_saturation
        nsw = (S - self.init_water_saturation) / (1 - self.init_water_saturation)
        nso = (S - self.init_water_saturation) / (1 - self.init_water_saturation - self.init_oleic_saturation)

        lambda_total = lambda_a + lambda_o
        f = lambda_a / lambda_total
        D = KK * lambda_o * f
        sigma_g = -10.001 / (G + 1) ** 2

        pc = (sigma * omega2 * phi ** 0.5) / (KK ** 0.5 * (1 - nso) ** (1 / omega1))
        pc_s = pc / (omega1 * (1 - nso))
        pc_g = (pc / sigma) * sigma_g + pc_s

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
