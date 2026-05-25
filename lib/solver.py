"""
Will hold the classes relevant for solving the elliptic equations and the transport equations

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

from surfactant import Surfactant
from water import Water
from polymer import Polymer
from grid import Grid
from enumerations import (
    ModelType,
    PolymerList,
    RelativePermeabilityFormula,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from Exceptions import SimulationCalcInputException
import numpy as np
import scipy as sp
from scipy.linalg import fractional_matrix_power
from scipy.sparse.linalg import bicgstab


class TransportEquationSolver():
    """
    Will hold all the relevant functions and parameters necessary for solving the transport equations
    """
    def __init__(
            self,
            grid : Grid,
            water : Water,
            surfactant : Surfactant,
            polymer : Polymer,
            pressure : np.ndarray,
            velocity : np.ndarray,
            permeability_matrix : np.ndarray,
            source_flow_magnitude: float
    ):
        """
        constructor for ``TransportEquationSolver``
        """
        self._water = water
        self._grid = grid
        self._surfactant = surfactant
        self._polymer = polymer

        self.pressure = pressure
        self.velocity = velocity
        self.permeability_matrix = permeability_matrix
        self.total_flow = source_flow_magnitude

    #Constant Parameter Definitions
    ## Flows
    _total_flow = None
    @property
    def total_flow(self):
        return self._total_flow
    @total_flow.setter
    def total_flow(self, value):
        self._total_flow = value
    @property
    def polymer_flow(self):
        return self.total_flow*self._water.concetration_scalar
    @property
    def surfactant_flow(self):
        return self.total_flow*self._surfactant.concentration
    ## Oil Viscosity
    @property
    def oil_viscosity(self):
        return self._water.miuo
    ## Elliptic pressure
    _pressure = None
    @property
    def pressure(self):
        return self._pressure
    @pressure.setter
    def pressure(self, value):
        self._pressure = value
    ## velocity
    _velocity = None
    @property
    def velocity(self):
        return self._velocity
    @velocity.setter
    def velocity(self, value):
        self._velocity = value
    ## Permeability Matrix
    _permeability_matrix = None
    @property
    def permeability_matrix(self):
        return self._permeability_matrix
    @permeability_matrix.setter
    def permeability_matrix(self, value):
        self._permeability_matrix = value
    ## Grid properties
    @property
    def m(self):
        return self._grid.m
    @property
    def n(self):
        return self._grid.n
    @property
    def dx(self):
        return self._grid.dx
    @property
    def dy(self):
        return self._grid.dy
    @property
    def x(self):
        return self._grid.x
    @property
    def y(self):
        return self._grid.y
    @property
    def dt_array(self):
        pass
    ## Porosity
    @property
    def porosity(self):
        return 1
    ## Capillary pressure parameters
    @property
    def omega(self):
        return [0.1, 0.4] #[Ω_1, Ω_2]
    ## Capillary Number
    @property
    def critical_capillary_aqueous_initial(self): #Nca0
        return SimulationConstants.Aqueous_Phase_Critical_Capillary_Num.value 
    @property
    def critical_capillary_oleic_initial(self): #Nco0
        return SimulationConstants.Oleic_Phase_Critical_Capillary_Num.value
    
    # Functions for parameter definitions (TODO: need to make sure that when current property value is None, respective calculations are automatically done when property is called)
    ## Surfactant concentration matrix
    @property
    def surfactant_concentration(self):
        return self._surfactant.concentration_matrix
    @surfactant_concentration.setter
    def surfactant_concentration(self, value):
        self._surfactant.concentration_matrix = value
    ## Polymer concentration matrix
    @property
    def polymer_concentration(self):
        return self._polymer.concentration_matrix
    @polymer_concentration.setter
    def polymer_concentration(self, value):
        self._polymer.concentration_matrix = value
    ## Water Saturation Matrix
    @property
    def water_saturation(self):
        return self._water.water_saturation
    @water_saturation.setter
    def water_saturation(self, value):
        self._water.water_saturation = value
    ## Aqueous Viscosity
    @property
    def aqueous_viscosity(self):
        return self._water.viscosity_array
    ## Residual saturations
    @property
    def swr_0(self):
        return SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
    _swr = None
    @property
    def swr(self):
        """
        swr: aqueous residual saturation

        Water that remains immobile after it has been displaced by a non-wetting phase (like oil or gas).
        """
        return self._swr
    @swr.setter
    def swr(self, value):
        self._swr = value
    _dswr_dg = None
    @property
    def dswr_dg(self):
        """
        derivative of the residual water saturation wrt surfactant concentration
        """
        return self._dswr_dg
    @dswr_dg.setter
    def dswr_dg(self, value):
        self._dswr_dg = value
    _nsw = None
    @property
    def nsw(self):
        """
        nsw: effective aqueous saturation
        """
        return self._nsw
    @nsw.setter
    def nsw(self, value):
        self._nsw = value
    _dnsw_dg = None
    @property
    def dnsw_dg(self):
        return self._dnsw_dg
    @dnsw_dg.setter
    def dnsw_dg(self, value):
        self._dnsw_dg = value
    
    @property
    def sor_0(self):
        return SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value
    _sor = None
    @property
    def sor(self):
        """
        sor: oleic residual saturation

        Oil that remains immobile after it has been displaced by a water-based phase 
        (this could be waterflooding or SP flooding).

        As the simulation is running, the residual oil saturation decreases as the interfacial tension (σ)
        decreases between wetting (water) and non-wetting (oil) phases.
        """
        return self._sor
    @sor.setter
    def sor(self, value):
        self._sor = value
    _dsor_dg = None
    @property
    def dsor_dg(self):
        """
        derivative of the residual oil saturation wrt surfactant concentration
        """
        return self._dsor_dg
    @dsor_dg.setter
    def dsor_dg(self, value):
        self._dsor_dg = value
    _nso = None
    @property
    def nso(self):
        """
        nso: effective oleic saturation
        """
        return self._nso
    @nso.setter
    def nso(self, value):
        self._nso = value
    _dnso_dg = None
    @property
    def dnso_dg(self):
        return self._dnso_dg
    @dnso_dg.setter
    def dnso_dg(self, value):
        self._dnso_dg = value

    def _compute_residual_saturations(self, sigma: np.ndarray):
        """
        wrapper function for running the ``self.water.compute_residual_saturation()`` method
        
        Compute swr, sor based on capillary numbers (came from compres.m MATLAB file)

        Args:
        -----
            sigma (np.ndarray): interfacial tension (IFT)

            u (np.ndarray): global pressure matrix.

            v (np.ndarray): velocity matrix.
        """
        [self.swr, self.sor] = self._water.compute_residual_saturation(sigma, self.pressure, self.velocity)

    def _compute_effective_saturations(
            self, w_sat_matrix : np.ndarray 
    ):
        assert (self.swr is not None) and (self.sor is not None),SimulationCalcInputException("SimulationInputException: residual saturations not computed")
 
        # effective water saturation
        self.nsw = (w_sat_matrix - self.swr) / (1 - self.swr)
        # effective oil saturation
        self.nso = (w_sat_matrix - self.swr) / (1 - self.swr - self.sor)
    
    def _derivative_effective_saturations(self, modified_water_saturation: np.ndarray | None = None):
        water_saturation_matrix = self.water_saturation
        if modified_water_saturation is not None:
            water_saturation_matrix = modified_water_saturation

        self.dnsw_dg = self.dswr_dg * (water_saturation_matrix - 1)/(1 - self.swr)**2
        self.dnso_dg = (self.dswr_dg * (water_saturation_matrix + self.sor - 1) + self.dsor_dg * (water_saturation_matrix - self.swr))/(1-self.swr-self.sor)**2

    def _derivative_residual_saturations(
            self, sigma, norm_nca, norm_nco
    ): #FIXME: Will need to update function to work with autodiff (v2.0)
        for j in range(self.n):
            for i in range(self.m):
                if norm_nca >= self.critical_capillary_aqueous_initial:
                    self.dswr_dg[j, i] = -(
                        self.swr_0
                        * 0.1534
                        * 10.001
                        * (
                            self.critical_capillary_aqueous_initial
                            ** 0.1534
                        )
                    ) / (
                        (
                            np.sqrt((self.pressure[j, i] ** 2) + (self.velocity[j, i] ** 2))
                            * self.aqueous_viscosity[j, i]
                        )
                        ** (0.1534)
                        * sigma[j, i] ** (0.8466)
                        * (self.surfactant_concentration[j, i] + 1) ** 2
                    )
                if norm_nco >= self.critical_capillary_oleic_initial:
                    self.dsor_dg[j, i] = -(
                        self.sor_0
                        * 0.5213
                        * 10.001
                        * (
                            self.critical_capillary_aqueous_initial
                            ** 0.5213
                        )
                    ) / (
                        (
                            np.sqrt((self.pressure[j, i] ** 2) + (self.velocity[j, i] ** 2))
                            * self.aqueous_viscosity[j, i]
                        )
                        ** (0.5213)
                        * sigma[j, i] ** (0.4787)
                        * (self.surfactant_concentration[j, i] + 1) ** 2
                    )
    
    ## Mobilities
    _lambda_a = None
    @property
    def lambda_a(self):
        return self._lambda_a
    @lambda_a.setter
    def lambda_a(self, value):
        self._lambda_a = value
    _lambda_o = None
    @property
    def lambda_o(self):
        return self._lambda_o
    @lambda_o.setter
    def lambda_o(self, value):
        self._lambda_o = value
    @property
    def lambda_total(self):
        return self.lambda_a + self.lambda_o
    @property
    def fractional_flow(self):
        return self.lambda_a / self.lambda_total
    @property
    def df_ds(self):
        return self.dkra_ds \
            * self.lambda_o \
            / (self.lambda_total**2 * self.aqueous_viscosity) \
            - self.dkro_ds \
            * self.lambda_a \
            / (self.lambda_total**2 * self.oil_viscosity)
    @property
    def df_dc(self):
        return (-1 * (self.lambda_o * self.lambda_a * self.oil_viscosity)) \
            / ((self.lambda_total**2) * self.aqueous_viscosity)
    @property
    def df_dg(self):
        return ( \
                (\
                    self.dkra_dg \
                    * self.lambda_o \
                )\
                / ((self.lambda_total**2) * self.aqueous_viscosity) \
            )\
            - (\
                (\
                    self.dkro_dg \
                    * self.lambda_a \
                )\
                / ((self.lambda_total**2) * self.aqueous_viscosity)\
            )\
    @property
    def D(self):
        return self.permeability_matrix*self.lambda_o*self.fractional_flow
    @property
    def dD_dg(self):
        return self.D * self.dpc_dg
    @property
    def dD_ds(self):
        return self.D * self.dpc_ds

    def _compute_lambda_a(
        self,
        rel_permeability_formula: RelativePermeabilityFormula,
        modified_water_saturation: np.ndarray | None = None,
    ):
        if(modified_water_saturation is not None):
            self.lambda_a = self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                True, 
                RelativePermeabilityFormula.CoreyTypeEquation.value, 
                modified_water_saturation
                )
        else:
            self.lambda_a = self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                True, 
                RelativePermeabilityFormula.CoreyTypeEquation.value 
                )

    def _compute_lambda_o(
        self,
        rel_permeability_formula: RelativePermeabilityFormula,
        modified_water_saturation: np.ndarray | None = None,
    ):
        if(modified_water_saturation is not None):
            self.lambda_o = self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                False, 
                RelativePermeabilityFormula.CoreyTypeEquation.value, 
                modified_water_saturation
                )
        else:
            self.lambda_o = self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                False, 
                RelativePermeabilityFormula.CoreyTypeEquation.value 
                )

    ## Capillary Number (ratio of viscous forces to surface tension forces)
    _nca = None
    @property
    def nca(self):
        return self._nca
    @nca.setter
    def nca(self, value):
        self._nca = value
    @property
    def norm_nca(self):
        return np.linalg.norm(self.nca)

    _nco = None
    @property
    def nco(self):
        return self._nco
    @nco.setter
    def nco(self, value):
        self._nco = value
    @property
    def norm_nco(self):
        return np.linalg.norm(self.nco)
    
    def _aqueous_capillary_number(self, sigma):
        self.nca = (
            np.sqrt(np.matmul(self.pressure, self.pressure) + np.matmul(self.velocity, self.velocity), dtype=np.complex128)
            * self.aqueous_viscosity.astype(np.complex128)
            / sigma.astype(np.complex128)
        )
        
    def _oleic_capillary_number(self, sigma):
        self.nco = (
            np.sqrt(np.matmul(self.pressure, self.pressure) + np.matmul(self.velocity, self.velocity), dtype=np.complex128)
            * SimulationConstants.Oil_Viscosity.value
            / sigma.astype(np.complex128)
        )

    ## Relative Permeability Derivatives wrt water saturation
    _dkra_ds = None
    @property
    def dkra_ds(self):
        return self._dkra_ds
    @dkra_ds.setter
    def dkra_ds(self, value):
        self._dkra_ds = value
    _dkra_dg = None
    @property
    def dkra_dg(self):
        return self._dkra_dg
    @dkra_dg.setter
    def dkra_dg(self, value):
        self._dkra_dg = value
    _dkro_ds = None
    @property
    def dkro_ds(self):
        return self._dkra_ds
    @dkro_ds.setter
    def dkro_ds(self, value):
        self._dkro_ds = value
    _dkro_dg = None
    @property
    def dkro_dg(self):
        return self._dkro_dg
    @dkro_dg.setter
    def dkro_dg(self, value):
        self._dkro_dg = value

    def _derivative_relative_permeabilities(self):
        self.dkra_dg = 2.5 * self.dswr_dg * (self.nsw**3 - self.nsw) \
            + (self.water_saturation - 1) \
            * (2.5 * self.swr * (3 * self.nsw**2 - 1) + 1) \
            * self.dnsw_dg \
            / (1 - self.swr) ** 2
        self.dkro_dg = 1 \
            - 5 * self.sor * self.nso \
            + (1 - self.nso) * (1 - 5 * self.nso * self.dsor_dg) \
            - (1 + 5 * self.sor - 10 * self.sor * self.nso) \
            * self.dnso_dg
        self.dkra_ds = 2.5 * self.swr * (3 * (self.nsw) ** 2 - 1) + 1
        self.dkro_ds = 10 * self.sor * self.nso - 5 * self.sor - 1

    ## Capillary pressure and its derivatives
    _pc = None
    @property
    def pc(self):
        return self._pc
    @pc.setter
    def pc(self, value):
        self._pc = value
    _dpc_ds = None
    @property
    def dpc_ds(self):
        return self._dpc_ds
    @dpc_ds.setter
    def dpc_ds(self, value):
        self._dpc_ds = value
    _dpc_dg = None
    @property
    def dpc_dg(self):
        return self._dpc_dg
    @dpc_dg.setter
    def dpc_dg(self, value):
        self._dpc_dg = value

    def _compute_capillary_pressure(self, sigma):
        self.pc = (
            sigma
            * self.omega[1]
            * np.sqrt(self.porosity)
        ) / (
            np.matmul(
                self.permeability_matrix ** (0.5),
                fractional_matrix_power(
                    1 - self.nso, 1 / self.omega[0]
                ),
            )
        )

    def _derivative_capillary_pressure(self, sigma, dsigma_dg):
        self.dpc_ds = self.pc / (self.omega[0]*(1-self.nso))
        self.dpc_dg = (self.pc/sigma)*dsigma_dg + self.dpc_ds

    def execute(self):
        """
        Primary function to execute Transport Equation calculations
        """
        
        # Initialize Parameters
        
        ## calculating interfacial tension
        sigma = self._surfactant.eval_IFT(self.surfactant_concentration)
        dsigma_dg = self._surfactant.eval_dIFT_dGamma(self.surfactant_concentration)
        
        ## calculating the residual saturations
        self._compute_residual_saturations(sigma)
        
        ## effective saturations
        self._compute_effective_saturations(self.water_saturation)
        
        ## mobility calculations 
        self._compute_lambda_a(RelativePermeabilityFormula.CoreyTypeEquation)
        self._compute_lambda_o(RelativePermeabilityFormula.CoreyTypeEquation)

        ## capillary number
        self._aqueous_capillary_number(sigma)
        self._oleic_capillary_number(sigma)

        # Water Saturation Computations
        wsat_old, wsat_modified, wsat_new = self._water_saturation_matrix_processing()
        self.water_saturation = wsat_new

        # Polymer Concentration Computations
        pconc_old, pconc_modified, pconc_new = self._polymer_concentration_matrix_processing(wsat_old)
        self.polymer_concentration = pconc_new


        # Surfactant Concentration Computation
        sconc_old, sconc_modified, sconc_new = self._surfactant_concentration_matrix_processing(wsat_old, wsat_modified)
        
        # Calculations for ROIP
        ocut = self.lambda_o[self.n-1, self.m-1] * self.total_flow / self.lambda_total[self.n-1,self.m-1]
        wcut = self.lambda_a[self.n-1, self.m-1] * self.total_flow / self.lambda_total[self.n-1,self.m-1]
        ROIP = 100*(np.sum(np.sum(1 - self.water_saturation)))/np.sum(np.ones((self.n*self.m,1)))
        
        return ocut, wcut, ROIP

    def _main_loop_computation(self, flag):
        """
        primary loop that will run for computations
        """
        idx = 1
        # setting matrices
        AAA = np.zeros((self.n * self.m, self.n * self.m))
        DDD = np.zeros((self.n * self.m, 1))
        fields = {
                1: self.water_saturation,
                2: self.polymer_concentration,
                3: self.surfactant_concentration
        }

        while (
            idx <= self.m * (self.n - 1) + 1
            and self.surfactant_concentration is not None
            and self.polymer_concentration is not None
        ):
            cnt = (idx - 1) // self.m  # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((self.n, self.m))
            AA = np.copy(BB)
            CC = np.copy(BB)
            DD = np.zeros((self.m, 1))


            for i in range(self.m):
                for j in range(self.n):
                    if i != j:
                        continue

                    field = fields[flag]

                    is_left = (i == 0)
                    is_right = (i == self.m - 1)

                    if is_left:
                        handler = self._leftmost_grid_calculations
                    elif is_right:
                        handler = self._rightmost_grid_calculations
                    else:
                        handler = self._interior_grid_calculations

                    handler(flag, idx, cnt, i, j, AA, BB, CC, DD, field)
            
            if cnt == 0:
                AAA[:self.n, : 2 * self.m] = np.hstack([BB, CC])
            elif cnt == self.n - 1:
                AAA[(self.m - 1) * self.n : self.m * self.n, (self.n - 2) * self.m : self.n * self.m] = np.hstack([AA, BB])
            else:
                AAA[cnt * self.n : (cnt + 1) * self.n, (cnt - 1) * self.m : (cnt + 2) * self.m] = np.hstack(
                    [AA, BB, CC]
                )

            DDD[cnt * self.m : (cnt + 1) * self.m] = DD
            idx += self.m
        
        TMM_flat, info = bicgstab(AAA, DDD, rtol=10 ** (-10), maxiter=600)
        if info != 0:
            import warnings
            warnings.warn(f"BiCGSTAB: convergence issue (info={info}) in water saturation solver")
        TMM_new = TMM_flat.reshape(self.m, self.n)

        return TMM_new
        

    def _water_saturation_matrix_processing(self):
        """
        run through computations for water saturation matrix
        """
        # redefine coordinates
        [xmod, ymod] = self._characteristic_coordinates(
            1,
            self.water_saturation,
            self.water_saturation,
        )
        
        #saving water saturation matrix and then modifying the current water saturation matrix for calcs
        water_saturation_old = np.copy(self.water_saturation)
        water_saturation_modified = self._matrix_reordering(np.copy(self.water_saturation), xmod, ymod)
        
        #computing the mobilities
        self._compute_lambda_a(RelativePermeabilityFormula.CoreyTypeEquation, water_saturation_modified)
        self._compute_lambda_o(RelativePermeabilityFormula.CoreyTypeEquation, water_saturation_modified)
        
        #computing derivatives of capillary pressure
        self._derivative_capillary_pressure(self._surfactant.eval_IFT, self._surfactant.eval_dIFT_dGamma)

        #executing main loop function
        water_saturation_new = self._main_loop_computation(1)

        return water_saturation_old, water_saturation_modified, water_saturation_new #[old matrix, modified matrix, new matrix]

    def _polymer_concentration_matrix_processing(self, water_saturation_old):
        """
        run through computations for the polymer concentration matrix
        """
        [xmod, ymod] = self._characteristic_coordinates(
            2,
            water_saturation_old,
            self.water_saturation,
        )

        #saving polymer matrix and then modifying the current polymer_concentration matrix for calcs
        polymer_concentration_old = np.copy(self.polymer_concentration)
        polymer_concentration_modified = self._matrix_reordering(np.copy(self.polymer_concentration), xmod, ymod)

        # executing main loop functions
        polymer_concentration_new = self._main_loop_computation(2)

        return polymer_concentration_old, polymer_concentration_modified, polymer_concentration_new # [old matrix, modified matrix, new matrix]

    def _surfactant_concentration_matrix_processing(self, water_saturation_old, water_saturation_modified):
        """
        run through computations for the surfactant concentration matrix
        """
        # redefine coordinates
        [xmod, ymod] = self._characteristic_coordinates(
            3,
            water_saturation_old,
            self.water_saturation,
        )

        #saving old surfactant concentration matrix and then modifying the current surfactant concentration matrix for calculations
        surfactant_concentration_old = np.copy(self.surfactant_concentration)
        surfactant_concentration_modified = self._matrix_reordering(np.copy(self.surfactant_concentration), xmod, ymod)

        #updating relevant coefficients using interpolated surfactant concentration matrix
        
        ##interfacial tension calculations
        sigma_modified = self._surfactant.eval_IFT(surfactant_concentration_modified)
        dsigma_dg_modified = self._surfactant.eval_dIFT_dGamma(surfactant_concentration_modified)

        ##compute residual saturations
        [swr, sor] = self._water.compute_residual_saturations(sigma_modified, self.pressure, self.velocity)
        self._compute_effective_saturations(water_saturation_modified)

        ##compute mobilities
        self._compute_lambda_a(RelativePermeabilityFormula.CoreyTypeEquation, water_saturation_modified)
        self._compute_lambda_o(RelativePermeabilityFormula.CoreyTypeEquation, water_saturation_modified)
        
        ##compute capillary pressues derivatives
        self._derivative_capillary_pressure(sigma_modified, dsigma_dg_modified)
        
        # running main loop functions
        surfactant_concentration_new = self._main_loop_computation(3)

        return surfactant_concentration_old, surfactant_concentration_modified, surfactant_concentration_new # [old matrix, new matrix]

    def _leftmost_grid_calculations(self, flag, row_index, cnt, i, j, AA, BB, CC, DD, TMM):
        """
        will conduct calculations related to the leftmost column of the grid at a particular row
        """
        match flag:
            case 1: #water saturaton matrix computations
                if row_index == 1:
                    DD[i] = (
                        (TMM[cnt, i] / self.dt_array[cnt, i])
                        + self.total_flow * (1 - self.fractional_flow[cnt, i])
                        + (
                            (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1]) / (self.dx**2)
                            + (self.dD_dg[cnt + 1, i] + self.dD_dg[cnt, i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i + 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt + 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt + 1, i]
                    )

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (self.dy**2)

                    BB[j, i] = (
                        1 / dt_array[cnt, i]
                        - (self.dD_ds[cnt+1, i] + self.dD_ds[cnt, i + 1]) / (self.dx**2)
                        - (self.dD_ds[cnt + 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )

                    BB[j, i + 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) * (self.dx**2)
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = (
                        (TMM[cnt, i] / self.dt_array[cnt, i])
                        + (
                            (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1]) / (self.dx**2)
                            + (self.dD_dg[cnt - 1, i] + self.dD_dg[cnt, i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i + 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt - 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt - 1, i]
                    )

                    AA[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt - 1, i]) / (self.dy**2)

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) / (self.dx**2)
                        - (self.dD_ds[cnt - 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )

                    BB[j, i + 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) / (self.dx**2)

                else:
                    DD[i] = (
                        (TMM[cnt, i] / self.dt_array[cnt, i])
                        - self.df_dc[cnt, i]
                        * (
                            self.velocity[cnt, i]
                            * (
                                self.polymer_concentration[cnt + 1, i]
                                - self.polymer_concentration[cnt, i]
                            )
                            / (2 * self.dy)
                        )
                        - self.df_dg[cnt, i]
                        * (
                            self.velocity[cnt, i]
                            * (
                                self.surfactant_concentration[cnt + 1, i]
                                - self.surfactant_concentration[cnt, i]
                            )
                            / (2 * self.dy)
                        )
                        + (
                            (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1]) / (self.dx**2)
                            + (
                                self.dD_dg[cnt - 1, i]
                                + 2 * self.dD_dg[cnt, i]
                                + self.dD_dg[cnt + 1, i]
                                )
                            / (2 * self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i + 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i + 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt + 1, i])
                        / (2 * self.dy**2)
                        * self.surfactant_concentration[cnt + 1, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt - 1, i])
                        / (2 * self.dy**2)
                        * self.surfactant_concentration[cnt - 1, i]
                    )

                    AA[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt - 1, i]) / (2 * self.dy**2)

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) / (self.dx**2)
                        - (
                            self.dD_ds[cnt - 1, i]
                            + 2 * self.dD_ds[cnt, i]
                            + self.dD_ds[cnt + 1, i]
                            )
                        / (2 * self.dy**2)
                    )
                    BB[j, i + 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) / (self.dx**2)

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (2 * self.dy**2)
            case 2: # polymer concentration matrix computations
                if row_index == 1: #bottom of grid
                    DD[i] = (
                        self.polymer_flow / self.water_saturation[cnt, i] + TMM[cnt, i] / self.dt_array[cnt, i]
                    )
                    BB[j, i] = 1 / self.dt_array[cnt, i] + self.total_flow / self.water_saturation[cnt, i]
            case 3: # surfactant concentration matrix computation
                F = self.D * self.dpc_dg / self.water_saturation
                if row_index == 1:
                    DD[i] = (
                        self.surfactant_flow / self.water_saturation[cnt, i] + TMM[cnt, i] / self.dt_array[cnt, i]
                    )
                    CC[j, i] = 2 * F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                        + self.total_flow / self.water_saturation[cnt, i]
                    )
                    BB[j, i + 1] = 2 * F[cnt][i] / (self.dx**2)
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = 2 * F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                        + self.total_flow / self.water_saturation[cnt, i]
                    )
                    BB[j, i + 1] = 2 * F[cnt, i] / (self.dx**2)
                else:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                    )
                    BB[j, i + 1] = 2 * F[cnt, i] / (self.dx**2)
                    CC[j, i] = F[cnt, i] / (self.dy**2)

    def _rightmost_calculations(self, flag, row_index, cnt, i, j, AA, BB, CC, DD, TMM):
        """
        will conduct calculatons related to the rightmost column of grid at a particular row
        """
        match flag:
            case 1: #water saturaton matrix computations
                if row_index == 1:
                    DD[i] = (
                        TMM[cnt, i] / self.dt_array[cnt, i]
                        + (
                            (self.dD_dg[cnt, i] + self.dD_dg[cnt][i - 1]) / (self.dx**2)
                            + (self.dD_dg[cnt + 1, i] + self.dD_dg[cnt, i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i - 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i - 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt + 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt + 1, i]
                    )

                    BB[j, i - 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dx**2)

                    BB[j, i] = (
                        1 / dt_array[cnt, i]
                        - (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dx**2)
                        - (self.dD_ds[cnt + 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (self.dy**2)
                    
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = (
                        (TMM[cnt, i] / self.dt_array[cnt, i])
                        + (
                            (self.dD_dg[cnt][i] + self.dD_dg[cnt][i - 1]) / (self.dx**2)
                            + (self.dD_dg[cnt - 1][i] + self.dD_dg[cnt][i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i - 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i - 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt - 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt - 1, i]
                    )

                    BB[j, i - 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dy**2)

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dx**2)
                        - (self.dD_ds[cnt - 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )

                    AA[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt - 1, i]) / (self.dy**2)
                else:
                    DD[i] = (
                        TMM[cnt, i] / self.dt_array[cnt, i]
                        - self.df_dc[cnt, i]
                        * (
                            self.velocity[cnt, i]
                            * (
                                self.polymer_concentration[cnt + 1, i]
                                - self.polymer_concentration[cnt, i]
                            )
                            / (2 * self.dy)
                        )
                        - self.df_dg[cnt, i]
                        * (
                            self.velocity[cnt, i]
                            * (
                                self.surfactant_concentration[cnt + 1, i]
                                - self.surfactant_concentration[cnt, i]
                            )
                            / (2 * self.dy)
                        )
                        + (
                            (self.dD_dg[cnt, i] + self.dD_dg[cnt, i - 1]) / (self.dx**2)
                            + (
                                self.dD_dg[cnt - 1, i]
                                + 2 * self.dD_dg[cnt, i]
                                + self.dD_dg[cnt + 1, i]
                            )
                            / (2 * self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt, i - 1])
                        / (self.dx**2)
                        * self.surfactant_concentration[cnt, i - 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt + 1, i])
                        / (2 * self.dy**2)
                        * self.surfactant_concentration[cnt + 1, i]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt - 1, i])
                        / (2 * self.dy**2)
                        * self.surfactant_concentration[cnt - 1, i]
                    )

                    AA[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt - 1, i]) / (2 * self.dy**2)

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dx**2)
                        - (
                            self.dD_ds[cnt - 1, i]
                            + 2 * self.dD_ds[cnt, i]
                            + self.dD_ds[cnt + 1, i]
                        )
                        / (2 * self.dy**2)
                    )
                    BB[j, i - 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i - 1]) / (self.dx**2)

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (2 * self.dy**2)

            case 2: # polymer concentration matrix computations
                if row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i] - self.total_flow * self.fractional_flow[cnt, i] / self.water_saturation[cnt, i]
                    )
            case 3: # surfactant concentration matrix computation
                F = self.D * self.dpc_dg / self.water_saturation
                if row_index == 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    CC[j, i] = 2 * F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                    )
                    BB[j, i - 1] = 2 * F[cnt, i] / (self.dx**2)
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = 2 * F[cnt, i] / (dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                        - ((self.total_flow * self.lambda_a[cnt][i]) / (self.lambda_total[cnt, i]))
                        / self.water_saturation[cnt, i]
                        + ((self.surfactant_flow * self.lambda_a[cnt][i]) / (self.lambda_total[cnt][i]))
                        / (self.water_saturation[cnt][i] * self.surfactant_concentration)
                    )
                    BB[j, i - 1] = 2 * F[cnt, i] / (self.dx**2)
                else:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = F[cnt][i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt][i]
                    )
                    BB[j, i - 1] = 2 * F[cnt, i] / (self.dx**2)
                    CC[j, i] = F[cnt, i] / (self.dy**2)

    def _interior_grid_calculations(self, flag, row_index, cnt, i, j, AA, BB, CC, DD, TMM):
        """
        will conduct calculations related to the interior column of grid at a particular row
        """
        match flag:
            case 1: #water saturaton matrix computations
                if row_index == 1:
                    DD[i] = (
                        TMM[cnt, i] / self.dt_array[cnt, i]
                        - self.df_dc[cnt, i]
                        * (
                            self.pressure[cnt, i]
                            * (
                                self.polymer_concentration[cnt, i + 1]
                                - self.polymer_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                        )
                        - self.df_dg[cnt, i]
                        * (
                            self.velocity[cnt, i]
                            * (
                                self.surfactant_concentration[cnt, i + 1]
                                - self.surfactant_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                        )
                        + (
                            (
                                self.dD_dg[cnt, i + 1]
                                + self.dD_dg[cnt, i - 1]
                                + 2 * self.dD_dg[cnt, i]
                            )
                            / (2 * self.dx**2)
                            + (self.dD_dg[cnt, i + 1] + self.dD_dg[cnt, i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i + 1] + self.dD_dg[cnt, i])
                        / (2 * self.dx**2)
                        * self.surfactant_concentration[cnt, i + 1]
                        - (self.dD_dg[cnt, i - 1] + self.dD_dg[cnt, i])
                        / (2 * self.dx**2)
                        * self.surfactant_concentration[cnt, i - 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt + 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt + 1, i]
                    )

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (
                            self.dD_ds[cnt, i + 1]
                            + self.dD_ds[cnt, i - 1]
                            + 2 * self.dD_ds[cnt, i]
                        )
                        / (2 * self.dx**2)
                        - (self.dD_ds[cnt + 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )
                    BB[j, i - 1] = (self.dD_ds[cnt, i - 1] + self.dD_ds[cnt, i]) / (
                        2 * self.dx**2
                    )
                    BB[j, i + 1] = (self.dD_ds[cnt, i + 1] + self.dD_ds[cnt, i]) / (
                        2 * self.dx**2
                    )

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (self.dy**2)
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = (
                        TMM[cnt, i] / self.dt_array[cnt, i]
                        - self.df_dc[cnt, i]
                        * (
                            self.pressure[cnt, i]
                            * (
                                self.polymer_concentration[cnt, i + 1]
                                - self.polymer_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                        )
                        - self.df_dg[cnt, i]
                        * (
                            self.pressure[cnt, i]
                            * (
                                self.surfactant_concentration[cnt, i + 1]
                                - self.surfactant_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                        )
                        + (
                            (
                                self.dD_dg[cnt, i + 1]
                                + self.dD_dg[cnt, i - 1]
                                + 2 * self.dD_dg[cnt, i]
                            )
                            / (2 * self.dx**2)
                            + (self.dD_dg[cnt, i + 1] + self.dD_dg[cnt, i]) / (self.dy**2)
                        )
                        * self.surfactant_concentration[cnt, i]
                        - (self.dD_dg[cnt, i + 1] + self.dD_dg[cnt, i])
                        / (2 * self.dx**2)
                        * self.surfactant_concentration[cnt, i + 1]
                        - (self.dD_dg[cnt, i - 1] + self.dD_dg[cnt, i])
                        / (2 * self.dx**2)
                        * self.surfactant_concentration[cnt, i - 1]
                        - (self.dD_dg[cnt, i] + self.dD_dg[cnt - 1, i])
                        / (self.dy**2)
                        * self.surfactant_concentration[cnt - 1, i]
                    )

                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - (
                            self.dD_ds[cnt, i + 1]
                            + self.dD_ds[cnt, i - 1]
                            + 2 * self.dD_ds[cnt, i]
                        )
                        / (2 * self.dx**2)
                        - (self.dD_ds[cnt - 1, i] + self.dD_ds[cnt, i]) / (self.dy**2)
                    )
                    BB[j, i - 1] = (self.dD_ds[cnt, i - 1] + self.dD_ds[cnt, i]) / (
                        2 * self.dx**2
                    )
                    BB[j, i + 1] = (self.dD_ds[cnt, i + 1] + self.dD_ds[cnt, i]) / (
                        2 * self.dx**2
                    )

                    AA[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt - 1, i]) / (self.dy**2)
                else:
                    DD[i] = (
                        (TMM[cnt, i] / self.dt_array[cnt, i])
                        - self.df_dc[cnt, i]
                        * (
                            self.pressure[cnt, i]
                            * (
                                self.polymer_concentration[cnt, i + 1]
                                - self.polymer_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                            + self.velocity[cnt, i]
                            * (
                                self.polymer_concentration[cnt + 1, i]
                                - self.polymer_concentration[cnt - 1, i]
                            )
                            / (2 * self.dy)
                        )
                        - self.df_dg[cnt, i]
                        * (
                            self.pressure[cnt, i]
                            * (
                                self.surfactant_concentration[cnt, i + 1]
                                - self.surfactant_concentration[cnt, i - 1]
                            )
                            / (2 * self.dx)
                            + self.velocity[cnt, i]
                            * (
                                self.surfactant_concentration[cnt + 1, i]
                                - self.surfactant_concentration[cnt - 1, i]
                            )
                            / (2 * self.dy)
                        )
                        - (
                            self.dD_dg[cnt, i + 1]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt, i + 1]
                                - self.surfactant_concentration[cnt, i]
                            )
                            - self.dD_dg[cnt, i - 1]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt, i - 1]
                                - self.surfactant_concentration[cnt, i]
                            )
                            + self.dD_dg[cnt, i - 1]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt, i - 1]
                                - self.surfactant_concentration[cnt, i + 1]
                            )
                            + self.dD_dg[cnt + 1, i]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt + 1, i]
                                - self.surfactant_concentration[cnt, i]
                            )
                            + self.dD_dg[cnt - 1, i]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt - 1, i]
                                - self.surfactant_concentration[cnt, i]
                            )
                            + self.dD_dg[cnt, i]
                            / (2 * self.dx**2)
                            * (
                                self.surfactant_concentration[cnt + 1, i]
                                - self.surfactant_concentration[cnt, i]
                            )
                        )
                    )
                    AA[j, i] = (self.dD_ds[cnt - 1, i] + self.dD_ds[cnt, i]) / (2 * self.dy**2)

                    CC[j, i] = (self.dD_ds[cnt, i] + self.dD_ds[cnt + 1, i]) / (2 * self.dy**2)

                    BB[j, i] = 1 / self.dt_array[cnt, i] - (
                        (1 / (2 * self.dx**2))
                        * (self.dD_ds[cnt, i] + 2 * self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1])
                        + (1 / (2 * self.dy**2))
                        * (
                            self.dD_ds[cnt - 1, i]
                            + 2 * self.dD_ds[cnt, i]
                            + self.dD_ds[cnt + 1, i]
                        )
                    )
                    BB[j, i + 1] = (self.dD_ds[cnt, i] + self.dD_ds[cnt, i + 1]) / (
                        2 * self.dx**2
                    )
                    BB[j, i - 1] = (self.dD_ds[cnt, i - 1] + self.dD_ds[cnt, i]) / (
                        2 * self.dx**2
                    )
            case 2: # polymer concentration matrix computations
                DD[i] = TMM[cnt][i] / self.dt_array[cnt, i]
                BB[j, i] = 1 / self.dt_array[cnt, i]
            case 3: # surfactant concentration matrix computation
                F = self.D * self.dpc_dg / self.water_saturation
                if row_index == 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    CC[j, i] = 2 * F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                    )
                    BB[j, i - 1] = F[cnt, i] / (self.dx**2)
                    BB[j, i + 1] = F[cnt, i] / (self.dx**2)
                elif row_index == (self.m) * (self.n - 1) + 1:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = 2 * F[cnt, i] / (self.dy**2)
                    BB[j, i + 1] = F[cnt, i] / (self.dx**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                    )
                    BB[j, i - 1] = F[cnt, i] / (self.dx**2)
                else:
                    DD[i] = TMM[cnt, i] / self.dt_array[cnt, i]
                    AA[j, i] = F[cnt, i] / (self.dy**2)
                    BB[j, i] = (
                        1 / self.dt_array[cnt, i]
                        - ((2 / (self.dx**2)) + (2 / (self.dy**2))) * F[cnt, i]
                    )
                    BB[j, i - 1] = F[cnt, i] / (self.dx**2)
                    BB[j, i + 1] = F[cnt, i] / (self.dx**2)
                    CC[j, i] = F[cnt, i] / (self.dy**2)

    def _matrix_reordering(self, transport_matrix, xmod, ymod):
        x1d = self.x[0, :]
        y1d = self.y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)

        # reorder Q if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            transport_matrix = transport_matrix[:, x_sort_idx]  # Sort columns of S
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            transport_matrix = transport_matrix[y_sort_idx, :]  # Sort rows of Q

        interp_func = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d), transport_matrix, method="linear", bounds_error=False, fill_value=None
        )

        query_points = np.stack([ymod.ravel(), xmod.ravel()], axis=-1)
        transport_matrix_modified = interp_func(query_points).reshape(xmod.shape)

        return transport_matrix_modified


    def _characteristic_coordinates(
        self,
        flag,
        old_water_saturation_matrix,
        new_water_saturation_matrix,
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
        assert self.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownWaterSaturationMatrix"
        )
        assert (
            self.surfactant_concentration is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownSurfactantConcentrationMatrix"
        )

        xjump = None
        yjump = None
        sold = old_water_saturation_matrix
        snew = new_water_saturation_matrix

        if flag == 1:
            xjump = self.x - self.df_ds * self.pressure * self.dt_array
            yjump = self.y - self.df_ds * self.velocity * self.dt_array
        elif flag == 2:
            # Calculate gradients
            sx, sy = self._get_gradient(sold)
            gx, gy = self._get_gradient(self.surfactant_concentration)

            xjump = (
                self.x
                - (
                    (self.fractional_flow / snew) * self.pressure
                    + (self.D * self.dpc_ds / snew) * sx
                    + (self.D * self.dpc_dg / snew) * gx
                )
                * self.dt_array
            )
            yjump = (
                self.y
                - (
                    (self.fractional_flow / snew) * self.velocity
                    + (self.D * self.dpc_ds / snew) * sy
                    + (self.D * self.dpc_dg / snew) * gy
                )
                * self.dt_array
            )
        elif flag == 3:
            sx, sy = self._get_gradient(sold)

            xjump = self.x - ((self.fractional_flow / snew) * self.pressure + (self.D * self.dpc_ds / snew) * sx) * self.dt_array
            yjump = self.y - ((self.fractional_flow / snew) * self.velocity + (self.D * self.dpc_ds / snew) * sy) * self.dt_array

        # Apply Neumann reflection conditions
        if xjump is None or yjump is None:
            raise SimulationCalcInputException(
                "SimulationInputException:UnknownXJumpYJumpMatrices"
            )

        xmod = np.copy(self.x)
        ymod = np.copy(self.y)

        for j in range(np.shape(self.y)[0]):
            for i in range(np.shape(self.x)[1]):
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

    def _get_gradient(self, vn):
        """
        (private method)

        Helper function to determine the gradients with respect to x and y dimensions
        
        Returns: (tuple[_Array[tuple[int, int], float64], NDArray[float64]])
        --------------------------------------------------------------------
            Tuple with px py which are numpy matrices that hold the gradient wrt to x and y dimensions
        """
        px = np.zeros((self.n + 1, self.m + 1))
        py = np.copy(px)

        for i in range(self.m + 1):
            for j in range(self.n + 1):
                if i != 0:
                    px[j, i] = (vn[j, i] - vn[j, i - 1]) / self.dx
                if i != self.m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i]) / self.dx
                if i != 0 and i != self.m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i - 1]) / (2 * self.dx)
                if j != 0:
                    py[j, i] = (vn[j, i] - vn[j - 1, i]) / self.dy
                if j != self.n:
                    py[j, i] = (vn[j + 1, i] - vn[j, i]) / self.dy
                if j != 0 and j != self.n:
                    py[j, i] = (vn[j + 1, i] - vn[j - 1, i]) / (2 * self.dy)

        return px, py


# class EllipticEquationSolver():
#     """
#     Will hold the relevant functions and parameters for solving the global pressure (u) and velocity (v)
#
#     (Not implemented for v1)
#     """
#     pass
