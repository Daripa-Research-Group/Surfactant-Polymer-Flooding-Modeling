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
from enumerations import *
from Exceptions import *
import numpy as np

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
        return self._total_flow*self._water.concetration_scalar
    @property
    def surfactant_flow(self):
        return self._total_flow*self._surfactant.concentration
    ## Water Saturation Matrix
    @property
    def water_saturation(self):
        return self._water.water_saturation
    ## Aqueous Viscosity
    @property
    def aqueous_viscosity(self):
        return self._water.viscosity_array
    ## Surfactant concentration matrix
    @property
    def surfactant_concentration(self):
        return self._surfactant.concentration
    ## Polymer concentration matrix
    @property
    def polymer_concentration(self):
        return self._polymer.concentration_matrix
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
        self._velocity = velocity
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
    
    # Functions for parameter definitions
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
        nsw: normalized aqueous saturation
        """
        return self._nsw
    @nsw.setter
    def nsw(self, value):
        self._nsw = value

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
        nso: normalized oleic saturation
        """
        return self._nso
    @nso.setter
    def nso(self, value):
        self._nso = value


    
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
        pass

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

    def _compute_lambda_a(
        self,
        rel_permeability_formula: RelativePermeabilityFormula,
        modified_water_saturation: np.ndarray | None = None,
    ):
        if(modified_water_saturation is not None):
            self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                True, 
                RelativePermeabilityFormula.CoreyTypeEquation.value, 
                modified_water_saturation
                )
        else:
            self._water.compute_mobility(
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
            self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                False, 
                RelativePermeabilityFormula.CoreyTypeEquation.value, 
                modified_water_saturation
                )
        else:
            self._water.compute_mobility(
                self.polymer_concentration, 
                self.sor, 
                self.swr, 
                False, 
                RelativePermeabilityFormula.CoreyTypeEquation.value 
                )

    ## Capillary Number (ratio of viscous forces to surface tension forces)
    def _aqueous_capillary_number(self):
        pass

    def _oleic_capillary_number(self):
        pass
    






    def execute(self):
        pass

    def _initialize_matrices(self):
        """
        initialize relevant matrices needed for computation
        """
        pass

    def _main_loop_computation(self):
        """
        primary loop that will run for computations
        """
        pass

    def _saturation_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the saturation matrix
        """
        pass

    def _surfactant_concentration_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the surfactant concentration matrix
        """
        pass

    def _polymer_concentration_matrix_preprocessing(self):
        """
        will conduct any preprocessing prior to computing the polymer concentration matrix
        """
        pass

    def _bottom_grid_calculations(self):
        """
        will conduct calculations related to the bottom of the grid
        """
        pass

    def _top_grid_calculations(self):
        """
        will conduct calculatons related to the top of the grid
        """
        pass

    def _interior_grid_calculations(self):
        """
        will conduct calculations related to the interior of the grid
        """
        pass


# class EllipticEquationSolver():
#     """
#     Will hold the relevant functions and parameters for solving the global pressure (u) and velocity (v)
#
#     (Not implemented for v1)
#     """
#     pass
