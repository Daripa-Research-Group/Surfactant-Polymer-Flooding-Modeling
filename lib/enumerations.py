"""
This python script contains the class definitions for the enumerations that are used within the simulation runs

The methods of this class were derived from the MATLAB Surfactant-Polymer Flooding Code developed by
Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""
from enum import Enum


class SimulationConstants(Enum):
    """
    Simulation constants
    (Taken from 2017 paper "Modeling and simulation of surfactant–polymer flooding using a new hybrid method")
    """

    Water_Viscosity = 1.26
    Water_Density = 1000  # kg/m^3
    Oil_Viscosity = 10

    Initial_Residual_Water_Saturation = 0.79
    Resid_Aqueous_Phase_Saturation_Initial = 0.1  # wetting phase
    Resid_Oleic_Phase_Saturation_Initial = 0.3  # non-wetting phase

    Aqueous_Phase_Critical_Capillary_Num = 10 ** (-5)
    Oleic_Phase_Critical_Capillary_Num = 10 ** (-5)

    Capillary_Pressure_Param_1 = 0.1  # omega_1
    Capillary_Pressure_Param_2 = 0.4  # omega_2

    Injection_Rate = 200
    Time_Step = 1 / 50
    Grid_Size = 29
    Source_Flow_Magnitude = 120000


class PolymerList(Enum):
    """
    List of Polymers that can be selected for the simulation runs
    """
    Xanthane = (1, 1500, [3.05428284, -0.27294817], [1.15410398e-04, 2.04937780e00])
    Schizophyllan = (2, 1300, [4.86265534, -0.41570227], [0.03647214, 1.32175949])
    No_Polymer = (3, 0, [0, 0], [0, 0])

    @property
    def Id(self):
        return self.value[0]

    @property
    def Density(self):
        return self.value[1]

    @property
    def n_coeff(self):
        return self.value[2]

    @property
    def e_coeff(self):
        return self.value[3]

    @classmethod
    def get_by_value(cls, value):
        member = next((member for member in cls if member.value[0] == value), None)
        return member


class SurfactantList(Enum):
    """
    List of Surfactants that can be selected for the simulation runs
    """
    Alkyl_Ether_Sulfate = [1, lambda GG: 10.001 / (GG + 1),lambda GG: (-10.001)/((GG + 1) ** 2)]
    No_Surfactant = [2, None, None]

    @property
    def Id(self):
        return self.value[0]

    @property
    def IFT_equation(self):
        return self.value[1]

    @property
    def derivative_IFT_equation(self):
        return self.value[2]

    @classmethod
    def get_by_value(cls, value):
        member = next((member for member in cls if member.value[0] == value), None)
        return member



class ModelType(Enum):
    """
    The simulation model types that can be selected
    """
    No_Shear_Thinning = 1
    Sourav_Implementation = 2
    Shear_Thinning_On = 3


class PlotType(Enum):
    """
    Selection of types of plots that can be created for the user
    """
    Saturation_Plot = 1
    Polymer_Concentration_Plot = 2
    Surfactant_Concentration_Plot = 3


class ResevoirGeometry(Enum):
    """
    Selection of the geometry of the resevoir for the simulation
    """
    Rectilinear = 1
    Quarter_Five_Spot = 2


class PermeabilityType(Enum):
    """
    Selection of the permeability profile for each of the simulation runs
    """
    Homogenous = 1
    Heterogenous = 2
