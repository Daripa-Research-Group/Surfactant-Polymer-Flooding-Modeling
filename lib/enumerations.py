from enum import Enum


class SimulationConstants(Enum):
    """
    Simulation constants
    (Taken from 2017 paper "Modeling and simulation of surfactant–polymer flooding using a new hybrid method")
    """

    Water_Viscosity = 1.26
    Water_Density = 1000  # kg/m^3
    Oil_Viscosity = 10

    Initial_Resident_Water_Saturation = 0.21
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


class SurfactantList(Enum):
    Alkyl_Ether_Sulfate = 1
    No_Surfactant = 2


class ModelType(Enum):
    No_Shear_Thinning = 1
    Sourav_Implementation = 2
    Shear_Thinning_On = 3


class PlotType(Enum):
    Saturation_Plot = 1
    Polymer_Concentration_Plot = 2
    Surfactant_Concentration_Plot = 3


class ResevoirGeometry(Enum):
    Rectilinear = 1
    Quarter_Five_Spot = 2


class PermeabilityType(Enum):
    Homogenous = 1
    Heterogenous = 2
