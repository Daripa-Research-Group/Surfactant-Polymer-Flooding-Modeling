import os

import numpy as np
from grid import Grid
from enumerations import (
    ModelType,
    PolymerList,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from Exceptions import SimulationCalcInputException, UserInputException
import simulation

MODEL = {
    "No_Shear_Thinning": 1,
    "Sourav_Implementation": 2,
    "Shear_Thinning": 3,
}

GEOMETRY = {
    "Rectilinear": 1,
    "Quarter Five Spot": 2,
}

PERMEABILITY = {
    "Homogeneous": 1,
    "Heterogeneous": 2,
}

POLYMER = {
    "Xanthane": 1,
    "Schizophyllan": 2,
    "No Polymer": 3,
}

SURFACTANT = {
    "Alkyl_Ether_Sulfate": 1,
    "No_Surfactant": 2,
}

# Making the simulation object:
user_dict = {
        "simulation_id": 1,
        "model_type": MODEL["Shear_Thinning"],
        "reservoir_geometry": GEOMETRY["Rectilinear"],
        "permeability": PERMEABILITY["Heterogeneous"],
        "polymer_type": POLYMER["Xanthane"],
        "polymer_concentration": 0.001,
        "surfactant_type": SURFACTANT["Alkyl_Ether_Sulfate"],
        "surfactant_concentration": 0,
}
sim1 = simulation.Simulation(user_input_dict=user_dict)

print(f"{sim1.phi}")
