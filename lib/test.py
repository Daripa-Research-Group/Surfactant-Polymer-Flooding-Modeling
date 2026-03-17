"""
Sandbox script for developer testing

Use this python script when implementing new features
"""

import os

import numpy as np
from .grid import Grid
from .enumerations import (
    ModelType,
    PolymerList,
    SurfactantList,
    PermeabilityType,
    ResevoirGeometry,
    SimulationConstants,
)
from .Exceptions import SimulationCalcInputException, UserInputException
import .simulation

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
    "model_type": MODEL["No_Shear_Thinning"],
    "reservoir_geometry": GEOMETRY["Rectilinear"],
    "permeability": PERMEABILITY["Homogeneous"],
    "polymer_type": POLYMER["Xanthane"],
    "polymer_concentration": 0.001,
    "surfactant_type": SURFACTANT["Alkyl_Ether_Sulfate"],
    "surfactant_concentration": 0,
}
# FIXME: discrepancy in polymer viscosity matrix. Checked the u, v, x, and y (x and y from 'Grid' class)
sim1 = simulation.Simulation(user_input_dict=user_dict)
sim1.run()
# print(f'polymer viscosity matrix{sim1.polymer.viscosity_matrix}')
# print(f'aqueous viscosity matrix{sim1.water.viscosity_array}')
