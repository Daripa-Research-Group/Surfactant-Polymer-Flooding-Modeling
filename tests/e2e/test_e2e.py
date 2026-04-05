import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../lib")))

import pytest
import numpy as np
from numpy.testing import assert_allclose
from lib.simulation import Simulation

MODEL = {
    "No_Shear_Thinning": 1,
    "Sourav_Implementation": 2,
    "Shear_Thinning": 3,
}

GEOMETRY = {
    "Rectilinear": 1,
    "Quarter_Five_Spot": 2,
}

PERMEABILITY = {
    "Homogeneous": 1,
    "Heterogeneous": 2,
}

POLYMER = {
    "Xanthane": 1,
    "Schizophyllan": 2,
    "No_Polymer": 3,
}

SURFACTANT = {
    "Alkyl_Ether_Sulfate": 1,
    "No_Surfactant": 2,
}

TRUE_VALUES_DIR = os.path.join(os.path.dirname(__file__), "true_values")

SIM_ID_TO_DIR = {
    1: "simulation_one",
    2: "simulation_two",
}


def load_coc(sim_id):
    path = os.path.join(TRUE_VALUES_DIR, SIM_ID_TO_DIR[sim_id], "COC.csv")
    return np.loadtxt(path).flatten()


@pytest.mark.e2e
@pytest.mark.parametrize(
    "simulation_id, model_type, reservoir_geometry, permeability, polymer_type, polymer_concentration, surfactant_type, surfactant_concentration",
    [
        (
            1,
            MODEL["Shear_Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Heterogeneous"],
            POLYMER["Xanthane"],
            0.001,
            SURFACTANT["No_Surfactant"],
            0,
        ),
        (
            2,
            MODEL["No_Shear_Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Heterogeneous"],
            POLYMER["Schizophyllan"],
            0.001,
            SURFACTANT["No_Surfactant"],
            0,
        ),
    ],
)
def test_e2e(
    simulation_id,
    model_type,
    reservoir_geometry,
    permeability,
    polymer_type,
    polymer_concentration,
    surfactant_type,
    surfactant_concentration,
):
    user_dict = {
        "simulation_id": simulation_id,
        "model_type": model_type,
        "reservoir_geometry": reservoir_geometry,
        "permeability": permeability,
        "polymer_type": polymer_type,
        "polymer_concentration": polymer_concentration,
        "surfactant_type": surfactant_type,
        "surfactant_concentration": surfactant_concentration,
    }

    sim = Simulation(user_input_dict=user_dict)
    sim.run()

    expected_coc = load_coc(simulation_id)

    print(f"sim.coc: {sim.COC[-1][-1]} and expected_coc: {expected_coc[-1]}")
    assert abs(sim.COC[-1] - expected_coc[-1]) < 1e-5, f"Final COC mismatch for simulation_id={simulation_id}"

    # assert_allclose(
    #     np.asarray(sim.COC).flatten()[-1],
    #     expected_coc[-1],
    #     rtol=1e-5,
    #     atol=1e-8,
    #     err_msg=f"Final COC mismatch for simulation_id={simulation_id}",
    # )
