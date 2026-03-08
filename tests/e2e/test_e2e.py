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


def load_true_values(sim_id):
    return {
        "COC": np.load(os.path.join(TRUE_VALUES_DIR, f"sim_{sim_id}_COC.npy")),
        "MFW": np.load(os.path.join(TRUE_VALUES_DIR, f"sim_{sim_id}_MFW.npy")),
    }


@pytest.mark.e2e
@pytest.mark.parametrize(
    "simulation_id, model_type, reservoir_geometry, permeability, polymer_type, polymer_concentration, surfactant_type, surfactant_concentration",
    [
        (
            1,
            MODEL["Shear_Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Homogeneous"],
            POLYMER["Xanthane"],
            0.001,
            SURFACTANT["No_Surfactant"],
            0,
        ),
        (
            2,
            MODEL["No_Shear_Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Homogeneous"],
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

    expected = load_true_values(simulation_id)

    assert_allclose(
        sim.COC,
        expected["COC"],
        rtol=1e-5,
        atol=1e-8,
        err_msg=f"COC mismatch for simulation_id={simulation_id}",
    )

    assert_allclose(
        np.array(sim.MFW),
        expected["MFW"],
        rtol=1e-5,
        atol=1e-8,
        err_msg=f"MFW mismatch for simulation_id={simulation_id}",
    )
