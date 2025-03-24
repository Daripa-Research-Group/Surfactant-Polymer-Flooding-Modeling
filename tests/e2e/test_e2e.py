import pytest
import numpy as np
from numpy.testing import assert_allclose
from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.surfactant import Surfactant
from lib.enumerations import (
    PolymerList,
    ModelType,
    PermeabilityType,
    PlotType,
    SurfactantList,
    ResevoirGeometry,
    SimulationConstants,
)

MODEL = {
    "No Shear Thinning": 1,
    "Sourav Implementation": 2,
    "Shear Thinning": 3,
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
    "Alkyl Ether Sulfate": 1,
    "No Surfactant": 2,
}


def load_true_values(sim_id):
    base = f"true_values/sim_{sim_id}_"
    return {
        "COC": np.load(base + "COC.npy"),
        "MFW": np.load(base + "MFW.npy"),
    }


@pytest.mark.e2e
@pytest.mark.parametrize(
    "simulation_id, model_type, reservoir_geometry, permeability, polymer_type, polymer_concentration, surfactant_type, surfactant_concentration",
    [
        (
            1,
            MODEL["Shear Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Homogeneous"],
            POLYMER["Xanthane"],
            0.001,
            SURFACTANT["No Surfactant"],
            0,
        ),
        (
            2,
            MODEL["No Shear Thinning"],
            GEOMETRY["Rectilinear"],
            PERMEABILITY["Homogeneous"],
            POLYMER["Schizophyllan"],
            0.001,
            SURFACTANT["No Surfactant"],
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
    plot_type = PlotType.Saturation_Plot  # TODO: MAKE DYNAMIC

    model_type = ModelType(model_type)
    reservoir_geometry = ResevoirGeometry(reservoir_geometry)
    permeability_flag = PermeabilityType(permeability)
    polymer_type = PolymerList.get_by_value(polymer_type)
    polymer_concentration = polymer_concentration
    surfactant_type = SurfactantList(surfactant_type)
    surfactant_concentration = surfactant_concentration

    polymer_obj = Polymer(
        name=polymer_type,
        initial_concentration=polymer_concentration,
        e_coeff=polymer_type.e_coeff,
        n_coeff=polymer_type.n_coeff,
    )

    surfactant_obj = Surfactant(
        name=surfactant_type,
        initial_concentration=surfactant_concentration,
        IFT_conc_equ=lambda GG: 10.001 / (GG + 1)
        - 0.001,  # TODO: make dynamic depending on the surfactant type
        derivative_IFT_conc_equ=lambda GG: (-10.001)
        / ((GG + 1) ** 2),  # TODO: make dynamic depending on the surfactant type
    )

    SOG = SimulationConstants.Grid_Size.value

    simulation = Simulation(
        sim_id=simulation_id,
        size_of_grid=SOG,
        polymer=polymer_obj,
        surfactant=surfactant_obj,
        resevoir_geometry=reservoir_geometry,
        permeability_flg=permeability_flag,
        mdl_id=model_type,
        plt_type=plot_type,
    )

    simulation_output = simulation.execute_simulation()
    expected_output = load_true_values(simulation_id)

    # compare simulation output with expected output
    assert_allclose(
        simulation_output["COC"],
        expected_output["COC"],
        rtol=1e-5,
        atol=1e-8,
        err_msg=f"COC mismatch in e2e test for simulation id = {simulation_id}",
    )
