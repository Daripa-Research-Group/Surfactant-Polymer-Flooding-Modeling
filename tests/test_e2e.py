import pytest
from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.surfactant import Surfactant
from lib.para import Box
from lib.enumerations import (
    SimulationConstants,
    PolymerList,
    ModelType,
    ReservoirGeometry,
    PermeabilityType,
    PlotType,
    SurfactantList,
)

# default polymer configurations
XANTHANE = Polymer(
    PolymerList.Xanthane, 0.001, PolymerList.e_coeff, PolymerList.n_coeff, None, None
)
SCHIZOPHYLLAN = Polymer(
    PolymerList.Schizophyllan,
    0.001,
    PolymerList.e_coeff,
    PolymerList.n_coeff,
    None,
    None,
)

# default surfactant configurations
ALKYL_ETHER_SULFATE = Surfactant(
    SurfactantList.Alkyl_Ether_Sulfate, 0.09, None, None, None
)

INITIAL_WATER_SATURATION = 0.79
box = Box()


@pytest.mark.integration
@pytest.mark.parametrize(
    "sim_id, size_of_grid, polymer, surfactant, init_water_saturation, reservoir_geometry, permeability_type, mesh_grid, model_type, plot_type",
    [
        (
            1,
            32,
            XANTHANE,
            ALKYL_ETHER_SULFATE,
            INITIAL_WATER_SATURATION,
            ReservoirGeometry.Rectilinear,
            PermeabilityType.Homogenous,
            box,
            ModelType.No_Shear_Thinning,
            PlotType.Saturation_Plot,
        ),
        # (2, 32, XANTHANE, ALKYL_ETHER_SULFATE, INITIAL_WATER_SATURATION, ReservoirGeometry.Rectilinear, PermeabilityType.Heterogenous, box, ModelType.No_Shear_Thinning, PlotType.Saturation_Plot),
    ],
)
def test_e2e(
    sim_id,
    size_of_grid,
    polymer,
    surfactant,
    init_water_saturation,
    reservoir_geometry,
    permeability_type,
    mesh_grid,
    model_type,
    plot_type,
):
    """Performs end-to-end simulations tests verifying validity of results"""

    simulation = Simulation(
        sim_id,
        size_of_grid,
        polymer,
        surfactant,
        init_water_saturation,
        reservoir_geometry,
        permeability_type,
        mesh_grid,
        model_type,
        plot_type,
    )


test_e2e(
    1,
    32,
    XANTHANE,
    ALKYL_ETHER_SULFATE,
    INITIAL_WATER_SATURATION,
    ReservoirGeometry.Rectilinear,
    PermeabilityType.Homogenous,
    box,
    ModelType.No_Shear_Thinning,
    PlotType.Saturation_Plot,
)
