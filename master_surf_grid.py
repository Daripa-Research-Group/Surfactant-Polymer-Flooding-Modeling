"""
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40

This code was derived from EOR repository developed by Sourav Dutta and Prabir Daripa

@author: Bhargav Akula Ramesh Kumar
"""

#### IMPORT STATEMENTS
import tkinter as tk
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
from user_input.gui import UserInputGUI
from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.surfactant import Surfactant
from lib.enumerations import (
    ModelType,
    PlotType,
    PermeabilityType,
    PolymerList,
    SurfactantList,
    ResevoirGeometry,
    SimulationConstants,
)


def sim_condition_initialization(simulation_ID: int, usr_input_dict: dict) -> dict:
    plot_type = PlotType.Saturation_Plot  # TODO: MAKE DYNAMIC

    model_type = ModelType(usr_input_dict["model_type"])
    assert model_type is not None, "Model Type not selected. Please try again"
    reservoir_geometry = ResevoirGeometry(usr_input_dict["reservoir_geometry"])
    assert reservoir_geometry is not None, "Resevoir Geometry  not selected. Please try again"
    permeability_flag = PermeabilityType(usr_input_dict["permeability"])
    assert permeability_flag is not None, "Permeability not properly selected. Please try again"
    polymer_type = PolymerList.get_by_value(usr_input_dict["polymer_type"])
    assert polymer_type is not None, "Polymer not properly selected. Please try again"
    polymer_concentration = usr_input_dict["polymer_concentration"]
    assert polymer_concentration is not None, "Polymer concentration not given. Please try again"
    surfactant_type = SurfactantList(usr_input_dict["surfactant_type"])
    assert surfactant_type is not None, "Surfactant not properly selected. Please try again"
    surfactant_concentration = usr_input_dict["surfactant_concentration"]
    assert surfactant_concentration is not None, "Surfactant concentration not given. Please try again"
    
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
        sim_id=simulation_ID,
        size_of_grid=SOG,
        polymer=polymer_obj,
        surfactant=surfactant_obj,
        resevoir_geometry=reservoir_geometry,
        permeability_flg=permeability_flag,
        mdl_id=model_type,
        plt_type=plot_type,
    )

    simulation_outputs = simulation.execute_simulation()

    return simulation_outputs


def main() -> None:
    root = tk.Tk()
    app = UserInputGUI(root)
    root.mainloop()

    user_input = app.get_input()
    for index, simulation in enumerate(user_input):
        sim_condition_initialization(index + 1, simulation)


if __name__ == "__main__":
    main()
