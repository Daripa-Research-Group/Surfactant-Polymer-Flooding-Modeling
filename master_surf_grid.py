'''
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40

This code was derived from EOR repository developed by Sourav Dutta and Prabir Daripa

@author: Bhargav Akula Ramesh Kumar
'''

#### IMPORT STATEMENTS
import numpy as np
import tkinter as tk
import seaborn as sb
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from user_input.gui import UserInputGUI
from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.surfactant import Surfactant
from lib.enumerations import ModelType, PlotType, PermeabilityType, PolymerList, SurfactantList, ResevoirGeometry, SimulationConstants

def sim_condition_initialization(simulation_ID : int, usr_input_dict : dict) -> dict:
    # dummy configurations for now
    # TODO: Update user input collection to adhere to class structure
    model_type = ModelType.No_Shear_Thinning
    plot_type = PlotType.Saturation_Plot
    permeability_flag = PermeabilityType.Homogenous
    polymer_type = PolymerList.Xanthane
    # model_type = ModelType(usr_input_dict['model_type'])  
    # plot_type = PlotType(usr_input_dict['plot_type'])  
    # permeability_flag = PermeabilityType(usr_input_dict['permeability_flag'])  
    # polymer_type = PolymerList(usr_input_dict['polymer']) 

    polymer = Polymer(
        name=polymer_type,
        initial_concentration=polymer_type.Density,
        e_coeff=polymer_type.e_coeff,
        n_coeff=polymer_type.n_coeff
    )

    surfactant = SurfactantList.No_Surfactant  # Defaulting to no surfactant
    surfactant_obj = Surfactant(
        name=surfactant,
        initial_concentration=0,  # No surfactant present
        IFT_conc_equ=lambda x: 0,  # Dummy lambda function
        derivative_IFT_conc_equ=lambda x: 0  # Dummy lambda function
    )

    reservoir_geometry = ResevoirGeometry.Rectilinear

    SOG = SimulationConstants.Grid_Size.value

    simulation = Simulation(
        sim_id=simulation_ID,
        size_of_grid=SOG,
        polymer=polymer,
        surfactant=surfactant_obj,
        resevoir_geometry=reservoir_geometry,
        permeability_flg=permeability_flag,
        mdl_id=model_type,
        plt_type=plot_type
    )
    
    simulation.execute_simulation()

    return simulation

def main() -> None:
    root = tk.Tk()
    app = UserInputGUI(root)
    root.mainloop()
    
    user_input = app.get_input()
    for index, simulation in enumerate(user_input):
        sim_condition_initialization(index + 1, simulation)



if(__name__ == "__main__"):
    main()

