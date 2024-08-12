'''
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40
'''

#import statements
import numpy as np
import matplotlib.pyplot as plt
#from get_phi_test import get_phi_test
#from KKdef import KKdef
#from s0c0 import s0c0
#from compvis import compvis


# declaring global variables to reduce computational time in passing through function calls
global c0_array, miuo, miuw, miup, swr0, sor0, dt, KK, s0, c0, g0, beta1, viscosityFlag, shearFlag, miup_array, polymerType, nsim, sog, k




"""
Runs the simulation

Parameters: None

Return: None
"""
def main():
    # clearing variables and setting all initial paramters
    plt.close("all")
    np.set_printoptions(precision=15, suppress=True) #adjusts print setting for improved readability
    
    #getting simulation parameters
    sim_params = simulation_parameters()
    pass



"""
Initiates the parameters for the simulation

Parameters: None

Return: Dict

"""
def simulation_parameters():
    '''
    #### python dictionary for defining mesh and grid points
    ############################################### ##############
    '''
    para_box = {
        "left" : 0,
        "right" : 1,
        "bottom" : 0,
        "top" : 1
    }

    #collecting user input for simulation
    print("1 = Original Model from JCP paper(Shear thinning off)")
    print("2 = Sourav's Model")
    print("3 = Dynamic Viscosity Model (Non-Newtonian Model(Shear thinning on)")
    model_type = int(input("Enter which computation model to use for this simulation: "))
    print("----------------------------------------------")

    figure_num = int(input("enter figure number being simulated here: "))
    print("----------------------------------------------")

    c0iter_val = float(input("Enter value for c0iter(0.001 as default): "))
    g0iter_val = float(input("Enter value for g0iter_val: "))
    print("----------------------------------------------")

    print("0 = Saturation(UU) Plot")
    print("1 = Polymer(CC) Plot")
    print("2 = Surfactant(GG) Plot")
    plot_type_input = int(input("Enter the type of plot desired: "))
    print("----------------------------------------------")

    print("0 = xanthane")
    print("1 = schizophyllan")
    polymer_input = int(input("Enter the desired polymer for the simulation: "))
    print("----------------------------------------------")

    print("1 = homogeneous with magnitude 1000")
    print("2 = continuous heterogeneous function")
    print("3 = impermeable block inclusion at the center")
    print("4 = impermeable block inclusion off the center")
    print("5 = Upper Ness from SPE10 model sections")
    print("6 = Tabert from SPE10 model sections")
    permeability_input = int(input("Enter the corresponding permeability flag for the simulation: "))
    print("----------------------------------------------")

    alpha_val = float(input("Enter value for alpha (range from [0.5 - 5] using 0.5 step increments):"))
    print("----------------------------------------------")

    #initializing key variables for simulation 
    
    #return will be a dictionary with input parameter information
    simulation_parameters_dictionary = {
            "figure_number" : figure_num,
            "figure_param" : para_box,
            "model_type" : model_type,
            "c0iter" : c0iter_val,
            "g0iter" : g0iter_val,
            "plot_type" : plot_type_input,
            "polymer" : polymer_input,
            "permeability_flag" : permeability_input,
            "alpha_value" : alpha_val
        }

    return simulation_parameters_dictionary



if __name__ == "__main__":
    main()

