'''
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40

This code was derived from EOR repository developed by Sourav Dutta and Prabir Daripa

@author: Bhargav Akula Ramesh Kumar
'''

#### IMPORT STATEMENTS
import numpy as np
# import matplotlib.pyplot as plt
import seaborn as sb
import math
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from lib.get_phi_test import get_phi_test
from lib.KKdef import KKdef
from lib.surfactant_polymer_conc_initial import inital_polymer_surfactant_concentration
from lib.Exceptions import OutOfRangeError
from lib.para import Para, Box
from lib.compvis import compvis


#initializing global variables
NSIM = 1
SOG = 29
K = 1
MASS_FLOW_MAGNITUDE = 12000 #represented as "src" in MATLAB version of code
TIME_MATRIX = np.zeros((1,3))
ITERATIONS = np.zeros((1,3))
CFL = 1
TIME_STOP = 500
TIME_START = 0




"""
Runs the simulation

Parameters: None

Return: None
"""
def main():
    # # clearing variables and setting all initial paramters
    # plt.close("all")
    # np.set_printoptions(precision=15, suppress=True) #adjusts print setting for improved readability
    # 
    # #getting simulation parameters
    # input_dict = sim_inputs()
    pass



"""
Allows the user to provide information on the type of simulation they want to run

Parameters: None

Return: Dict

"""
def usr_inputs():
    #need to implement feature to take in file with user inputs instead if user wants to do multiple runs

    #collecting user input for simulation
    print("1 = Original Model from JCP paper(Shear thinning off)")
    print("2 = Sourav's Model")
    print("3 = Dynamic Viscosity Model (Non-Newtonian Model(Shear thinning on)")
    while True:
        try:
            model_type = int(input("Enter which computation model to use for this simulation: "))
            if not (1 <= model_type <= 3):
                raise OutOfRangeError(model_type)
            break

        except OutOfRangeError as e:
              print(e)
    print("----------------------------------------------")

    figure_num = int(input("enter figure number being simulated here: "))
    print("----------------------------------------------")

    c0iter_val = float(input("Enter value for c0iter(0.001 as default): "))
    g0iter_val = float(input("Enter value for g0iter_val: "))
    print("----------------------------------------------")

    print("0 = Saturation(UU) Plot")
    print("1 = Polymer(CC) Plot")
    print("2 = Surfactant(GG) Plot")
    while True:
        try:
            plot_type_input = int(input("Enter the type of plot desired: "))
            if not(0 <= plot_type_input <= 2):
                raise OutOfRangeError(plot_type_input)
            break
        except OutOfRangeError as e:
            print(e)
    print("----------------------------------------------")

    print("0 = xanthane")
    print("1 = schizophyllan")
    while True:
        try:    
            polymer_input = int(input("Enter the desired polymer for the simulation: "))
            if not(0 <= polymer_input <= 1):
                raise OutOfRangeError(polymer_input)
            break
        except OutOfRangeError as e:
            print(e)
    print("----------------------------------------------")

    print("1 = homogeneous with magnitude 1000")
    print("2 = continuous heterogeneous function")
    print("3 = impermeable block inclusion at the center")
    print("4 = impermeable block inclusion off the center")
    print("5 = Upper Ness from SPE10 model sections")
    print("6 = Tabert from SPE10 model sections")
    while True:
        try:
            permeability_input = int(input("Enter the corresponding permeability flag for the simulation: "))
            if not(1 <= permeability_input <= 6):
                raise OutOfRangeError(permeability_input)
            break
        except OutOfRangeError as e:
            print(e)
    print("----------------------------------------------")
    while True:
        try:
            alpha_val = float(input("Enter value for alpha (range from [0.5 - 5] using 0.5 step increments):"))
            if not(0.5 <= alpha_val <= 5):
                raise OutOfRangeError(alpha_val)
            break
        except OutOfRangeError as e:
            print(e)
    print("----------------------------------------------")

    
    #return will be a dictionary with input parameter information
    usr_input_dict = {
            "figure_number" : figure_num,
            "model_type" : model_type,
            "c0iter" : c0iter_val,
            "g0iter" : g0iter_val,
            "plot_type" : plot_type_input,
            "polymer" : polymer_input,
            "permeability_flag" : permeability_input,
            "alpha_value" : alpha_val,
            "nsim" : 1,
            'Todd_Longstaff_mixing_parameter' : 0,
            'shear_flag' : 0
            #"size_of_grid" : [SOG, SOG, SOG, SOG] <-- use when implementing multiple runs
        }

    return usr_input_dict


def sim_condition_initialization(usr_input_dict):
    #### Initializing grid points
    mesh_dimensions = Box()
    mesh_dimensions.m = SOG
    mesh_dimensions.n = SOG
    mesh_dimensions.calculate_spacing

    [x,y] = np.meshgrid(
            range(mesh_dimensions.left, mesh_dimensions.right, round(mesh_dimensions.dx)), 
            range(mesh_dimensions.bottom, mesh_dimensions.top, round(mesh_dimensions.dy)))
    
    #### matrix containing the evaluation of the distance based level set function at each of the grid points
    perm_input = usr_input_dict['permeability_input']
    phi_test = get_phi_test(mesh_dimensions, perm_input)
    
    #### Defining the rhs of the eliptical system:
    src_matrix = np.zeros((mesh_dimensions.n + 1, mesh_dimensions.m + 1))
    
    ### Setting Permeability State:
    if(perm_input == 1 or perm_input == 2):
        src_matrix[:,0] = MASS_FLOW_MAGNITUDE
        src_matrix[:,mesh_dimensions.m+1] = -1*MASS_FLOW_MAGNITUDE
    elif(perm_input == 5):
        src_matrix[0,0] = MASS_FLOW_MAGNITUDE
        src_matrix[mesh_dimensions.n+1, mesh_dimensions.m+1] = -1*MASS_FLOW_MAGNITUDE
    else:
        print('Improper permeability input provided')
        exit(1)
    
    initalized_param_dict = {
            "np_mesh_grid" : {"x" : x, "y" : y},
            "dimensions" : mesh_dimensions,
            "phi_initial" : phi_test,
            "source_matrix" : src_matrix,
            "permeability_flag" : perm_input
    }

    return initalized_param_dict


def sim_auto_runs(usr_input_dict, initalized_param_dict):
    """This will contain the code for enabling automatic runs"""
    #### Initialize RP variables
    s0 = 0.79
    c0 = usr_input_dict['c0iter']
    g0 = usr_input_dict['g0iter']
    c0_array = c0 * np.ones((SOG+1, SOG+1))
    [UU, CC, GG] = inital_polymer_surfactant_concentration(
            initalized_param_dict['dimensions'], 
            initalized_param_dict['phi_initial'], 
            s0, c0, g0, 1)
    interface = np.zeros((60,1))
    
    #### Determining the viscosity of surfactant, water, and polymer
    viscosity_oil = 10 #miuo in MATLAB code
    viscosity_water = 1.26 #miuw in MATLAB code
    
    #determining the viscosity of polymer
    beta_1 = 15000
    viscosity_polymer = viscosity_water*(1+beta_1*c0)
    viscosity_polymer_array = viscosity_polymer*np.ones((SOG+1, SOG+1))

    #### Initializing the residual saturations before critical capillary num
    resid_sat_oil = 0.1
    resid_sat_water = 0.3
    
    #### Initializing Flags for simulation
    viscosity_flg = usr_input_dict['model_type']
    polymer_flg = usr_input_dict['polymer_input']
    TLmxixingflg = usr_input_dict['Todd_Longstaff_mixing_parameter']
    shear_flag = usr_input_dict['shear_flag']

    #### Initializing time parameters and matrices for loop
    mesh_grid = initalized_param_dict['dimensions']
    t = TIME_START
    tcal = TIME_START
    tstop = TIME_STOP
    dt = CFL * mesh_grid.dx / MASS_FLOW_MAGNITUDE
    u = np.zeros((mesh_grid.n + 1, mesh_grid.m +1))
    v = u
    coc_matrix = np.zeros((NSIM, 2000))
    prod_rate_matrix = np.zeros((NSIM, math.floor(tstop / dt)))
    croip_matrix = np.zeros((NSIM, math.floor(tstop/dt)))
    mu_save = np.zeros()
    shear_save = np.zeros()
    mass_flow_total = 0;
    sum_UU = 0
    
    while(t < tstop and UU[mesh_grid.n + 1,mesh_grid.m + 1] <= 0.70):
        #Updating the source total:
        mass_flow_total += MASS_FLOW_MAGNITUDE
        t += dt
        inner_iterator = 0
        epsilon = 10
        
        #solving elliptic equations
        alpha = usr_input_dict['alpha_value']
        sigma = 10.001 / ((alpha * GG) + 1) - 0.001
        
        #Determining the aqueous soln visocity as a function of polymer
        compvis_params = {
                "model_type" : usr_input_dict['model_type'],
                "polymer_type" : usr_input_dict['polymer_input'],
                "beta_1" : beta_1,
                "c0" : c0,
                "c0_array" : c0_array
            }
        viscosity_dict = {
                "water" : viscosity_water,
                "oil" : viscosity_oil,
                "polymer" : viscosity_polymer,
                "polymer_array" : viscosity_polymer_array
                }
        [viscosity_aqu, shear] = compvis(CC, u, v, initalized_param_dict['np_mesh_grid']['x'], initalized_param_dict['np_mesh_grid']['y'], compvis_params, viscosity_dict)
        


    

def sim_visualization():
    
    pass

if(__name__ == "__main__"):
    main()

