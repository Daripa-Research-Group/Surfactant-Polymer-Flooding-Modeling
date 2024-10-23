'''
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40

This code was derived from EOR repository developed by Sourav Dutta and Prabir Daripa

@author: Bhargav Akula Ramesh Kumar
'''

#### IMPORT STATEMENTS
import numpy as np
import tkinter as tk
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
from user_input.json_input_parser import get_user_input
from user_input.gui import UserInputGUI


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

def sim_condition_initialization(usr_input_dict) -> dict:
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
        
        
def main() -> None:
    root = tk.Tk()
    app = UserInputGUI(root)
    root.mainloop()
    
    user_input = app.get_input()
    for simulation in user_input:
        sim_condition_initialization(simulation)



if(__name__ == "__main__"):
    main()

