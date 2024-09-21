import numpy as np

# function to compute viscosity of injected
# displacing phase containing polymer
# miuo = Displaced phase viscosity, c = polymer concentration
#

### INITIALIZING GLOBAL PARAMETERS HERE!
from .. import master_surf_grid 
from . import Exceptions

def compvis(c, U, V, X, Y, params, viscosity):
    ### Initializing variables:
    gamma_dot = np.zeros_like(c)
    viscosityFlag = params['model_type']
    vis_water = viscosity['water']
    vis_oil = viscosity['oil']
    vis_polymer = viscosity['polymer']
    polymer_type = params['polymer_type']
    beta1 = params['beta1']
    c0 = params['c0']
    c0_array = params['c0_array']

    if ( viscosityFlag == 1 ):
        # Newtonian Model (NO SHEAR THINNING INVOLVED => MODEL TYPE #1)
        n = np.shape(( c,1 ))
        m = np.shape(( c,2 ))
        if c0 == 0:
            vis_aqueous = vis_water * np.ones((n, m))
        else:
            vis_aqueous = vis_water * (1 + beta1 * c)
    elif ( viscosityFlag == 2 ):
        # Sourav's Implementation (MODEL TYPE #2)
        n = np.shape((c,1))
        if(c0 == 0):
            vis_aqueous = vis_water*np.ones(n)
        else:
            vis_aqueous = vis_oil*(0.5+c)
    elif( viscosityFlag == 3 ):
        # Dynamic Viscosity (SHEAR THINNING ON => MODEL TYPE #3)
        rho_water = 1000 # kg/m^3
        rho_xanthane = 1500 #kg/m^3
        rho_schizophyllan = 1300 #kg/m^3
        
        if(polymer_type == 0):
            #Xanthane polymer
            w1 = rho_xanthane * c
            w10 = rho_xanthane * c0_array
        elif(polymer_type == 1):
            #Schizophyllan polymer
            w1 = rho_schizophyllan * c
            w10 = rho_schizophyllan * c0_array
        else:
            print("<will raise exception here>")
            exit(1)
        
        w2 = rho_water * (1-c)
        w20 = rho_water * (1- c0_array)

        wppm0 = (w10 / (w10 + w20)) * 1*10^(6)

        #determining the epsilon and n values for the power law equation:

        



    pass

