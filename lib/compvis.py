import numpy as np

# function to compute viscosity of injected
# displacing phase containing polymer
# miuo = Displaced phase viscosity, c = polymer concentration
#

# from .. import master_surf_grid 
# need to also import exceptions here!

def compvis(c, U, V, X, Y, params, viscosity):
    ### Initializing variables:
    gamma_dot = np.zeros_like(c)
    viscosityFlag = params['model_type']
    vis_water = viscosity['water']
    vis_oil = viscosity['oil']
    vis_polymer = viscosity['polymer']
    vis_polymer_array = viscosity['polymer_array']
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
        wppm = (w1 / (w1 + w2))* (10**6)
        #determining the epsilon and n values for the power law equation:
        epsilon_coefficients = {
            "xanthane" : [1.15410398e-04, 2.04937780e+00],
            'schizophyllan' : [0.03647214, 1.32175949]
        }
        n_coefficients = {
            "xanthane" : [3.05428284, -0.27294817],
            "schizophyllan" : [ 4.86265534, -0.41570227]
        }

        if(polymer_type == 0):
            e_coeff = epsilon_coefficients['xanthane']
            n_coeff = n_coefficients['xanthane']
        else:
            e_coeff = epsilon_coefficients['schizophyllan']
            n_coeff = n_coefficients['schizophyllan']

            
        e_power_value = pow((e_coeff[0]*wppm0), e_coeff[1])
        n_power_value = min(pow((n_coeff[0]*wppm0), n_coeff[1]))
        e_vector = pow((e_coeff[0]*wppm), e_coeff[1])
        n_vector = min(pow((n_coeff[0]*wppm0), n_coeff[1]))
        
        n = np.shape((c,1))
        m = np.shape((c,2))

        vis_aqueous = vis_water * np.ones((n,m))
        
        dList = []
        dList[0] = divergence(X,V)
        dList[1] = divergence(Y,U)
        dList[2] = divergence(X,U)
        dList[3] = divergence(Y,V)

        pi_D = np.abs(-0.25* ( (dList[0] + dList[1])**2 ) + (dList[2] * dList[3])) 
        
        #Updating the polymer viscosity matrix
        for ii in range(n):
            for jj in range(m):
                if c[ii, jj] > 0:
                    gamma_dot[ii, jj] = 2 * np.sqrt(pi_D[ii, jj])
                    if gamma_dot[ii, jj] != 0:
                        vis_polymer_array[ii, jj] = e_power_value[ii, jj] * (gamma_dot[ii, jj] ** (n_power_value[ii, jj] - 1))
                        vis_aqueous[ii, jj] = e_vector[ii, jj] * (gamma_dot[ii, jj] ** (n_vector[ii, jj] - 1))

                        # Applying constraints
                        if vis_aqueous[ii, jj] < vis_water:
                            vis_aqueous[ii, jj] = vis_water
                        if vis_aqueous[ii, jj] > 100:
                            vis_aqueous[ii, jj] = 100
                        if vis_polymer_array[ii, jj] < vis_water:
                            vis_polymer_array[ii, jj] = vis_water
                        if vis_polymer_array[ii, jj] > 100:
                            vis_polymer_array[ii, jj] = 100


def divergence(F1, F2):
    """
    Calculate the divergence of a 2D vector field.
    
    Parameters:
    F1, F2 : 2D numpy arrays
        Components of the vector field
    
    Returns:
    div : 2D numpy array
        Divergence of the vector field
    """
    return np.gradient(F1, axis=1) + np.gradient(F2, axis=0)



