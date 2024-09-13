import numpy as np

# function to compute viscosity of injected
# displacing phase containing polymer
# miuo = Displaced phase viscosity, c = polymer concentration
#

def compvis(c, U, V, X, Y):
    global viscosityFlag
    global miuw, miuo, beta1, c0, miup_array, c0_array, polymerType
    
    gamma_dot = np.zeros_like(c)
    n, m = c.shape
    
    if viscosityFlag == 1:
        if c0 == 0:
            miua = miuw * np.ones((n, m))
        else:
            miua = miuw * (1 + beta1 * c)  #miuo*(0.5+c);#
    elif viscosityFlag == 2:
        # Sourav's implementation
        if c0 == 0:
            miua = miuw * np.ones(n)
        else:
            miua = miuo * (0.5 + c)
    elif viscosityFlag == 3:
        # Dynamic Viscosity Model (Non-Newtonian Model)
        rho_water = 1000
        rho_xanthane = 1500
        
        w1 = rho_xanthane * c
        w2 = rho_water * (1 - c)
        wppm = (w1 / (w1 + w2)) * 1e6

        w10 = rho_xanthane * c0_array
        w20 = rho_water * (1 - c0_array)
        wppm0 = (w10 / (w10 + w20)) * 1e6     
#         [ 4.86265534 -0.41570227] (Schizo – n)
#         [0.03647214 1.32175949] (Schizo – epsilon)
#                 
#         [ 3.05428284 -0.27294817] (Xanthane – n)
#         [1.15410398e-04 2.04937780e+00] (Xanthane – epsilon)   

    if polymerType:
        epsCoeff = [0.03647214, 1.32175949]
        nCoeff = [4.86265534, -0.41570227]
    else:
        epsCoeff = [1.15410398e-04, 2.04937780e+00]
        nCoeff = [3.05428284, -0.27294817]

    epsilon0 = epsCoeff[0] * np.power(wppm0, epsCoeff[1])
    power_n0 = np.minimum(nCoeff[0] * np.power(wppm0, nCoeff[1]), 1)
    epsilon = epsCoeff[0] * np.power(wppm, epsCoeff[1])
    power_n = np.minimum(nCoeff[0] * np.power(wppm, nCoeff[1]), 1)

#         epsilon0(:) = 250.*ones;
#         power_n0(:) = 0.5.*ones;
#         epsilon(:) = 250.*ones;
#         power_n(:) = 0.5.*ones;

    miua = miuw * np.ones((n, m))
    a1 = np.sum(np.gradient(X, axis=(0, 1), axis=(1, 2)))
    a2 = np.sum(np.gradient(Y, axis=(0, 1), axis=(1, 2)))
    a3 = np.sum(np.gradient(X, axis=(1, 2), axis=(0, 1)))
    a4 = np.sum(np.gradient(Y, axis=(1, 2), axis=(0, 1)))
    pi_D = np.abs(-0.25 * ((a1 + a2) ** 2) + a3 * a4)
            
    for ii in range(n):
        for jj in range(m):
            if c[ii, jj] > 0:
                gamma_dot[ii, jj] = 2 * np.sqrt(pi_D[ii, jj])
                if gamma_dot[ii, jj] == 0:
                    pass
                else:
                    miup_array[ii, jj] = epsilon0[ii, jj] * (gamma_dot[ii, jj] ** (power_n0[ii, jj] - 1))
                    miua[ii, jj] = epsilon[ii, jj] * (gamma_dot[ii, jj] ** (power_n[ii, jj] - 1))
                    if miua[ii, jj] < miuw:
                        miua[ii, jj] = miuw
                    if miua[ii, jj] > 100:
                        miua[ii, jj] = 100
                    if miup_array[ii, jj] < miuw:
                        miup_array[ii, jj] = miuw
                    if miup_array[ii, jj] > 100:
                        miup_array[ii, jj] = 100
                        
    gammaMax = np.max(gamma_dot)
