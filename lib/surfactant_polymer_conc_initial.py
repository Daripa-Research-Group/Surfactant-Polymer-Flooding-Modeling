import numpy as np

# function to initialize s,c,g in the domain
# flag = 0 is no surfactantimplementation#
# flag = 1 is with surfactant
# 1-s0 = initial residual saturation 
# c0 = concentration of polymer in injected mixture
# g0 = concentration of surfactant in injected mixture

### Vectorized implementation
def inital_polymer_surfactant_concentration(para_box, phi, s_0, c_0, g_0, flag):
    s0 = np.zeros((para_box["n"] + 1, para_box["m"] + 1))
    c0 = np.copy(s0)
    D = (phi > 1e-10) + (np.abs(phi) < 1e-10)
    if flag == 0:
        g0 = []
        s0 = (~D) + D * (1 - s_0)
        c0 = (~D) * c_0
    elif flag == 1:
        s0 = (~D) + D * (1 - s_0)
        c0 = (~D) * c_0
        g0 = (~D) * g_0
        
    return s0, c0, g0
        
