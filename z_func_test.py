import numpy as np

## Specifying the initial position of the water front
# A function describing the initial position of 
# the water front in the shape of a circular arc 
# $$ z(x,y) = x^2+y^2-0.015 $$ 
# This can take array input
## 

def z_func_test(x, y, permeability_input):
    init_front_hs = 0.1
    # out=(x-0.15)^2+(y-0.15)^2 -0.5*(x+y-0.35)^2;



    # # perturbed initial saturation front for special fingering simulations
    # out = x.^2 + y.^2 - 0.015*(1+0.1*sin(18*atan(y./x)))^2;

    #setting default value to out
    out = 0
    
    #homogenous
    if permeability_input == 1:
        out = y - init_front_hs + 0.01 * (np.cos(80 * np.pi * x))
    elif permeability_input == 2:
        out = y - init_front_hs ## Rectilinear Homogenous
    elif permeability_input == 5:
        out=(x)**2+(y)**2-0.015; # Normal unperturbed initial saturation front 
    
    return out