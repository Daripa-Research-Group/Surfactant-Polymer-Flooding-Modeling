'''
Determine which subdomain the grid points belong to
- Evaluate the level set function on the grid points
- Negative output signifies the domain $$\Omega^-$$.
- Positive output signifies the domain $$\Omega^+$$.

This code was derived from MATLAB

@author: Carlos Acosta Carpio
'''

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
import numpy as np
from z_func_test import z_func_test
from src.para import Box

def get_phi_test(para_box : Box, permeability_input : int):
    m = para_box.m
    n = para_box.n
    
    dx = para_box.dx
    dy = para_box.dy
    
    left = para_box.left
    bottom = para_box.bottom
    
    # ------Loop Implementation
    # get_phi=zeros(n+1,m+1);
    # for ii=1:m+1
    #     for jj=1:n+1
    # 
    #             get_phi(jj,ii)=z_func_test(left+(ii-1)*dx,bottom+(jj-1)*dy);
    # 
    #     end
    # end
    #-------------------
    
    # --- Vectorized implementation
    jj, ii = np.meshgrid(np.arange(1, n + 2), np.arange(1, m + 2))
    get_phi = z_func_test(left+(ii - 1) * dx, bottom+(jj - 1)*dy, permeability_input)
    
    return get_phi
