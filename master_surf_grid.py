'''
This is the main program for a full flooding simulationdisp
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40
'''
import numpy as np
import matplotlib.pyplot as plt
from get_phi_test import get_phi_test
from KKdef import KKdef
from s0c0 import s0c0
from compvis import compvis

# clearing variables and setting all initial paramters
plt.close("all")
np.set_printoptions(precision=15, suppress=True) #adjusts print setting for improved readability

# declaring global variables to reduce computational time in passing through function calls
global c0_array, miuo, miuw, miup, swr0, sor0, dt, KK, s0, c0, g0, beta1, viscosityFlag, shearFlag, miup_array, polymerType

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
model_type = input("Enter which computation model to use for this simulation: ")
print("----------------------------------------------")
figure_num = input("enter figure number being simulated here: ")
c0iter_val = input("Enter value for c0iter(0.001 as default): ")
g0iter_val = input("Enter value for g0iter_val: ")
print("----------------------------------------------")
print("0 = Saturation(UU) Plot")
print("1 = Polymer(CC) Plot")
print("2 = Surfactant(GG) Plot")
plot_type_input = input("Enter the type of plot desired: ")
print("----------------------------------------------")
print("0 = xanthane")
print("1 = schizophyllan")
polymer_input = input("Enter the desired polymer for the simulation: ")
print("----------------------------------------------")
print("1 = homogeneous with magnitude 1000")
print("2 = continuous heterogeneous function")
print("3 = impermeable block inclusion at the center")
print("4 = impermeable block inclusion off the center")
print("5 = Upper Ness from SPE10 model sections")
print("6 = Tabert from SPE10 model sections")
permeability_input = input("Enter the corresponding permeability flag for the simulation: ")

'''
#### Initialization    ########################################
### nsim = number of runs
### sizeofgrid = mesh size for each run
### c0iter = injection concentration of polymer for each run
### g0iter = injection concentration of surfactant for each run
###############################################################
'''
nsim = 1
sog = 29
sizeofgrid = [sog, sog, sog, sog]
c0iter = [c0iter_val, 0, 0, 0.001] # 0.0002 (300 wppm) 0.0006 (900 wpppm) 0.001 (1500 wpppm)
g0iter = [g0iter_val, 0, 0, 0] # [0.01 0.01 0.01]
k = 1

'''
### preallocating variables that will store cpu time taken for 
### each run and number of iterations/time-steps of each run
'''
time = [0, 0, 0]
iterations = [0, 0 , 0]


#running through simulations
for counter in range(1, nsim + 1):
    start_time = time.time() # starting timing for the simulation
    
    
    #setting x and y grid points in the para.box dictionary
    para_box["m"] = sizeofgrid[counter - 1]
    para_box["n"] = sizeofgrid[counter - 1]
    #calculate dx and dy based on the updated m and n values
    para_box["dx"] = (para_box["right"] - para_box["left"]) / para_box["m"]
    para_box["dy"] = (para_box["top"] - para_box["bottom"]) / para_box["n"]
    #creating grid points and grids with meshgrid function
    x, y = np.meshgrid(np.arange(para_box['left'], para_box['right'] + para_box['dx'], para_box['dx']),
                       np.arange(para_box['bottom'], para_box['top'] + para_box['dy'], para_box['dy']))
    
    #Matrix containing the evaluation of the distance based level set function at each of the grid points
    phi_test = get_phi_test(para_box, permeability_input)

    #Defining the right hand of elliptic system - source terms
    f = np.zeros(para_box["n"] + 1, para_box["m"] + 1)
    MFW = 0
    interX_save = 0
    
    src = 120000; #50000; #5000; #2; #200 for QFS    # magnitude of mass flow rate at source (1.2 for heterogeneous)1.2 to 30 for shear cases
    #rate = 0
    
    ## setting permeability state
    if permeability_input == 1 or permeability_input == 2:
        f[ : 0] = src # intensity of injection well = src
        f[ : para_box["m"] + 1] = -src # intensity of production well = -src
    elif permeability_input == 5:
        f[0:0] = src # intensity of injection well = src
        f[para_box["n"] + 1 : para_box["m"] + 1] = -src # intensity of production well = -src
        
    # defining permeability matrix
    permeabilityFlag = permeability_input
    KKdef(counter, sizeofgrid, x, y, permeabilityFlag)

    #     KK = 100*KK;
    #-------------------------------------------------------------------
    
    
    ########### Preparing automatic runs  ###########################
    
    # initialize RP variables----------------------------------------------
    # UU=zeros(para.box.n+1,para.box.m+1);
    s0 = 0.79;  # initial residual water saturation inside the reservoir = 1-s0
    c0 = c0iter[counter - 1] # 0.1
    g0 = g0iter[counter - 1] # 0.005
    c0_array = c0 * np.ones((sog + 1, sog + 1))
    UU, CC, GG = s0c0(para_box, phi_test, s0, c0, g0, 1)
    interface = np.zeros(60)
    ###
    miuo=10 #2e-3;#12.6;  ### 0.95
    miuw=1.26 #8.9e-4;#1.26;#1.26; #0.0430;   ### 0.095   # Effective aq. viscosity = miuo*(0.5+c)
    beta1= 15000 #2;
    miup = miuw * (1 + beta1 * c0) ## Assumption
    miup_array = miup * np.ones((sog+1,sog+1))
    #alpha=0.6;  ### 0.6
    # define initial residual saturations before critical capillary number 
    swr0 = 0.1 # actual initial residual saturation is higher than this value in practice
    sor0 = 0.3
    
    viscosityFlag = model_type #Flag to enable different viscosity models
        # 3 = Dynamic Viscosity Model (Non-Newtonian model) SHEAR THINNING
        # ON
        # 2 = Sourav's Model
        # 1 = Original model from JCP paper SHEAR THINNING OFF

    polymerType = polymer_input; #Flag to initialize with a particular polymer
        # 0 = Xanthane
        # 1 = Schizophyllan

    TLmixingFlag = 0; #Flag to enable Todd-Longstaff mixing parameter
    # 1 = TL mixing model ON
    # 0 = TL mixing model OFF
   
    shearFlag = 0; #Flag to enable shear effects with ad-hoc model from Eclipse code
    # 1 = Shear effects ON
    # 0 = Shear effects OFF
    
    t = 0
    # cfl = 0.2
    u = np.zeros((para_box["n"] + 1, para_box["m"] + 1))
    v = np.zeros((para_box["n"] + 1, para_box["m"] + 1))
    
    #     u_old = zeros(para.box.n+1,para.box.m+1); v_old=u_old; u=u_old; v=v_old;
    CFL = 1
    dt = CFL * para_box["dx"] / src
    if permeability_input == 5: # 100 only for quarter five spot
        dt = CFL * para_box["dx"] / src * 100
    #dt = 1/25
    tstop = 500
    COC = np.zeros((nsim, 2000)) # floor(tstop/dt)
    ProdRate = np.zeros((nsim, int(tstop / dt)))
    CROIP = np.zeros((nsim, int(tstop / dt)))
    miuaSave = np.zeros((0,0))
    shearSave = np.zeros((0,0))
    concSave = np.zeros((0,0))
    src_total = 0
    sumUU = 0
    
    while t < tstop and UU[para_box["n"] + 1, para_box["m"] + 1] <= 0.70:
      #      src = src + rate*tcal*0.22941; #200 for QFS    # magnitude of mass flow rate at source (1.2 for heterogeneous)1.2 to 30 for shear cases
        src_total = src_total + src
        
        t += dt
        innerIter = 0
        epsilon = 10
        #Solve elliptic equation.................................
#        while(epsilon > 1e-4) 
        # IFT as a function of surf concentration
        # $$ \sigma = \frac{10.001}{\Gamma +1} - 0.001 $$
        sigma = 10.001 / (GG + 1) - 0.001 #THIS IS SURFACE TENSION
        
        # aq soln viscosity as a function of polymer $$ \mu_a =
        # \mu_o(0.5+c) $$
        miua, shear = compvis(CC, u, v, x, y)  
        if viscosityFlag == 1:
            miup = np.max(miua[0 : ])
        miup_array = miup * np.ones((sog + 1, sog + 1))
        if TLmixingFlag:
            u_w_eff, m_mu = TLmixing(miua, CC)
            miua = mu_w_eff           
        
