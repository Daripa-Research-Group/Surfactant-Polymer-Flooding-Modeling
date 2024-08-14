import numpy as np

def KKdef(counter, sizeofgrid, x, y, flag):
    global KK
    
    
    if flag == 1:
        ###########################################
        # Defining a homogeneous permeability field
        ###########################################
        #
        Kmax = 100 # JCP 2017 Daripa and Dutta used 1000
    elif flag == 2:
        ###############################################################
        # Defining a continuous heterogeneous permeability function
        ###############################################################
 #        KK = 1000*sin(4*x) + 1000*cos(5*y) + 100*exp((x-.5)^2+(y-.5)^2);        
        Kmax = 100 # JCP 2017 Daripa and Dutta used 1000
        KK = Kmax * (0.5 * (1 - 10**(-7)) * (np.sin(6 * np.pi * np.cos(x)) * np.cos(4 * np.pi * np.sin(3 * y)) - 1) + 1)
    elif flag == 3:
        #############################################################
        # Defining a permeability function with an impermeable block
        #############################################################        
        bn = sizeofgrid[counter] + 1
        KK = 3000 * np.ones((bn, bn))
        KK[bn // 2 - bn // 8:bn // 2 + bn // 8, bn // 2 - bn // 8:bn // 2 + bn // 8] = 3
    elif flag == 4:
        #############################################################
        # Defining a permeability function with an impermeable block
        ############################################################# 
        bn = sizeofgrid[counter] + 1
        KK = 3000 * np.ones((bn, bn))
        KK[3 * bn // 4 - bn // 12:3 * bn // 4 + bn // 12, 2 * bn // 3 - bn // 12:2 * bn // 3 + bn // 12] = 3
        KK[bn // 3 - bn // 10:bn // 3 + bn // 10, bn // 3 - bn // 10:bn // 3 + bn // 10] = 3
    elif flag == 5:
        ################################################################
        # Simulating with Upper Ness permeability from SPE10
        ################################################################        
        KK = np.load('KK30Ness.npy')
        print('Upper Ness formation permeability loaded')
    elif flag == 6:
        ################################################################
        # Simulating with Tabert permeability from SPE10
        ################################################################
        KK = np.load('KK30Tabert.npy')
        print('Tarbert formation permeability loaded')
    
#     ############ Plotting the permeability field ##################
#     figure(100)
#     surf(x,y,KK);
#     export_fig(sprintf('KK#dx#d.pdf',para.box.m+1, para.box.n+1),'-opengl');
#     
#     pause
#     #---------------------------------------------------------------



#     #---------------------------------------------------------------
#     ########################################################
#     # Defining a random matrix for heterogeneous coefficient
#     ########################################################
#     #
#     if counter == 1
#         load('KK16.mat');
#     elseif counter == 2
#         load('KK16.mat');
#         KK1 = interp2(KK,1); KK = KK1; clear KK1;
#     else 
#         load('KK16.mat')
#         KK2 = interp2(KK,2); KK = KK2; clear KK2;
#     end
#     #---------------------------------------------------------------

    

    