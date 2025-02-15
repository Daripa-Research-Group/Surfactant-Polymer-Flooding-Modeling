"""
This python script contains the class definition for running simulations
"""

import sys
import os
import scipy as sp
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

#Relevant imports
from lib.Exceptions import SimulationCalcInputException
from lib.para import Box
from lib.enumerations import SimulationConstants, PolymerList, ModelType, ResevoirGeometry, PermeabilityType


class Simulation:
    """
    Class is used to generate an instance of a simulation which will run the SP-flooding model based on the given parameters provided by the user
    """
    def __init__(self, sim_id, size_of_grid, polymer, surfactant, init_water_saturation, resevoir_geometry, permeability_flg, mesh_grid, mdl_id, plt_type):
        """
        creates instance of the simulation class which will enable for calculating changes in system parameters at every time-step

        :param sim_id: Simulation number
        :type sim_id: int

        :param size_of_grid: size of mesh
        :type size_of_grid: int

        :param polymer: Polymer object used in SP-flooding run
        :type polymer: Polymer

        :param surfactant: Surfactant object used in SP-flooding run (can also be non if surfactant concentration = 0)
        :type surfactant: Surfactant, None

        :param init_water_saturation: Initial Water Saturation (scalar quantitiy)
        :type init_water_saturation: float

        :param resevoir_geometry: Type of resevoir geometry (is it a rectilinear or quarter-five-spot geometry)
        :type resevoir_geometry: enum 'ResevoirGeometry'

        :param permeability_flg: Homogenous vs. Heterogenous porosity in resevoir
        :type permeability_flg: enum 'PermeabilityType'

        :param mesh_grid: mesh_grid used in the SP-flooding run
        :type mesh_grid: Box

        :param mdl_id: the model id for the simulation (whether shear thinning is on or off)
        :type mdl_id: enum 'ModelType'

        :param plt_type: the plot type outputted by the program for the simulation run
        :type plt_type: enum 'PlotType'
        """
        
        self.sim_id = sim_id
        self.sog = size_of_grid

        #Instances of the polymer and surfactant objects
        self.polymer = polymer
        self.surfactant = surfactant

        #Simulaiton Properties
        self.resevoir_geometry = resevoir_geometry
        self.permeability_flg = permeability_flg
        self.mesh = mesh_grid
        self.init_water_saturation_scalar = init_water_saturation
        
        #Model and plotting flags
        self.mdl_id = mdl_id #Model type
        self.plt_type = plt_type #types of plots to generate

    ### DEPENDENT VARIABLES OF SIMULATION CLASS:
    _phi_ = None
    @property
    def phi(self):
        if(self._phi_ is None):
            self._phi_ = self.get_phi_value()
        return self._phi_

    _water_saturation_vector_form_ = None
    @property
    def water_saturation(self): #vector version
        return self._water_saturation_vector_form_

    @water_saturation.setter
    def water_saturation(self, value):
        self._water_saturation_vector_form_ = value

    _aqueous_viscosity_ = None
    @property
    def aqueous_viscosity(self):
        return self._aqueous_viscosity_
    
    @aqueous_viscosity.setter
    def aqueous_viscosity(self, value):
        self._aqueous_viscosity_ = value

    _oleic_mobility_ = None
    @property
    def oleic_mobility(self):
        return self._oleic_mobility_

    @oleic_mobility.setter
    def oleic_mobility(self, value):
        self._oleic_mobility_ = value

    _aqueous_mobility_ = None
    @property
    def aqueous_mobility(self):
        return self._aqueous_mobility_

    @aqueous_mobility.setter
    def aqueous_mobility(self, value):
        self._aqueous_mobility_ = value

    @property
    def sigma(self): #TODO: Need to update to make sure that i calculate the concentration using the lambda function within the surfactant object
        # print(type(self.surfactant.vec_concentration))
        return self.surfactant.IFT_conc_equ(np.array(self.surfactant.vec_concentration))

    @property
    def derivative_sigma(self):
        return self.surfactant.derivative_IFT_conc_equ(np.array(self.surfactant.vec_concentration))


    _is_surfactant_ = None
    @property
    def is_surfactant(self):
        """
        This function differentiates whether the model is SP-Flooding or just polymer flooding
        """
        if(self._is_surfactant_ is None):
            if(self.surfactant.concentration == 0):
                self._is_surfactant_ = False
            else:
                self._is_surfactant_ = True
        return self._is_surfactant_


    def get_phi_value(self):
        """
        function to initialize s,c,g in the domain
        flag = 0 is no surfactantimplementation#
        flag = 1 is with surfactant
        1-s0 = initial residual saturation 
        c0 = concentration of polymer in injected mixture
        g0 = concentration of surfactant in injected mixture

        :return: phi value in vector form
        :rtype: np.meshgrid
        """
        try:
            #getting values from mesh:
            if(isinstance(self.mesh, Box)):
                m = self.mesh.m
                n = self.mesh.n
                
                dx = self.mesh.dx
                dy = self.mesh.dy
                
                left = self.mesh.left
                bottom = self.mesh.bottom
            else:
                raise SimulationCalcInputException("SimulationError | get_phi_value | mesh values not provided")
            
            # --- Vectorized implementation
            jj, ii = np.meshgrid(np.arange(1, n + 2), np.arange(1, m + 2))
            phi_vec = self.z_func_test(left + (ii - 1) * dx, bottom + (jj - 1) * dy)
            return phi_vec
        except Exception as e:
            print(e)
            exit(1)


    def z_func_test(self, x, y):
        """
        Specifying the initial position of the water front
        A function describing the initial position of 
        the water front in the shape of a circular arc 
        $$ z(x,y) = x^2+y^2-0.015 $$ 
        This can take array input

        :param x: x dimension value
        :type x: float

        :param y: y dimension value
        :type y: float


        :return: returns water front position
        :rtype: np.array
        """
        init_front_hs = 0.1

        #setting default value to out
        out = 0
        
        #homogenous
        if ( self.permeability_flg == PermeabilityType.Homogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            out = y - init_front_hs + 0.01 * (np.cos(80 * np.pi * x))
        elif ( self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            out = y - init_front_hs ## Rectilinear Homogenous
        elif ( self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Quarter_Five_Spot ):
            out = ( (x)**2 ) + ( (y)**2 ) - 0.015 # Normal unperturbed initial saturation front 
        
        return out


    def initial_concentration_matrix(self):
        """
        function to initialize s,c,g in the domain
        flag = 0 is no surfactantimplementation#
        flag = 1 is with surfactant
        1-s0 = initial residual saturation 

        Vectorized implementation
        
        :return: List object with the vectorized implementation of the water saturation (s0), polymer concentration (c0), and surfactant concentration (g0)
        :rtype: List
        """
        try:
            if(self.polymer is not None and self.surfactant is not None and self.init_water_saturation_scalar is not None):
                # retrieving scalar quantities of concentration for surfactant and polymer
                # scalar quantity of initial water saturation (fraction of pore space filled with water)
                print("reaches here")
                s_0 = self.init_water_saturation_scalar
                c_0 = self.polymer.initial_concentration
                g_0 = self.surfactant.concentration
                self.water_saturation = np.zeros((self.mesh.m + 1, self.mesh.n + 1))
                self.polymer.vec_concentration = np.copy(self.water_saturation)
            else:
                raise SimulationCalcInputException("SimulationError: Did not provide intial water saturation or Polymer and/or Surfactant were not intialized. Please try again...")

            if(self.phi is not None):
                D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
            else:
                raise SimulationCalcInputException("SimulationError: phi value not calculated!")

            self.water_saturation = (~D) + D * (1 - s_0)
            self.polymer.vec_concentration = (~D) * c_0
            self.surfactant.vec_concentration = (~D) * g_0
            
            return [ self.water_saturation, self.polymer.vec_concentration, self.surfactant.vec_concentration ]
        
        except Exception as e:
            print(e)
            exit(1)

    
    def compvis(self, U, V, X, Y, beta1, c0_array):
        """
        function to compute viscosity of injected
        displacing phase containing polymer
        """
        ### Initializing variables:
        try:
            if(self.polymer is not None and self.surfactant is not None):
                gamma_dot = np.zeros_like(self.polymer.vec_concentration)
                print("reaches here")
                vis_water = SimulationConstants.Water_Viscosity.value #viscosity['water']
                vis_oil = SimulationConstants.Oil_Viscosity.value #viscosity['oil']
                polymer_obj = self.polymer
            else:
                raise SimulationCalcInputException("SimulationError: Surfactant and/or Polymer Not Initialized")

            if (self.mdl_id == ModelType.No_Shear_Thinning):
                # Newtonian Model (NO SHEAR THINNING INVOLVED => MODEL TYPE #1)
                print("shape of polymer vec concentration:", np.shape(polymer_obj.vec_concentration))
                n = np.shape(polymer_obj.vec_concentration)[0]
                m = np.shape(polymer_obj.vec_concentration)[1]
                if polymer_obj.vec_concentration.all() == 0:
                    self.aqueous_viscosity = vis_water * np.ones((n, m))
                    print('aqueous viscosity: ', self.aqueous_viscosity)
                else:
                    self.aqueous_viscosity = vis_water * (1 + beta1 * polymer_obj.vec_concentration)
            elif (self.mdl_id == ModelType.Sourav_Implementation):
                # Sourav's Implementation (MODEL TYPE #2)
                n = np.shape((polymer_obj.vec_concentration,1))
                if(polymer_obj.initial_concentration == 0):
                    self.aqueous_viscosity = vis_water * np.ones(n)
                else:
                    self.aqueous_viscosity = vis_oil * (0.5 + polymer_obj.vec_concentration)
            elif(self.mdl_id == ModelType.Shear_Thinning_On):
                # Dynamic Viscosity (SHEAR THINNING ON => MODEL TYPE #3)
                rho_water = SimulationConstants.Water_Density.value # kg/m^3
                rho_xanthane = PolymerList.Xanthane.Density #kg/m^3
                rho_schizophyllan = PolymerList.Schizophyllan.Density #kg/m^3
                
                if(polymer_obj.name == PolymerList.Xanthane.Id):
                    #Xanthane polymer
                    w1 = rho_xanthane * polymer_obj.vec_concentration
                    w10 = rho_xanthane * c0_array
                elif(polymer_obj.name == PolymerList.Schizophyllan.Id):
                    #Schizophyllan polymer
                    w1 = rho_schizophyllan * polymer_obj.vec_concentration
                    w10 = rho_schizophyllan * c0_array 
                else:
                    raise SimulationCalcInputException("SimulationError: Polymer Not Part of 'PolymerTypes' Enumerations List")
                
                w2 = rho_water * (1-polymer_obj.vec_concentration)
                w20 = rho_water * (1- c0_array)

                wppm0 = (w10 / (w10 + w20)) * 1*10^(6)
                wppm = (w1 / (w1 + w2))* (10**6)

                #determining the epsilon and n values for the power law equation:
                e_coeff = polymer_obj.e_coeff
                n_coeff = polymer_obj.n_coeff
                    
                e_power_value = pow((e_coeff[0]*wppm0), e_coeff[1])
                n_power_value = min(pow((n_coeff[0]*wppm0), n_coeff[1]))
                e_vector = pow((e_coeff[0]*wppm), e_coeff[1])
                n_vector = min(pow((n_coeff[0]*wppm0), n_coeff[1]))
                
                n = np.shape((polymer_obj.vec_concentration,1))
                m = np.shape((polymer_obj.vec_concentration,2))

                self.aqueous_viscosity = vis_water * np.ones((n,m))
                
                dList = []
                dList[0] = self.divergence(X,V)
                dList[1] = self.divergence(Y,U)
                dList[2] = self.divergence(X,U)
                dList[3] = self.divergence(Y,V)

                pi_D = np.abs(-0.25* ( (dList[0] + dList[1])**2 ) + (dList[2] * dList[3])) 
                
                #Updating the polymer viscosity matrix
                for ii in range(n):
                    for jj in range(m):
                        if polymer_obj.vec_concentration[ii, jj] > 0:
                            gamma_dot[ii, jj] = 2 * np.sqrt(pi_D[ii, jj])
                            if gamma_dot[ii, jj] != 0:
                                polymer_obj.vec_concentration[ii, jj] = e_power_value[ii, jj] * (gamma_dot[ii, jj] ** (n_power_value[ii, jj] - 1))
                                self.aqueous_viscosity[ii, jj] = e_vector[ii, jj] * (gamma_dot[ii, jj] ** (n_vector[ii, jj] - 1))

                                # Applying constraints
                                if self.aqueous_viscosity[ii, jj] < vis_water:
                                    self.aqueous_viscosity[ii, jj] = vis_water
                                if self.aqueous_viscosity[ii, jj] > 100:
                                    self.aqueous_viscosity[ii, jj] = 100
                                if polymer_obj.vec_concentration[ii, jj] < vis_water:
                                    polymer_obj.vec_concentration[ii, jj] = vis_water
                                if polymer_obj.vec_concentration[ii, jj] > 100:
                                    polymer_obj.vec_concentration[ii, jj] = 100
                
            else:
                raise SimulationCalcInputException("SimulationError: Model Type not part of 'ModelType' Enumerations list")
            return [self.aqueous_viscosity, gamma_dot]
        except Exception as e:
            print(e)
            exit(1)


    def compmob(self, sor, swr, flag):
        """
        function to compute mobility

        :param sor: residual oil saturation at IFT sigma (matrix)
        :type sor: float

        :param swr: residual water saturation at IFT sigma (matrix)
        :type swr: float
        
        :param sor0: residual oil saturation at IFT sigma0 (constant)
        :type sor0: float

        :param swr0: residual water saturation at IFT sigma0 (constant)
        :type swr0: float

        :param flag: denotes which phase 0=oleic and 1=aqueous
        :type flag: int

        :return: mobility of aqueous or oleic phases
        :rtype: np.array
        """
        oil_viscosity = SimulationConstants.Oil_Viscosity.value
        swr0 = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
        sor0 = SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value
        if(self.is_surfactant == 0 and self.water_saturation is not None): # (without surfactant)
            # Normalized saturations of water and oil at IFT sigma0
            nsw0 = (self.water_saturation - swr0) / (1 - swr0)
            nso0 = (self.water_saturation - swr0) / (1 - swr0 - sor0)

            # Corey type relative permeability in the absence of surfactant
            krw0 = nsw0**3.5
            kro0 = ((1 - nso0)**2) * (1 - nso0**1.5)

            if(flag == 0): #calculating mobility in oleic phase
                self.oleic_mobility = kro0/oil_viscosity
                return self.oleic_mobility
            elif(flag == 1): #calculating mobility in aqueous phase
                self.aqueous_mobility = krw0/self.aqueous_viscosity
                return self.aqueous_mobility
            else:
                raise SimulationCalcInputException("SimulationError: Flag unknown. Flag can only be a 0 (for oleic mobility calculation) and 1 (for aqueous mobility calculation). Please try again...")
        elif(self.is_surfactant == 1 and self.water_saturation is not None): # (with surfactant)
            # normalized saturations of water and oil at IFT sigma
            nsw = (self.water_saturation - swr)/(1 - swr)
            nso = (self.water_saturation - swr)/(1 - swr - sor)

            # rel perm in presence of surfactant
            krw = nsw * (2.5 * swr * ((nsw**2) - 1) + 1)
            kro = (1 - nso) * (1 - 5 * sor * nso)

            if(flag == 0): #calculating mobility in oleic phase
                self.oleic_mobility = kro/oil_viscosity
                return self.oleic_mobility
            elif(flag == 1): #calculating mobilitty of the aqueous phase
                self.aqueous_mobility = krw/self.aqueous_viscosity
                return self.aqueous_mobility
            else:
                raise SimulationCalcInputException("SimulationError: Flag unknown. Flag can only be a 0 (for oleic mobility calculation) and 1 (for aqueous mobility calculation). Please try again...")
        else:
            raise SimulationCalcInputException("SimulationError: Need proper statement of is_surfactant (0 for Polymer Flooding Simulation & 1 for Surfactant-Polymer Flooding Simulation). Please try again...")

    def compres(self, u, v):
        """
        function to compute residual saturations as a function of surfactant
        concentration via capillary number variation. Hence it varies with change
        in surf concentration and velocity evolution. Must be recomputed at every
        time step.
        """

        # define critical capillary numbers 
        # ie $$N_c $$ at which $$s_{ro}$$ and $$ s_{ra}$$ begin to decrease
        Nco0 = ( 1.44 )*( 10**(-4) )  # Values from Amafuele Handy 1982
        Nca0 = ( 1.44 )*( 10**(-4) )  # these two do not have to be the same
        
        miuo = SimulationConstants.Oil_Viscosity.value
        swr0 = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
        sor0 = SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value

        # compute capillary number
        nca = np.sqrt(u**2 + v**2) * self.aqueous_viscosity / self.sigma
        nco = np.sqrt(u**2 + v**2) * miuo / self.sigma
        Nca = np.linalg.norm(nca)  # compute 2-norm (largest singular value)
        Nco = np.linalg.norm(nco)

        # define residual saturations as functions of capillary numbers
        if Nco < Nco0:
            sor = sor0
        else:
            sor = sor0 * (Nco0 / Nco)**0.5213

        if Nca < Nca0:
            swr = swr0
        else:
            swr = swr0 * (Nca0 / Nca)**0.1534

        return [swr, sor]

    def setTri(self):
        """
        Setting up triangulations for the FEM grid:

        U = cell array with each element = array of vertices of Upper Triangle of the rectangular cell 
        L = cell array with each element = array of vertices of Lower Triangle of the rectangular cell
        At every point (i,j), U{i,j} & L{i,j} are cells with coordinates of vertices of the two triangles 
        obtained by bisecting the rectangle starting at (i,j). The bisection line goes from NW to SE.
         
        """
        U = np.empty((self.mesh.m, self.mesh.n), dtype=object)
        L = np.empty((self.mesh.m, self.mesh.n), dtype=object)

        for j in range(0,self.mesh.m):
            for k in range(0, self.mesh.n):
                x1 = self.mesh.left + j * self.mesh.dx
                y1 = self.mesh.bottom + k * self.mesh.dy
                x2 = self.mesh.left + (j + 1) * self.mesh.dx
                y2 = y1
                x3 = x1
                y3 = self.mesh.bottom + (k + 1) * self.mesh.dy
                x4 = x2
                y4 = y3

                l = {
                    'x': [x1, x2, x3],
                    'y': [y1, y2, y3]
                }

                u = {
                    'x': [x4, x3, x2],
                    'y': [y4, y3, y2]
                }

                U[j, k] = u
                L[j, k] = l
        return [U, L]


    def set_grid(self,U, L, beta):
        m = self.mesh.m
        n = self.mesh.n

        out = [[None for _ in range(n+1)] for _ in range(m+1)]
        

        for j in range(m+1):
            for l in range(n+1):
                if j == 0 and l != 0 and l != n:
                    t1 = self.weak(L[j][l],[1, 0, 0], beta)
                    t2 = [0, 0, 0, 0]
                    t3 = [0, 0, 0, 0]
                    t4 = [0, 0, 0, 0]
                    t5 = self.weak(L[j][l-1],[0, 0, 1], beta)
                    t6 = self.weak(U[j][l-1], [0, 1, 0], beta)
                elif j == m and l != 0 and l != n:
                    t1 = [0, 0, 0, 0]
                    t2 = self.weak(U[j-1][l],[0, 0, 1], beta)
                    t3 = self.weak(L[j-1][l], [0, 1, 0], beta)
                    t4 = self.weak(U[j-1][l-1], [1, 0, 0], beta)
                    t5 = [0, 0, 0, 0]
                    t6 = [0, 0, 0, 0]
                elif j != 0 and j != m and l == 0:
                    t1 = self.weak(L[j][l],[1, 0, 0], beta)
                    t2 = self.weak(U[j-1][l], [0, 0, 1], beta)
                    t3 = self.weak(L[j-1][l], [0, 1, 0], beta)
                    t4 = [0, 0, 0, 0]
                    t5 = [0, 0, 0, 0]
                    t6 = [0, 0, 0, 0]
                elif j != 0 and j != m and l == n:
                    t1 = [0, 0, 0, 0]
                    t2 = [0, 0, 0, 0]
                    t3 = [0, 0, 0, 0]
                    t4 = self.weak(U[j-1][l-1],[1, 0, 0], beta)
                    t5 = self.weak(L[j][l-1],[0, 0, 1], beta)
                    t6 = self.weak(U[j][l-1], [0, 1, 0], beta)
                elif j == 0 and l == 0:
                    t1 = self.weak(L[j][l],[1, 0, 0], beta)
                    t2 = [0, 0, 0, 0]
                    t3 = [0, 0, 0, 0]
                    t4 = [0, 0, 0, 0]
                    t5 = [0, 0, 0, 0]
                    t6 = [0, 0, 0, 0]
                elif j == 0 and l == n:
                    t1 = [0, 0, 0, 0]
                    t2 = [0, 0, 0, 0]
                    t3 = [0, 0, 0, 0]
                    t4 = [0, 0, 0, 0]
                    t5 = self.weak(L[j][l-1], [0, 0, 1], beta)
                    t6 = self.weak(U[j][l-1], [0, 1, 0], beta)
                elif j == m and l == 0:
                    t1 = [0, 0, 0, 0]
                    t2 = self.weak(U[j-1][l], [0, 0, 1], beta)
                    t3 = self.weak(L[j-1][l], [0, 1, 0], beta)
                    t4 = [0, 0, 0, 0]
                    t5 = [0, 0, 0, 0]
                    t6 = [0, 0, 0, 0]
                elif j == m and l == n:
                    t1 = [0, 0, 0, 0]
                    t2 = [0, 0, 0, 0]
                    t3 = [0, 0, 0, 0]
                    t4 = self.weak(U[j-1][l-1], [1, 0, 0], beta)
                    t5 = [0, 0, 0, 0]
                    t6 = [0, 0, 0, 0]
                else:
                    t1 = self.weak(L[j][l], [1, 0, 0], beta)
                    t2 = self.weak(U[j-1][l], [0, 0, 1], beta)
                    t3 = self.weak(L[j-1][l], [0, 1, 0], beta)
                    t4 = self.weak(U[j-1][l-1], [1, 0, 0], beta)
                    t5 = self.weak(L[j][l-1], [0, 0, 1], beta)
                    t6 = self.weak(U[j][l-1], [0, 1, 0], beta)

                grid = {
                    'c': t1[0] + t2[2] + t3[1] + t4[0] + t5[2] + t6[1],
                    'w': t3[0] + t4[1],
                    's': t4[2] + t5[0],
                    'n': t1[2] + t2[0],
                    'e': t1[1] + t6[0],
                    'nw': t2[1] + t3[2],
                    'se': t5[1] + t6[2],
                    'const': t1[3] + t2[3] + t3[3] + t4[3] + t5[3] + t6[3]
                }
                
                out[j][l] = grid

        return out

    def weak(self, T, v, beta):
        b1 = self.beta_func(T.x[0],T.y[0], beta)
        b2 = self.beta_func(T.x[1],T.y[1], beta)
        b3 = self.beta_func(T.x[2],T.y[2], beta)
        b_avg = (b1 + b2 + b3) / 3
        
        s = self.polyarea(T.x, T.y)

        # Create and manipulate matrix M
        M = np.vstack((T.x, T.y, [1, 1, 1])).T
        M_inv = np.linalg.inv(M)
        M = M_inv[:2, :]  # Extract the first two rows of M_inv

        # Calculate vdiff and inte
        vdiff = np.dot(M, v)
        inte = np.dot(vdiff.T, beta * np.dot(M, s))

        # Output result
        out = [inte, 0] 

        return out

    def beta_func(self, x, y, beta):
        dx = self.mesh.dx
        dy = self.mesh.dy
        
        left = self.mesh.left
        bottom = self.mesh.bottom
        
        nn = np.round(( x - left )/dx) + 1
        mm = np.round(( y - bottom )/dy) + 1
        
        out = beta[nn][mm]
        return out

    def setRightHand(self, src_matrix, U, L):
        m = self.mesh.m
        n = self.mesh.n
        
        rh = np.zeros((m+1) * (n+1))
        
        for j in range(1, m + 2):
            for l in range(1, n + 2):
                id = j + (l - 1) * (m + 1) - 1  # Adjust for Python 0-based indexing
                
                if j == 1 and l != 1 and l != (n + 1):
                    t1 = self.fInt(L[j][l], src_matrix, [1, 0, 0])
                    t2 = t3 = t4 = 0
                    t5 = self.fInt(L[j][l - 1], src_matrix, [0, 0, 1])
                    t6 = self.fInt(U[j][l - 1], src_matrix, [0, 1, 0])
                
                elif j == (m + 1) and l != 1 and l != (n + 1):
                    t1 = 0
                    t2 = self.fInt(U[j - 1][l], src_matrix, [0, 0, 1])
                    t3 = self.fInt(L[j - 1][l], src_matrix, [0, 1, 0])
                    t4 = self.fInt(U[j - 1][l - 1], src_matrix, [1, 0, 0])
                    t5 = t6 = 0
                
                elif j != 1 and j != (m + 1) and l == 1:
                    t1 = self.fInt(L[j][l], src_matrix, [1, 0, 0])
                    t2 = self.fInt(U[j - 1][l], src_matrix, [0, 0, 1])
                    t3 = self.fInt(L[j - 1][l], src_matrix, [0, 1, 0])
                    t4 = t5 = t6 = 0
                
                elif j != 1 and j != (m + 1) and l == (n + 1):
                    t1 = t2 = t3 = 0
                    t4 = self.fInt(U[j - 1][l - 1], src_matrix, [1, 0, 0])
                    t5 = self.fInt(L[j][l - 1], src_matrix, [0, 0, 1])
                    t6 = self.fInt(U[j][l - 1], src_matrix, [0, 1, 0])
                
                elif j == 1 and l == 1:
                    t1 = self.fInt(L[j][l], src_matrix, [1, 0, 0])
                    t2 = t3 = t4 = t5 = t6 = 0
                
                elif j == 1 and l == (n + 1):
                    t1 = t2 = t3 = t4 = 0
                    t5 = self.fInt(L[j][l - 1], src_matrix, [0, 0, 1])
                    t6 = self.fInt(U[j][l - 1], src_matrix, [0, 1, 0])
                
                elif j == (m + 1) and l == 1:
                    t1 = 0
                    t2 = self.fInt(U[j - 1][l], src_matrix, [0, 0, 1])
                    t3 = self.fInt(L[j - 1][l], src_matrix, [0, 1, 0])
                    t4 = t5 = t6 = 0
                
                elif j == (m + 1) and l == (n + 1):
                    t1 = t2 = t3 = 0
                    t4 = self.fInt(U[j - 1][l - 1], src_matrix, [1, 0, 0])
                    t5 = t6 = 0
                
                else:
                    t1 = self.fInt(L[j][l], src_matrix, [1, 0, 0])
                    t2 = self.fInt(U[j - 1][l], src_matrix, [0, 0, 1])
                    t3 = self.fInt(L[j - 1][l], src_matrix, [0, 1, 0])
                    t4 = self.fInt(U[j - 1][l - 1], src_matrix, [1, 0, 0])
                    t5 = self.fInt(L[j][l - 1], src_matrix, [0, 0, 1])
                    t6 = self.fInt(U[j][l - 1], src_matrix, [0, 1, 0])
                
                rh[id] = t1 + t2 + t3 + t4 + t5 + t6
        
        return rh

    def fInt(self, T, src_matrix, v):
        f0 = self.f_func(T.x[0], T.y[0],src_matrix)
        f1 = self.f_func(T.x[1], T.y[1],src_matrix)
        f2 = self.f_func(T.x[2], T.y[2],src_matrix)
        
        s = self.polyarea(T.x, T.y)
        
        f_avg = (f0 + f1 + f2) / 3
        v_avg = (v[0], v[1] + v[2]) / 3

        f3 = (f1 + f2) / 2
        f4 = (f0 + f2) / 2
        f5 = (f0 + f1) / 2

        v3 = (v[1] + v[2]) / 2
        v4 = (v[2] + v[0]) / 2
        v5 = (v[0] + v[1]) / 2

        out = (f3*v3 + f4*v4 + f5*v5 + f_avg*v_avg)*s / 4

        return out


    def f_func(self, x, y, src_matrix):
        dx = self.mesh.dx
        dy = self.mesh.dy
        
        left = self.mesh.left
        bottom = self.mesh.bottom
        
        nn = np.round(( x - left )/dx) + 1
        mm = np.round(( y - bottom )/dy) + 1
        
        out = src_matrix[nn][mm]
        return out

    def setA(self, grid):
        m = self.mesh.m
        n = self.mesh.n
        
        A = np.zeros(((m + 1) * (n + 1) * 7, 3))
        list_index = 0
        
        for j in range(1, m + 2):
            for l in range(1, n + 2):
                a = grid[j][l]
                id_ = j + (l - 1) * (m + 1)
                
                # Center
                A[list_index, :] = [id_, id_, a['c']]
                list_index += 1
                
                # West
                if j != 1:
                    A[list_index, :] = [id_, id_ - 1, a['w']]
                    list_index += 1
                
                # Northwest
                if j != 1 and l != (n + 1):
                    A[list_index, :] = [id_, id_ + m, a['nw']]
                    list_index += 1
                
                # North
                if l != (n + 1):
                    A[list_index, :] = [id_, id_ + m + 1, a['n']]
                    list_index += 1
                
                # East
                if j != (m + 1):
                    A[list_index, :] = [id_, id_ + 1, a['e']]
                    list_index += 1
                
                # South
                if l != 1:
                    A[list_index, :] = [id_, id_ - m - 1, a['s']]
                    list_index += 1
                
                # Southeast
                if j != (m + 1) and l != 1:
                    A[list_index, :] = [id_, id_ - m, a['se']]
                    list_index += 1
        
        # Trim the array to remove unused rows
        A = A[:list_index, :]
        return A

    def setB(self, grid, rh):
        m = self.mesh.m
        n = self.mesh.n
        
        B = np.zeros((m + 1) * (n + 1))
        
        for j in range(1, m + 2):
            for l in range(1, n + 2):
                a = grid[j][l]
                id_ = j + (l - 1) * (m + 1)
                B[id_ - 1] = a['const']  # Adjust for zero-indexing
        
        B = rh - B
        return B

    def get_u_val(self, A, B):
        """
        This method is a helper functtion to formulate the mesh
        """
        maximum_iterations = 300
        out = bicgstab(A, B, maxiter=maximum_iterations)
        
        return out

    def get_vn_val(self, u):
        """
        This is a helper function to formulate mesh for simulation object
        """
        m = self.mesh.m
        n = self.mesh.n

        vn = np.zeros((n+1, m+1))
        for ii in range(m+1):
            for jj in range(n+1):
                vn[ii, jj] = u[(jj-1)*(m+1)+ii]
        
        return vn

    def saturation_equ_solver(self, dt, u, v):
        """
        -- Solving Saturation Equations --
        code to compute solution of saturation,concentration and 
        surfactant equations by Modified Method of Characteristics
        using explicit formulation (Yuan Yi-Rang 1993) and implicit finite
        difference method
        """
        g1 = self.init_water_saturation_scalar
        g2 = self.init_water_saturation_scalar * self.polymer.initial_concentration
        g3 = self.init_water_saturation_scalar * self.surfactant.concentration

        if(self.water_saturation is not None):
            m = np.size(self.water_saturation, 1)
            n = np.size(self.water_saturation, 0)
        else:
            raise SimulationCalcInputException("SimulationError: Water Saturation matrix not initialized...")
        
        dx = self.mesh.dx
        dy = self.mesh.dy

        Q = self.water_saturation

        dt_array = dt*np.ones((n,m))
        
        #Defining const. parameters for Pc
        omega1 = 0.1
        omega2 = 0.4

        phi = 1 #porosity
        [x, y] = np.meshgrid(
                np.arange(self.mesh.left, self.mesh.right + self.mesh.dx, self.mesh.dx), 
                np.arange(self.mesh.bottom, self.mesh.top + self.mesh.dy, self.mesh.dy))
        

        #Define critical capillary numbers
        Nco0 = 10**(-5)
        Nca0 = 10**(-5)

        ###PARAMETER DEFINITION
        
        #recompute the residual saturations using (n+1)th time velocities & define the normalized saturrations of water and oil at IFT sigma as:
            # $$ \bar{s} = \frac{s-s_{ra}}{1-s_{ra}} $$
            #
            # $$ \tilde{s} = \frac{s-s_{ra}}{1-s_{ra}-s_{ro}} $$
        [swr, sor] = self.compres(u, v)
        nsw = (Q - swr) / (1-swr)
        nso = (Q - swr) / (1-swr-sor)
        
        #recompute mobilities
        lambda_a = self.compmob(sor, swr, 1)
        lambda_o = self.compmob(sor, swr, 0)
        lambda_total = lambda_a + lambda_o

        #recompute fractional flows (lambda_a / lambda_total)
        f = lambda_a / lambda_total
        
        #Getting KK and Kmax from KK_def() function
        [Kmax, KK] = self.KK_def(x, y)

        D = KK*lambda_o*f

        #Calculating derivative of IFT with respect to surfactant concentration
        derivative_sigma_g = self.derivative_sigma
        
        #compute the capillary number
        


    def KK_def(self, x, y):
        if(self.permeability_flg == PermeabilityType.Homogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            # Represents a homogenous rectilinear model
            Kmax = 1000
            KK = Kmax*np.ones(self.sog+1)
        elif(self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            Kmax = 100
            KK = Kmax*( 0.5*(1-10^(-7))*(np.sin(6*np.pi*np.cos(x))*np.cos(4*np.pi*np.sin(3*y))-1)+1)
        elif(self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Quarter_Five_Spot):
            # need to load the KK30Tabert.mat file... need to use the scipy.io.loadmat() method
            [Kmax, KK] = sp.io.loadmat('./Resources/KK30Tabert.mat')

        return [Kmax, KK]



    def get_gradient(self, vn):
        m = self.mesh.m
        n = self.mesh.n

        dx = self.mesh.dx
        dy = self.mesh.dy

        px = np.zeros((n+1, m+1))
        py = px
        
        for i in range(m + 2):
            for j in range(n + 2):
                if i != 0:
                    px[j, i] = (vn[j, i] - vn[j, i - 1]) / dx
                if i != m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i]) / dx
                if i != 0 and i != m:
                    px[j, i] = (vn[j, i + 1] - vn[j, i - 1]) / (2 * dx)
                if j != 0:
                    py[j, i] = (vn[j, i] - vn[j - 1, i]) / dy
                if j != n:
                    py[j, i] = (vn[j + 1, i] - vn[j, i]) / dy
                if j != 0 and j != n:
                    py[j, i] = (vn[j + 1, i] - vn[j - 1, i]) / (2 * dy)
        return [px, py]
    
    def calculate_shear_effects(self):
        """
        Will calculate the change in the aqueous viscosity due to polymer shear thinning (if flag is turned on)
        """
        pass

    ### HELPER FUNCTIONS TO THE SIMULATION CLASS 
    def divergence(self,F1, F2):
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

    def polyarea(self, x, y):
        """
        Calculate the area of a polygon using the Shoelace formula.
        The vertices are defined by the x and y coordinates.
        
        Parameters:
        x (list or array): x-coordinates of the polygon vertices
        y (list or array): y-coordinates of the polygon vertices
        
        Returns:
        float: Area of the polygon
        """
        return 0.5 * abs(sum(x[i] * y[i + 1] - y[i] * x[i + 1] for i in range(-1, len(x) - 1)))

    
