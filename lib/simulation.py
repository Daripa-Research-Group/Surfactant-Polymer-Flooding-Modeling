"""
This python script contains the class definition for running simulations
"""

import sys
import os
import scipy as sp
from scipy.sparse.linalg import bicgstab
from scipy.sparse import coo_matrix
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

#Relevant imports
from lib.Exceptions import SimulationCalcInputException
from lib.para import Box
from lib.enumerations import SimulationConstants, PolymerList, ModelType, ResevoirGeometry, PermeabilityType
from lib.surfactant import Surfactant
from lib.polymer import Polymer

class Simulation:
    """
    Class is used to generate an instance of a simulation which will run the SP-flooding model based on the given parameters provided by the user
    """
    def __init__(self, sim_id, size_of_grid, polymer : Polymer, surfactant : Surfactant, resevoir_geometry, permeability_flg, mdl_id, plt_type, init_water_saturation = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value, init_oleic_saturation = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value):
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

        :param init_oleic_saturation: Initial oleic phase saturation (scalar quantity)
        :type init_oleic_saturation: float

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
        self.init_water_saturation_scalar = init_water_saturation
        self.init_oleic_saturation_scalar = init_oleic_saturation
        
        #Model and plotting flags
        self.mdl_id = mdl_id #Model type (shear thinning on or off type model)
        self.plt_type = plt_type #types of plots to generate

    ### DEPENDENT VARIABLES OF SIMULATION CLASS:

    #TODO: Need to add oleic viscosity and oleic saturation properties!!!
    _mesh_ = None
    @property
    def mesh(self):
        if(self._mesh_ is None):
            self._mesh_ = Box()
            self._mesh_.calculate_spacing
        return self._mesh_

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

    def execute_simulation(self):
        """
        This method will run the simulation under the parameters provided by the user. This will be the only method called within the 'main' function
        """
        
        #initialize the mesh
        self.mesh
        [x, y] = np.meshgrid(
            np.arange(self.mesh.left, self.mesh.right + self.mesh.dx, self.mesh.dx), 
            np.arange(self.mesh.bottom, self.mesh.top + self.mesh.dy, self.mesh.dy))

        #running the 'get_phi_value' method:
        phi_test = self.get_phi_value()

        #defining rhs of elliptic system (source terms)
        f = np.zeros((self.mesh.n+1,self.mesh.m+1))
        MFW = 0
        iterX_save = 0

        # Magnitude of flowrate at source
        mag_source_flow = SimulationConstants.Source_Flow_Magnitude.value

        # Setting the permeability state
        bool_rectilinear_homogenous = self.permeability_flg == PermeabilityType.Homogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear
        bool_rectilinear_heterogenous = self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear
        bool_quarterfivespot_heterogenous = self.permeability_flg == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Quarter_Five_Spot
 
        if(bool_rectilinear_homogenous or bool_rectilinear_heterogenous):
            f[:,0] = mag_source_flow #intensity of injection well
            f[:,self.mesh.m] = -1*mag_source_flow #intensity of production well
        elif(bool_quarterfivespot_heterogenous):
            f[0,0] = mag_source_flow #intensity of injection well
            f[self.mesh.n, self.mesh.m] = -1*mag_source_flow #intensity of production well

        # Developing the permeability matrix
        Kinfo = self.KK_def(x,y)
        KK = None
        Kmax = None
        if(Kinfo is not None):
            Kmax = Kinfo[0]
            KK = Kinfo[1]
        else:
            raise SimulationCalcInputException("SimulationError: 'KK' and 'Kmax' variables not instantiated...")


        # Initializing the water saturation, polymer concentration, and surfactant concentration matrix
        c0_array = self.polymer.initial_concentration*np.ones((SimulationConstants.Grid_Size.value+1,SimulationConstants.Grid_Size.value+1))
        [water_sat_matrix, polymer_matrix, surfactant_matrix] = self.initial_concentration_matrix()
        interface = np.zeros((60,1))

        # Determining viscosities of oil, water, and polymer
        viscosity_oil = SimulationConstants.Oil_Viscosity.value
        viscosity_water = SimulationConstants.Water_Viscosity.value
        beta_1 = 15000
        self.polymer.viscosity_scalar = viscosity_water*(1*beta_1*polymer_matrix)
        self.polymer.viscosity_matrix = self.polymer.viscosity_scalar*np.ones((SimulationConstants.Grid_Size.value+1,SimulationConstants.Grid_Size.value+1))

        # Defining parameters that need to be updated during each iteration of the while loop
        t = 0
        t_cal = 0
        dt = self.mesh.dx/mag_source_flow
        if(bool_quarterfivespot_heterogenous):
            dt = (self.mesh.dx/mag_source_flow)*100
        tf = 500
        tSave = 0
        u = np.zeros((self.mesh.n+1, self.mesh.m+1))
        v = u
        COC = np.zeros((1, 2000)) #cumulative oil captured
        ProdRate = np.zeros((1, np.floor(tf/dt))) #rate of production
        CROIP = np.zeros((1, np.floor(tf/dt))) #cummulative residual oil in place
        viscosity_aqueous_save = 0
        shear_force_save = 0
        concentration_save = 0
        source_flow_magnitude_total = 0
        sum_of_saturation_matrix = 0
        innerIter_save = []

        #while loop to update parameters during each iteration
            #while t is < tf and water isn't starting to show up at production well, keep iterating:
        if(self.water_saturation is not None and self.aqueous_viscosity is not None):
            while(t < tf and self.water_saturation[self.mesh.n+1, self.mesh.m+1] <=0.70):
                #updating source flow magnitude:
                source_flow_magnitude_total += mag_source_flow

                #updating time index:
                t += dt
                innerIter = 0
                epsilon = 10

                #computing viscosities:
                shear_force = self.compvis(u, v, x, y, beta_1, c0_array)[1]

                if(self.mdl_id == ModelType.Shear_Thinning_On):
                    self.polymer.viscosity_scalar = max(self.aqueous_viscosity[1, :])
                self.polymer.viscosity_matrix = self.polymer.viscosity_scalar*np.ones((SimulationConstants.Grid_Size.value+1,SimulationConstants.Grid_Size.value+1))
                
                #updating tSave
                tSave += t_cal

                #determining size of shear force matrix
                [m_shear, n_shear] = np.shape(shear_force)

                #compute residual saturations
                [swr, sor] = self.compres(u, v)

                #compute mobilities
                lambda_a = self.compmob(swr=swr, sor=sor, flag=1)
                lambda_o = self.compmob(swr=swr, sor=sor, flag=0)
                lambda_total = lambda_a + lambda_o

                #update 'beta' value for polymer flow:
                beta = KK*lambda_total

                #updating grid
                [U, L] = self.setTri()
                grid_size = self.set_grid(U, L, beta)
                rh = self.setRightHand(f, U, L)
                A = self.setA(grid_size)
                A = coo_matrix(A)
                B = self.setB(grid_size, rh)
                uOld = u
                vOld = v
                u = self.get_u_val(A, B)
                vn = self.get_vn_val(u)
                [px, py] = self.get_gradient(vn)
                u = -1*beta*px
                v = -1*beta*py
                if(innerIter > 1000):
                    break
                innerIter+=1
                innerIter_save[t_cal+1] = innerIter

                #Solving transport problem:


                #Save relevant results in each iteration for plotting










        pass


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

    def compres(self, u, v, sigma_mod=None):
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
        if(sigma_mod is None):
            nca = np.sqrt(u**2 + v**2) * self.aqueous_viscosity / self.sigma
            nco = np.sqrt(u**2 + v**2) * miuo / self.sigma
            Nca = np.linalg.norm(nca)  # compute 2-norm (largest singular value)
            Nco = np.linalg.norm(nco)
        else:
            nca = np.sqrt(u**2 + v**2) * self.aqueous_viscosity / sigma_mod
            nco = np.sqrt(u**2 + v**2) * miuo / sigma_mod
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
        miuo = SimulationConstants.Oil_Viscosity.value


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
        norm_nco0 = 10**(-5)
        norm_nca0 = 10**(-5)

        ###PARAMETER DEFINITION
        
        #recompute the residual saturations using (n+1)th time velocities & define the normalized saturations of water and oil at IFT sigma as:
            # $$ \bar{s} = \frac{s-s_{ra}}{1-s_{ra}} $$
            #
            # $$ \tilde{s} = \frac{s-s_{ra}}{1-s_{ra}-s_{ro}} $$
        swr0 = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
        sor0 = SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value
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
        KK_info = self.KK_def(x, y)
        if(KK_info is not None):
            Kmax = KK_info[0]
            KK = KK_info[1]
        else:
            raise SimulationCalcInputException("SimulationError: KK and Kmax not calculated properly...")

        D = KK*lambda_o*f

        #Calculating derivative of IFT with respect to surfactant concentration
        sigma_g = self.derivative_sigma
        
        #compute the capillary number
        nca = np.sqrt((u**2)+(v**2))*self.aqueous_viscosity/self.sigma
        nco = np.sqrt((u**2)+(v**2))*miuo/self.sigma
        norm_nca = np.linalg.norm(nca, ord=2)
        norm_nco = np.linalg.norm(nco, ord=2)

        #recalculating derivative of residual saturation with respect to surfactant concentration
        swr_g = np.zeros((n,m))
        sor_g = swr_g

        if(self.aqueous_viscosity is not None and self.surfactant.vec_concentration is not None):
            for i in range(n): #traversing column
                for j in range(m): #traversing row
                    if(norm_nca >= norm_nca0):
                        swr_g[j][i] = - (swr0*0.1534*10.001*norm_nca0**0.1534)/(np.sqrt(u[j][i]**2 + v[j][i]**2)*self.aqueous_viscosity[j][i]**(0.1534)*self.sigma[j][i]**(0.8466)*( self.surfactant.vec_concentration[j][i] + 1 )**2) 
                    elif(norm_nco >= norm_nco0):
                        swr_g[j][i] = - (sor0*0.5213*10.001*norm_nco0**0.5213)/(np.sqrt(u[j][i]**2 + v[j][i]**2)*self.aqueous_viscosity[j][i]**(0.5213)*self.sigma[j][i]**(0.4787)*( self.surfactant.vec_concentration[j][i] + 1 )**2) 
        
        # Determining the derivatives of normalized saturations with respect to surfactant concentration
        nsw_g = swr_g * (Q - 1)/ ((1-swr)**2)
        nso_g = (swr_g*(Q+sor-1)+sor_g*(Q-swr))/(1-swr-sor)**2

        # Determining relative permeability with respect to saturation
        kra_s = 2.5*swr*(3*(nsw)**2-1)+1
        kro_s = 10*sor*nso-5*sor-1

        # Determining derivative of relative permeability with respect to surfactant concentration
        kra_g = 2.5*swr_g*(nsw**2-nsw) + (Q - 1)*(2.5*swr*(3*(nsw**2)-1)+1)*nsw_g / ((1-swr)**2)
        kro_g = 1-5*sor*nso+(1-nso)*(105*nso*sor_g)-(1+5*sor-10*sor*nso)*nso_g

        # Determining derivative of fraction flow function with respect to saturation, concentration, and surfactant concentration
        f_s = kra_s*lambda_o/((lambda_total**2)*self.aqueous_viscosity) - kro_s*lambda_a/(( lambda_total**2 )*miuo)

        # Determining capillary pressure and its derivatives
        pc = (self.sigma*omega2*self.phi**(0.5))/(KK**(0.5)*(1-nso)**(1/omega1))
        pc_s = pc/(omega1*(1-nso))# derivative of capillary pressure with respect to saturation
        pc_g = (pc/self.sigma)*sigma_g + pc_s

        #Solving for saturation
            #Will use characteristic equations and finite difference discretization

        [xmod, ymod] = self.eval_Xsurf_neumann(
                                            flag=1,
                                            x=x,
                                            y=y,
                                            s=Q,
                                            snew=Q,
                                            g=self.surfactant.vec_concentration,
                                            f=f,
                                            f_s=f_s,
                                            D=D,
                                            pc_s=pc_s,
                                            pc_g=pc_g,
                                            u=u,
                                            v=v,
                                            dt=dt
                                        )
        #performing 2-D interpolation using the scipy.interpolate package
        interp = sp.interpolate.RegularGridInterpolator((x,y), Q)
        Qmod = interp((xmod, ymod))

        #Recalculate the normalized saturation of water and oil
        nso = (Qmod - swr)/(1-swr-sor)

        #Updating coefficients with interpolated saturations
        lambda_a = self.compmob(sor,swr, 1)
        lambda_o = self.compmob(sor, swr, 0)
        lambda_total = lambda_a + lambda_o
        f = lambda_a/lambda_total #fractional flow of the aqueous phase
        [Kmax, KK] = self.KK_def(x, y)
        D = KK*lambda_o*f
        pc_s = pc/(omega1*(1-nso))# derivative of capillary pressure with respect to saturation
        pc_g = (pc/self.sigma)*sigma_g + pc_s
        f_c = (lambda_o*lambda_a*miuo)/(( lambda_total**2 )*self.aqueous_viscosity)
        f_g = (kra_g*lambda_o)/(( lambda_total**2 )*self.aqueous_viscosity) - (kro_g*lambda_a)/(( lambda_total**2 )*miuo)

        #intermediate calculation of param 'D'
        D_g = D*pc_g
        D_s = D*pc_s

        idx = 1
        AAA = np.zeros((n*m))
        DDD = np.zeros_like((n*m,1))
        
        while(idx <= (m)*(n-1)+1 and self.surfactant.vec_concentration is not None and self.polymer.vec_concentration is not None):
            cnt = (idx - 1) / m # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n,m))
            AA = BB
            CC = BB
            DD = np.zeros((m,1))

            #'cnt+1' in matlab is 'cnt' in python as matlab indexes from 1 but python indexes from 0
            for i in range(m):
                for j in range(n):
                    if(idx == 1):
                        if(i == 0): #first/left column
                            DD[i] = (Qmod[cnt][i]/dt[cnt][i]) + g1*(1-f[cnt][i]) \
                                    + ((D_g[cnt][i]+D_g[cnt][i+1])/(dx**2) + (D_g[cnt+1][i] + D_g[cnt+1][i])/(dx**1) )*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i] + D_g[cnt][i+1])/(dx**2)*self.surfactant.vec_concentration[cnt][i+1] \
                                    - (D_g[cnt][i] + D_g[cnt+2][i])/(dy**2)*self.surfactant.vec_concentration[cnt+2][i]

                            CC[j][i] = (D_s[cnt][i] + D_s[cnt+1][i])/(dy**2)

                            BB[j][i] = 1/dt[cnt][i] - (D_s[cnt+1][i] + D_s[cnt][i+1])/(dx**2) \
                                    - (D_s[cnt+1][i] + D_s[cnt][i])/(dy**2)

                            BB[j][i+1] = (D_s[cnt][i] + D_s[cnt][i+1])*(dx**2)
                        elif(i==m): #last/rightmost column 
                            DD[i] = Qmod[cnt][i]/dt[cnt][i] \
                                    + ((D_g[cnt][i] + D_g[cnt][i-1])/(dx**2) + (D_g[cnt+1][i] + D_g[cnt][i])/(dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i] + D_g[cnt][i-1])/(dx**2)*self.surfactant.vec_concentration[cnt][i-1] \
                                    - (D_g[cnt][i]+ D_g[cnt+1][i])/(dy**2)*self.surfactant.vec_concentration[cnt+1][i]
                            
                            BB[j][i-1] = (D_s[cnt][i] + D_s[cnt][i-1])/(dx**2)

                            BB[i][i] = 1/dt[cnt][i] - (D_s[cnt][i] + D_s[cnt][i-1])/(dx**2) \
                                    - (D_s[cnt+1][i] + D_s[cnt][i])/(dy**2)

                            CC[j][i] = (D_s[cnt][i] + D_s[cnt+1][i])/(dy**2)
                        else:
                            DD[i] = Qmod[cnt][i]/dt[cnt][i] \
                                    - f_c[cnt][i]*(u[cnt][i]*(self.polymer.vec_concentration[cnt][i+1] - self.polymer.vec_concentration[cnt][i-1])/(2*dx)) \
                                    - f_g[cnt][i]*(u[cnt][i]*(self.surfactant.vec_concentration[cnt][i+1] - self.surfactant.vec_concentration[cnt][i-1])/(2*dx)) \
                                    + ((D_g[cnt][i+1]+D_g[cnt][i-1] + 2*D_g[cnt][i])/(2*dx**2)+(D_g[cnt][i+1]+D_g[cnt][i])/(dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i+1]+D_g[cnt][i])/(2*dx**2)*self.surfactant.vec_concentration[cnt][i+1] \
                                    - (D_g[cnt][i-1]+D_g[cnt][i])/(2*dx**2)*self.surfactant.vec_concentration[cnt][i-1] \
                                    - (D_g[cnt][i]+D_g[cnt+1][i])/(dy**2)*self.surfactant.vec_concentration[cnt+1][i]

                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i+1]+D_s[cnt][i-1]+2*D_s[cnt][i])/(2*dx**2)-(D_s[cnt+1][i]+D_s[cnt][i])/(dy**2)
                            BB[j][i-1] = (D_s[cnt][i-1]+D_s[cnt][i])/(2*dx**2)
                            BB[j][i+1] = (D_s[cnt][i+1]+D_s[cnt][i])/(2*dx**2)

                            CC[j][i] = (D_s[cnt][i]+D_s[cnt+1][i])/(dy**2)
                    elif(idx == (m)*(n-1)+1): #topmost row of grid
                        if(i == 0): #leftmost column
                            DD[i] = (Qmod[cnt][i]/dt[cnt][i]) \
                                    + ((D_g[cnt][i] + D_g[cnt][i+1])/(dx**2) + (D_g[cnt-1][i]+D_g[cnt][i])/(dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i]+D_g[cnt][i+1])/(dx**2)*self.surfactant.vec_concentration[cnt][i+1]\
                                    - (D_g[cnt][i] + D_g[cnt-1][i])/(dy**2)*self.surfactant.vec_concentration[cnt-1][i]

                            AA[j][i] = (D_s[cnt][i]+D_s[cnt-1][i])/(dy**2)
                            
                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i]+D_s[cnt][i+1])/(dx**2)\
                                    -(D_s[cnt-1][i]+D_s[cnt][i])/(dy**2)
                            
                            BB[j][i+1] = (D_s[cnt][i]+D_s[cnt][i+1])/(dx**2)
                        elif(i == m - 1): #rightmost column
                            DD[i] = (Qmod[cnt][i]/dt[cnt][i])\
                                    + ((D_g[cnt][i] + D_g[cnt][i-1])/(dx**2) + (D_g[cnt-1][i]+D_g[cnt][i])/(dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i]+D_g[cnt][i-1])/(dx**2)*self.surfactant.vec_concentration[cnt][i-1]\
                                    - (D_g[cnt][i] + D_g[cnt-1][i])/(dy**2)*self.surfactant.vec_concentration[cnt-1][i]

                            BB[j][i-1] = (D_s[cnt][i]+D_s[cnt][i-1])/(dy**2)

                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i]+D_s[cnt][i-1])/(dx**2)\
                                    -(D_s[cnt-1][i]+D_s[cnt][i])/(dy**2)

                            AA[j][i] = (D_s[cnt][i]+D_s[cnt-1][i])/(dy**2)
                        else:
                            DD[i] = Qmod[cnt][i]/dt[cnt][i] \
                                    - f_c[cnt][i]*(u[cnt][i]*(self.polymer.vec_concentration[cnt][i+1] - self.polymer.vec_concentration[cnt][i-1])/(2*dx)) \
                                    - f_g[cnt][i]*(u[cnt][i]*(self.surfactant.vec_concentration[cnt][i+1] - self.surfactant.vec_concentration[cnt][i-1])/(2*dx)) \
                                    + ((D_g[cnt][i+1]+D_g[cnt][i-1] + 2*D_g[cnt][i])/(2*dx**2)+(D_g[cnt][i+1]+D_g[cnt][i])/(dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i+1]+D_g[cnt][i])/(2*dx**2)*self.surfactant.vec_concentration[cnt][i+1] \
                                    - (D_g[cnt][i-1]+D_g[cnt][i])/(2*dx**2)*self.surfactant.vec_concentration[cnt][i-1] \
                                    - (D_g[cnt][i]+D_g[cnt-1][i])/(dy**2)*self.surfactant.vec_concentration[cnt-1][i]

                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i+1]+D_s[cnt][i-1]+2*D_s[cnt][i])/(2*dx**2)-(D_s[cnt-1][i]+D_s[cnt][i])/(dy**2)
                            BB[j][i-1] = (D_s[cnt][i-1]+D_s[cnt][i])/(2*dx**2)
                            BB[j][i+1] = (D_s[cnt][i+1]+D_s[cnt][i])/(2*dx**2)

                            AA[j][i] = (D_s[cnt][i]+D_s[cnt-1][i])/(dy**2)
                    else:
                        if(i == 0):
                            DD[i] = Qmod[cnt][i]/dt[cnt][i] \
                                    - f_c[cnt][i]*(v[cnt][i]*(self.polymer.vec_concentration[cnt+1][i] - self.polymer.vec_concentration[cnt][i])/(2*dy)) \
                                    - f_g[cnt][i]*(v[cnt][i]*(self.surfactant.vec_concentration[cnt+1][i] - self.surfactant.vec_concentration[cnt][i])/(2*dy)) \
                                    + ((D_g[cnt][i]+D_g[cnt][i+1])/(dx**2)+(D_g[cnt-1][i]+2*D_g[cnt][i]+D_g[cnt+1][i])/(2*dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i]+D_g[cnt][i+1])/(dx**2)*self.surfactant.vec_concentration[cnt][i+1] \
                                    - (D_g[cnt][i]+D_g[cnt+1][i])/(2*dy**2)*self.surfactant.vec_concentration[cnt+1][i] \
                                    - (D_g[cnt][i]+D_g[cnt-1][i])/(2*dy**2)*self.surfactant.vec_concentration[cnt-1][i]

                            AA[j][i] = (D_s[cnt][i]+D_s[cnt-1][i])/(2*dy**2)

                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i]+D_s[cnt][i+1])/(dx**2)-(D_s[cnt-1][i]+2*D_s[cnt][i]+D_s[cnt+1][i])/(2*dy**2)
                            BB[j][i+1] = (D_s[cnt][i+1]+D_s[cnt][i])/(dx**2)

                            CC[j][i] = (D_s[cnt][i]+D_s[cnt+1][i])/(2*dy**2)
                        elif(i == m - 1):
                            DD[i] = Qmod[cnt][i]/dt[cnt][i] \
                                    - f_c[cnt][i]*(v[cnt][i]*(self.polymer.vec_concentration[cnt+1][i] - self.polymer.vec_concentration[cnt][i])/(2*dy)) \
                                    - f_g[cnt][i]*(v[cnt][i]*(self.surfactant.vec_concentration[cnt+1][i] - self.surfactant.vec_concentration[cnt][i])/(2*dy)) \
                                    + ((D_g[cnt][i]+D_g[cnt][i-1])/(dx**2)+(D_g[cnt-1][i]+2*D_g[cnt][i]+D_g[cnt+1][i])/(2*dy**2))*self.surfactant.vec_concentration[cnt][i] \
                                    - (D_g[cnt][i]+D_g[cnt][i-1])/(dx**2)*self.surfactant.vec_concentration[cnt][i-1] \
                                    - (D_g[cnt][i]+D_g[cnt+1][i])/(2*dy**2)*self.surfactant.vec_concentration[cnt+1][i] \
                                    - (D_g[cnt][i]+D_g[cnt-1][i])/(2*dy**2)*self.surfactant.vec_concentration[cnt-1][i]

                            AA[j][i] = (D_s[cnt][i]+D_s[cnt-1][i])/(2*dy**2)

                            BB[j][i] = 1/dt[cnt][i]-(D_s[cnt][i]+D_s[cnt][i-1])/(dx**2)-(D_s[cnt-1][i]+2*D_s[cnt][i]+D_s[cnt+1][i])/(2*dy**2)
                            BB[j][i-1] = (D_s[cnt][i]+D_s[cnt][i-1])/(dx**2)

                            CC[j][i] = (D_s[cnt][i]+D_s[cnt+1][i])/(2*dy**2)
                        else:
                            DD[i] = (Qmod[cnt][i]/dt[cnt][i]) \
                                    - f_c[cnt][i]*(u[cnt][i]*(self.polymer.vec_concentration[cnt][i+1]-self.polymer.vec_concentration[cnt][i-1])/(2*dx) \
                                    +v[cnt][i]*(self.polymer.vec_concentration[cnt+1][i]-self.polymer.vec_concentration[cnt-1][i])/(2*dy)) \
                                    - f_g[cnt][i]*(u[cnt][i]*(self.surfactant.vec_concentration[cnt][i+1] - self.surfactant.vec_concentration[cnt][i-1])/(2*dx) \
                                    +v[cnt][i]*(self.surfactant.vec_concentration[cnt+1][i]-self.surfactant.vec_concentration[cnt-1][i])/(2*dy)) \
                                    - (D_g[cnt][i+1]/(2*dx**2)*(self.surfactant.vec_concentration[cnt][i+1]-self.surfactant.vec_concentration[cnt][i]) \
                                    -D_g[cnt][i-1]/(2*dx**2)*(self.surfactant.vec_concentration[cnt][i-1]-self.surfactant.vec_concentration[cnt][i]) \
                                    +D_g[cnt][i-1]/(2*dx**2)*(self.surfactant.vec_concentration[cnt][i-1]-self.surfactant.vec_concentration[cnt][i+1]) \
                                    +D_g[cnt+1][i]/(2*dx**2)*(self.surfactant.vec_concentration[cnt+1][i]-self.surfactant.vec_concentration[cnt][i]) \
                                    +D_g[cnt-1][i]/(2*dx**2)*(self.surfactant.vec_concentration[cnt-1][i]-self.surfactant.vec_concentration[cnt][i]) \
                                    +D_g[cnt][i]/(2*dx**2)*(self.surfactant.vec_concentration[cnt+1][i]-self.surfactant.vec_concentration[cnt][i]))

                            AA[j][i] = (D_s[cnt-1][i]+D_s[cnt][i])/(2*dy**2)

                            CC[j][i] = (D_s[cnt][i]+D_s[cnt+1][i])(2*dy**2)

                            BB[j][i] = 1/dt[cnt][i]-((1/(2*dx**2))*(D_s[cnt][i]+2*D_s[cnt][i]+D_s[cnt][i+1])+(1/(2*dy**2))*(D_s[cnt-1][i]+2*D_s[cnt][i]+D_s[cnt+1][i]))
                            BB[j][i+1] = (D_s[cnt][i]+D_s[cnt][i+1])/(2*dx**2)
                            BB[j][i-1] = (D_s[cnt][i-1]+D_s[cnt][i])/(2*dx**2)
            if cnt == 0:
                AAA[:n, :2*m] = np.hstack([BB, CC])
            elif cnt == n - 1:
                AAA[(m-1)*n:m*n, (n-2)*m:n*m] = np.hstack([AA, BB])
            else:
                AAA[cnt*n:(cnt+1)*n, (cnt-1)*m:(cnt+2)*m] = np.hstack([AA, BB, CC])

            DDD[cnt*m:(cnt+1)*m] = DD
            idx += m
        
        #Entering the bicgstab for the saturation calculations
            # bicgstab (Biconjugate Gradient Stabilized) - Iterative algorithm to solve large, sparse, and non-symmetric linear systems of the form Ax = b

        Qnew_flat, info = bicgstab(AAA,DDD,rtol=10**(-10),maxiter=600)
        Qnew = Qnew_flat = Qnew_flat.reshape(m, n)

        Qnew[Qnew > 1] = 1

        Smax = np.max(Qnew)
        Smin = np.min(Qnew) 

        #Solving for concentration of polymer in resevoir
        
        [xmod2, ymod2] = self.eval_Xsurf_neumann(
                                            flag=2,
                                            x=x,
                                            y=y,
                                            s=Q,
                                            snew=Qnew,
                                            g=self.surfactant.vec_concentration,
                                            f=f,
                                            f_s=f_s,
                                            D=D,
                                            pc_s=pc_s,
                                            pc_g=pc_g,
                                            u=u,
                                            v=v,
                                            dt=dt
                                        )
        
        interp = sp.interpolate.RegularGridInterpolator((x,y), self.polymer.vec_concentration)
        Cmod = interp((xmod2, ymod2))

        idx = 1
        AAA = np.zeros((n*m))
        DDD = np.zeros_like((n*m,1))

        while(idx <=(m)*(n-1)+1):
            cnt = (idx - 1) / m # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n,m))
            AA = BB
            CC = BB
            DD = np.zeros((m,1))
            for i in range(m):
                for j in range(n):
                    if(j == i):
                        if(idx == 1): #lowermost row of grid
                            if(i == 1): #leftmost point (source)
                                DD[i] = g2/Qnew[cnt][i]+Cmod[cnt][i]/dt[cnt][i]
                                BB[j][i] = 1/dt[cnt][i] + g1/Qnew[cnt][i]
                            else:
                                DD[i] = Cmod[cnt][i] / dt[cnt][i]
                                BB[j][i] = 1/dt[cnt][i]
                        elif(idx == (m)*(n-1)+1):
                            if(i == m - 1):
                                DD[i] = Cmod[cnt][i] / dt[cnt][i]
                                BB[j][i] = 1/dt[cnt][i] - g1*f[cnt][i]/Qnew[cnt][i]
                            else:
                                DD[i] = Cmod[cnt][i]/dt[cnt][i]
                                BB[j][i] = 1/dt[cnt][i]
                        else:
                            DD[i] = Cmod[cnt][i]/dt[cnt][i]
                            BB[j][i] = 1/dt[cnt][i]

            if cnt == 0:
                AAA[0:n, 0:2*m] = np.hstack([BB, CC])
            elif cnt == n-1:
                AAA[(m-1)*n:m*n, (n-2)*m:n*m] = np.hstack([AA, BB])
            else:
                AAA[cnt*n:(cnt+1)*n, (cnt-1)*m:(cnt+2)*m] = np.hstack([AA, BB, CC])

            DDD[cnt*m:(cnt+1)*m] = DD

            idx += m

        Cnew_flat, info = bicgstab(AAA, DDD, rtol=10**(-10), maxiter=600)
        Cnew = Cnew_flat.reshape(m,n)


        [xmod2, ymod2] = self.eval_Xsurf_neumann(
                                            flag=3,
                                            x=x,
                                            y=y,
                                            s=Q,
                                            snew=Qnew,
                                            g=self.surfactant.vec_concentration,
                                            f=f,
                                            f_s=f_s,
                                            D=D,
                                            pc_s=pc_s,
                                            pc_g=pc_g,
                                            u=u,
                                            v=v,
                                            dt=dt
                                        )
        interp = sp.interpolate.RegularGridInterpolator((x,y), self.surfactant.vec_concentration)
        Gmod = interp((xmod2, ymod2))
        
        #Updating coefficients using interpolated surfactant concentration
        sigma_mod = self.surfactant.IFT_conc_equ(Gmod)
        sigma_g_mod = self.surfactant.derivative_IFT_conc_equ(Gmod)
        [swr, sor] = self.compres(u,v,sigma_mod)
        nso = (Qmod-swr)/(1-swr-sor)
        lambda_a = self.compmob(sor, swr, 1)
        lambda_o = self.compmob(sor, swr, 0)
        lambda_total = lambda_a + lambda_o
        f = lambda_a/lambda_total
        KK_info = self.KK_def(x, y)
        if(KK_info is not None):
            Kmax = KK_info[0]
            KK = KK_info[1]
        else:
            raise SimulationCalcInputException("SimulationError: KK and Kmax not calculated properly...")
        D = KK*lambda_o*f
        pc_s = pc/(omega1*(1-nso))
        pc_g = (pc/sigma_mod)*sigma_g_mod + pc_s

        #intermediate parameters for code:
        F = D*pc_g/Qnew
        idx = 1
        AAA = np.zeros(n*m)
        DDD = np.zeros(n*m)[0]

        while(idx <= (m)*(n-1)+1):
            cnt = (idx - 1) / m # cnt = 0, 1, 2, ... for idx = 1, m+1, 2m+1, 3m+1, ...
            BB = np.zeros((n,m))
            AA = BB
            CC = BB
            DD = np.zeros((m,1))
            for i in range(m):
                for j in range(n):
                    if(j == i):
                        if(idx == 1):
                            if(i == 0):
                                DD[i] = g3/Qnew[cnt][i] + Gmod[cnt][i]/dt[cnt][i]
                                CC[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i] + g1/Qnew[cnt][i]
                                BB[j][i+1] = 2*F[cnt][i]/(dx**2)
                            elif(i == m - 1): # Bottom right point
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                CC[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i-1] = 2*F[cnt][i]/(dx**2)
                            else:
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                CC[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i-1] = F[cnt][i]/(dx**2)
                                BB[j][i+1] = 2*F[cnt][i]/(dx**2)
                        elif(idx == (m)*(n-1)+1):
                            if(i == 0):
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i] + g1/Qnew[cnt][i]
                                BB[j][i+1] = 2*F[cnt][i]/(dx**2)
                            elif(i == m-1):
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i] \
                                        - (g1*lambda_a[cnt][i])/[lambda_total[cnt][i]])/Qnew[cnt][i] \
                                        + (g3*lambda_a[cnt][i])/[lambda_total[cnt][i]])/(Qnew[cnt][i]*self.surfactant.concentration)
                                BB[j][i-1] = 2*F[cnt][i]/(dx**2)
                            else:
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = 2*F[cnt][i]/(dy**2)
                                BB[j][i+1] = F[cnt][i]/(dx**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i-1] = F[cnt][i]/(dx**2)
                        else:
                            if(i == 0):
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i+1] = 2*F[cnt][i]/(dx**2)
                                CC[j][i] = F[cnt][i]/(dy**2)
                            elif(i == m - 1):
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i-1] = 2*F[cnt][i]/(dx**2)
                                CC[j][i] = F[cnt][i]/(dy**2)
                            else:
                                DD[i] = Gmod[cnt][i]/dt[cnt][i]
                                AA[j][i] = F[cnt][i]/(dy**2)
                                BB[j][i] = 1/dt[cnt][i] - ((2/(dx**2))+(2/(dy**2)))*F[cnt][i]
                                BB[j][i-1] = 2*F[cnt][i]/(dx**2)
                                BB[j][i+1] = 2*F[cnt][i]/(dx**2)
                                CC[j][i] = F[cnt][i]/(dy**2)
            if cnt == 0:
                AAA[:n, :2*m] = np.hstack([BB, CC])
            elif cnt == n - 1:
                AAA[(m-1)*n:m*n, (n-2)*m:n*m] = np.hstack([AA, BB])
            else:
                AAA[cnt*n:(cnt+1)*n, (cnt-1)*m:(cnt+2)*m] = np.hstack([AA, BB, CC])

            DDD[cnt*m:(cnt+1)*m] = DD
            idx += m

        Gnew_flat, info = bicgstab(AAA,DDD,rtol=10**(-10),maxiter=600)
        Gnew = Gnew_flat.reshape(m, n)

        self.water_saturation = Qnew
        self.polymer.vec_concentration = Cnew
        self.surfactant.vec_concentration = Gnew

        ocut = lambda_o[n][m]*self.init_water_saturation_scalar/lambda_total[n][m] # volume of oil recovered in production well
        wcut = lambda_a[n][m]*self.init_water_saturation_scalar/lambda_total[n][m] # volume of water recovered in production well
        roip = 100*np.sum(np.sum(1-Qnew))/sum(np.ones((n*m,1))) # Oil still in place as a percentage of volume fraction

        return_dict = {
                'water_saturation' : self.water_saturation,
                'polymer_vec_saturation' : self.polymer.vec_concentration,
                'surfactant_vec_saturation' : self.surfactant.vec_concentration,
                'prod_oil_vol' : ocut,
                'prod_water_vol' : wcut,
                'ROIP' : roip
                }
        
        return return_dict


    def eval_Xsurf_neumann(self, flag, x, y, s, snew, g, f, f_s, D, pc_s, pc_g, u, v, dt):
        x_jump = None
        y_jump = None
        if(flag == 1):
            x_jump = x - f_s * u * dt
            y_jump = y - f_s * v * dt
        elif(flag == 2):
            [sx, sy] = self.get_gradient(s)
            [gx, gy] = self.get_gradient(g)

            x_jump = x - ((f/snew)*u + (D*pc_s/snew)*sx + (D*pc_g/snew)*gx)*dt
            y_jump = y - ((f/snew)*v + (D*pc_s/snew)*sy + (D*pc_g/snew)*gy)*dt
        elif(flag == 3):
            [sx, sy] = self.get_gradient(s)

            x_jump = x - ((f/snew)*u + (D*pc_s/snew)*sx)*dt
            y_jump = y - ((f/snew)*v + (D*pc_s/snew)*sy)*dt
        
        xmod = x
        ymod = y
        if(x_jump is not None and y_jump is not None):
            for j in range(np.shape(y)[0]):
                for i in range(np.shape(x)[1]):
                    if(x_jump[j][i] <= 1 and y_jump[j][i] <= 1):
                        xmod[j][i] = abs(x_jump[j][i])
                        ymod[j][i] = abs(y_jump[j][i])
                    elif(x_jump[j][i] > 1 and y_jump[j][i] <= 1):
                        xmod[j][i] = 2 - x_jump[j][i]
                        ymod[j][i] = abs(y_jump[j][i])
                    elif(x_jump[j][i] <= 1 and y_jump[j][i] > 1):
                        xmod[j][i] = abs(x_jump[j][i])
                        ymod[j][i] = 2 - y_jump[j][i]
                    elif(x_jump[j][i] > 1 and y_jump[j][i] > 1):
                        xmod[j][i] = 2 - x_jump[j][i]
                        ymod[j][i] = 2 - y_jump[j][i]
        
        return [xmod, ymod]
                        
    def KK_def(self, x, y):
        Kmax = None
        KK = None
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
        
        if(Kmax is not None and KK is not None):
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
    
    # def calculate_shear_effects(self):
    #     """
    #     Will calculate the change in the aqueous viscosity due to polymer shear thinning (if flag is turned on)
    #     """
    #     pass

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

    
