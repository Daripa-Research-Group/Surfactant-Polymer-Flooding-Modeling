"""
This python script contains the class definition for running simulations
"""

#Relevant imports
from Exceptions import SimulationCalcInputException
from para import Box
import numpy as np
from enumerations import SimulationConstants, PolymerList, ModelType, ResevoirGeometry, PermeabilityType


class Simulation:
    """
    Class is used to generate an instance of a simulation which will run the SP-flooding model based on the given parameters provided by the user
    """
    def __init__(self, sim_id, polymer, surfactant, init_water_saturation, resevoir_geometry, permeability_flg, mesh_grid, is_surfactant, mdl_id, plt_type):
        """
        creates instance of the simulation class which will enable for calculating changes in system parameters at every time-step

        :param sim_id: Simulation number
        :type sim_id: int

        :param polymer: Polymer object used in SP-flooding run
        :type polymer: Polymer

        :param surfactant: Surfactant object used in SP-flooding run
        :type surfactant: Surfactant

        :param init_water_saturation: Initial Water Saturation
        :type init_water_saturation: float

        :param resevoir_geometry: Type of resevoir geometry (is it a rectilinear or quarter-five-spot geometry)
        :type resevoir_geometry: enum 'ResevoirGeometry'

        :param permeability_flg: Homogenous vs. Heterogenous porosity in resevoir
        :type permeability_flg: enum 'PermeabilityType'

        :param mesh_grid: mesh_grid used in the SP-flooding run
        :type mesh_grid: Box

        :param is_surfactant: boolean that states whether there is surfactant in the system or not
        :type is_surfactant: bool

        :param mdl_id: the model id for the simulation
        :type mdl_id: enum 'ModelType'

        :param plt_type: the plot type outputted by the program for the simulation run
        :type plt_type: enum 'PlotType'
        """
        
        #Polymer and Surfactant Properties
        self._polymer_ = None
        self._surfactant_ = None
        self.polymer = polymer
        self.surfactant = surfactant

        #User Inputs
        self._resevoir_geometry_ = None
        self._permeability_flag_ = None
        self._mesh_ = Box()
        self.resevoir_geometry = resevoir_geometry
        self.permeability_flg = permeability_flg
        self.mesh = mesh_grid

        #simulation properties
        self._water_saturation_ = 0
        self._init_water_saturation_scalar_ = init_water_saturation
        self._aqueous_viscosity_ = 0
        self._oleic_mobility_ = 0
        self._aqueous_mobility_ = 0 
        self._phi_ = 0 
        self._sigma_ = 0 # interfacial tension
        
        #General Parameters in Simulation
        self.sim_id = sim_id
        self.is_surfactant = is_surfactant
        self.mdl_id = mdl_id
        self.plt_type = plt_type


    @property
    def polymer(self):
        return self._polymer_

    @polymer.setter
    def polymer(self, value):
        self._polymer_ = value
        
    @property
    def surfactant(self):
        return self._surfactant_

    @surfactant.setter
    def surfactant(self, value):
        self._surfactant_ = value

    @property
    def resevoir_geometry(self):
        return self._resevoir_geometry_

    @resevoir_geometry.setter
    def resevoir_geometry(self, value):
        self._resevoir_geometry_ = value

    @property
    def permeability_flg(self):
        return self._permeability_flag_

    @permeability_flg.setter
    def permeability_flg(self, value):
        self._permeability_flag_ = value

    @property
    def mesh(self):
        return self._mesh_

    @mesh.setter
    def mesh(self, value):
        self._mesh_ = value

    @property
    def phi(self):
        self._phi_ = self.get_phi_value()
        return self._phi_

    @property
    def water_saturation(self):
        return self._water_saturation_

    @water_saturation.setter
    def water_saturation(self, value):
        self._water_saturation_ = value

    @property
    def init_water_saturation_scalar(self):
        return self._init_water_saturation_scalar_

    @init_water_saturation_scalar.setter
    def init_water_saturation_scalar(self, value):
        self._init_water_saturation_scalar_ = value

    @property
    def aqueous_viscosity(self):
        return self._aqueous_viscosity_
    
    @aqueous_viscosity.setter
    def aqueous_viscosity(self, value):
        self._aqueous_viscosity_ = value

    @property
    def oleic_mobility(self):
        return self._oleic_mobility_

    @oleic_mobility.setter
    def oleic_mobility(self, value):
        self._oleic_mobility_ = value

    @property
    def aqueous_mobility(self):
        return self._aqueous_mobility_

    @aqueous_mobility.setter
    def aqueous_mobility(self, value):
        self._aqueous_mobility_ = value

    @property
    def sigma(self):
        return self._sigma_

    @sigma.setter
    def sigma(self, value):
        self._sigma_ = value


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
            if(self.mesh is Box()):
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
            if(self.polymer is not None and self.surfactant is not None and self.water_saturation is not None):
                s_0 = self.init_water_saturation_scalar
                c_0 = self.polymer.concentration
                g_0 = self.surfactant.concentration
                self.water_saturation = np.zeros((self.mesh.m + 1, self.mesh.n + 1))
                self.polymer.vec_concentration = np.copy(self.water_saturation)
            else:
                raise SimulationCalcInputException("SimulationError: Did not provide intial water saturation or Polymer and/or Surfactant were not intialized. Please try again...")

            if(self.phi is not None):
                D = (self.phi > 1e-10) + (np.abs(self.phi) < 1e-10)
            else:
                raise SimulationCalcInputException("SimulationError: phi value not calculated!")

            if self.is_surfactant == 0 and self.polymer is not None and self.surfactant is not None:
                self.surfactant.vec_concentration = []
                self.water_saturation = (~D) + D * (1 - s_0)
                self.polymer.vec_concentration = (~D) * c_0
            elif self.is_surfactant == 1 and self.polymer is not None and self.surfactant is not None:
                self.water_saturation = (~D) + D * (1 - s_0)
                self.polymer.vec_concentration = (~D) * c_0
                self.surfactant.vec_concentration = (~D) * g_0
            else:
                raise SimulationCalcInputException("SimulationError: Surfactant and/or polymer not initialized")
            
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
                vis_water = SimulationConstants.Water_Viscosity.value #viscosity['water']
                vis_oil = SimulationConstants.Oil_Viscosity.value #viscosity['oil']
                # vis_polymer_array = self.polymer.visosity #viscosity['polymer_array']
                # polymer_type = self.polymer.name #params['polymer_type']
                polymer_obj = self.polymer
                # beta1 =  beta1 #params['beta1']
            else:
                raise SimulationCalcInputException("SimulationError: Surfactant and/or Polymer Not Initialized")

            if (self.mdl_id == ModelType.No_Shear_Thinning):
                # Newtonian Model (NO SHEAR THINNING INVOLVED => MODEL TYPE #1)
                n = np.shape(( polymer_obj.vec_concentration,1 ))
                m = np.shape(( polymer_obj.vec_concentration,2 ))
                if polymer_obj.vec_concentration == 0:
                    self.aqueous_viscosity = vis_water * np.ones((n, m))
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
                e_coeff = polymer_obj.epsilon_coefficients
                n_coeff = polymer_obj.n_coefficients
                    
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
        if(self.is_surfactant == 0): # (without surfactant)
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
        elif(self.is_surfactant == 1): # (with surfactant)
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

    def setGrid(self):
        
        pass

    def setRHS(self):
        pass

    def setA(self):
        pass

    def setB(self):
        pass

    def transport_solver(self):
        pass


    def __str__(self):
        return ""

    def export_data(self):
        pass


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
