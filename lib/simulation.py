"""
This python script contains the class definition for running simulations

This Python code has been derived from the MATLAB Surfactant-Polymer Flooding Simulation 
developed by Sourav Dutta and Rohit Mishra.

@author: Bhargav Akula Ramesh Kumar and Carlos Acosta Caripo

"""

import os

#Relevant imports
from lib.Exceptions import SimulationCalcInputException
from lib.para import Box
from lib.enumerations import SimulationConstants, PolymerList, ModelType, ResevoirGeometry, PermeabilityType, PlotType
from lib.polymer import Polymer
from lib.surfactant import Surfactant


class Simulation:
    """
    Simulation class to run SP-flooding simulations based on MATLAB translation.
    """
    def __init__(self, sim_id : int, size_of_grid : int, polymer : Polymer, surfactant : Surfactant, init_water_saturation : float, resevoir_geometry : ResevoirGeometry, permeability_type : PermeabilityType, mesh_grid : Box, model_type : ModelType, plot_type : PlotType):
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

        :param permeability_type: Homogenous vs. Heterogenous porosity in resevoir
        :type permeability_type: enum 'PermeabilityType'

        :param mesh_grid: mesh_grid used in the SP-flooding run
        :type mesh_grid: Box

        :param model_type: the model id for the simulation (whether shear thinning is on or off)
        :type model_type: enum 'ModelType'

        :param plot_type: the plot type outputted by the program for the simulation run
        :type plot_type: enum 'PlotType'
        """
        ## Performs checks on the user input dictionary passed in:
        try:
            model_type = ModelType(user_input_dict["model_type"])
            if model_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadModelTypeAssignment", user_input_dict
            )

        try:
            reservoir_geometry = ResevoirGeometry(user_input_dict["reservoir_geometry"])
            if reservoir_geometry is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadReservoirGeometryAssignment", user_input_dict
            )

        #Simulaiton Properties
        self.resevoir_geometry = resevoir_geometry
        self.permeability_type = permeability_type
        self.mesh = mesh_grid
        self.init_water_saturation_scalar = init_water_saturation
        
        #Model and plotting flags
        self.model_type = model_type
        self.plot_type = plot_type #types of plots to generate

        try:
            polymer_type = PolymerList.get_by_value(user_input_dict["polymer_type"])
            if polymer_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadPolymerTypeAssignment", user_input_dict
            )

        try:
            polymer_concentration = user_input_dict["polymer_concentration"]
            if polymer_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadPolymerConcentrationAssignment", user_input_dict
            )

        try:
            surfactant_type = SurfactantList.get_by_value(
                user_input_dict["surfactant_type"]
            )
            if surfactant_type is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadSurfactantTypeAssignment", user_input_dict
            )

        try:
            surfactant_concentration = user_input_dict["surfactant_concentration"]
            if surfactant_concentration is None:
                raise ValueError
        except (KeyError, ValueError, TypeError):
            raise UserInputException(
                "UserInputError:BadSurfactantConcentrationAssignment", user_input_dict
            )

        ## Instantiates the required simulation properties:

        # initializing properties that hold simulation constants
        self.grid_size = SimulationConstants.Grid_Size.value
        self.source_flow_magnitude = SimulationConstants.Source_Flow_Magnitude.value

        # Initializing Simulation Flags:
        self.permeability_flag = permeability_flag
        self.reservoir_geometry = reservoir_geometry
        self.model_type = model_type
        self.relative_permeability_formula = (
            RelativePermeabilityFormula.AmaefuleHandEquation
        )  # WILL NEED TO UPDATE TO INCLUDE IN GUI

        # Initializing sim properties
        self.mesh, self.FE_mesh = self._create_mesh()
        grid_shape = (self.mesh.n, self.mesh.m)
        self.x, self.y = self.mesh.get_meshgrid
        self.phi = None  # FIXME: Should check if this functionality code
        self.KK = None  # Permeability tensor
        self.time_step = None
        self.u, self.v = self._initialize_pressure_and_velocity()
        self._initialize_simulation()  # will initialize phi, KK, and time_step
        if (
            self.phi is None or self.KK is None or self.time_step is None
        ):  # Raise Exception if not properly initialized...
            raise SimulationCalcInputException(
                "SimulationCalcInputError:BadInitialSimulationPropertiesCalculation"
            )

        # Initalizing Polymer Object
        self.polymer = Polymer(
            name=polymer_type,
            e_coeff=polymer_type.e_coeff,
            n_coeff=polymer_type.n_coeff,
            rho=polymer_type.Density,
            concentration_scalar=polymer_concentration,
            phi=self.phi,
        )
        self.polymer.initialize(grid_shape=grid_shape)

        # Initializing Surfactant Object
        self.surfactant = Surfactant(
            name=surfactant_type,
            initial_concentration=surfactant_concentration,
            IFT_equation=surfactant_type.IFT_equation,
            derivative_IFT_equation=surfactant_type.derivative_IFT_equation,
            phi=self.phi,
        )
        self.surfactant.initialize()

        # Initializing Water Object
        self.water = Water(
            init_water_saturation=SimulationConstants.Initial_Residual_Water_Saturation.value,
            init_aqueous_saturation=SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value,
            init_oleic_saturation=SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value,
            miuw=SimulationConstants.Water_Viscosity.value,
            miuo=SimulationConstants.Oil_Viscosity.value,
            phi=self.phi,
        )
        self.water.initialize(grid_shape=grid_shape)

        # Properties for Exporting Simulation Results
        self.COC = np.zeros((1, 2000))  # Cumulative Oil Recovered
        self.miuaTcal = np.zeros((1, 2000))  # Total Aqueous Viscosity
        self.lambdaTcal = np.zeros((1, 2000))  # Total Mobility (λ_o + λ_a)

        ##The following properties require memmaps:
        self.ProdRate, self.CROIP = (
            self._initialize_memmap_properties()
        )  # ProdRate (Production Rate) / CROIP (Cummulative Remaining Oil In Place)
        self.MFW = []
        self.integrated_inlet_flow = 0  # "src_total" in the MATLAB version of the code

    # Property of Simulation Class
    _grid_size = None
    @property
    def grid_size(self):
        """
        grid_size (float): the dimensions of the square grid
        """
        return self._grid_size
    @grid_size.setter
    def grid_size(self, value):
        self._grid_size = value

    _source_flow_magnitude = None
    @property
    def source_flow_magnitude(self):
        """
        source_flow_magnitude (float): flow rate at injection site
        """
        return self._source_flow_magnitude
    @source_flow_magnitude.setter
    def source_flow_magnitude(self, value):
        self._source_flow_magnitude = value

    _permeability_flag = None
    @property
    def permeability_flag(self):
        """
        permeability_flag (enum 'PermeabilityType'): sets the permeability field based on enum ``PermeabilityType``
        """
        return self._permeability_flag
    @permeability_flag.setter
    def permeability_flag(self, value):
        self._permeability_flag = value

    _reservoir_geometry = None
    @property
    def reservoir_geometry(self):
        """
        reservoir_geometry (enum 'ResevoirGeometry'): sets the reservoir geometry based on the enum ``ResevoirGeometry``
        """
        return self._reservoir_geometry
    @reservoir_geometry.setter
    def reservoir_geometry(self, value):
        self._reservoir_geometry = value

    _model_type = None
    @property
    def model_type(self):
        """
        model_type (enum 'ModelType'): sets the type of simulation being run based on the enum ``ModelType``
        """
        return self._model_type
    @model_type.setter
    def model_type(self, value):
        self._model_type = value

    _relative_permeability_formula = None
    @property
    def relative_permeability_formula(self):
        """
        relative_permeability_formula (enum 'RelativePermeabilityFormula'): sets the type of permeability formula being used in the simulation, based on the enum ``RelativePermeabilityFormula``
        """
        return self._relative_permeability_formula
    @relative_permeability_formula.setter
    def relative_permeability_formula(self, value):
        self._relative_permeability_formula = value

    _phi = None
    @property
    def phi(self):
        """
        phi (np.ndarray): porosity matrix
        """
        return self._phi
    @phi.setter
    def phi(self, value):
        self._phi = value

    _KK = None
    @property
    def KK(self):
       """
       KK (np.ndarray): the permeability matrix
       """
       return self._KK
    @KK.setter
    def KK(self, value):
       self._KK = value

    _time_step = None
    @property
    def time_step(self):
        """
        time_step (float): The Δt
        """
        return self._time_step
    @time_step.setter
    def time_step(self, value):
        self._time_step = value

    _polymer = None
    @property
    def polymer(self):
        """
        polymer (Polymer): Holds the ``Polymer`` object
        """
        return self._polymer
    @polymer.setter
    def polymer(self, value):
        self._polymer = value

    _surfactant = None
    @property
    def surfactant(self):
        """
        surfactant (Surfactant): Holds the ``Surfactant`` object
        """
        return self._surfactant
    @surfactant.setter
    def surfactant(self, value):
        self._surfactant = value

    _water = None
    @property
    def water(self):
        """
        water (Water): Holds the ``Water`` object
        """
        return self._water
    @water.setter
    def water(self, value):
        self._water = value

    _COC = None
    @property
    def COC(self):
        """
        COC (np.ndarray): An array that holds the cummulative oil captured
        """
        return self._COC
    @COC.setter
    def COC(self, value):
        self._COC = value

    _miuaTcal = None
    @property
    def miuaTcal(self):
        """
        miuTcal (np.ndarray): An array that caputres the change in the total aqueous viscosity over time
        """
        return self._miuaTcal
    @miuaTcal.setter
    def miuaTcal(self, value):
        self._miuaTcal = value

    _lambdaTcal = None
    @property
    def lambdaTcal(self):
        """
        lambdaTcal (np.ndarray): Array that holds the change in the total mobility (λ_T = λ_a + λ_o)
        """
        return self._lambdaTcal
    @lambdaTcal.setter
    def lambdaTcal(self, value):
        self._lambdaTcal = value

    _MFW = None
    @property
    def MFW(self):
        """
        MFW (np.ndarray): Array that holds change in the MFW (mean finger width)
        """
        return self._MFW
    @MFW.setter
    def MFW(self, value):
        self._MFW = value

    _integrated_inlet_flow = None
    @property
    def integrated_inlet_flow(self):
        """
        integrated_inlet_flow (float): Basically the integrating the source flow rate over time
        """
        return self._integrated_inlet_flow
    @integrated_inlet_flow.setter
    def integrated_inlet_flow(self, value):
        self._integrated_inlet_flow = value

    _source_prod_flow = None
    @property
    def source_prod_flow(self):
        """
        source_prod_flow (np.ndarray): The matrix with the source & and production well flow rates

        assuming that the source flow = production well flow (flow magnitudes are the same!)
        """
        # setting permeability state
        if self._source_prod_flow is None:
            self._source_prod_flow = np.zeros((self.mesh.n + 1, self.mesh.m + 1))
            bool_Homogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Homogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Quarter_Five_Spot = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (
                self.reservoir_geometry.value
                == ResevoirGeometry.Quarter_Five_Spot.value
            )
            if bool_Homogenous_and_Rectilinear or bool_Heterogenous_and_Rectilinear:
                self._source_prod_flow[:, 0] = (
                    self.source_flow_magnitude
                )  # intensity of injection well = src
                self._source_prod_flow[:, -1] = (
                    -1 * self.source_flow_magnitude
                )  # intensity of production well = -src
            elif bool_Heterogenous_and_Quarter_Five_Spot:  # Quarter-Five Spot
                self._source_prod_flow[0, 0] = (
                    self.source_flow_magnitude
                )  # Intensity of injection well = src
                self._source_prod_flow[-1, -1] = (
                    -1 * self.source_flow_magnitude
                )  # Intensity of production well = -src

        return self._source_prod_flow

    _scenario_flag = None
    @property
    def scenario_flag(self):
        """
        Determines the scenario based on the chosen reservoir geometry and permeability.
        Returns the integer value that represents a type of scenario run
        """
        if self._scenario_flag is None:
            bool_Homogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Homogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Rectilinear = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
            bool_Heterogenous_and_Quarter_Five_Spot = (
                self.permeability_flag.value == PermeabilityType.Heterogenous.value
            ) and (
                self.reservoir_geometry.value
                == ResevoirGeometry.Quarter_Five_Spot.value
            )

            if bool_Homogenous_and_Rectilinear:
                self._scenario_flag = 1
            elif bool_Heterogenous_and_Rectilinear:
                self._scenario_flag = 2
            elif bool_Heterogenous_and_Quarter_Five_Spot:
                self._scenario_flag = 3
            else:
                raise SimulationCalcInputException(
                    "SimulationCalcInputError:InvalidSimulationCase"
                )

        return self._scenario_flag

    # Private Methods of the Simulation Class
    def _initialize_pressure_and_velocity(self):
        """
        (private method)

        Initializing global pressure ('u') and velocity matrices ('v')
        Will use the ``n`` and ``m`` properties from ``Grid`` Class for initialization

        Returns: (tuple[np.ndarray, np.ndarray])
        ----------------------------------------
            The global pressure matrix (index 0) and velocity matrix (index 1)
        """
        u = np.zeros((self.mesh.n + 1, self.mesh.m + 1), dtype=np.complex128)
        v = np.zeros((self.mesh.n + 1, self.mesh.m + 1), dtype=np.complex128)

        return u, v

    def _compute_pressure_and_velocity_matrices(self, sparsed_A, B, beta):
        """
        (private method)

        dependent property to calculate the pressure matrix (``u``)
        and velocity matrix (``v``). Will rely on functions in the ``FEMesh`` class.
        
        Returns: (list[np.ndarray])
        ---------------------------
            List of update pressure matrix and velocity matrix => [u,v]
        """
        max_iterations = 1000
        new_u, convergence_flag = bicgstab(
            sparsed_A, B, rtol=1e-10, atol=0, maxiter=max_iterations
        )  # new_u is of shape (900, )
        if convergence_flag != 0:
            import warnings
            warnings.warn(f"BiCGSTAB: convergence issue (info={convergence_flag}) in pressure solver")

        #setting default value to out
        out = 0
        
        #homogenous
        if ( self.permeability_type == PermeabilityType.Homogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            out = y - init_front_hs + 0.01 * (np.cos(80 * np.pi * x))
        elif ( self.permeability_type == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            out = y - init_front_hs ## Rectilinear Homogenous
        elif ( self.permeability_type == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Quarter_Five_Spot ):
            out = ( (x)**2 ) + ( (y)**2 ) - 0.015 # Normal unperturbed initial saturation front 
        
        return out

        px, py = self._get_gradient(new_v)
        new_u = (-1 * beta) * px
        new_v = (-1 * beta) * py

        return new_u, new_v

    def _get_gradient(self, vn):
        """
        (private method)

        Helper function to determine the gradients with respect to x and y dimensions
        
        Returns: (tuple[_Array[tuple[int, int], float64], NDArray[float64]])
        --------------------------------------------------------------------
            Tuple with px py which are numpy matrices that hold the gradient wrt to x and y dimensions
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

            if (self.model_type == ModelType.No_Shear_Thinning):
                # Newtonian Model (NO SHEAR THINNING INVOLVED => MODEL TYPE #1)
                print("shape of polymer vec concentration:", np.shape(polymer_obj.vec_concentration))
                n = np.shape(polymer_obj.vec_concentration)[0]
                m = np.shape(polymer_obj.vec_concentration)[1]
                if polymer_obj.vec_concentration.all() == 0:
                    self.aqueous_viscosity = vis_water * np.ones((n, m))
                    print('aqueous viscosity: ', self.aqueous_viscosity)
                else:
                    self.aqueous_viscosity = vis_water * (1 + beta1 * polymer_obj.vec_concentration)
            elif (self.model_type == ModelType.Sourav_Implementation):
                # Sourav's Implementation (MODEL TYPE #2)
                n = np.shape((polymer_obj.vec_concentration,1))
                if(polymer_obj.initial_concentration == 0):
                    self.aqueous_viscosity = vis_water * np.ones(n)
                else:
                    self.aqueous_viscosity = vis_oil * (0.5 + polymer_obj.vec_concentration)
            elif(self.model_type == ModelType.Shear_Thinning_On):
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
        if(self.permeability_type == PermeabilityType.Homogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            # Represents a homogenous rectilinear model
            Kmax = 1000
            KK = Kmax*np.ones(self.sog+1)
        elif(self.permeability_type == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Rectilinear):
            Kmax = 100
            KK = Kmax*( 0.5*(1-10^(-7))*(np.sin(6*np.pi*np.cos(x))*np.cos(4*np.pi*np.sin(3*y))-1)+1)
        elif(self.permeability_type == PermeabilityType.Heterogenous and self.resevoir_geometry == ResevoirGeometry.Quarter_Five_Spot):
            # need to load the KK30Tabert.mat file... need to use the scipy.io.loadmat() method
            [Kmax, KK] = sp.io.loadmat('./Resources/KK30Tabert.mat')

        return [Kmax, KK]



    def get_gradient(self, vn):
        m = self.mesh.m
        n = self.mesh.n

        dx = self.mesh.dx
        dy = self.mesh.dy

        px = np.zeros((n + 1, m + 1))
        py = np.copy(px)

        for i in range(m + 1):
            for j in range(n + 1):
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

        return px, py

    def _initialize_memmap_properties(self):
        """
        (private method)

        Will initialize the ``ProdRate`` and ``CROIP`` properties.
        Using memmaps to allow window's users to run program.
        
        Returns: (tuple[np.ndarray, np.ndarray])
        ----------------------------------------
            Initialized ``ProdRate`` and ``CROIP`` properties
        """
        os.makedirs("memmaps", exist_ok=True)
        tf = 500
        dt = self.mesh.dx / self.source_flow_magnitude
        timestamps = int(np.floor(tf / dt))
        CROIP = np.memmap(
            "memmaps/CROIP.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )
        ProdRate = np.memmap(
            "memmaps/ProdRate.dat", dtype="float64", mode="w+", shape=(1, timestamps)
        )

        return ProdRate, CROIP

    def _create_mesh(self):
        """
        (private method)

        Returns: (Tuple[Grid, FEMesh])
        -----------------------------
            Initialized FD and FE mesh
        """
        FD_mesh = Grid(self.grid_size, self.grid_size)

        FE_mesh = FEMesh(self.grid_size, self.grid_size)

        return FD_mesh, FE_mesh

    def _initialize_simulation(self):
        """
        (private method)

        Sets up initial reservoir fields, permeability, and time step.

        Returns: (None)
        ---------------
            Initialized properties of simulation. Required in ``__init__`` function
        """
        self.phi = self._compute_phi()

        # Initialize permeability matrix
        self.KK = self._compute_permeability()

        # Calculate time step
        self.time_step = self.mesh.dx / self.source_flow_magnitude
        if (
            self.permeability_flag == PermeabilityType.Heterogenous
            and self.reservoir_geometry == ResevoirGeometry.Quarter_Five_Spot
        ):
            self.time_step *= 100

    def _compute_phi(self):
        """
        (private method)

        Returns: (np.ndarray)
        ---------------------
            Computes and returns the level set function phi at each grid point. 
            (Equivalent to MATLAB get_phi_test function.)
        """
        m = self.mesh.m
        n = self.mesh.n
        dx = self.mesh.dx
        dy = self.mesh.dy
        left = self.mesh.left
        bottom = self.mesh.bottom

        jj, ii = np.meshgrid(np.arange(1, n + 2), np.arange(1, m + 2))
        x_coords = left + (ii - 1) * dx
        y_coords = bottom + (jj - 1) * dy

        phi = self._z_func_test(x_coords, y_coords)
        return phi

    def _z_func_test(self, x, y):
        """
        (private method)

        Args:
        -----
            x (np.ndarry): x-dimension coordinates 

            y (np.ndarray): y-dimension coordinate points
        
        Returns: 
        --------
            Compute and returns the initial position of the water front.
            (Equivalent to MATLAB z_func_test.)

        """
        init_front_hs = 0.1
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        if bool_Homogenous_and_Rectilinear:
            # Homogeneous
            out = y - init_front_hs + 0.01 * np.cos(80 * np.pi * x)
        elif bool_Heterogenous_and_Rectilinear:
            # Rectilinear Heterogeneous
            out = y - init_front_hs
        elif bool_Heterogenous_and_Quarter_Five_Spot:
            # Quarter five spot
            out = np.square(x) + np.square(y) - 0.015
        else:
            out = np.zeros_like(x)

        return out

    def _compute_permeability(self):
        """
        (private method)
        
        Returns:
        -------
            Compute and returns the permeability matrix KK based on the ``scenario_flag``.
            (Equivalent to MATLAB KKdef function.)
        """
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        m = self.mesh.m
        KK = None
        if bool_Homogenous_and_Rectilinear:
            # Homogeneous
            Kmax = 100
            KK = Kmax * np.ones((m + 1, m + 1))
        elif bool_Heterogenous_and_Rectilinear:
            # Heterogeneous
            Kmax = 100
            KK = Kmax * (
                0.5
                * (1 - 10 ** (-7))
                * (
                    np.sin(6 * np.pi * np.cos(self.x))
                    * np.cos(4 * np.pi * np.sin(3 * self.y))
                    - 1
                )
                + 1
            )
        elif bool_Heterogenous_and_Quarter_Five_Spot:
            # Load Upper Ness formation (SPE10)
            mat_data = loadmat(
                "./lib/Resources/KK30Ness.mat"
            )  # FIXME: when using master_surf_grid need to change this path
            if "KK" not in mat_data:
                raise SimulationCalcInputException(
                    "SimulationInputException: KK matrix not found in KK30Ness.mat file."
                )
            KK = mat_data["KK"]
        return KK

    def _transport_equation_solver(self, dt):
        """
        (Private method)

        This method of the ``Simulation`` class will solve the transport equations
        for the surfactant concentration, polymer concentration, and water saturation

        Based on the 'nmmoc_surf_mod_neumann.m' function within the MATLAB version

        Note:
            dg - derivative with respect to surfactant concentration
            ds - derivative with respect to water saturation
            dc - derivative with respect to polymer concentration

        Returns: (tuple[float, float, float])
        ----------------------------------------------------
            tuple[ocut, wcut, ROIP]
            ocut - volume of oil in the production well
            wcut - volume of water in the production well
            ROIP - residual oil in place (as a volume fraction)
        """
        # Initialize constant parameters
        const_parameters = {}
        const_parameters["inlet_total_flow"] = self.source_flow_magnitude  # G1
        const_parameters["inlet_polymer_flow"] = (
            self.polymer.concetration_scalar * self.source_flow_magnitude
        )  # G2
        const_parameters["inlet_surfactant_flow"] = (
            self.surfactant.concentration * self.source_flow_magnitude
        )  # G3
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UndefinedWaterSaturationMatrix"
        )
        n = np.shape(self.water.water_saturation)[0]
        m = np.shape(self.water.water_saturation)[1]
        const_parameters["FD_grid_constants"] = {
            "n": n,
            "m": m,
            "dx": self.mesh.dx,
            "dy": self.mesh.dy,
            "dt": dt,
            "dt_matrix": dt * np.ones((n, m)),
            "x": self.mesh.x,
            "y": self.mesh.y,
        }
        ## parameters around Pc calculations
        const_parameters["Pc_constants"] = {
            "omega1": SimulationConstants.Capillary_Pressure_Param_1.value,
            "omega2": SimulationConstants.Capillary_Pressure_Param_2.value,
        }
        ## parameters around capillary number for residual saturation calcs
        const_parameters["resid_saturation_constants"] = {
            "Nco0": 10 ** (-5),
            "Nca0": 10 ** (-5),
        }
        ## porosity parameter
        const_parameters["porosity"] = 1

        ## Permeability
        const_parameters["KK"] = self.KK
        const_parameters["relative_permeability_formula"] = (
            self.relative_permeability_formula
        )

        # Initialize variable parameters
        varying_parameters = {}
        ## residual water and oil saturations
        [swr, sor] = self.water.compute_residual_saturations(
            sigma=self.surfactant.eval_IFT, u=self.u, v=self.v
        )
        varying_parameters["swr"] = swr
        varying_parameters["sor"] = sor
        nsw = (self.water.water_saturation - swr) / (1 - swr)
        nso = (self.water.water_saturation - swr) / (1 - swr - sor)
        varying_parameters["nsw"] = nsw
        varying_parameters["nso"] = nso
        ## recomputing mobilities
        assert (
            self.polymer.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPolymerConcentrationMatrix"
        )
        lambda_a = self.water.compute_mobility(
            c=self.polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=True,
            rel_permeability_formula=self.relative_permeability_formula,
        )
        lambda_o = self.water.compute_mobility(
            c=self.polymer.concentration_matrix,
            sor=float(sor),
            swr=float(swr),
            aqueous=False,
            rel_permeability_formula=self.relative_permeability_formula,
        )
        lambda_total = lambda_a + lambda_o
        varying_parameters["mobility_parameters"] = {
            "lambda_a": lambda_a,
            "lambda_o": lambda_o,
            "lambda_total": lambda_total,
        }
        ## fractional flow calculations
        f = lambda_a / lambda_total
        assert self.KK is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPermeabilityTensor"
        )
        D = self.KK * lambda_o * f
        varying_parameters["fractional_flow_parameters"] = {"f": f, "D": D}
        ## calculating dσ/dΓ
        varying_parameters["IFT_and_derivative"] = {
            "IFT": self.surfactant.eval_IFT,
            "dIFT_dg": self.surfactant.eval_dIFT_dGamma,
        }
        ## compute capillary number
        assert self.water.viscosity_array is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownAqueousViscosityMatrix"
        )
        nca = (
            np.sqrt(np.matmul(self.u, self.u) + np.matmul(self.v, self.v), dtype=np.complex128)
            * self.water.viscosity_array.astype(np.complex128)
            / self.surfactant.eval_IFT.astype(np.complex128)
        )
        nco = (
            np.sqrt(np.matmul(self.u, self.u) + np.matmul(self.v, self.v), dtype=np.complex128)
            * SimulationConstants.Oil_Viscosity.value
            / self.surfactant.eval_IFT
        )
        norm_nca = np.linalg.norm(nca)  # 2-norm of nca
        norm_nco = np.linalg.norm(nco)  # 2-norm of nco
        ## compute derivatives of residual saturations with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        dswr_dg = np.zeros((n, m))
        dsor_dg = np.zeros((n, m))
        swr_0 = SimulationConstants.Resid_Aqueous_Phase_Saturation_Initial.value
        sor_0 = SimulationConstants.Resid_Oleic_Phase_Saturation_Initial.value
        assert (
            self.surfactant.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownSurfactantConcentrationMatrix"
        )
        for j in range(n):
            for i in range(m):
                if norm_nca >= const_parameters["resid_saturation_constants"]["Nca0"]:
                    dswr_dg[j, i] = -(
                        swr_0
                        * 0.1534
                        * 10.001
                        * (
                            const_parameters["resid_saturation_constants"]["Nca0"]
                            ** 0.1534
                        )
                    ) / (
                        (
                            np.sqrt((self.u[j, i] ** 2) + (self.v[j, i] ** 2))
                            * self.water.viscosity_array[j, i]
                        )
                        ** (0.1534)
                        * self.surfactant.eval_IFT[j, i] ** (0.8466)
                        * (self.surfactant.concentration_matrix[j, i] + 1) ** 2
                    )
                if norm_nco >= const_parameters["resid_saturation_constants"]["Nco0"]:
                    dsor_dg[j, i] = -(
                        sor_0
                        * 0.5213
                        * 10.001
                        * (
                            const_parameters["resid_saturation_constants"]["Nca0"]
                            ** 0.5213
                        )
                    ) / (
                        (
                            np.sqrt((self.u[j, i] ** 2) + (self.v[j, i] ** 2))
                            * self.water.viscosity_array[j, i]
                        )
                        ** (0.5213)
                        * self.surfactant.eval_IFT[j, i] ** (0.4787)
                        * (self.surfactant.concentration_matrix[j, i] + 1) ** 2
                    )
        varying_parameters["resid_saturation_derivatives"] = {
            "dswr_dg": dswr_dg,
            "dsor_dg": dsor_dg,
        }
        ## compute derivatives of normalized saturation with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["normalized_saturation_derivatives"] = {
            "dnsw_dg": dswr_dg * (self.water.water_saturation - 1) / (1 - swr) ** 2,
            "dnso_dg": (
                dswr_dg * (self.water.water_saturation + sor - 1)
                + dsor_dg * (self.water.water_saturation - swr)
            )
            / (1 - swr - sor) ** 2,
        }
        ## computing relative permeability with respect to surfactant concentration (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["relative_permeability_derivatives"] = {
            "dkra_dg": 2.5 * dswr_dg * (nsw**3 - nsw)
            + (self.water.water_saturation - 1)
            * (2.5 * swr * (3 * nsw**2 - 1) + 1)
            * varying_parameters["normalized_saturation_derivatives"]["dnsw_dg"]
            / (1 - swr) ** 2,
            "dkro_dg": 1
            - 5 * sor * nso
            + (1 - nso) * (1 - 5 * nso * dsor_dg)
            - (1 + 5 * sor - 10 * sor * nso)
            * varying_parameters["normalized_saturation_derivatives"]["dnso_dg"],
            "dkra_ds": 2.5 * swr * (3 * (nsw) ** 2 - 1) + 1,
            "dkro_ds": 10 * sor * nso - 5 * sor - 1,
        }
        ## computing capillary pressure derivatives with respect to concentrations and saturations (FIXME: Need to update when inplementing autodiff!)
        pc = (
            self.surfactant.eval_IFT
            * const_parameters["Pc_constants"]["omega2"]
            * np.sqrt(const_parameters["porosity"])
        ) / (
            np.matmul(
                self.KK ** (0.5),
                fractional_matrix_power(
                    1 - nso, 1 / const_parameters["Pc_constants"]["omega1"]
                ),
            )
        )
        dpc_ds = pc / (const_parameters["Pc_constants"]["omega1"] * (1 - nso))
        dpc_dg = (
            pc / self.surfactant.eval_IFT
        ) * self.surfactant.eval_dIFT_dGamma + dpc_ds
        varying_parameters["capillary_pressure_and_derivatives"] = {
            "pc": pc,
            "dpc_ds": dpc_ds,
            "dpc_dg": dpc_dg,
        }
        ## computing fractional flow derivatives with respect to concentrations and saturations (FIXME: Need to update when inplementing autodiff!)
        varying_parameters["fractional_flow_derivatives"] = {
            "df_ds": varying_parameters["relative_permeability_derivatives"]["dkra_ds"]
            * lambda_o
            / (lambda_total**2 * self.water.viscosity_array)
            - varying_parameters["relative_permeability_derivatives"]["dkro_ds"]
            * lambda_a
            / (lambda_total**2 * self.water.miuo),
            "df_dc": (-1 * (lambda_o * lambda_a * self.water.miuo))
            / ((lambda_total**2) * self.water.viscosity_array),
            "df_dg": (
                (
                    varying_parameters["relative_permeability_derivatives"]["dkra_dg"]
                    * lambda_o
                )
                / ((lambda_total**2) * self.water.viscosity_array)
            )
            - (
                (
                    varying_parameters["relative_permeability_derivatives"]["dkro_dg"]
                    * lambda_a
                )
                / ((lambda_total**2) * self.water.viscosity_array)
            ),
            "dD_dg": D
            * varying_parameters["capillary_pressure_and_derivatives"]["dpc_dg"],
            "dD_ds": D
            * varying_parameters["capillary_pressure_and_derivatives"]["dpc_ds"],
        }
        # Update Water Saturation Matrix
        ## Calculate ``xmod`` and ``ymod``
        [xmod, ymod] = self._characteristic_coordinates(
            1,
            self.water.water_saturation,
            self.water.water_saturation,
            const_parameters,
            varying_parameters,
        )

        ## Pass in parameters into ``compute_water_saturation`` method of the ``Water`` class
        Q_old = np.copy(self.water.water_saturation)
        Qmod, varying_parameters = self.water.compute_water_saturation(
            grid=self.mesh,
            surfactant=self.surfactant,
            polymer=self.polymer,
            u=self.u,
            v=self.v,
            xmod=xmod,
            ymod=ymod,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
        )

        # Update the Polymer Concentration Matrix
        [xmod, ymod] = self._characteristic_coordinates(
            2, Q_old, self.water.water_saturation, const_parameters, varying_parameters
        )
        C_old = np.copy(self.polymer.concentration_matrix)
        varying_parameters = self.polymer.compute_concentration(
            grid=self.mesh,
            water_sat=self.water.water_saturation,
            u=self.u,
            v=self.v,
            xmod=xmod,
            ymod=ymod,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
        )

        # Update the Surfactant Concentration matrix
        [xmod, ymod] = self._characteristic_coordinates(
            3, Q_old, self.water.water_saturation, const_parameters, varying_parameters
        )
        G_old = np.copy(self.surfactant.concentration_matrix)
        G = np.copy(self.surfactant.concentration_matrix)
        x1d = self.mesh.x[0, :]
        y1d = self.mesh.y[:, 0]
        x_sorted = np.all(np.diff(x1d) > 0)
        y_sorted = np.all(np.diff(y1d) > 0)

        # reorder surfactant.concentration_matrix if a dimension isn't sorted
        if not x_sorted:
            x_sort_idx = np.argsort(x1d)
            x1d = x1d[x_sort_idx]
            G = G[:, x_sort_idx]  # Sort columns of surfactant.concentration_matrix
        if not y_sorted:
            y_sort_idx = np.argsort(y1d)
            y1d = y1d[y_sort_idx]
            G = G[y_sort_idx, :]  # Sort rows of surfactant.concentration_matrix

        interp = sp.interpolate.RegularGridInterpolator(
            (y1d, x1d),
            G,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )
        query_points = np.stack([ymod.ravel(), xmod.ravel()], axis=-1)
        Gmod = interp(query_points).reshape(xmod.shape)

        # Updating coefficients using interpolated surfactant concentration
        assert self.surfactant.IFT_conc_equ is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownIFTEquation"
        )
        assert (
            self.surfactant.derivative_IFT_conc_equ is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownDerivativeIFTEquation"
        )
        sigma_mod = self.surfactant.IFT_conc_equ(Gmod)
        sigma_g_mod = self.surfactant.derivative_IFT_conc_equ(Gmod)
        [swr, sor] = self.water.compute_residual_saturations(
            sigma=sigma_mod, u=self.u, v=self.v
        )
        lambda_a = self.water.compute_mobility(
            c=C_old,
            sor=float(sor),
            swr=float(swr),
            aqueous=True,
            rel_permeability_formula=self.relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_o = self.water.compute_mobility(
            c=C_old,
            sor=float(sor),
            swr=float(swr),
            aqueous=False,
            rel_permeability_formula=self.relative_permeability_formula,
            modified_water_saturation=Qmod,
        )
        lambda_total = lambda_a + lambda_o
        varying_parameters["mobility_parameters"] = {
            "lambda_a": lambda_a,
            "lambda_o": lambda_o,
            "lambda_total": lambda_total,
        }
        f = lambda_a / lambda_total
        assert self.KK is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownPermeabilityTensor"
        )
        D = self.KK * lambda_o * f
        varying_parameters["fractional_flow_parameters"] = {"f": f, "D": D}
        dpc_ds = pc / (const_parameters["Pc_constants"]["omega1"] * (1 - nso))
        dpc_dg = (pc / sigma_mod) * sigma_g_mod + dpc_ds
        varying_parameters["capillary_pressure_and_derivatives"] = {
            "pc": pc,
            "dpc_ds": dpc_ds,
            "dpc_dg": dpc_dg,
        }
        F = D * dpc_dg / self.water.water_saturation
        varying_parameters = self.surfactant.compute_concentration(
            grid=self.mesh,
            water_sat=self.water.water_saturation,
            const_parameters=const_parameters,
            varying_parameters=varying_parameters,
            F=F,
            Gmod=Gmod,
        )

        # calculate the Oil capture, water captured, and residual oil in place
        ocut = lambda_o[n-1, m-1] * self.source_flow_magnitude / lambda_total[n-1,m-1]
        wcut = lambda_a[n-1, m-1] * self.source_flow_magnitude / lambda_total[n-1,m-1]
        ROIP = 100*(np.sum(np.sum(1 - self.water.water_saturation)))/np.sum(np.ones((n*m,1)))
        return ocut, wcut, ROIP

    def _characteristic_coordinates(
        self,
        flag,
        old_water_saturation_matrix,
        new_water_saturation_matrix,
        const_parameters,
        varying_parameters,
    ):
        """
        (private method)

        Compute redefined characteristic coordinates (xmod, ymod) according to Neumann boundary conditions.

        will be a helper function to the ``self._transport_equation_solver()`` method.

        Args:
        ------
            flag (int): scenario flag for the simulation the user wants to run

            old_water_saturation_matrix (np.ndarray): water saturation matrix from previous iteration

            new_water_saturation_matrix (np.ndarray): water saturation from current iteration

            const_parameters (dict): constant parameters to help with calculations

            varying_parameters (dict): parameters that vary but assist with calculations for water saturation, polymer concentration, and surfactant concentration

        Returns: (tuple[np.ndarray, np.ndarray])
        ---------------------------------------
            This method returns ``xmod`` and ``ymod``, which are the modified characteristic coordinates according to the Neumann boundary conditions
        """
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:UnknownWaterSaturationMatrix"
        )
        assert (
            self.surfactant.concentration_matrix is not None
        ), SimulationCalcInputException(
            "SimulationCalcInputError:UnknownSurfactantConcentrationMatrix"
        )
        x, y = (
            const_parameters["FD_grid_constants"]["x"],
            const_parameters["FD_grid_constants"]["y"],
        )
        dt_matrix = const_parameters["FD_grid_constants"]["dt_matrix"]
        xjump = None
        yjump = None
        f = varying_parameters["fractional_flow_parameters"]["f"]
        f_s = varying_parameters["fractional_flow_derivatives"]["df_ds"]
        D = varying_parameters["fractional_flow_parameters"]["D"]
        pc_s = varying_parameters["capillary_pressure_and_derivatives"]["dpc_ds"]
        pc_g = varying_parameters["capillary_pressure_and_derivatives"]["dpc_dg"]
        sold = old_water_saturation_matrix
        snew = new_water_saturation_matrix

        if flag == 1:
            xjump = x - f_s * self.u * dt_matrix
            yjump = y - f_s * self.v * dt_matrix
        elif flag == 2:
            # Calculate gradients
            sx, sy = self._get_gradient(sold)
            gx, gy = self._get_gradient(self.surfactant.concentration_matrix)

            xjump = (
                x
                - (
                    (f / snew) * self.u
                    + (D * pc_s / snew) * sx
                    + (D * pc_g / snew) * gx
                )
                * dt_matrix
            )
            yjump = (
                y
                - (
                    (f / snew) * self.v
                    + (D * pc_s / snew) * sy
                    + (D * pc_g / snew) * gy
                )
                * dt_matrix
            )
        elif flag == 3:
            sx, sy = self._get_gradient(sold)

            xjump = x - ((f / snew) * self.u + (D * pc_s / snew) * sx) * dt_matrix
            yjump = y - ((f / snew) * self.v + (D * pc_s / snew) * sy) * dt_matrix

        # Apply Neumann reflection conditions
        if xjump is None or yjump is None:
            raise SimulationCalcInputException(
                "SimulationInputException:UnknownXJumpYJumpMatrices"
            )

        xmod = np.copy(x)
        ymod = np.copy(y)

        for j in range(np.shape(y)[0]):
            for i in range(np.shape(x)[1]):
                if xjump[j, i] <= 1 and yjump[j, i] <= 1:
                    xmod[j, i] = np.abs(xjump[j, i])
                    ymod[j, i] = np.abs(yjump[j, i])
                elif xjump[j, i] > 1 and yjump[j, i] <= 1:
                    xmod[j, i] = 2 - xjump[j, i]
                    ymod[j, i] = np.abs(yjump[j, i])
                elif xjump[j, i] <= 1 and yjump[j, i] > 1:
                    xmod[j, i] = np.abs(xjump[j, i])
                    ymod[j, i] = 2 - yjump[j, i]
                elif xjump[j, i] > 1 and yjump[j, i] > 1:
                    xmod[j, i] = 2 - xjump[j, i]
                    ymod[j, i] = 2 - yjump[j, i]

        return xmod, ymod

    def _export_results(self):
        """
        (private method)

        Saves simulation sim_results to CSV files.
        """
        import os

        os.makedirs("sim_results", exist_ok=True)

        np.savetxt(f"sim_results/COC_scenario_{self.scenario_flag}.csv", self.COC, delimiter=",")
        np.savetxt(f"sim_results/MFW_scenario_{self.scenario_flag}.csv", self.MFW, delimiter=",")
        np.savetxt(f"sim_results/CROIP_scenario_{self.scenario_flag}.csv", self.CROIP, delimiter=",")
        np.savetxt(f"sim_results/ProdRate_scenario_{self.scenario_flag}.csv", self.ProdRate, delimiter=",")
        # if hasattr(self, "lambdaTcal"):
        #     np.savetxt(
        #         "sim_results/lambdaTcal.csv", np.array(self.lambdaTcal), delimiter=","
        #     )
        # if hasattr(self, "miuaTcal"):
        #     np.savetxt(
        #         "sim_results/miuaTcal.csv", np.array(self.miuaTcal), delimiter=","
        #     )

        print("Simulation sim_results exported to /sim_results/ folder.")
        
    def _compute_MFW(self, UU):
        """
        (private function)

        Args:
        -----
            UU (np.ndarray): water saturation matrix

        Returns: (np.array)
        -------
            list of values which are the mean finger width during each iteration.
            Note: Only will run under the Rectilinear Homogenous and Rectilinear Heterogeneous simulation scenarios
        """
        # post processing of finger width 
        interface = np.zeros((29, 1))
        mean_UU_save = np.zeros((29, 29))
        store_UU = 0
        
        # find average concentration for each Y level
        iter_x = 0
        iter_y = 0
        iter_counter = 0
        
        while iter_y < 29:
            while iter_x < 29:
                
                if UU[iter_y, iter_x] > 0.21 and UU[iter_y, iter_x] < 0.99:
                    iter_counter += 1
                    store_UU += UU[iter_y, iter_x]
                    
                iter_x += 1
                
            mean_UU_save.T.flat[iter_y] = store_UU / iter_counter
            iter_counter = 0
            store_UU = 0
            iter_y += 1
            iter_x = 0
        
        iter_y = 0
        iter_x = 0
        while iter_y < 29:
            while iter_x < 29:
                if UU[iter_y, iter_x] < mean_UU_save.T.flat[iter_y]:
                    interface[iter_y, 0] = iter_x
                    break
                iter_x += 1
                
            iter_x = 0
            iter_y += 1
        
        # find location of interface front
        iter_x = 0
        iter_y = 0
        while iter_x < 29:
            if np.mean(UU[:, iter_x]) < np.mean(mean_UU_save[:, 0]):
                break
            iter_x += 1
        
        # find presence of saturation along the mixing layer
        rows = UU.shape[0]
        cols = UU.shape[1]
        check_concentration = np.zeros((rows, cols))
        for ii in range(rows):
            for jj in range(cols):
                if UU[ii, jj] < mean_UU_save[0, iter_y]:
                    check_concentration[ii, jj] = 1
                else:
                    check_concentration[ii, jj] = 0
        
        last = 0
        counter = 0
        mean_finger_width = 0
        total_concentration = 0
        jj = iter_x - 1
        for i in range(1):
            for ii in range(cols):
                if check_concentration[ii, jj] == 1:
                    new_last = 1
                    total_concentration += 1
                else:
                    new_last = 0
                
                if (new_last == 1 and last == 0) or (new_last == 0 and last == 1):
                    counter += 1
                last = new_last
        
            old_MFW = mean_finger_width
            mean_finger_width = 2 * total_concentration / (counter + 1)
            mean_finger_width = np.max(mean_finger_width, old_MFW)
            jj += 1
            
        return interface, mean_finger_width, iter_x

    # Public Method of Simulation Class
    def run(self):
        """
        Executes simulation loop.
        
        Raises:
        ------
            SimulationCalcInputException: If there is an issue with the execution of a calculation during runtime
        """
        assert self.water.water_saturation is not None, SimulationCalcInputException(
            "SimulationCalcInputError:WaterSaturationMatrixUnavailable"
        )
        bool_Homogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Homogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Rectilinear = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (self.reservoir_geometry.value == ResevoirGeometry.Rectilinear.value)
        bool_Heterogenous_and_Quarter_Five_Spot = (
            self.permeability_flag.value == PermeabilityType.Heterogenous.value
        ) and (
            self.reservoir_geometry.value == ResevoirGeometry.Quarter_Five_Spot.value
        )
        try:
            ## STEP 1: initializing start time, end time, and time step
            t = 0
            t_cal = 0
            t_stop = 500
            dt = self.mesh.dx / self.source_flow_magnitude
            # dt value if running under Quarter Five Spot & Heterogeneous scenario:
            if bool_Heterogenous_and_Quarter_Five_Spot:
                dt *= 100

            ## STEP 2: Initiating the primary 'while' loop that will keep running until water shows up in production well
            while(t < t_stop and self.water.water_saturation[self.mesh.n, self.mesh.m] <= 0.70):
                print(f'{t},{self.water.water_saturation[self.mesh.n, self.mesh.m]}')
                # while t < 1:
                ## STEP 2.1: Increment time and amount of feed used:
                self.integrated_inlet_flow += self.source_flow_magnitude
                t += dt
                ## STEP 2.2: Compute viscosities:
                if (
                    self.model_type.value == ModelType.No_Shear_Thinning.value
                ):  # if No Polymer Shear Thinning
                    self.water.compute_viscosity(
                        grid=self.mesh,
                        model_type=self.model_type,
                        polymer=self.polymer,
                        u=self.u,
                        v=self.v,
                    )
                    self.polymer.compute_viscosity(
                        grid=self.mesh,
                        u=self.u,
                        v=self.v,
                        model_type=self.model_type,
                        aqueous_viscosity=self.water.viscosity_array,
                    )
                elif self.model_type.value == ModelType.Shear_Thinning_On.value:
                    self.polymer.compute_viscosity(
                        grid=self.mesh,
                        u=self.u,
                        v=self.v,
                        model_type=self.model_type,
                        aqueous_viscosity=None,
                    )
                    self.water.compute_viscosity(
                        grid=self.mesh,
                        model_type=self.model_type,
                        polymer=self.polymer,
                        u=self.u,
                        v=self.v,
                    )
                ## STEP 2.2: Computing Residual Saturation:
                assert (
                    self.surfactant.IFT_conc_equ is not None
                ), SimulationCalcInputException(
                    "SimulationCalcInputError:SurfactantIFTEquationUnavailable"
                )
                interfacial_tension_matrix = self.surfactant.IFT_conc_equ(
                    self.surfactant.concentration_matrix
                )
                [resid_water_saturation, resid_oleic_saturation] = (
                    self.water.compute_residual_saturations(
                        sigma=interfacial_tension_matrix, u=self.u, v=self.v
                    )
                )
                ## STEP 2.3: Compute mobilities:
                assert (
                    self.polymer.concentration_matrix is not None
                ), SimulationCalcInputException(
                    "SimulationCalcInputError:PolymerConcentrationMatrixUnavailable"
                )
                aqueous_mobility = self.water.compute_mobility(
                    c=self.polymer.concentration_matrix,
                    sor=float(resid_oleic_saturation),
                    swr=float(resid_water_saturation),
                    aqueous=True,
                    rel_permeability_formula=self.relative_permeability_formula,
                )
                oleic_mobility = self.water.compute_mobility(
                    c=self.polymer.concentration_matrix,
                    sor=float(resid_oleic_saturation),
                    swr=float(resid_water_saturation),
                    aqueous=False,
                    rel_permeability_formula=self.relative_permeability_formula,
                )
                total_mobility = aqueous_mobility + oleic_mobility
                assert self.KK is not None, SimulationCalcInputException(
                    "SimulationCalcInputError:PermeabilityTensorUnavailable"
                )
                beta = self.KK * total_mobility

                ## STEP 2.4: Calculating Global Pressure and velocity
                ### STEP 2.4.1: setting FEM Mesh
                self.FE_mesh.set_triangulation()
                self.FE_mesh.set_FE_meshgrid(beta)
                self.FE_mesh.set_right_hand(self.source_prod_flow)
                self.FE_mesh.get_A_B_matrices()
                ### STEP 2.4.2: updating the pressure & velocity matrices
                u_old = self.u  # storing old pressure matrix
                v_old = self.v  # storing old velocity matrix
                self.u, self.v = self._compute_pressure_and_velocity_matrices(
                    self.FE_mesh.sparsed_A, self.FE_mesh.B, beta
                )

                ## STEP 2.5: Solving Transport Equations
                ocut, wcut, ROIP = self._transport_equation_solver(dt)
                # self._transport_equation_solver(dt)
                ## Step 2.6: MFW post processing (excluding QFS)
                if (self.scenario_flag != 3): # FIXME: compute_MFW currently operates for rectilinear geometries. Implement MFW computation for QFS
                    interface, MFW_val, _ = self._compute_MFW(self.water.water_saturation)
                    self.MFW.append(MFW_val)

                ## STEP 2.7: Updating the cummulative oil captured, Production rate, and the residual oil in place 
                # arrays for exporting to CSV files
                if (t_cal == 0):
                    self.COC[0,t_cal] = ocut
                else:
                    self.COC[0,t_cal] = self.COC[0,t_cal - 1] + ocut

                self.ProdRate[0,t_cal] = ocut/dt
                self.CROIP[0,t_cal] = ROIP

                t_cal += 1

            self._export_results()

        except Exception as e:
            print(e)
