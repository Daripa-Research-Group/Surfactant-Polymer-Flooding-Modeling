"""
This python script contains the class definition for running simulations
"""

#Relevant imports
from surfactant import Surfactant
from polymer import Polymer
from Exceptions import OutOfRangeError, SimulationCalcInputException
from para import Box

import numpy as np



class Simulation:
    def __init__(self, sim_id : int, polymer : Polymer, surfactant : Surfactant, resevoir_geometry : str, permeability_flg : str, mesh_grid : Box):
        """
        creates instance of the simulation class which will enable for calculating changes in system parameters at every time-step

        :param sim_id: Simulation number
        :type sim_id: int

        :param polymer: Polymer object used in SP-flooding run
        :type polymer: Polymer

        :param surfactant: Surfactant object used in SP-flooding run
        :type surfactant: Surfactant

        :param resevoir_geometry: Type of resevoir geometry (is it a rectilinear or quarter-five-spot geometry)
        :type resevoir_geometry: str

        :param permeability_flg: Homogenous vs. Heterogenous porosity in resevoir
        :type permeability_flg: str

        :param mesh_grid: mesh_grid used in the SP-flooding run
        :type mesh_grid: Box
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

        #General Parameters in Simulation
        self._phi_ = None 
        self.sim_id = sim_id


    #PROPERTIES OF SIMULATION CLASS
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


    def z_func_test(self, x, y):
        """
        Specifying the initial position of the water front
        A function describing the initial position of 
        the water front in the shape of a circular arc 
        $$ z(x,y) = x^2+y^2-0.015 $$ 
        This can take array input

        :return: returns water front position
        """
        init_front_hs = 0.1
        # out=(x-0.15)^2+(y-0.15)^2 -0.5*(x+y-0.35)^2;
        # # perturbed initial saturation front for special fingering simulations
        # out = x.^2 + y.^2 - 0.015*(1+0.1*sin(18*atan(y./x)))^2;

        #setting default value to out
        out = 0
        
        #homogenous
        if ( self.permeability_flg == "HOMOGENOUS" and self.resevoir_geometry == "RECTILINEAR" ):
            out = y - init_front_hs + 0.01 * (np.cos(80 * np.pi * x))
        elif ( self.permeability_flg == "HETEROGENOUS" and self.resevoir_geometry == "RECTILINEAR" ):
            out = y - init_front_hs ## Rectilinear Homogenous
        elif ( self.permeability_flg == "HETEROGENOUS" and self.resevoir_geometry == "QUARTER_FIVE_SPOT" ):
            out = ( (x)**2 ) + ( (y)**2 ) - 0.015 # Normal unperturbed initial saturation front 
        
        return out
