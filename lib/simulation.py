"""
This python script contains the class definition for running simulations
"""

from .surfactant import Surfactant
from .polymer import Polymer
from . import Exceptions
from .para import Box

class Simulation:
    def __init__(self, sim_id, polymer, surfactant, resevoir_geometry, permeability_flg, mesh_grid):
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
        
        # PROTECTED MEMBERS
        self._polymer_ = None
        self._surfactant_ = None
        self._resevoir_geometry_ = None
        self._permeability_flag_ = None
        self._mesh_ = None
        
        #PUBLIC MEMBERS
        self.sim_id = sim_id
        self.polymer = polymer
        self.surfactant = surfactant
        self.resevoir_geometry = resevoir_geometry
        self.permeability_flg = permeability_flg
        self.mesh = mesh_grid

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

