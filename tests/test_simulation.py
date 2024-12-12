"""
This file should hold the test cases for the 'simulation.py' class
"""

from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.enumerations import PolymerList

def test_initializing_simulation():
    """
    This pytest will initialize the simulation object

    This test will generate a simulation object for homogenous porous media, 
    rectilinear resevoir geometry, with 0 wppm Surfactant and 
    0.001 wppm injection concentration of Xanthane.
    """

    sim_id = 1
    polymer_object = Polymer(PolymerList.Xanthane, 0.001, PolymerList.Xanthane.e_coeff, PolymerList.Xanthane.n_coeff) 
    
    pass

def test_get_phi_value():
    pass

def test_z_func_calculation():
    pass

def test_initialize_concentration():
    pass

def test_compute_viscosity():
    pass

def test_compute_resid_saturations():
    pass

def test_compute_mobility():
    pass

def test_solving_saturation_equations():
    pass
