"""
This file should hold the test cases for the 'simulation.py' class

@author: Bhargav Akula Ramesh Kumar
"""
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.para import Box
from lib.simulation import Simulation
from lib.polymer import Polymer
from lib.surfactant import Surfactant
from lib.enumerations import PolymerList, SurfactantList, ResevoirGeometry, PermeabilityType, PlotType, ModelType

#pytest import for unit tests
import pytest

def initializing_simulation():
    """
    This function will initialize the simulation object

    This test will generate a simulation object for homogenous porous media, 
    rectilinear resevoir geometry, with 0 wppm Surfactant and 
    0.001 wppm injection concentration of Xanthane.

    :return: Test simulation object
    :rtype: Simulation
    """
    
    #simulation id
    sim_id = 1

    #initializing polymer object
    polymer_object = Polymer(PolymerList.Xanthane, 0.001, PolymerList.Xanthane.e_coeff, PolymerList.Xanthane.n_coeff)
    
    #initializing surfactant object
    lambda_IFT_equation = lambda GG: 10.001 / (GG + 1) - 0.001
    surfactant_object = Surfactant(SurfactantList.Alkyl_Ether_Sulfate, 0, lambda_IFT_equation)

    #initializing initial water saturation
    initial_water_saturation = 0.79

    #initializing resevoir geometry
    resevoir_geometry = ResevoirGeometry.Rectilinear

    #initializing permeability
    permeability_type = PermeabilityType.Homogenous

    #initializing mesh_grid
    mesh_grid = Box()

    #initialiing whether surfactant will be present
    is_surfactant = True

    #plot type (Polymer, Surfactant, or Saturation Plots)
    plot_type = PlotType.Saturation_Plot

    #model type
    mdl_type = ModelType.No_Shear_Thinning 

    sim_obj_test = Simulation(
            sim_id= sim_id,
            polymer= polymer_object,
            surfactant= surfactant_object,
            init_water_saturation= initial_water_saturation,
            resevoir_geometry= resevoir_geometry,
            permeability_flg= permeability_type,
            mesh_grid= mesh_grid,
            mdl_id= mdl_type,
            plt_type= plot_type,
            is_surfactant= is_surfactant
            )
    
    return sim_obj_test


def test_get_phi_value():
    """
    This function will test calculation of phi value. This test also goes through testing the z_func

    """
    #initializing simulation object
    test_sim_object = initializing_simulation()
     
    #getting the phi value
    phi_vec = test_sim_object.get_phi_value()
    print("The phi vector is: ", phi_vec)


def test_initialize_concentration():
    """
    This function will be use test the initialization of surfactant and polymer concentrations
    at the injection site
    """
    #intializing simulation object
    test_sim_object = initializing_simulation()
    
    #determining intial concentration matrix
    concentration_matrix = test_sim_object.initial_concentration_matrix()

    print(concentration_matrix)

    pass

def test_compute_viscosity():
    pass

def test_compute_resid_saturations():
    pass

def test_compute_mobility():
    pass

def test_solving_saturation_equations():
    pass


if __name__ == "__main__":
    test_initialize_concentration()
