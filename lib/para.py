"""
This class hold the para and box classes which will be used to formulate the mesh grid

@author: Bhargav Akula Ramesh Kumar
"""

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))

from lib.enumerations import (
    SimulationConstants,
)

class Para:
    def __init__(self):
        """Initializing default values for generating mesh grid for simulation"""
        self.m = SimulationConstants.Grid_Size.value   # Number of points in x-direction
        self.n =  SimulationConstants.Grid_Size.value # Number of points in y-direction
        self.dx = None  # Spacing in x-direction (to be calculated later)
        self.dy = None  # Spacing in y-direction (to be calculated later)


class Box(Para):
    def __init__(self):
        super().__init__()  # Initialize the parent class (Para)
        self.left = 0  # x_min
        self.right = 1  # x_max
        self.bottom = 0  # y_min
        self.top = 1  # y_max
        self.calculate_spacing

    @property
    def calculate_spacing(self):
        # Calculate dx and dy based on the box dimensions
        self.dx = (self.right - self.left) / (self.m)
        self.dy = (self.top - self.bottom) / (self.n)
