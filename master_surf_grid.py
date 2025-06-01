"""
This is the main program for a full flooding simulation
    -- Grid sizes used 15x15, 20x20, 30x30, 40x40

This code was derived from EOR repository developed by Sourav Dutta, Rohit Mishra, and Prabir Daripa

Code Refactoring to the Python programming language was done by Bhargav Akula Ramesh Kumar and Carlos Acosta Caripo

@author: Bhargav Akula Ramesh Kumar, Carlos Acosta Caripo
"""

#### IMPORT STATEMENTS
import tkinter as tk

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
from user_input.gui import UserInputGUI
from lib.simulation import Simulation

def main() -> None:
    root = tk.Tk()
    app = UserInputGUI(root)
    root.mainloop()

    user_input = app.get_input()
    for index, simulation in enumerate(user_input):
        # Will need to pass the dictionary into the simulation class
        sim_object = Simulation(user_input_dict=simulation)
        sim_object.run()

if __name__ == "__main__":
    main()
