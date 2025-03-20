from typing import List
import json

################################################################
# This utility handles user inputs and settings configuration. #
################################################################


def get_user_input(file_name: str) -> List[dict]:
    simulations: List[dict] = []
    with open(file_name, "r") as input_file:
        input_simulations = json.load(input_file)
        for simulation in input_simulations:

            simulations.append(simulation)
    return simulations
