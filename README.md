# Description
Surfactant-Polymer (SP) flooding is an enhanced-oil-recovery (EOR) process which involves injecting a polymer and surfactant-laden aqueous phase into oil reservoirs to improve oil recovery through:
 - Reduction of interfacial tension (IFT) to decrease residual oil saturation
 - Increasing the viscosity of the displacing front thereby improving the overall sweep efficiency. 

This repository contains the code and data used in our research to model and analyze the effects of SP flooding, with a  particular focus on the shear-thinning behavior of polymers. 

## Purpose
This research aims to view the effects of surfactant in EOR while accounting for the shear-thinning behavior of polymers. 

## Key Features:
 - Simulation models for SP flooding with shear-thinning effects
 - Analysis of polymer and surfactant interactions in the EOR process
 - Data visualization and interpretation of cumulative oil recovery
 - Comparison of different polymer, geometric, and permeability field conditions

## Usage
 - Execute simulation scripts to explore various SP flooding scenarios allowing for the manipulation of various parameters such as:
    - Geometry
    - Porous Media
    - Surfactant Concentration
    - Initial Polymer Concentration (IPC)
    - Injection Rate (IR)
 - Analyze and visualize results to gain insights into SP flooding dynamics

## Requirements
 - ```Python 3.13.1``` (with dependencies listed in ```requirments.txt```)
 
### Installing Dependencies
1. cd to the location of the repository within your computer
2. Next create a virtual environment using the following command: ```python -m venv my_project_env```
3. Activate the virtual environment
    - **Windows**: ```.\my_project_env\Scripts\activate```
    - **MacOS/Linux**: ```source my_project_env/bin/activate```
4. Install the dependencies using the following command: ```pip install -r requirements.txt```
 

## Acknowledgments
### Principal Investigator:
 - [Professor Prabir Daripa](https://www.math.tamu.edu/directory/formalpg.php?user=daripa) - Texas A&M University, Department of Mathematics
 
### Previous Works: 
This research builds upon the following key works: 

[1] Prabir Daripa and R. Mishra, “Modeling shear thinning polymer flooding using a dynamic viscosity model,” Physics of Fluids, vol. 35, no. 4, Apr. 2023, doi: https://doi.org/10.1063/5.0145061.
‌

[2] Prabir Daripa and S. Dutta, “Modeling and simulation of surfactant–polymer flooding using a new hybrid method,” Journal of Computational Physics, vol. 335, pp. 249–282, Apr. 2017, doi: https://doi.org/10.1016/j.jcp.2017.01.038.
‌

[3] Prabir Daripa and S. Dutta, “On the convergence analysis of a hybrid numerical method for multicomponent transport in porous media,” Applied Numerical Mathematics, vol. 146, pp. 199–220, Dec. 2019, doi: https://doi.org/10.1016/j.apnum.2019.07.009.
‌

[4] P. Daripa and S. Dutta, "DFEM-MMOC based EOR code in MATLAB," GitHub repository, GitHub, 2020. [Online]. Available: https://github.com/daripa8371/EOR. 
