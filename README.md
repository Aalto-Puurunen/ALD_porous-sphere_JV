# ALD_porous-sphere_JV
Version: 0.3

Diffusion-reaction model for ALD on porous spheres

# Project description
This Python script provides the solution of the one-dimensional diffusion equation with surface reaction for volumetric reactant number density and surface coverage as function of time along the radial coordinate of a porous sphere (Eq. 1 of Heikkinen et al., Phys. Chem. Chem. Phys., 2024, 26, 7580, DOI: [10.1039/d3cp05639b](https://doi.org/10.1039/D3CP05639B), for spheres). For the surface reaction, the model uses Langmuir adsorption and includes a desorption term (as in Eq. 16 of Ylilammi et al., J. Appl. Phys. 123, 205301, 2018, DOI: [10.1063/1.5028178](https://doi.org/10.1063/1.5028178)). The partial differential equations for diffusion and reaction are solved numerically. For the calculations, the effective diffusion coefficient is assumed to be constant along the porous sphere. This script can be used in a wide range of diffusion regimes (Kn number from Kn<<1 to Kn>>1). The script was written by Dr. Jorge A. Velasco, by request of Prof. Riikka L. Puurunen (Catalysis Group, Aalto University) with funding partially provided by the Genesis EU project (Chips JU, Horizon Europe), and the GreenAro project (BussinesFinland).  

# Usage
The parameters for the simulation are entered in the accompanying file “Parameters_ALD_porous-spheres.xlsx”. The simulation is performed by running “ALD_porous-spheres_JV_v03.py”. There’s no need to modify the Python script nor the name of the parameters file. Once the simulation ends, the script creates a simulation results file (the given running code/name in the parameters file is shown in the result’s file name: “given running code”_simulation_results_v03.xlsx). The simulation results file includes input parameters, calculated values, final profiles for volumetric reactant number density, reactant pressure, surface coverage, and fraction of coated volume along the radial distance of the sphere. Calculated values include: diffusion coefficients, Knudsen number, Thiele modulus, weight percent of deposited metal, penetration depth and slope at half coverage.     

# Citing 
Please cite as:
J.A. Velasco and R. L. Puurunen, ALD_porous-sphere_JV – Diffusion-reaction model for ALD on porous spheres, (2025), Github repository, [https://github.com/Aalto-Puurunen/ALD_porous-sphere_JV](https://github.com/Aalto-Puurunen/ALD_porous-sphere_JV).

# Copyright and license
MIT License
Copyright 2025 (c) Jorge A. Velasco and Riikka Puurunen, Aalto University

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
