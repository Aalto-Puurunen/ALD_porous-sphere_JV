"""
ALD_PorouSSphere_JV
Version: 0.3

Diffusion-reaction model for ALD on porous spheres. 

Solution of of the one-dimensional diffusion equation with surface reaction for volumetric reactant number density and 
surface coverage as function of time along the radial coordinate of a 
porous sphere (Eq. 1 of Heikkinen et al., Phys. Chem. Chem. Phys.,2024, 26, 7580, DOI: 10.1039/d3cp05639b, for spheres). 
For the surface reaction, the model uses Langmuir adsorption and includes a desorption term (as in Eq. 16 of Ylilammi et al., J. Appl. Phys. 123, 205301, 2018, DOI: 10.1063/1.5028178). 
The partial differential equations for diffusion and reaction are solved numerically.
For the calculations, the effective diffusion coefficient is assumed to be constant along the porous sphere. 
This script can be used in a wide range of diffusion regimes (Kn number from Kn<<1 to Kn>>1).

DO NOT modify this script. Use/modify parameters in the file “Parameters_ALD_porous-spheres.xlsx”

Created on August 19, 2025
@author: Jorge A. Velasco, by request of Prof. Riikka L. Puurunen (Catalysis Group, Aalto University). 
Partially funded by Genesis EU project (Chips JU, Horizon Europe), and the GreenAro project (BussinesFinland). 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd
from openpyxl import load_workbook
import os

# ---- Loading the Excel workbook ----
wb = load_workbook("Parameters_ALD_porous-spheres.xlsx", data_only=True)  # data_only=True reads the value, not formula
# Selecting the active workshee
ws = wb["Parameters"]
# Reading parameter values from specific cells
Runcode = ws["B1"].value
Rp = ws["B3"].value     # radius of the sphere (m)
Por = ws["B4"].value    # Porosity of the sphere(-)
Tau = ws["B5"].value    # Tortuosity (-)
dpore = ws["B6"].value  # Pore diameter (m)
SA = ws["B7"].value     # Specific surface area (m2/g)
VP = ws["B8"].value     # Pore volume (m3/g)
DenS = ws["B9"].value   # Skeletal density of the porous sphere (kg/m3)
c = ws["B10"].value     # Initial sticking probability (-)
q = ws["B11"].value     # Adsorption capacity (m-2) 
Pd = ws["B12"].value    # Desorption probability (s-1)
t_pulse = ws["B13"].value  # Pulse time (s) 
T = ws["B14"].value     # Temperature (K) 
pA0 = ws["B15"].value   # Partial pressure of reactant at z=0 (Pa) 
pI = ws["B16"].value    # Partial pressure of inert gas (Pa)
MA = ws["B17"].value    # Molar mass of reactant(kg/mol) 
MI = ws["B18"].value    # Molar mass of inert gas (kg/mol)
dA = ws["B19"].value    # Diameter of reactant (m) 
dI = ws["B20"].value    # Diameter of inert gas (m)
MM = ws["B21"].value    # Atomic mass of metal (kg/mol)
MMOx = ws["B22"].value  # Molar mass of metal oxide (kg/mol)
N = ws["B23"].value     # number of spatial grid points in r axis
Ntime = ws["B24"].value # Number of time grid points
ScriptN = os.path.basename(__file__)

R = 8.31446  # Gas constant J/(mol*K); 
NAv = 6.022E23 # Avogadro's number (mol-1);   
kb = 1.38064852E-23; # Boltzmann constant (m2kg/s2K)

nA0 = NAv*pA0/(R*T)

# ---- Calculated parameters using set values ----
vA = np.sqrt((8*R*T)/(np.pi*MA)) # Average speed of A molecule (m/s)
DKn = (1/3)*vA*dpore             #
zA = np.pi/4*(dA+dI)**2*((8*R*T/np.pi*(1/MA + 1/MI))**(1/2))*pI*NAv/(R*T)+np.pi*(dA)**2*((16*R*T/(np.pi*MA))**(1/2))*pA0*NAv/(R*T)
DA = 3*np.pi*(vA**2)/(16*zA)  # Gas-phase diffusion constant of A molecules (m2/s)
D1 = (1/(1/DA + 1/DKn)) 
D = (Por/Tau)*D1 # Based on Evans in https://pubs.acs.org/doi/pdf/10.1021/ie00055a004
s_AP = SA/VP    # ratio of specific surface area and pore volume (m-1)

# Knudsen number calculations
sAA = np.pi*((dA/2)+(dA/2))**2
sAB = np.pi*((dA/2)+(dI/2))**2
term1 = (2**0.5)*sAA*pA0
term2 = ((1 + (MA/MI))**0.5)*pI*sAB
MFP = kb*T/(term1 + term2)
Knudsen = MFP/dpore
# Thiele modulus calculations
Alph = (Rp**2)*c*vA/(dpore*D)
Thiele = (Alph)**0.5

# Rate constants 
k1 = (1/4)*s_AP*vA*c       # rate constant for nA
k2 = (1/4)*(1/q)*vA*c      # rate constant for theta

# ---- Solving PDE's ----
# Spatial discretization
r = np.linspace(0, Rp, N)      # radial positions
dr = r[1] - r[0]

# Initial conditions
nA_init = np.zeros(N)          # nA(r, 0)
theta_init = np.zeros(N)       # theta(r, 0)
y0 = np.concatenate([nA_init, theta_init])

# Function to compute derivatives
def system(t, y):
    nA = y[:N]
    theta = y[N:]
    
    dnA_dt = np.zeros(N)
    dtheta_dt = np.zeros(N)

    # Second derivative with spherical term
    for i in range(1, N-1):
        d2nA = (nA[i+1] - 2*nA[i] + nA[i-1]) / dr**2
        dnA_dr = (nA[i+1] - nA[i-1]) / (2*dr)
        dnA_dt[i] = D * (d2nA + (2/r[i]) * dnA_dr) - k1 * nA[i] * (1 - theta[i])

    # Center (r = 0): use symmetric forward difference
    d2nA_0 = 2 * (nA[1] - nA[0]) / dr**2
    #dnA_dr_0 = (nA[1] - nA[0]) / dr
    dnA_dt[0] = D*d2nA_0 - k1*nA[0]*(1 - theta[0])

    # Surface (r = Rp): Dirichlet BC
    nA[-1] = nA0
    dnA_dt[-1] = 0  # hold boundary fixed

    # Theta equation (kinetics including a desorption term)
    dtheta_dt = k2 * nA * (1 - theta) - (theta * Pd)

    return np.concatenate([dnA_dt, dtheta_dt])

# Time span
t_span = (0, t_pulse)
t_eval = np.linspace(*t_span, Ntime)

# Solve the PDEs
sol = solve_ivp(system, t_span, y0, t_eval=t_eval, method='BDF')

# Extract results
nA_sol = sol.y[:N, :]
theta_sol = sol.y[N:, :]

# ---- Estimation of metal loading-----
# fraction of coated volume as function of "r"
FracVol = 1-(r/Rp)**3
from scipy.integrate import simpson
thetafinal = theta_sol[:, -1]
Frac = simpson(thetafinal[::-1], x=FracVol[::-1])
Vsphere = (4/3)*np.pi*Rp**3 # Volume of sphere (m3)
msphere = (1-Por)*Vsphere*DenS # Mass of sphere (kg)
Metalmass = Frac*(q*MM*(SA*1000)/NAv)* msphere # 
MOxmass = Metalmass*MMOx/MM
wtper = Metalmass*100/(msphere+MOxmass)
wtfull = (q*MM*(SA*1000)/NAv)*msphere*100/(msphere+MOxmass)
#print(Frac, wtfull, wtper) 

# ---- Extracting penetration depth at theta=0.5 and slope ---- 
# Find indices where Theta crosses 0.5
crossing_indices = np.where(np.diff(np.sign(thetafinal - 0.5)))[0]
if len(crossing_indices) == 0:
    print("No crossing found near Theta = 0.5")
    x_interp = 0
    Theta2 = 1
    Theta1 = 0
    x2=1
    x1=0
else:
    # Use the first crossing 
    i = crossing_indices[0]
    # Surrounding values
    Theta1, x1 = thetafinal[i], r[i]
    Theta2, x2 = thetafinal[i + 1], r[i + 1]
    # Linear interpolation for x at Theta = 0.5
    x_interp = x1 + (0.5 - Theta1) * (x2 - x1) / (Theta2 - Theta1)
    #print(f"Theta1 = {Theta1}, x1 = {x1}")
    #print(f"Theta2 = {Theta2}, x2 = {x2}")
    #print(f"Interpolated x at Theta = 0.5: x_interp = {x_interp}")
PD50 = (Rp-x_interp)  # in (m)
PD50rel = (Rp-x_interp)/Rp
slopePD50rel = abs(Rp * (Theta2 - Theta1) / (x2 - x1))
ThetaDiff = abs(Theta1 - Theta2)*100

#------- Temporal evolution profile plots -----------
# Plot of evolution of partial pressure, pA
X, Tgrid = np.meshgrid((Rp-r)*1e6, sol.t) # Rp-r values in micrometers
plt.figure(figsize=(10, 5))
plt.contourf(X, Tgrid,(R*T/NAv)*nA_sol.T, levels=50, cmap='viridis')
plt.colorbar(label='Pressure pA(r,t)')
plt.xlabel('d=Rp-r (μm)')
plt.ylabel('Time (s)')
plt.title('Temporal evolution of partial pressure of A')
plt.show()

# Plot evolution of surface coverage
plt.figure(figsize=(10, 5))
plt.contourf(X, Tgrid, theta_sol.T, levels=50, cmap='inferno')
plt.colorbar(label='Surface coverage, theta (x,t)')
plt.xlabel('d=Rp-r (μm)')
plt.ylabel('Time (s)')
plt.title('Temporal evolution of surface coverage')
plt.show()

#---- Final profile plots ----

# Plot final profiles
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot((Rp-r)*1e6, nA_sol[:, -1], label='nA at t_final')
plt.xlabel('d=Rp-r (µm)')
plt.ylabel('nA (1/m^3)')
plt.title('Final nA profile')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot((Rp-r)*1e6, theta_sol[:, -1], label='theta at t_final', color='orange')
plt.xlabel('d=Rp-r (µm)')
plt.ylabel('theta (-)')
plt.title('Final θ profile')
plt.grid(True)

# Plot final profiles
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot((Rp-r)*1e6, (R*T/NAv)*nA_sol[:, -1], label='PA at t_final')
plt.xlabel('d=Rp-r (µm)')
plt.ylabel('pA (Pa)')
plt.title('Final pA profile')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot((Rp-r)*1e6, theta_sol[:, -1], label='theta at t_final', color='orange')
plt.xlabel('d=Rp-r (µm)')
plt.ylabel('theta (-)')
plt.title('Final θ profile')
plt.grid(True)

# Plot final profiles
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(FracVol, (R*T/NAv)*nA_sol[:, -1], label='PA at t_final')
plt.xlabel('Coated Volume / Total Volume  (-)')
plt.ylabel('pA (Pa)')
plt.title('Final pA profile')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(FracVol, theta_sol[:, -1], label='theta at t_final', color='orange')
plt.xlabel('Coated Volume / Total Volume')
plt.ylabel('theta (-)')
plt.title('Final θ profile')
plt.grid(True)

plt.tight_layout()
plt.show()

# --- Creating a dataframe for exporting---
# Input parameters 
param_data = {
    "Parameter": ["Running Code", "Script's name & version:", "Radius of the sphere, R_p (m)", "Porosity of the sphere, ε (-)", "Tortuosity of the sphere, τ (-)", "Pore diameter d_pore (m)", "Specific surface area, S (m2/g)", "Pore volume, V_pore (m3/g)", "Skeletal density of sphere, ρ_S (kg/m3)", "Sticking coeficient, c (-)", "Adsorption Capacity, q (m-2)", "Desorption probability, P_d (1/s)", "Pulse time, t_pulse (s)", "Temperature, T (K)", "Partial pressure of reactant at surface, p_A0 (Pa)", "Partial pressure of inert gas, p_I (Pa)", "Molar mass of reactant A, M_A (kg/mol)", "Molar mass of inert, M_I (kg/mol)","Molecular diameter of reactant A, d_A (m)","Molecular diameter of inert, d_I (m)", "Atomic mass of metal, M_M (kg/mol)", "Molar mass of metal oxide deposited, M_MOx (kg/mol)"],
    "Value": [Runcode, ScriptN, Rp, Por, Tau, dpore, SA, VP, DenS, c, q, Pd, t_pulse, T, pA0, pI, MA, MI, dA, dI, MM, MMOx]
}
param_df = pd.DataFrame(param_data)

# Calculated parameter values
calcvalues_data = {
    "Calculated Values": ["Volumetric reactant number density, n_A0 (1/m3)", "Thermal velocity of molecule A, v_A (m/s)", "Knudsen diffusion coefficient, D_Kn (m2/s)", "Molecular diffusion coefficient, D_A (m2/s)", "Effective diffusion coefficient, D_eff (m2/s)", "Aspect ratio as R_p/d_pore, AR (-)", "Knudsen number, Kn (-)", "Thiele modulus, h_T", "Exposure, P_A0*t (Pa-s)", "Metal weight percent, w (%)", "Vol.fraction coated, V_coated/V_sphere (%)", "Max. (possible) metal weight percent, w_max (%)", "Penetration depth (Rp-r) at theta=0.5, d_50% (m)", "Normalized Penetration depth (R_p-r/R_p) at theta=0.5, (d/R_p)_50% (-)", "Slope 'theta/(R_p-r/R_p)' at theta=0.5, (-)", "Diff. between the two data points (y-axis) around theta=0.5 (%)", "Spatial points", "Temporal points"],
    "Value": [nA0, vA, DKn, DA, D, Rp/dpore, Knudsen, Thiele, pA0*t_pulse, wtper, Frac*100, wtfull, PD50, PD50rel, slopePD50rel, ThetaDiff, N, Ntime]
}
calcvalues_df = pd.DataFrame(calcvalues_data)

# Profiles data (Rp-r, pA, theta, nA, )
data_df = pd.DataFrame({
    "Rp-r (µm)": (Rp-r)*1e6,
    "pA (Pa)": (R*T/NAv)*nA_sol[:, -1],
    "theta (-)": theta_sol[:, -1],
    "nA (m-3)": nA_sol[:, -1],
    "Vcoated/Vsphere (-)": FracVol
})

# --- Exporting to Excel ---
run_label = Runcode
filename1 = f"{run_label}_simulation_results_v03.xlsx"
with pd.ExcelWriter(filename1, engine="openpyxl") as writer:
    param_df.to_excel(writer, sheet_name="Results", index=False, startrow=0, startcol=0)
    calcvalues_df.to_excel(writer, sheet_name="Results", index=False, startrow=24, startcol=0)
    data_df.to_excel(writer, sheet_name="Results", index=False, startrow=0, startcol=4)

# ---- Figure of cross-sectional section ---------------

# Creating a 2D grid (X, Y) to represent the circular cross-section
grid_size = 300  # resolution of the image
x = np.linspace(-Rp, Rp, grid_size)
y = np.linspace(-Rp, Rp, grid_size)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)  # radial distance from center
# Interpolate theta values onto the 2D grid
# Values outside the sphere (R > Rp) are masked
from scipy.interpolate import interp1d
theta_interp = interp1d(r, theta_sol[:, -1], bounds_error=False, fill_value=0)
Z = theta_interp(R)
# Masking values outside the sphere
mask = R > Rp
Z_masked = np.ma.array(Z, mask=mask)
# Colormap going from black to yellow 
from matplotlib.colors import LinearSegmentedColormap
colors = [(0, 0, 0), (1, 1, 0)]  # black to yellow (RGB)
custom_cmap = LinearSegmentedColormap.from_list("black_yellow", colors)
# --- Plot ---
fig, ax = plt.subplots(figsize=(6, 5))
cmp = ax.pcolormesh(X * 1e6, Y * 1e6, Z_masked, shading='auto', cmap=custom_cmap)
# Set colorbar
cmp.set_clim(0, 1)  # <-- set color limits here
cb = plt.colorbar(cmp, label='θ (theta)')
# Set dark gray for the outside of the sphere
cmp.cmap.set_bad(color='dimgray')  # outside region
# Final plot settings
ax.set_aspect('equal')
ax.set_xlabel('r: distance from centre of sphere [μm]')
ax.set_ylabel('r: distance from centre of sphere [μm]')
ax.set_title('Radial Profile of θ in a Sphere')
plt.tight_layout()
plt.show()

# --- Final message ---
print(f"Excel with results saved as {run_label}_simulation_results_v03.xlsx")

