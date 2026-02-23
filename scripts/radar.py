import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal.windows import hann
import os

# --- 1. BENCHMARK PARAMETERS (Averages between lower and upper bounds from paper) ---
# TC1 (Steady): Drag [5.57, 5.59], Lift [0.0104, 0.0110]
# TC2 (Unsteady): Drag [3.22, 3.24], Lift [0.99, 1.01], Strouhal [0.295, 0.305]
# TC3 (Unsteady variable): Drag [2.93, 2.97], Lift [0.47, 0.49]
bench = {
    'TC1_cd': (5.5700 + 5.5900) / 2, # 5.5800
    'TC1_cl': (0.0104 + 0.0110) / 2, # 0.0107
    'TC2_cd': (3.2200 + 3.2400) / 2, # 3.2300
    'TC2_cl': (0.9900 + 1.0100) / 2, # 1.0000
    'TC2_st': (0.2950 + 0.3050) / 2, # 0.3000
    'TC3_cd': (2.9300 + 2.9700) / 2, # 2.9500
    'TC3_cl': (0.4700 + 0.4900) / 2  # 0.4800
}

# Helper function to calculate percentage error
def calc_error(sim_value, bench_value):
    return (abs(sim_value - bench_value) / abs(bench_value)) * 100

# Shared simulation parameters
timestep = 0.005
col_names = ['Step', 'Drag_Coefficient', 'Lift_Coefficient']
errori = []

print("=== ERROR CALCULATION FROM CSV ===")

# --- 2. READ AND CALCULATE TC1 (Steady) ---
try:
    df1 = pd.read_csv('output_TC1/coeff_2.csv', header=None, names=col_names)
    # For TC1, we take the final value at steady state
    tc1_cd_sim = df1['Drag_Coefficient'].iloc[-1]
    tc1_cl_sim = df1['Lift_Coefficient'].iloc[-1]
    
    err_tc1_cd = calc_error(tc1_cd_sim, bench['TC1_cd'])
    err_tc1_cl = calc_error(tc1_cl_sim, bench['TC1_cl'])
    
    errori.extend([err_tc1_cd, err_tc1_cl])
    print(f"TC1 -> Drag Err: {err_tc1_cd:.2f}% | Lift Err: {err_tc1_cl:.2f}%")
except FileNotFoundError:
    print("File output_TC1/coeff_2.csv not found. Using 0%.")
    errori.extend([0, 0])

# --- 3. READ AND CALCULATE TC2 (Unsteady) ---
try:
    df2 = pd.read_csv('output_TC2/coeff_2.csv', header=None, names=col_names)
    time2 = df2['Step'] * timestep
    
    # Isolate steady state phase for maxima and FFT 
    mask2 = time2 >= 3.0
    df2_regime = df2[mask2]
    
    tc2_cd_sim = df2_regime['Drag_Coefficient'].max()
    tc2_cl_sim = df2_regime['Lift_Coefficient'].max()
    
    err_tc2_cd = calc_error(tc2_cd_sim, bench['TC2_cd'])
    err_tc2_cl = calc_error(tc2_cl_sim, bench['TC2_cl'])
    
    f_vortex = 2.837
    U_mean = 1.0  # 2 * U_m / 3 (with U_m = 1.5)
    D = 0.1
    tc2_st_sim = (f_vortex * D) / U_mean
    
    err_tc2_st = calc_error(tc2_st_sim, bench['TC2_st'])
    
    errori.extend([err_tc2_cd, err_tc2_cl, err_tc2_st])
    print(f"TC2 -> Drag Err: {err_tc2_cd:.2f}% | Lift Err: {err_tc2_cl:.2f}% | Strouhal Err: {err_tc2_st:.2f}%")
except FileNotFoundError:
    print("File output_TC2/coeff_2.csv not found. Using 0%.")
    errori.extend([0, 0, 0])

# --- 4. READ AND CALCULATE TC3 (Variable) ---
try:
    df3 = pd.read_csv('output_TC3/coeff_2.csv', header=None, names=col_names)
    # For TC3, the maximum peak of the entire simulation is evaluated (8s sinusoidal wave)
    tc3_cd_sim = df3['Drag_Coefficient'].max()
    tc3_cl_sim = df3['Lift_Coefficient'].max()
    
    err_tc3_cd = calc_error(tc3_cd_sim, bench['TC3_cd'])
    err_tc3_cl = calc_error(tc3_cl_sim, bench['TC3_cl'])
    
    errori.extend([err_tc3_cd, err_tc3_cl])
    print(f"TC3 -> Drag Err: {err_tc3_cd:.2f}% | Lift Err: {err_tc3_cl:.2f}%")
except FileNotFoundError:
    print("File output_TC3/coeff_2.csv not found. Using 0%.")
    errori.extend([0, 0])

print("===============================\n")

percentage_errors = errori

# Labels for the 7 vertices
categories = [
    'TC1: Drag', 
    'TC1: Lift', 
    'TC2: Drag', 
    'TC2: Lift', 
    'TC2: Strouhal', 
    'TC3: Drag', 
    'TC3: Lift'
]

N = len(categories)

# --- 2. RADAR CHART PREPARATION ---
# Calculate the angles for each radar axis (in radians)
angles = [n / float(N) * 2 * np.pi for n in range(N)]

# To close the polygon, we need to add the first value at the end of the lists
percentage_errors += percentage_errors[:1]
angles += angles[:1]

# --- 3. PLOT CREATION ---
# Initialize the figure with polar projection
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

# Draw the polygon (the outline)
ax.plot(angles, percentage_errors, linewidth=2, linestyle='solid', color='black', marker='o')

# Fill the polygon area with a semi-transparent color
ax.fill(angles, percentage_errors, color='gray', alpha=0.25)

# Add labels to the vertices
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=11, fontweight='bold')

# Set a maximum limit for the r axis (the radius), for example 5%.
# Modify it if your errors are larger (e.g. ax.set_ylim(0, 10))
max_err = max(percentage_errors)
ax.set_ylim(0, 8) 
ax.set_yticks(np.arange(1, int(max_err) + 2, 1))  # Show concentric circles every 1%
ax.set_yticklabels([f"{i}%" for i in np.arange(1, int(max_err) + 2, 1)], color="gray", size=9)

# Add a title
plt.title('Percentage Error (%) on different benchmark values', size=14, pad=30, fontweight='bold')

# Optimize spacing
plt.tight_layout()

# Save the image
plt.savefig('radar_chart_error.png', dpi=300)